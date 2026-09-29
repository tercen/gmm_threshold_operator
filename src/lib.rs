//! gmm_threshold_operator — a positivity threshold per marker from a two-component Gaussian
//! mixture, the rule of `state_analysis.py` in the CYTOSHRINK pipeline.
//!
//! Rows are channels, columns are cells, y is the value. The result is one row per cell with
//! its SOM node and its metacluster, plus a table describing the map.
//!
//! **The shape of the problem is a transpose.** The crosstab arrives as scattered
//! `(.ri, .ci, .y)` triples and a map needs each cell's whole vector across channels. There is
//! no order to rely on — R's own client scatters by index rather than assuming one — so the
//! operator gathers the matrix itself, `n_cells × n_channels` of `f32`, and the memory model
//! declares that cost rather than hiding it. `f32` because the R pipeline is `f32` anyway
//! (flowCore stores expressions as 32-bit floats), so the second half of the mantissa was never
//! real.
//!
//! One pass to gather, one map, one assignment, one write.
pub mod context;
pub mod gmm;
pub mod input;
pub mod output;
pub mod pagecache;
pub mod progress;
pub mod props;
pub mod tson;
pub mod upload;

use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result, anyhow};
use tercen_rs::context::ContextBase;
use tercen_rs::{DevContext, TercenClient};

use progress::Reporter;
use tson::TsonWriter;

const CHUNK: usize = 1_000_000;

pub fn init_tracing() {
    let filter = tracing_subscriber::EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info"));
    let _ = tracing_subscriber::fmt().with_env_filter(filter).try_init();
}

pub fn require_env(name: &str) -> Result<String> {
    std::env::var(name).map_err(|_| anyhow!("{name} is not set"))
}

pub async fn run(task_id: &str) -> Result<()> {
    tracing::info!("gmm_threshold_operator starting (task_id={task_id})");
    let client = build_client().await?;
    let ctx = context::from_task_id(client, task_id).await?;
    execute(
        &ctx,
        Mode::Production {
            task_id: task_id.to_string(),
        },
    )
    .await
}

pub async fn run_dev(workflow_id: &str, step_id: &str) -> Result<()> {
    tracing::info!("gmm_threshold_operator starting in dev mode ({workflow_id} / {step_id})");
    let client = build_client().await?;
    let ctx = DevContext::from_workflow_step(client, workflow_id, step_id)
        .await
        .map_err(|e| anyhow!("load workflow {workflow_id} / step {step_id}: {e}"))?;
    execute(
        &ctx,
        Mode::Dev {
            workflow_id: workflow_id.to_string(),
            step_id: step_id.to_string(),
        },
    )
    .await
}

enum Mode {
    Production {
        task_id: String,
    },
    Dev {
        workflow_id: String,
        step_id: String,
    },
}

async fn build_client() -> Result<Arc<TercenClient>> {
    let client = TercenClient::from_env()
        .await
        .map_err(|e| anyhow!("connect to Tercen: {e}"))?;
    tracing::info!("connected to Tercen");
    Ok(Arc::new(client))
}

async fn execute(ctx: &ContextBase, mode: Mode) -> Result<()> {
    let t_start = Instant::now();
    tracing::info!(
        workflow = ctx.workflow_id(),
        step = ctx.step_id(),
        namespace = ctx.namespace(),
        "context loaded"
    );
    let rep = match &mode {
        Mode::Production { task_id } => Reporter::spawn(Arc::clone(ctx.client()), task_id.clone()),
        Mode::Dev { .. } => Reporter::silent(),
    };
    let s = props::read(ctx)?;
    tracing::info!(?s, "properties");

    let n_values = input::cell_count(ctx).await?;
    let markers = input::row_labels(ctx).await?;
    let p = markers.len();
    if p == 0 {
        anyhow::bail!("the projection has no rows: put the markers on rows");
    }
    let n_cells = n_values / p;
    tracing::info!(n_values, n_cells, markers = p, "projection");

    rep.at(0, "Reading the crosstab");
    let data = gather(ctx, n_values, n_cells, p, &rep).await?;

    // One seeded subsample of cells, shared by every marker, as state_analysis.py draws it.
    let sub: Vec<usize> = if s.subsample > 0 && n_cells > s.subsample {
        draw(n_cells, s.subsample, s.seed)
    } else {
        (0..n_cells).collect()
    };
    rep.at(progress::FIT.0, "Fitting the mixtures");
    let t = Instant::now();
    let mut threshold = Vec::with_capacity(p);
    let mut mode_s = Vec::with_capacity(p);
    let (mut sep, mut m0, mut s0, mut m1, mut s1, mut w0, mut pct) =
        (vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
    let mut iters = Vec::with_capacity(p);
    for r in 0..p {
        let col = &data[r * n_cells..(r + 1) * n_cells];
        let x: Vec<f64> = sub.iter().map(|&i| col[i] as f64).collect();
        let f = gmm::fit(&x, s.max_iter);
        let th = gmm::threshold(f, s.min_separation, s.tail_sd);
        let positive = col.iter().filter(|v| (**v as f64) > th.threshold).count();
        threshold.push(th.threshold);
        mode_s.push(if th.bimodal { "bimodal" } else { "unimodal" }.to_string());
        sep.push(th.separation);
        m0.push(f.mean[0]);
        s0.push(f.sd[0]);
        m1.push(f.mean[1]);
        s1.push(f.sd[1]);
        w0.push(f.weight[0]);
        pct.push(100.0 * positive as f64 / n_cells.max(1) as f64);
        iters.push(f.iterations as i32);
        tracing::info!(marker = %markers[r], threshold = th.threshold, bimodal = th.bimodal, separation = th.separation, "fitted");
    }
    let n_bimodal = mode_s.iter().filter(|m| *m == "bimodal").count();
    tracing::info!(
        secs = format!("{:.1}", t.elapsed().as_secs_f64()),
        n_bimodal,
        "mixtures fitted"
    );
    rep.info(format!(
        "GMM thresholds: {p} markers on {} of {n_cells} cells; {n_bimodal} bimodal, {} unimodal",
        sub.len(),
        p - n_bimodal
    ));

    rep.at(progress::WRITE.0, "Writing the result");
    let work_root = std::env::temp_dir().join(format!(
        "gmm_threshold_op_{}_{}",
        ctx.workflow_id(),
        ctx.step_id()
    ));
    std::fs::create_dir_all(&work_root)
        .with_context(|| format!("create {}", work_root.display()))?;
    let _guard = TempDirGuard(work_root.clone());
    let result_path = work_root.join("result.tson");
    {
        let f = std::fs::File::create(&result_path)
            .with_context(|| format!("create {}", result_path.display()))?;
        let w = std::io::BufWriter::with_capacity(4 << 20, pagecache::Releasing::new(f, 256 << 20));
        let mut w = TsonWriter::new(w)?;
        let ns = ctx.namespace();
        let name = |c: &str| format!("{ns}.{c}");
        let cols = [
            output::Col::F64(&name("threshold"), threshold),
            output::Col::Str(&name("mode"), mode_s),
            output::Col::F64(&name("separation"), sep),
            output::Col::F64(&name("mean_low"), m0),
            output::Col::F64(&name("sd_low"), s0),
            output::Col::F64(&name("mean_high"), m1),
            output::Col::F64(&name("sd_high"), s1),
            output::Col::F64(&name("weight_low"), w0),
            output::Col::F64(&name("pct_positive"), pct),
            output::Col::I32(&name("cells_fitted"), vec![sub.len() as i32; p]),
            output::Col::I32(&name("em_iterations"), iters),
        ];
        output::write_row_table(&mut w, &table_name(ctx), p, &cols)?;
        output::write_footer(&mut w)?;
    }
    let bytes = std::fs::metadata(&result_path)?.len();
    tracing::info!(bytes, "result written");
    pagecache::release_path(&result_path);

    rep.at(progress::UPLOAD.0, "Uploading the result");
    match mode {
        Mode::Production { task_id } => {
            upload::save_production(ctx, &task_id, &result_path, &rep).await?
        }
        Mode::Dev {
            workflow_id,
            step_id,
        } => {
            let saved = upload::save_dev(ctx, &workflow_id, &step_id, &result_path).await?;
            tracing::info!(task_id = saved.task_id, "dev result saved");
        }
    }
    rep.at(100, "Done");
    tracing::info!(
        total_secs = format!("{:.1}", t_start.elapsed().as_secs_f64()),
        peak_rss_kb = peak_rss_kb().unwrap_or(0),
        "done"
    );
    Ok(())
}

/// `k` distinct indices below `n`, by a seeded shuffle — the same cells for every marker.
fn draw(n: usize, k: usize, seed: u64) -> Vec<usize> {
    use rand::SeedableRng;
    use rand::seq::SliceRandom;
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
    let mut idx: Vec<usize> = (0..n).collect();
    idx.shuffle(&mut rng);
    idx.truncate(k);
    idx.sort_unstable();
    idx
}

/// Gather the crosstab into a cell-by-channel matrix.
///
/// Scattered by `(.ri, .ci)`, so no ordering is assumed.
///
/// **`f32` is a parity decision, not only a memory one.** The R operator does
/// `flowCore::flowFrame(as.matrix(data))`, and flowCore stores expressions as 32-bit floats, so
/// every value it clusters has already been rounded — by up to 2.4e-7. Rounding here too is what
/// makes the two agree; keeping full precision would be more accurate and would not be FlowSOM.
/// It also halves the peak, which is the part the memory model books.
async fn gather(
    ctx: &ContextBase,
    n_values: usize,
    n_cells: usize,
    p: usize,
    rep: &Reporter,
) -> Result<Vec<f32>> {
    let mut data = vec![f32::NAN; n_cells * p];
    let mut seen = 0usize;
    let mut out_of_range = 0usize;
    input::for_each_chunk(ctx, &[".ri", ".ci", ".y"], n_values, CHUNK, |c| {
        for k in 0..c.y.len() {
            let (ri, ci) = (c.ri[k] as usize, c.ci[k] as usize);
            if ri >= p || ci >= n_cells {
                out_of_range += 1;
                continue;
            }
            data[ci + ri * n_cells] = c.y[k] as f32;
        }
        seen += c.len();
        rep.at(
            progress::band(progress::READ, seen, n_values.max(1)),
            format!("Read {seen} of {n_values}"),
        );
        Ok(())
    })
    .await?;
    if out_of_range > 0 {
        anyhow::bail!(
            "{out_of_range} values fell outside the {n_cells} x {p} crosstab the schema \
             described; the projection changed under the run"
        );
    }
    let missing = data.iter().filter(|v| v.is_nan()).count();
    if missing > 0 {
        anyhow::bail!(
            "{missing} of {} cell-channel pairs have no value. a mixture needs a complete matrix: \
             every cell must have every channel. Filter the projection, or fill the gaps \
             upstream.",
            data.len()
        );
    }
    Ok(data)
}

struct TempDirGuard(std::path::PathBuf);
impl Drop for TempDirGuard {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn peak_rss_kb() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1)?.parse().ok())
        })
}

fn table_name(ctx: &ContextBase) -> String {
    format!("{}_{}", ctx.step_id(), ctx.qt_hash())
}
