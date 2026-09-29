//! Operator properties.
//!
//! The rule is `state_analysis.py`'s (Jamie's pipeline, section 14): per marker, fit a
//! two-component Gaussian mixture to a seeded subsample; if the two means are at least
//! `min_separation` apart the marker is bimodal and the threshold is their midpoint, else it is
//! unimodal and the threshold is the lower mean plus `tail_sd` standard deviations of the lower
//! component.
use anyhow::Result;
use tercen_rs::PropertyReader;
use tercen_rs::context::ContextBase;

#[derive(Debug, Clone)]
pub struct Settings {
    /// Cells drawn (without replacement, seeded) to fit each marker's mixture; 0 = all.
    pub subsample: usize,
    pub seed: u64,
    /// Below this distance between the two means the marker is treated as unimodal.
    pub min_separation: f64,
    /// Unimodal threshold: `mean_low + tail_sd * sd_low`.
    pub tail_sd: f64,
    /// EM iterations at most.
    pub max_iter: usize,
    /// EM starts per marker; the best log-likelihood wins (sklearn's `n_init`).
    pub n_init: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            subsample: 500_000,
            seed: 5024,
            min_separation: 0.5,
            tail_sd: 2.0,
            max_iter: 100,
            n_init: 1,
        }
    }
}

pub fn read(ctx: &ContextBase) -> Result<Settings> {
    let pr = PropertyReader::from_operator_settings(ctx.operator_settings());
    let d = Settings::default();
    let num = |name: &str, dflt: f64| -> Result<f64> {
        let raw = pr.get_string(name, &dflt.to_string());
        raw.trim()
            .parse::<f64>()
            .map_err(|_| anyhow::anyhow!("property '{name}' is not a number: {raw:?}"))
    };
    let s = Settings {
        subsample: num("subsample", d.subsample as f64)?.max(0.0) as usize,
        seed: num("seed", d.seed as f64)?.max(0.0) as u64,
        min_separation: num("min_separation", d.min_separation)?,
        tail_sd: num("tail_sd", d.tail_sd)?,
        max_iter: num("max_iter", d.max_iter as f64)?.max(1.0) as usize,
        n_init: num("n_init", d.n_init as f64)?.max(1.0) as usize,
    };
    Ok(s)
}
