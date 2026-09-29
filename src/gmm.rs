//! One-dimensional two-component Gaussian mixture by EM, and the positivity threshold
//! `state_analysis.py` derives from it.
//!
//! sklearn's `GaussianMixture(2, random_state=seed, n_init=1)` initialises from k-means and runs
//! EM to a log-likelihood tolerance of 1e-3 with a covariance floor of 1e-6. This does the same
//! (k-means++ style seeding from the two extreme quantiles, then Lloyd, then EM) so the fitted
//! means and standard deviations agree to the precision the threshold rule needs, not bit for
//! bit: the rule only asks whether the means are at least `min_separation` apart, and then for
//! their midpoint or `mean_low + tail_sd * sd_low`.

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fit {
    /// Sorted: component 0 is the lower mean.
    pub mean: [f64; 2],
    pub sd: [f64; 2],
    pub weight: [f64; 2],
    pub iterations: usize,
    pub log_likelihood: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Threshold {
    pub fit: Fit,
    pub threshold: f64,
    /// `mean[1] - mean[0]`.
    pub separation: f64,
    pub bimodal: bool,
}

const REG_COVAR: f64 = 1e-6;
const TOL: f64 = 1e-3;

/// Fit the two-component mixture. `x` must have at least two distinct values.
pub fn fit(x: &[f64], max_iter: usize) -> Fit {
    let n = x.len().max(1) as f64;
    // k-means init from the 10th and 90th percentiles, then Lloyd iterations.
    let mut sorted = x.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let q = |p: f64| sorted[((sorted.len() - 1) as f64 * p).round() as usize];
    let mut c = [q(0.1), q(0.9)];
    if c[0] == c[1] {
        c = [sorted[0], sorted[sorted.len() - 1]];
    }
    for _ in 0..20 {
        let (mut s0, mut n0, mut s1, mut n1) = (0.0, 0.0, 0.0, 0.0);
        for &v in x {
            if (v - c[0]).abs() <= (v - c[1]).abs() {
                s0 += v;
                n0 += 1.0;
            } else {
                s1 += v;
                n1 += 1.0;
            }
        }
        let next = [
            if n0 > 0.0 { s0 / n0 } else { c[0] },
            if n1 > 0.0 { s1 / n1 } else { c[1] },
        ];
        if next == c {
            break;
        }
        c = next;
    }
    // Initial responsibilities from the hard assignment.
    let mut mean = c;
    let mut var = [0.0f64; 2];
    let mut w = [0.0f64; 2];
    for &v in x {
        let k = if (v - c[0]).abs() <= (v - c[1]).abs() {
            0
        } else {
            1
        };
        var[k] += (v - mean[k]).powi(2);
        w[k] += 1.0;
    }
    for k in 0..2 {
        var[k] = if w[k] > 0.0 { var[k] / w[k] } else { 1.0 } + REG_COVAR;
        w[k] = (w[k] / n).max(1e-12);
    }
    let mut prev_ll = f64::NEG_INFINITY;
    let mut iterations = 0;
    let mut ll = prev_ll;
    let mut resp = vec![0.0f64; x.len()];
    while iterations < max_iter {
        iterations += 1;
        // E step: responsibility of component 1 for each point, and the log-likelihood.
        ll = 0.0;
        for (i, &v) in x.iter().enumerate() {
            let l0 = w[0].ln()
                - 0.5 * (2.0 * std::f64::consts::PI * var[0]).ln()
                - (v - mean[0]).powi(2) / (2.0 * var[0]);
            let l1 = w[1].ln()
                - 0.5 * (2.0 * std::f64::consts::PI * var[1]).ln()
                - (v - mean[1]).powi(2) / (2.0 * var[1]);
            let m = l0.max(l1);
            let lse = m + ((l0 - m).exp() + (l1 - m).exp()).ln();
            ll += lse;
            resp[i] = (l1 - lse).exp();
        }
        ll /= n;
        // M step.
        let mut s = [0.0f64; 2];
        let mut nk = [0.0f64; 2];
        for (i, &v) in x.iter().enumerate() {
            let r1 = resp[i];
            s[0] += (1.0 - r1) * v;
            s[1] += r1 * v;
            nk[0] += 1.0 - r1;
            nk[1] += r1;
        }
        for k in 0..2 {
            nk[k] = nk[k].max(1e-12);
            mean[k] = s[k] / nk[k];
        }
        let mut ss = [0.0f64; 2];
        for (i, &v) in x.iter().enumerate() {
            let r1 = resp[i];
            ss[0] += (1.0 - r1) * (v - mean[0]).powi(2);
            ss[1] += r1 * (v - mean[1]).powi(2);
        }
        for k in 0..2 {
            var[k] = ss[k] / nk[k] + REG_COVAR;
            w[k] = nk[k] / n;
        }
        if (ll - prev_ll).abs() < TOL {
            break;
        }
        prev_ll = ll;
    }
    let mut out = Fit {
        mean,
        sd: [var[0].sqrt(), var[1].sqrt()],
        weight: w,
        iterations,
        log_likelihood: ll,
    };
    if out.mean[0] > out.mean[1] {
        out.mean.swap(0, 1);
        out.sd.swap(0, 1);
        out.weight.swap(0, 1);
    }
    out
}

/// The threshold rule.
pub fn threshold(fit: Fit, min_separation: f64, tail_sd: f64) -> Threshold {
    let separation = fit.mean[1] - fit.mean[0];
    let bimodal = separation >= min_separation;
    let threshold = if bimodal {
        (fit.mean[0] + fit.mean[1]) / 2.0
    } else {
        fit.mean[0] + tail_sd * fit.sd[0]
    };
    Threshold {
        fit,
        threshold,
        separation,
        bimodal,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_distr::{Distribution, Normal};

    fn draw(rng: &mut rand_chacha::ChaCha8Rng, mean: f64, sd: f64, n: usize) -> Vec<f64> {
        let d = Normal::new(mean, sd).unwrap();
        (0..n).map(|_| d.sample(rng)).collect()
    }

    #[test]
    fn recovers_a_clear_bimodal_mixture() {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(1);
        let mut x = draw(&mut rng, 0.2, 0.3, 30_000);
        x.extend(draw(&mut rng, 2.5, 0.5, 10_000));
        let f = fit(&x, 100);
        assert!((f.mean[0] - 0.2).abs() < 0.03, "{f:?}");
        assert!((f.mean[1] - 2.5).abs() < 0.05, "{f:?}");
        assert!(
            (f.sd[0] - 0.3).abs() < 0.03 && (f.sd[1] - 0.5).abs() < 0.05,
            "{f:?}"
        );
        assert!((f.weight[0] - 0.75).abs() < 0.02, "{f:?}");
        let t = threshold(f, 0.5, 2.0);
        assert!(t.bimodal);
        assert!((t.threshold - 1.35).abs() < 0.05, "{t:?}");
    }

    #[test]
    fn unimodal_data_uses_the_tail_rule() {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(2);
        let x = draw(&mut rng, 0.5, 0.2, 20_000);
        let f = fit(&x, 100);
        let t = threshold(f, 0.5, 2.0);
        assert!(!t.bimodal, "{t:?}");
        // both components sit inside one hump; the threshold is the lower mean + 2 sd
        assert!(t.threshold > 0.5 && t.threshold < 1.2, "{t:?}");
    }

    #[test]
    fn components_are_sorted_by_mean() {
        let x: Vec<f64> = (0..1000)
            .map(|i| {
                if i % 4 == 0 {
                    3.0 + (i % 7) as f64 * 0.01
                } else {
                    (i % 5) as f64 * 0.01
                }
            })
            .collect();
        let f = fit(&x, 50);
        assert!(f.mean[0] < f.mean[1]);
    }
}
