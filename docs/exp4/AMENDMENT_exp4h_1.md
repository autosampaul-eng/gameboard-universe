# Experiment 4h: Amendment 1 (October 10, 2026), before the rerun
**What changed:** the particle-only tracker only.

## Why
- In the first 4h run, diamond control 0c failed: the particle-only speed read 0.11% below the whole-box speed.
- **Cause:** the tracker used the density of the grid-mapped node totals |S|². The conserved quantity is the link probability Σ|φ|². The ratio between the two differs from one wavelength to another within the particle. So the |S|² centre moves at a differently weighted average speed.

## The fix
- Before the inverse transform, the forward-cone spectrum is multiplied by √w(k) / W_cic(k):
  - w(k) is the link probability per unit |S|² of the one-branch plane-wave mode (the same ansatz used to build packets), averaged over 20,000 sample nodes;
  - W_cic(k) is the cloud-in-cell window.
- So the tracked density is the particle's conserved probability, the same quantity the whole-box centre uses.

## Checked on the 216k-derived test networks before the rerun
- **Test diamond (2×2×2 tiled from the 216k diamond), 0.8c:**
  - corrected 0.7910 / 0.7869 against whole-box 0.7910 / 0.7869 ([100] / [1,2,3]): agreement better than 0.01%;
  - the old tracker reproduces its 0.11% low reading (0.7902 / 0.7859).
- **Test glass:** the correction shifts the particle-only speed by +0.08%.
- **At 0.5c (σ = 40 Å)** the tracker is still biased, because the cone/shell trims part of the particle. This amendment does not fix that. The 0.5c speeds stay unscored, as in the original design.

## Not changed
- Runs, criteria (including 0c and its reference values), the scattering measure s(t), and Γ.
- The original tracker is also recorded in every job (`coh_orig`), for comparison.
