# Experiment 4j: Results — separating wavelength from speed
Gameboard Universe. Run October 12, 2026: 148 GitHub jobs in one workflow run (3 h 46 min). Design: `DESIGN_exp4j.md`, registered before any run. The scoring script, committed before any run, is used unchanged.

## Short answer
- **The glass's extra slowdown fades as about 1/λ^2.2 at fixed speed**: faster than an ordinary lattice's 1/λ², with no sign of a floor.
- **At a fixed wavelength, slower (heavier) particles are slowed more**, by about (γβ)^−0.5.
- **Scattering depends on wavelength only.** The mean free path grows as λ^3.67 with no bend over 17.8–50.3 bonds, and is the same at all three speeds.
- **One setup check failed**, so by the registered rule Q1–Q3 are labeled unreliable. A post hoc check shows the failure was a mistake in how I set up the comparison, not a box-size effect (see below).

## Registered verdict

| # | Criterion | Result |
|---|---|---|
| 0a | Code check and regression job | **PASS**: identical to the earlier engineering test (difference 0) on both the old and new code paths |
| 0b | Diamond controls (21 cells) | **PASS**: 21 of 21. No scattering; particle-only and whole-box speeds equal |
| 0c | Box check vs 4i / 4i-R at the three η = 0.15 cells | **FAIL**: all 6 scattering comparisons pass (≤ 1.3 SE), and all 3 s3 speed comparisons pass (≤ 1.9 SE). All 3 s4 speed comparisons fail (6.4, 7.9 and 9.0 SE). Rule: Q1–Q3 are labeled unreliable |
| **Q1** | How fast the excess shortfall fades (fixed speed) | **q = 2.200 ± 0.013**, χ²/dof 1.15. **"Fades faster than an ordinary lattice"** (unreliable label per 0c) |
| Q2 | Does the scattering law bend? | p = 3.673 ± 0.017 at 29.9 bonds; b = 0.026 ± 0.048; χ²/dof 0.48. **"No bend detected."** Slope at 42.3 bonds is 3.65 ± 0.05, so 4 is excluded (unreliable label per 0c) |
| Q3 | Speed at fixed wavelength | Shortfall: s = −0.507 ± 0.012, **"depends on speed"**. Scattering: s′ = −0.024 ± 0.014, **"no speed dependence detected"** (unreliable label per 0c) |
| — | Clocks within 3%; total conserved | **PASS**: −1.93% to −0.10%; ≤ 2.3×10⁻¹¹ |

## Why the box check failed (post hoc)
- **For s4, the registered comparison paired runs that were not the same.**
  - 4i-R kept 4i's wavelengths on s4 by raising the launch speed: 0.8045c, 0.6917c and 0.5607c, at η = 0.15.
  - 4j keeps the launch speed and raises η by 1.6% instead (0.8000c, 0.6860c, 0.5547c at η = 0.1524).
  - So the s4 settled speeds were bound to differ, by about the launch difference (0.5–1.1%). I should have seen this when designing the check.
- **Compared fairly** (settled speed as a fraction of each run's own launch speed), all six cells agree:

| Cell | Glass | 4j (3×3×3) | 4i / 4i-R (2×2×2) | Difference |
|---|---|---|---|---|
| 0.800c, 17.8 bonds | s3 | 0.93002 | 0.93187 | 1.9 SE |
| 0.800c, 17.8 bonds | s4 | 0.92790 | 0.93008 | 1.9 SE |
| 0.686c, 25.2 bonds | s3 | 0.95950 | 0.95841 | 1.1 SE |
| 0.686c, 25.2 bonds | s4 | 0.95795 | 0.95800 | 0.0 SE |
| 0.555c, 35.6 bonds | s3 | 0.97562 | 0.97467 | 0.7 SE |
| 0.555c, 35.6 bonds | s4 | 0.97411 | 0.97331 | 0.7 SE |

- **The scattering comparisons, which don't depend on this, all passed.**
- **My reading:** the bigger box did not change the results. The registered verdict still stands as FAIL, and Q1–Q3 carry the "unreliable" label.

## The numbers (glass, mean of 6 runs per cell)

| Wavelength | δ at 0.800c | δ at 0.686c | δ at 0.555c | Mean free path (all speeds) |
|---|---|---|---|---|
| 17.8 bonds | 6.24 ± 0.08% | 7.26 ± 0.07% | 8.83 ± 0.07% | 570–585 bonds |
| 21.2 | 4.22 ± 0.04% | 4.89 ± 0.05% | 6.01 ± 0.06% | 1,069–1,090 |
| 25.2 | 2.91 ± 0.04% | 3.37 ± 0.05% | 4.16 ± 0.06% | 2,069–2,099 |
| 29.9 | 1.99 ± 0.04% | 2.28 ± 0.05% | 2.84 ± 0.07% | 3,899–3,942 |
| 35.6 | 1.36 ± 0.03% | 1.56 ± 0.05% | 1.94 ± 0.07% | 7,290–7,455 |
| 42.3 | 0.93 ± 0.03% | 1.05 ± 0.03% | 1.33 ± 0.06% | 14,044–14,092 |
| 50.3 | 0.64 ± 0.02% | 0.72 ± 0.02% | 0.91 ± 0.04% | 24,659–25,353 |

- **δ** is the shortfall relative to the computed packet speed: the registered "excess shortfall", with the packet-shape part removed.
- **Raw shortfalls** relative to launch are 0.6–0.9 percentage points higher: from 9.4% at 0.555c and 17.8 bonds down to 1.5% at 50.3 bonds.
- **Scattering at 25.2 bonds matches 4i:** 1.742–1.768×10⁻⁴ per Å here, vs 1.742×10⁻⁴ for 4i.

## Checks of how solid Q1 is (post hoc)
- **Step by step, the slope of δ is steady:** 2.12–2.28 at every step and every speed, with no drift toward 2 at the long end.
  - A fit with an ordinary k² term plus a k⁴ correction fits worse (χ²/dof 2.06 vs 1.15).
- **The exponent depends somewhat on the reference packet speed.**
  - Using the packet's actual tracker-weighted spectrum instead of the registered Gaussian gives q = 2.12 ± 0.01.
  - Shifting the reference by ±0.1% of c moves q between 2.11 and 2.30.
  - The diamond crystal's δ goes slightly negative (−0.06%) at its longest wavelength at 0.8c, which suggests the reference is uncertain by about 0.1%.
  - So the honest statement is **q ≈ 2.1–2.3**: at least as fast as an ordinary lattice, and probably a little faster.
- **No floor:** a fit with a constant floor gives B = +0.006 ± 0.034%, consistent with zero.
- **The two glasses agree:** s4's δ averages 5% above s3's.
- **Drag over the run:** −0.12% per 100 ticks on average (halves −0.20%). This repeats the small residual slowdown seen before.

## Reading
- **The slowdown is mostly the glass's disorder, not lattice-ness.** The diamond crystal, with the same packets, is slowed 5–30× less at similar wavenumbers (0.17–2.4% at 0.555c, falling as about 1/λ^2.4).
- **At fixed speed, both wavelength and mass shrink together** when everything is scaled up toward real particle sizes, and δ falls as about 1/λ^2.2. In 4i, mass was held fixed while wavelength changed, which is why it looked like λ^−1.4.
  - Written in terms of k and η separately, the fit is δ ∝ k^1.7 η^0.5.
- **What this does and doesn't say:**
  - On this range (2.8× in wavelength), the slowdown fades at least as fast as an ordinary lattice's.
  - Taken at face value, a k^2.2 slowdown would be far too small to see at real particle sizes.
  - But the range is short, so this is a local slope, not a law at 10¹⁶ spacings.
- **Scattering:** 1/ℓ ∝ k^3.67 over 17.8–50.3 bonds, with no bend and no dependence on speed or mass.
  - The rising step-by-step slopes in 4i (3.4 → 4.1) were not confirmed. With 6 runs per cell, the slope is flat at about 3.67 across the whole range.
  - It is still not the Rayleigh value of 4.

## Files
- **Code:** `run_4j.py`, `run_big.py` and `run_4h2.py` (per-run η), `score_4j.py` (registered), `fig4j.py`.
- **Results:** `res_4j/` (148 jobs; also `results/exp4j/` on GitHub), `score_4j.json`, `exp4j_results.png`.
