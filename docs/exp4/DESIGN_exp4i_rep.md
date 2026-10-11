# Experiment 4i-R: independent-glass replication (pre-registered)
Gameboard Universe. Written October 10, 2026, before the new glass was built.

## Question
Experiment 4i found, on one glass sample (s3):
- mean free path ℓ from 561 to 7,324 bonds over wavelengths of 17.8–35.6 bonds;
- a single power law, 1/ℓ ∝ k₀^3.74 ± 0.04;
- directions agreeing except at the longest wavelength.

Are these properties of this kind of glass, or of that one sample?

## The new sample
- **Glass s4:** built by exactly the same recipe as s3, with only the random seed changed.
  - Command: `make_raw_big.py 60 4`: a random 4-regular network of 1,728,000 nodes in a 326 Å box, relaxed with the same settings.
  - Seeds 1–3 are already used; s4 has never been run.
- **Built on GitHub** (workflow `build_glass4`), stored in `networks/`.
- **Before any particle run**, its quality numbers are reported against s3: bond rms, angle mean and rms, S(111) and energy per atom.

## Its own speed limit c
- **Measured** by the 4c method (diag_c3), on s4 untiled:
  - the node totals' wave equation with η = 0.05;
  - c is read from the frequency of the longest wave fitting the box along [100], the definition that gave s3's 1.5495.
- **Also reported:** the n = 2, n = 3 and [1,2,3] waves.

## Runs (31 GitHub jobs; code exactly as in 4i, `run_4i.py`)
- **Glass s4, tiled 2×2×2** (652 Å box): 5 points × 2 directions × 3 starts = 30 runs, 600 ticks each.
- **Same wavelengths and widths as 4i.** k₀ and σ are identical to 4i point by point.
  - The launch speed is the one that gives that k₀ under s4's measured c: γv = k₀ c_s4 / η.
  - If c_s4 = c_s3, the launch speeds equal 4i's.
- **Pipeline check:** 1 diamond job (point 5, [100]), which must reproduce 4i's diamond point-5 result to 10⁻¹².
  - 4i's five diamond controls ran on the same deterministic network and code, so they are not repeated.
- **Measures and fit** exactly as in 4i: 1/τ from −ln f over ticks 100–600; ℓ = v_p τ; weighted power-law fit of ln(1/ℓ) against ln k₀.

## Pre-registered criteria
0. **Setup:**
   - (a) the pipeline check reproduces 4i to 10⁻¹²;
   - (b) c_s4 has a direction spread ≤ 0.5% ([100] n = 1 vs [1,2,3]), as for s3;
   - (c) s4's quality numbers are reported (not scored).

   If (a) fails, nothing else is judged.
1. **Falls with wavelength:** 1/ℓ is lower at each step to a longer wavelength, by more than 2 SE, in all 4 steps.
2. **Single power law:** χ²/dof < 3.
3. **Exponent replicates:** |p_s4 − 3.74| < 2 × the combined SE.
4. **Size replicates:** the fitted 1/ℓ at the middle wavelength (point 3, k₀ = 0.09127) is within 10% of 4i's.
   - Also reported, not scored: point-by-point ratios and their combined-SE tests.
   - Starts within one sample don't capture sample-to-sample variation, so a 10% band is used rather than SE.
5. **Directions:** [100] and [1,2,3] agree within 2 SE at each point.
   - This tests whether 4i's point-5 mismatch belongs to the s3 sample.
6. **Clocks:** within 3% of relativity (using s4's c) in every glass run. **Total** conserved to 10⁻⁹.

## Reported, not scored
- Settled particle speed vs launch at each point, compared with 4i's shortfall pattern (6.8% shrinking to 2.5%).
- Drag, halves, the share cut off by the cone/shell.
- A two-sample combined fit, giving the exponent and how much the size varies between samples.

## Limits stated in advance
- Two samples is the minimum for a replication. It shows whether s3 was unusual, not the full spread between samples.
- Both samples share the recipe: the same box, density and relaxation. A different recipe (for example, smoothed glass) is a separate question.

Any outcome is reported as found.
