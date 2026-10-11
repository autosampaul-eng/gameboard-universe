# Experiment 4k: is the direction difference a finite-sample effect? (pre-registered)
Gameboard Universe. Written October 11, 2026, before the new glass was built.

## Question
- In 4i (glass s3) and 4i-R (glass s4), scattering along [100] and [1,2,3] differed by up to 28%. The sign differed between the two samples, and the size grew with wavelength.
- **Proposed explanation (post hoc in 4i-R):** each sample's own random anisotropy. Both boxes are one 326 Å glass repeated 2×2×2, so they contain only one independent 326 Å piece of glass.
- **If so,** a glass with more independent material should be less lopsided. For a random effect averaging over independent volume, the spread should fall as 1/√volume.

## The new sample (s5)
- **Same recipe as s3/s4,** seed 5, but built directly at 13,824,000 nodes in a 652 Å box (ncell 120), not tiled.
  - That is 8× the independent volume of s3/s4.
  - **Prediction if the effect is finite-sample:** the direction spread is 1/√8 ≈ 0.35 of s3/s4's.
- **Engineering changes,** none of which changes the recipe:
  - The builder's pair list is held in a compact form (memory). It is checked to give bit-identical networks to the original builder on a small test size before use.
  - The relaxation is split across GitHub jobs by pausing and resuming its exact state. It is checked to be bit-identical to an uninterrupted run before use.
  - Repair-pass cap 400 (as in 4i-R Amendment 1).
- **Its quality numbers** are reported against s3/s4: energy, bond and angle rms, and the residual force (must be below 10⁻³, as for s3/s4).
- **Its speed limit c** is measured by the 4c method at the same wavenumber that defined s3's c: k = 2π/326 Å, which is the n = 2 wave along [100] in the 652 Å box.
  - Also reported: n = 1, 3, 4 and [1,2,3].

## Runs (code as in 4i; `run_4k.py`)
- **Wavelength points:** 4i points 1 (17.8 bonds, σ = 40) and 5 (35.6 bonds, σ = 80). Same k₀ and σ as 4i; launch speed from each glass's own c.
- **Directions (6):** [100], [010], [001], [1,2,3], [3,1,2] and [2,3,1].
- **Starts:** 2 per direction, the first two of 4f–4i, in units of each network's base box.
- **Glass s5 (untiled, 652 Å):** 2 points × 6 directions × 2 starts = 24 runs.
- **Glass s3 and s4 (tiled 2×2×2, 652 Å):** [100] and [1,2,3] reuse 4i and 4i-R (starts 1 and 2). The other 4 directions are new: 2 × 4 × 2 = 16 runs each, 32 in all.
- **Total:** 56 GitHub jobs. Measures as in 4i: 1/τ over ticks 100–600; ℓ = v_p τ.

## Measure
- **Per network and point:** m_d = mean 1/ℓ over the 2 starts for each of the 6 directions, and its start-to-start SE, se_d.
- **Direction spread:**

\[ A = \frac{\sqrt{\max\left(\mathrm{var}_d(m_d) - \overline{se_d^2},\ 0\right)}}{\overline{m_d}} \]

  - That is the spread across directions, with start-to-start noise removed.
- **Tiled reference:**

\[ A_{tiled} = \sqrt{\left(A_{s3}^2 + A_{s4}^2\right)/2} \]

- **Ratio:** R = A_s5 / A_tiled.

## Pre-registered criteria
0. **Setup:**
   - (a) the lean builder and the resumed relaxation are bit-identical to the originals on a test size;
   - (b) s5's residual force is below 10⁻³;
   - (c) s5's c direction spread is ≤ 0.5%.

   If any fails, Q is not judged.
- **Q (primary, point 5, the long wavelength, where the effect is largest):**
  - **R ≤ 0.6:** "shrinks as a finite-sample effect" (expected 0.35);
  - **R ≥ 0.85:** "does not shrink", meaning the direction difference is not explained by sample size;
  - **otherwise:** "inconclusive".
- **Q1 (secondary, point 1):** the same rule.
- **S (same law):** direction-averaged 1/ℓ on s5 is within 10% of the s3/s4 average, at both points.
- **Clocks:** within 3% in every glass run. **Total** conserved to 10⁻⁹.

## Reported, not scored
- A bootstrap 68% range for R (resampling starts within each direction).
- Per-direction values for all three glasses.
- Whether s5's direction pattern correlates with s3's or s4's (it shouldn't, if the anisotropy is random).
- Settled-speed shortfall on s5.

## Limits stated in advance
- **Six directions** give a rough estimate of a spread. The relative uncertainty of A from 6 directions is about 30%. The thresholds leave room for this: the expected 0.35 vs the cut at 0.6.
- **One untiled sample.** A second would be needed to measure the spread of A itself.
- **The 1/√volume prediction** assumes the anisotropy comes from short-range-correlated random structure.

Any outcome is reported as found.

## Amendment 1 (October 11, 2026, before any s5 relaxation or particle run)
- **What happened:** the first s5 build (GitHub run 38130936935) stopped at the 400-pass repair cap with 2 nodes still one link short. These were the last 2 of 1,254,447 deficient nodes after the greedy pass, and they stayed stuck from about pass 100 to pass 399.
  - The builder's repair moves a missing link randomly, so two lone gaps far apart in a 13.8M-node box rarely meet.
- **Change:** if gaps remain after the repair passes, a final guided step (`bigglass.guided_endgame`) pairs them up.
  - It uses the repair's own move: remove link x–y, add u–x, and the gap passes to y.
  - Half the steps choose x among u's neighbours closest to the partner gap; half are the repair's ordinary random choice.
  - Same constraints as the repair: no duplicate links and no rings of 4 or fewer.
  - It touches only the links along one short path. Networks that the repair completes are unchanged, so the bit-identity test of 0a still holds.
- **Tested before use:**
  - Two gaps 18 Å apart (6k-node box) closed in 54 steps; two gaps 78 Å apart (110k-node box) closed in 219 steps.
  - Both results are symmetric, with no duplicate links and no ring of 4 or fewer anywhere.
- **Chosen instead of a new seed:** the 4i-R s4 build also ended with 2 stuck nodes at 80 passes. At 8× the size, a stuck pair is likely for any seed.
- **Seed and everything else unchanged.** The number of endgame steps will be reported.
