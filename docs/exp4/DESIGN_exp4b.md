# Experiment 4b: does smoothing the glass reduce the short-wavelength clock error? (pre-registered)
Gameboard Universe. Written October 8, 2026, before any 4b runs. Rule, particles and measurements are exactly as in Experiment 4, including Amendment 1.

## Question
- In Experiment 4, raw glass clocks failed at 0.7c and 0.8c.
- The error followed err% ≈ a·(k·b)², with a = −63 on raw glass (s3) and a = −5.5 on diamond (b = 2.35 Å).
- Does smoothing the bond angles (30° → 20° spread) make |a| smaller?

## Networks (all 1,728,000 nodes, box 325.6 Å)
1. **Smoothed, tiled:** `s1_final` (216k nodes) repeated 2×2×2. The same glass passed Experiment 2.
2. **Raw, tiled:** `s1_raw`, the unsmoothed parent of the same sample, repeated 2×2×2 the same way.
   - This is the paired control: the only difference between 1 and 2 is the smoothing.
3. **Raw, untiled:** `s3`, already run in Experiment 4. Used only to check for tiling artifacts.

## Runs
- Same as Experiment 4: η = 0.15, packets 40 Å wide, 250 ticks.
- Speeds: rest, 0.3c, 0.5c, 0.7c, 0.8c.
- Directions: [100] and [1,2,3].
- The classical check is not repeated.

## Pre-registered checks
1. **Tiling check:** |a| for raw tiled agrees with raw untiled (−63) within a factor of 1.5. If not, the tiled results are reported as unreliable.
2. **Main question:** compare a_smooth with a_raw,tiled.
   - **"Smoothing helps":** |a_smooth| ≤ 0.5 × |a_raw,tiled|.
   - **"Smoothing doesn't matter":** the two agree within ±25%.
   - **Anything in between:** inconclusive.
3. **Experiment 4 criteria, reported per network** (unchanged):
   - clock within 3% of relativity at each speed;
   - speed change between halves < 1%;
   - directions agree within 1%;
   - total conserved to 10⁻⁹;
   - scattered fraction ≤ 1%.
4. **Fit method:** a = least-squares slope of err% against (k·b)² through the origin, over all 8 moving runs. The scatter is reported.

## Limits stated in advance
- Tiling makes the randomness repeat every 163 Å. A wave can then pick up a slight ordered-crystal character, which might lower |a| artificially. Check 1 is there to catch this.
- One parent sample.

Any outcome is reported as found.

## Amendment 1 (logged after the networks loaded, before any 4b particle result was seen)
- **The raw parent differs from the raw glass in Experiment 4.** `s1_raw` is the network before any relaxation: bond length 2.70 ± 0.50 Å. The Experiment 4 raw glass, `s3`, was relaxed: 2.35 ± 0.27 Å. The smoothed `s1_final` measures 2.385 ± 0.22 Å.
  - So raw tiled vs smoothed tiled compares unrelaxed against relaxed-and-smoothed. It does not isolate bond angles alone.
  - So check 1 (tiling) compares two different kinds of network. If it fails, the failure can't be pinned on tiling.
- **What will be reported:**
  - The **primary comparison is smoothed tiled vs s3** (relaxed raw, untiled). Both are relaxed, and they differ in angle spread (20° vs 30°) and in tiling.
  - **Raw tiled vs smoothed tiled** is reported as secondary.
  - The thresholds are unchanged: ≤ 0.5× means "helps", within ±25% means "doesn't matter".
- **Bond length in the fit:** b = 2.35 Å as registered. A second fit using each network's own mean bond is also reported.

## Correction to Amendment 1 (logged after the runs)
- Amendment 1 gave s3's bond length as 2.35 ± 0.27 Å. That was wrong: I took the nominal 2.35 Å as the mean, and the 0.27 in s3's stored metadata is a different quantity.
- Measured directly, s3's mean bond is **2.735 Å**, about the same as s1_raw's 2.704 Å.
- So the raw tiled control is a closer match to s3 than Amendment 1 said. Check 1 (tiling) is a valid test after all.
- Both comparisons are reported. Neither threshold was changed.
