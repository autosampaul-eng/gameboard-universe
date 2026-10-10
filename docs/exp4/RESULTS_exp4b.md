# Experiment 4b: Results — does smoothing the glass help?
Gameboard Universe. Run October 8, 2026. Design: `DESIGN_exp4b.md`, with Amendment 1 and a correction.

> **Correction (Oct 8, Experiment 4c):** the glass speed limit c used here came from a smooth-average formula. The glass's measured c is about 3.5–4% lower. With the measured c, every glass clock in this experiment is within 1.7% of relativity, so all pass. The "(k·b)² law" and "glass 9× worse than crystal" conclusions came from the wrong c and are withdrawn. The registered verdicts below stand as recorded. See `RESULTS_exp4c.md`.

## Verdict
- **Smoothing doesn't matter.** Smoothing the bond angles (30° → 20° spread) does not reduce the clock error at high speed. The pre-registered "doesn't matter" outcome is met under every comparison.
- **Tiling check: PASS.** Repeating a 216k-node glass 2×2×2 gives the same error coefficient as an independent 1.73M-node glass.
- **All three glasses follow one error law.** Measured against each network's own bond length, the coefficient is about −50, roughly 9× worse than the diamond crystal's −5.5.
- So the error comes from random connections (which node links to which), not from bond angles. Experiment 3b found the same for drag.

## Error coefficient a (clock error % ≈ a × (k·b)²)

| Network | Mean bond (Å) | c (Å/tick) | a, b = 2.35 Å (registered) | a, own mean bond | Scatter |
|---|---|---|---|---|---|
| Raw glass s3, untiled | 2.735 | 1.607 | −63.3 | −46.7 | 0.14% |
| Raw glass s1, tiled | 2.704 | 1.588 | −63.2 | −47.7 | 0.17% |
| Smoothed glass s1, tiled | 2.385 | 1.383 | −53.0 | −51.4 | 0.24% |
| Diamond crystal | 2.350 | 1.357 | −5.5 | −5.5 | 0.13% |

## Pre-registered checks
1. **Tiling:** raw tiled ÷ raw untiled = 1.00. The limit was a factor of 1.5, so this passes.
2. **Main question:**
   - Primary comparison, smoothed vs s3: the ratio is 0.84 with the registered b, or 1.10 with each network's own bond. Both fall within ±25%, so **"doesn't matter"**.
   - Secondary comparison, smoothed vs raw tiled: 0.84 or 1.08. Also **"doesn't matter"**.
3. **Experiment 4 criteria on the smoothed glass:**
   - **Clocks:** pass at 0.3c (−0.3%) and 0.5c (−1.2%). Fail at 0.7c (−3.6%, −3.4%) and 0.8c (−6.4%, −5.6%).
   - **Speed change between halves:** fails at 0.7c and 0.8c (+1.7% to +7.1%). The packets speed up, rather than slowing as drag would make them.
   - **Directions:** at 0.8c, the two directions' speeds differ by 2.2% (fail).
   - **Total:** conserved to 2×10⁻¹⁴. Scattered fraction ≤ 0.02%.
4. **Raw tiled glass:** almost identical to s3. Clocks fail at 0.7c (−3.3% and −3.0%) and 0.8c (−5.7% and −5.3%).

## What this means
- **Smoothing does help in one other place.** It helped the shortest-path light clock in Experiment 2. For wave particles under the GBU's own rule, it does nothing.
- **The issue is random connectivity itself.** The crystal's ordered connectivity gives an error 9× smaller at the same wavelength.
- **For the theory:** a glass substrate stays consistent with relativity only in the long-wavelength limit. The deviation grows as (spacing ÷ wavelength)². At any physical spacing anywhere near the Planck scale, this is far too small to measure. It is still a real difference between glass and crystal: a quadratic dispersion term with a coefficient set by the disorder.
- **Glass vs crystal, where things stand:**
  - The crystal is better here, by 9×.
  - At long wavelengths both are direction-independent. In Stage A, c varies by 0.04% across directions on diamond, 0.07–0.09% on s3 and the smoothed glass, and 0.38% on raw tiled glass.
  - Glass's argued advantage is that it has no preferred axes at any scale. That wasn't tested at short wavelengths here, where a crystal's cubic axes would show.
  - Neither is ruled out.

## Correction to Experiment 4
- Experiment 4 used b = 2.35 Å for the glass. Its measured mean bond is 2.735 Å.
- Against its own bond length, the raw-glass coefficient is −47, not −63. The pass/fail results are unchanged.

## Next (to be pre-registered)
- **Smaller mass η,** e.g. 0.075, on raw glass s3. At the same speed the wavelength doubles, so the error should fall about 4×. This tests the (k·b)² law directly.
- **A partly ordered network** (a glass with a controlled amount of crystal order), to see whether a falls smoothly between −50 and −5.5.

## Files
- `run_4b.log`: the run log.
- `res_smooth_tiled.json`, `res_raw_tiled.json`: all results.
- `fits_4b.json`: the fits.
- `fig4b.py`, `exp4b_results.png`: the figure.
