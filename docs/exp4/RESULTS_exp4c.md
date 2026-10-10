# Experiment 4c: Results — lighter and heavier particles on raw glass
Gameboard Universe. Run October 8, 2026. Design: `DESIGN_exp4c.md`, with a table correction logged before any moving-particle result.

## Headline
- **The pre-registered verdict is "speed-driven", not the wavelength law.**
  - Q = 0.79. The wavelength law predicted 0.25; a pure speed effect predicts 1, and anything ≥ 0.75 counts as speed-driven.
  - The heavy-particle check also fails: heavy errors were about 3× smaller than the wavelength law predicts.
- **Cause, found by a post-hoc diagnostic (not pre-registered):** the speed limit c used for the glass was wrong.
  - I had computed c with a smooth-average formula. That formula is exact for a crystal but overestimates the speed in a glass.
  - Measured directly, the raw glass's real long-wave c is 1.549 Å/tick, 3.5% below the formula's 1.606.
  - Using a c that is too high makes the relativity prediction too high, and the gap grows with speed. That gives exactly the speed-driven error seen here.
- **With the measured c, every glass clock passes.**
  - Raw, raw tiled and smoothed glass, all at 0.3c–0.8c: within 1.7% (limit 3%).
  - Diamond: within 0.8%.
- **This overturns the main glass conclusions of Experiments 4 and 4b:**
  - glass failing at 0.7c and 0.8c;
  - the a·(k·b)² "wavelength law";
  - glass being 9× worse than crystal.
  - All three came from the wrong c.

## The diagnostic: measuring c directly
- **Method:** run the node totals' own wave equation with a tiny mass (η = 0.05), using waves that fit the box exactly. c is read from the wave frequency, with no particle tracking involved.
- **Control:** diamond gives 1.0000 × the formula at every wave tested, as it should for a crystal.
- **First attempt (discarded):** it used waves that don't fit the box, and failed the diamond control (0.96–1.20). The results above come from the corrected method.

| Network | Formula c | Measured c | Ratio | Direction spread |
|---|---|---|---|---|
| Raw glass s3 | 1.606 | 1.549 | 0.965 | ≤ 0.4% |
| Raw glass s1, tiled | 1.585 | 1.529 | 0.965 | 0.3% |
| Smoothed glass s1, tiled | 1.382 | 1.325 | 0.959 | 0.2% |
| Diamond | 1.356 | 1.356 | 1.0000 | 0.04% |

## Clock errors recomputed with the measured c
- Same runs, same clock readings; only c changes.

| Network (η = 0.15) | 0.3c | 0.5c | 0.7c | 0.8c |
|---|---|---|---|---|
| Raw glass s3 | +0.1/+0.3% | −0.1/+0.1% | −0.6/−0.4% | −1.5/−1.4% |
| Raw s1, tiled | +0.1/+0.3% | −0.1/+0.2% | −0.7/−0.4% | −1.7/−1.3% |
| Smoothed, tiled | 0.0/+0.2% | −0.1/0.0% | −0.6/−0.5% | −1.6/−1.2% |
| Diamond | 0.0/+0.2% | −0.1/0.0% | −0.3/−0.2% | −0.5/−0.8% |

- **Experiment 4c's own runs, with the measured c:**
  - Light (η = 0.075): −0.2% to +1.9%, plus +2.9% at 0.3c, the point flagged in advance as least reliable.
  - Reference (η = 0.15): −0.2% to −1.4%.
  - Heavy (η = 0.30): −0.5% to 0.0%.
  - All pass the 3% limit.
- **What is left:** at 0.8c the glasses run about 1.2–1.7% slow and diamond about 0.5–0.8%.
  - That leftover may be a real effect of disorder, roughly 2× the crystal's, not 9×.
  - It is too small to pin down here; the light particles are too noisy.
- **Still failing, whatever c is used:**
  - the speed change between halves: glass packets speed up by 1–3%, and by up to 7% on smoothed glass at 0.8c;
  - some direction speed agreements.
  - These look like particle-centre tracking problems: motionless link modes, and the spread of wavelengths in a packet. They are not settled.

## Honest accounting
- **The registered verdicts stand as recorded:**
  - Experiment 4, glass: FAIL at 0.7c and 0.8c.
  - Experiment 4b: "smoothing doesn't matter".
  - Experiment 4c: "speed-driven".
- **The reanalysis is post-hoc.** It is credible for three reasons:
  - c was measured without looking at any clock data;
  - the diamond control is exact;
  - the same correction fixes all three glasses at once, with no adjustable numbers.
- **This was a design error on my part.** Stage A in Experiment 4 computed c from a formula rather than measuring it. The diamond runs couldn't catch it, because the formula is exact for a crystal.
- **Still true:** smoothing doesn't matter. With the measured c, smoothed and raw glass agree.

## What this means for the theory
- Under the GBU's own rule, raw glass gives relativistic time dilation within 1.7% up to 0.8c, in every direction tested. That is as good as needed at this size.
- The glass vs crystal gap at high speed is now small and unresolved. Glass keeps its argued advantage of having no preferred directions.
- **The glass's real speed limit isn't given by a simple average over its links.** It has to be measured. That matters for any future comparison.

## Next (to be pre-registered)
- **Find and fix the speed-tracking problem** (packets speeding up between halves). Then re-test speed stability and direction agreement, which are the remaining failed criteria.
- **Two-particle or interaction tests.**
- **Bigger glass,** to pin down the leftover 0.8c effect.

## Files
- `run_4c.py`, `run_4c.log`, `res_4c.json`: the registered runs.
- `diag_c3.py`, `diag_c3.json`, `diag_c4.py`, `diag_c4.json`: the c measurement.
- `diag_c2.py`: the discarded first attempt.
- `recalc_ctrue.json`: the recomputed errors.
- `fig4c.py`, `exp4c_results.png`: the figure.
