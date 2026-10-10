# Experiment 4d: Results — speed stability and direction agreement

> **Note added Oct 9, 2026 (after Experiment 4g):** the speed measurements here use the whole-box centre. Experiment 4g found that on glass this centre is pulled around by a scattered cloud whose effect depends on box size. The start-up transient, overshoot and direction gap reported here are most likely that tracking artifact. The clock results are not affected. See `RESULTS_exp4g.md`.
Gameboard Universe. Run October 8, 2026. Design: `DESIGN_exp4d.md`, written before any registered run. Raw glass s3 and diamond, 1.73M nodes each, judged with the measured speed limits.

## Verdict
- **At 0.5c: PASS, with one strict-rule miss.**
  - No drag on [100]: −0.01 ± 0.09% per 100 ticks.
  - On [1,2,3] the trend is +0.41 ± 0.12% per 100 ticks. That's well inside the 1% limit, but 3.4 standard errors from zero, so under the rules it counts as a fail.
  - Directions agree: −0.17 ± 0.28% (pass).
- **At 0.8c: FAIL.** The glass particles don't slow down; they **speed up**, by +3.5% per 100 ticks on [100] and +4.9% on [1,2,3].
  - The two directions' average speeds differ by 1.48% (fail).
- **The wobble really comes from the glass: PASS.**
  - Diamond's wobble is 0.005%.
  - The glass's wobble shrinks as the particle grows: W(50 Å) ÷ W(30 Å) = 0.60 (limit ≤ 0.8).
- **Clocks: PASS.** All 24 glass runs are within 2.1% of relativity.
- **Total:** conserved to 2×10⁻¹⁴ (pass).

## What the 0.8c failure actually is
- **It isn't drag or random wobble. It's a start-up transient.**
  - Every 0.8c glass run starts at about 0.72c.
  - It climbs between ticks 60 and 120.
  - It then levels off near 0.75–0.77c.
  - The same curve appears from all 4 starting points, so it is a property of the particle, not of where it happens to be.
- **It shrinks for bigger particles:** the trend is +4.5% at 30 Å, +3.5% at 40 Å and +2.8% at 50 Å.
- **Diamond shows nothing:** 0.000% trend, with speeds flat to 4 decimal places.
- **Likely cause:** I build the particle as a clean wave. A clean wave is a natural travelling state of the crystal, but not of the glass. On glass, the particle has to settle into the glass's own travelling wave, which carries it slightly faster. At 0.5c the wave is long and the glass looks smooth to it, so there is almost nothing to settle.
- **The direction failure at 0.8c comes from the same transient.** It settles more slowly along [1,2,3].
  - Post hoc, using only the last three windows (ticks 140–180), the directions agree: −0.27 ± 0.56%.
  - That comparison wasn't pre-registered. It is only a lead.

## Numbers (glass, mean ± standard error over 4 starting points)

| Case | Speed v/c | Trend D (% per 100 ticks) | Wobble W | Clock error |
|---|---|---|---|---|
| 0.5c [100], 40 Å | 0.4795 ± 0.0010 | −0.01 ± 0.09 | 0.07% | −0.18 to −0.04% |
| 0.5c [1,2,3], 40 Å | 0.4804 ± 0.0009 | +0.41 ± 0.12 | 0.12% | +0.05 to +0.15% |
| 0.8c [100], 30 Å | 0.7417 | +4.49 ± 0.65 | 2.15% | −2.0 to −1.0% |
| 0.8c [100], 40 Å | 0.7485 ± 0.0013 | +3.48 ± 0.47 | 1.69% | −2.1 to −1.4% |
| 0.8c [100], 50 Å | 0.7508 | +2.84 ± 0.32 | 1.29% | −1.9 to −1.5% |
| 0.8c [1,2,3], 40 Å | 0.7376 ± 0.0012 | +4.91 ± 0.37 | 0.80% | −1.5 to −1.1% |
| Diamond 0.8c, both directions | 0.7911 / 0.7868 | 0.00 | 0.005% | −0.52 / −0.80% |

- At 0.8c the wobble W partly comes from fitting a straight line to an S-shaped rise, so it overstates the patch-to-patch wobble.

## What this settles and what it doesn't
- **Settled:**
  - The earlier "speed creep" is not a flaw in the tracking: three tracking methods agree.
  - It's also not drag: nothing slows down at any speed.
  - At 0.5c, glass particles keep a steady speed, the same in both directions, to within a few tenths of a percent.
  - Clocks follow relativity everywhere.
- **Not settled:** whether a fast particle on glass, once it has settled, keeps a constant speed that's the same in all directions. These runs were too short to see much after the transient.

## Next (to be pre-registered)
- **Experiment 4e: longer runs at 0.8c,** about 400 ticks, judging speed only after the transient (from tick 150 on), from the same 4 starting points and both directions.
  - Pass rules: the same as here.
  - Long paths cross the 326 Å box more than once. A wrap-safe tracking check comes first.
- **Alternatively, build "pre-settled" particles,** by letting a packet settle and then re-centring it, so the transient never happens.

## Files
- `diag_track.py`: the debug runs on the 216k-node test networks.
- `run_4d.py`, `run_4d.log`, `res_4d.json`: the registered runs.
- `fig4d.py`, `exp4d_results.png`: the figure.
