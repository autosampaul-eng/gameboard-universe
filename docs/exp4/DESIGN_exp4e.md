# Experiment 4e: do fast particles on glass hold a steady speed once settled? (pre-registered)
Gameboard Universe. Written October 8, 2026, before any 4e runs. Rule, clock and tracking are as in Experiment 4d, with the measured c.

## Question
- In Experiment 4d, glass particles at 0.8c started at about 0.72c, climbed between ticks 60 and 120, and levelled off near 0.76c.
- The runs ended at tick 200, too soon to judge what happens after settling.
- Once settled, is the speed steady, with no drag, and the same in both directions?

## Runs
- **Glass:** raw glass s3, 1.73M nodes, η = 0.15, v/c = 0.8 (launch setting, measured c = 1.5495 Å/tick), **400 ticks**.
  - **Main set:** σ = 40 Å, directions [100] and [1,2,3], the same 4 starting points as in 4d. 8 runs.
  - **Size set:** σ = 50 Å, [100], the same 4 starting points. 4 runs. These show whether bigger particles settle faster.
- **Control:** diamond, σ = 40 Å, 0.8c, [100] and [1,2,3], 1 starting point each. 2 runs.
- **Wrapping:** over 400 ticks a particle travels about 500 Å, more than the 326 Å box, so it wraps around once.
  - Tracking uses nearest-image coordinates around the current centre, which handles the wrap.
  - The diamond control checks it: a crystal particle should keep a flat speed through the wrap.

## Measures
- Speed in windows of 40 ticks, stepping 20 ticks.
- **Settled period:** ticks 150 to 400.
- **Drag D:** slope of window speed vs time over the settled period, in % per 100 ticks.
- **Settled speed:** a straight-line fit of position vs time over the settled period.
- **Halves check:** settled speed over ticks 150–275 vs 275–400.
- **Clock:** the carrier-frequency clock, fitted over the settled period, using the settled speed.
- **Settling time:** the centre of the first window whose speed is within 0.5% of the settled speed, reported but not scored.
- Results are given as mean ± SE over the starting points.

## Pre-registered pass criteria
1. **Wrap-safety control (scored first):** on diamond, |D| < 0.1% per 100 ticks and wobble W < 0.1%. If this fails, the glass results are reported as unreliable.
2. **No drag after settling:** for each glass direction (σ = 40 Å), |mean D| < 1% per 100 ticks, **and** within 2 SE of zero.
3. **Steady:** for each glass direction, the halves check differs by < 1%.
4. **Directions agree:** the settled speeds for [100] and [1,2,3] differ by < 1%, **and** by < 2 SE of the difference.
5. **Clocks:** every glass run is within 3% of relativity over the settled period.
6. **Total:** conserved to 10⁻⁹.

If criteria 1–4 pass, the 0.8c failures in Experiment 4d are put down to a start-up transient. That would mean fast particles on glass, once settled, move steadily and the same in every direction. If any of 2–4 fails, the result is reported as a real effect of the glass.

## Limits stated in advance
- One glass sample. The starting points share the same glass, and the particles cross the same box more than once.
- Free particles spread over time. By tick 400 the packet's tails may reach halfway round the box; the diamond control is there to catch any effect of that.

Any outcome is reported as found.
