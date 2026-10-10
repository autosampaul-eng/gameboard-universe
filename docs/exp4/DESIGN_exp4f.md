# Experiment 4f: settled speed of fast glass particles over a long path (pre-registered)
Gameboard Universe. Written October 9, 2026, before any 4f runs. Rule, clock, measured c and tracking are as in Experiment 4e.

## Question
- In Experiment 4e, scoring started at tick 150. On [1,2,3] that was before its speed overshoot had ended.
- With scoring after tick 250 and a long scored period, is a settled 0.8c glass particle:
  - free of drag;
  - steady between halves;
  - the same speed in both directions?

## Runs
- **Glass:** raw glass s3, η = 0.15, launch v/c = 0.8 (measured c = 1.5495 Å/tick), σ = 40 Å, **600 ticks**.
  - Directions [100] and [1,2,3].
  - **6 starting points each:** the 4 from Experiments 4d and 4e, plus (0.75, 0.75, 0.75) and (0.5, 0.25, 0.5), as fractions of the box. 12 runs.
- **Control:** diamond, the same settings, [100] and [1,2,3], 1 starting point each. 2 runs.

## Measures
- **Scored period:** ticks 250 to 600 (350 ticks).
- **Drag D:** slope of window speed (40-tick windows, stepping 20) vs time over the scored period, in % per 100 ticks.
- **Settled speed:** a straight-line fit of position vs time over ticks 250–600.
- **Halves check:** ticks 250–425 vs 425–600.
- **Clock:** the carrier-frequency clock over ticks 250–600, using the settled speed.
- **Tracking:**
  - **Registered:** the whole box, nearest-image around the current centre (as in 4d and 4e).
  - **Spreading check:** by tick 600 the packets may reach halfway round the box. Settled speeds are also computed with the centre limited to within 0.45 × box of the current centre. If the two methods differ by more than 0.3% for any run, that run is flagged, and the glass results are reported as possibly affected by wrap-around.

## Pre-registered pass criteria
1. **Control (scored first):** diamond |D| < 0.1% per 100 ticks and wobble < 0.1%. Diamond shows no settling, so its scored period is the same.
2. **No drag:** for each glass direction, |mean D| < 1% per 100 ticks, **and** within 2 SE of zero.
3. **Steady:** for each glass direction, the mean halves difference < 1%.
4. **Directions agree:** the settled speeds for [100] and [1,2,3] differ by < 1%, **and** by < 2 SE of the difference.
5. **Clocks:** every glass run is within 3% of relativity over the scored period.
6. **Total:** conserved to 10⁻⁹.

- **If 1–4 pass:** fast particles on this glass, once settled, move steadily and the same in both tested directions. The 4d and 4e failures are put down to the start-up transient.
- **If 4 fails while 2–3 pass:** the result is reported as a real direction difference in this glass sample at this wavelength (about 48 Å, or 18 bonds). Telling that apart from a quirk of this one sample needs a second glass.

## Limits stated in advance
- One glass sample. Over 600 ticks the particles cross the 326 Å box about twice, so they pass through the same glass more than once.
- Packets spread out, so their tails may wrap around the box. The spreading check is there to catch this.

Any outcome is reported as found.
