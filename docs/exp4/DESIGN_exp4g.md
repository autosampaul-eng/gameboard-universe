# Experiment 4g: the direction gap in a box big enough to avoid wrap-around (pre-registered)
Gameboard Universe. Written October 9, 2026, before any 4g runs. Rule, clock, measured c and measures are as in Experiment 4f. Only the box changes.

## Question
- In Experiment 4f the packets spread to about 135 Å, close to half the 326 Å box. Their tails wrapped around and spoiled the tracking, so the diamond control failed.
- Across Experiments 4e and 4f, fast glass particles along [1,2,3] were about 1% faster than along [100], but no clean test has shown it.
- **The question:** with wrap-around removed, do settled 0.8c glass particles have no drag, a steady speed, and the same speed in both directions?

## Networks
- **Glass:** raw glass s3 (1.73M nodes, 326 Å) tiled 2×2×2, giving **13.8M nodes in a 652 Å box**.
  - Tiling copies the same glass, so the structure repeats every 326 Å. The particles meet exactly the glass they crossed in 4f; only the self-overlap is gone.
  - Experiment 4b's tiling check confirmed that tiling preserves the network.
- **Control:** diamond (1.73M nodes) tiled 2×2×2, giving 13.8M nodes. This is the same crystal as before, in the same box size as the glass.
- **Engineering shortcuts (they don't change the physics):**
  - The tiled network's slot-pairing table and bond vectors are copies of the original's, so they are reused rather than recomputed, to fit in 8 GB of memory.
  - Tracking is moved into compiled code. It is checked against the 4f tracking code on the 216k tiled test glass before any registered run.

## Runs
- η = 0.15, launch v/c = 0.8 (measured c: glass 1.5495, diamond 1.3562 Å/tick), σ = 40 Å, **600 ticks**.
- **Glass:** directions [100] and [1,2,3], with the **same 6 starting positions as 4f**, at the same distances in Å. 12 runs.
- **Diamond:** [100] and [1,2,3], 1 starting point each. 2 runs.
- **Estimated time:** about 30–40 minutes per run, about 8 hours in total.

## Measures (identical to 4f)
- **Scored period:** ticks 250–600.
- **Drag D:** trend of 40-tick window speeds (stepping 20), in % per 100 ticks.
- **Settled speed:** a position-vs-time fit over ticks 250–600.
- **Halves:** ticks 250–425 vs 425–600.
- **Clock:** carrier frequency − k·v, over ticks 250–600.
- **Tracking:** whole box, nearest-image (registered), plus the limited-window check (within 0.45 × box of the centre).

## Pre-registered pass criteria
0. **Spreading check:** every run's two tracking methods agree to within 0.3%. If any run fails, it is reported and the glass results are marked unreliable.
1. **Control:** diamond |D| < 0.1% per 100 ticks and wobble < 0.1%, in both directions.
2. **No drag:** for each glass direction, |mean D| < 1% per 100 ticks and within 2 SE of zero.
3. **Steady:** for each glass direction, the mean halves difference < 1%.
4. **Directions agree:** the settled speeds for [100] and [1,2,3] differ by < 1% and by < 2 SE.
5. **Clocks:** every glass run within 3% of relativity.
6. **Total:** conserved to 10⁻⁹.

## Pre-registered reading
- **If 0–1 pass and 4 fails:** the direction gap is real for this glass at this wavelength (about 48 Å, or 18 bonds), not a tracking artifact. Whether it is general to glass or a quirk of this sample is left to an independent glass.
- **If 0–4 all pass:** fast particles on this glass are steady and the same in both tested directions, and the earlier gap was a box effect.
- **Comparison with 4f:** settled speeds from the same starting points are compared directly. A change of more than 0.3% would show how much wrap-around distorted 4f.

## Limits stated in advance
- One glass sample, repeated by tiling. Particles cross the same 326 Å of glass about twice.
- The tiled box can't show differences between glass samples. Only an independent glass can.

Any outcome is reported as found.
