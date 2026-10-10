# Experiment 4d: speed stability and direction agreement, done properly (pre-registered)
Gameboard Universe. Written October 8, 2026, before any 4d runs. Rule and clock are as in Experiment 4 (Amendment 1). The speed limit is the **measured** c from Experiment 4c.

## What the debug runs showed (216k-node test glass, not a registered network)
- In Experiments 4, 4b and 4c, the glass particles appeared to speed up between the two halves of a run. Three tracking methods gave the same result:
  - the original window;
  - the whole box, with no window;
  - node totals instead of link amplitudes.
  - So it isn't a tracking artifact.
- Only about 1.4% of the particle stays behind at the start, which is too little to matter.
- Over a longer run, the speed on glass **rises and falls**, from 0.684c to 0.729c and back to 0.682c in windows of 30 ticks. It doesn't drift one way. The same run on diamond is flat to 0.01%.
- **Interpretation:** a 40 Å particle samples only a few patches of glass at a time, and each patch has its own local wave speed. The particle's speed wobbles as it passes through them.
  - Splitting a run into halves can't tell this wobble from drag.
  - A single path can't tell it from a real difference between directions.

## The fixed method (registered)
- **Tracking:** centre of |amplitude|² over the whole box, using nearest-image coordinates around the current centre, with no window.
- **Speed vs time:** speed fitted in windows of 40 ticks, stepping 20 ticks, from tick 20 to tick 200.
- **Drag:** D = slope of a straight-line fit to window speed vs time, expressed as % change per 100 ticks.
- **Wobble:** W = the rms spread of the window speeds around that fit, in %.
- **Averaging over the glass:** each case is run from **4 starting points** spread through the box. Results are given as mean ± standard error (SE) over the 4.
- **Particle speeds** are set and judged with the measured c: s3 1.5495 Å/tick, diamond 1.3562 Å/tick.

## Runs
- **Network:** raw glass s3 (1.73M nodes). Diamond is the control.
- **Particles:** η = 0.15, 200 ticks.
- **Main set (glass):**
  - σ = 40 Å;
  - v/c = 0.5 and 0.8;
  - directions [100] and [1,2,3];
  - 4 starting points each.
  - 16 runs.
- **Size set (glass):** v/c = 0.8, [100], σ = 30 Å and 50 Å, 4 starting points each. 8 runs.
- **Control (diamond):** σ = 40 Å, v/c = 0.8, both directions, 2 starting points each. 4 runs.

## Pre-registered pass criteria
1. **No drag (glass, main set):** for every speed and direction, |mean D| < 1% per 100 ticks, **and** within 2 SE of zero.
2. **Directions agree (glass, main set):** at each speed, the mean speeds for [100] and [1,2,3] differ by < 1%, **and** by < 2 SE of the difference.
3. **The wobble comes from the glass, not the method:**
   - (a) diamond W < 0.3%;
   - (b) glass W falls as the particle grows, with W(50 Å) ÷ W(30 Å) ≤ 0.8. A random-patch picture predicts about (30/50)^1.5 ≈ 0.46.
4. **Clocks (glass, all runs):** within 3% of relativity, using the measured c and each run's mean speed.
5. **Total:** conserved to 10⁻⁹.

If 1 and 2 pass, the speed-stability and direction-agreement failures from Experiments 4 to 4c are put down to path wobble, not physics. If either fails, the result is reported as a real effect.

## Limits stated in advance
- One glass sample. The 4 starting points share the same glass.
- The paths are 150–250 Å long, so each particle crosses only a few patches. The standard errors will be fairly large.
- The packet-size range is limited by the box: at σ = 50 Å the particle fills a good fraction of it.

Any outcome is reported as found.
