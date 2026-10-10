# Experiment 4: Results — the GBU's own update rule (version 2, corrected speed limit)
Gameboard Universe. Run October 8, 2026. Design: `DESIGN_exp4.md`, including Amendment 1, which was logged before any registered run.

> **What changed in version 2:**
> - Version 1 judged the glass against a speed limit c from a smooth-average formula: 1.607 Å/tick.
> - Experiment 4c measured the glass's real long-wave c directly: 1.5495 Å/tick, 3.6% lower. The method was checked on diamond, where it gives exactly the formula value, as it must for a crystal.
> - This version re-scores the same runs with the measured c.
> - The registered verdict (using the formula c) is kept below, as recorded.
> - The "(k·b)² law" and "glass 9× worse than crystal" conclusions in version 1 came from the wrong c and are withdrawn.

## What was run
- **Rule:** selection with memory. Each node holds 4 complex amplitudes, one per link (direction-valued, phase-carrying). The same mixing step runs everywhere, then every amplitude moves one link. One free number: η = 0.15 per tick (the mass).
- **Networks:** both 1,728,000 nodes, box 325.6 Å.
  - Diamond crystal.
  - Raw glass (`glass_n1728000_s3_raw`).
- **Particles:** wave packets 40 Å wide, run for 250 ticks.
  - Speeds: rest, 0.3c, 0.5c, 0.7c, 0.8c.
  - Two directions: [100] and [1,2,3].
- **Classical check:** the same rule with plain probabilities instead of amplitudes.

## Summary (corrected)
- **Clocks: PASS on both networks.** Every moving clock is within 1.5% of relativity up to 0.8c, in both directions (limit 3%).
  - Diamond: within 0.8%.
  - Glass: within 1.5%.
- **Conservation:** the total was conserved to 2×10⁻¹⁴ (limit 10⁻⁹). Nothing scattered out of the packets.
- **Classical check:** no motion and no clock. Phase is needed.
- **Speed stability and direction agreement are not settled here.**
  - On glass the speed changed 1.4–3.1% between halves at 0.3c and 0.8c.
  - **Experiment 4d re-tested this properly** (4 starting points, speed tracked over time):
    - **At 0.5c:** steady speed, and the two directions agree to 0.2%.
    - **At 0.8c:** a start-up transient. The particle starts at about 0.72c, settles near 0.76c by about tick 120, and does not slow down. It is absent on diamond and smaller for bigger particles.
    - Under the registered rules, 0.8c still fails speed stability and direction agreement. Whether the settled speed is steady is the next test.
    - See `RESULTS_exp4d.md`.

## Clock error vs relativity, with the measured c

| Speed | Diamond [100] / [1,2,3] | Raw glass [100] / [1,2,3] |
|---|---|---|
| 0.3c | +0.02% / +0.18% | +0.13% / +0.26% |
| 0.5c | −0.07% / +0.04% | −0.05% / +0.13% |
| 0.7c | −0.25% / −0.22% | −0.56% / −0.44% |
| 0.8c | −0.52% / −0.80% | −1.48% / −1.43% |

- **Measured speed limits:**
  - Diamond: 1.3562 Å/tick, the same as the formula.
  - Glass: 1.5495 Å/tick, against the formula's 1.607.
  - Both are the same in every direction to within 0.4%.
- **The leftover at 0.8c:** glass runs about 1.5% slow and diamond about 0.5–0.8%. Glass is roughly 2× the crystal, not 9×. Whether this leftover is a real effect of disorder is still open.

## Registered verdict, as recorded (formula c)
- **Diamond: PASS on clocks.**
  - Two speed readings at 0.3c missed their 1% limits by a hair: 1.02% between halves, and 1.07% between directions.
- **Glass: PASS at 0.3c and 0.5c, FAIL at 0.7c and 0.8c.** Clocks read 3.0–3.1% and 5.3–5.5% slow.
  - Speed crept up 1.4–3.1% at 0.3c and 0.8c (fail).
  - Directions agreed within 1% (pass).
- **Why the glass clock fails were wrong:** a speed limit set too high makes the relativity prediction too high, and the gap grows with speed. The measured c removes the failures without any adjustable numbers. The same correction fixes three different glasses (Experiment 4c).

## What this shows and doesn't show
- **Shows:** a GBU-native rule produces relativistic time dilation on both a crystal and a raw random glass, to within 1.5% up to 0.8c. No relativity was put in by hand: direction-valued amplitudes on 4 links, one fixed mixing step, and one-link-per-tick motion are enough.
- **Lesson:** on a disordered network the speed limit isn't a simple average over the links. It has to be measured.
- **Doesn't show:**
  - These are free particles: no binding, no interactions, no gravity.
  - Length contraction can't be tested on free packets.
  - One glass sample.
- **Caveats:**
  - **Amendment 1:** the way the clock is read was changed after the 216k-node debug runs and before the registered runs. The reason was a known spreading effect that fails a perfect crystal too.
  - **The c correction** came after seeing results (Experiment 4c). It is credible because c was measured without looking at any clock data, the diamond control is exact, and no adjustable numbers are involved.

## Withdrawn from version 1
- **"Glass error ≈ −63·(k·b)²" and "glass is 9× worse than crystal":** both came from the wrong c.
- **"b = 2.35 Å" for the glass:** its mean bond is 2.735 Å.

## Files
- `exp4.py`: the rule and the measurements.
- `run_main.py`: the registered run.
- `res_diamond.json`, `res_glass.json`: all results.
- `check_recurrence.py`: proof that the node totals obey the relativistic wave recurrence.
- `exp4_results.png`: version 1 figure (formula c).
- `exp4_results_v2.png`: version 2 figure (measured c).
- `diag_c3.py`, `diag_c3.json`: the c measurement.
