# Experiment 4e: Results — fast particles on glass after settling

> **Note added Oct 9, 2026 (after Experiment 4g):** the speed measurements here use the whole-box centre. Experiment 4g found that on glass this centre is pulled around by a scattered cloud whose effect depends on box size. The start-up transient, overshoot and direction gap reported here are most likely that tracking artifact. The clock results are not affected. See `RESULTS_exp4g.md`.
Gameboard Universe. Run October 8, 2026. Design: `DESIGN_exp4e.md`, written before any 4e run. Raw glass s3 and diamond, 1.73M nodes, 0.8c launch, 400 ticks, with the measured c.

## Verdict
- **Wrap-safety control: PASS.** Diamond particles keep a flat speed through the wrap: trend ≤ 0.06% per 100 ticks, wobble ≤ 0.02%. So the tracking is reliable over 400 ticks.
- **Direction [100]: PASS.**
  - Trend −0.17 ± 0.48% per 100 ticks (no drag).
  - Halves differ by −0.36%.
- **Direction [1,2,3]: FAIL by a hair on two rules.**
  - Trend −0.96 ± 0.22% per 100 ticks: inside the 1% limit, but 4.5 standard errors from zero.
  - Halves differ by −1.02% (limit 1%).
- **Directions agree: FAIL.** The settled speeds differ by 0.88%, which is inside 1%, but 3.8 standard errors from zero.
- **Clocks: PASS.** All 12 glass runs are within 2.0% of relativity, using each run's settled speed.
- **Total:** conserved to 4×10⁻¹⁴ (pass).

## What the speed actually does

| Phase | Ticks | Speed | Notes |
|---|---|---|---|
| Start | 0–80 | dips to 0.72c | the clean wave I build doesn't suit the glass |
| Rise | 80–120 | climbs | |
| Overshoot | about 120 on [100], about 180 on [1,2,3] | about 0.765c | |
| Settled | from about 220–250 | around 0.745c | wobbles ±1% as the particle crosses patches of glass |

- **The overshoot comes later on [1,2,3].** So its scored period, which starts at tick 150, includes the fall from the overshoot. That explains why [1,2,3] shows a slight slowdown, and probably the direction gap too.
- **Post hoc only,** using just the period after tick 250:
  - settled speeds 0.7420 ([100]) and 0.7471 ([1,2,3]);
  - the directions differ by 0.69 ± 0.43%, which is inside 1% and 1.6 standard errors.
  - This is a lead, not a pass.
- **Bigger particles (50 Å):** the same rise, overshoot and settle, with no drag (−0.16 ± 0.61%) and a settled speed of 0.7468c.
- **Experiment 4d's "late speed" of about 0.756c was the overshoot,** not the settled value.

## Numbers (mean ± SE over 4 starting points)

| Case | Settled v/c | Trend D (% per 100 ticks) | Halves | Clock error |
|---|---|---|---|---|
| Glass [100], 40 Å | 0.7426 ± 0.0010 | −0.17 ± 0.48 | −0.36% | −1.85 to −1.49% |
| Glass [1,2,3], 40 Å | 0.7492 ± 0.0013 | −0.96 ± 0.22 | −1.02% | −1.72 to −1.50% |
| Glass [100], 50 Å | 0.7468 ± 0.0006 | −0.16 ± 0.61 | −0.55% | −1.95 to −1.70% |
| Diamond [100] / [1,2,3] | 0.7910 / 0.7871 | 0.00 / +0.06 | 0.00 / +0.07% | −0.53 / −0.80% |

## What this settles and what it doesn't
- **Settled:**
  - Glass doesn't drag fast particles. Along [100] the settled speed holds steady over 250 ticks, and the 50 Å particles behave the same.
  - Throughout, the clocks follow relativity at the particle's actual speed.
  - Diamond shows none of these effects.
- **Not settled:**
  - Whether the [1,2,3] slowdown and the 0.7–0.9% direction gap are real.
  - Both may just be the overshoot ending later on [1,2,3], since the scored period started too early for it.
  - A real direction difference at this short wavelength (about 48 Å, or 18 bonds) in a single glass sample is also possible. Stage A found the long-wave speed limit the same in every direction to within 0.4%.
- **Another open point:** the settled speed, about 0.745c, is about 7% below what a relativistic particle with this wave would have (0.8c). The clock still matches relativity at the speed actually travelled. So on glass, the relationship between a particle's wave and its speed is shifted, but time dilation still follows the speed actually travelled.

## Next (to be pre-registered)
- **Start the scored period after tick 250,** and run to tick 600 so it covers 350 ticks. On the 326 Å box that means wrapping about twice; the diamond control covers this.
- **Or skip the transient altogether** with pre-settled particles: run a packet for 250 ticks, then use that settled state as the starting particle.
- **A second glass sample,** to tell a real direction effect from a quirk of this sample. A bigger machine would help here.

## Files
- `run_4e.py`, `run_4e.log`, `res_4e.json`: the registered runs.
- `fig4e.py`, `exp4e_results.png`: the figure.
