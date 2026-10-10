# Experiment 4f: Results — 600-tick runs at 0.8c, scored after tick 250

> **Note added Oct 9, 2026 (after Experiment 4g):** the speed measurements here use the whole-box centre. Experiment 4g found that on glass this centre is pulled around by a scattered cloud whose effect depends on box size. The start-up transient, overshoot and direction gap reported here are most likely that tracking artifact. The clock results are not affected. See `RESULTS_exp4g.md`.
Gameboard Universe. Run October 9, 2026. Design: `DESIGN_exp4f.md`, written before any 4f run. Raw glass s3 and diamond, 1.73M nodes, 326 Å box, 0.8c launch, σ = 40 Å, 6 starting points per glass direction.

## Verdict: the control failed, so the glass results are unreliable (as pre-registered)
- **Control: FAIL.**
  - Diamond [100] is flat: trend 0.00% per 100 ticks, wobble 0.00%.
  - Diamond [1,2,3] drifts: trend +0.28% per 100 ticks and wobble 0.11%, against limits of 0.1%. The drift starts after about tick 400.
- **Spreading check: FAIL.** 9 of 12 glass runs give settled speeds that differ by more than 0.3% between the two tracking methods; the biggest difference is 1.5%.
  - By tick 600 the particles have spread to about 130–140 Å, close to half the box width (163 Å).
  - Their tails wrap around the box and pull on the measured centre.
- **Glass criteria, scored as registered but unreliable:**

| Criterion | [100] | [1,2,3] | Result |
|---|---|---|---|
| No drag (< 1% per 100 ticks and < 2 SE) | −0.51 ± 0.19 (2.7 SE) | +0.44 ± 0.39 (1.1 SE) | [100] fail, [1,2,3] pass |
| Steady (halves < 1%) | −0.81 ± 0.50% | +1.07 ± 0.46% | [100] pass, [1,2,3] fail |
| Directions agree (< 1% and < 2 SE) | 0.7422 vs 0.7509 | difference −1.16 ± 0.18% (6.4 SE) | fail |
| Clocks (within 3%) | −2.11 to −1.63% | −1.66 to −1.21% | pass |
| Total conserved | 5×10⁻¹⁴ | | pass |

## What holds up despite the tracking problem
- **The direction gap is consistent, with [1,2,3] faster in every measurement:**

| Measurement | Gap |
|---|---|
| 4e, ticks 150–400 | 0.88% |
| 4e, ticks 250–400, post hoc | 0.69% |
| 4f, whole-box tracking | 1.16% |
| 4f, limited-window tracking | 1.00% (4.0 SE) |

  - Diamond has its own small gap at this wavelength: 0.7910 vs 0.7890, which is 0.25%.
  - The glass gap is about 4× bigger than that.
- **The drag and halves trends have no consistent sign.** They flip between directions and starting points. That points to wobble and wrap-around error, not a steady slowdown.
- **The settled speed matches 4e:** about 0.742c ([100]) and 0.751c ([1,2,3]), compared with 0.743c and 0.749c in 4e.
- **Clocks follow relativity at the speed actually travelled,** within 2.1%, in every run.

## What it means
- **The method has reached its limit on this box.** A 326 Å box can't hold a 0.8c particle for 600 ticks, because it spreads until it wraps onto itself. The diamond control caught this, as it was designed to.
- **The question is still open:** is there a direction difference of about 1% for fast, short-wavelength (18-bond) particles in this glass? Every measurement so far points to yes, but none passes a clean test.
  - If real, it may be a quirk of this one glass sample or a general property of glass at short wavelengths.
  - Stage A found the long-wave speed limit the same in every direction to within 0.4%.

## Options for a clean test (none pre-registered yet)
1. **A bigger box, made locally:** tile the s3 glass 2×2×2 to 13.8M nodes and a 652 Å box. The particle would then never reach halfway round the box in 600 ticks.
   - The glass would repeat every 326 Å, which is the same glass these particles already cross.
   - It would be tight on this machine's 8 GB of memory, at about 30–40 minutes per run.
2. **Pre-settled particles:** start from a settled state, so scoring can begin at tick 0 and stop well before the tails wrap (about tick 350).
3. **A second, independent glass** of the same size, to tell a quirk of this sample from a general glass property. This is the clearest case for a bigger machine.

## Files
- `run_4f.py`, `run_4f.log`, `res_4f.json`: the registered runs.
- `fig4f.py`, `exp4f_results.png`: the figure.
