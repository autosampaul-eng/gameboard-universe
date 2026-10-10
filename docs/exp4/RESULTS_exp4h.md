# Experiment 4h: Results — glass scatters fast particles continuously
Gameboard Universe. Run October 10, 2026, as 25 parallel jobs on GitHub (repository `autosampaul-eng/gameboard-universe`, workflow `exp4h`). Design: `DESIGN_exp4h.md`, written before any 4h run.
- **Method check beforehand:** on the 216k test networks. Diamond showed exactly zero scattering. The test glass showed steady scattering, with a smooth particle-only speed. Those were seen before the registered runs and changed nothing in the design.

## Registered verdict

| # | Criterion | Result |
|---|---|---|
| 0a | Diamond: \|Δs\| < 0.5 points at all ticks | **PASS** (0.00 at every tick, all 4 runs) |
| 0b | Diamond: \|Γ\| < 0.2 per 100 ticks | **PASS** (0.000) |
| 0c | Diamond particle-only speed within 0.1% of 4g | **FAIL by a hair:** −0.115% ([100]), −0.111% ([1,2,3]) |

- **Control 0 failed** (on 0c), so, as pre-registered, **the method is reported as not working, and Q1 and Q2 are not judged.**
- The scattering half of the control (0a, 0b) passed perfectly. The failure is a fixed 0.11% offset in the particle-only speed, the same in both directions with no trend. It most likely comes from measuring the grid density of the node totals rather than the link density used in 4g.
- **What Q1 and Q2 would have given,** shown for the record but not counted as verdicts:

| Criterion | [100] | [1,2,3] | Would be |
|---|---|---|---|
| Q1: late scattering rate Γ (points per 100 ticks) | 5.44 ± 0.07 (81 SE) | 5.41 ± 0.09 (64 SE) | **"Continuous"** in both directions |
| Q2: drag D (% per 100 ticks) | −0.19 ± 0.09 (2.0 SE) | −0.17 ± 0.06 (2.6 SE) | [100] pass, [1,2,3] fail (tiny but non-zero) |
| Q2: halves | −0.31% | −0.26% | pass |
| Q2: directions | 0.7442 vs 0.7446 | gap −0.05 ± 0.10% (0.5 SE) | pass |
| Clocks (within 3%) | −2.07 to −1.75% | −1.85 to −1.53% | pass |
| Total conserved | 3×10⁻¹² | | pass |

## What the runs show

### Glass scatters fast particles steadily, from the first tick
- **Rate:** at 0.8c the scattered share of the wave grows almost linearly, by about 5.4 percentage points per 100 ticks.
  - It reaches about 17% by tick 250 and 35–38% by tick 600.
  - It is the same in both directions and from all 6 starting points (spread about ±0.2 points).
- **Not a set-up loss:** the scattering doesn't level off. The "start-up loss" at tick 250 is just the same steady rate running from tick 0.
- **Not a box effect:** the 326 Å box gives Γ = 5.53 ± 0.10 vs 5.52 ± 0.11 for the matching big-box runs, and the same start-up loss (17.3 vs 17.2 points).
- **Mean free path ≈ 2,100 Å.** That is about 780 bonds, or about 44 of the particle's own wavelengths (48.7 Å).
- **Diamond scatters nothing,** at either speed.

### Speed dependence (reported, not scored)
- **At 0.5c** (wavelength about 112 Å, or 41 bonds), scattering is 16–28× slower: Γ = 0.34 ± 0.02 ([100]) and 0.19 ± 0.06 ([1,2,3]).
  - Mean free path about 22,000–59,000 Å.
- **Per unit distance,** the rate falls about 13× when the wavenumber halves (k₀ ratio 2.31). That is roughly a k³ dependence.
  - Post hoc, from two speeds only. Not a measured law.

### With the scattered waves removed, the particle moves steadily
- The particle-only speed is 0.744–0.746c from every starting point and in both directions. The window speeds change by less than 1% over 600 ticks.
- The big dip and overshoot seen in 4d–4g appear only in the whole-box centre, which mixes particle and cloud. That confirms the 4g reading: they were tracking artifacts.
- **The direction gap is gone** for the particle itself: −0.05 ± 0.10%.
- **A slight slowdown remains,** about −0.18% per 100 ticks, statistically clear on [1,2,3].
  - It could be small-angle scattering staying inside the 30° cone, or a real effect. These runs can't tell.
- **The clocks follow relativity** at the particle's own speed, within 2.1%.

### Measurement limit found at 0.5c (not scored)
- At 0.5c the particle's own spread of directions is wider than the cone/shell allows. The baseline s(0) is 6.5–7% on glass and 3% on diamond.
- So the filter trims part of the particle itself. That biases the particle-only speed along [1,2,3]: diamond at 0.5c shows +0.54% against the whole-box speed.
- The 0.5c "direction gap" (1.25%) is therefore not trustworthy. The 0.5c scattering rates are less affected, since Δs removes the fixed baseline.

## What it means
- **The earlier puzzles are explained.** Moving through this glass, a fast particle steadily sheds part of itself into waves heading off in other directions. The part that remains moves steadily, at the same speed in both directions, and its clock follows relativity.
- **This is a real physical property of the glass network, not a tracking problem,** and an important one:
  - It is the same in every direction and box, and diamond shows none.
  - A particle that loses coherence while moving through empty space, at a rate set by its wavelength, would make the network's rest frame detectable.
- **What decides whether the model survives** is how the mean free path grows as the wavelength gets longer compared to the bond length.
  - Real particles are 10¹⁶ or more node spacings across, so even a steep power law could put the effect far beyond any measurement.
  - With two speeds, the exponent (about 3) isn't pinned down. Extrapolating 15 orders of magnitude from it would not be defensible.

## Next (to be pre-registered)
- **Experiment 4i, scattering vs wavelength:**
  - Γ at 4–5 wavelengths (about 18 to 80 bonds), using a wider σ so the slow particles fit inside the cone and shell.
  - Fit the exponent, with errors, on this glass.
  - Then repeat on an independent glass sample, for which the GitHub machines are well suited.
- **Fix the 0.11% tracker offset:** check particle-only tracking with the link density on diamond, so that control 0c can pass in future.

## Files
- `run_4h.py`: the runner (also on GitHub, `code/run_4h.py`).
- `res_4h/`: 25 job results and logs, also under `results/exp4h/` on GitHub.
- `test_4h_*.json`: the method checks on the 216k networks.
- `fig4h.py`, `exp4h_results.png`: the figure.
