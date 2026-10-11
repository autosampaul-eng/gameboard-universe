# Experiment 4i-R: Results — independent-glass replication
Gameboard Universe. Run October 10–11, 2026, on GitHub: the s4 build (workflow `build_glass4`), then 31 jobs (`exp4i_rep`, 34 minutes). Design: `DESIGN_exp4i_rep.md`, including Amendment 1 (the repair-pass cap), both written before any s4 particle run.

## The new sample
- **Glass s4:** same recipe as s3, seed 4. It needed 101 repair passes (s3: 46), under the cap raised by Amendment 1.
- **Both are fully relaxed:** max residual force 8.0×10⁻⁴ (s4) vs 9.6×10⁻⁴ (s3), below the 10⁻³ stopping rule.

| | s3 | s4 |
|---|---|---|
| Energy per atom | 6.79 | 7.59 (+12%) |
| Bond-length rms | 0.275 | 0.292 |
| Bond-angle rms | 30.4° | 30.8° |
| Measured c ([100], longest wave) | 1.5495 | 1.5740 (+1.6%) |
| c, [100] vs [1,2,3] | ≤ 0.4% | 0.31% |

S(111) is quoted by the builder at a single k-vector, so it is noisy (0.075 vs 0.40). Both are far below a crystal's value (about N); not used.

## Registered verdict

| # | Criterion | Result |
|---|---|---|
| 0a | Pipeline check reproduces 4i | **PASS**: identical to 3×10⁻¹⁴ |
| 0b | c direction spread ≤ 0.5% | **PASS**: 0.31% |
| 1 | Scattering falls at each step (> 2 SE) | **PASS**: 15.3, 10.9, 10.1 and 9.3 SE |
| 2 | Single power law (χ²/dof < 3) | **PASS**: 1.09 |
| 3 | Exponent replicates | **PASS**: p = 3.68 ± 0.07 vs 3.74 ± 0.04 (0.8 SE apart) |
| 4 | Size replicates (fitted 1/ℓ at point 3 within 10%) | **PASS**: ratio 0.999 |
| 5 | Directions agree at each point (2 SE) | **FAIL at all 5 points**: 3.3, 2.6, 2.2, 2.4 and 6.1 SE |
| 6 | Clocks within 3%; total conserved | **PASS**: −2.21% to −0.21%; 3×10⁻¹² |

## The numbers (s4; mean of 6 runs each)

| Wavelength | Launch | Mean free path | vs s3 | [100] vs [1,2,3] | Settled speed (shortfall) |
|---|---|---|---|---|---|
| 17.8 bonds | 0.804c | 580 bonds | 0.97× | +8% | 0.7482c (7.0%) |
| 21.2 bonds | 0.751c | 1,042 bonds | 1.01× | +14% | 0.7110c (5.4%) |
| 25.2 bonds | 0.692c | 1,986 bonds | 1.06× | +15% | 0.6626c (4.2%) |
| 29.9 bonds | 0.627c | 3,876 bonds | 1.05× | +17% | 0.6067c (3.3%) |
| 35.6 bonds | 0.561c | 7,815 bonds | 0.94× | +28% | 0.5458c (2.7%) |

- **Point by point** the mean free paths agree with s3 within 1.3 SE at every wavelength.
- **Joint fit of both samples** (shared exponent, separate sizes): p = 3.73 ± 0.03, χ²/dof 1.10, size ratio s4/s3 = 0.99.

## Reading
- **The scattering law belongs to the glass recipe, not to one sample.** Same exponent, same size, same point-by-point values.
  - This holds even though s4 is measurably more strained (+12% energy) and has a 1.6% higher c.
  - The combined result is 1/ℓ ∝ k₀^3.73 ± 0.03 over wavelengths of 18–36 bonds.
- **Settled speed shortfall replicates too:** 7.0% → 2.7% here vs 6.8% → 2.5% on s3.
- **Direction agreement fails, in a telling way** (post hoc reading):
  - On s4, [100] scatters more than [1,2,3] at every wavelength, by 8% rising to 28%.
  - On s3, the difference went the other way at long wavelengths (−9% at point 5) and was near zero at short ones.
  - So the sign is set by the sample, and the size grows with wavelength in both.
  - That is what a finite sample's own random anisotropy would do. Each box is one 326 Å glass repeated 2×2×2, and the longer waves feel the larger-scale structure, of which a single box contains fewer independent pieces.
  - It points to a sample-size effect rather than a preferred direction in the glass. But that is an interpretation, not a test; it could be tested with larger untiled samples.

## Next (to be pre-registered)
- **Direction test:** a larger untiled glass (for example, 13.8M nodes built directly, not tiled), to see whether the direction difference shrinks as the sample grows. The lean code makes this feasible on GitHub.
- **Wavelength range:** 4j on a 3×3×3 box, up to about 53 bonds.
- **Speed vs wavelength:** vary η at fixed speed.

## Files
- `run_4i_rep.py`, `build_glass4.py`; `res_4i_rep/` (31 jobs; also `results/exp4i_rep/` on GitHub).
- Glass s4 in `networks/` on GitHub, with its measured c (`glass_n1728000_s4_c.json`) and build log (`results/glass4/build.log`).
- `fit_4i_rep.json`, `fig4i_rep.py`, `exp4i_rep_results.png`.
