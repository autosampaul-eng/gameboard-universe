# Experiment 4h rerun (Amendment 1): Results
Gameboard Universe. Run October 10, 2026, as 25 GitHub jobs (workflow `exp4h2`). Design: `DESIGN_exp4h.md`. Change: `AMENDMENT_exp4h_1.md`, a probability-weighted particle-only tracker, logged and checked on the 216k-derived test networks before the rerun.

## Registered verdict

| # | Criterion | Result |
|---|---|---|
| 0a | Diamond \|Δs\| < 0.5 points | **PASS** (0.000) |
| 0b | Diamond \|Γ\| < 0.2 per 100 ticks | **PASS** (0.0000) |
| 0c | Diamond particle-only speed within 0.1% of 4g | **PASS:** 0.79114 / 0.78664, identical to 4g to 5 digits (was −0.11% before the fix) |
| Q1 | Scattering timing (glass, 0.8c) | **"Continuous"**: Γ = 5.44 ± 0.07 ([100], 81 SE) and 5.41 ± 0.09 ([1,2,3], 64 SE) per 100 ticks |
| Q2 | No drag (< 1% and within 2 SE) | [100] −0.18 ± 0.10 (1.85 SE) **pass**; [1,2,3] −0.17 ± 0.07 (2.5 SE) **fail** |
| Q2 | Steady (halves < 1%) | **PASS**: −0.29% / −0.25% |
| Q2 | Directions agree | **PASS**: 0.7450 vs 0.7453, gap −0.05 ± 0.10% (0.5 SE) |
| 3 | Clocks within 3% | **PASS**: −2.09% to −0.10% |
| 4 | Total conserved | **PASS**: 3×10⁻¹² |

- **Q1:** glass scatters fast particles continuously. A mean free path of about 2,100 Å (780 bonds, 44 wavelengths), the same in both directions and in both box sizes (326 Å box: 5.53 ± 0.10).
- **Q2:** the particle itself moves at the same speed in both directions and holds steady between halves. A very slight slowdown (about 0.18% per 100 ticks) is statistically clear along [1,2,3]. So the strict no-drag rule fails there.
- **Scattering numbers are unchanged** from the first 4h run (identical to 10⁻¹⁴), because only the tracker changed. Particle-only speeds rose by about 0.1% (0.7442 → 0.7450).

## Not scored
- **At 0.5c the tracker is still biased** by cone trimming: diamond +0.40% / +0.90% against whole-box. This was expected; see the amendment. So the 0.5c direction gap (1.23%) is not meaningful. Experiment 4i's scale-matched particles address this.
- **0.5c scattering:** Γ = 0.34 ± 0.02 ([100]) and 0.19 ± 0.06 ([1,2,3]), unchanged.

## Files
- `run_4h2.py`; `res_4h2/` (25 job results and logs; also `results/exp4h2/` on GitHub).
