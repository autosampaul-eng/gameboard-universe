# Experiment 4i: Results — scattering vs wavelength
Gameboard Universe. Run October 10, 2026, as 35 GitHub jobs (workflow `exp4i`, 50 minutes). Design: `DESIGN_exp4i.md`, written before any 4i run and before the 4h rerun results were seen. Tracker: 4h Amendment 1.

## Registered verdict

| # | Criterion | Result |
|---|---|---|
| 0 | Diamond control at all 5 wavelengths | **PASS**: Δs = 0.0000; \|1/τ\| ≤ 3×10⁻¹³ per tick; particle-only speed equals whole-box speed to ≤ 0.001% |
| 1 | Scattering falls at each step to a longer wavelength (> 2 SE) | **PASS**: steps of 20.1, 15.9, 14.1 and 15.7 SE |
| 2 | Single power law (χ²/dof < 3) | **PASS**: χ²/dof = 1.28 |
| 3 | Exponent p (1/ℓ ∝ k₀^p) | **p = 3.74 ± 0.04.** None of the reference values 2, 3 or 4 is within 2 SE; 4 is 6.4 SE away |
| 4 | Same in every direction (each point within 2 SE) | **FAIL at point 5:** 1.0, 0.1, 0.4 and 1.1 SE at points 1–4; 3.8 SE at point 5 ([1,2,3] scatters 10% more) |
| 5 | Clocks within 3%; total conserved | **PASS**: −2.04% to −0.21%; 3×10⁻¹² |

## The numbers

| Wavelength | Launch speed | Mean free path ℓ | In wavelengths | Settled particle speed |
|---|---|---|---|---|
| 17.8 bonds (48.7 Å) | 0.800c | 1,533 Å = 561 bonds | 31 | 0.7455c |
| 21.2 bonds (57.9 Å) | 0.746c | 2,870 Å = 1,049 bonds | 50 | 0.7067c |
| 25.2 bonds (68.8 Å) | 0.686c | 5,741 Å = 2,099 bonds | 83 | 0.6575c |
| 29.9 bonds (81.9 Å) | 0.621c | 11,086 Å = 4,053 bonds | 135 | 0.6014c |
| 35.6 bonds (97.4 Å) | 0.555c | 20,031 Å = 7,324 bonds | 206 | 0.5406c |

- **Each value** is the mean of 6 glass runs: 3 starts × 2 directions. The SE of 1/ℓ is 2–3%.
- **Doubling the wavelength makes the mean free path 13× longer.** That is 1/ℓ ∝ k₀^3.74, or ℓ/λ ∝ λ^2.74.
- **Per tick** rather than per distance, the exponent is 4.20 ± 0.04 (χ²/dof 0.99). The two differ because the slower particles also cover less distance per tick.
- **Step-by-step exponents** are 3.62, 4.00, 3.80 and 3.41. There's no clear trend; the last step is the lowest.

## Reading
- **Long-wavelength particles are much better protected.** This is the behaviour of waves scattering off disorder much smaller than themselves. Each 19% increase in wavelength cuts the scattering rate roughly in half.
- **The exponent, 3.74, is close to the Rayleigh value of 4 but measurably below it.** Two possible reasons, neither tested:
  - speed and wavelength change together in these runs;
  - 18–36 bonds may not yet be fully in the long-wavelength regime.
- **Diamond stays perfectly clean at every wavelength**, so the scattering belongs to the glass's disorder, not to the rule or the measurement.
- **The direction mismatch at the longest wavelength is small** (10%). It may be a statistical-sample effect of a single glass. At that wavelength the particle covers more of the box, so it samples fewer independent regions. Not resolved here.

## Not scored
- **The settled speed is below the launch speed** at every wavelength. The shortfall shrinks steadily with wavelength: 6.8%, 5.3%, 4.2%, 3.2% and 2.5%.
  - A launch speed is set from the crystal-like ansatz with the measured glass c. So the shortfall says the glass's own speed-vs-wavelength relation differs from the ansatz, and the difference gets smaller as the particle gets larger.
  - This is a property worth its own test.
- **Drag** is about −0.11 to −0.16% per 100 ticks at all points. Halves are within 0.35%.
- **The fraction of the particle cut off by the cone/shell** is 0.06–0.27% for every point. Scale-matching worked; compare 3–7% at 0.5c in 4h.
- **Comparison with 4h:**
  - The 4h figure of about 2,100 Å came from the late linear rate (ticks 250–600).
  - The 4i definition (exponential decay, ticks 100–600) gives 1,533 Å for the same particle. The scattering itself is identical (Γ = 5.46 vs 5.44).

## Limits
- **One glass sample, tiled 2×2×2. Wavelength range 2×.** p = 3.74 is a local slope over 18–36 bonds, not a law.
  - A real particle at 10¹⁶ spacings is 15 orders of magnitude beyond this range, so no prediction at that scale is claimed.
  - What can be said: if the slope stayed near 3.7 at long wavelengths, ℓ/λ would grow by about 2.7 orders of magnitude per decade of wavelength.
- **Speed and wavelength are not separated.**

## Next (to be pre-registered)
- **Replicate on an independent glass sample.** This checks p and the point-5 direction mismatch.
- **Separate speed from wavelength:** vary η at fixed speed so k₀ changes with v held constant.
- **Extend the range** with a 3×3×3 tiled box (46.7M nodes). Scaling from the 4.7 GB peak at 13.8M nodes gives about 16 GB, right at the GitHub runner's limit. So the code's memory use would need to be cut first.

## Files
- `run_4i.py`; `res_4i/` (35 job results and logs; also `results/exp4i/` on GitHub); `fit_4i.json`; `fig4i.py`; `exp4i_results.png`.
