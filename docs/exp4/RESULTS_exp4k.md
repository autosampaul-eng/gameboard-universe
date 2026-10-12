# Experiment 4k: Results — direction spread in an untiled glass
Gameboard Universe. Run October 11–12, 2026, on GitHub. Pre-registration and its three amendments: `DESIGN_exp4k.md`. The scoring script (`score_4k.py`) was committed before any results arrived.

## Short answer
At the long wavelength, where 4i and 4i-R found the large direction differences, the untiled glass shows no direction spread beyond run-to-run noise. The tiled glasses show 4–7%. This fits the finite-sample explanation.

At the short wavelength, every glass, tiled or not, shows only a small spread (2–3.5%), about the size of the noise. There it does not shrink, and there was little to shrink.

## The new glass (s5)
- **Same recipe as s3/s4,** seed 5, built directly at 13,824,000 nodes in a 652 Å box. That is 8× the independent material of the tiled boxes.
- **Three amendments, all logged before any particle run:**
  1. **A guided final repair step** closed the last 2 stuck gaps (644 moves) after the 400-pass repair stopped.
  2. **A 0.1 Å per-step move cap** in the relaxation.
     - Without it, the relaxation exploded at step 110–120: a few atoms overshot and the stiff energy ran away.
     - Control: s4 rebuilt from its seed with the cap has identical links and E/atom 7.592908, vs 7.592909 for the published s4.
  3. **Continuing the same relaxation past 5,000 steps.** It converged at step 7,242.

| | s3 | s4 | s5 (untiled) |
|---|---|---|---|
| Energy per atom | 6.79 | 7.59 | 7.37 |
| Bond-length rms | 0.275 | 0.292 | 0.287 |
| Bond-angle rms | 30.4° | 30.8° | 30.7° |
| Largest residual force | 9.6×10⁻⁴ | 8.0×10⁻⁴ | 3.4×10⁻⁴ |
| Speed limit c | 1.5495 | 1.5740 | 1.5659 |
| c spread, [100] vs [1,2,3] | ≤ 0.4% | 0.31% | 0.35% |

One caveat on the 0.35% c spread: the [1,2,3] wave available in this box has a 3.7× higher wavenumber than the [100] wave. Some of the 0.35% is therefore ordinary dispersion rather than direction. The true direction spread is no larger than 0.35%.

## Registered verdict

| # | Criterion | Result |
|---|---|---|
| 0a | Lean builder and resumed relaxation bit-identical to the originals | **PASS** |
| 0b | s5 largest force < 10⁻³ | **PASS**: 3.4×10⁻⁴ |
| 0c | s5 c direction spread ≤ 0.5% | **PASS**: 0.35% |
| **Q** | Long wavelength (35.6 bonds): R = A_s5 / A_tiled | **R = 0.00** (bootstrap 68%: 0.00–0.54). **"Shrinks as a finite-sample effect"** |
| Q1 | Short wavelength (17.8 bonds) | **R = 1.28** (bootstrap 68%: 0.73–1.58). **"Does not shrink"** |
| S | Same law: s5 average within 10% of s3/s4 | **PASS**: 1.005 (long), 0.982 (short) |
| — | Clocks within 3%; total conserved | **PASS**: −2.22% to −0.19%; ≤ 3×10⁻¹² |

## The numbers
A is the direction spread with run-to-run noise removed (the registered measure). The raw spread and the noise per direction are shown alongside it.

| Glass | Long wavelength: A | raw spread | noise per direction | Short wavelength: A | raw spread | noise per direction |
|---|---|---|---|---|---|---|
| s3, tiled | 4.4% | 5.0% | 2.4% | 0.0% | 2.1% | 2.3% |
| s4, tiled | 6.5% | 7.0% | 2.5% | 2.8% | 3.5% | 2.1% |
| s5, untiled | 0.0% | 1.7% | 3.2% | 2.5% | 3.2% | 2.0% |

- **Long wavelength, raw spreads:** untiled vs tiled is 1.7% vs 6.1% (root-mean-square of s3 and s4), a ratio of 0.28.
  - The 1/√8 prediction was 0.35.
- **Direction patterns are unrelated between glasses** (correlations −0.75 to +0.61, with no consistent sign). That is what random, sample-specific anisotropy looks like.
- **Mean free paths agree** across all three glasses (long: s5 ℓ = 7,452 bonds vs 7,180 for s3 and 7,822 for s4). The speed shortfall repeats as well: 2.4% long, 7.0% short.

## Reading
- **Long wavelength:** the 10–28% direction differences in 4i/4i-R came from finite tiled samples.
  - Each tiled box repeats one 326 Å piece of glass, so its random lopsidedness at long wavelengths is copied, not averaged.
  - With 8× the independent glass, the spread falls to within noise, by about the predicted amount.
  - No preferred direction in the glass is needed to explain them.
- **Short wavelength (post hoc reading):** all three glasses are at about the noise level (2–3.5% raw, against about 2% noise).
  - The registered "does not shrink" is a correct reading of the rule. In practice it means the measurement can't resolve any short-wavelength spread with 2 starts per direction.
  - Any real spread at 17.8 bonds is a few percent at most, in every sample.
  - The 4i-R short-wavelength failure (+8% at 3.3 SE, [100] vs [1,2,3]) is about 2.5× the noise level here, so it is unlikely to be pure noise, but it does not reappear clearly in the 6-direction view.
- **Limits:**
  - One untiled sample, six directions and two starts.
  - The long-wavelength R has a wide range (0–0.54).
  - A second untiled glass, or more starts, would tighten it.

## Files
- **Code:** `bigglass.py` (lean builder, guided final repair, resumable relaxation with move cap), `glass5.py`, `run_4k.py`, `score_4k.py`, `fig4k.py`.
- **Logs and checks:** `results/glass5/` (build, diagnostic, s4 control, continuation, finish logs).
- **Results:** `res_4k/` (56 jobs), `score_4k.json`, `exp4k_results.png`.
- **The s5 network itself:** the [glass-s5 release on GitHub](https://github.com/autosampaul-eng/gameboard-universe/releases/tag/glass-s5). Its quality and c are in `networks/glass_n13824000_s5_*.json`.
