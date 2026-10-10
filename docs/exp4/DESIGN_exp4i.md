# Experiment 4i: how glass scattering depends on wavelength (pre-registered)
Gameboard Universe. Written October 10, 2026, before any 4i runs and before the Experiment 4h rerun results were seen. Rule (η = 0.15), packet construction, measured c, the 652 Å tiled networks, the scattering measure and the corrected particle-only tracker (4h Amendment 1) are as in Experiment 4h.

## Question
- In Experiment 4h, glass scattered a 0.8c particle (wavelength 17.8 bonds) steadily, with a mean free path of about 2,100 Å. At 0.5c it scattered 16–28× less.
- **How does the mean free path ℓ grow with wavelength?** This is the number that decides whether the effect could be invisible for real particles, which are about 10¹⁶ node spacings across or more.

## Scale-matched particles
- To compare wavelengths fairly, every particle has the **same shape relative to its wavelength.** The packet width σ is scaled with the wavelength, so its spread of directions is the same for all (σ_k/k₀ = 0.097).
- So the 30° cone and the ±50% shell hold the same share of every particle. This fixes the 0.5c trimming problem found in 4h.
- Wavelength is changed through the speed, at fixed η. Five wavelengths, each 2^(1/4) ≈ 1.19× longer than the last, spanning a factor of 2:

| Point | Speed v/c | k₀ (1/Å) | Wavelength | σ (Å) |
|---|---|---|---|---|
| 1 | 0.8000 | 0.12907 | 48.7 Å = 17.8 bonds | 40.0 |
| 2 | 0.7463 | 0.10854 | 57.9 Å = 21.2 bonds | 47.6 |
| 3 | 0.6860 | 0.09127 | 68.8 Å = 25.2 bonds | 56.6 |
| 4 | 0.6213 | 0.07675 | 81.9 Å = 29.9 bonds | 67.3 |
| 5 | 0.5547 | 0.06454 | 97.4 Å = 35.6 bonds | 80.0 |

- **Why the range stops at 2×:** σ = 80 Å is the widest particle that fits in the 652 Å box. At the box edge its initial probability is below 10⁻⁴. Longer wavelengths need a bigger box than these machines can hold.
- **A trade-off stated in advance:** because wavelength is changed through speed, speed and wavelength change together. The fits treat the rate per unit distance as a function of k₀.

## Runs (35 GitHub jobs)
- **Glass:** 600 ticks; at each of the 5 points, [100] and [1,2,3], with the first 3 starting points from 4f–4h. 30 runs.
- **Diamond control:** at each of the 5 points, [100], first starting point. 5 runs.

## Measures
- **Coherent fraction:** f(t) = 1 − Δs(t), with Δs as in 4h (shell ±50% of k₀, forward cone 30°).
- **Decay rate 1/τ:** slope of −ln f(t) against t, fitted over ticks 100–600.
  - The fit starts at tick 100 so the packet has settled; it uses the full run after that.
- **Mean free path:** ℓ = v_p τ, with v_p the particle-only settled speed (ticks 250–600) in Å per tick.
- **Also reported:** the linear rate Γ from 4h (slope of Δs over ticks 250–600).
- Each point's 1/ℓ is the mean over its 6 runs ± SE.
- **Power-law fit:** weighted least squares of ln(1/ℓ) against ln k₀ over the 5 points, giving the exponent p (1/ℓ ∝ k₀^p) ± SE, and χ²/dof.

## Pre-registered criteria
0. **Diamond control** at each point:
   - |Δs| < 0.5 points at all ticks, and |1/τ| < 2×10⁻⁵ per tick;
   - corrected particle-only speed within 0.1% of the same run's whole-box speed.

   A point whose control fails is left out of the fits, and that is reported. If 2 or more points fail, the experiment is reported as not working.
1. **Falls with wavelength:** 1/ℓ is lower at each step to a longer wavelength, by more than 2 SE of the difference, in all 4 steps.
2. **Single power law:** χ²/dof < 3 (3 degrees of freedom). Otherwise, "not a single power law over this range".
3. **Exponent:** p ± SE is reported, along with which reference values (2, 3, 4) lie within 2 SE.
   - 4 is the value for waves scattering off disorder much smaller than the wavelength (Rayleigh-type).
   - This is a reading, not a pass/fail.
4. **Same in every direction:** at each point, [100] and [1,2,3] give 1/ℓ within 2 SE of each other.
5. **Clocks:** every glass run within 3% of relativity, using the particle-only speed. **Total** conserved to 10⁻⁹.

## Reported, not scored
- The particle-only speed at each point vs the launch speed. In 4h, a 0.8c launch settled at 0.744c.
- Drag and halves of the particle-only speed.
- Off-shell fraction.
- Comparison with the 4h 0.5c result (not scale-matched).

## Limits stated in advance
- One glass sample, tiled. Wavelength range only 2× (17.8–35.6 bonds). An exponent from it is a local slope, not a law valid at 10¹⁶ spacings.
- With fixed η, slower particles are also closer to the mass scale (v/c 0.55 at the longest wavelength). The rate per distance may depend on speed as well as wavelength, and these runs can't separate the two.
- Re-scattering of scattered waves back into the cone is ignored in the decay fit. Its effect is largest for the fastest point (up to about 6.7% of the scattered share, the cone's fraction of all directions).

Any outcome is reported as found.
