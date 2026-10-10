# Experiment 4h: how much of a fast particle does the glass scatter, and does it keep scattering? (pre-registered)
Gameboard Universe. Written October 10, 2026, before any 4h runs. Rule, packet construction, measured c, clock and the 652 Å tiled networks are as in Experiment 4g. What's new is how the particle is measured.

## Question
- In Experiment 4g, glass particles spread twice as far as crystal ones (rms about 210 Å vs 105 Å). The whole-box centre was pulled around by what looked like a scattered cloud.
- **Q1:** Is part of the particle being scattered out of its travelling wave? If so, does the scattering stop after the start (a set-up loss), or keep going (a mean free path, which would make the glass's rest frame special)?
- **Q2:** With the cloud removed, does the particle itself move steadily, with no drag and the same speed in both directions?

## How the particle and the cloud are separated
At every 5th tick:
1. **Map onto a grid.** The node totals S_i (the sum of a node's 4 link amplitudes) are spread onto a regular grid with cloud-in-cell weighting. Grid spacing: 652 Å / 128 = 5.09 Å for the big box; 326 Å / 64 for the small box.
2. **Fourier transform.** This gives the wave's spectrum of directions and wavelengths.
3. **The shell:** all waves with wavenumber within ±50% of the particle's carrier wavenumber k₀.
   - Scattering by a static network keeps the frequency, so scattered waves stay near this shell but change direction.
4. **The forward cone:** the part of the shell within 30° of the carrier direction.
   - The particle's own spread of directions is about 5° (1/(2σk₀) with σ = 40 Å), so 30° keeps the whole particle.
5. **Scattered fraction s(t) = 1 − (forward-cone weight ÷ shell weight).**
   - Δs(t) = s(t) − s(0) removes the baseline from the grid mapping and from the glass's position noise.
   - Also reported: the fraction of all weight that is off the shell.
6. **The particle's position:** the forward-cone part is transformed back onto the grid, and its density centre is tracked (nearest-image, as before). This is the particle with the cloud removed.
   - The whole-box centre from 4g is also recorded, for comparison.

The method is checked on the 216k test networks before any registered run.

## Runs (25 GitHub jobs, one run each)
- η = 0.15, σ = 40 Å, 600 ticks. Starting positions are the same in Å as in 4f and 4g.
- **A. Glass, 652 Å box, 0.8c:** [100] and [1,2,3], 6 starting points each. 12 runs.
- **B. Glass, 652 Å box, 0.5c:** [100] and [1,2,3], the first 3 starting points. 6 runs.
- **C. Diamond, 652 Å box:** 0.8c and 0.5c, [100] and [1,2,3]. 4 runs.
- **D. Glass, original 326 Å box, 0.8c:** [100], the first 3 starting points, as a box-size check. 3 runs.

## Measures
- **Late scattering rate Γ:** slope of Δs(t) over ticks 250–600, in percentage points of the shell weight per 100 ticks.
- **Start-up loss:** Δs at tick 250.
- **Mean free path:** if Γ > 0, then ℓ ≈ v × 100 ÷ (Γ/100), in Å, reported for reference.
- **Particle-only speed measures,** as in 4f/4g but from the particle's position:
  - settled speed (ticks 250–600);
  - drag D (trend of 40-tick window speeds, % per 100 ticks);
  - halves (ticks 250–425 vs 425–600).
- **Clock:** carrier frequency − k·v as before, with v the particle-only settled speed.
- Results are given as mean ± SE over starting points.

## Pre-registered criteria
0. **Diamond control:**
   - |Δs(t)| < 0.5 percentage points at every recorded tick;
   - |Γ| < 0.2 per 100 ticks;
   - particle-only settled speed within 0.1% of the 4g whole-box diamond speed (0.7911 / 0.7866 at 0.8c).

   If this fails, the method is reported as not working, and Q1 and Q2 are not judged.
1. **Q1, scattering timing** (glass, 652 Å, 0.8c, each direction):
   - **"Start-up only":** Γ < 0.5 per 100 ticks and within 2 SE of zero, in both directions.
   - **"Continuous":** Γ > 0.5 per 100 ticks and more than 3 SE above zero, in either direction.
   - **Otherwise:** inconclusive.

   At this sensitivity, "start-up only" means any continuing loss is slower than about 0.5% per 100 ticks, a mean free path above about 25,000 Å (about 9,000 bonds). That is a limit of these runs, not proof of zero.
2. **Q2, particle motion** (glass, 652 Å, 0.8c):
   - no drag: |mean D| < 1% per 100 ticks and within 2 SE, each direction;
   - steady: mean halves < 1%, each direction;
   - directions agree: settled speeds differ by < 1% and by < 2 SE.
3. **Clocks:** every glass run within 3% of relativity.
4. **Total:** conserved to 10⁻⁹.

## Reported, not scored
- **0.5c vs 0.8c:** start-up loss and Γ.
- **Box check:** small-box Γ and start-up loss vs the big box (runs D vs the matching runs in A).
- **Off-shell fraction** over time.
- **Whole-box vs particle-only speeds.**

## Limits stated in advance
- One glass sample, tiled. The 30° cone and ±50% shell are fixed choices. A different choice would shift the baseline s(0), but should not create a trend.
- The scattering measure depends on the grid mapping. Δs removes the fixed part, but not any part that changes with the particle's position in the glass. The diamond control and the spread of results across starting points show how large that is.

Any outcome is reported as found.
