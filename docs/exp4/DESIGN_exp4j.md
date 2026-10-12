# Experiment 4j: separating wavelength from speed (pre-registered)
Gameboard Universe. Drafted October 12, 2026; approved by Paul and registered the same day, before any 4j run. **The changes from the reviewed draft are listed in the "Registration notes" section at the end; all were made before any run.**

Everything not mentioned here is as in Experiment 4i: the update rule, packet construction, particle-only tracker, scattering measure, starting points, 600 ticks, and the decay fit over ticks 100–600.

## Why this experiment
In 4i and 4i-R, a particle's settled speed came out below its launch speed: 7.0% short at 17.8 bonds, 2.7% short at 35.6 bonds. Two things are unknown:

1. **Does it fade fast enough at long wavelengths?**
   - Every lattice slows short waves (dispersion). For an ordinary lattice the slowdown fades as 1/λ² (∝ k², with k = 2π/λ). That is harmless at real particle sizes, which are around 10¹⁶ node spacings or more.
   - If it fades more slowly, particles of different energy would travel at measurably different speeds, which real observations rule out very tightly.
   - Taken at face value, the 4i numbers fall as about k^1.4.
2. **Is it wavelength or speed?** In 4i, wavelength was changed by changing speed (η fixed at 0.15), so longer waves were also slower (0.80c → 0.55c). The two effects can't be separated from 4i.

The same confound affects the scattering law. The step-by-step slopes in 4i rose from 3.4 to 4.1 toward long wavelengths, and part of that could be speed rather than wavelength.

**The fix:** the wavenumber is set by k₀ = γβη/c, so changing the mass parameter η lets us change wavelength at a fixed speed, and speed at a fixed wavelength.

## A correction found while drafting (stated before any 4j run)
- **Part of the "shortfall" is just packet shape.** A packet contains a spread of directions (σ_k/k₀ = 0.097 in every run), and its sideways parts don't move forward at full speed. So even in a perfect continuum, its center moves slightly slower than the launch speed.
- **Computed** for a continuum, by averaging the true relativistic group velocity over the packet's spread of wavenumbers: 0.92% at 0.80c down to 0.58% at 0.55c.
- **This roughly matches the diamond controls in 4i,** which were 1.0–1.1% short with zero scattering.
- **So 4j measures the excess shortfall** against this computed packet speed, not against the launch speed:

\[ \delta = 1 - \frac{v_p}{v_{pkt}} \]

  - v_p is the particle-only settled speed (ticks 250–600, as in 4i).
  - v_pkt is the packet-averaged continuum speed. It is computed from each glass's own measured c, the run's η and k₀, and the packet's spread. It is computed numerically and stored with each job.
  - δ = 0 would mean the glass behaves exactly like relativity for that packet.

## Grid
- **Speeds:** 0.800c, 0.686c and 0.555c (the 4i points 1, 3 and 5 speeds).
- **Wavelengths:** 7 values, each 2^(1/4) longer than the last: the 5 from 4i plus 2 longer ones, for a span of 2.83×.
- **Packet width** σ ∝ wavelength, as in 4i.
- **η** is set to give each speed at each wavelength. Values below are for s3 (c = 1.5495); s4 uses its own c (1.5740), so its η is 1.6% higher at the same point.

| j | k₀ (1/Å) | Wavelength (bonds) | σ (Å) | η at 0.800c | η at 0.686c | η at 0.555c |
|---|---|---|---|---|---|---|
| 1 | 0.12907 | 17.8 | 40.0 | **0.150** | 0.212 | 0.300 |
| 2 | 0.10853 | 21.2 | 47.6 | 0.126 | 0.178 | 0.252 |
| 3 | 0.09127 | 25.2 | 56.6 | 0.106 | **0.150** | 0.212 |
| 4 | 0.07675 | 29.9 | 67.3 | 0.089 | 0.126 | 0.178 |
| 5 | 0.06454 | 35.6 | 80.0 | 0.075 | 0.106 | **0.150** |
| 6 | 0.05427 | 42.3 | 95.1 | 0.063 | 0.089 | 0.126 |
| 7 | 0.04563 | 50.3 | 113.1 | 0.053 | 0.075 | 0.106 |

The three bold cells (η = 0.150) repeat 4i points 1, 3 and 5 exactly. They serve as a check that the bigger box changes nothing.

## Box and runs
- **Box:** every run uses the 3×3×3 tiled box (977 Å, 46.7M nodes) with the memory-lean code. The longest packets (σ = 113 Å) do not fit in the 652 Å box under 4i's rule (initial probability below 10⁻⁴ at the box edge).
- **Glass runs:** at each of the 21 grid points, 6 runs: glasses s3 and s4 × 3 directions.
  - The directions are [100] (start 1), [1,2,3] (start 2) and [3,1,2] (start 3).
  - That makes 126 runs.
- **Diamond control:** one run per grid point, [100] and start 1, for 21 runs.
- **Regression job:** diamond at η = 0.15 and 0.555c in the 3×3×3 box, run through the unchanged default code path. It must reproduce the earlier engineering test (particle-only speed 0.5494640716399202c) to 10⁻¹⁰.
  - The diamond grid cell with the same settings (0.555c, j = 5) runs through the new per-run-η path and must match it to 10⁻¹⁰ too.
- **Total: 148 GitHub jobs** in one workflow run (`exp4j.yml`, `run_4j.py`). About 15–20 minutes each.

## Code change (done and checked before registering)
- **η is now a per-run setting** in the lean runner (`run_big.one`) and in the tracker weighting (`run_4h2.tracker_weight`). Both default to the fixed 0.15, so earlier experiments are unaffected.
  - The clock error is recomputed with the run's own η.
- **Check: passed.** With η = 0.15 passed explicitly, the changed code gave bit-identical results to the previous code on the 216k test networks (2×2×2, 300 ticks, 4i point 5), for one glass run and one diamond run.
  - All 48 recorded fields were identical; only the timing fields differed.
- **Sanity runs at the extreme η values** (test networks only, 300 ticks):
  - Diamond at η = 0.300 gave no scattering, particle-only and whole-box speeds equal, clock −0.61%, total conserved to 10⁻¹³.
  - Test glass at η = 0.053 ran normally (clock +1.36%, total conserved to 10⁻¹⁴). Its 326 Å box is too small for that packet, so its speeds mean nothing.
- **Packet-speed reference (tabulated before any run):**
  - v_pkt is the continuum group velocity averaged by Gauss–Hermite quadrature over a 3D Gaussian wavenumber distribution of std 1/(2σ) per component about k₀.
  - Because σ·k₀ is the same in every glass cell, it depends only on speed.

| Speed | v_pkt / launch, s3 and s4 (all 7 wavelengths) | diamond (reference) |
|---|---|---|
| 0.800c | −0.918% | −0.704% |
| 0.686c | −0.785% | −0.602% |
| 0.555c | −0.583% | −0.448% |

  - Diamond differs because its c is lower, so the same speed and η give it a larger k₀ (wavelengths of 42.6–120.5 Å, vs 48.7–137.6 Å for the glasses) and a narrower relative spread.
- **A second version, v_pkt_spec, is recorded but not scored.** It averages over the actual initial packet's tracker-weighted spectrum.
  - It includes the tracker's own weighting, and on the test networks it differed from v_pkt by about 1%.
  - It is reported alongside δ as a check on how much the reference choice matters.

## Measures (per grid point; mean ± SE over the 6 glass runs)
- **δ:** the excess speed shortfall, defined above.
- **1/ℓ:** the scattering rate per distance, as in 4i.
- **Joint fits over all 21 points,** weighted least squares, one for each (x = ln(k₀/k_ref), k_ref = 0.07675 Å⁻¹, the middle wavelength of 29.9 bonds):

\[ \delta = A \,(k_0/k_{ref})^{q}\,(\gamma\beta)^{s} \quad \text{(fitted on the plain values of } \delta\text{)} \]

\[ \ln(1/\ell) = a' + p\,x + b\,x^2 + s' \ln(\gamma\beta) \]

- **q** is how fast the excess shortfall fades with wavelength at fixed speed; ordinary dispersion gives 2.
- **s and s′** measure the effect of speed at fixed wavelength.
- **p and b** give the scattering slope at 29.9 bonds and whether it bends. The local slope is p + 2b·x, so b < 0 means steeper at long wavelengths.
- All fits report χ²/dof.

## Pre-registered criteria
0. **Setup:**
   - **(a) Code:** the code check and the regression job pass.
   - **(b) Diamond controls,** at each grid point (as in 4i):
     - |Δs| < 0.5 points at all ticks;
     - |1/τ| < 2×10⁻⁵ per tick;
     - particle-only speed within 0.1% of the whole-box speed.

     A failing point is dropped from the fits and reported. If 4 or more points fail, the experiment is reported as not working.
   - **(c) Box check:** at the three η = 0.15 points, 1/ℓ and v_p agree with 4i (s3 runs) and 4i-R (s4 runs) within 2 combined SE.

     If any fails, the box size matters, and Q1–Q3 are reported but labeled unreliable.
- **Q1 (main): how fast does the excess shortfall fade?** Uses q from the joint fit (18 degrees of freedom).
  - q within 2 SE of 2, with χ²/dof < 3: **"ordinary lattice dispersion"**.
  - q below 2 by more than 2 SE: **"fades more slowly than an ordinary lattice"**. This would be a problem for the model at physical scales, and would need a dedicated follow-up.
  - q above 2 by more than 2 SE: **"fades faster than an ordinary lattice"**.
  - χ²/dof ≥ 3: **"not a single power law"**, reported with the slope at each speed.
- **Q2: does the scattering law bend?**
  - b more than 2 SE below 0: **"steepens toward long wavelengths"**. Reported with the local slope at 42.3 bonds (j = 6, the middle of j = 5–7) and whether 4 lies within 2 SE of it.
  - Otherwise: **"no bend detected"**, reported with p.
- **Q3: does speed matter at a fixed wavelength?** For δ (s) and for scattering (s′) separately:
  - |s| > 2 SE: **"depends on speed"**, reported with its size;
  - otherwise: **"no speed dependence detected"**.

  No direction is predicted.
- **Clocks:** every glass run within 3% of relativity (using the particle-only speed and the run's η). **Total** conserved to 10⁻⁹.

## Reported, not scored
- δ and 1/ℓ at every grid point, by glass and by direction.
- Separate fits for each speed (q, p).
- Whether δ looks like it levels off at a nonzero floor (a fit of δ = A k^q + B). The range is too short to score this, but it matters.
- **Diamond's own δ and its slope q,** as a reference: what an ordered crystal does with the same packets.
- Drag (speed change over the run), halves, and off-shell fraction.

## Limits stated in advance
- **Tiled glass.** 4k showed that 2×2×2 tiling produces sample-specific direction differences of 4–7% at long wavelengths, and 3×3×3 tiling repeats the same 326 Å piece 27 times. Averaging 3 directions × 2 glasses reduces this but doesn't remove it.
  - An untiled 977 Å glass (46.7M nodes) can't be built on GitHub's free runners: the builder's sort step alone needs about 18 GB.
- **Range is 2.83× in wavelength.** Every exponent is still a local slope, not a law valid at 10¹⁶ spacings.
- **Large η approaches the lattice scale.** At η = 0.30 the particle's Compton wavelength is about 12 bonds, so the slowest/shortest cell is the most lattice-affected. That is part of what Q3 measures.
- **Low scattering at the long end.** At j = 7 and 0.555c, only about 0.6% of the packet is expected to scatter over the run. Run-to-run noise there is expected to be about 10–15% per run, about 5% for the 6-run mean.
- **v_pkt assumes a perfect continuum** with the glass's long-wavelength c. Any difference from it, including the glass's own small change of c with wavelength (0.1% between k = 0.01 and 0.04), counts as excess shortfall. That is intended: δ measures the total departure from relativity.

Any outcome is reported as found.

## Registration notes (changes from the reviewed draft, all before any run)
1. **Fit of δ on plain values** instead of ln δ. The sanity runs showed δ can be near zero or negative in a cell, where a log fit breaks. With δ well above its noise, the two fits agree.
2. **Wavenumbers measured from the middle wavelength** (k_ref = 0.07675, 29.9 bonds), so p is the slope at the middle. The long-end slope is evaluated at 42.3 bonds.
3. **Regression check made explicit:** one job uses the default path and the matching diamond cell uses the new path. Both must match the earlier engineering test to 10⁻¹⁰.
4. **One workflow run** of 148 jobs instead of two batches.
5. **Scoring code:** `score_4j.py` is committed with this registration, before any run.
   - It was dry-run on made-up data and recovered the planted values: q = 2.00 (planted 2), s = 0.50 (0.5), p = 3.71 (3.7), b = −0.25 (−0.3) and s′ = 0.00 (0).
6. **Known small box effect (stated before the box check):** in diamond, the 3×3×3 box gave a particle-only speed 0.03% above the 2×2×2 box at 0.555c (0.54946c vs 0.54928c). If the glass v_p box check fails by a similar amount, that is the likely cause; it still counts as a fail under criterion 0c.
