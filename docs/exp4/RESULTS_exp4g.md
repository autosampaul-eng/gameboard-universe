# Experiment 4g: Results — 0.8c particles in a 652 Å tiled box
Gameboard Universe. Run October 9, 2026. Design: `DESIGN_exp4g.md`, written before any 4g run.
- **Networks:** raw glass s3 tiled 2×2×2 (13.8M nodes, 652 Å box), and diamond tiled the same way.
- **Settings:** 600 ticks, 6 glass starting points per direction, at the same positions in Å as 4f.
- **Code check:** the new code reproduced the old code exactly on the 216k test glass. Packets and tracking agreed to 10⁻¹³.
- **Run times:** about 31 minutes per glass run; peak memory 5.3 GB.

## Registered verdict

| # | Criterion | Result |
|---|---|---|
| 0 | Spreading check (tracking methods agree within 0.3%) | **FAIL**: all 12 glass runs differ by 0.5–2.1%. Diamond agrees to 0.001%. |
| 1 | Diamond control | **PASS**: drift and wobble below 0.0002% in both directions |
| 2 | No drag | **FAIL**: [100] −1.58 ± 0.52, [1,2,3] −1.17 ± 0.49 (% per 100 ticks) |
| 3 | Steady (halves) | **FAIL**: [100] −3.06%, [1,2,3] −2.48% |
| 4 | Directions agree | **FAIL**: 0.7336 vs 0.7568, a gap of 3.07% (6.3 SE) |
| 5 | Clocks within 3% | **PASS**: −2.36% to −1.25% |
| 6 | Total conserved | **PASS**: 6.5×10⁻¹² |

- Criterion 0 failed, so, as pre-registered, the glass speed results are **unreliable**.
- **Comparison with 4f:** from the same starting points, settled speeds changed by up to 1.9%, well over the 0.3% threshold. So the 326 Å box distorted 4f.
- The diamond result is unchanged from the small box (0.7911 vs 0.7910), which confirms the code.

## What the bigger box revealed (post hoc, not pre-registered)
The fault is in the tracking method, not the box size.

- **The glass particle spreads twice as far as the crystal one.** By tick 600 its spread (rms) is 199–217 Å on glass vs 105 Å on diamond. In the 326 Å box this was hidden at about 135 Å, because the box can't hold more.
  - The likely cause: the glass scatters part of the particle into a slowly spreading cloud. Roughly 10–35% of the total, judging from the extra spread; that is a rough estimate.
  - Part of the extra spread may instead be faster spreading of the particle itself. The saved data can't separate the two.
- **The cloud drags the measured centre.** The whole-box centre is the average over particle and cloud together.
  - **Spreading cloud:** while the cloud spreads around the starting point, it holds the centre back, and the measured speed dips.
  - **Box-filling cloud:** once the cloud reaches half the box width and wraps, it fills the box evenly. It then stops pulling, the centre catches up, and the speed overshoots.
- **This predicts the timing, and the timing matches:**

| Box | Cloud reaches half-box (≈ half-box ÷ c) | Measured speed peak |
|---|---|---|
| 326 Å | about tick 105 | ticks 120–180 |
| 652 Å | about tick 210 | about tick 280 ([100]), 320–360 ([1,2,3]) |

  - The dip is deeper in the big box (0.66c vs 0.72c), because the cloud can grow larger before it wraps.
  - The two boxes already disagree at ticks 50–150, before the particle can have reached any box edge. A real particle effect would not depend on box size that early.
- **So the start-up transient and settling reported in 4d and 4e were most likely this tracking artifact, not particle behaviour.** The same probably goes for most of the direction gap:
  - **Whole-box centre:** gap 3.07% (big box), 1.16% (small box). It depends heavily on the box.
  - **Limited-window centre:** gap 0.53 ± 0.27% (2.0 SE, big box). It is smaller and not clearly real.
- **What survives:** the clocks follow relativity at the measured speed, within 2.4%, in every run and in both boxes.
  - The clock comes from the particle's wave phase, which the cloud barely affects.

## What this means for the speed and direction questions
- The registered drag, steadiness and direction verdicts from 4d–4g all rest on the whole-box centre. On glass that centre is contaminated by the scattered cloud, so those verdicts don't measure the particle cleanly.
- **The real physical question is now:** how much of a fast particle does the glass scatter, and how fast?
  - If a steady fraction scatters every tick, the glass gives particles a mean free path. A real particle would then lose coherence just by moving through empty space.
  - That would mark the glass's rest frame as special, which experiments already constrain tightly.
  - If instead the scattering happens only at the start (the built packet not matching the glass's own travelling state), it's a set-up issue.

## Next (to be pre-registered)
**Experiment 4h:** measure the particle and the cloud separately.
- **Particle position:** the centre of density within a sphere around the peak, about 2.5σ wide, following the peak.
- **Coherent fraction:** the size of the wave-matched sum |Σ S e^{−ik·r}| over time. This is a direct measure of how much of the particle is still a coherent travelling wave.
- **Scattering:** the fraction outside the sphere over time. Does it level off (a start-up loss) or keep growing (an ongoing loss)?
- **Comparisons:** diamond as the zero-scattering control; 0.5c vs 0.8c; the small and big boxes.

## Files
- `run_4g.py`, `test_4g.py`, `run_4g.log`, `res_4g.json`: the registered runs and code check.
- `fig4g.py`, `exp4g_results.png`: the figure.
