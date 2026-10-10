# Memory-lean code and the 3×3×3 box: engineering test
Gameboard Universe, October 10, 2026. Engineering only: no registered glass was used, and no physics verdicts are made here.

## What changed (`lean.py`, `run_big.py`)

| Item | Old (run_4g / run_4h2) | Lean |
|---|---|---|
| Neighbour sums A·v | Sparse matrix copy of the network (56 B/node) | Computed straight from the link table (0 B) |
| Reverse-link table | Full tiled copy (16 B/node) | Base network only, indexed by i mod N₀ |
| Tiled link table | Built as 64-bit, then converted | Built directly as 32-bit |
| Packet build | SciPy CG plus many full-size temporaries | Hand-written CG and Chebyshev, updated in place (same stopping rule) |

- Mathematics and measurements are unchanged.
- The packet build is also about 5× faster.

## Same answers as the old code
Checked on the 216k test glass, tiled 2×2×2, at 4i point 2, 100 ticks:
- scattered share: same to 6×10⁻¹⁴;
- speeds: same to 2×10⁻¹⁵;
- Γ, clock error and 1/τ: same to 10⁻¹³ or better.

## Memory

| Nodes | Old peak | Lean peak |
|---|---|---|
| 13.8M (652 Å box) | 4.7 GB (GitHub smoke test) | 2.8–2.9 GB |
| 46.7M (977 Å box, 3×3×3) | about 16 GB (estimated; would not fit) | **8.7–9.2 GB** (measured on GitHub; the machine has 16 GB) |

## 3×3×3 test runs on GitHub (600 ticks each)

| Run | Nodes | Build | Per tick | Total | Peak | Result |
|---|---|---|---|---|---|---|
| Diamond 3×3×3, 4i point 5 (σ = 80) | 46.7M | 225 s | 0.37 s | 16 min | 9.24 GB | No scattering (max \|Δs\| 0.000); particle-only = whole-box = 0.54946c; clock −0.15% |
| 216k test glass 6×6×6, point 5 | 46.7M | 384 s | 0.59 s | 20 min | 8.71 GB | 1/τ = 2.724×10⁻⁵, Γ = 0.262, particle-only 0.54122c |
| 216k test glass 4×4×4, point 5 (comparison) | 13.8M | 71 s | 0.10 s | 4 min | 2.90 GB | 1/τ = 2.749×10⁻⁵, Γ = 0.262, particle-only 0.54083c |
| 216k test glass 6×6×6, point 1 (σ = 40) | 46.7M | 358 s | 0.57 s | 20 min | 8.71 GB | 1/τ = 6.60×10⁻⁴, Γ = 4.91, particle-only 0.74547c |

- **Same glass, bigger box:** the scattering rate is unchanged (−0.9% in 1/τ, Γ identical) and the speed shifts by +0.07%. So the 652 Å box was not distorting the widest (σ = 80) particle's scattering.
- **Diamond at point 5 reads 0.03% faster** in the 977 Å box than in the 652 Å box used in 4i (0.54946 vs 0.54928). That is within 4i's 0.1% control limit. It suggests a small box-size effect on speed for the widest particles; the scattering rate is unaffected.
- **Room for growth:** at about 190 bytes per node at peak, about 70M nodes fit in 16 GB with a safety margin. A 4×4×4 tiling (110.6M) would not fit without further cuts, such as computing positions on the fly.

## What this makes possible
- A pre-registered Experiment 4j on the main glass, tiled 3×3×3 (977 Å box): wavelengths up to about 53 bonds (σ ≈ 120 Å), extending 4i's range from 2× to 3×.
- About 20–25 minutes per run on GitHub.

## Files
`lean.py`, `run_big.py`, `.github/workflows/bigtest.yml`; results in `results/bigtest/` on GitHub.
