# Gameboard Universe

Toy models for a discrete-lattice ("gameboard") picture of space: every node has 4 links, and particles are waves hopping node to node under one local update rule. The question is how much of special relativity (a speed limit, moving clocks running slow) comes out of the rules alone, on a crystal (diamond) and on a random network (glass).

Working practice: every experiment is **pre-registered** (the `DESIGN_*.md` written before the runs), debugging uses only the small 216k networks, amendments are logged, and verdicts are reported as found. Analyses done after the fact are labelled post hoc.

## Layout
- `code/`: the update rule and analysis (`exp4.py`), the glass builder (`crn.py`, `build.py`), tiling (`exp2.py`), and experiment runners.
- `networks/`: the glass networks (positions in Å; neighbour lists, 4 per node).
  - `glass_n1728000_s3_raw`: the main glass, 1.73M nodes.
  - `glass_n216000_s1_*`: the 216k debug glasses.
- `docs/exp4/`: designs and results for Experiments 4–4g.
- `.github/workflows/`: jobs that run the code on GitHub's machines.

## Status (October 2026)
- Moving clocks follow relativity within about 2–2.5% up to 0.8c, on both diamond and glass, using each network's measured speed limit.
- On glass, a fast particle loses part of itself to a scattered cloud. That spoils the simple "centre of everything" position measurement. Experiment 4h (next) measures the scattering directly.
