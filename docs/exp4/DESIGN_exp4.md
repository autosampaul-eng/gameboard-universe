# Experiment 4: The GBU's own update rule (pre-registered)
Gameboard Universe. Written October 8, 2026, before any Experiment 4 runs.

## Paul's decisions this rule is built on (Oct 8)
- **Update rule:** selection with memory, accepted as a fair reading of Axiom 7.
- **Velocity value:** a direction.
- **Conserved:** both the total value and the energy.
- **Phase:** node values may carry a phase (amplitudes, not plain probabilities).

## The rule (GBU-native, no borrowed physics)
- **State:** each node holds 4 complex amplitudes, one on each of its links. The amplitude on a link is the part of the node's value that is heading out along that link. That is the "velocity value as a direction".
- **Each tick, two steps:**
  1. **Selection (mixing):** at each node, the 4 amplitudes are re-split among the 4 directions by the same fixed rule everywhere:
     - new amplitude on link a = (cos η / 2) × (sum of all 4) − e^(−iη) × (old amplitude on a)
     - η is the only free number. It sets how strongly a node's value turns rather than carries straight on.
  2. **Move:** every amplitude moves one link along its direction and arrives at the neighbour.
- **Conservation:** the total Σ|amplitude|² is conserved exactly, because the mixing is unitary. Energy (the walk's frequency content) is also conserved exactly, because the rule is the same every tick.
- **Speed limit:** nothing can move more than one link per tick.
- **Proved before running** (`check_recurrence.py`, error 10⁻¹⁵): the node totals S obey S(t+1) + S(t−1) = (cos η / 2) × (sum of the neighbours' S) exactly.
  - On large scales that is the relativistic wave equation for a particle of mass η (in units of 1/tick), with a speed c set by the network.
  - So on smooth large scales the rule is relativistic by construction. **What this experiment tests is whether the real networks (raw glass vs diamond crystal) spoil it, using the GBU's own rule for the first time.**

## Networks
- **Glass:** `glass_n1728000_s3_raw`. 1,728,000 nodes, 4 links each, box 325.6 Å, 30° angle spread. Unrepeated.
- **Crystal:** diamond with 60×60×60 cells. 1,728,000 nodes, same bond length, box 325.6 Å.

## Stage A: speed limit c and its direction dependence (geometry only)
- c is computed from the network for long waves, c = lim θ/k.
- 13 directions, at k = 0.02, 0.04, 0.064 Å⁻¹.
- **PASS A:** c agrees across directions within 1% at k ≤ 0.064 Å⁻¹ (wavelength ≥ 98 Å).

## Stage B: the particle at rest
- **Particle:** a wave packet of the node totals, Gaussian width 40 Å, mass η = 0.15 per tick, built on the slow-moving branch.
- **Clock:** the phase at the packet centre turns at the rest rate.
- **PASS B:** rest clock rate = η within 0.5%.

## Stage C: moving particles
- **Launch:** packets at v/c = 0.3, 0.5, 0.7, 0.8. Carrier k = γmv/c², with m = η, c from Stage A, and γ = 1/√(1 − v²/c²).
- **Directions:** [100] and a generic direction [1,2,3]/√14.
- **Run:** 250 ticks.
- **Measured:** speed v from the centre of |amplitude|², and the clock (the rate the phase turns at the moving centre, with the carrier removed).
- **PASS C:**
  - (a) moving clock ÷ rest clock is within 3% of √(1 − v²/c²), using the measured v and c, at every speed;
  - (b) speed changes by less than 1% between the first and second halves of the run (no drag);
  - (c) the two directions agree within 1% in speed and clock;
  - (d) the total is conserved to 10⁻⁹;
  - (e) on the glass, at most 1% of the total is scattered out of the packet over the run.

## Stage D: the classical check (4a)
- **Setup:** the same rule with plain probabilities instead of amplitudes, using the squared weights of the mixing step and the same starting packet.
- **Prediction:** no clock (there is no phase to turn), and the packet spreads out with no lasting forward motion.
- This checks the code and confirms why the phase is needed. It is reported, not scored.

## Limits stated in advance
- These are free particles, so they spread out slowly, as real quantum particles do. A particle that holds itself together needs interactions, which come later.
- Length contraction can't be tested on free packets, because their shape is set by how they are made. Clocks and drag can be tested.
- One sample of each network.

Any outcome is reported as found.

## Amendment 1 (Oct 8): measurement fixes found in debug runs on the 216k-node test networks, before any registered run
Pass criteria are unchanged. Only how the particle is built and how its clock is read have changed.

1. **Building the particle.**
   - The local recipe for setting each node's 4 link amplitudes was not exact. It mixed in the opposite-frequency branch, and the measured speed and clock came out wrong. Error in the node totals: 0.19 on diamond, 0.29 on glass.
   - New method: the next-tick node totals are made exactly one-branch, S1 = e^{iΩ(A)} S0 with cos Ω = cos η · A/4. This is evaluated with a Chebyshev series on the network itself. The link amplitudes are then solved exactly by conjugate gradients.
   - Errors in the debug runs: 1e-5 or smaller (diamond), 1e-10 (glass).
2. **Reading the clock.**
   - The phase at the centre of a spreading packet also turns because of the spreading itself (the Gouy phase, about 3c²/(4ησ²)). That made the rest clock read 5.6% fast on a σ = 25 Å packet in a perfect crystal, so the centre-phase reading fails Stage B whatever network it runs on.
   - **Primary clock (new):** the frequency ω of the packet's carrier wave, measured from the phase of its k-component, minus k·v using the measured speed. For an exact relativistic particle this equals η/γ.
   - At rest this reading is exact by construction, so **PASS B becomes a code check, not a test.**
   - **Secondary clock (the original centre-phase reading):** still reported, as moving ÷ rest.
3. **Speed.** Still measured from the centre of |amplitude|², as registered. Two known biases are reported, not corrected:
   - Part of the link amplitude sits in modes that don't move (link patterns whose node totals are zero), which pulls the measured centre back a few percent.
   - A packet built from a spread of wavelengths moves at their average speed, not at the speed of the central wavelength.
   - Criterion C(a) already uses the measured v, so neither bias affects the pass/fail test.
