# Experiment 4c: is the glass clock error set by wavelength? (pre-registered)
Gameboard Universe. Written October 8, 2026, before any 4c runs. Rule and measurement are as in Experiment 4, with Amendment 1.

## Question
- On every glass, the clock error fitted err% ≈ a·(k·b)², with a ≈ −47 for raw glass s3 (b = own mean bond, 2.735 Å).
- But in Experiments 4 and 4b, k and the speed v rose together, since k = γηv/c². So the error might be driven by wavelength (the lattice bending the wave law) or by speed itself.
- Changing the mass η separates the two. At the same speed:
  - halving η halves k, so a wavelength effect drops **4×**;
  - doubling η doubles k, so a wavelength effect rises **4×**;
  - a pure speed effect wouldn't change.

## Network
- Raw glass s3 (`glass_n1728000_s3_raw`), the same network as in Experiment 4.

## Runs
- **Packets:** width σ = 50 Å, 120 ticks, speed fitted from tick 20 on.
  - This changes from Experiment 4 (40 Å, 250 ticks) because the light particle spreads faster. It keeps the packet inside the 326 Å box.
- **Set L (light), η = 0.075:**
  - rest;
  - 0.3c, 0.5c, 0.7c, 0.8c;
  - directions [100] and [1,2,3].
- **Set R (reference), η = 0.15:** 0.7c and 0.8c, both directions. Same packet width and run length as set L.
- **Set H (heavy), η = 0.30:** 0.3c and 0.5c, both directions. Higher speeds would put k·b far outside the range already tested.

## Predictions from the wavelength law, a = −47

| Set | v/c | k·b | Predicted error |
|---|---|---|---|
| L | 0.7 | 0.125 | −0.74% |
| L | 0.8 | 0.170 | −1.36% |
| R | 0.7 | 0.250 | −2.94% |
| R | 0.8 | 0.340 | −5.45% |
| H | 0.3 | 0.161 | −1.2% |
| H | 0.5 | 0.295 | −4.08% |

(Calculated from k = γηv/c² with c = 1.607 Å/tick. The exact values are recomputed in the analysis.)

## Pre-registered verdict
- **Main measure:** the ratio Q = (L error at 0.8c) ÷ (R error at 0.8c), averaged over the two directions. Wavelength law: Q = 0.25. Pure speed effect: Q = 1.
  - **"Wavelength law confirmed":** Q ≤ 0.40, **and** a fitted to set L is within ±30% of a fitted to set R.
  - **"Speed-driven":** Q ≥ 0.75.
  - **Anything in between:** inconclusive.
- **Heavy check:** the set H errors are within ±30% of a_R·(k·b)². Reported as pass or fail.
- **Experiment 4 criteria:** also reported per run. These are free particles, so the speed-change and direction criteria are noisier over 120 ticks.

## Limits stated in advance
- At η = 0.075 and 0.3c, the wavelength (about 430 Å) is longer than the box. The clock is still readable from a localised packet, but this point is the least reliable.
- Set L's errors are small (≤ 1.4%), near the measurement noise of about 0.15–0.2%. Q relies mainly on the 0.8c point.

Any outcome is reported as found.

## Correction (logged after only the rest run had finished)
- The 0.7c rows of the prediction table were hand-calculation slips. They have been replaced with the script's exact values. The verdict rules are unchanged.
