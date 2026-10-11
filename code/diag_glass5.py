"""Diagnose the s5 relaxation failure: rebuild s5 exactly, then report bond lengths, initial energy, and the first FIRE steps in detail."""
import numpy as np, crn, bigglass as B, time
from numba import njit
ncell, seed = 120, 5; a = 4 * crn.D / np.sqrt(3); Lb = ncell * a; n = 8 * ncell ** 3
log = lambda s: print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)
pos, nb = B.random_network_lean(n, Lb, seed, passes=400)
@njit
def bonds(pos, nb, Lb):
    n = pos.shape[0]; out = np.empty((n, 4))
    for i in range(n):
        for k in range(4):
            j = nb[i, k]; s = 0.0
            for c in range(3):
                d = pos[j, c] - pos[i, c]; d -= Lb * np.round(d / Lb); s += d * d
            out[i, k] = np.sqrt(s)
    return out
@njit
def per_atom(pos, nb, Lb):
    n = pos.shape[0]; g = np.zeros_like(pos); e = np.empty(n)
    for i in range(n): e[i] = crn.atom_eg(i, pos, nb, Lb, g)
    return e, g
bl = bonds(pos, nb, Lb); log(f'bonds: mean {bl.mean():.3f} max {bl.max():.3f}  >2D {(bl > 2*crn.D).sum()}  >9A {(bl > 9).sum()}')
log(f'positions: min {pos.min():.3f} max {pos.max():.3f} Lb {Lb:.3f}; nb min {nb.min()} max {nb.max()} dtype {nb.dtype}')
e, g = per_atom(pos, nb, Lb); f = np.sqrt((g ** 2).sum(1)); o = np.argsort(-e)
log(f'initial E/atom {e.mean():.4f}  max atom E {e.max():.3e}  fmax {np.abs(g).max():.3e}  atoms with E>1e3: {(e > 1e3).sum()}')
for i in o[:8]: log(f'  atom {i} E {e[i]:.3e} |f| {f[i]:.3e} nb {nb[i].tolist()} bonds {np.round(bl[i],3).tolist()}')
for q in [50, 90, 99, 99.9, 99.99]: log(f'  E percentile {q}: {np.percentile(e, q):.3f}')
v = np.zeros_like(pos); st = B.new_state(); gg = np.zeros_like(pos)
for k in range(30):
    B.fire_chunk(pos, nb, Lb, gg, v, st, 1e-3, 5000, 10)
    e, g = per_atom(pos, nb, Lb); bl = bonds(pos, nb, Lb)
    log(f'it {int(st[3])} E/atom {st[4]/n:.4f} fmax {st[5]:.3e} dt {st[0]:.4f} | bond max {bl.max():.2f} >9A {(bl > 9).sum()} | atom E max {e.max():.3e} at {int(e.argmax())}, E>1e3: {(e > 1e3).sum()}')
