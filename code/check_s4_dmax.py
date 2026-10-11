"""Amendment 2 check: rebuild s4 (ncell 60, seed 4) with the lean builder and the 0.1 A move cap; compare with the published s4."""
import numpy as np, json, crn, bigglass as B, time
ncell, seed = 60, 4; a = 4 * crn.D / np.sqrt(3); Lb = ncell * a; n = 8 * ncell ** 3
t = time.time(); pos, nb = B.random_network_lean(n, Lb, seed, passes=400)
nb0 = np.load('../networks/glass_n1728000_s4_raw_nb.npy').astype(np.int64); p0 = np.load('../networks/glass_n1728000_s4_raw_pos.npy')
print('links identical to published s4:', np.array_equal(nb, nb0), flush=True)
v = np.zeros_like(pos); st = B.new_state(); B.relax(pos, nb, Lb, v, st, dmax=0.1)
q = crn.quality(pos, nb, Lb)
out = dict(links_identical=bool(np.array_equal(nb, nb0)), positions_identical=bool(np.array_equal(pos, p0)), max_pos_diff=float(np.abs(pos - p0).max()),
           fire_iterations=int(st[3]), converged=float(st[6]), clamped_atom_steps=int(st[7]), E_per_atom=q['E_per_atom'], published_E_per_atom=json.load(open('../networks/glass_n1728000_s4_raw.json'))['E_per_atom'], secs=time.time() - t)
print(json.dumps(out), flush=True); json.dump(out, open('../results/glass5/check_s4_dmax.json', 'w'))
