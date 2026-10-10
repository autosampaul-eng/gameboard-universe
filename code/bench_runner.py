"""Machine check (no physics results): load the 13.8M tiled glass, time ticks of the update rule, report memory."""
import time, resource, numpy as np
from run_4g import load_tiled, adj
from exp4 import step
t0 = time.time(); pos, nb, rev, d0, LB, LB0, bip = load_tiled('glass'); A = adj(nb)
print(f'loaded N={len(nb)} box={LB:.1f} A in {time.time()-t0:.0f}s', flush=True)
phi = (np.random.default_rng(0).standard_normal((len(nb), 4)) + 0j); out = np.zeros_like(phi)
ce = np.cos(0.15) / 2; em = np.exp(-0.15j)
step(phi, nb, rev, ce, em, out)          # compile
t1 = time.time()
for _ in range(20): step(phi, nb, rev, ce, em, out); phi, out = out, phi
print(f'tick time at 13.8M nodes: {(time.time()-t1)/20:.3f}s', flush=True)
x = np.random.default_rng(1).standard_normal(len(nb)) + 0j; t1 = time.time()
for _ in range(20): x = A @ x
print(f'sparse matvec time: {(time.time()-t1)/20:.3f}s', flush=True)
print(f'peak memory: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f} MB', flush=True)
