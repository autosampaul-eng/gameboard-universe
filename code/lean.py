"""Memory-lean versions of the tiling, packet build and update (same mathematics as run_4g / exp4).
- No sparse matrix: A@v is computed directly from the link table (A_ij = links i->j).
- Reverse-slot table kept at base size and indexed by i % N0 (tiling preserves slot structure).
- Tiled link table built directly as int32.
- Packet build: in-place Chebyshev recursion and a hand-written conjugate-gradient solver (same stopping rule as
  scipy's cg: ||r|| <= rtol * ||b||), with fused kernels so no full-size temporaries are created."""
import numpy as np, gc, time
from numba import njit, prange

@njit(parallel=True, cache=True)
def tile_nb(nb0, shift, T, NB):
    n = nb0.shape[0]
    for q in prange(T ** 3):
        cx = q % T; cy = (q // T) % T; cz = q // (T * T)
        for i in range(n):
            for a in range(4):
                tx = (cx + shift[i, a, 0]) % T; ty = (cy + shift[i, a, 1]) % T; tz = (cz + shift[i, a, 2]) % T
                NB[q * n + i, a] = nb0[i, a] + n * (tx + T * ty + T * T * tz)

@njit(parallel=True, cache=True)
def tile_pos(p0, Lb, T, P):
    n = p0.shape[0]
    for q in prange(T ** 3):
        cx = q % T; cy = (q // T) % T; cz = q // (T * T)
        for i in range(n):
            P[q * n + i, 0] = p0[i, 0] + Lb * cx; P[q * n + i, 1] = p0[i, 1] + Lb * cy; P[q * n + i, 2] = p0[i, 2] + Lb * cz

@njit(parallel=True, cache=True)
def check_rev(NB, rev0):
    N = NB.shape[0]; N0 = rev0.shape[0]; bad = 0
    for i in prange(N):
        for a in range(4):
            if NB[NB[i, a], rev0[i % N0, a]] != i: bad += 1
    return bad

def tile(p0, nb0, rev0, Lb, T):
    p0 = p0 % Lb
    r = p0[nb0] - p0[:, None, :]; r -= Lb * np.round(r / Lb)
    shift = np.round((p0[:, None, :] + r - p0[nb0]) / Lb).astype(np.int64)
    N = len(p0) * T ** 3
    NB = np.empty((N, 4), np.int32); tile_nb(nb0.astype(np.int64), shift, T, NB)
    P = np.empty((N, 3)); tile_pos(p0, Lb, T, P)
    bad = check_rev(NB, rev0.astype(np.int32)); assert bad == 0, 'tiled reverse-slot table inconsistent'
    return P, NB, Lb * T

@njit(parallel=True, cache=True)
def nbsum(v, nb, out):
    for i in prange(nb.shape[0]):
        s = 0j
        for a in range(4): s += v[nb[i, a]]
        out[i] = s

@njit(parallel=True, cache=True)
def cheb_step(T1, T0, nb, c, acc):
    """T0 <- 2*(A/4)T1 - T0 ; acc += c*T0 (in place)."""
    for i in prange(nb.shape[0]):
        s = 0j
        for a in range(4): s += T1[nb[i, a]]
        t = 0.5 * s - T0[i]; T0[i] = t; acc[i] += c * t

@njit(parallel=True, cache=True)
def mop(p, tmp, nb, Ap):
    """Ap = (16 - A A) p using tmp; returns <p, Ap> (real for this Hermitian operator)."""
    N = nb.shape[0]
    for i in prange(N):
        s = 0j
        for a in range(4): s += p[nb[i, a]]
        tmp[i] = s
    d = 0.0
    for i in prange(N):
        s = 0j
        for a in range(4): s += tmp[nb[i, a]]
        z = 16 * p[i] - s; Ap[i] = z; d += (p[i].conjugate() * z).real
    return d

@njit(parallel=True, cache=True)
def upd_xr(x, r, p, Ap, al):
    rs = 0.0
    for i in prange(x.shape[0]):
        x[i] += al * p[i]; z = r[i] - al * Ap[i]; r[i] = z; rs += z.real * z.real + z.imag * z.imag
    return rs

@njit(parallel=True, cache=True)
def upd_p(p, r, be):
    for i in prange(p.shape[0]): p[i] = r[i] + be * p[i]

@njit(parallel=True, cache=True)
def norm2(v):
    s = 0.0
    for i in prange(v.shape[0]): s += v[i].real * v[i].real + v[i].imag * v[i].imag
    return s

def cg(nb, b, x, rtol=1e-9, maxiter=5000):
    """Solve (16 - A A) x = b in place (x holds the start). b is overwritten with the residual."""
    tmp = np.empty_like(b); Ap = np.empty_like(b)
    mop(x, tmp, nb, Ap); b -= Ap; r = b                     # r = b - M x
    bn = np.sqrt(norm2(r + Ap))                             # ||b|| (one transient temp)
    p = r.copy(); rs = norm2(r); it = 0
    while np.sqrt(rs) > rtol * bn and it < maxiter:
        d = mop(p, tmp, nb, Ap); al = rs / d
        rs2 = upd_xr(x, r, p, Ap, al); upd_p(p, r, rs2 / rs); rs = rs2; it += 1
    del tmp, Ap, p
    return it if np.sqrt(rs) <= rtol * bn else -it

@njit(parallel=True, cache=True)
def link_from_xy(x, y, nb, out):
    for i in prange(nb.shape[0]):
        for a in range(4): out[i, a] = x[i] + y[nb[i, a]]

@njit(parallel=True, cache=True)
def gauss_S0(pos, c0, LB, kv, sigma, S0):
    for i in prange(pos.shape[0]):
        rx = pos[i, 0] - c0[0]; rx -= LB * np.round(rx / LB)
        ry = pos[i, 1] - c0[1]; ry -= LB * np.round(ry / LB)
        rz = pos[i, 2] - c0[2]; rz -= LB * np.round(rz / LB)
        ph = kv[0] * rx + kv[1] * ry + kv[2] * rz
        S0[i] = np.exp(-(rx * rx + ry * ry + rz * rz) / (4 * sigma * sigma)) * (np.cos(ph) + 1j * np.sin(ph))

@njit(parallel=True, cache=True)
def out_err(phi, S0):
    m = 0.0
    for i in prange(phi.shape[0]):
        v = abs(phi[i, 0] + phi[i, 1] + phi[i, 2] + phi[i, 3] - S0[i])
        m = max(m, v)
    return m

@njit(cache=True)
def in_err(phi, nb, X):
    N = nb.shape[0]; acc = np.zeros(N, np.complex128)
    for i in range(N):
        for a in range(4): acc[nb[i, a]] += phi[i, a]
    m = 0.0
    for i in range(N):
        v = abs(acc[i] - X[i])
        if v > m: m = v
    return m

def build_packet(pos, nb, d0, LB, c0, kvec, w, sigma, eta, mu_of_k, bip=None):
    """Same construction as run_4g.build_packet (vvec is unused there at t=0), memory-lean."""
    N = len(nb); S0 = np.empty(N, np.complex128); gauss_S0(pos, c0, LB, kvec, sigma, S0)
    ce = np.cos(eta); nC = 400
    xs = np.cos(np.pi * (np.arange(nC) + 0.5) / nC); fv = np.sqrt(1 - ce ** 2 * xs ** 2)
    coef = np.array([2.0 / nC * (fv * np.cos(j * np.pi * (np.arange(nC) + 0.5) / nC)).sum() for j in range(nC)]); coef[0] /= 2
    nuse = int(np.nonzero(np.abs(coef) > 1e-13)[0].max()) + 1
    T0 = S0.copy(); T1 = np.empty_like(S0); nbsum(S0, nb, T1); T1 *= 0.25
    acc = coef[0] * T0; acc += coef[1] * T1
    for j in range(2, nuse):
        cheb_step(T1, T0, nb, coef[j], acc); T0, T1 = T1, T0
    AS0 = T0; nbsum(S0, nb, AS0); del T1; gc.collect()          # reuse buffer
    # X = e^{i eta} (ce/2 AS0 - S1), S1 = ce/4 AS0 + i acc  ->  X = e^{i eta} (ce/4 AS0 - i acc)
    acc *= -1j; acc += (ce / 4) * AS0; acc *= np.exp(1j * eta); X = acc; del AS0; gc.collect()
    rhs = np.empty_like(S0); nbsum(X, nb, rhs); rhs *= -1; _add(rhs, S0, 4.0)
    if bip is not None:
        m = (bip @ rhs) / len(bip); _axpy_real(rhs, bip, -m)
    lam = 4 - mu_of_k(d0, kvec); Xc = np.exp(1j * eta) * (np.cos(eta) / 2 * lam - np.exp(1j * w)); den = 16 - lam ** 2
    p = 1 / 8 if den < 1e-9 else (4 - lam * Xc) / den
    x = S0 * p
    info = cg(nb, rhs, x); del rhs; gc.collect()
    y = np.empty_like(S0); nbsum(x, nb, y); y *= -1; _add(y, X, 1.0); y *= 0.25      # y = (X - A x)/4
    phi = np.empty((N, 4), np.complex128); link_from_xy(x, y, nb, phi); del x, y; gc.collect()
    e_out = float(out_err(phi, S0)); e_in = float(in_err(phi, nb, X))
    del S0, X; gc.collect()
    return phi, e_in, e_out, info

@njit(parallel=True, cache=True)
def _add(a, b, c):
    for i in prange(a.shape[0]): a[i] += c * b[i]

@njit(parallel=True, cache=True)
def _axpy_real(a, bip, m):
    for i in prange(a.shape[0]): a[i] += m * bip[i]

@njit(parallel=True, cache=True)
def step(phi, nb, rev0, ce, em, out):
    N = phi.shape[0]; N0 = rev0.shape[0]
    for i in prange(N):
        s = phi[i, 0] + phi[i, 1] + phi[i, 2] + phi[i, 3]; ii = i % N0
        for a in range(4):
            out[nb[i, a], rev0[ii, a]] = ce * s - em * phi[i, a]
