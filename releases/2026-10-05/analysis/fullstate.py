#!/usr/bin/env python3
"""Full-state error decomposition for compressed finite-temperature candidates,
evaluated on the saved IDMate Hamiltonians (H2) and preceding-map frames (U1).

Metric: embedded weighted norm  d(A,B) = sqrt(2 * sum_k w_k ||A_k - B_k||_F^2)
for per-spin (half) density matrices, identical to the published d_full.
"""
import gzip, json, math, os, sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
# Saved Hamiltonian records (<material>/H2_U1.json.gz and input.json), available from the authors. Set
# IDMATE_H2_RUNS to their directory; the default is ./h2_records next to this script.
RUNS = os.environ.get('IDMATE_H2_RUNS') or os.path.join(_HERE, 'h2_records')

def read_json(p):
    op = gzip.open if p.endswith('.gz') else open
    with op(p, 'rt') as s:
        return json.load(s)

def matrix(v):
    r, c = int(v['rows']), int(v['cols'])
    a = np.asarray(v['complex_re_im'], float)
    return (a[:, 0] + 1j * a[:, 1]).reshape((r, c), order='F')

def fermi(x):
    # stable Fermi function 1/(1+e^x)
    y = np.exp(-np.abs(x))
    return np.where(x >= 0, y / (1 + y), 1 / (1 + y))

def solve_mu(spectra, weights, tau, target):
    """common mu such that sum_k w_k sum_i f((e-mu)/tau) = target (per spin)."""
    lo = min(e.min() for e in spectra) - 50 * tau - 1.0
    hi = max(e.max() for e in spectra) + 50 * tau + 1.0
    def count(mu):
        return sum(w * fermi((e - mu) / tau).sum() for w, e in zip(weights, spectra))
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if count(mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-15 * max(1.0, abs(mid)):
            break
    mu = 0.5 * (lo + hi)
    return mu, count(mu) - target

def fdd(x, y, tau):
    """first divided difference of f(e)=1/(1+exp(e/tau)), elementwise; x,y shifted by mu."""
    x = np.asarray(x)[:, None]; y = np.asarray(y)[None, :]
    d = x - y
    fx, fy = fermi(x / tau), fermi(y / tau)
    small = np.abs(d) < 1e-9 * tau
    with np.errstate(divide='ignore', invalid='ignore'):
        out = np.where(small, -fx * (1 - fx) / tau, (fx - fy) / np.where(small, 1.0, d))
    return out

def load(material):
    path = os.path.join(RUNS, material)
    inp = read_json(os.path.join(path, 'input.json'))
    hist = read_json(os.path.join(path, 'H2_U1.json.gz'))
    hs = [matrix(h) for h in hist['H2']]
    old = hist['U1_records_residuals_against_H2_not_H1']
    U1 = [matrix(p['U']) for p in old]
    f1 = [np.asarray(p['occupations_f'], float) for p in old]
    w = np.array([p['weight'] for p in inp['kpoints']], float)
    return dict(inp=inp, hs=hs, U1=U1, f1=f1, w=w, tau=float(inp['tau_ha']),
                N=float(inp['target_electrons']), old=old)

def canonical_full(hs, w, tau, Nspin):
    spec, frames = [], []
    for h in hs:
        e, u = np.linalg.eigh(h)
        spec.append(e); frames.append(u)
    mu, dN = solve_mu(spec, w, tau, Nspin)
    D = [(u * fermi((e - mu) / tau)).dot(u.conj().T) for e, u in zip(spec, frames)]
    return D, mu, spec, frames

def dist(A, B, w):
    return math.sqrt(2.0 * sum(wk * np.linalg.norm(a - b) ** 2 for wk, a, b in zip(w, A, B)))

def analyse(data, extra=2, tau=None, ref=None):
    hs, U1, f1, w = data['hs'], data['U1'], data['f1'], data['w']
    tau = data['tau'] if tau is None else tau
    Nspin = data['N'] / 2.0
    if ref is None:
        ref = canonical_full(hs, w, tau, Nspin)
    Dstar, mu_star, spec, frames = ref
    # retained frames (historical rule: last occupied f>1e-12 plus `extra` bands)
    P, HP_eigs, HP_vecs, Qs, QHQ_eigs, QHQ_vecs, E_blocks = [], [], [], [], [], [], []
    for h, u, f in zip(hs, U1, f1):
        occ = np.flatnonzero(f > 1e-12)
        m = min(int(occ[-1]) + 1 + extra, u.shape[1])
        p = u[:, :m]
        hp = p.conj().T @ h @ p
        th, v = np.linalg.eigh(hp)
        pr = p @ v                                  # Ritz vectors
        n = h.shape[0]
        # complement basis
        q = u[:, m:]                                # U1 is unitary -> complement frame
        hq = q.conj().T @ h @ q
        lam, z = (np.linalg.eigh(hq) if q.shape[1] else (np.zeros(0), np.zeros((0, 0))))
        qr = q @ z if q.shape[1] else q
        coup = qr.conj().T @ h @ pr                 # (n-m) x m  coupling in eigenbases
        P.append(pr); HP_eigs.append(th); Qs.append(qr); QHQ_eigs.append(lam); E_blocks.append(coup)
    return analyse_blocks(hs, w, tau, Nspin, P, HP_eigs, Qs, QHQ_eigs, E_blocks, Dstar, mu_star, extra=extra)

def analyse_blocks(hs, w, tau, Nspin, P, HP_eigs, Qs, QHQ_eigs, E_blocks, Dstar, mu_star, extra=None):
    """All estimators for a Ritz candidate given per-k Ritz vectors P (Ritz values HP_eigs), complement eigenvectors Qs
    (eigenvalues QHQ_eigs of the compressed complement Hamiltonian), coupling blocks E = Q^dag H P in those eigenbases,
    and the reference D* (with mu_star). Shared by the saved-Hamiltonian sweeps and the in-loop rerun (inloop.py)."""
    # Ritz candidate with common mu, zero complement
    mu_c, _ = solve_mu(HP_eigs, w, tau, Nspin)
    Dc = [(pr * fermi((th - mu_c) / tau)).dot(pr.conj().T) for pr, th in zip(P, HP_eigs)]
    # canonical state of block-diagonal H0 (common mu over both blocks)
    allspec = [np.concatenate([th, lam]) for th, lam in zip(HP_eigs, QHQ_eigs)]
    mu0, _ = solve_mu(allspec, w, tau, Nspin)
    D0 = [(pr * fermi((th - mu0) / tau)).dot(pr.conj().T) + ((qr * fermi((lam - mu0) / tau)).dot(qr.conj().T) if qr.shape[1] else 0)
          for pr, th, qr, lam in zip(P, HP_eigs, Qs, QHQ_eigs)]
    # in-subspace and complement parts of ||Dc - D0||
    eP = math.sqrt(2.0 * sum(wk * np.linalg.norm(fermi((th - mu_c) / tau) - fermi((th - mu0) / tau)) ** 2
                              for wk, th in zip(w, HP_eigs)))
    eQ = math.sqrt(2.0 * sum(wk * np.linalg.norm(fermi((lam - mu0) / tau)) ** 2 for wk, lam in zip(w, QHQ_eigs)))
    # coupling norms
    c = math.sqrt(sum(wk * np.linalg.norm(E) ** 2 for wk, E in zip(w, E_blocks)))   # sqrt(sum w ||QHP||^2)
    bound_cpl = 2.0 * c / (4.0 * tau)          # embedded metric: sqrt2 (embedding) * sqrt2 (two blocks)
    # first-order Daleckii-Krein estimate (embedded)
    dk2 = 0.0
    for wk, th, lam, E in zip(w, HP_eigs, QHQ_eigs, E_blocks):
        if E.size == 0:
            continue
        g = fdd(lam - mu0, th - mu0, tau)       # (n-m) x m
        dk2 += wk * 2.0 * np.sum(np.abs(g) ** 2 * np.abs(E) ** 2)
    d_dk = math.sqrt(2.0 * dk2)
    # computable complement-occupation bound C_hat >= C and the chemical-potential interval [mu_hat, mu_c]
    # that contains mu_0 (mu_hat solves the Ritz-only count for N - C_hat)
    Chat = sum(wk * lam.size * float(fermi((lam.min() - mu_c) / tau)) for wk, lam in zip(w, QHQ_eigs) if lam.size)
    if Chat < Nspin * (1 - 1e-12):
        mu_hat, _ = solve_mu(HP_eigs, w, tau, Nspin - Chat)
    else:
        mu_hat = None
    mu_lo = mu_hat if mu_hat is not None else min(th.min() for th in HP_eigs) - 40 * tau
    mus = [mu_c] if mu_c - mu_lo < 1e-12 else list(np.linspace(mu_lo, mu_c, 11))
    # cheap estimate: per-Ritz residual norms and lowest complement level only; the divided-difference
    # weight is the supremum over lambda >= lambda_min and over mu in [mu_hat, mu_c] (grid supremum)
    ch2 = 0.0
    for wk, th, lam, E in zip(w, HP_eigs, QHQ_eigs, E_blocks):
        if E.size == 0:
            continue
        rho = np.linalg.norm(E, axis=0)          # ||Q H p_i|| per Ritz vector
        lmin = lam.min()
        g = np.zeros(th.size)
        for mu in mus:
            ygrid = (lmin - mu) + np.concatenate([np.linspace(0, 60 * tau, 601), np.geomspace(60 * tau + 1e-3, 1e3, 200)])
            g = np.maximum(g, np.max(np.abs(fdd(th - mu, ygrid, tau)), axis=1))
        ch2 += wk * 2.0 * np.sum(g ** 2 * rho ** 2)
    d_cheap = math.sqrt(2.0 * ch2)
    # computable second-order near-Fermi estimate from the same data (rho_ki, lam_min,k):
    # |delta theta_ki| <~ min(rho_ki, rho_ki^2 / (lam_min,k - theta_ki)); shift by this maximum and re-solve mu
    THs = []
    for th, lam, E in zip(HP_eigs, QHQ_eigs, E_blocks):
        if E.size == 0:
            THs.append(th); continue
        rho = np.linalg.norm(E, axis=0); gap = lam.min() - th
        # Kato-Temple-type cap: a Ritz value moves by at most its residual norm rho_ki
        THs.append(th - np.where(gap > 0, np.minimum(rho, rho ** 2 / np.maximum(gap, 1e-300)), 0.0))
    mu_s, _ = solve_mu(THs, w, tau, Nspin)
    eP2_est = math.sqrt(2.0 * sum(wk * np.linalg.norm(fermi((th - mu_c) / tau) - fermi((ts - mu_s) / tau)) ** 2
                                  for wk, th, ts in zip(w, HP_eigs, THs)))
    # exact orthogonal blocks of Dc - D* (embedded): P block, Q block, off-diagonal block
    sPP = sQQ = sOFF = 0.0
    for wk, pr, qr, dc, ds in zip(w, P, Qs, Dc, Dstar):
        A = dc - ds
        Bpr = pr.conj().T @ A @ pr; Bqq = qr.conj().T @ A @ qr; Bqp = qr.conj().T @ A @ pr
        sPP += wk * np.linalg.norm(Bpr) ** 2; sQQ += wk * np.linalg.norm(Bqq) ** 2; sOFF += wk * 2 * np.linalg.norm(Bqp) ** 2
    blk_P, blk_Q, blk_off = math.sqrt(2 * sPP), math.sqrt(2 * sQQ), math.sqrt(2 * sOFF)
    # second-order level shifts of Ritz and complement levels, re-solved common mu
    th2, lam2 = [], []
    for th, lam, E in zip(HP_eigs, QHQ_eigs, E_blocks):
        if E.size == 0:
            th2.append(th); lam2.append(lam); continue
        den = lam[:, None] - th[None, :]
        den = np.where(np.abs(den) < 1e-12, 1e-12, den)
        W = np.abs(E) ** 2 / den
        th2.append(th - W.sum(axis=0)); lam2.append(lam + W.sum(axis=1))
    allspec2 = [np.concatenate([a, b]) for a, b in zip(th2, lam2)]
    mu2, _ = solve_mu(allspec2, w, tau, Nspin)
    eP2 = math.sqrt(2.0 * sum(wk * np.linalg.norm(fermi((th - mu_c) / tau) - fermi((t2 - mu2) / tau)) ** 2
                              for wk, th, t2 in zip(w, HP_eigs, th2)))
    eQ2 = math.sqrt(2.0 * sum(wk * np.linalg.norm(fermi((l2 - mu2) / tau)) ** 2 for wk, l2 in zip(w, lam2) if l2.size))
    # practical complement estimate from the lowest complement level only:
    # e_Q^2 = 2 sum_k w_k sum_j f(lam_j - mu)^2 <= 2 sum_k w_k n_Q,k f(lam_min,k - mu)^2   (mu -> candidate mu_c)
    eq2 = 0.0
    for wk, lam in zip(w, QHQ_eigs):
        if lam.size:
            eq2 += wk * lam.size * float(fermi((lam.min() - mu_c) / tau)) ** 2
    eQ_est = math.sqrt(2.0 * eq2)
    # in-subspace estimate from the complement occupation count (trace-defect argument)
    # in-subspace estimate: mu_hat <= mu0 <= mu_c and occupations are monotone in mu, so
    # ||f(th-mu_c) - f(th-mu_hat)|| >= e_P; the mu_hat -> -inf limit gives ||f(th-mu_c)|| (always valid)
    eP_cap = math.sqrt(2.0 * sum(wk * np.linalg.norm(fermi((th - mu_c) / tau)) ** 2 for wk, th in zip(w, HP_eigs)))
    eP_est_crude = math.sqrt(2.0) * Chat / math.sqrt(min(w))      # sum w d^2 <= (C/w_min) C
    if mu_hat is not None:
        eP_est = math.sqrt(2.0 * sum(wk * np.linalg.norm(fermi((th - mu_c) / tau) - fermi((th - mu_hat) / tau)) ** 2
                                     for wk, th in zip(w, HP_eigs)))
        eP_est = min(eP_est, eP_cap)
    else:
        mu_hat, eP_est = float('-inf'), eP_cap
    indicator = math.sqrt(eP_est ** 2 + eQ_est ** 2 + d_cheap ** 2)
    indicator2 = math.sqrt(indicator ** 2 + eP2_est ** 2)        # I_1 plus second-order near-Fermi estimate
    # exact coupling-induced change ||D*(H0) - D*(H)||
    d_cpl_exact = dist(D0, Dstar, w)
    d_full = dist(Dc, Dstar, w)
    rigorous = math.sqrt(eP ** 2 + eQ ** 2) + bound_cpl
    gapQ = min(lam.min() for lam in QHQ_eigs if lam.size) - mu0 if any(lam.size for lam in QHQ_eigs) else float('nan')
    # largest Ritz value retained relative to mu
    topP = max(th.max() for th in HP_eigs) - mu0
    return dict(extra=extra, tau=tau, m=[p.shape[1] for p in P], n=[h.shape[0] for h in hs],
                d_full=d_full, eP=eP, eQ=eQ, c_w=c, bound_cpl=bound_cpl, rigorous=rigorous,
                d_dk=d_dk, d_cheap=d_cheap, eQ_est=eQ_est, eP_est=eP_est, eP_est_crude=eP_est_crude, rigorous_computable=math.sqrt(eP_est ** 2 + eQ_est ** 2) + bound_cpl, eP2_est=eP2_est, indicator2=indicator2, C_hat=Chat, mu_hat=mu_hat, indicator=indicator, blk_P=blk_P, blk_Q=blk_Q, blk_off=blk_off, eP2=eP2, eQ2=eQ2, pyth2=math.sqrt(eP2 ** 2 + eQ2 ** 2 + d_dk ** 2), d_cpl_exact=d_cpl_exact, mu_c=mu_c, mu0=mu0, mu_star=mu_star,
                gapQ=gapQ, topP=topP)

if __name__ == '__main__':
    for mat in ['si', 'graphene', 'al']:
        data = load(mat)
        r = analyse(data, extra=2)
        print(mat, 'tau', data['tau'], 'n', r['n'][:3], 'm', r['m'][:3])
        for k in ['d_full', 'd_cpl_exact', 'd_dk', 'd_cheap', 'bound_cpl', 'rigorous', 'eP', 'eQ', 'c_w', 'gapQ', 'topP']:
            print('   %-12s %.6e' % (k, r[k]))

# ---------------------------------------------------------------- densities
def grid_phases(data, k):
    inp = data['inp']
    nx, ny, nz = [int(x) for x in inp['grid_dims']]
    ids = np.arange(nx * ny * nz)
    frac = np.column_stack(((ids % nx + 0.5) / nx, ((ids // nx) % ny + 0.5) / ny, (ids // (nx * ny) + 0.5) / nz))
    idx = np.asarray(data['old'][k]['plane_wave_indices'], float)
    kf = np.asarray(data['old'][k]['k_frac'], float)
    return np.exp(2j * np.pi * frac.dot((idx + kf).T))      # npts x n

def density(data, Ds, phases=None):
    """rho(r) = (2/Omega) sum_k w_k phi^dag D phi  on the cell-centred grid (replay convention)."""
    vol = float(data['inp']['volume_bohr3'])
    rho = 0.0
    for k, (wk, D) in enumerate(zip(data['w'], Ds)):
        Ph = phases[k] if phases is not None else grid_phases(data, k)
        val = np.real(np.sum(np.conj(Ph) * (Ph @ D.T), axis=1))
        rho = rho + wk * 2.0 * val / vol
    return rho

def eps_rho(data, rho_a, rho_b):
    vol = float(data['inp']['volume_bohr3']); N = float(data['N'])
    return float(vol * np.linalg.norm(rho_a - rho_b) / math.sqrt(rho_a.size) / N)
