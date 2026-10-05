"""Split the exact diagonal blocks of Dc - D* into diagonal/off-diagonal parts (Ritz/complement eigenbases)
and compare with the second-order formula
<p_i|D2|p_i> = f'(th_i)(dth_i - mu2) - sum_j [f(th_i) - f(lam_j)] |E_ji|^2/(lam_j - th_i)^2 ."""
import math, numpy as np, fullstate as fs
def blocks(mat, b):
    d = fs.load(mat); tau = d['tau']; w = d['w']; Ns = d['N'] / 2
    Dstar, mu_s, _, _ = fs.canonical_full(d['hs'], w, tau, Ns)
    P, TH, Q, LAM, E = [], [], [], [], []
    for h, u, f in zip(d['hs'], d['U1'], d['f1']):
        occ = np.flatnonzero(f > 1e-12); m = int(occ[-1]) + 1 + b
        p = u[:, :m]; th, v = np.linalg.eigh(p.conj().T @ h @ p); pr = p @ v
        q = u[:, m:]; lam, z = np.linalg.eigh(q.conj().T @ h @ q); qr = q @ z
        P.append(pr); TH.append(th); Q.append(qr); LAM.append(lam); E.append(qr.conj().T @ h @ pr)
    mu_c, _ = fs.solve_mu(TH, w, tau, Ns)
    mu0, _ = fs.solve_mu([np.r_[a, c] for a, c in zip(TH, LAM)], w, tau, Ns)
    # second-order chemical potential from re-solving with shifted levels
    th2 = [th - np.sum(np.abs(Ek) ** 2 / (lam[:, None] - th[None, :]), axis=0) for th, lam, Ek in zip(TH, LAM, E)]
    lam2 = [lam + np.sum(np.abs(Ek) ** 2 / (lam[:, None] - th[None, :]), axis=1) for th, lam, Ek in zip(TH, LAM, E)]
    mu2tot, _ = fs.solve_mu([np.r_[a, c] for a, c in zip(th2, lam2)], w, tau, Ns)
    sPd = sPo = sQd = sQo = s_shift = s_full = 0.0
    for wk, pr, th, qr, lam, Ek, Ds in zip(w, P, TH, Q, LAM, E, Dstar):
        Dc = (pr * fs.fermi((th - mu_c) / tau)) @ pr.conj().T
        A = Dc - Ds
        Bp = pr.conj().T @ A @ pr; Bq = qr.conj().T @ A @ qr
        sPd += wk * np.sum(np.abs(np.diag(Bp)) ** 2); sPo += wk * (np.linalg.norm(Bp) ** 2 - np.sum(np.abs(np.diag(Bp)) ** 2))
        sQd += wk * np.sum(np.abs(np.diag(Bq)) ** 2); sQo += wk * (np.linalg.norm(Bq) ** 2 - np.sum(np.abs(np.diag(Bq)) ** 2))
        x = (th - mu0) / tau; f = fs.fermi(x); fp = -f * (1 - f) / tau
        dth = -np.sum(np.abs(Ek) ** 2 / (lam[:, None] - th[None, :]), axis=0)
        fl = fs.fermi((lam - mu0) / tau)
        transfer = -np.sum((f[None, :] - fl[:, None]) * np.abs(Ek) ** 2 / (lam[:, None] - th[None, :]) ** 2, axis=0)
        shift = fp * (dth - (mu2tot - mu0))
        # Dc - D* diagonal (P block) ~ (f_c - f0) - (shift + transfer);  f_c = f0 here
        s_shift += wk * np.sum(shift ** 2); s_full += wk * np.sum((shift + transfer) ** 2)
    e = lambda s: math.sqrt(2 * s)            # embedded
    print(mat, 'b=%d: P diag %.3e off %.3e | Q diag %.3e off %.3e | 2nd-order P-diag: shift-only %.3e, full %.3e' % (
        b, e(sPd), e(sPo), e(sQd), e(sQo), e(s_shift), e(s_full)))
blocks('si', 2); blocks('graphene', 2); blocks('graphene', 16); blocks('al', 2)
