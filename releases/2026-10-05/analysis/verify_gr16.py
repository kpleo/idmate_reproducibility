"""Verify the graphene b=16 second-order near-Fermi mechanism quoted in the SI."""
import math, numpy as np, fullstate as fs
d = fs.load('graphene'); tau = d['tau']; w = d['w']; Ns = d['N'] / 2
Dstar, mu_s, spec, frames = fs.canonical_full(d['hs'], w, tau, Ns)
P, TH, Q, LAM, E = [], [], [], [], []
for h, u, f in zip(d['hs'], d['U1'], d['f1']):
    occ = np.flatnonzero(f > 1e-12); m = int(occ[-1]) + 1 + 16
    p = u[:, :m]; th, v = np.linalg.eigh(p.conj().T @ h @ p); pr = p @ v
    q = u[:, m:]; lam, z = np.linalg.eigh(q.conj().T @ h @ q); qr = q @ z
    P.append(pr); TH.append(th); Q.append(qr); LAM.append(lam); E.append(qr.conj().T @ h @ pr)
mu_c, _ = fs.solve_mu(TH, w, tau, Ns)
mu0, _ = fs.solve_mu([np.r_[a, b] for a, b in zip(TH, LAM)], w, tau, Ns)
print('mu_c-mu0 %.2e  mu*-mu0 %.3e' % (mu_c - mu0, mu_s - mu0))
tot_excess = 0; per_k = []
for k, (pr, th, qr, lam, Ek, Ds) in enumerate(zip(P, TH, Q, LAM, E, Dstar)):
    Dc = (pr * fs.fermi((th - mu_c) / tau)) @ pr.conj().T
    A = Dc - Ds
    full2 = 2 * w[k] * np.linalg.norm(A) ** 2
    g = fs.fdd(lam - mu0, th - mu0, tau)
    dk2 = 2 * w[k] * 2 * np.sum(np.abs(g) ** 2 * np.abs(Ek) ** 2)
    per_k.append(full2 - dk2)
    # close pairs
    gaps = np.abs(lam[:, None] - th[None, :])
    j, i = np.unravel_index(np.argmin(gaps), gaps.shape)
    # near-Fermi occupation errors in the Ritz basis
    Bpp = pr.conj().T @ A @ pr
    near = np.flatnonzero(np.abs(th - mu0) < 0.05)
    dth = -np.sum(np.abs(Ek) ** 2 / (lam[:, None] - th[None, :]), axis=0)
    print('k%d excess %.3e | closest pair gap/tau %.3f at theta-mu0 %.3f coupling %.1e | near-Fermi: th-mu0 %s diag err %s dtheta %s' % (
        k, full2 - dk2, gaps[j, i] / tau, th[i] - mu0, abs(Ek[j, i]),
        np.round(th[near] - mu0, 4), np.round(np.real(np.diag(Bpp))[near], 4), np.round(dth[near], 5)))
per_k = np.array(per_k); print('total excess %.3e; share of k4 %.3f' % (per_k.sum(), per_k[4] / per_k.sum()))
