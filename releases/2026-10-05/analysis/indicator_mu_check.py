"""Compare the indicator with divided differences at mu_c (used) and at mu_0 (needs complement spectrum)."""
import math, json, numpy as np, fullstate as fs
def cheap(data, extra, tau, at):
    hs, U1, f1, w = data['hs'], data['U1'], data['f1'], data['w']; Ns = data['N'] / 2
    TH, LAM, EB = [], [], []
    for h, u, f in zip(hs, U1, f1):
        occ = np.flatnonzero(f > 1e-12); m = min(int(occ[-1]) + 1 + extra, u.shape[1])
        p = u[:, :m]; th, v = np.linalg.eigh(p.conj().T @ h @ p); pr = p @ v
        q = u[:, m:]; lam, z = np.linalg.eigh(q.conj().T @ h @ q); qr = q @ z
        TH.append(th); LAM.append(lam); EB.append(qr.conj().T @ h @ pr)
    mu_c, _ = fs.solve_mu(TH, w, tau, Ns)
    mu0, _ = fs.solve_mu([np.r_[a, b] for a, b in zip(TH, LAM)], w, tau, Ns)
    mu = mu_c if at == 'c' else mu0
    ch2 = 0.0
    for wk, th, lam, E in zip(w, TH, LAM, EB):
        rho = np.linalg.norm(E, axis=0); lmin = lam.min()
        yg = (lmin - mu) + np.concatenate([np.linspace(0, 60 * tau, 601), np.geomspace(60 * tau + 1e-3, 1e3, 200)])
        g = np.max(np.abs(fs.fdd(th - mu, yg, tau)), axis=1)
        ch2 += wk * 2.0 * np.sum(g ** 2 * rho ** 2)
    return math.sqrt(2 * ch2)
o = json.load(open('decomposition_results.json'))
worst = 0; worst_case = None
for mat in ['si', 'graphene', 'al']:
    data = fs.load(mat)
    cases = [(r['extra'], data['tau'], r) for r in o[mat]['extra_sweep']] + [(2, r['tau'], r) for r in o[mat]['tau_sweep'] if not r['duplicate_of_buffer_sweep']]
    for b, t, r in cases:
        c0 = cheap(data, b, t, '0')
        I0 = math.sqrt(r['eP_est'] ** 2 + r['eQ_est'] ** 2 + c0 ** 2)
        rel = abs(I0 - r['indicator']) / r['indicator']
        if rel > worst: worst, worst_case = rel, (mat, b, t, r['indicator'] / r['d_full'], I0 / r['d_full'])
print('max relative change of I1 (mu_0 vs mu_c): %.4f' % worst, worst_case)
