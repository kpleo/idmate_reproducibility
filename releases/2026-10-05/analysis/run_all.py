import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fullstate as fs

def candidate_matrices(data, extra, tau):
    """Return Dc (Ritz candidate), D* (full), D0* (block-diagonal canonical) lists for density work."""
    hs, U1, f1, w = data['hs'], data['U1'], data['f1'], data['w']
    Nspin = data['N'] / 2
    Dstar, mu_star, spec, frames = fs.canonical_full(hs, w, tau, Nspin)
    P, TH = [], []
    for h, u, f in zip(hs, U1, f1):
        occ = np.flatnonzero(f > 1e-12)
        m = min(int(occ[-1]) + 1 + extra, u.shape[1])
        p = u[:, :m]; th, v = np.linalg.eigh(p.conj().T @ h @ p); P.append(p @ v); TH.append(th)
    mu_c, _ = fs.solve_mu(TH, w, tau, Nspin)
    Dc = [(pr * fs.fermi((th - mu_c) / tau)).dot(pr.conj().T) for pr, th in zip(P, TH)]
    return Dc, Dstar

out = {}
for mat in ['si', 'graphene', 'al']:
    data = fs.load(mat)
    phases = [fs.grid_phases(data, k) for k in range(len(data['hs']))]
    nmin = min(h.shape[0] for h in data['hs'])
    occmax = max(int(np.flatnonzero(f > 1e-12)[-1]) + 1 for f in data['f1'])
    tau0 = data['tau']
    ref0 = fs.canonical_full(data['hs'], data['w'], tau0, data['N'] / 2)
    rho_star0 = fs.density(data, ref0[0], phases)
    rows = []
    for e in range(0, nmin - occmax):          # every buffer size that leaves a nonempty complement
        r = fs.analyse(data, extra=e, ref=ref0)
        Dc, _ = candidate_matrices(data, e, tau0)
        r['eps_rho'] = fs.eps_rho(data, fs.density(data, Dc, phases), rho_star0)
        r['pyth'] = math.sqrt(r['eP'] ** 2 + r['eQ'] ** 2 + r['d_dk'] ** 2)
        r['duplicate_of_buffer_sweep'] = False
        rows.append(r)
    trows = []
    for t in np.geomspace(0.002, 0.2, 15):
        t = float(t)
        r = fs.analyse(data, extra=2, tau=t)
        Dc, Dstar = candidate_matrices(data, 2, t)
        r['eps_rho'] = fs.eps_rho(data, fs.density(data, Dc, phases), fs.density(data, Dstar, phases))
        r['pyth'] = math.sqrt(r['eP'] ** 2 + r['eQ'] ** 2 + r['d_dk'] ** 2)
        r['duplicate_of_buffer_sweep'] = abs(t - tau0) < 1e-12 * tau0
        trows.append(r)
    out[mat] = {'tau': tau0, 'N': data['N'], 'extra_sweep': rows, 'tau_sweep': trows,
                'nk': len(data['hs']), 'w_min': float(np.min(data['w']))}
    for r in rows[::5]:
        print(mat, 'e=%3d d=%.4e pyth=%.4e rel=%.2e cheap/d=%.2f bound/d=%.1f eps_rho=%.4f' % (
            r['extra'], r['d_full'], r['pyth'], abs(r['pyth'] - r['d_full']) / r['d_full'], r['d_cheap'] / r['d_full'], r['rigorous'] / r['d_full'], r['eps_rho']))
    for r in trows:
        print(mat, 'tau=%.4f d=%.4e pyth=%.4e rel=%.2e eQ=%.2e dk=%.2e bound/d=%.2f eps_rho=%.4f' % (
            r['tau'], r['d_full'], r['pyth'], abs(r['pyth'] - r['d_full']) / r['d_full'], r['eQ'], r['d_dk'], r['rigorous'] / r['d_full'], r['eps_rho']))
json.dump(out, open(os.path.join(HERE, 'decomposition_results.json'), 'w'), indent=1, default=float)

# rigorous-bound check over everything
viol = 0; n = 0; worst = 0
for mat in out:
    for r in out[mat]['extra_sweep'] + out[mat]['tau_sweep']:
        n += 1; worst = max(worst, r['d_full'] / r['rigorous'])
        if r['d_full'] > r['rigorous'] * (1 + 1e-12): viol += 1
print('rigorous bound checks:', n, 'violations', viol, 'max d/bound', worst)
print('distinct candidates:', sum(1 for m in out for r in out[m]['extra_sweep'] + out[m]['tau_sweep'] if not r['duplicate_of_buffer_sweep']))
