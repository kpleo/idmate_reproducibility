"""Export the k-resolved block matrices used in Fig. 1(a,b): the saved Si Hamiltonian H2
written in the basis [Ritz vectors of P | eigenvectors of QHQ], and the full-state error
D_c - D* in the same basis.  Writes derived/si_block_k{K}.npz (small, no private data)."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fullstate as fs

def blocks(mat='si', k=0, extra=2):
    d = fs.load(mat)
    tau, Nspin = d['tau'], d['N'] / 2
    Dstar, mu, spec, frames = fs.canonical_full(d['hs'], d['w'], tau, Nspin)
    TH, PR, QR = [], [], []
    for h, u, f in zip(d['hs'], d['U1'], d['f1']):
        occ = np.flatnonzero(f > 1e-12); m = int(occ[-1]) + 1 + extra
        p = u[:, :m]; th, v = np.linalg.eigh(p.conj().T @ h @ p)
        q = u[:, m:]; lam, z = np.linalg.eigh(q.conj().T @ h @ q)
        TH.append(th); PR.append(p @ v); QR.append(q @ z)
    mu_c, _ = fs.solve_mu(TH, d['w'], tau, Nspin)
    h = d['hs'][k]; B = np.hstack([PR[k], QR[k]])
    Hb = B.conj().T @ h @ B
    Dc = (PR[k] * fs.fermi((TH[k] - mu_c) / tau)) @ PR[k].conj().T
    dD = B.conj().T @ (Dc - Dstar[k]) @ B
    return dict(H=Hb, dD=dD, m=PR[k].shape[1], n=h.shape[0], mu=mu_c, theta=TH[k], tau=tau)

if __name__ == '__main__':
    for k in (0,):
        b = blocks('si', k)
        np.savez_compressed(f'derived/si_block_k{k}.npz', H_abs=np.abs(b['H']), dD_abs=np.abs(b['dD']),
                            m=b['m'], n=b['n'], mu=b['mu'], theta=b['theta'], tau=b['tau'])
        P = slice(0, b['m']); Q = slice(b['m'], b['n'])
        print('k', k, 'm', b['m'], 'n', b['n'],
              '|dD| blocks: PP %.2e QQ %.2e PQ %.2e' % (np.linalg.norm(b['dD'][P, P]), np.linalg.norm(b['dD'][Q, Q]),
                                                         np.linalg.norm(b['dD'][P, Q])),
              '|H| offdiag %.3e' % np.linalg.norm(b['H'][Q, P]))
