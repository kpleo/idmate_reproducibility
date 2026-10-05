#!/usr/bin/env python3
"""In-loop candidates of the 2026-10-05 rerun with saved Hamiltonians.

For every candidate presented to the screen in the SCF tests (reference-scoring arm), the instrumented harness (IDMate commit
559dfc1 plus the test-only hdump patch) wrote the screen's full-cutoff H_k, the trial orbitals and occupations, the
exact band frame the auditor measures against, the weights, tau and N_e. This script

1. reproduces the recorded auditor full-state distance from the dumped frames;
2. recomputes the reference state D* from H_k and checks it against the dumped exact frame;
3. for the Rayleigh--Ritz families (low cutoff, oracle subspace) evaluates the exact block identity, the first- and
   second-order decomposition, the indicator I1 and the rigorous/computable bounds with the same estimators as the
   saved-Hamiltonian study (fullstate.analyse_blocks);
4. writes inloop_results.json.

Input layout: <RUNS>/<leg>/hdump/i<it>_<kind>[_rN]/{meta.json, kXX_{H,P,X}.bin[.xz|.gz]} with RUNS taken from
IDMATE_INLOOP_RUNS (default: ../inloop_rerun_2026-10-05/runs relative to this folder)."""
import gzip, json, lzma, math, os, sys, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import fullstate as fs

RUNS = os.environ.get('IDMATE_INLOOP_RUNS') or os.path.join(HERE, '..', 'inloop_rerun_2026-10-05', 'runs')
RITZ_KINDS = ('low_cutoff_continuation', 'converged_subspace_reuse')


def _bin(d, v):
    p = os.path.join(d, v['file'])
    if os.path.exists(p):
        raw = open(p, 'rb').read()
    elif os.path.exists(p + '.xz'):
        raw = lzma.open(p + '.xz', 'rb').read()
    else:
        raw = gzip.open(p + '.gz', 'rb').read()
    a = np.frombuffer(raw, dtype='<c16')
    assert a.size == int(v['rows']) * int(v['cols']), (p, a.size, v)
    return a.reshape((int(v['rows']), int(v['cols'])), order='F')


def load_candidate(d):
    meta = json.load(open(os.path.join(d, 'meta.json')))
    ks = []
    for p in meta['kpoints']:
        ks.append(dict(H=_bin(d, p['H']), P=_bin(d, p['proposal']['vectors']),
                       fP=np.asarray(p['proposal']['occupations'], float),
                       eP=np.asarray(p['proposal']['eigenvalues_ha'], float),
                       win=np.asarray(p['proposal']['window_spectrum_ha'], float),
                       X=_bin(d, p['exact']['vectors']), fX=np.asarray(p['exact']['occupations'], float),
                       eX=np.asarray(p['exact']['eigenvalues_ha'], float), w=float(p['weight'])))
    return meta, ks


def code_density(V, f):
    """The auditor's complex_spectral_density: per-spin density, occupations below 1e-12 dropped."""
    g = np.where(f < 1e-12, 0.0, f) / 2.0
    return (V * g) @ V.conj().T


def analyse_candidate(meta, ks):
    tau = float(meta['temperature_ha'])
    Ns = float(meta['electron_count']) / 2.0
    w = np.array([k['w'] for k in ks])
    hs = [k['H'] for k in ks]
    scr, aud = meta['screen'], meta['auditor']
    out = dict(iteration=int(meta['iteration']), kind=meta['kind'], decision=scr['decision'],
               abstain_reason=scr['abstain_reason'], window_bound=scr['distance_bound'],
               window_bands=scr['window_bands'], d_window=aud['true_distance_window'],
               d_full_recorded=aud['true_distance_full'], tau=tau, N=2 * Ns, nk=len(ks),
               n=[int(h.shape[0]) for h in hs], m=[int(k['P'].shape[1]) for k in ks])
    # 1. the recorded auditor distance, recomputed from the dumped frames
    out['d_full_reproduced'] = math.sqrt(sum(k['w'] * 2.0 * np.linalg.norm(code_density(k['P'], k['fP'])
                                                                             - code_density(k['X'], k['fX'])) ** 2
                                             for k in ks))
    # 2. reference state from H_k (dense eigh, common mu) against the dumped exact frame
    ref = fs.canonical_full(hs, w, tau, Ns)
    Dstar, mu_star = ref[0], ref[1]
    out['ref_vs_exact_frame'] = math.sqrt(sum(wk * 2.0 * np.linalg.norm(code_density(k['X'], k['fX']) - Ds) ** 2
                                              for wk, k, Ds in zip(w, ks, Dstar)))
    out['mu_star'] = float(mu_star)
    out['mu_exact_recorded'] = float(meta['exact_chemical_potential_ha'])
    out['herm_max'] = float(max(np.abs(h - h.conj().T).max() for h in hs))
    if meta['kind'] not in RITZ_KINDS:
        # coarse-k upsample: the exact frame with transferred occupations (no complement); its whole error is the
        # occupation transfer, i.e. the retained block.
        out['ritz'] = False
        return out
    # 3. Rayleigh--Ritz structure of the trial subspace
    out['ritz'] = True
    orth, offdiag, eigdev = 0.0, 0.0, 0.0
    P, TH, Qs, LAM, E = [], [], [], [], []
    for k in ks:
        p, h = k['P'], k['H']
        m = p.shape[1]
        orth = max(orth, float(np.abs(p.conj().T @ p - np.eye(m)).max()))
        hp = p.conj().T @ h @ p
        hp = 0.5 * (hp + hp.conj().T)
        offdiag = max(offdiag, float(np.linalg.norm(hp - np.diag(np.diag(hp))) / max(np.linalg.norm(hp), 1e-300)))
        th, v = np.linalg.eigh(hp)
        eigdev = max(eigdev, float(np.abs(np.sort(k['eP']) - th).max()))
        pr = p @ v
        qfull, _ = np.linalg.qr(p, mode='complete')
        q = qfull[:, m:]
        hq = q.conj().T @ h @ q
        hq = 0.5 * (hq + hq.conj().T)
        lam, z = np.linalg.eigh(hq)
        qr = q @ z
        P.append(pr); TH.append(th); Qs.append(qr); LAM.append(lam); E.append(qr.conj().T @ h @ pr)
    out.update(orthonormality=orth, ritz_offdiag_rel=offdiag, ritz_vs_recorded_eigs=eigdev)
    r = fs.analyse_blocks(hs, w, tau, Ns, P, TH, Qs, LAM, E, Dstar, mu_star)
    # the analysed candidate (Ritz vectors, common mu_c from the Ritz values) must be the recorded one
    out['candidate_reconstruction'] = abs(r['d_full'] - out['d_full_reproduced'])
    out['mu_c_recomputed'] = float(r['mu_c'])
    out['mu_c_recorded'] = float(meta['proposal_chemical_potential_ha'])
    keep = ('d_full', 'eP', 'eQ', 'c_w', 'bound_cpl', 'rigorous', 'rigorous_computable', 'd_dk', 'd_cheap', 'eQ_est',
            'eP_est', 'eP2_est', 'indicator', 'indicator2', 'C_hat', 'blk_P', 'blk_Q', 'blk_off', 'eP2', 'eQ2', 'pyth2',
            'd_cpl_exact', 'mu_c', 'mu0', 'gapQ', 'topP')
    for key in keep:
        v = r[key]
        out[key] = float(v) if v is not None else None
    out['pyth'] = math.sqrt(r['eP'] ** 2 + r['eQ'] ** 2 + r['d_dk'] ** 2)
    out['identity_residual'] = abs(math.sqrt(r['blk_P'] ** 2 + r['blk_Q'] ** 2 + r['blk_off'] ** 2) - r['d_full'])
    return out


def leg_label(leg):
    """Map a harness leg name to the paper's run id (material, family, insertion iteration)."""
    mat = 'Si' if 'si_diamond' in leg else ('graphene' if 'graphene' in leg else 'Al')
    if leg.startswith('s8_'):
        return mat, 'multiple'
    code = leg.rsplit('_', 1)[1]
    fam = {'a': 'low_cutoff', 'b': 'target_subspace', 'c': 'coarse_k'}[code[0]]
    return mat, f'single_{fam}_{code[1:]}'


def main():
    rows = []
    for leg_dir in sorted(glob.glob(os.path.join(RUNS, '*'))):
        leg = os.path.basename(leg_dir)
        for cdir in sorted(glob.glob(os.path.join(leg_dir, 'hdump', 'i*'))):
            meta, ks = load_candidate(cdir)
            r = analyse_candidate(meta, ks)
            mat, run = leg_label(leg)
            r.update(leg=leg, material=mat, run_id=f'{mat}_{run}', dump=os.path.basename(cdir))
            rows.append(r)
            print('%-46s i%02d %-25s %-7s d=%.6e rec=%.1e%s' % (
                leg, r['iteration'], r['kind'], r['decision'], r['d_full_recorded'],
                abs(r['d_full_reproduced'] - r['d_full_recorded']) / max(r['d_full_recorded'], 1e-300),
                '' if not r['ritz'] else '  I1/d=%.3f pyth dev=%.2e bound/d=%.1f' % (
                    r['indicator'] / r['d_full'], abs(r['pyth'] - r['d_full']) / r['d_full'], r['rigorous'] / r['d_full'])),
                flush=True)
    json.dump(rows, open(os.path.join(HERE, 'inloop_results.json'), 'w'), indent=1)
    print('candidates:', len(rows))


if __name__ == '__main__':
    main()
