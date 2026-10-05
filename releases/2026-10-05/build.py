"""Rebuild every figure and Supplemental table of this release from the derived inputs, then run the numerical
tests and the input-identity check. Outputs are written under build/ (not distributed)."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
FIGURES = ['fig1_blocks.py', 'fig2_inloop.py', 'fig3_decomposition.py', 'figS1_blocks.py', 'figS2_silicon.py',
           'figS3S4_auxiliary.py']


def run(cmd, cwd):
    done = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return done.returncode, (done.stdout + done.stderr).strip().splitlines()[-1:] if (done.stdout + done.stderr).strip() else []


def main():
    from matplotlib import font_manager
    if 'Liberation Sans' not in {f.name for f in font_manager.fontManager.ttflist}:
        sys.exit('Liberation Sans is required for the figures (text and math); install it and rerun.')
    summary = {'figures': {}, 'tables': None, 'tests': None, 'inputs': None}
    for script in FIGURES:
        rc, tail = run([sys.executable, script], os.path.join(HERE, 'figures'))
        summary['figures'][script] = 'ok' if rc == 0 else f'failed: {tail}'
    rc, _ = run([sys.executable, os.path.join('tables', 'make_tables.py')], HERE)
    summary['tables'] = 'rows identical to the manuscript' if rc == 0 else 'differences, see build/tables'
    rc, tail = run([sys.executable, '-m', 'unittest', 'discover', '-s', 'tests'], HERE)
    summary['tests'] = 'passed' if rc == 0 else f'failed: {tail}'
    rc, tail = run([sys.executable, 'verify.py'], HERE)
    summary['inputs'] = 'manifest verified' if rc == 0 else f'failed: {tail}'
    os.makedirs(os.path.join(HERE, 'build'), exist_ok=True)
    with open(os.path.join(HERE, 'build', 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))
    ok = all(v == 'ok' for v in summary['figures'].values()) and summary['tests'] == 'passed' \
        and summary['inputs'] == 'manifest verified' and summary['tables'].startswith('rows identical')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
