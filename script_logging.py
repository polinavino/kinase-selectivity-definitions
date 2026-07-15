"""Per-script bookkeeping shared by the analysis scripts.

Add this as the first line of an analysis script::

    import script_logging; script_logging.capture(__file__)

``capture`` does two things on every run, so the repository stays in sync with
the code without any manual steps:

1. Mirrors the script's stdout to ``outputs/<name>.txt`` -- a committed text
   file holding its printed results, so the numbers quoted in the paper can be
   checked without rerunning anything. These files change only when a script
   (or its data) changes.

2. Propagates any manuscript figure the script (re)generates into the two
   locations the paper builds from: ``paper/figures/<name>.png`` (the LaTeX
   include) and ``paper/submission_figures/Figure<N>.png`` (the journal upload).
   Only figures actually rewritten during the run are copied, so rerunning one
   script never disturbs another's figure.
"""
import os
import sys
import atexit
import shutil

# Manuscript figure basename (as written next to the scripts) -> journal figure
# number. Order follows the figures' appearance in the compiled paper.
FIGURE_MAP = {
    "metric_fidelity.png": "Figure1.png",
    "cross_dataset_correlations.png": "Figure2.png",
    "instability_by_family.png": "Figure3.png",
    "binding_profiles.png": "Figure4.png",
    "panel_size_stability.png": "Figure5.png",
    "candidate_panel_convergence.png": "Figure6.png",
}


def _mtime(path):
    try:
        return os.path.getmtime(path)
    except OSError:
        return None


def _mirror_figures(repo, before):
    fig_dir = os.path.join(repo, "paper", "figures")
    sub_dir = os.path.join(repo, "paper", "submission_figures")
    if not os.path.isdir(os.path.join(repo, "paper")):
        return
    for name, journal_name in FIGURE_MAP.items():
        src = os.path.join(repo, name)
        now = _mtime(src)
        # Only copy figures this run actually (re)wrote.
        if now is None or now == before.get(name):
            continue
        os.makedirs(fig_dir, exist_ok=True)
        os.makedirs(sub_dir, exist_ok=True)
        shutil.copyfile(src, os.path.join(fig_dir, name))
        shutil.copyfile(src, os.path.join(sub_dir, journal_name))


def capture(script_path):
    repo = os.path.dirname(os.path.abspath(script_path))

    # (2) snapshot figure mtimes before the script runs
    before = {name: _mtime(os.path.join(repo, name)) for name in FIGURE_MAP}

    # (1) tee stdout to outputs/<name>.txt
    name = os.path.splitext(os.path.basename(script_path))[0]
    outdir = os.path.join(repo, "outputs")
    os.makedirs(outdir, exist_ok=True)
    fh = open(os.path.join(outdir, name + ".txt"), "w")
    real = sys.stdout

    class _Tee:
        def write(self, data):
            real.write(data)
            if not fh.closed:
                fh.write(data)
            return len(data)

        def flush(self):
            real.flush()
            if not fh.closed:
                fh.flush()

    sys.stdout = _Tee()

    def _finish():
        # Restore the real stream first so output emitted during interpreter
        # shutdown does not hit the closed log file.
        sys.stdout = real
        fh.flush()
        fh.close()
        _mirror_figures(repo, before)

    atexit.register(_finish)
