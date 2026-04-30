"""Regenerate tutorial Links/Oracle artifacts so they are loadable under pandas >=2.0.

Background
----------
The original tutorial artifacts (`links_louvain_v20220406.celloracle.links`,
`Paul_etal_v20220406.celloracle.oracle`) were saved with pandas 1.5 and
anndata <=0.10. Their pickle streams reference internal classes that were
removed in pandas 2.0 (`pandas.core.indexes.numeric`) and changed in
anndata 0.12 (`AnnData.__setstate__` no longer accepts the old key set).

This regeneration runs in two phases across two conda envs:

  Phase 1: in a pandas <2.0 / anndata <=0.10 env (e.g. test310)
    - Load Links; strip Int64Index from filtered_links; dump the new file.
    - Load Oracle; export oracle.adata as standalone h5ad and the remaining
      attributes (which contain no pandas internal-class references) as a
      portable pickle.

  Phase 2: in a pandas >=2.0 / anndata >=0.12 env (e.g. test312)
    - Read the h5ad with the new anndata, unpickle the remaining attributes,
      reattach, and re-dump as the new .celloracle.oracle artifact.

Run as:
    conda activate test310 && python regen_tutorial_artifacts.py phase1
    conda activate test312 && python regen_tutorial_artifacts.py phase2
"""

import os
import pickle
import sys

import pandas as pd

import celloracle as co

# Output filenames must match what `load_data.py` looks up under version "0.20.0".
LINKS_OUT_NAME = "links_louvain_v20260430.celloracle.links"
ORACLE_OUT_NAME = "Paul_etal_v20260430.celloracle.oracle"

OUT_DIR = os.path.join(
    os.path.dirname(os.path.abspath(co.__file__)),
    "data",
    "tutorial_data",
)
os.makedirs(OUT_DIR, exist_ok=True)
LINKS_OUT = os.path.join(OUT_DIR, LINKS_OUT_NAME)
ORACLE_OUT = os.path.join(OUT_DIR, ORACLE_OUT_NAME)

# Intermediate files written by phase 1 and consumed by phase 2.
ADATA_TMP = os.path.join(OUT_DIR, "_oracle_adata_v20260430.h5ad")
STATE_TMP = os.path.join(OUT_DIR, "_oracle_state_v20260430.pkl")


def _is_pandas1():
    return int(pd.__version__.split(".")[0]) < 2


def phase1():
    assert _is_pandas1(), (
        f"phase1 must run with pandas <2.0 (got {pd.__version__})"
    )
    from celloracle.utility.hdf5_processing import dump_hdf5

    print("[phase1] Regenerating Links artifact...")
    links = co.data.load_tutorial_links_object()
    for k, df in list(links.filtered_links.items()):
        # Int64Index from pandas 1.x -> RangeIndex (no semantic loss; the
        # original integer values were just original-row pointers from links_dict).
        links.filtered_links[k] = df.reset_index(drop=True)
    print(f"[phase1]   writing {LINKS_OUT}")
    dump_hdf5(links, LINKS_OUT)

    print("[phase1] Exporting Oracle adata + non-adata state...")
    oracle = co.data.load_tutorial_oracle_object()

    print(f"[phase1]   writing {ADATA_TMP}")
    oracle.adata.write_h5ad(ADATA_TMP, compression="gzip")

    state = {k: v for k, v in oracle.__dict__.items() if k != "adata"}
    print(f"[phase1]   writing {STATE_TMP}")
    with open(STATE_TMP, "wb") as f:
        pickle.dump(state, f, protocol=4)

    print("[phase1] done. Run phase2 in a pandas>=2.0 env next.")


def phase2():
    assert not _is_pandas1(), (
        f"phase2 must run with pandas>=2.0 (got {pd.__version__})"
    )
    import anndata
    from celloracle.trajectory.oracle_core import Oracle
    from celloracle.utility.hdf5_processing import dump_hdf5

    print("[phase2] Reading exported adata + state...")
    adata = anndata.read_h5ad(ADATA_TMP)
    with open(STATE_TMP, "rb") as f:
        state = pickle.load(f)

    oracle = object.__new__(Oracle)
    for k, v in state.items():
        setattr(oracle, k, v)
    oracle.adata = adata

    print(f"[phase2] writing {ORACLE_OUT}")
    dump_hdf5(oracle, ORACLE_OUT)

    print("[phase2] cleaning intermediate files")
    for p in (ADATA_TMP, STATE_TMP):
        try:
            os.remove(p)
        except FileNotFoundError:
            pass

    print("[phase2] done.")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in {"phase1", "phase2"}:
        print(__doc__)
        sys.exit(2)
    if sys.argv[1] == "phase1":
        phase1()
    else:
        phase2()
