"""
Filter the adult BraTS-GLI cohort to GBM-like cases (ET > 100 voxels) and write
the surviving subject IDs to cc_results_gbm/gbm_cases.txt.
Run before cc_registration_gbm.py.
"""
import nibabel as nib
import numpy as np
import os
from multiprocessing import Pool, cpu_count

BASE    = os.path.dirname(os.path.abspath(__file__))
DATA    = os.path.join(BASE, "data")

_GLI_ROOT = os.path.join(DATA, "BraTS-GLI")
def _resolve_gli_dir(root):
    if not os.path.isdir(root):
        return root
    entries = [os.path.join(root, e) for e in os.listdir(root)]
    subdirs = [e for e in entries if os.path.isdir(e)]
    if any(os.path.basename(d).startswith("BraTS-GLI-") for d in subdirs):
        return root
    for d in subdirs:
        for sub in [d] + [os.path.join(d, x) for x in os.listdir(d)
                          if os.path.isdir(os.path.join(d, x))]:
            if any(f.startswith("BraTS-GLI-") for f in os.listdir(sub)):
                return sub
    return root

GLI_DIR   = _resolve_gli_dir(_GLI_ROOT)
OUT_DIR   = os.path.join(BASE, "cc_results_gbm")
OUT_FILE  = os.path.join(OUT_DIR, "gbm_cases.txt")
ET_THRESH = 100

os.makedirs(OUT_DIR, exist_ok=True)


def get_volumes(sid):
    seg_path = os.path.join(GLI_DIR, sid, f"{sid}-seg.nii.gz")
    if not os.path.exists(seg_path):
        return sid, 0, 0
    seg = nib.load(seg_path).get_fdata()
    return sid, int((seg > 0).sum()), int((seg == 3).sum())


if __name__ == "__main__":
    cases = sorted(d for d in os.listdir(GLI_DIR)
                   if os.path.isdir(os.path.join(GLI_DIR, d)))
    print(f"GLI_DIR = {GLI_DIR}")
    print(f"Scanning {len(cases)} GLI cases for ET volume...")

    with Pool(max(1, cpu_count() - 2)) as p:
        results = p.map(get_volumes, cases)

    total    = len(results)
    et_zero  = sum(1 for _, _, et in results if et == 0)
    et_gt    = sum(1 for _, _, et in results if et > ET_THRESH)

    print(f"Total: {total}")
    print(f"ET = 0:        {et_zero} ({100*et_zero/total:.1f}%) — likely LGG")
    print(f"ET > {ET_THRESH}:      {et_gt} ({100*et_gt/total:.1f}%) — GBM proxy")

    gbm_cases = [sid for sid, _, et in results if et > ET_THRESH]
    with open(OUT_FILE, "w") as f:
        f.write("\n".join(gbm_cases))
    print(f"\nGBM-likely list saved: {len(gbm_cases)} cases -> {OUT_FILE}")
    print("SCAN_COMPLETE")
