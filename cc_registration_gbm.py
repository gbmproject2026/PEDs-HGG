"""
CC involvement pipeline for adult BraTS-GLI (GBM cases).
Reads the GBM-like subset (ET > 100) from cc_results_gbm/gbm_cases.txt,
which scan_gli_et.py generates. Run scan_gli_et.py first.
"""
import ants
import nibabel as nib
import numpy as np
import os, json, time
from multiprocessing import Pool, cpu_count

BASE        = os.path.dirname(os.path.abspath(__file__))
DATA        = os.path.join(BASE, "data")

_GLI_ROOT   = os.path.join(DATA, "BraTS-GLI")
def _resolve_gli_dir(root):
    if not os.path.isdir(root):
        return root
    subdirs = [os.path.join(root, e) for e in os.listdir(root)
               if os.path.isdir(os.path.join(root, e))]
    if any(os.path.basename(d).startswith("BraTS-GLI-") for d in subdirs):
        return root
    for d in subdirs:
        for sub in [d] + [os.path.join(d, x) for x in os.listdir(d)
                          if os.path.isdir(os.path.join(d, x))]:
            if any(f.startswith("BraTS-GLI-") for f in os.listdir(sub)):
                return sub
    return root

GLI_DIR     = _resolve_gli_dir(_GLI_ROOT)
ATLAS_DIR   = os.path.join(BASE, "atlas")
RESULTS_DIR = os.path.join(BASE, "cc_results_gbm")
CASES_FILE  = os.path.join(RESULTS_DIR, "gbm_cases.txt")
MNI_PATH    = ants.get_ants_data('mni')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(os.path.join(RESULTS_DIR, "registrations"), exist_ok=True)

CC_MASKS = {
    "whole":    ants.image_read(f"{ATLAS_DIR}/CC_mask.nii.gz").numpy() > 0,
    "genu":     ants.image_read(f"{ATLAS_DIR}/CC_genu.nii.gz").numpy() > 0,
    "body":     ants.image_read(f"{ATLAS_DIR}/CC_body.nii.gz").numpy() > 0,
    "splenium": ants.image_read(f"{ATLAS_DIR}/CC_splenium.nii.gz").numpy() > 0,
}
CC_VOL = {k: int(v.sum()) for k, v in CC_MASKS.items()}
MNI    = ants.image_read(MNI_PATH)


def register_case(subject_id):
    out_json = os.path.join(RESULTS_DIR, f"{subject_id}_cc.json")
    if os.path.exists(out_json):
        with open(out_json) as f:
            return json.load(f)

    t1c_path = os.path.join(GLI_DIR, subject_id, f"{subject_id}-t1c.nii.gz")
    seg_path = os.path.join(GLI_DIR, subject_id, f"{subject_id}-seg.nii.gz")
    if not os.path.exists(t1c_path) or not os.path.exists(seg_path):
        return {"subject": subject_id, "error": "missing files"}

    try:
        t0 = time.time()

        t1c = ants.image_read(t1c_path)
        seg = ants.image_read(seg_path)

        tumor_bin = ants.threshold_image(seg, low_thresh=0.5, high_thresh=4.5, inval=1, outval=0)
        tumor_dil = ants.morphology(tumor_bin, operation="dilate", radius=5, mtype="binary", shape="ball")
        cost_mask = ants.threshold_image(tumor_dil, low_thresh=0.5, inval=0, outval=1)

        reg = ants.registration(
            fixed=MNI, moving=t1c,
            type_of_transform="SyNRA",
            mask=cost_mask,
            reg_iterations=(40, 20, 10),
            random_seed=42,
            verbose=False,
        )

        seg_mni = ants.apply_transforms(
            fixed=MNI, moving=seg,
            transformlist=reg["fwdtransforms"],
            interpolator="nearestNeighbor",
        )
        seg_np = seg_mni.numpy()

        reg_out = os.path.join(RESULTS_DIR, "registrations", f"{subject_id}_seg_MNI.nii.gz")
        ants.image_write(seg_mni, reg_out)

        wt = seg_np > 0
        et = seg_np == 3
        tc = (seg_np == 1) | (seg_np == 3)
        mid = seg_np.shape[0] // 2
        cc_w = CC_MASKS["whole"]

        r = {
            "subject": subject_id,
            "wt_vol": int(wt.sum()),
            "et_vol": int(et.sum()),
            "elapsed_s": round(time.time() - t0, 1),
            "error": None,
        }
        for mname, marr in CC_MASKS.items():
            for lname, larr in [("wt", wt), ("et", et), ("tc", tc)]:
                ov = int((larr & marr).sum())
                r[f"cc_{mname}_{lname}_overlap"] = ov
                r[f"cc_{mname}_{lname}_involved"] = ov > 0
                r[f"cc_{mname}_{lname}_frac"]    = round(ov / CC_VOL[mname], 4)

        r["butterfly_et"] = bool((et & cc_w)[:mid].sum() > 5 and (et & cc_w)[mid:].sum() > 5)
        r["butterfly_wt"] = bool((wt & cc_w)[:mid].sum() > 20 and (wt & cc_w)[mid:].sum() > 20)

        with open(out_json, "w") as f:
            json.dump(r, f)

        print(f"  {subject_id}: CC_wt={r['cc_whole_wt_involved']} CC_et={r['cc_whole_et_involved']} butterfly_ET={r['butterfly_et']}  [{r['elapsed_s']}s]")
        return r

    except Exception as e:
        err = {"subject": subject_id, "error": str(e)}
        with open(out_json, "w") as f:
            json.dump(err, f)
        print(f"  {subject_id}: ERROR {e}")
        return err


if __name__ == "__main__":
    if not os.path.exists(CASES_FILE):
        raise SystemExit(f"Missing {CASES_FILE}. Run scan_gli_et.py first.")

    with open(CASES_FILE) as f:
        gbm_ids = [l.strip() for l in f if l.strip()]

    print(f"GLI_DIR = {GLI_DIR}")
    print(f"Running CC pipeline on {len(gbm_ids)} GBM-likely cases")
    print(f"CPUs available: {cpu_count()}, using {max(1, cpu_count()-2)} workers\n")

    n_workers = max(1, cpu_count() - 2)
    with Pool(n_workers) as pool:
        all_results = pool.map(register_case, gbm_ids)

    valid = [r for r in all_results if not r.get("error")]
    cc_wt = sum(1 for r in valid if r.get("cc_whole_wt_involved"))
    cc_et = sum(1 for r in valid if r.get("cc_whole_et_involved"))
    bf_et = sum(1 for r in valid if r.get("butterfly_et"))
    n = max(1, len(valid))

    print(f"\n{'='*55}")
    print(f"ADULT GBM CC RESULTS — {len(valid)} cases")
    print(f"{'='*55}")
    print(f"CC involvement — whole tumor:     {cc_wt}/{len(valid)} ({100*cc_wt/n:.1f}%)")
    print(f"CC involvement — ET only:         {cc_et}/{len(valid)} ({100*cc_et/n:.1f}%)")
    print(f"Butterfly pattern (ET bilateral): {bf_et}/{len(valid)} ({100*bf_et/n:.1f}%)")

    import pandas as pd
    pd.DataFrame(valid).to_csv(f"{RESULTS_DIR}/cc_results_gbm.csv", index=False)
    print(f"\nSaved -> {RESULTS_DIR}/cc_results_gbm.csv")
