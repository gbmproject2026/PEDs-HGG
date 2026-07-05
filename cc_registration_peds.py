"""
Corpus Callosum involvement analysis pipeline for pHGG (BraTS-PEDs).
Registers each case to ANTs MNI space, computes CC overlap metrics.
"""
import ants
import nibabel as nib
import numpy as np
import csv, os, json, time
from multiprocessing import Pool, cpu_count

BASE         = os.path.dirname(os.path.abspath(__file__))
DATA         = os.path.join(BASE, "data")

_IMG_ROOT    = os.path.join(DATA, "BraTS-PEDs-v1")
TRAIN_DIR    = os.path.join(_IMG_ROOT, "Training") \
               if os.path.isdir(os.path.join(_IMG_ROOT, "Training")) else _IMG_ROOT

ATLAS_DIR    = os.path.join(BASE, "atlas")
RESULTS_DIR  = os.path.join(BASE, "cc_results_peds")

METADATA_TSV = os.path.join(DATA, "BraTS-PEDs_metadata.tsv")
CBTN_TSV     = os.path.join(DATA, "histologies.tsv")
COHORT_COL   = "BraTS2025_cohort"

MNI_PATH     = ants.get_ants_data('mni')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(os.path.join(RESULTS_DIR, "registrations"), exist_ok=True)
CC_MASKS = {
    "whole":    ants.image_read(f"{ATLAS_DIR}/CC_mask.nii.gz").numpy() > 0,
    "genu":     ants.image_read(f"{ATLAS_DIR}/CC_genu.nii.gz").numpy() > 0,
    "body":     ants.image_read(f"{ATLAS_DIR}/CC_body.nii.gz").numpy() > 0,
    "splenium": ants.image_read(f"{ATLAS_DIR}/CC_splenium.nii.gz").numpy() > 0,
}
CC_VOL = {k: v.sum() for k, v in CC_MASKS.items()}
MNI    = ants.image_read(MNI_PATH)


def get_phgg_ids():
    brats_meta = {}
    with open(METADATA_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            brats_meta[row["BraTS-SubjectID"]] = row

    cbtn_map = {}
    with open(CBTN_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cid = row["cohort_participant_id"].strip()
            if cid and cid not in cbtn_map:
                cbtn_map[cid] = row["pathology_diagnosis"]

    ids = []
    for sid, row in brats_meta.items():
        if row[COHORT_COL] != "Training":
            continue
        src = row["Source"]
        is_phgg = (src == "DFCI-BCH-BWH-PEDs-HGG") or \
                  (src == "CBTN" and "High-grade" in cbtn_map.get(row["MappingID"], ""))
        if is_phgg:
            ids.append(sid)
    return sorted(ids)


def register_case(subject_id):
    out_json = os.path.join(RESULTS_DIR, f"{subject_id}_cc.json")
    if os.path.exists(out_json):
        with open(out_json) as f:
            return json.load(f)

    seg_path = os.path.join(TRAIN_DIR, subject_id, f"{subject_id}-seg.nii.gz")
    t1c_path = os.path.join(TRAIN_DIR, subject_id, f"{subject_id}-t1c.nii.gz")

    if not os.path.exists(seg_path) or not os.path.exists(t1c_path):
        return {"subject": subject_id, "error": "missing files"}

    try:
        t0 = time.time()

        t1c = ants.image_read(t1c_path)
        seg = ants.image_read(seg_path)

        # Cost function mask: dilate tumor, invert → exclude from registration
        tumor_bin = ants.threshold_image(seg, low_thresh=0.5, high_thresh=4.5, inval=1, outval=0)
        tumor_dil = ants.morphology(tumor_bin, operation="dilate", radius=5, mtype="binary", shape="ball")
        cost_mask = ants.threshold_image(tumor_dil, low_thresh=0.5, inval=0, outval=1)  # invert

        reg = ants.registration(
            fixed=MNI,
            moving=t1c,
            type_of_transform="SyNRA",
            mask=cost_mask,
            reg_iterations=(40, 20, 10),
            random_seed=42,
            verbose=False,
        )

        # Warp segmentation to MNI space — nearest neighbor for labels
        seg_mni = ants.apply_transforms(
            fixed=MNI,
            moving=seg,
            transformlist=reg["fwdtransforms"],
            interpolator="nearestNeighbor",
        )
        seg_np = seg_mni.numpy()

        reg_out = os.path.join(RESULTS_DIR, "registrations", f"{subject_id}_seg_MNI.nii.gz")
        ants.image_write(seg_mni, reg_out)

        wt  = seg_np > 0
        et  = seg_np == 3
        tc  = (seg_np == 1) | (seg_np == 3)

        result = {
            "subject":  subject_id,
            "wt_vol":   int(wt.sum()),
            "et_vol":   int(et.sum()),
            "elapsed_s": round(time.time() - t0, 1),
            "error":    None,
        }

        for mask_name, cc_arr in CC_MASKS.items():
            for label_name, tumor_arr in [("wt", wt), ("et", et), ("tc", tc)]:
                overlap = int((tumor_arr & cc_arr).sum())
                result[f"cc_{mask_name}_{label_name}_overlap"] = overlap
                result[f"cc_{mask_name}_{label_name}_involved"] = overlap > 0
                result[f"cc_{mask_name}_{label_name}_frac"] = round(overlap / CC_VOL[mask_name], 4)

        # Butterfly: bilateral involvement in CC (left x<mid, right x>=mid)
        mid = seg_np.shape[0] // 2
        cc_whole = CC_MASKS["whole"]
        et_cc = et & cc_whole
        result["butterfly_et"] = bool(et_cc[:mid].sum() > 5 and et_cc[mid:].sum() > 5)
        result["butterfly_wt"] = bool(
            (wt & cc_whole)[:mid].sum() > 20 and (wt & cc_whole)[mid:].sum() > 20
        )

        with open(out_json, "w") as f:
            json.dump(result, f, indent=2)

        print(f"  {subject_id}: CC_wt={result['cc_whole_wt_involved']} "
              f"CC_et={result['cc_whole_et_involved']} "
              f"butterfly_ET={result['butterfly_et']}  [{result['elapsed_s']}s]")
        return result

    except Exception as e:
        err = {"subject": subject_id, "error": str(e)}
        with open(out_json, "w") as f:
            json.dump(err, f)
        print(f"  {subject_id}: ERROR {e}")
        return err


if __name__ == "__main__":
    print(f"TRAIN_DIR  = {TRAIN_DIR}")
    print(f"ATLAS_DIR  = {ATLAS_DIR}")
    print(f"RESULTS    = {RESULTS_DIR}\n")

    phgg_ids = get_phgg_ids()
    print(f"Running CC registration pipeline on {len(phgg_ids)} pHGG cases")
    print(f"CPUs available: {cpu_count()}")

    TEST_ONLY = False
    ids_to_run = phgg_ids[:3] if TEST_ONLY else phgg_ids

    n_workers = max(1, cpu_count() - 2)
    print(f"Using {n_workers} parallel workers\n")

    with Pool(n_workers) as pool:
        all_results = pool.map(register_case, ids_to_run)

    valid = [r for r in all_results if not r.get("error")]
    cc_wt = sum(1 for r in valid if r.get("cc_whole_wt_involved"))
    cc_et = sum(1 for r in valid if r.get("cc_whole_et_involved"))
    butterfly = sum(1 for r in valid if r.get("butterfly_et"))
    n = max(1, len(valid))

    print(f"\n{'='*50}")
    print(f"Results: {len(valid)}/{len(all_results)} cases completed")
    print(f"CC involvement (whole tumor): {cc_wt}/{len(valid)} ({100*cc_wt/n:.1f}%)")
    print(f"CC involvement (ET only):     {cc_et}/{len(valid)} ({100*cc_et/n:.1f}%)")
    print(f"Butterfly pattern (ET):       {butterfly}/{len(valid)} ({100*butterfly/n:.1f}%)")

    import pandas as pd
    pd.DataFrame(valid).to_csv(os.path.join(RESULTS_DIR, "cc_results_phgg.csv"), index=False)
    print(f"\nResults saved to {RESULTS_DIR}/cc_results_phgg.csv")
