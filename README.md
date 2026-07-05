# Corpus Callosum Involvement — Registration Pipeline (pHGG vs Adult GBM)

Registers pediatric high-grade glioma (pHGG, BraTS-PEDs) and adult glioblastoma
(GBM, BraTS-GLI) MRI cases to MNI space and measures tumor overlap with the
corpus callosum (CC) and its sub-regions (genu, body, splenium).

This repo covers the registration stage only: from data download through
per-case CC-overlap metrics and a results CSV. Statistics and figures are done
separately.

## Pipeline

For each subject the pipeline:

1. Registers the T1c image to the ANTs MNI152 template (SyNRA), with a
   tumor-excluded cost-function mask so the tumor does not bias alignment.
2. Warps the tumor segmentation into MNI space (nearest neighbor).
3. Computes overlap of WT / ET / TC with the whole CC and each sub-region, plus
   bilateral "butterfly" flags.

Output per cohort: one `<subject>_cc.json` per case and a combined
`cc_results_*.csv`.

## Layout

```
PEDs-HGG/
├── README.md
├── requirements.txt
├── scan_gli_et.py              # filter GLI to GBM-like cases (ET > 100)
├── download_gbm_synapse.py     # download adult GBM (BraTS-GLI) from Synapse
├── cc_registration_peds.py     # register pHGG -> cc_results_peds/
├── cc_registration_gbm.py      # register GBM  -> cc_results_gbm/
├── atlas/                      # CC masks (committed)
│   ├── CC_mask.nii.gz
│   ├── CC_genu.nii.gz
│   ├── CC_body.nii.gz
│   └── CC_splenium.nii.gz
└── data/
    ├── histologies.tsv             # OpenPedCan v15 (pHGG cohort selection)
    ├── BraTS-PEDs_metadata.tsv     # BraTS-PEDs metadata
    ├── BraTS-PEDs-v1/  (download separately, ~33 GB)
    └── BraTS-GLI/      (download separately, ~12 GB)
```

The two imaging folders are not tracked in git. Download them as below.

## 1. Environment

Python 3.11 is recommended (best antspyx wheel coverage; 3.9 often forces a slow
source build).

```bash
cd PEDs-HGG
python3.11 -m venv venv_ants
source venv_ants/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

Check:

```bash
python -c "import ants, nibabel, numpy, pandas; print('all OK')"
```

If antspyx tries to build from source, force a prebuilt wheel:

```bash
pip install antspyx==0.6.3 --only-binary=:all:
```

If no wheel exists for your Python, use Python 3.11.

## 2. Data

### Pediatric HGG — BraTS-PEDs (TCIA)

https://www.cancerimagingarchive.net/collection/brats-peds/

Download the Images archive (~33 GB, NIfTI) and extract so case folders land
under:

```
data/BraTS-PEDs-v1/Training/BraTS-PED-00001-000/BraTS-PED-00001-000-t1c.nii.gz
                                                 BraTS-PED-00001-000-seg.nii.gz
```

`data/BraTS-PEDs_metadata.tsv` and `data/histologies.tsv` are already included
and are used for cohort selection.

### Adult GBM — BraTS-GLI (Synapse)

Needs a free Synapse account and a Personal Access Token (synapse.org → Account
Settings → Personal Access Tokens).

```bash
export SYNAPSE_AUTH_TOKEN=your_token_here
python download_gbm_synapse.py
```

Or run without the env var to log in interactively. Extract so GLI case folders
resolve under `data/BraTS-GLI/` (a nested challenge folder is handled).

Only `-t1c.nii.gz` and `-seg.nii.gz` are read per case.

## 3. Run

Activate the venv first: `source venv_ants/bin/activate`

Pediatric HGG:

```bash
python cc_registration_peds.py
```

Prints `Running CC registration pipeline on N pHGG cases` (N should be 104).
Writes `cc_results_peds/*_cc.json`, warped segs, and `cc_results_phgg.csv`. Set
`TEST_ONLY = True` near the bottom to run only 3 cases first.

Adult GBM (two steps, in order):

```bash
python scan_gli_et.py            # step 1: define the GBM cohort
python cc_registration_gbm.py    # step 2: register that cohort
```

Step 1 — `scan_gli_et.py` selects which GLI cases count as GBM. The raw BraTS-GLI
download contains ~1,251 cases of mixed grade. This script reads every
segmentation, counts enhancing-tumor (ET) voxels, and keeps only cases with
ET > 100 (enhancement used as a GBM proxy; cases with no/near-zero enhancement are
likely low-grade and dropped). It writes the surviving subject IDs (~1,213) to
`cc_results_gbm/gbm_cases.txt`. It runs once, is fast (no registration), and must
be run before step 2 — `cc_registration_gbm.py` reads `gbm_cases.txt` and will
exit with an error if it is missing.

Step 2 — `cc_registration_gbm.py` registers each case in `gbm_cases.txt` and
writes `cc_results_gbm/*_cc.json`, warped segs, and `cc_results_gbm.csv`.

(The pediatric run needs no such filter — its cohort comes from the metadata and
histology labels, not from an enhancement threshold.)

Both scripts skip any case that already has a `_cc.json`, so an interrupted run
resumes on re-run. Delete error JSONs first so they retry:

```bash
python -c "import json,glob,os;[os.remove(f) for f in glob.glob('cc_results_peds/*_cc.json') if json.load(open(f)).get('error')]"
```

## Registration output (descriptive)

CC involvement counts from the registration stage, before any statistical
analysis. WT = whole tumor, TC = tumor core, ET = enhancing tumor.

| Metric | pHGG (n=104) | Adult GBM (n=1213) |
|--------|--------------|--------------------|
| CC whole — WT | 62 (59.6%) | 623 (51.4%) |
| CC whole — TC | 34 (32.7%) | 341 (28.1%) |
| CC whole — ET | 16 (15.4%) | 327 (27.0%) |
| CC genu — TC | 27 (26.0%) | 179 (14.8%) |
| CC genu — WT | 57 (54.8%) | 467 (38.5%) |
| CC body — TC | 17 (16.3%) | 216 (17.8%) |
| CC splenium — WT | 10 (9.6%) | 275 (22.7%) |
| CC splenium — TC | 5 (4.8%) | 173 (14.3%) |
| Butterfly WT | 45 (43.3%) | 422 (34.8%) |
| Butterfly ET | 6 (5.8%) | 208 (17.1%) |

Counts are cases with any overlap in that region. ET metrics reflect that a large
share of pHGG is non-enhancing (et_vol = 0), unlike the enhancement-filtered GBM
cohort. Significance testing is done in the separate analysis stage.

## Registration to MNI152 Space

Sample cases warped to the MNI152 template. Red = warped tumor (whole tumor) in
MNI space; cyan = corpus callosum contour. Each panel shows sagittal, coronal,
and axial views at the midline. These span CC-involved and non-involved cases so
the computed overlap flags can be checked against anatomy.

| Case | CC whole (WT) | CC genu (TC) | Butterfly (WT) |
|------|---------------|--------------|----------------|
| BraTS-PED-00033-000 | yes | yes | yes |
| BraTS-PED-00001-000 | yes | yes | yes |
| BraTS-PED-00043-000 | yes | no  | yes |
| BraTS-PED-00128-000 | no  | no  | no |
| BraTS-PED-00151-000 | no  | no  | no |

![QC BraTS-PED-00033-000](mni-mapping-sample/QC_BraTS-PED-00033-000.png)
![QC BraTS-PED-00001-000](mni-mapping-sample/QC_BraTS-PED-00001-000.png)
![QC BraTS-PED-00043-000](mni-mapping-sample/QC_BraTS-PED-00043-000.png)
![QC BraTS-PED-00128-000](mni-mapping-sample/QC_BraTS-PED-00128-000.png)
![QC BraTS-PED-00151-000](mni-mapping-sample/QC_BraTS-PED-00151-000.png)

## Notes

- Each case takes about 4–8 min; scripts use (CPU cores − 2) workers.
- ANTs writes large temporary warp fields to `$TMPDIR`. Keep several GB free, or
  set `export TMPDIR=/path/with/space`.
- The pHGG run takes roughly an hour on a modern multi-core Mac. The GBM run
  (~1213 cases) is much longer.
- Outputs are small (~5 MB for 104 pHGG cases); the imaging inputs are the large
  part.

## Data sources

- `atlas/CC_*.nii.gz` — from the JHU-ICBM DTI-81 White-Matter Labels atlas
  (labels 3=genu, 4=body, 5=splenium), MNI152 space. References: Mori et al.,
  NeuroImage 2008; Wakana et al., NeuroImage 2007.
- `data/histologies.tsv` — OpenPedCan v15 release (D3b Center), used to classify
  CBTN cases as high-grade for cohort selection.
- `data/BraTS-PEDs_metadata.tsv` — BraTS-PEDs metadata (Source, cohort split,
  MappingID).
