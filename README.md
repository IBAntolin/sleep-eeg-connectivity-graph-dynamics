# sleep-eeg-connectivity-graph-dynamics

Exploratory study of sleep EEG functional connectivity across temporal scales:
graph-theoretical analysis of stage-mean wPLI networks and self-supervised
GATv2 autoencoder learning on 30-second delta-band connectivity graphs to assess
within-stage temporal stability.

## Phase 1: EEG processing, connectivity construction, and quality control

This repository begins with a reproducible preprocessing and quality-control
workflow that transforms overnight ANPHY-Sleep EEG recordings into
frequency-specific functional-connectivity graphs.

The purpose of this first phase is to create an auditable and reusable dataset
for later graph-theoretical, temporal-dynamics, statistical, and
graph-neural-network analyses. It does not yet test final population-level
hypotheses about sleep-stage differences or report definitive neurophysiological
findings.

Each retained 30-second EEG epoch is represented as a weighted graph:

- **Nodes:** 83 scalp EEG channels.
- **Edges:** weighted phase-lag index (wPLI) values between pairs of EEG
  channels.
- **Frequency bands:** delta, theta, alpha, sigma, and beta.
- **Epoch labels:** scored sleep stage.
- **Audit metadata:** original epoch index, EEG channel names, and number of
  channels interpolated in that epoch.

---

## Phase 1 workflow

```text
ANPHY-Sleep source files
│
├── EDF recording
│   ├── 93 channels recorded at 1000 Hz
│   ├── 83 scalp EEG channels
│   └── 10 auxiliary PSG channels
│
├── Sleep-scoring TXT file
│   └── sleep-stage labels and epoch timing
│
├── Artifact matrix (*.mat)
│   └── clean/bad channel flag for each epoch and electrode
│
└── Electrode-position file (*.pos)
    └── 3D EEG electrode coordinates
        │
        ▼
Memory-efficient per-subject processing
│
├── Stream one 30-second epoch at a time from the EDF
├── Select scalp EEG channels and exclude auxiliary PSG channels
├── Assign each epoch a sleep-stage label
├── Downsample 1000 Hz → 200 Hz
├── Regress chin-EMG-linked components from EEG channels
├── Drop epochs with too many artifact-flagged channels
├── Interpolate a limited number of bad channels when appropriate
└── Immediately discard the time-domain epoch after graph construction
        │
        ▼
Connectivity graph construction
│
├── Delta: 0.5–4 Hz
├── Theta: 4–8 Hz
├── Alpha: 8–12 Hz
├── Sigma: 12–16 Hz
└── Beta: 16–30 Hz
        │
        ▼
One processed HDF5 file per subject
│
├── One 83 × 83 wPLI matrix per retained epoch and band
├── Sleep-stage labels
├── Original epoch indices
├── Channel names
└── Interpolation audit metadata
        │
        ▼
Exploratory quality-control and sensitivity checks
├── HDF5 structure and numerical-range validation
├── Montage/channel-order validation
├── Connectivity visualization
├── Peripheral-channel / possible EMG sensitivity analysis
├── Frontal-channel / possible EOG sensitivity analysis
└── Interpolation-related connectivity checks
```

---

## Data inputs

Phase 1 uses the following ANPHY-Sleep source files:

| File type | Role in this pipeline |
|---|---|
| `.edf` | Raw overnight EEG and auxiliary polysomnographic recordings |
| Sleep-scoring `.txt` | Sleep-stage labels and start times for epochs |
| `*_artndxn.mat` | Per-epoch, per-channel artifact matrix |
| `.pos` | Co-registered 3D scalp-electrode coordinates |
| Subject-information `.csv` | Participant metadata for potential later analyses |

The EDF recordings contain 93 channels sampled at 1000 Hz:

```text
83 scalp EEG channels:
Used as graph nodes.

10 auxiliary PSG channels:
ChEMG1, ChEMG2,
RLEG-, RLEG+,
LLEG-, LLEG+,
ECG1, ECG2,
EOG1, EOG2
```

The auxiliary channels are excluded from the 83-node connectivity graph. The
chin-EMG channels are used as regressors in the per-epoch EMG-related
correction.

---

## Processing steps

### EDF streaming and epoch selection

The pipeline processes one 30-second segment at a time instead of loading a
complete overnight recording into memory. This keeps memory usage manageable
and allows processing of high-density recordings on ordinary hardware.

Each epoch is linked to a sleep-stage label from the scoring file. The
`light_on_off` marker is excluded by default because it is not a conventional
sleep stage.

### Downsampling

The raw 1000 Hz data are downsampled to 200 Hz using polyphase resampling with
anti-alias filtering.

```text
Original sampling rate: 1000 Hz
Processed sampling rate:  200 Hz
Epoch duration:            30 seconds
Samples per processed epoch: 6,000 per channel
```

The 200 Hz rate is sufficient for the current delta-to-beta analysis while
substantially reducing computational cost.

### Chin-EMG regression

For every epoch, the pipeline can regress EEG variance associated with the
two chin-EMG channels (`ChEMG1`, `ChEMG2`) from all EEG channels. This is done
before interpolation and connectivity estimation.

The aim is to reduce possible muscle-linked signal contributions, particularly
because muscle activity can differ across sleep stages and may otherwise affect
higher-frequency connectivity estimates.

### Artifact handling

The artifact matrix contains one flag per channel per epoch:

```text
1 = clean
0 = artifact-flagged / bad
```

The default policy is:

```text
If more than 15% of EEG channels are flagged bad:
    drop the full epoch.

If 15% or fewer are flagged bad:
    retain the epoch and interpolate flagged channels using
    MNE spherical-spline interpolation and the electrode coordinates.
```

The number of interpolated channels is saved for every retained epoch. This
provides an audit trail for later sensitivity analyses.

### wPLI graph construction

For every retained epoch, the pipeline calculates a wPLI matrix in five
frequency bands:

| Band | Frequency range |
|---|---:|
| Delta | 0.5–4 Hz |
| Theta | 4–8 Hz |
| Alpha | 8–12 Hz |
| Sigma | 12–16 Hz |
| Beta | 16–30 Hz |

Each matrix describes pairwise phase-lagged connectivity between the 83 scalp
EEG channels:

```text
Matrix shape: 83 × 83
Value range:  0 to 1
Diagonal:     0
```

The pipeline calculates graphs immediately after preprocessing each epoch, then
discards the time-domain epoch. Only the graph representation and associated
metadata are retained.

---

## Processed outputs

One HDF5 file is written per subject:

```text
graphs_output/
├── EPCTL01.h5
├── EPCTL02.h5
├── ...
└── EPCTL29.h5
```

Each HDF5 file contains:

| Item | Description |
|---|---|
| `subject_id` | Subject identifier |
| `n_epochs` | Number of retained epochs |
| `bands` | Frequency bands stored in the file |
| `max_bad_fraction` | Bad-channel threshold used during processing |
| `epoch_idx` | Original scoring/artifact-matrix index for each retained epoch |
| `stage` | Sleep-stage label for each retained epoch |
| `channel_names` | EEG channel name corresponding to every matrix row and column |
| `n_bad_channels` | Number of channels interpolated in each retained epoch |
| `delta` | Delta-band wPLI matrices |
| `theta` | Theta-band wPLI matrices |
| `alpha` | Alpha-band wPLI matrices |
| `sigma` | Sigma-band wPLI matrices |
| `beta` | Beta-band wPLI matrices |

Each band dataset has the following layout:

```text
(number of retained epochs, 83, 83)
```

For example:

```python
result["delta"][epoch_number, i, j]
```

is the delta-band wPLI value between:

```python
result["channel_names"][i]
```

and:

```python
result["channel_names"][j]
```

for one retained epoch.

---

## Files in this phase

| File | Purpose |
|---|---|
| `download_data.py` | Downloads and extracts ANPHY-Sleep source files from OSF, including subject archives, artifact matrices, electrode positions, and subject information |
| `preprocessing.py` | Main processing module: EDF streaming, epoch labeling, downsampling, chin-EMG regression, artifact handling, interpolation, HDF5 output, and HDF5 loading |
| `connectivity.py` | Defines frequency bands and calculates per-epoch wPLI connectivity matrices |
| `sanity_check.py` | Checks HDF5 structure, expected dimensions, stage labels, numerical wPLI ranges, graph diagonals, and interpolation metadata |
| `peripheral_channel_check.py` | Tests whether peripheral/outer-ring channels are disproportionately represented among strong connectivity edges |
| `eog_check.py` | Tests whether frontal and eye-adjacent electrodes are over-represented among strong low-frequency edges, including stage-aware checks |
| `check_max_wpli_source.py` | Performs interpolation-related diagnostic checks, including maximum-edge inspection, two-way and three-way edge comparisons, bad-channel summaries, and distance-based checks |
| `main_1_2.ipynb` | Main notebook for initial EDF inspection, processing orchestration, and generation of subject-level HDF5 graph files |
| `main_checks.ipynb` | Notebook for exploratory visualization, processed-data inspection, and artifact/interpolation sensitivity checks |
| `guideline.ipynb` | Notes the project aim and the intended architecture for later phases |

---

## Quality-control checks

The goal of these checks is not to establish final neuroscience results. Their
role is to verify that the graph dataset is sensible, interpretable, and
methodologically transparent before later analysis.

### HDF5 sanity checks

`sanity_check.py` confirms that every processed subject file has:

- The expected datasets and attributes.
- 83 EEG channels.
- Matching numbers of stage labels and epoch indices.
- Consistent graph shapes across frequency bands.
- Finite sampled wPLI values within the expected 0–1 range.
- Approximately zero diagonals.
- Readable stage distributions.
- Evidence that artifact filtering/interpolation metadata were recorded.

Example:

```python
from sanity_check import sanity_check

sanity_check("graphs_output/EPCTL01.h5")
```

### Montage and graph visualization checks

`main_checks.ipynb` includes an exploratory single-epoch scalp plot of the
strongest wPLI edges across the five frequency bands.

This plot is used to:

- Confirm that channel names and electrode positions are aligned.
- Confirm that the generated graphs are visually interpretable.
- Inspect the names and spatial locations of the strongest connections.
- Identify unexpected spatial patterns that should be investigated
  systematically.

This is an illustrative visualization only. It does not provide group-level
evidence for sleep-stage effects or artifact contamination.

### Subject-level descriptive plots

The same notebook includes:

- Violin plots of epoch-level mean wPLI grouped by sleep stage and band.
- A temporal plot of mean wPLI across epochs together with a hypnogram.
- A heatmap of mean wPLI by frequency band and sleep stage.

These figures are useful for identifying temporal discontinuities, unusual
periods, stage-associated trends, and future hypotheses. They are descriptive
subject-level checks, not population-level statistical analyses.

### Peripheral-channel / possible EMG sensitivity check

`peripheral_channel_check.py` investigates whether the strongest connectivity
edges are disproportionately associated with outer-ring, temporal, zygomatic,
supraorbital, or jaw-adjacent electrodes.

It distinguishes:

```text
both_central:
    neither endpoint is peripheral

one_peripheral:
    exactly one endpoint is peripheral

both_peripheral:
    both endpoints are peripheral
```

This is particularly relevant for beta-band connectivity because residual
muscle activity may be stronger at peripheral electrodes and can vary with
sleep stage.

A later beta-band finding should therefore be repeated after excluding
peripheral channels. Agreement between full-montage and peripheral-excluded
results would increase confidence; disappearance of an effect would indicate
that the original result should be interpreted cautiously.

### Frontal-channel / possible EOG sensitivity check

`eog_check.py` investigates whether frontal, eye-adjacent channels are
over-represented among the strongest delta- and theta-band edges.

The check can be run across subjects or within individual subjects while
comparing REM against non-REM epochs. Its purpose is to identify whether
low-frequency connectivity patterns might be influenced by residual
eye-movement or blink-related activity.

### Interpolation sensitivity checks

`check_max_wpli_source.py` investigates whether bad-channel interpolation could
systematically influence wPLI values.

The available checks are:

| Check | Purpose |
|---|---|
| Maximum-edge inspection | Identifies the global maximum wPLI edge and checks whether either channel was interpolated in that epoch |
| Two-way edge comparison | Compares edges touching an interpolated channel against edges with two clean endpoints |
| Three-way edge comparison | Separates `both_clean`, `one_interpolated`, and `both_interpolated` edges |
| Bad-channel summary | Describes how many channels were flagged bad per epoch for selected subjects |
| Distance-based check | Tests whether wPLI involving one interpolated channel is stronger for nearby scalp electrodes |
| Channel-order plausibility check | Ranks channels by bad-channel fraction to assess whether artifact-matrix column order appears plausible |

These are sensitivity analyses and diagnostic checks, not definitive proof of
an interpolation effect. The artifact matrices do not contain channel-name
metadata, so their column order is assumed to match the channel order used in
the processed graph files.

---

## Installation

Install the required Python packages:

```bash
pip install numpy pandas scipy h5py pyedflib mne requests tqdm matplotlib
```

The main packages are used as follows:

| Package | Purpose |
|---|---|
| `numpy` | Numerical arrays and connectivity matrices |
| `pandas` | Sleep-scoring tables and metadata |
| `scipy` | Resampling, filtering, and Hilbert transform |
| `h5py` | HDF5 connectivity files and MATLAB v7.3 artifact matrices |
| `pyedflib` | Memory-efficient EDF access |
| `mne` | Bad-channel interpolation and montage handling |
| `requests` | Data download from OSF |
| `tqdm` | Progress bars during processing |
| `matplotlib` | Quality-control and exploratory plots |

---

## How to run Phase 1

### 1. Download shared data files

```python
from download_data import download_shared_files

DATA_ROOT = "anphy_sleep_data"

download_shared_files(DATA_ROOT)
```

This downloads the artifact matrices, electrode-position file, and subject
information.

### 2. Download one subject

Start with one subject as a small test:

```python
from download_data import download_subject

download_subject(1, DATA_ROOT)
```

### 3. Inspect the EDF layout

```python
import pyedflib

with pyedflib.EdfReader("anphy_sleep_data/EPCTL01/EPCTL01.edf") as f:
    labels = f.getSignalLabels()

aux_channel_names = {
    "ChEMG1", "ChEMG2",
    "RLEG-", "RLEG+",
    "LLEG-", "LLEG+",
    "ECG1", "ECG2",
    "EOG1", "EOG2",
}

eeg = [ch for ch in labels if ch not in aux_channel_names]
aux = [ch for ch in labels if ch in aux_channel_names]

print(f"EEG channels: {len(eeg)}")
print(f"Auxiliary channels: {len(aux)}")
print("Auxiliary channels:", aux)
```

Expected result:

```text
EEG channels: 83
Auxiliary channels: 10
```

### 4. Process a short test run

Use `max_epochs` for a short smoke test before processing a full night:

```python
from pathlib import Path

from download_data import (
    get_subject_artifact_path,
    get_subject_edf_path,
    get_subject_scoring_path,
)
from preprocessing import process_subject

DATA_ROOT = Path("anphy_sleep_data")
OUTPUT_DIR = Path("graphs_output")

subject_num = 1
subject_id = f"EPCTL{subject_num:02d}"

edf_path = get_subject_edf_path(subject_num, DATA_ROOT)
scoring_path = get_subject_scoring_path(subject_num, DATA_ROOT)
artifact_path = get_subject_artifact_path(subject_num, DATA_ROOT)
pos_path = DATA_ROOT / "co-registered_average_positions.pos"

output_path = process_subject(
    subject_id=subject_id,
    edf_path=edf_path,
    scoring_path=scoring_path,
    artifact_path=artifact_path,
    output_dir=OUTPUT_DIR,
    pos_path=pos_path,
    max_epochs=10,
)

print(f"Created: {output_path}")
```

After the test passes, remove:

```python
max_epochs=10
```

to process the full recording.

### 5. Validate processed output

```python
from sanity_check import sanity_check

sanity_check("graphs_output/EPCTL01.h5")
```

### 6. Run exploratory and quality-control checks

Open:

```text
main_checks.ipynb
```

This notebook contains the scalp visualizations, temporal/stage summaries, and
the artifact-related sensitivity checks described above.

---

## Data policy

Raw and generated participant-level data should remain local and should not be
committed to the Git repository.

Recommended `.gitignore` entries:

```gitignore
# Raw and generated EEG data
anphy_sleep_data/
graphs_output/
*.edf
*.h5
*.mat
*.zip

# Jupyter
.ipynb_checkpoints/

# Python cache and environments
__pycache__/
*.py[cod]
.venv/
venv/
.env

# Editor and operating-system files
.vscode/
.DS_Store
```

The repository should contain the code, notebooks, documentation, and
dependency information required to reproduce the outputs. A future frozen data
release, if needed, should be archived separately through an appropriate
research-data repository rather than stored in Git.

---

## Phase 1 outcome

Phase 1 produces a structured collection of subject-level, frequency-specific
wPLI connectivity graphs with sleep-stage labels and preprocessing audit
metadata.

The completed output is suitable for the next project phase, including:

- Subject-aware stage-mean connectivity aggregation.
- Weighted graph-theoretical feature extraction.
- Temporal network-dynamics analysis.
- Statistical comparison across sleep stages.
- Self-supervised graph representation learning on 30-second connectivity
  graphs.

The present phase does not yet treat the exploratory plots or artifact checks as
final population-level neuroscience findings. Their role is to establish a
transparent, inspectable, and reproducible foundation for the analyses that
follow.
