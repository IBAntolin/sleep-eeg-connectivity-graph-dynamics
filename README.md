# sleep-eeg-connectivity-graph-dynamics

Exploratory study of sleep EEG functional connectivity across temporal scales:
graph-theoretical analysis of stage-mean wPLI networks and self-supervised
GATv2 autoencoder learning on 30-second delta-band connectivity graphs to assess
within-stage temporal stability.

## Introduction

Sleep is not a uniform state. Across a night, the brain transitions between
wakefulness, light non-REM sleep, deep non-REM sleep, and REM sleep. These
states differ in neural oscillations, large-scale coordination, sensory
responsiveness, physiological regulation, and information processing.

This project investigates sleep-stage transitions using a network-neuroscience framework. 
Rather than examining EEG channels independently, it represents the scalp EEG
recording as a functional network in which electrodes are nodes and phase-based
statistical relationships between pairs of electrodes are weighted edges. This
makes it possible to ask how functional-network organization differs across
sleep stages and how it evolves over time within a nominally stable stage.

The initial connectivity representation is based on weighted phase-lag index
(wPLI). For each scored 30-second EEG epoch, wPLI estimates the consistency of
non-zero-lag phase relationships between pairs of scalp EEG channels. Each
epoch is therefore transformed into a weighted connectivity matrix and can be
treated as a functional brain graph rather than only as a collection of
separate channel time series.

```text
One scored 30-second EEG epoch
        ↓ 
Band-specific wPLI connectivity estimation
        ↓
One weighted EEG functional-connectivity graph
        ↓
Static graph analysis and temporal graph-representation analysis
```

The project combines two complementary analytical time scales.

At the **stage-mean scale**, all valid 30-second graphs belonging to the same
participant, sleep stage, and frequency band are averaged to produce one
stage-representative functional-connectivity graph. Classical graph-theoretical
metrics are then used to describe network integration, local segregation,
community-like structure, and threshold robustness across Wake, N1, N2, N3,
and REM sleep.

At the **epoch-resolved scale**, individual 30-second delta-band graphs are
analysed with a self-supervised GATv2 graph autoencoder. The model learns a
compact representation of graph structure without using sleep stages as direct
training labels. Distances between embeddings of adjacent 30-second graphs are
then used to quantify within-stage temporal stability.

```text
Stage-mean graph analysis:
    How does average functional-network topology differ across sleep stages?

Epoch-resolved GATv2 graph analysis:
    How stable is the learned functional-network representation from one
    30-second epoch to the next within the same sleep stage?
```

## Research questions

This exploratory study addresses the following questions:

- How do frequency-specific EEG functional-connectivity networks differ across
  Wake, N1, N2, N3, and REM sleep?
- Do stage-mean networks differ in graph-theoretical measures of integration,
  local segregation, efficiency, characteristic path length, modularity, or
  community structure?
- Are descriptive static-network patterns robust to reasonable graph-density
  choices?
- How stable are 30-second delta-band functional-connectivity graphs within a
  nominally homogeneous sleep stage?
- Can a self-supervised GATv2 autoencoder learn graph representations that
  capture meaningful within-stage temporal variation?
- Do static graph-theoretical measures and learned graph embeddings provide
  complementary descriptions of sleep-network dynamics?

## Project workflow

```text
ANPHY-Sleep recordings and annotations
│
├── Overnight high-density EEG recordings
├── Sleep-stage annotations
├── Channel-level artifact matrices
└── Electrode-position information
        │
        ▼
1. EEG processing and connectivity construction
│
├── Stream 30-second EEG epochs from EDF recordings
├── Align epochs with sleep-stage labels
├── Downsample 1000 Hz → 200 Hz
├── Regress chin-EMG-linked components
├── Drop heavily artifact-contaminated epochs
├── Interpolate a limited number of bad channels
└── Compute wPLI matrices in delta, theta, alpha, sigma, and beta bands
        │
        ▼
2. Quality control and artifact sensitivity analysis
│
├── Validate HDF5 output structure and wPLI value ranges
├── Check montage and channel-order alignment
├── Inspect representative connectivity graphs
├── Check peripheral-channel enrichment and possible EMG sensitivity
├── Check frontal-channel enrichment and possible EOG sensitivity
└── Check possible interpolation-related effects
        │
        ▼
3. Classical stage-mean graph-theoretical analysis
│
├── Average graphs within subject × stage × band
├── Calculate weighted graph metrics
├── Evaluate unthresholded weighted networks
├── Evaluate thresholded network robustness
├── Use connected thresholds for path-based integration metrics
└── Compare descriptive within-subject patterns across stages
        │
        ▼
4. Phase B: self-supervised graph representation learning
│
├── Construct thresholded individual delta-band graphs
├── Train a GATv2 graph autoencoder using participant-level data splits
├── Learn a graph embedding for each 30-second epoch
├── Compare embeddings of adjacent stable-stage epochs
├── Summarize within-stage stability at the participant level
└── Perform exploratory repeated-measures permutation testing
        │
        ▼
Interpretation of static and temporal sleep-network organization
```

## Example connectivity representation

The figure below illustrates the graph representation used throughout the
project. It shows the strongest within-band wPLI connections from one
30-second N1 epoch of one of the subjects.

![Illustrative multi-band scalp connectivity graph from one 30-second N1 epoch
of participant EPCTL01. Nodes are scalp EEG electrodes. Lines represent the
strongest 4% of wPLI edges within each frequency band. A shared colour scale
indicates absolute wPLI magnitude across delta, theta, alpha, sigma, and beta.
This single-epoch visualization demonstrates the graph representation used in
the project and is not a group-level result.](figures/scalp_graph.png)

Only the strongest 4% of edges are shown separately within each frequency band
to keep the scalp networks readable, and the common colour scale allows the
strength of the displayed edges to be compared across bands.

## Current project scope

The repository currently contains an end-to-end exploratory workflow from raw
overnight EEG recordings to static and temporal analyses of functional
connectivity.

### Implemented components

- Memory-efficient processing of overnight high-density EEG recordings.
- Frequency-specific wPLI graph construction for delta, theta, alpha, sigma,
  and beta bands.
- HDF5 storage of participant-level epoch graphs, sleep-stage labels, channel
  names, original epoch indices, and interpolation audit information.
- Data-quality, montage, EMG-related peripheral-channel, EOG-related frontal
  channel, and interpolation sensitivity checks.
- Stage-mean graph construction for each participant × sleep stage × frequency
  band combination.
- Unthresholded weighted graph-theoretical analysis.
- Proportional-threshold robustness analysis using sparse and fully connected
  graph-density families.
- Static graph metrics including mean wPLI, node strength, weighted clustering,
  weighted global efficiency, characteristic path length, Louvain modularity,
  and community count.
- A beta-band peripheral–peripheral edge-exclusion sensitivity specification.
- Self-supervised GATv2 graph-autoencoder training on individual delta-band
  connectivity graphs.
- Participant-level within-stage temporal-stability analysis based on distances
  between adjacent graph embeddings.
- Exploratory repeated-measures Friedman permutation analysis of stage-related
  differences in learned within-stage stability.

### Project outputs

The workflow produces three main categories of output:

```text
1. Epoch-level connectivity graphs
   One wPLI matrix per retained 30-second epoch, band, and participant.

2. Stage-mean graphs and static graph metrics
   One stage-mean graph and one set of graph-theoretical measures per
   participant × sleep stage × frequency band combination.

3. Learned graph representations and within-stage stability measures
   One graph embedding per selected delta-band epoch, plus participant-level
   summaries of embedding-distance stability within each sleep stage.
```

### Interpretation principles

This is an exploratory research project. The repository aims to make the
analysis reproducible, inspectable, and scientifically cautious.

- Raw epoch pairs are not treated as independent participant observations.
- Stage-mean static graphs are treated as repeated measurements within
  participants.
- Participant-level summaries are used before cohort-level descriptions or
  repeated-measures tests.
- Graph visualizations are interpreted as exploratory and diagnostic unless
  supported by formal subject-aware statistical analysis.
- Artifact-related findings are treated as sensitivity checks rather than
  automatic evidence that a connectivity pattern is biological or artifactual.
- The beta peripheral-edge exclusion is a narrow robustness specification and
  is kept separate from the full-beta analysis.
- GATv2 embedding distance is a learned representation-space measure; it is
  not a direct physiological distance.
- The current statistical analyses are exploratory and should be followed by
  clearly specified pairwise tests, effect sizes, confidence intervals, and
  appropriate multiple-comparison correction where relevant.

> **Scope note:** the repository now implements the complete exploratory
> workflow: preprocessing and quality control, classical static
> graph-theoretical analysis, and self-supervised graph-representation learning
> for within-stage temporal stability. The results remain exploratory and are
> not presented as definitive clinical or neurophysiological conclusions.

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

## Computational environment

This workflow was developed and run on a Windows laptop using Python.
The current preprocessing and classical graph-analysis stages are CPU-based and
do not require a GPU. Exact package requirements are listed in
`requirements.txt`.

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

## Data source

This project uses the **ANPHY-Sleep** dataset: an open overnight
high-density scalp EEG dataset from healthy adults.

- **Dataset on OSF:** [ANPHY-Sleep / OSF project](https://doi.org/10.17605/OSF.IO/R26FH)
- **Dataset paper:** Wei, X. et al. (2024). *ANPHY-Sleep: an Open Sleep
  Database from Healthy Adults Using High-Density Scalp
  Electroencephalogram*. Scientific Data, 11, 896.
  [https://doi.org/10.1038/s41597-024-03722-1](https://doi.org/10.1038/s41597-024-03722-1)

The dataset contains overnight polysomnographic recordings from 29 healthy
adults, including 83-channel high-density scalp EEG together with EOG, EMG,
ECG, electrode-position information, and sleep-scoring annotations.

The raw dataset and derived participant-level connectivity files are not
included in this repository. Users should obtain the original data directly
from the ANPHY-Sleep OSF project and should follow the dataset's terms,
citation requirements, and any applicable ethical or data-use guidance.

The scripts in this repository are designed to download and process the source
data locally:

```python
from download_data import download_shared_files, download_subject

DATA_ROOT = "anphy_sleep_data"

download_shared_files(DATA_ROOT)
download_subject(1, DATA_ROOT)
```
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

---

## Phase 2

## Static graph-theoretical analysis

After preprocessing, each participant has one wPLI connectivity graph for every
retained 30-second epoch and frequency band. This stage of the project moves
from epoch-level graphs to **stage-mean functional-connectivity networks** and
calculates classical static graph-theoretical measures.

The goal is to describe how the topology of EEG functional-connectivity
networks differs across sleep stages and frequency bands, while checking that
observed patterns are not dependent on one arbitrary graph-threshold choice.

```text
Processed epoch-level wPLI graphs
        ↓
Group epochs within each subject by sleep stage and frequency band
        ↓
Average wPLI adjacency matrices element by element
        ↓
One stage-mean graph per subject × stage × band
        ↓
Calculate static graph-theoretical metrics
        ↓
Evaluate threshold sensitivity and prepare data for later statistics
```

This analysis remains exploratory and descriptive at its current stage. The
metric tables and figures are intended to identify robust candidate patterns
for formal subject-aware statistical testing; they do not alone establish
population-level sleep-stage effects.

### Stage-mean graph construction

For every participant, sleep stage, and frequency band, all valid epoch-level
wPLI matrices are averaged element by element.

```text
stage_mean_graph[subject, band, stage]
    = average of all valid 30-second wPLI matrices
      from that subject, band, and stage
```

This produces one weighted, undirected 83 × 83 functional-connectivity graph
for each available combination of:

```text
subject × sleep stage × frequency band
```

The number of retained epochs contributing to each stage-mean graph is stored
as metadata. This is important because stage duration and the number of valid
epochs differ between participants.

The stage-mean graphs are saved in separate HDF5 files, allowing later
statistical, visualization, and machine-learning analyses to reuse exactly the
same graph representation without repeating aggregation.

```text
graphs_output/
└── stage_mean_wpli_graphs_<beta_mode>/
    ├── EPCTL01_stage_means.h5
    ├── EPCTL02_stage_means.h5
    └── ...
```

### Unthresholded weighted analysis

The initial static-metric analysis uses the complete stage-mean weighted wPLI
networks.

```text
All positive wPLI edges are retained.
No proportional threshold is applied.
No absolute wPLI threshold is applied.
The graph diagonal is set to zero.
```

The following metrics are calculated for every subject × stage × band graph:

| Metric | Meaning |
|---|---|
| `mean_wpli` | Average wPLI across all unique channel pairs; overall functional-connectivity strength |
| `mean_node_strength` | Average sum of weighted connections incident on a node |
| `mean_weighted_clustering` | Local weighted neighbourhood connectivity; a descriptive measure related to local segregation |
| `weighted_global_efficiency` | Global integration through short weighted paths |
| `weighted_characteristic_path_length` | Mean shortest weighted path length; lower values correspond to shorter routes |
| `modularity_louvain` | Strength of weighted community-like organization identified with the Louvain method |
| `n_louvain_communities` | Number of Louvain communities found |
| `n_nodes` | Number of graph nodes; expected to be 83 |
| `n_edges` | Number of positive weighted edges |
| `edge_density` | Fraction of possible edges present |
| `n_connected_components` | Number of disconnected components |
| `largest_component_size` | Size of the largest connected component |
| `minimum_nonzero_wpli` | Smallest retained positive wPLI value |

For shortest-path measures, wPLI connection strength is converted to graph
distance:

```text
edge distance = 1 / wPLI
```

A stronger wPLI connection therefore corresponds to a shorter functional route
through the network.

Louvain community detection uses a fixed random seed to make results
reproducible across reruns.

### Beta-band sensitivity specification

Earlier quality-control checks found a beta-specific concern: peripheral
channels were over-represented among the strongest beta-band edges, with the
pattern concentrated in edges connecting two peripheral electrodes. This may
reflect residual muscle-related signal contributions at outer-ring, temporal,
zygomatic, or jaw-adjacent electrodes.

The primary sensitivity version of the beta-band analysis retains all 83 EEG
nodes but removes only edges with **both** endpoints in the predefined
peripheral-channel set.

```text
Total possible undirected edges:       3,403
Both-peripheral beta edges excluded:     231
Fraction of possible edges excluded:    6.8%
Edges still available:                 93.2%
```

This is a deliberately narrow sensitivity correction. It does not assume that
all peripheral beta connectivity is artifactual. Instead, it tests whether
stage-related beta-network patterns remain when the most specifically
identified potential contamination source is excluded.

The selected beta-edge mode is saved in the output HDF5 files and metric
tables:

```text
full
    All positive beta-band edges are retained.

exclude_both_peripheral
    All 83 nodes remain, but beta edges connecting two peripheral channels are
    set to zero.
```

Results from these two versions should be treated as separate analyses and
should not be combined in the same statistical comparison.

### Thresholded robustness analysis

Thresholding is used as a robustness analysis of static network topology.

A proportional threshold retains the same percentage of strongest positive wPLI
edges in each graph:

```text
10% density:
    retain the strongest 10% of possible undirected edges.

85% density:
    retain the strongest 85% of possible undirected edges.
```

This gives every subject × stage × band graph comparable edge density,
regardless of its overall raw wPLI magnitude.

Two threshold families are evaluated because sparse and connected graphs allow
different questions to be addressed.

| Threshold family | Densities | Primary purpose | Metrics interpreted |
|---|---:|---|---|
| Sparse sensitivity | 10%, 20%, 30% | Inspect topology among strongest edges | Weighted clustering, Louvain modularity, community count, component quality control |
| Connected robustness | 85%, 90% | Test integration patterns while preserving full connectivity | Global efficiency, characteristic path length, clustering, modularity |

At sparse densities, many graphs become disconnected. This is expected when
only the strongest edges are retained. Therefore, shortest-path-based measures
are not calculated or interpreted for the sparse threshold family.

At the connected thresholds, all 700 available subject × stage × band graphs
were fully connected:

```text
85% density:
    700 / 700 graphs fully connected
    all graphs contain one connected component of 83 nodes

90% density:
    700 / 700 graphs fully connected
    all graphs contain one connected component of 83 nodes
```

This makes weighted global efficiency and weighted characteristic path length
well-defined for the connected-threshold analysis.  

### Thresholded graph metrics

For thresholded networks, graph connectedness is recorded explicitly:

| Field | Meaning |
|---|---|
| `requested_density` | Proportion of possible edges targeted by the threshold |
| `edge_density` | Observed retained edge density |
| `n_connected_components` | Number of disconnected graph components |
| `largest_component_size` | Number of nodes in the largest component |
| `fully_connected` | Whether all 83 nodes are part of one connected graph |

The primary connected-threshold metrics are:

| Metric | Interpretation |
|---|---|
| `weighted_global_efficiency` | Higher values indicate more efficient communication through strong, short weighted routes |
| `weighted_characteristic_path_length` | Lower values indicate shorter average routes between nodes |
| `mean_weighted_clustering` | Describes local weighted neighbourhood structure |
| `modularity_louvain` | Describes weighted community-like organization |

### Exploratory connected-threshold results

The figure below shows participant-level weighted global efficiency across sleep
stages for the five frequency bands. The upper row retains 85% of possible
edges and the lower row retains 90%.

Grey lines connect repeated measurements from the same participant. Boxplots
and points show the distribution of participant-level stage-mean network
values. All graphs were fully connected at both thresholds, so weighted global
efficiency is well-defined for every participant × stage × band graph.

![Participant-level weighted global efficiency across sleep stages and
frequency bands. Rows show proportional thresholds retaining 85% and 90% of
edges; grey lines connect repeated measurements from the same participant,
while boxplots and points summarize participant-level values. Beta-band graphs
use the peripheral–peripheral edge-exclusion sensitivity specification. This
figure is descriptive and exploratory; it is not a formal statistical
test.](figures/global_eff_2.png)

### Current descriptive findings

The connected-threshold analysis shows qualitatively similar sleep-stage
patterns at both 85% and 90% graph density. This supports robustness of the
descriptive integration patterns to the choice between these two connected
thresholds.  

The main candidate patterns for later formal testing are:

- **Sigma:** networks appear relatively more integrated during N2 and N3 and
  less integrated during REM.
- **Delta:** integration appears relatively higher during Wake and N2 and lower
  during N1 and REM.
- **Alpha:** REM often shows lower integration, while Wake and N3 appear
  relatively more integrated; participant variability is substantial.
- **Theta:** integration appears comparatively stable across sleep stages.
- **Beta:** integration tends to decrease from Wake toward N2/N3 and partially
  recover in REM, but between-participant variation is large.

The principal descriptive observation is that sigma-band networks appear more
integrated during N2/N3 and less integrated during REM. These observations are
not yet inferential conclusions and require subject-aware statistical tests
with appropriate control of repeated measures and multiple comparisons.  

Weighted clustering and Louvain modularity are retained as complementary,
supplementary topology measures. They can provide useful information about
local segregation and community-like organization, but they are more sensitive
to graph density, threshold selection, and community-detection choices than the
primary connected-network integration measures.

### Visualizations

The static graph-analysis notebooks generate participant-level descriptive
figures that combine:

- One trajectory per participant across sleep stages.
- Boxplots summarizing the participant distribution per stage.
- Individual participant values.
- Separate panels for each frequency band.
- Separate rows for 85% and 90% connected thresholds.

These figures are useful for showing individual heterogeneity and for checking
whether apparent stage-related changes are consistent across participants and
across the two connected densities.

Wake-referenced heatmaps are also generated for selected metrics. In these
figures, each cell represents the change from the same participant's Wake
value:

```text
change from Wake
    = metric value for a given stage
      minus
      metric value for Wake in the same subject and band
```

Positive values indicate an increase relative to Wake; negative values indicate
a decrease relative to Wake.

> The figures are exploratory and descriptive. They should not be read as
> group-level statistical evidence until the next analysis stage applies
> appropriate within-subject statistical modelling.

### Generated tables

The thresholded analysis writes separate CSV tables for sparse and connected
threshold families:

```text
threshold_metrics_sparse_10_20_30_exclude_both_peripheral_beta.csv

threshold_metrics_connected_85_90_exclude_both_peripheral_beta.csv
```

Each row corresponds to:

```text
subject × frequency band × sleep stage × threshold density
```

The tables include the subject ID, sleep stage, frequency band, number of
epochs used to form the stage-mean graph, beta-edge mode, graph-connectivity
quality-control fields, and static graph metrics.

These tables are the primary numerical outputs for the next stage:
subject-aware statistical comparison of sleep-stage effects.

### Interpretation notes

- Each subject × stage × band graph is one repeated-measures observation; it is
  not an independent epoch-level sample.
- Later group-level analysis must account for repeated observations within
  participants.
- The number of epochs averaged for a stage graph should be considered in
  sensitivity analyses, because stage coverage differs between participants.
- Sparse graphs may be disconnected; do not interpret path-length or global
  efficiency values from disconnected sparse networks.
- Agreement between 85% and 90% density supports robustness to those connected
  thresholds, but does not prove robustness to every possible graph
  construction choice.
- Do not mix full-beta and peripheral-edge-excluded beta results in one
  analysis.
- Current results are descriptive. Formal statistical inference remains a
  subsequent step.

  ---

## Phase 3

## Phase B: self-supervised graph representation learning and within-stage dynamics

The final project component examines sleep EEG connectivity at the temporal
resolution of individual 30-second epochs.

The earlier static graph-theoretical analysis aggregates all valid epochs from a
participant, sleep stage, and frequency band into one stage-mean network. That
approach is useful for describing average network topology, but it cannot show
whether the connectivity graph remains stable or changes substantially from one
30-second epoch to the next within the same scored sleep stage.

Phase B therefore uses a self-supervised GATv2 graph autoencoder to learn a
compact representation of individual delta-band functional-connectivity graphs.
The learned graph embeddings are then used to quantify short-timescale
within-stage temporal stability.

```text
30-second delta-band wPLI connectivity graph
        ↓
Thresholded graph representation
        ↓
Self-supervised GATv2 graph autoencoder
        ↓
Learned low-dimensional graph embedding
        ↓
Compare embeddings of adjacent 30-second epochs
        ↓
Quantify within-stage temporal stability
```

### Aim

The main question is:

> Does the temporal stability of the learned delta-band functional-connectivity
> representation differ across Wake, N1, N2, N3, and REM sleep?

The analysis does not treat conventional sleep stages as perfectly homogeneous.
Instead, it asks whether consecutive 30-second graphs assigned to the same
stage remain similarly represented by the learned graph model.

This provides a complementary perspective to stage-mean graph metrics:

```text
Static graph analysis:
    How does average network topology differ between sleep stages?

Phase B graph-embedding analysis:
    How stable is the network representation from one 30-second epoch
    to the next within the same sleep stage?
```

### Input graphs

Phase B uses individual 30-second delta-band wPLI connectivity matrices
created during preprocessing.

The delta band is defined as:

```text
0.5–4 Hz
```

Each epoch is represented as an 83-node weighted functional-connectivity graph:

```text
Nodes:
    83 scalp EEG channels

Edges:
    delta-band wPLI connectivity values between electrode pairs

Graph source:
    subject-level processed HDF5 files created by the preprocessing pipeline
```

The graph representation uses proportional edge retention. The trained
autoencoder configuration records the selected edge density and epoch sampling
stride in its output filenames and metadata.

```text
Example configuration:
    Epoch stride: 10
    Retained graph density: 20%
```

The selected graph-construction configuration should be kept fixed when
interpreting embeddings and within-stage distances.

### Data partitioning

The autoencoder uses participant-level data partitioning.

```text
Training participants:
    Used to optimize model parameters.

Validation participants:
    Used for model selection and selection of the best validation epoch.

Held-out test participants:
    Not used during fitting, validation, or model selection.
    Used only for an independent descriptive evaluation.
```

Participant-level splitting is essential because multiple epochs from the same
participant are strongly related. Splitting individual epochs across training
and validation/test sets would allow participant-specific connectivity
properties to leak into evaluation.

The model partition is used for machine-learning evaluation only. It is not a
biological category and is displayed in plots only to show the origin of each
participant's embedding values.

### Self-supervised GATv2 autoencoder

The model is a graph-attention autoencoder based on GATv2 layers.

The autoencoder is trained without sleep-stage labels as direct prediction
targets. Instead, it learns to reconstruct the observed graph edge structure:

```text
Input:
    Thresholded delta-band connectivity graph

Encoder:
    GATv2 message-passing layers produce node-level latent representations

Graph embedding:
    Node-level latent representations are pooled into one fixed-length
    embedding per 30-second graph

Decoder:
    Estimates whether an edge should be present between pairs of nodes

Training objective:
    Distinguish retained positive edges from sampled non-edges
```

The model is trained with binary edge-reconstruction loss. It receives retained
graph edges as positive examples and sampled non-edges as negative examples.

A successful model should assign higher predicted probabilities to retained
positive edges than to sampled non-edges. Training and validation reconstruction
losses, together with positive-versus-negative edge probability separation, are
tracked to assess training quality and select the best validation epoch.

### Training quality control

The notebook produces a training-history figure with three panels:

| Panel | What it evaluates |
|---|---|
| Reconstruction loss | Training and validation binary edge-reconstruction loss across training epochs |
| Positive/negative edge probability | Whether retained graph edges receive higher predicted probabilities than sampled non-edges |
| Probability gap | Mean retained-edge probability minus mean non-edge probability |

The best model state is selected using the minimum validation reconstruction
loss.

The training-history figure is a model-quality-control result. It verifies that
the autoencoder learns meaningful graph structure, but it is not itself a
sleep-stage finding.

### Graph embeddings

After training, the encoder produces one fixed-length embedding for each input
30-second connectivity graph.

```text
One graph at time t
        ↓
Encoder
        ↓
One low-dimensional graph embedding at time t
```

The distance between two graph embeddings measures how much the learned graph
representation changes between the corresponding connectivity epochs.

```text
embedding distance =
    Euclidean distance between two graph embeddings
```

Interpretation:

```text
Smaller embedding distance:
    Consecutive graphs have more similar learned network representations.
    The within-stage network state is more temporally stable.

Larger embedding distance:
    Consecutive graphs differ more in the learned network representation.
    The within-stage network state changes more from one 30-second epoch
    to the next.
```

Embedding distance is a learned representation-space measure. It should not be
interpreted as a direct physiological distance or as a replacement for
classical graph-theoretical metrics.

### Within-stage stability analysis

The stability analysis is restricted to adjacent graph pairs with the same
scored sleep-stage label.

```text
Included stable pair:
    Wake → Wake
    N1   → N1
    N2   → N2
    N3   → N3
    REM  → REM

Excluded transition pair:
    Wake → N1
    N1   → N2
    N2   → N3
    N3   → REM
    or any other stage transition
```

Only exact adjacent 30-second pairs are included:

```text
stored epoch gap = 1
time gap = 0.5 minutes
```

For every valid stable adjacent pair:

```text
embedding_distance =
    distance between the graph embedding at time t
    and the graph embedding at time t - 1
```

The analysis then summarizes pair-level distances within each participant and
stage:

```text
participant-stage stability value =
    median embedding distance across that participant's valid stable
    adjacent 30-second pairs in the stage
```

The median is used because it is less sensitive to occasional large
epoch-to-epoch graph changes than a simple mean.

This participant-level summary is the unit used for cohort-level descriptive
statistics and exploratory tests. Raw counts of stable 30-second pairs are not
treated as independent observations.

### Exploratory results

The notebook produces participant-level and cohort-level descriptive summaries
of within-stage graph-embedding stability.

The analysis includes:

- Participant-level median embedding-distance profiles across stages.
- A participant-by-stage heatmap of median embedding distances.
- Cohort summaries that give every participant equal weight within a stage.
- Bootstrap 95% confidence intervals around the mean and median of
  participant-level values.
- A small held-out-test-participant summary for independent descriptive
  inspection.

The available descriptive output indicates that the learned delta-band graph
representation has different within-stage stability profiles across the scored
sleep stages. In the current held-out test-participant summary, the median
embedding distances were approximately:

| Stable stage | Held-out participants | Median of participant medians |
|---|---:|---:|
| Wake | 3 | 0.243 |
| N1 | 3 | 0.277 |
| N2 | 3 | 0.275 |
| N3 | 3 | 0.294 |
| REM | 3 | 0.296 |

Lower values correspond to greater temporal stability of the learned graph
representation. The held-out summary is intentionally descriptive because it
contains only three participants.    

### Exploratory permutation check

An exploratory repeated-measures Friedman analysis tests whether
participant-level within-stage stability differs across the five sleep stages.

The observed Friedman chi-square statistic was: 20.714.

A permutation null distribution is generated by repeatedly shuffling stage
labels within participants and recalculating the Friedman statistic. This
preserves the participant-level repeated-measures structure while breaking the
association between a participant's stability values and their stage labels.

The observed statistic lies beyond the visible permutation null distribution,
supporting an exploratory stage-related difference in within-stage embedding
stability. The result indicates that at least one stage differs from at
least one other stage; it does not identify the responsible stage pair(s) or
their direction.   

> The exact permutation p-value should be reported together with the number of
> permutations used. Use the finite-sample corrected estimate:
>
> ```text
> permutation p-value =
> (number of permuted statistics greater than or equal to the observed statistic + 1)
> divided by
> (number of permutations + 1)
> ```
>
> Do not report a finite permutation result as `p = 0`.

### Outputs

Generated Phase B results are saved locally and excluded from Git version
control. These may include:

```text
phase_b_gatv2_autoencoder/
├── model checkpoints
├── training history
├── graph embeddings
├── edge-reconstruction metrics
├── adjacent-pair embedding-distance table
├── participant-level within-stage stability table
├── participant-by-stage stability matrix
├── bootstrap summary table
├── held-out test-subject summary
└── exploratory stability figures
```

Typical saved outputs include:

```text
<band>_autoencoder_training_history.png
stable_adjacent_30second_embedding_pairs.csv
subject_level_within_stage_stability.csv
participant_by_stage_stability_matrix.csv
within_stage_stability_descriptive_summary.csv
within_stage_stability_bootstrap_summary.csv
heldout_test_within_stage_stability_summary.csv
within_stage_stability_combined.png
```

The code and notebook required to recreate these outputs are included in the
repository. Large learned embeddings, model checkpoints, generated tables, and
training figures remain local data products and should not be committed to the
main repository.

### Interpretation and limitations

This is an exploratory proof-of-concept analysis.

- Embedding distances describe differences in a learned model representation;
  they are not direct measures of physiological distance.
- The autoencoder is self-supervised and does not use sleep stage as a training
  label, but learned representations may still reflect subject, recording,
  preprocessing, graph-construction, or signal-quality characteristics.
- Participant-level summaries are used to avoid treating repeated 30-second
  pairs as independent observations.
- Only pairs with the same stage label are included in the primary stability
  analysis; stage transitions are deliberately excluded from the outcome.
- A Friedman test is an general comparison. Pairwise within-participant
  follow-up tests, effect sizes, confidence intervals, and correction for
  multiple comparisons are needed to identify which stages differ.
- The held-out test subset provides independent descriptive context but is too
  small for strong standalone inference.
- Results should be interpreted together with the earlier EEG preprocessing,
  artifact-sensitivity, and graph-construction checks.

### Phase B conclusion

Phase B demonstrates a self-supervised graph-learning workflow for representing
individual 30-second delta-band EEG functional-connectivity graphs and
quantifying their short-timescale stability within conventionally scored sleep
stages.

The analysis provides exploratory evidence that the learned graph
representation does not have identical within-stage stability across all sleep
stages. This complements the stage-mean static network analysis by adding a
temporal, epoch-resolved description of functional-connectivity dynamics.

Further work should focus on preregistered or clearly specified follow-up
comparisons, participant-aware effect-size estimation, robustness to graph
density and model hyperparameters, comparison with non-learned graph-distance
measures, and interpretation of which connectivity features drive learned
embedding differences.

---

## Overall conclusions

This project establishes an end-to-end workflow for studying sleep EEG
functional connectivity as a sequence of weighted brain-network graphs.

```text
Raw overnight EEG and annotations
        ↓
Artifact-aware preprocessing
        ↓
Per-epoch wPLI connectivity graphs
        ↓
Stage-mean graph-theoretical analysis
        +
Epoch-resolved GATv2 graph-embedding analysis
        ↓
Complementary static and temporal views of sleep-network organization
```

The workflow includes:

- Memory-efficient processing of overnight high-density EEG recordings.
- Frequency-specific wPLI graphs for 83 scalp EEG channels in delta, theta,
  alpha, sigma, and beta bands.
- Subject-level HDF5 outputs containing graph matrices, sleep-stage labels,
  original epoch indices, channel names, and interpolation metadata.
- Quality-control and sensitivity checks for montage alignment, HDF5 outputs,
  interpolation, peripheral-channel effects, and frontal-channel effects.
- Stage-mean functional networks and classical weighted graph metrics.
- Threshold robustness checks for static network measures.
- A self-supervised GATv2 autoencoder for learning representations of
  individual 30-second delta-band graphs.
- Participant-level analysis of within-stage temporal stability using
  distances between adjacent graph embeddings.

Together, these components provide a reusable framework for studying sleep EEG
networks at both stage-mean and 30-second time scales.

### Static network analysis

The static analysis created one stage-mean connectivity graph for every
participant, sleep stage, and frequency band.

All 700 participant × stage × band graphs remained fully connected at both 85%
and 90% retained graph density. Global efficiency and characteristic path
length were almost unchanged between these two connected thresholds, supporting
the robustness of the descriptive integration patterns to this threshold choice.

The descriptive results suggest that:

- Sigma-band networks are relatively more integrated during N2 and N3 and less
  integrated during REM.
- Delta-band integration is relatively higher during Wake and N2 and lower
  during N1 and REM.
- Alpha-band integration may be lower during REM.
- Theta-band integration appears comparatively stable across stages.
- Beta-band integration may decrease toward N2/N3 and partly recover in REM,
  although participant-to-participant variation is substantial.

Beta-band results should be interpreted cautiously. A peripheral-channel
sensitivity analysis motivated a narrow robustness check that excludes only
peripheral–peripheral beta edges while keeping all 83 EEG nodes.

### Within-stage graph dynamics

The GATv2 autoencoder analysis adds an epoch-resolved perspective by learning
a graph embedding for each selected 30-second delta-band connectivity graph.

```text
Lower embedding distance:
    More similar learned graph representations between adjacent epochs.

Higher embedding distance:
    Greater change in the learned graph representation between adjacent epochs.
```

Within-stage stability was calculated only from adjacent 30-second epoch pairs
with the same scored stage label. Distances were summarized within each
participant and stage before cohort-level descriptions, so individual epoch
pairs were not treated as independent observations.

The results show substantial participant-specific variation in within-stage
delta-band graph stability. Wake had the lowest descriptive cohort estimate and
REM the highest, but differences between stages were small relative to
individual variation and bootstrap confidence intervals overlapped.

The current analysis therefore supports graph embeddings as a useful
descriptive tool for examining short-timescale sleep-network dynamics. It does
not yet support a strong cohort-wide conclusion that one sleep stage is always
more stable than another.

An exploratory Friedman permutation analysis produced an observed statistic of
20.714 outside the visible permutation null distribution. This suggests that
participant-level stability profiles are not fully identical across stages.
However, this general result does not identify which stage pairs differ and
should be followed by paired, participant-aware comparisons.

### Interpretation

The two analysis branches answer complementary questions:

```text
Static graph metrics:
    How does average network organization differ between stages?

Graph embeddings:
    How much does network structure change from one 30-second epoch to the
    next within the same stage?
```

A sleep stage can have a characteristic average connectivity pattern while also
showing meaningful short-timescale variation within individual participants.

All results remain exploratory and should not be treated as definitive
physiological or clinical conclusions without confirmatory analysis,
robustness checks, and replication.

---

## Future work

### Statistical analysis

- Perform paired, participant-aware follow-up comparisons between sleep stages.
- Report effect sizes and confidence intervals alongside p-values.
- Apply multiple-comparison correction across graph metrics, frequency bands,
  and stage comparisons.
- Define clear inclusion rules for participants with incomplete stage coverage.
- Use repeated-measures models or participant-level permutation tests.

### Method extensions

- Complete and evaluate minimum-spanning-tree and graph-entropy analyses.
- Add node-level, regional, hub, and community-consensus analyses.
- Compare learned embedding distances with classical graph metrics and
  non-learned graph-distance measures.
- Analyse stage transitions separately from stable within-stage pairs.
- Study graph-embedding trajectories across the full night.
- Use explainability methods to investigate which connections contribute most
  to learned graph representations.

## Final perspective

The main contribution of this project is the creation of an auditable workflow
linking high-density sleep EEG preprocessing, connectivity analysis, classical
network science, and self-supervised graph learning.

It provides a practical foundation for studying both average sleep-stage
network organization and short-timescale changes within sleep stages.
