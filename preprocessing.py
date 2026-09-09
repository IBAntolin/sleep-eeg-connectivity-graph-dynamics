"""
Phase 1 -- Data Loading & Preprocessing (ANPHY-Sleep)

Pipeline: stream one epoch at a time from the EDF (never the full night
into memory) -> label with sleep stage -> drop epochs with more than
MAX_BAD_FRACTION of channels flagged bad; downsample 1000 Hz -> 200 Hz
(downsample_signal, below) -> regress chin EMG out of every EEG channel,
per epoch (regress_out_emg, "Step 3b" below) -> interpolate bad channels
within the MAX_BAD_FRACTION threshold, from already EMG-clean reference
channels (see "Step 3" below) -> immediately build per-band wPLI graphs
(connectivity.py), discarding the signal right away -> accumulate the
small graphs in memory for the whole subject -> write everything to one
output_dir/<subject_id>.h5 file (graphs + epoch labels + interpolation
audit trail together, see process_subject()). Filtering-wise, the pipeline
applies: chin-EMG regression (per-epoch OLS against the chin EMG
channel(s) -- chosen over true sample-by-sample LMS, since muscle-tone
contamination changes on a per-epoch, not per-sample, timescale, and OLS
is fully vectorizable), bad-channel interpolation (spherical spline, via
MNE -- chosen over full-signal ICA reconstruction specifically because it
doesn't touch the phase relationships of the *good* channels, which wPLI
depends on), and the anti-aliasing filter that scipy.signal.resample_poly
applies internally as part of downsampling -- a different problem (safe
sample-rate reduction) from the per-band Butterworth filters in
connectivity.py (delta/theta/alpha/sigma/beta isolation for wPLI), so
none of these overlap or duplicate each other's work.

Confirmed from real files:
  - .pos montage: index/name/x/y/z rows, cm units
  - artifact matrix: *_artndxn.mat, MATLAB v7.3/HDF5, (n_epochs, n_channels)
    binary, 1=clean
  - scoring TXT: tab/whitespace-delimited, no header,
    [stage_raw, start_time_sec, duration_sec]; codes W/N1/N2/N3/R/L
    (L = light on/off marker, kept as its own label, not folded into
    Wake or dropped automatically -- filter it out downstream if you
    don't want it as a modeled stage)
  - EDF: 93 channels, all at 1000 Hz. 83 are EEG (matching the .pos
    file), the other 10 are PSG: ChEMG1, ChEMG2, RLEG-, RLEG+, LLEG-,
    LLEG+, ECG1, EOG1, ECG2, EOG2. A full night at 83 channels x 1000 Hz
    is ~18 GB in float64 -- too large to preload, hence the streaming
    design below (mne.io.read_raw_edf tried to preload anyway even with
    preload=False and crashed; pyedflib's readSignal(chn, start, n)
    reads only the requested slice, which is what we use instead).

Requires: numpy, pandas, scipy, h5py (artifact .mat), pyedflib (EDF),
mne (bad-channel interpolation only -- used on small in-memory chunks,
not for loading the EDF itself).
"""

import gc
import re
from fractions import Fraction
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyedflib
from scipy.signal import resample_poly

from connectivity import build_epoch_graphs

# ---------------------------------------------------------------------------
# Config (move to config.py once this stabilizes)
# ---------------------------------------------------------------------------

STAGE_MAP = {
    "W": "Wake",
    "L": "light_on_off",  # experiment light on/off marker, not a sleep stage
    "N1": "N1",
    "N2": "N2",
    "N3": "N3",
    "R": "REM",
}

# Confirmed non-EEG channel names present in the EDF alongside the 83 EEG
# channels -- everything else in the file is treated as EEG.
NON_EEG_CHANNELS = {
    "ChEMG1", "ChEMG2", "RLEG-", "RLEG+", "LLEG-", "LLEG+",
    "ECG1", "EOG1", "ECG2", "EOG2",
}

# Chin EMG channels used as the regressor in regress_out_emg() to remove
# muscle-tone contamination from the EEG before interpolation/graph
# building. Both are used together as a multi-regressor design (see
# get_emg_channel_indices() / regress_out_emg()).
EMG_CHANNEL_NAMES = ["ChEMG1", "ChEMG2"]

EPOCH_LEN_SEC = 30

# Downsample target: raw EDF is 1000 Hz, far more than delta-beta (0.5-30 Hz)
# needs. 200 Hz keeps a conservative >3x Nyquist margin above beta's 30 Hz
# edge while cutting per-epoch sample count (and downstream filtfilt/wPLI
# compute cost) by 5x. See downsample_signal() below.
TARGET_SFREQ = 200.0

# Bad-channel interpolation threshold: an epoch with more than this
# fraction of channels flagged bad is dropped entirely (interpolating too
# many channels from too few good reference points stops being a
# trustworthy reconstruction); at or below it, the bad channels are
# spherical-spline interpolated (via MNE) from the good ones and the
# epoch is kept. See interpolate_bad_channels() / iter_clean_epochs().
MAX_BAD_FRACTION = 0.15

# Subject 08 was removed from the dataset -- skip it everywhere.
SKIP_SUBJECT_NUMBERS = {8}
N_TOTAL_SUBJECTS = 29


def list_subject_ids(skip: set[int] = SKIP_SUBJECT_NUMBERS,
                      n_total: int = N_TOTAL_SUBJECTS) -> list[str]:
    """EPCTL01, EPCTL02, ... EPCTL29, excluding any subject numbers in `skip`."""
    return [f"EPCTL{i:02d}" for i in range(1, n_total + 1) if i not in skip]


# ---------------------------------------------------------------------------
# Step 0 -- discovery (run this first, on one subject, before trusting the rest)
# ---------------------------------------------------------------------------

def discover_subject_files(subject_dir: str | Path) -> dict:
    """
    List everything in a subject's folder and print a quick preview of each
    file type, so we can confirm real filenames/formats before finalizing
    the loaders below.
    """
    subject_dir = Path(subject_dir)
    found = {"edf": [], "scoring_txt": [], "artifact": [], "other": []}

    for f in sorted(subject_dir.rglob("*")):
        if f.is_dir():
            continue
        suffix = f.suffix.lower()
        if suffix == ".edf":
            found["edf"].append(f)
        elif suffix == ".txt" and "scor" in f.stem.lower():
            found["scoring_txt"].append(f)
        elif "artifact" in f.stem.lower() or "artifact" in str(f.parent).lower():
            found["artifact"].append(f)
        else:
            found["other"].append(f)

    print(f"--- {subject_dir} ---")
    for category, files in found.items():
        print(f"\n[{category}] ({len(files)} files)")
        for f in files[:5]:
            print(f"  {f.name}  ({f.stat().st_size / 1024:.1f} KB)")
            if f.suffix.lower() in (".txt", ".csv"):
                try:
                    with open(f, "r", errors="replace") as fh:
                        preview = [next(fh) for _ in range(3)]
                    print("    preview:", [line.strip() for line in preview])
                except (StopIteration, OSError):
                    pass

    return found


# ---------------------------------------------------------------------------
# Step 1 -- electrode positions
# ---------------------------------------------------------------------------

def parse_pos_file(pos_path: str | Path) -> dict[str, np.ndarray]:
    """
    Parse the ANPHY-Sleep .pos file: first line = channel count, then one
    row per channel: index<TAB>channel_name<TAB>x<TAB>y<TAB>z. Coordinates
    are in cm; converted to meters here (MNE/most downstream tools expect
    meters for head coordinates).
    """
    pos_path = Path(pos_path)
    ch_pos = {}
    with open(pos_path, "r") as f:
        n_channels = int(f.readline().strip())
        for _ in range(n_channels):
            parts = f.readline().split()
            _, name, x, y, z = parts
            ch_pos[name] = np.array([float(x), float(y), float(z)]) / 100.0  # cm -> m
    return ch_pos


# ---------------------------------------------------------------------------
# Step 2 -- sleep scoring + artifact matrix -> per-epoch labels
# ---------------------------------------------------------------------------

def load_scoring(scoring_path: str | Path) -> pd.DataFrame:
    """
    Parse the sleep scoring TXT. Confirmed format: whitespace/tab-delimited,
    no header, three columns: stage_raw, start_time (sec), duration_sec.
    """
    df = pd.read_csv(
        scoring_path,
        sep=r"\s+",
        header=None,
        names=["stage_raw", "start_time", "duration_sec"],
    )
    df["stage"] = df["stage_raw"].map(STAGE_MAP)
    return df


def load_artifact_matrix(artifact_path: str | Path) -> np.ndarray:
    """
    Load the artifact matrix from the ANPHY-Sleep *_artndxn.mat file.
    Confirmed format: MATLAB v7.3 (HDF5) containing an 'artndxn' dataset,
    shape (n_epochs, n_channels), binary -- 1 = clean, 0 = artifact, per
    channel per epoch. Returned as a boolean array.
    """
    with h5py.File(artifact_path, "r") as f:
        artndxn = np.array(f["artndxn"])
    return artndxn.astype(bool)


def epoch_clean_mask(artifact_matrix: np.ndarray, mode: str = "all") -> np.ndarray:
    """
    Reduce a (n_epochs, n_channels) clean/artifact matrix to a single
    per-epoch boolean mask.

    mode="all": epoch kept only if every channel is clean (strict; safest
                for network/connectivity analyses where one noisy channel
                can distort many edges).
    mode="majority": epoch kept if more than half the channels are clean
                      (looser -- use only if "all" drops too many epochs).
    """
    if mode == "all":
        return artifact_matrix.all(axis=1)
    if mode == "majority":
        return artifact_matrix.mean(axis=1) > 0.5
    raise ValueError(f"Unknown mode: {mode!r}")


# ---------------------------------------------------------------------------
# Step 3 -- bad-channel interpolation
# ---------------------------------------------------------------------------
#
# The strict "drop the whole epoch if any channel is dirty" policy
# (epoch_clean_mask(mode="all")) was found to drop very different
# fractions of the night across subjects (as low as ~8% for some, as high
# as ~50% for others) -- since it only takes one bad channel in 83 to lose
# an entire epoch. Interpolating a small number of bad channels per epoch
# (spherical spline, via MNE, using electrode positions from the .pos
# file) recovers most of those epochs without touching the phase
# relationships of the other, actually-clean channels -- important since
# wPLI depends specifically on inter-channel phase, which full-signal
# methods like ICA reconstruction can risk subtly altering. Epochs with
# more than MAX_BAD_FRACTION of channels bad are still dropped outright --
# interpolating too many channels from too few good reference points
# stops being a trustworthy reconstruction.
#
# NOTE / assumption: this assumes the artifact matrix's channel columns
# are in the same order as the .pos file / EEG channel list. This wasn't
# independently cross-checked against real column headers (the .mat file
# has no channel-name metadata, just a plain (n_epochs, n_channels)
# matrix) -- if you want to double check, compare each channel's overall
# bad-fraction (artifact_matrix.mean(axis=0)) against known-noisy
# electrode positions (e.g. temporal/mastoid channels are typically
# noisier than central ones) as a plausibility check.


def normalize_channel_name(name: str) -> str:
    """
    Shared channel-name normalization used to match EDF channel names
    against the .pos file's names despite the two files disagreeing on
    capitalization for at least some channels (e.g. EDF "FZ" vs. .pos
    "Fz") and the EDF carrying reference/average/CAR suffixes the .pos
    file doesn't. Lowercases, strips a trailing -ref/_ref/-avg/_avg/-car/
    _car suffix if present, then strips any leading/trailing -, +, _.
    Used by both build_montage_info() (EDF ch_names -> .pos positions) and
    get_eeg_channel_indices() (EDF channel list -> .pos anomaly check), so
    the two stay consistent instead of drifting apart.
    """
    name = name.lower().strip()
    for suffix in ["-ref", "_ref", "-avg", "_avg", "-car", "_car"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    name = re.sub(r'^[\-\+_]+|[\-\+_]+$', '', name)
    return name


def build_montage_info(ch_names: list[str], pos_path: str | Path, sfreq: float):
    """
    Build an mne.Info object with a montage attached, once per subject --
    reused for every epoch's interpolation rather than rebuilt each time,
    since electrode geometry doesn't change epoch to epoch. Requires mne
    (only used here, on small in-memory chunks -- unrelated to the earlier
    full-night mne.io.read_raw_edf memory problem).

    Matches ch_names (from the EDF) against the .pos file's channel names
    via normalize_channel_name() (case-insensitive, suffix-stripped).
    Raises a clear ValueError (not a bare KeyError) listing every EDF
    channel name that still can't be found in the .pos file at all, even
    after normalization.
    """
    import mne

    ch_pos = parse_pos_file(pos_path)
    ch_pos_norm = {normalize_channel_name(name): pos for name, pos in ch_pos.items()}

    missing = []
    matched_pos = {}

    for name in ch_names:
        key = normalize_channel_name(name)
        if key in ch_pos_norm:
            matched_pos[name] = ch_pos_norm[key]
        else:
            missing.append(name)

    if missing:
        raise ValueError(
            f"{len(missing)} channel(s) from the EDF have no match in the "
            f".pos file (even after normalization): {missing}.\n"
            f".pos file channels: {sorted(ch_pos.keys())}"
        )

    montage = mne.channels.make_dig_montage(
        ch_pos=matched_pos, coord_frame="head"
    )
    info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types="eeg")
    info.set_montage(montage)
    return info


def interpolate_bad_channels(signal: np.ndarray, bad_channel_mask: np.ndarray,
                              montage_info) -> np.ndarray:
    """
    Spherical-spline interpolate the channels flagged True in
    bad_channel_mask, reconstructing them from the good channels via MNE.
    signal: (n_channels, n_samples), same channel order as montage_info's
    ch_names. Returns a corrected (n_channels, n_samples) float32 array;
    if bad_channel_mask is all False, returns signal unchanged (no MNE
    call, to avoid the overhead on the majority of epochs that have zero
    bad channels).
    """
    if not bad_channel_mask.any():
        return signal

    import mne

    raw = mne.io.RawArray(signal.astype(np.float64), montage_info, verbose="error")
    raw.info["bads"] = [name for name, bad in zip(montage_info["ch_names"], bad_channel_mask) if bad]
    raw.interpolate_bads(reset_bads=True, verbose="error")
    return raw.get_data().astype(np.float32)


# ---------------------------------------------------------------------------
# Step 4 -- stream epochs directly from the EDF (this replaces loading the
# whole raw recording into memory)
# ---------------------------------------------------------------------------

def get_eeg_channel_indices(
    edf_labels: list[str],
    pos_path: str | Path | None = None,
) -> tuple[list[int], list[str]]:
    """
    Split the EDF's channel list into EEG channel (index, name) pairs.

    Base filter: exclude anything in NON_EEG_CHANNELS (the 10 known PSG
    channels). If pos_path is given, every surviving channel is also
    cross-checked against the .pos file's channel names (via
    normalize_channel_name(), same matching used by build_montage_info):

    - Dropped: a channel whose normalized name has no match at all in the
      .pos file. Printed so the anomaly is visible, not silent.

    - Collision: normalize_channel_name() strips trailing -/+/_ and
      -ref/-avg/-car suffixes, so an anomalous duplicate channel with a
      trailing character (e.g. EPCTL06's stray "C1-" alongside the
      genuine "C1") normalizes to the *same* target as the real channel
      instead of failing to match -- so a plain "was it found" check
      would silently accept both, letting the anomalous duplicate through
      as a second, spurious "C1". Resolved by keeping only the exact
      match (raw name, lowercased and stripped of whitespace only, equal
      to the .pos name) and dropping the other(s); if no candidate is an
      unambiguous exact match, raises ValueError rather than guessing.

    If pos_path is not given, only the NON_EEG_CHANNELS filter applies
    (no cross-check, no collision detection) -- fine for a quick
    discovery pass, but pos_path should be supplied whenever the result
    feeds interpolation/graph building downstream.
    """
    eeg = [(i, name) for i, name in enumerate(edf_labels) if name not in NON_EEG_CHANNELS]

    if pos_path is not None:
        pos_names_norm = {normalize_channel_name(n) for n in parse_pos_file(pos_path)}
        matched = [(i, name) for i, name in eeg if normalize_channel_name(name) in pos_names_norm]
        dropped = [name for i, name in eeg if normalize_channel_name(name) not in pos_names_norm]
        if dropped:
            print(f"get_eeg_channel_indices: dropping {len(dropped)} channel(s) "
                  f"not found in .pos file (anomalous, not genuine EEG): {dropped}")

        by_target: dict[str, list[tuple[int, str]]] = {}
        for i, name in matched:
            by_target.setdefault(normalize_channel_name(name), []).append((i, name))

        resolved = []
        unresolved_collisions = []
        for target, candidates in by_target.items():
            if len(candidates) == 1:
                resolved.append(candidates[0])
                continue
            exact = [(i, name) for i, name in candidates if name.lower().strip() == target]
            if len(exact) == 1:
                resolved.append(exact[0])
                dupes = [name for _, name in candidates if name != exact[0][1]]
                print(f"get_eeg_channel_indices: name collision at {target!r} -- "
                      f"kept exact match {exact[0][1]!r}, dropped anomalous "
                      f"duplicate(s) {dupes}")
            else:
                unresolved_collisions.append((target, [name for _, name in candidates]))

        if unresolved_collisions:
            raise ValueError(
                f"Ambiguous channel name collision(s) after normalization -- "
                f"multiple EDF channels map to the same .pos electrode and "
                f"none is an unambiguous exact match: {unresolved_collisions}. "
                f"Resolve manually before proceeding (e.g. inspect the raw "
                f"EDF channel list and decide which one is genuine)."
            )

        eeg = sorted(resolved, key=lambda pair: pair[0])  # restore original EDF order

    indices, names = zip(*eeg) if eeg else ([], [])
    return list(indices), list(names)


def get_emg_channel_indices(
    edf_labels: list[str],
    emg_names: list[str] = EMG_CHANNEL_NAMES,
) -> tuple[list[int], list[str]]:
    """
    Locate the chin EMG reference channel(s) (EMG_CHANNEL_NAMES, default
    ChEMG1 + ChEMG2) within the EDF's full channel list -- these are the
    regressors used in regress_out_emg() to remove muscle-tone
    contamination from the EEG. Raises ValueError if any expected EMG
    channel is missing from the EDF, rather than silently regressing
    against fewer reference channels than intended.
    """
    found = [(i, name) for i, name in enumerate(edf_labels) if name in emg_names]
    missing = sorted(set(emg_names) - {name for _, name in found})
    if missing:
        raise ValueError(f"Expected EMG channel(s) not found in EDF: {missing}")
    indices, names = zip(*found)
    return list(indices), list(names)


def downsample_signal(
    signal: np.ndarray,
    orig_sfreq: float,
    target_sfreq: float = TARGET_SFREQ,
) -> tuple[np.ndarray, float]:
    """
    Downsample a (n_channels, n_samples) signal from orig_sfreq to
    target_sfreq using polyphase filtering (scipy.signal.resample_poly),
    which applies its own anti-aliasing low-pass filter as part of the
    decimation -- unlike naive slicing (signal[:, ::step]), this doesn't
    fold high-frequency content (EMG, movement artifact) back down into
    the delta-beta range as aliasing distortion.

    orig_sfreq/target_sfreq is converted to a small integer up/down ratio
    (e.g. 1000/200 -> up=1, down=5) since resample_poly requires integers;
    limit_denominator guards against float sfreq values (e.g. 999.99999)
    producing a huge, slow ratio.

    Returns (downsampled_signal, actual_new_sfreq) -- actual_new_sfreq is
    computed from the exact up/down ratio applied, which may differ
    infinitesimally from target_sfreq if orig_sfreq/target_sfreq wasn't
    an exact ratio; always use the returned value as the epoch's sfreq
    from here on, not target_sfreq itself.
    """
    ratio = Fraction(target_sfreq / orig_sfreq).limit_denominator(1000)
    up, down = ratio.numerator, ratio.denominator

    downsampled = resample_poly(signal, up, down, axis=-1).astype(np.float32)
    actual_new_sfreq = orig_sfreq * up / down
    return downsampled, actual_new_sfreq


# ---------------------------------------------------------------------------
# Step 3b -- per-epoch chin-EMG regression (runs before interpolation, see
# iter_clean_epochs -- interpolation reconstructs bad channels from good
# ones, so those good channels need to already be EMG-clean or the
# interpolated channel would inherit the contamination)
# ---------------------------------------------------------------------------

def regress_out_emg(eeg_signal: np.ndarray, emg_signal: np.ndarray) -> np.ndarray:
    """
    Remove EMG-linked variance from every EEG channel via ordinary least
    squares regression against the chin EMG reference channel(s), fit
    fresh per epoch rather than once for the whole night -- this is the
    "per-epoch regression" alternative to true sample-by-sample LMS:
    it lets the correction strength adapt as muscle tone changes across
    the night (the thing LMS would track), at the timescale that
    contamination actually changes on (epoch-to-epoch, not
    sample-to-sample), using a method that's fully vectorizable and cheap.

    eeg_signal: (n_eeg_channels, n_samples)
    emg_signal: (n_emg_channels, n_samples) -- same epoch, same sfreq,
    already downsampled consistently with eeg_signal.

    All n_eeg_channels are fit in a single np.linalg.lstsq call (one
    shared design matrix -- the EMG reference channel(s) plus an
    intercept column -- against all EEG channels as multiple right-hand
    sides at once), not a per-channel Python loop.

    Only the EMG-driven component of each channel's fit is subtracted;
    each channel's own intercept (its DC level within the epoch) is left
    untouched, since that offset isn't part of the EMG contamination.

    Returns a corrected (n_eeg_channels, n_samples) float32 array, same
    shape/dtype as eeg_signal.
    """
    n_samples = eeg_signal.shape[1]
    design = np.vstack([emg_signal.astype(np.float64), np.ones(n_samples)]).T  # (n_samples, n_emg+1)

    coeffs, *_ = np.linalg.lstsq(design, eeg_signal.T.astype(np.float64), rcond=None)  # (n_emg+1, n_eeg_channels)
    fitted = design @ coeffs  # (n_samples, n_eeg_channels)

    intercept = coeffs[-1, :]  # (n_eeg_channels,) -- kept, not subtracted
    emg_component = (fitted - intercept).T  # (n_eeg_channels, n_samples)

    corrected = eeg_signal.astype(np.float64) - emg_component
    return corrected.astype(np.float32)


def iter_clean_epochs(
    edf_path: str | Path,
    scoring_df: pd.DataFrame,
    artifact_matrix: np.ndarray | None = None,
    epoch_len: int = EPOCH_LEN_SEC,
    stages_to_drop: tuple[str, ...] = ("light_on_off",),
    target_sfreq: float | None = TARGET_SFREQ,
    max_epochs: int | None = None,
    pos_path: str | Path | None = None,
    max_bad_fraction: float = MAX_BAD_FRACTION,
    emg_correct: bool = True,
    emg_names: list[str] = EMG_CHANNEL_NAMES,
):
    """
    Generator: yields one (epoch_index, stage, channel_names, signal,
    sfreq, n_bad_channels) tuple at a time, reading only that epoch's
    samples from disk via pyedflib -- the full night is never loaded into
    memory. signal has shape (n_eeg_channels, n_samples_per_epoch), dtype
    float32, unfiltered by any band filter (band filtering happens
    per-band in connectivity.build_epoch_graphs).

    Per-epoch order of operations: read EEG (+ EMG, if emg_correct) ->
    downsample both -> regress out EMG (regress_out_emg(), if emg_correct)
    -> interpolate any bad channels (interpolate_bad_channels(), if any
    are flagged) -> yield. EMG regression runs before interpolation
    deliberately: interpolation reconstructs bad channels from the good
    ones, so those good channels need to already be EMG-clean or the
    interpolated channel would inherit the contamination.

    artifact_matrix: the raw (n_epochs, n_channels) boolean array from
    load_artifact_matrix() (True = clean), not a pre-reduced per-epoch
    mask. Per epoch: if the fraction of bad channels exceeds
    max_bad_fraction, the epoch is dropped entirely; otherwise, any bad
    channels present are interpolated (via interpolate_bad_channels(),
    requires pos_path) and the epoch is kept. Pass artifact_matrix=None to
    skip artifact-based filtering/interpolation entirely.

    pos_path: required whenever artifact_matrix is given, to build the
    montage used for interpolation (see build_montage_info()) -- raises
    ValueError if artifact_matrix is provided without it. Also passed
    through to get_eeg_channel_indices() (whenever given) so anomalous
    non-EEG channels absent from the .pos file get excluded up front --
    see get_eeg_channel_indices().

    emg_correct: if True (default), regress the chin EMG channel(s)
    (emg_names, default EMG_CHANNEL_NAMES) out of every EEG channel each
    epoch via regress_out_emg() -- see that function's docstring for why
    per-epoch OLS regression was chosen over true sample-by-sample LMS.
    Set False to skip EMG correction entirely (e.g. for debugging, or
    comparing corrected vs. uncorrected graphs).

    If target_sfreq is not None (default: module-level TARGET_SFREQ,
    200 Hz), each epoch (EEG and, if emg_correct, EMG alike) is
    downsampled via downsample_signal() right after reading it from the
    EDF and before yielding -- so every downstream consumer
    (connectivity.build_epoch_graphs) works on the downsampled signal and
    the correct (downsampled) sfreq automatically, and EMG regression
    operates on EEG/EMG at matching sample rates. Pass target_sfreq=None
    to skip downsampling and yield the raw EDF signal (e.g. 1000 Hz)
    unchanged, if ever needed for debugging.

    max_epochs: if set, stop after this many *clean* (kept) epochs have
    been yielded -- for a quick smoke test on a handful of epochs before
    committing to a full ~900-epoch subject.

    Epochs with stage in `stages_to_drop` are also skipped. Rows in
    scoring_df with no matching artifact_matrix entry are skipped too
    (length mismatch safety).
    """
    if artifact_matrix is not None and pos_path is None:
        raise ValueError(
            "pos_path is required when artifact_matrix is given, to build "
            "the montage used for bad-channel interpolation."
        )

    reader = pyedflib.EdfReader(str(edf_path))
    try:
        edf_labels = reader.getSignalLabels()
        eeg_indices, eeg_names = get_eeg_channel_indices(edf_labels, pos_path=pos_path)
        sfreq = reader.getSampleFrequency(eeg_indices[0])
        n_samples_per_epoch = int(epoch_len * sfreq)

        emg_indices = []
        if emg_correct:
            emg_indices, _ = get_emg_channel_indices(edf_labels, emg_names=emg_names)

        montage_info = None  # built lazily below, once, at the target sfreq

        bad = None
        bad_fraction = None
        if artifact_matrix is not None:
            bad = ~artifact_matrix  # True = artifact/bad
            bad_fraction = bad.mean(axis=1)

        n_yielded = 0
        for i, row in scoring_df.iterrows():
            if max_epochs is not None and n_yielded >= max_epochs:
                break

            if row["stage"] in stages_to_drop or pd.isna(row["stage"]):
                continue

            n_bad_channels = 0
            if bad is not None:
                if i >= len(bad):
                    continue
                if bad_fraction[i] > max_bad_fraction:
                    continue  # too many bad channels to trust interpolation -- drop
                n_bad_channels = int(bad[i].sum())

            start_sample = int(row["start_time"] * sfreq)
            signal = np.stack([
                reader.readSignal(ch, start=start_sample, n=n_samples_per_epoch)
                for ch in eeg_indices
            ]).astype(np.float32)

            emg_signal = None
            if emg_correct:
                emg_signal = np.stack([
                    reader.readSignal(ch, start=start_sample, n=n_samples_per_epoch)
                    for ch in emg_indices
                ]).astype(np.float32)

            epoch_sfreq = sfreq
            if target_sfreq is not None:
                signal, epoch_sfreq = downsample_signal(signal, sfreq, target_sfreq)
                if emg_correct:
                    emg_signal, _ = downsample_signal(emg_signal, sfreq, target_sfreq)

            if emg_correct:
                signal = regress_out_emg(signal, emg_signal)
                del emg_signal

            if n_bad_channels > 0:
                if montage_info is None:
                    montage_info = build_montage_info(eeg_names, pos_path, epoch_sfreq)
                signal = interpolate_bad_channels(signal, bad[i], montage_info)

            n_yielded += 1
            yield i, row["stage"], eeg_names, signal, epoch_sfreq, n_bad_channels
    finally:
        reader.close()


# ---------------------------------------------------------------------------
# Step 5 -- full per-subject pipeline: stream, filter, build graphs, save
# only the small graphs (merged Phase 1 + Phase 2 -- the filtered signal
# itself is never written to disk, only the resulting wPLI matrices)
# ---------------------------------------------------------------------------

def process_subject(
    subject_id: str,
    edf_path: str | Path,
    scoring_path: str | Path,
    artifact_path: str | Path | None,
    output_dir: str | Path,
    pos_path: str | Path | None = None,
    max_bad_fraction: float = MAX_BAD_FRACTION,
    max_epochs: int | None = None,
    show_progress: bool = True,
    emg_correct: bool = True,
) -> Path:
    """
    Run Phase 1 + Phase 2 for one subject in a single streaming pass: for
    each epoch, read + downsample it (never more than one epoch's raw
    signal in memory), regress out chin-EMG contamination (if
    emg_correct, see regress_out_emg()), interpolate any bad channels
    within max_bad_fraction (dropping the epoch outright if more than
    that are bad -- see iter_clean_epochs), immediately build its
    per-band wPLI graphs, and accumulate just those small graphs (not the
    signal) in memory for the whole subject -- at ~137 KB/epoch even ~900
    epochs is only ~120 MB, fine to hold at once. Once the subject is done,
    everything is written to a single HDF5 file: output_dir/<subject_id>.h5,
    containing one dataset per band, shape (n_kept_epochs, n_channels,
    n_channels), plus "epoch_idx", "stage", and "n_bad_channels" datasets
    (the last recording how many channels were interpolated in that row,
    0 for epochs that needed no correction) -- no separate CSV needed, the
    metadata travels with the data in the same file.

    pos_path: electrode position file, required if artifact_path is given
    (needed to build the interpolation montage). Not required if
    artifact_path is None (no artifact-based filtering/interpolation at
    all in that case).

    max_bad_fraction: see MAX_BAD_FRACTION / iter_clean_epochs.

    max_epochs: cap the number of kept epochs processed (see
    iter_clean_epochs) -- for a quick smoke test on e.g. 10 epochs before
    committing to a full ~900-epoch subject. None (default) processes the
    whole night.

    show_progress: print a live per-epoch progress bar (via tqdm) showing
    epochs/sec and elapsed time -- useful for checking the pipeline isn't
    stuck or unexpectedly slow before letting it run unattended on more
    subjects.

    emg_correct: if True (default), regress the chin EMG channel(s) out
    of every EEG channel each epoch before interpolation/graph building
    -- see iter_clean_epochs() / regress_out_emg(). Set False to skip.

    One file per subject (instead of one .npz per epoch) keeps the file
    count manageable across all subjects (~28 files instead of ~25,000),
    while gzip compression inside HDF5 keeps size comparable to the old
    per-epoch .npz files. Use read_subject_graphs() to load it back.

    Returns the path to the written .h5 file.
    """
    subject_num = int("".join(filter(str.isdigit, subject_id)))
    if subject_num in SKIP_SUBJECT_NUMBERS:
        raise ValueError(
            f"{subject_id} is in SKIP_SUBJECT_NUMBERS (removed from the "
            f"dataset) -- not processing it."
        )

    scoring_df = load_scoring(scoring_path)
    artifact_matrix = None
    if artifact_path:
        artifact_matrix = load_artifact_matrix(artifact_path)
        if pos_path is None:
            raise ValueError(
                "pos_path is required when artifact_path is given, to build "
                "the montage used for bad-channel interpolation."
            )

    graphs_by_band: dict[str, list[np.ndarray]] = {}
    epoch_indices: list[int] = []
    stages: list[str] = []
    n_bad_channels_list: list[int] = []
    band_names: list[str] | None = None
    channel_names: list[str] | None = None

    epoch_iter = iter_clean_epochs(
        edf_path, scoring_df, artifact_matrix,
        max_epochs=max_epochs, pos_path=pos_path, max_bad_fraction=max_bad_fraction,
        emg_correct=emg_correct,
    )
    if show_progress:
        try:
            from tqdm.auto import tqdm
            epoch_iter = tqdm(
                epoch_iter,
                total=max_epochs,
                desc=f"{subject_id} epochs",
                unit="epoch",
            )
        except ImportError:
            print("tqdm not installed (pip install tqdm) -- running without a progress bar.")

    for epoch_idx, stage, ch_names, signal, sfreq, n_bad_channels in epoch_iter:
        if channel_names is None:
            channel_names = ch_names  # same EEG channel order for every epoch of this subject

        graphs = build_epoch_graphs(signal, sfreq)

        if band_names is None:
            band_names = list(graphs.keys())
            graphs_by_band = {band: [] for band in band_names}

        for band in band_names:
            graphs_by_band[band].append(graphs[band])

        epoch_indices.append(epoch_idx)
        stages.append(stage)
        n_bad_channels_list.append(n_bad_channels)
        del signal, graphs

    gc.collect()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{subject_id}.h5"

    with h5py.File(output_path, "w") as f:
        f.attrs["subject_id"] = subject_id
        f.attrs["n_epochs"] = len(epoch_indices)
        f.attrs["bands"] = band_names or []
        f.attrs["max_bad_fraction"] = max_bad_fraction

        for band in band_names or []:
            stacked = np.stack(graphs_by_band[band]).astype(np.float32)  # (n_epochs, 83, 83)
            f.create_dataset(band, data=stacked, compression="gzip", compression_opts=4)

        f.create_dataset("epoch_idx", data=np.array(epoch_indices, dtype=np.int32))
        f.create_dataset("stage", data=np.array(stages, dtype=h5py.string_dtype()))
        f.create_dataset("n_bad_channels", data=np.array(n_bad_channels_list, dtype=np.int16))
        # row/column i of every band's adjacency matrix corresponds to
        # channel_names[i] -- same order for every epoch of this subject,
        # so a matrix's rows/cols can be mapped to electrode names (and,
        # via parse_pos_file, to 3D positions) without touching the EDF
        # again -- needed both for plotting and as GNN node identity/
        # position features later.
        f.create_dataset("channel_names", data=np.array(channel_names or [], dtype=h5py.string_dtype()))

    return output_path


def read_subject_graphs(h5_path: str | Path) -> dict:
    """
    Load a subject's .h5 file back into memory: returns a dict with keys
    "subject_id" (str), "epoch_idx" (int array), "stage" (str array),
    "n_bad_channels" (int array -- how many channels were interpolated in
    that epoch, 0 if none), "channel_names" (str array, length
    n_channels -- row/col i of every band's matrices corresponds to
    channel_names[i]), and one key per band (float32 array, shape
    (n_epochs, n_channels, n_channels)).
    """
    with h5py.File(h5_path, "r") as f:
        result = {
            "subject_id": f.attrs["subject_id"],
            "epoch_idx": f["epoch_idx"][:],
            "stage": np.array([s.decode() if isinstance(s, bytes) else s
                                for s in f["stage"][:]]),
            "n_bad_channels": f["n_bad_channels"][:] if "n_bad_channels" in f else None,
            "channel_names": np.array([c.decode() if isinstance(c, bytes) else c
                                        for c in f["channel_names"][:]]),
        }
        for band in f.attrs["bands"]:
            result[band] = f[band][:]
    return result


def build_global_manifest(output_dir: str | Path) -> pd.DataFrame:
    """
    Scan output_dir for all *.h5 subject files and build one manifest
    DataFrame (subject_id, epoch_idx, stage, n_bad_channels) across all of
    them -- a lightweight, easy-to-pandas-load index for Phase 3/4,
    without needing to open every subject's HDF5 file just to see what
    stages/interpolation activity it has. n_bad_channels is NaN for older
    files written before that dataset was added.
    """
    output_dir = Path(output_dir)
    rows = []
    for h5_path in sorted(output_dir.glob("*.h5")):
        with h5py.File(h5_path, "r") as f:
            subject_id = f.attrs["subject_id"]
            n_bad = f["n_bad_channels"][:] if "n_bad_channels" in f else None
            for idx, (epoch_idx, stage) in enumerate(zip(f["epoch_idx"][:], f["stage"][:])):
                stage_str = stage.decode() if isinstance(stage, bytes) else stage
                rows.append({
                    "subject_id": subject_id,
                    "epoch_idx": int(epoch_idx),
                    "stage": stage_str,
                    "n_bad_channels": int(n_bad[idx]) if n_bad is not None else None,
                    "filepath": str(h5_path),
                })
    return pd.DataFrame(rows)