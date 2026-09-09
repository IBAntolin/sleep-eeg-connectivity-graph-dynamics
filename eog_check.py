"""
Checks whether frontal, eye-adjacent channels (the ones most prone to
picking up EOG -- eye movement/blink -- contamination) are over-
represented among the strongest wPLI edges in delta/theta, the low-
frequency bands where eye artifact classically shows up (unlike the
EMG/muscle pattern found earlier, which concentrated in beta and involved
temporal/outer-ring channels instead).

Mirrors the logic already validated in peripheral_channel_check.py (2-way
"touches a flagged channel" enrichment check), generalized to a different
channel set and different bands, kept as a separate module rather than
modifying the already-validated peripheral check.

Two levels of check, because they answer different questions:

  run_frontal_check()          -- pooled across all subjects/epochs. Answers
                                   "does the dataset as a whole need EOG
                                   correction?" Good for the pipeline-level
                                   decision, but can hide two things that
                                   matter specifically for EOG:
                                     - a few outlier subjects (poor
                                       electrode contact near the eyes,
                                       unusually blink-heavy nights) whose
                                       signal gets diluted by everyone else
                                     - stage-dependence: eye movements
                                       concentrate in REM, so mixing REM
                                       and NREM epochs into one average can
                                       wash out a real, expected pattern

  run_frontal_check_subject()  -- one subject, broken down by sleep stage.
                                   Answers "does THIS subject show frontal
                                   enrichment, and specifically in the
                                   stages where EOG contamination is
                                   actually expected (REM), rather than
                                   uniformly across the night?" Run this on
                                   individual subjects (e.g. ones that look
                                   borderline in the pooled check, or just
                                   a handful spot-checked at random) before
                                   concluding the pooled "no effect" result
                                   is the full picture.

Usage:
    from eog_check import run_frontal_check, run_frontal_check_subject

    run_frontal_check(Path("graphs_output"), bands=("delta", "theta"))
    run_frontal_check_subject(Path("graphs_output/EPCTL05.h5"),
                               bands=("delta", "theta"))
"""

from pathlib import Path

import h5py
import numpy as np

# Frontopolar / eye-adjacent channels: the classic EOG-prone sites,
# directly above or immediately adjacent to the eyes. Matched
# case-insensitively below, since this dataset's own channel names aren't
# consistently capitalized even within the same list (e.g. "Fp1" vs
# "FPZ" vs "AFZ" all appear as different case styles for related
# channels) -- same lesson learned from every other naming mismatch in
# this project so far.
FRONTAL_CHANNELS = {"FP1", "FP2", "FPZ", "AF7", "AF3", "AFZ", "AF4", "AF8"}

# Stages grouped as REM vs everything else, since eye-movement density is
# the thing that differs -- not each individual NREM substage. Adjust the
# exact stage label strings here if this dataset's scoring_df uses
# different names (check f["stage"][:] once to confirm).
REM_STAGE_NAMES = {"REM", "R", "rem"}


def _is_frontal(name: str) -> bool:
    return name.upper() in FRONTAL_CHANNELS


def _expected_baseline_fraction(ch_names: list[str]) -> float:
    """Chance-level fraction of all edges touching at least one frontal channel."""
    n = len(ch_names)
    is_f = np.array([_is_frontal(name) for name in ch_names])
    n_f, n_other = is_f.sum(), n - is_f.sum()
    total = n * (n - 1) / 2
    other_only = n_other * (n_other - 1) / 2
    return 1 - (other_only / total)


def _top_edge_frontal_fraction(con: np.ndarray, ch_names: list[str], top_percent: float) -> float:
    iu_i, iu_j = np.triu_indices_from(con, k=1)
    vals = con[iu_i, iu_j]
    thresh = np.percentile(vals, 100 - top_percent)
    keep = vals >= thresh

    is_f = np.array([_is_frontal(name) for name in ch_names])
    touches_frontal = is_f[iu_i] | is_f[iu_j]
    kept = touches_frontal[keep]
    return kept.mean() if len(kept) else np.nan


def check_one_subject(h5_path: str | Path, band: str, top_percent: float, epoch_stride: int):
    with h5py.File(h5_path, "r") as f:
        if band not in f.attrs["bands"]:
            return None, None
        ch_names = [c.decode() if isinstance(c, bytes) else c for c in f["channel_names"][:]]
        n_epochs = f.attrs["n_epochs"]

        fractions = []
        for e in range(0, n_epochs, epoch_stride):
            con = f[band][e]
            frac = _top_edge_frontal_fraction(con, ch_names, top_percent)
            if not np.isnan(frac):
                fractions.append(frac)

    baseline = _expected_baseline_fraction(ch_names)
    return (np.mean(fractions) if fractions else np.nan), baseline


def run_frontal_check(output_dir: str | Path, bands: tuple[str, ...] = ("delta", "theta"),
                       top_percent: float = 4, epoch_stride: int = 10) -> None:
    """
    Across every subject's .h5 file, for each band in `bands`: what
    fraction of the top X% strongest edges touch a frontal/eye-adjacent
    channel, vs. the chance baseline. Run this BEFORE deciding whether an
    EOG correction is worth adding -- if delta/theta show no enrichment
    (expected, since eye artifact wasn't hypothesized to concentrate
    there in the same way EMG did in beta), that's real evidence EOG
    correction isn't needed, not just an assumption.

    This is a POOLED check (all subjects, all epochs of a given stage mix,
    averaged together) -- see module docstring for why a per-subject,
    per-stage check (run_frontal_check_subject) is also worth running
    before trusting a "no effect" result here.
    """
    output_dir = Path(output_dir)
    h5_files = sorted(output_dir.glob("*.h5"))
    if not h5_files:
        print(f"No .h5 files found in {output_dir}")
        return

    print(f"Checking {len(h5_files)} subject(s), top {top_percent}% edges, "
          f"every {epoch_stride}th epoch, bands={bands}...\n")

    print(f"{'Band':<8} {'Baseline':>10} {'Observed':>10} {'Enrichment':>12}")
    print("-" * 44)

    for band in bands:
        fracs, baseline = [], None
        for h5_path in h5_files:
            frac, subj_baseline = check_one_subject(h5_path, band, top_percent, epoch_stride)
            if frac is None or np.isnan(frac):
                continue
            baseline = subj_baseline
            fracs.append(frac)

        if not fracs:
            print(f"{band:<8} (no data)")
            continue

        observed = np.mean(fracs)
        enrichment = observed / baseline if baseline > 0 else np.nan
        flag = "  <-- [FLAG]" if enrichment > 1.5 else ""
        print(f"{band:<8} {baseline:>9.1%} {observed:>9.1%} {enrichment:>11.2f}x{flag}")

    print("\nA ratio near 1.0x means frontal channels show up in the strongest edges "
          "about as often as chance predicts -- no evidence of EOG concentration, "
          "supporting leaving EOG uncorrected. A ratio well above 1.0x would be real "
          "evidence worth investigating before deciding.")


def check_one_subject_by_stage(h5_path: str | Path, band: str, top_percent: float,
                                epoch_stride: int):
    """
    Like check_one_subject, but returns two separate fractions: one
    computed only over REM epochs, one over everything else. Returns
    (rem_fraction, non_rem_fraction, baseline) -- any of the three can be
    NaN if that stage has no epochs at all for this subject/band.
    """
    with h5py.File(h5_path, "r") as f:
        if band not in f.attrs["bands"]:
            return None, None, None
        ch_names = [c.decode() if isinstance(c, bytes) else c for c in f["channel_names"][:]]
        n_epochs = f.attrs["n_epochs"]
        stages = np.array([s.decode() if isinstance(s, bytes) else s for s in f["stage"][:]])

        rem_fracs, non_rem_fracs = [], []
        for e in range(0, n_epochs, epoch_stride):
            con = f[band][e]
            frac = _top_edge_frontal_fraction(con, ch_names, top_percent)
            if np.isnan(frac):
                continue
            if stages[e] in REM_STAGE_NAMES:
                rem_fracs.append(frac)
            else:
                non_rem_fracs.append(frac)

    baseline = _expected_baseline_fraction(ch_names)
    rem_mean = np.mean(rem_fracs) if rem_fracs else np.nan
    non_rem_mean = np.mean(non_rem_fracs) if non_rem_fracs else np.nan
    return rem_mean, non_rem_mean, baseline


def run_frontal_check_subject(h5_path: str | Path, bands: tuple[str, ...] = ("delta", "theta"),
                               top_percent: float = 4, epoch_stride: int = 5) -> None:
    """
    Single-subject version of run_frontal_check, broken down by REM vs
    non-REM (see module docstring for why: eye-movement density differs
    sharply by stage, so a subject-wide average can hide a real,
    REM-concentrated effect -- or, just as informative, confirm ITS
    ABSENCE even in REM, which is stronger evidence than an absence in a
    stage-blind average).

    epoch_stride defaults smaller than run_frontal_check's (5 vs 10) since
    a single subject has far fewer epochs to draw from per stage,
    especially REM (typically a minority of the night) -- denser sampling
    keeps the REM-epoch count from being too small to trust.

    A subject with too few REM epochs at this stride will print "n/a" for
    that column rather than a misleading single-epoch estimate.
    """
    h5_path = Path(h5_path)
    with h5py.File(h5_path, "r") as f:
        subject_id = f.attrs["subject_id"]

    print(f"Subject: {subject_id}  ({h5_path.name})")
    print(f"top {top_percent}% edges, every {epoch_stride}th epoch, bands={bands}\n")

    print(f"{'Band':<8} {'Baseline':>10} {'REM obs':>10} {'REM enr':>9} "
          f"{'nonREM obs':>11} {'nonREM enr':>11}")
    print("-" * 64)

    for band in bands:
        rem_frac, non_rem_frac, baseline = check_one_subject_by_stage(
            h5_path, band, top_percent, epoch_stride
        )
        if baseline is None:
            print(f"{band:<8} (band not present in this file)")
            continue

        rem_enr = rem_frac / baseline if baseline > 0 and not np.isnan(rem_frac) else np.nan
        non_rem_enr = (non_rem_frac / baseline
                        if baseline > 0 and not np.isnan(non_rem_frac) else np.nan)

        rem_str = f"{rem_frac:>9.1%}" if not np.isnan(rem_frac) else "      n/a"
        non_rem_str = f"{non_rem_frac:>10.1%}" if not np.isnan(non_rem_frac) else "       n/a"
        rem_enr_str = f"{rem_enr:>8.2f}x" if not np.isnan(rem_enr) else "     n/a"
        non_rem_enr_str = f"{non_rem_enr:>10.2f}x" if not np.isnan(non_rem_enr) else "       n/a"

        flag = ""
        if not np.isnan(rem_enr) and rem_enr > 1.5:
            flag = "  <-- [FLAG: REM]"
        elif not np.isnan(non_rem_enr) and non_rem_enr > 1.5:
            flag = "  <-- [FLAG: non-REM]"

        print(f"{band:<8} {baseline:>9.1%} {rem_str} {rem_enr_str} "
              f"{non_rem_str} {non_rem_enr_str}{flag}")

    print("\nCompare the REM and non-REM columns, not just the overall level: real EOG "
          "contamination should show up mainly (or only) in REM, since that's where eye "
          "movement density is high. Enrichment that's equally high (or equally absent) "
          "in both REM and non-REM is less consistent with an eye-movement explanation, "
          "and worth a second look -- it may be a different artifact, or a real "
          "neurogenic pattern rather than an artifact at all.")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
        if path.is_dir():
            run_frontal_check(path)
        else:
            run_frontal_check_subject(path)
    else:
        run_frontal_check(Path("graphs_output"))