"""
Checks whether peripheral/outer-ring channels (temporal, jaw-adjacent,
zygomatic, supraorbital sites -- the ones most prone to picking up muscle
(EMG) or eye (EOG) artifact rather than clean cortical signal) are
over-represented among the STRONGEST wPLI edges, and whether that
over-representation is worse in higher frequency bands (alpha/sigma/beta,
where EMG contamination classically shows up) than in lower bands
(delta/theta) -- across many epochs and subjects, not just one example
epoch.

Method: for each band, for each sampled epoch, take the top X% of edges
(same selection used for plotting) and compute what fraction of them
touch at least one peripheral channel. Compare that against the BASELINE
fraction you'd expect just from how many peripheral channels there are,
if peripheral-ness had no relationship to edge strength at all. A gap
between "observed in top edges" and "expected by chance" is the signal to
look for -- and if that gap is much bigger in alpha/sigma/beta than in
delta/theta, that's evidence of a real, systematic higher-band artifact
pattern rather than a one-epoch fluke.

Usage:
    from peripheral_channel_check import run_check, compare_peripheral_pair_types

    # Chance-enrichment check (top-edge fraction touching a peripheral channel)
    run_check(Path("graphs_output"), top_percent=4, epoch_stride=10)

    # Three-way both_central / one_peripheral / both_peripheral comparison,
    # across all five bands at once -- this is the one to rerun after any
    # EMG correction change, to confirm it actually fixed the affected
    # band(s) without disturbing the others
    compare_peripheral_pair_types(Path("graphs_output"), epoch_stride=10)
"""

from pathlib import Path

import h5py
import numpy as np

# Outer-ring / peripheral channels in this montage: extended 10-10/10-5
# sites beyond the standard 10-20 ring (F9-F12, FT9-FT12, T9/T10,
# TP9-TP12, P9-P12), zygomatic (ZY1/ZY2, near jaw/ear), and supraorbital
# (SO1/SO2, near the eyes) -- all sites where EMG (muscle) or EOG (eye)
# contamination is classically more likely than at central/standard sites.
PERIPHERAL_CHANNELS = {
    "F9", "F10", "F11", "F12", "FT9", "FT10", "FT11", "FT12", "T9", "T10",
    "TP9", "TP10", "TP11", "TP12", "P9", "P10", "P11", "P12", "ZY1", "ZY2", "SO1", "SO2",
}


def exclude_peripheral_channels(
    con: np.ndarray, ch_names: list[str], peripheral: set[str] = PERIPHERAL_CHANNELS,
) -> tuple[np.ndarray, list[str]]:
    """
    Return a smaller version of a wPLI adjacency matrix with all
    peripheral/outer-ring channels removed -- both their rows and columns
    dropped, leaving only central-channel-to-central-channel edges.

    This is the lighter alternative chosen over full ICA-based muscle-
    component removal: rather than reconstructing the signal, re-run the
    same graph metric or stage comparison twice -- once on the full graph,
    once on this peripheral-excluded version -- specifically for BETA band
    (where run_check() found a real, graded 1.36x enrichment of peripheral
    channels among the strongest edges, consistent with EMG/muscle
    contamination, which is itself a concern because muscle tone changes
    systematically across sleep stages -- near-normal in Wake, reduced in
    NREM, atonic in REM -- so it could masquerade as a stage-related
    network finding that isn't really one). If a beta-band stage
    difference holds up in both versions, that's evidence it's a real
    network effect and not an artifact of muscle tone dropping out toward
    REM. If it disappears once peripheral channels are excluded, that
    tells you the original result shouldn't be trusted. Delta/theta/sigma
    didn't show this pattern (0.94x, 0.94x, 1.10x respectively, from a
    28-subject check) and don't need this treatment.

    Example:
        con_full = f["beta"][epoch]
        con_central, ch_names_central = exclude_peripheral_channels(con_full, channel_names)
        # compute your graph metric on both con_full and con_central,
        # compare whether the stage-comparison conclusion changes
    """
    keep_mask = np.array([name not in peripheral for name in ch_names])
    con_filtered = con[np.ix_(keep_mask, keep_mask)]
    ch_names_filtered = [name for name, keep in zip(ch_names, keep_mask) if keep]
    return con_filtered, ch_names_filtered


def _split_edges_by_peripheral(con: np.ndarray, ch_names: list[str],
                                peripheral: set[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Split one epoch's edges into three groups by how many of the two
    channels are peripheral: 0 (both central), 1 (one peripheral, one
    central), or 2 (both peripheral). Returns the wPLI values for each
    group as three separate arrays.
    """
    iu_i, iu_j = np.triu_indices_from(con, k=1)
    is_periph = np.array([name in peripheral for name in ch_names])
    n_periph_in_pair = is_periph[iu_i].astype(int) + is_periph[iu_j].astype(int)
    vals = con[iu_i, iu_j]

    return (
        vals[n_periph_in_pair == 0],  # both central
        vals[n_periph_in_pair == 1],  # one peripheral
        vals[n_periph_in_pair == 2],  # both peripheral
    )


ALL_BANDS = ("delta", "theta", "alpha", "sigma", "beta")


def compare_peripheral_pair_types(
    output_dir: str | Path, bands: tuple[str, ...] = ALL_BANDS, epoch_stride: int = 10,
) -> dict:
    """
    Splits every edge, across every subject's .h5 file, into three groups
    -- both channels central, exactly one peripheral, both peripheral --
    and compares their wPLI distributions for each band in `bands`
    (defaults to all five: delta/theta/alpha/sigma/beta). This is more
    specific than the earlier chance-enrichment check: it tells us WHERE
    within the "touches a peripheral channel" category the effect
    concentrates, and -- checked across the whole spectrum rather than
    just one band at a time -- whether a correction (e.g. the per-band
    EMG regression fix) actually resolved the problem everywhere it
    mattered, not just in the one band originally flagged.

    Why the both-peripheral distinction matters for choosing a
    correction: if the inflation is concentrated specifically in
    both-peripheral edges (two artifact-prone channels picking up the
    same shared muscle activity), a correction only needs to target that
    small group -- 231 of 3403 possible edges (~6.8%) for this
    83-channel montage (22 peripheral channels) -- rather than every edge
    touching a peripheral channel (1573 of 3403, ~46%, the "one or more
    peripheral" group used in the earlier chance check). That's a much
    smaller, more defensible footprint: most peripheral-channel
    information (their connections to central/brain channels) would be
    left completely untouched.

    Running all five bands together (rather than one band vs. a single
    reference band) is what lets you see the full picture in one pass:
    e.g. confirming the effect really was concentrated in beta before a
    fix, and checking after a fix that beta improved without silently
    introducing a new problem in a band that was previously clean
    (sigma, alpha) -- something a two-band comparison can't reveal.
    """
    output_dir = Path(output_dir)
    h5_files = sorted(output_dir.glob("*.h5"))
    if not h5_files:
        print(f"No .h5 files found in {output_dir}")
        return {}

    bands_to_check = list(bands)
    results = {b: {"both_central": [], "one_peripheral": [], "both_peripheral": []} for b in bands_to_check}

    for h5_path in h5_files:
        with h5py.File(h5_path, "r") as f:
            ch_names = [c.decode() if isinstance(c, bytes) else c for c in f["channel_names"][:]]
            n_epochs = f.attrs["n_epochs"]
            for b in bands_to_check:
                if b not in f:
                    continue
                for e in range(0, n_epochs, epoch_stride):
                    both_c, one_p, both_p = _split_edges_by_peripheral(f[b][e], ch_names, PERIPHERAL_CHANNELS)
                    results[b]["both_central"].append(both_c)
                    results[b]["one_peripheral"].append(one_p)
                    results[b]["both_peripheral"].append(both_p)

    print(f"Checked {len(h5_files)} subject(s), every {epoch_stride}th epoch\n")

    summary = {}
    for b in bands_to_check:
        print(f"--- {b} band: three-way edge comparison ---")
        group_summary = {}
        for group_name in ("both_central", "one_peripheral", "both_peripheral"):
            vals = np.concatenate(results[b][group_name]) if results[b][group_name] else np.array([])
            group_summary[group_name] = vals
            if len(vals):
                print(f"  {group_name:16s} n={len(vals):9d}  mean={vals.mean():.4f}  "
                      f"median={np.median(vals):.4f}")
            else:
                print(f"  {group_name:16s} (none)")

        if len(group_summary["both_central"]) and len(group_summary["both_peripheral"]):
            diff_both_periph = group_summary["both_peripheral"].mean() - group_summary["both_central"].mean()
            diff_one_periph = (group_summary["one_peripheral"].mean() - group_summary["both_central"].mean()
                                if len(group_summary["one_peripheral"]) else float("nan"))
            print(f"  both_peripheral vs both_central: {diff_both_periph:+.4f}")
            print(f"  one_peripheral  vs both_central: {diff_one_periph:+.4f}")
            if abs(diff_both_periph) > abs(diff_one_periph) * 1.5:
                print("  [PATTERN] effect concentrates in both-peripheral edges -- a correction "
                      "targeting only that small group (~6.8% of edges) would likely be sufficient.")
            elif abs(diff_one_periph) > 0.03:
                print("  [PATTERN] effect is spread across one-peripheral edges too, not just "
                      "both-peripheral -- a narrowly-targeted correction may not be enough.")
            else:
                print("  [OK] no meaningful effect in either peripheral group for this band.")
        print()

        summary[b] = group_summary

    return summary


def _expected_baseline_fraction(ch_names: list[str], peripheral: set[str]) -> float:
    """
    Fraction of ALL possible edges that touch at least one peripheral
    channel, if peripheral-ness had no relationship to which edges are
    strong -- i.e. the null/chance expectation to compare "top edges"
    against.
    """
    n = len(ch_names)
    is_peripheral = np.array([name in peripheral for name in ch_names])
    n_periph = is_peripheral.sum()
    n_central = n - n_periph
    total_pairs = n * (n - 1) / 2
    central_only_pairs = n_central * (n_central - 1) / 2
    return 1 - (central_only_pairs / total_pairs)


def _top_edge_peripheral_fraction(con: np.ndarray, ch_names: list[str],
                                   peripheral: set[str], top_percent: float) -> float:
    """Fraction of the top X% strongest edges in this epoch/band that
    touch at least one peripheral channel."""
    iu_i, iu_j = np.triu_indices_from(con, k=1)
    vals = con[iu_i, iu_j]
    thresh = np.percentile(vals, 100 - top_percent)
    keep = vals >= thresh

    is_peripheral = np.array([name in peripheral for name in ch_names])
    touches_peripheral = is_peripheral[iu_i] | is_peripheral[iu_j]

    kept_touches_peripheral = touches_peripheral[keep]
    return kept_touches_peripheral.mean() if len(kept_touches_peripheral) else np.nan


def check_one_subject(h5_path: str | Path, top_percent: float = 4,
                       epoch_stride: int = 10) -> dict:
    """
    For one subject's .h5 file: for each band, sample every `epoch_stride`
    epochs (not all, for speed) and compute the mean top-edge peripheral
    fraction across those sampled epochs. Returns {band: mean_fraction}.
    """
    with h5py.File(h5_path, "r") as f:
        bands = list(f.attrs["bands"])
        ch_names = [c.decode() if isinstance(c, bytes) else c for c in f["channel_names"][:]]
        n_epochs = f.attrs["n_epochs"]

        results = {}
        for band in bands:
            fractions = []
            for e in range(0, n_epochs, epoch_stride):
                con = f[band][e]
                frac = _top_edge_peripheral_fraction(con, ch_names, PERIPHERAL_CHANNELS, top_percent)
                if not np.isnan(frac):
                    fractions.append(frac)
            results[band] = np.mean(fractions) if fractions else np.nan

    baseline = _expected_baseline_fraction(ch_names, PERIPHERAL_CHANNELS)
    return results, baseline


def run_check(output_dir: str | Path, top_percent: float = 4, epoch_stride: int = 10) -> None:
    """
    Run check_one_subject() across every *.h5 file in output_dir and print
    a per-band summary: mean observed top-edge peripheral fraction
    (averaged across subjects and sampled epochs) vs. the chance baseline,
    and the enrichment ratio between them.
    """
    output_dir = Path(output_dir)
    h5_files = sorted(output_dir.glob("*.h5"))
    if not h5_files:
        print(f"No .h5 files found in {output_dir}")
        return

    print(f"Checking {len(h5_files)} subject(s), top {top_percent}% edges, "
          f"every {epoch_stride}th epoch...\n")

    all_band_results: dict[str, list[float]] = {}
    baseline = None

    for h5_path in h5_files:
        results, subj_baseline = check_one_subject(h5_path, top_percent, epoch_stride)
        baseline = subj_baseline  # same for every subject (same 83 channels)
        for band, frac in results.items():
            all_band_results.setdefault(band, []).append(frac)

    print(f"Baseline (chance) fraction of edges touching a peripheral channel: {baseline:.1%}\n")
    print(f"{'Band':<8} {'Observed (top edges)':>22} {'Enrichment vs. chance':>24}")
    print("-" * 56)
    for band, fracs in all_band_results.items():
        mean_frac = np.nanmean(fracs)
        enrichment = mean_frac / baseline if baseline > 0 else np.nan
        flag = "  <-- [FLAG] much higher than chance" if enrichment > 1.5 else ""
        print(f"{band:<8} {mean_frac:>21.1%} {enrichment:>23.2f}x{flag}")

    print("\nA ratio near 1.0x means peripheral channels show up in the strongest edges "
          "about as often as pure chance would predict -- no artifact concentration. "
          "A ratio well above 1.0x, especially concentrated in alpha/sigma/beta rather "
          "than delta/theta, supports the muscle/eye-artifact concentration hypothesis "
          "from the single-epoch plot.")


if __name__ == "__main__":
    import sys
    output_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("graphs_output")
    run_check(output_dir)
    print()
    compare_peripheral_pair_types(output_dir)