"""
Find which epoch and channel pair produced the global max wPLI value for a
given subject/band, and check whether either channel was flagged bad
(interpolated) in that specific epoch -- to tell apart a genuine high-
synchronization finding from a possible interpolation artifact (an
interpolated channel is, by construction, a blend of its neighbors, which
can inflate wPLI between it and those neighbors).

Also runs a broader, less anecdotal check: compares wPLI statistics for
edges that touch an interpolated channel against edges that don't, across
ALL epochs (not just the single global max) -- a systematic pattern here
matters more than one example.

CAVEAT: this relies on an unverified assumption -- that the artifact
matrix's channel columns are in the same order as channel_names (i.e. the
EEG channel order used everywhere else in the pipeline). This was never
independently confirmed against real column headers (the .mat file has no
channel-name metadata). A rough plausibility check is included: channels
known to be typically noisier in real recordings (temporal/mastoid sites
like T9/T10/SO1/SO2) should show a higher bad-fraction than central sites
if the alignment is correct -- if that pattern is reversed or absent,
the column order may not actually match, and both the interpolation
itself and this check's conclusions would need re-examining.

Usage:
    from check_max_wpli_source import find_max_wpli_source, compare_interpolated_vs_clean_edges
    find_max_wpli_source("EPCTL06.h5", "EPCTL06_artndxn.mat", band="delta")
    compare_interpolated_vs_clean_edges("EPCTL06.h5", "EPCTL06_artndxn.mat", band="delta")
"""

from pathlib import Path

import h5py
import numpy as np

from preprocessing import load_artifact_matrix, parse_pos_file


def _load(h5_path):
    with h5py.File(h5_path, "r") as f:
        channel_names = [c.decode() if isinstance(c, bytes) else c for c in f["channel_names"][:]]
        epoch_idx = f["epoch_idx"][:]
        bands = list(f.attrs["bands"])
        band_data = {band: f[band][:] for band in bands}
    return channel_names, epoch_idx, band_data


def find_max_wpli_source(h5_path: str | Path, artifact_path: str | Path, band: str = "delta"):
    channel_names, epoch_idx, band_data = _load(h5_path)
    graphs = band_data[band]  # (n_epochs, n_channels, n_channels)

    flat_argmax = np.argmax(graphs)
    e, i, j = np.unravel_index(flat_argmax, graphs.shape)
    max_val = graphs[e, i, j]
    orig_epoch_idx = int(epoch_idx[e])
    ch_i, ch_j = channel_names[i], channel_names[j]

    print(f"Max {band} wPLI = {max_val:.4f}")
    print(f"  stored epoch position: {e} (original scoring-file epoch index: {orig_epoch_idx})")
    print(f"  channel pair: {ch_i} <-> {ch_j}")

    artifact_matrix = load_artifact_matrix(artifact_path)
    if orig_epoch_idx >= len(artifact_matrix):
        print("  could not cross-reference -- epoch index out of artifact matrix range")
        return e, i, j, max_val

    bad_row = ~artifact_matrix[orig_epoch_idx]  # True = bad; SAME column-order assumption as the pipeline
    i_bad = bool(bad_row[i]) if i < len(bad_row) else None
    j_bad = bool(bad_row[j]) if j < len(bad_row) else None
    print(f"  {ch_i} flagged bad in this epoch? {i_bad}")
    print(f"  {ch_j} flagged bad in this epoch? {j_bad}")

    if i_bad or j_bad:
        print("  [FLAG] one of the channels in the max-value pair was interpolated in this "
              "epoch -- this specific value may be an interpolation artifact rather than a "
              "genuine synchronization finding.")
    else:
        print("  [OK] neither channel was flagged bad in this epoch -- this max value comes "
              "from an unmodified signal pair, consistent with a genuine finding (e.g. delta "
              "synchronization is well documented as high during deep sleep).")

    return e, i, j, max_val


def compare_interpolated_vs_clean_edges(h5_path: str | Path, artifact_path: str | Path,
                                         band: str = "delta") -> dict:
    """
    Systematic version of the check above: across every epoch, split edges
    into "touches an interpolated channel" vs "both channels clean", and
    compare their wPLI value distributions. A large, consistent gap here
    (interpolated edges running noticeably higher) would indicate a
    systematic inflation pattern worth correcting for in later analysis,
    not just a one-off in the single max value.
    """
    channel_names, epoch_idx, band_data = _load(h5_path)
    graphs = band_data[band]  # (n_epochs, n_channels, n_channels)
    n_epochs, n_channels, _ = graphs.shape

    artifact_matrix = load_artifact_matrix(artifact_path)

    touches_interp_vals = []
    clean_vals = []

    for e in range(n_epochs):
        orig_idx = int(epoch_idx[e])
        if orig_idx >= len(artifact_matrix):
            continue
        bad = ~artifact_matrix[orig_idx]
        if not bad.any():
            # no interpolation this epoch -- every edge here is "clean", sample sparsely
            # to keep this from dwarfing the comparison
            iu = np.triu_indices(n_channels, k=1)
            vals = graphs[e][iu]
            clean_vals.append(vals[::37])  # sparse sample, arbitrary stride
            continue

        iu_i, iu_j = np.triu_indices(n_channels, k=1)
        edge_touches_bad = bad[iu_i] | bad[iu_j]
        vals = graphs[e][iu_i, iu_j]
        touches_interp_vals.append(vals[edge_touches_bad])
        clean_vals.append(vals[~edge_touches_bad][::37])

    touches_interp_vals = np.concatenate(touches_interp_vals) if touches_interp_vals else np.array([])
    clean_vals = np.concatenate(clean_vals) if clean_vals else np.array([])

    print(f"--- {band} band: interpolated-touching edges vs. clean edges ---")
    print(f"n edges touching an interpolated channel: {len(touches_interp_vals)}")
    print(f"  mean={touches_interp_vals.mean():.4f}  median={np.median(touches_interp_vals):.4f}  "
          f"max={touches_interp_vals.max():.4f}" if len(touches_interp_vals) else "  (none)")
    print(f"n clean edges (sampled): {len(clean_vals)}")
    print(f"  mean={clean_vals.mean():.4f}  median={np.median(clean_vals):.4f}  "
          f"max={clean_vals.max():.4f}" if len(clean_vals) else "  (none)")

    if len(touches_interp_vals) and len(clean_vals):
        diff = touches_interp_vals.mean() - clean_vals.mean()
        print(f"mean difference (interpolated-touching minus clean): {diff:+.4f}")
        if abs(diff) > 0.05:
            print("[FLAG] noticeable systematic difference -- worth investigating whether "
                  "interpolation is inflating (or deflating) wPLI at those edges.")
        else:
            print("[OK] no large systematic difference between the two groups.")

    return {
        "touches_interpolated": touches_interp_vals,
        "clean": clean_vals,
    }


def compare_interpolated_edges_three_way(h5_path: str | Path, artifact_path: str | Path,
                                          band: str = "delta") -> dict:
    """
    Refined version of compare_interpolated_vs_clean_edges() that splits
    edges into THREE groups instead of two:
      - both_clean:        neither channel interpolated in this epoch
      - one_interpolated:  exactly one of the two channels interpolated
      - both_interpolated: both channels interpolated in this epoch

    Rationale: the original two-way check (touches-interpolated vs. clean)
    can't distinguish "interpolated channel talking to a genuine neighbor"
    from "two interpolated channels -- both reconstructed from overlapping
    sets of nearby good channels -- talking to each other". The latter is
    a distinct, more concerning failure mode: two channels built from
    overlapping source signals would be expected to look artificially
    synchronized with each other specifically, more so than either looks
    with an untouched channel. If inflation is concentrated in
    both_interpolated, that supports this mechanism; if it's spread evenly
    across one_interpolated and both_interpolated, that points elsewhere.
    """
    channel_names, epoch_idx, band_data = _load(h5_path)
    graphs = band_data[band]  # (n_epochs, n_channels, n_channels)
    n_epochs, n_channels, _ = graphs.shape

    artifact_matrix = load_artifact_matrix(artifact_path)

    both_clean_vals = []
    one_interp_vals = []
    both_interp_vals = []

    iu_i, iu_j = np.triu_indices(n_channels, k=1)

    for e in range(n_epochs):
        orig_idx = int(epoch_idx[e])
        if orig_idx >= len(artifact_matrix):
            continue
        bad = ~artifact_matrix[orig_idx]
        vals = graphs[e][iu_i, iu_j]

        if not bad.any():
            # no interpolation this epoch -- every edge is both_clean, sample
            # sparsely to keep this from dwarfing the comparison
            both_clean_vals.append(vals[::37])
            continue

        i_bad = bad[iu_i]
        j_bad = bad[iu_j]
        n_bad_endpoints = i_bad.astype(int) + j_bad.astype(int)  # 0, 1, or 2

        both_clean_vals.append(vals[n_bad_endpoints == 0][::37])
        one_interp_vals.append(vals[n_bad_endpoints == 1])
        both_interp_vals.append(vals[n_bad_endpoints == 2])

    groups = {
        "both_clean": np.concatenate(both_clean_vals) if both_clean_vals else np.array([]),
        "one_interpolated": np.concatenate(one_interp_vals) if one_interp_vals else np.array([]),
        "both_interpolated": np.concatenate(both_interp_vals) if both_interp_vals else np.array([]),
    }

    print(f"--- {band} band: three-way edge comparison ---")
    for name, vals in groups.items():
        if len(vals):
            print(f"{name:18s} n={len(vals):6d}  mean={vals.mean():.4f}  "
                  f"median={np.median(vals):.4f}  max={vals.max():.4f}")
        else:
            print(f"{name:18s} (none)")

    baseline = groups["both_clean"]
    if len(baseline):
        for name in ("one_interpolated", "both_interpolated"):
            vals = groups[name]
            if len(vals):
                diff = vals.mean() - baseline.mean()
                flag = " [FLAG]" if abs(diff) > 0.05 else ""
                print(f"mean difference ({name} minus both_clean): {diff:+.4f}{flag}")

    if len(groups["one_interpolated"]) and len(groups["both_interpolated"]):
        gap = groups["both_interpolated"].mean() - groups["one_interpolated"].mean()
        print(f"mean difference (both_interpolated minus one_interpolated): {gap:+.4f}")
        if gap > 0.03:
            print("[FLAG] inflation concentrates specifically in both-interpolated edges -- "
                  "consistent with the overlapping-neighbor-reconstruction hypothesis "
                  "(two channels rebuilt from shared nearby good channels looking "
                  "artificially synchronized with each other).")
        else:
            print("[INFO] inflation is not clearly concentrated in both-interpolated edges "
                  "specifically -- one_interpolated and both_interpolated look similar, so "
                  "the overlapping-neighbor hypothesis is not strongly supported by this alone.")

    return groups


def per_epoch_bad_channel_summary(artifact_path: str | Path, subject_label: str = "") -> dict:
    """
    Summarizes, across all epochs in a subject's artifact matrix, how many
    channels were flagged bad per epoch. Used to check whether subjects
    showing the interpolation-inflation flag (e.g. in
    compare_interpolated_edges_three_way) tend to have higher bad-channel
    counts per epoch than subjects that don't -- which would support the
    hypothesis that inflation stems from overlapping-neighbor
    reconstruction when many channels are bad at once, not just 1-3.
    """
    artifact_matrix = load_artifact_matrix(artifact_path)
    n_bad_per_epoch = (~artifact_matrix).sum(axis=1)

    label = f" ({subject_label})" if subject_label else ""
    print(f"--- per-epoch bad-channel counts{label} ---")
    print(f"n_epochs: {len(n_bad_per_epoch)}")
    print(f"mean bad channels/epoch: {n_bad_per_epoch.mean():.2f}")
    print(f"median: {np.median(n_bad_per_epoch):.1f}  max: {n_bad_per_epoch.max()}")

    n_channels = artifact_matrix.shape[1]
    thresholds = [1, 3, 5, 10, int(round(0.15 * n_channels)), int(round(0.20 * n_channels))]
    thresholds = sorted(set(t for t in thresholds if t > 0))
    for t in thresholds:
        frac = (n_bad_per_epoch >= t).mean()
        print(f"epochs with >= {t} bad channels: {frac:.1%}")

    return {
        "n_bad_per_epoch": n_bad_per_epoch,
        "mean": float(n_bad_per_epoch.mean()),
        "median": float(np.median(n_bad_per_epoch)),
        "max": int(n_bad_per_epoch.max()),
    }


def check_interpolation_distance_effect(h5_path: str | Path, artifact_path: str | Path,
                                         pos_path: str | Path, band: str = "delta",
                                         n_bins: int = 4) -> dict:
    """
    Tests the mechanism directly implicated by interpolate_bads(): spherical-
    spline interpolation reconstructs a bad channel as a WEIGHTED
    combination of good channels, with weights that decay with physical
    distance (channels close to the interpolated site contribute more).
    If that's what's driving the inflation seen in one_interpolated edges,
    the inflation should be concentrated on edges to NEARBY good channels
    and fade out for good channels far from the interpolated site.

    Restricted to edges with exactly one interpolated endpoint (the group
    that showed the consistent, robust inflation across subjects), since
    "distance to the interpolated channel" is only unambiguous when there's
    a single interpolated channel in the pair.

    NOTE: this uses physical distance as a proxy for actual spline weight
    (the true per-channel weights MNE computes internally aren't exposed
    through the public interpolate_bads() API). Distance is a reasonable
    proxy since spherical-spline weights are known to decay smoothly with
    it, but it isn't the literal weight MNE used -- treat a positive
    finding here as supportive, not as exact confirmation of MNE's
    internal weighting.
    """
    channel_names, epoch_idx, band_data = _load(h5_path)
    graphs = band_data[band]  # (n_epochs, n_channels, n_channels)
    n_epochs, n_channels, _ = graphs.shape

    # match channel_names against the .pos file the same way build_montage_info
    # does (case-insensitive, strips common ref/avg/car suffixes) -- done
    # locally here to avoid pulling in mne just to read positions
    ch_pos = parse_pos_file(pos_path)

    def normalize(name: str) -> str:
        name = name.lower().strip()
        for suffix in ["-ref", "_ref", "-avg", "_avg", "-car", "_car"]:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
        return name.strip("-+_")

    ch_pos_norm = {normalize(k): v for k, v in ch_pos.items()}
    missing = [name for name in channel_names if normalize(name) not in ch_pos_norm]
    if missing:
        raise ValueError(
            f"{len(missing)} channel(s) from the h5 have no match in the .pos "
            f"file (even after normalization): {missing}"
        )
    positions = np.array([ch_pos_norm[normalize(name)] for name in channel_names])  # (n_channels, 3)

    # pairwise Euclidean distance between all channels (proxy for spline weight)
    diffs = positions[:, None, :] - positions[None, :, :]
    dist_matrix = np.linalg.norm(diffs, axis=-1)  # (n_channels, n_channels)

    artifact_matrix = load_artifact_matrix(artifact_path)
    iu_i, iu_j = np.triu_indices(n_channels, k=1)

    edge_vals = []
    edge_dists = []

    for e in range(n_epochs):
        orig_idx = int(epoch_idx[e])
        if orig_idx >= len(artifact_matrix):
            continue
        bad = ~artifact_matrix[orig_idx]
        if not bad.any():
            continue

        i_bad = bad[iu_i]
        j_bad = bad[iu_j]
        n_bad_endpoints = i_bad.astype(int) + j_bad.astype(int)
        one_interp = n_bad_endpoints == 1
        if not one_interp.any():
            continue

        vals = graphs[e][iu_i, iu_j][one_interp]
        dists = dist_matrix[iu_i, iu_j][one_interp]
        edge_vals.append(vals)
        edge_dists.append(dists)

    edge_vals = np.concatenate(edge_vals) if edge_vals else np.array([])
    edge_dists = np.concatenate(edge_dists) if edge_dists else np.array([])

    print(f"--- {band} band: interpolation-distance effect (one_interpolated edges only) ---")
    print(f"n edges: {len(edge_vals)}")
    if len(edge_vals) < 2:
        print("  not enough edges to analyze")
        return {"values": edge_vals, "distances": edge_dists}

    corr = np.corrcoef(edge_dists, edge_vals)[0, 1]
    print(f"correlation (distance to interpolated channel vs. wPLI): r={corr:.3f}")
    if corr < -0.05:
        print("[SUPPORTIVE] negative correlation -- wPLI tends to be higher for edges "
              "closer to the interpolated channel, consistent with spline-weight-driven "
              "inflation (nearby channels contribute more to the reconstruction, so they "
              "end up looking more synchronized with it).")
    else:
        print("[NOT SUPPORTIVE] no clear negative correlation -- distance to the "
              "interpolated channel doesn't explain the inflation by itself.")

    # bin by distance quantile and compare means directly (easier to read than
    # a bare correlation coefficient, and robust to a non-linear relationship)
    quantile_edges = np.quantile(edge_dists, np.linspace(0, 1, n_bins + 1))
    print(f"\nBy distance bin (nearest to farthest):")
    bin_means = []
    for b in range(n_bins):
        lo, hi = quantile_edges[b], quantile_edges[b + 1]
        mask = (edge_dists >= lo) & (edge_dists <= hi if b == n_bins - 1 else edge_dists < hi)
        vals_in_bin = edge_vals[mask]
        if len(vals_in_bin):
            bin_means.append(vals_in_bin.mean())
            print(f"  bin {b+1}/{n_bins}  dist [{lo:.3f}, {hi:.3f})m  n={len(vals_in_bin):6d}  "
                  f"mean wPLI={vals_in_bin.mean():.4f}")
        else:
            bin_means.append(np.nan)
            print(f"  bin {b+1}/{n_bins}  dist [{lo:.3f}, {hi:.3f})m  (none)")

    if len(bin_means) >= 2 and not np.isnan(bin_means[0]) and not np.isnan(bin_means[-1]):
        drop = bin_means[0] - bin_means[-1]
        print(f"\nnearest-bin minus farthest-bin mean wPLI: {drop:+.4f}")

    return {"values": edge_vals, "distances": edge_dists, "correlation": corr, "bin_means": bin_means}


def channel_bad_fraction_plausibility_check(artifact_path: str | Path, channel_names: list[str]) -> None:
    """
    Rough plausibility check for whether the artifact matrix's column
    order actually matches channel_names: prints each channel's overall
    bad-fraction across the night, sorted worst to best. Channels at
    typically noisier real-world sites (temporal, mastoid -- e.g. T9, T10,
    SO1, SO2 in this dataset's naming) showing up disproportionately near
    the top would be consistent with correct alignment; if central/frontal
    channels dominate the "worst" list instead, or the ranking looks
    physiologically implausible, the column order may not match.
    """
    artifact_matrix = load_artifact_matrix(artifact_path)
    bad_fraction = (~artifact_matrix).mean(axis=0)

    n = min(len(bad_fraction), len(channel_names))
    pairs = sorted(zip(channel_names[:n], bad_fraction[:n]), key=lambda x: -x[1])

    print("Channels ranked by overall bad-fraction (worst first):")
    for name, frac in pairs[:15]:
        print(f"  {name:8s} {frac:.1%}")