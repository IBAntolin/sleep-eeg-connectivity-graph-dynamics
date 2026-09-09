"""
Quick sanity check for one subject's processed .h5 file -- confirms the
structure matches what process_subject() is supposed to write: consistent
shapes across bands, matching epoch counts, sane wPLI value ranges, and a
readable stage-label breakdown.

Usage:
    python sanity_check.py path/to/EPCTL01.h5
or from a notebook:
    from sanity_check import sanity_check
    sanity_check("anphy_sleep_data/processed/EPCTL01.h5")
"""

from pathlib import Path

import h5py
import numpy as np


def sanity_check(h5_path: str | Path, expected_n_channels: int = 83) -> bool:
    h5_path = Path(h5_path)
    ok = True

    with h5py.File(h5_path, "r") as f:
        print(f"--- {h5_path.name} ---")
        subject_id = f.attrs.get("subject_id", "?")
        bands = list(f.attrs.get("bands", []))
        n_epochs_attr = f.attrs.get("n_epochs", None)
        print(f"subject_id: {subject_id}")
        print(f"bands recorded in attrs: {bands}")

        # 1. required datasets present
        required = {"epoch_idx", "stage", "channel_names", *bands}
        missing = required - set(f.keys())
        if missing:
            print(f"[FAIL] missing datasets: {missing}")
            ok = False
        else:
            print("[OK] all expected datasets present")

        if not bands:
            print("[FAIL] no bands recorded -- nothing further to check")
            return False

        # 2. channel_names length
        n_channels = f["channel_names"].shape[0]
        print(f"channel_names: {n_channels} channels")
        if n_channels != expected_n_channels:
            print(f"[FAIL] expected {expected_n_channels} channels, got {n_channels}")
            ok = False
        else:
            print(f"[OK] channel count matches expected ({expected_n_channels})")

        # 3. epoch_idx / stage lengths match each other
        n_epoch_idx = f["epoch_idx"].shape[0]
        n_stage = f["stage"].shape[0]
        print(f"epoch_idx: {n_epoch_idx} entries, stage: {n_stage} entries")
        if n_epoch_idx != n_stage:
            print(f"[FAIL] epoch_idx ({n_epoch_idx}) and stage ({n_stage}) length mismatch")
            ok = False
        if n_epochs_attr is not None and n_epochs_attr != n_epoch_idx:
            print(f"[FAIL] attrs['n_epochs'] ({n_epochs_attr}) != len(epoch_idx) ({n_epoch_idx})")
            ok = False

        # 4. every band's shape matches (n_epochs, n_channels, n_channels)
        # and matches epoch_idx's length
        shapes = {}
        for band in bands:
            shape = f[band].shape
            shapes[band] = shape
            print(f"{band}: shape {shape}")

        expected_shape = (n_epoch_idx, n_channels, n_channels)
        for band, shape in shapes.items():
            if shape != expected_shape:
                print(f"[FAIL] {band} shape {shape} != expected {expected_shape}")
                ok = False

        if len(set(shapes.values())) == 1:
            print("[OK] all bands have identical, consistent shapes")
        else:
            print("[FAIL] bands have inconsistent shapes with each other")
            ok = False

        # 5. wPLI value sanity: should be in [0, 1], not NaN, not all-zero
        for band in bands:
            arr = f[band][:5] if f[band].shape[0] > 5 else f[band][:]  # sample a few epochs
            n_nan = np.isnan(arr).sum()
            below0 = (arr < 0).sum()
            above1 = (arr > 1).sum()
            nonzero_frac = np.count_nonzero(arr) / arr.size
            print(f"{band}: sampled range [{arr.min():.4f}, {arr.max():.4f}], "
                  f"NaNs={n_nan}, <0={below0}, >1={above1}, nonzero={nonzero_frac:.1%}")
            if n_nan > 0 or below0 > 0 or above1 > 0:
                print(f"[FAIL] {band} has out-of-range or NaN values")
                ok = False
            if nonzero_frac < 0.01:
                print(f"[WARN] {band} is almost entirely zero -- check wPLI computation")

        # 6. stage label breakdown
        stages = np.array([s.decode() if isinstance(s, bytes) else s for s in f["stage"][:]])
        unique, counts = np.unique(stages, return_counts=True)
        print("Stage breakdown:")
        for stage, count in zip(unique, counts):
            print(f"  {stage:10s} {count:5d} epochs")

        # 6b. epoch_idx gap check + interpolation activity check.
        #
        # Under the current pipeline (interpolate bad channels within
        # MAX_BAD_FRACTION, only drop epochs beyond that), it's entirely
        # normal for a subject to have ZERO dropped epochs if none of
        # their artifact-flagged epochs exceeded the threshold -- so gaps
        # in epoch_idx are no longer a reliable signal that artifact
        # handling is working. n_bad_channels is: if it's present and has
        # some nonzero entries, that confirms the artifact matrix was
        # read and interpolation actually happened, independent of
        # whether any epoch was dropped outright.
        epoch_idx = f["epoch_idx"][:]
        if len(epoch_idx) > 1:
            span = int(epoch_idx.max() - epoch_idx.min() + 1)
            n_stored = len(epoch_idx)
            n_dropped = span - n_stored
            print(f"epoch_idx range: {epoch_idx.min()}-{epoch_idx.max()} "
                  f"(span {span}, {n_stored} stored, {n_dropped} dropped)")

        if "n_bad_channels" in f:
            n_bad = f["n_bad_channels"][:]
            n_epochs_with_interp = int(np.count_nonzero(n_bad))
            total_channels_interp = int(n_bad.sum())
            print(f"n_bad_channels: {n_epochs_with_interp}/{len(n_bad)} epochs had "
                  f"1+ interpolated channel(s), {total_channels_interp} channel-epochs "
                  f"interpolated in total (max in one epoch: {n_bad.max() if len(n_bad) else 0})")
            if n_epochs_with_interp == 0 and n_dropped == 0:
                print("[INFO] no interpolation activity and no drops -- either this "
                      "subject's recording was fully clean, or artifact_path wasn't "
                      "actually passed to process_subject(). Worth a quick check if "
                      "unexpected.")
            else:
                print("[OK] artifact matrix is being used (interpolation and/or drops occurred)")
        else:
            print("[INFO] no n_bad_channels dataset -- file predates the interpolation "
                  "pipeline, or artifact_path wasn't passed. epoch_idx gaps are the only "
                  "available signal for this file.")
            if len(epoch_idx) > 1 and n_dropped == 0:
                print("[WARN] epoch_idx is perfectly consecutive and there's no "
                      "n_bad_channels dataset to cross-check -- can't confirm whether "
                      "artifact filtering ran for this (older-format) file.")

        # 7. diagonal should be ~0 (self-connections set to 0 in wpli_matrix)
        for band in bands:
            sample = f[band][0]
            diag_max = np.abs(np.diag(sample)).max()
            if diag_max > 1e-6:
                print(f"[WARN] {band} epoch 0 diagonal not ~0 (max abs {diag_max:.2e})")

    print(f"\n{'PASSED' if ok else 'FAILED'} sanity check for {h5_path.name}")
    return ok


if __name__ == "__main__":
    import sys
    sanity_check(sys.argv[1])