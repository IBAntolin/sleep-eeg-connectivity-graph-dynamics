"""
Phase 2 -- Graph Construction (wPLI, per band, per epoch)

Builds one 83x83 weighted phase-lag-index (wPLI) connectivity matrix per
frequency band from a single epoch's already-filtered EEG signal. Designed
to be called per epoch, inline, right after preprocessing.iter_clean_epochs
yields a filtered signal -- so the raw filtered signal never has to be
written to disk, only the small resulting graphs.

Band edges are a standard sleep-EEG convention; adjust BANDS if your
analysis calls for different cutoffs.
"""

import numpy as np
from scipy.signal import butter, filtfilt, hilbert

BANDS = {
    "delta": (0.5, 4.0),
    "theta": (4.0, 8.0),
    "alpha": (8.0, 12.0),
    "sigma": (12.0, 16.0),
    "beta": (16.0, 30.0),
}


def band_filter(signal: np.ndarray, sfreq: float, low: float, high: float,
                 order: int = 4) -> np.ndarray:
    """Zero-phase Butterworth band-pass for one frequency band, along the last axis."""
    nyq = sfreq / 2
    b, a = butter(order, [low / nyq, high / nyq], btype="band")
    return filtfilt(b, a, signal, axis=-1)


def wpli_matrix(analytic_signal: np.ndarray) -> np.ndarray:
    """
    Weighted phase-lag index across all channel pairs.

    analytic_signal: complex array, shape (n_channels, n_samples), e.g. the
    output of scipy.signal.hilbert applied along the time axis.

    Computed row-by-row rather than as one (n_channels, n_channels,
    n_samples) broadcast -- the full broadcast would need several GB for
    83 channels x 30,000 samples; row-by-row keeps peak memory to a single
    (n_channels, n_samples) array (a few MB) at a time.
    """
    n_channels = analytic_signal.shape[0]
    wpli = np.zeros((n_channels, n_channels), dtype=np.float32)

    for i in range(n_channels):
        csd = analytic_signal[i][None, :] * np.conj(analytic_signal)  # (n_channels, n_samples)
        imag_csd = csd.imag
        numerator = np.abs(imag_csd.mean(axis=1))
        denominator = np.abs(imag_csd).mean(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            row = np.where(denominator > 0, numerator / denominator, 0.0)
        wpli[i] = row

    np.fill_diagonal(wpli, 0.0)  # self-connections undefined (imag part is 0), set to 0
    return wpli


def build_epoch_graphs(signal: np.ndarray, sfreq: float,
                        bands: dict = BANDS) -> dict[str, np.ndarray]:
    """
    Build one wPLI graph per frequency band for a single epoch's signal.

    signal: real array, shape (n_channels, n_samples) -- already broadband
    filtered (e.g. by preprocessing.iter_clean_epochs).

    Returns: {band_name: (n_channels, n_channels) float32 wPLI matrix}
    """
    graphs = {}
    for band_name, (low, high) in bands.items():
        band_signal = band_filter(signal, sfreq, low, high)
        analytic = hilbert(band_signal, axis=-1).astype(np.complex64)
        graphs[band_name] = wpli_matrix(analytic)
    return graphs
