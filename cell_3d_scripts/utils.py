import logging
import operator
import re
from collections import defaultdict
from collections.abc import Callable

import numpy as np
from brainglobe_utils.cells.cells import Cell
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

_cell_filter_pat = re.compile(
    r"^(log10:)?"
    r"([\w.]+)"
    r"([<>=!]{1,2})"
    r"(p|peak|mean|gaussianS|gaussian_p|gaussian_mean|biGaussianS|biGaussianTrough)?"
    r"(-?[0-9]*\.?[0-9]*)"
    r"(,-?[0-9]*\.?[0-9]*)?$"
)

# region, outside of which, it's definitely a bad value
MEASURES = {
    "intensity": {"log_x": True, "fit_gaussian": True},
    "paor_intensity_total": {"log_x": True, "fit_gaussian": True},
    "min_intensity": {"log_x": True},
    # denominator is at least 1
    "intensity_rel": {"log_x": True, "fit_gaussian": True},
    # denominator is at least 1
    "paor_mean_intensity": {"log_x": True, "fit_gaussian": True},
    "paor_mean_intensity_rel": {"log_x": True, "fit_gaussian": True},
    "r_xy_um": {"domain": (1e-9, None)},
    "r_z_um": {"domain": (1e-9, None)},
    "r_xy_um_max_std": {"domain": (1e-9, None), "log_x": True},
    "r_z_um_max_std": {"domain": (1e-9, None), "log_x": True},
    "paor_cuboid_um3": {"log_x": True, "domain": (1, None), "fit_gaussian": True},
    "paor_volume_um3": {"log_x": True, "domain": (1, None), "fit_gaussian": True},
    "paor_d1_um5": {"log_x": True, "fit_gaussian": True, "domain": (1, None)},
    "paor_d2_um5": {"log_x": True, "fit_gaussian": True, "domain": (1, None)},
    "paor_d3_um5": {"log_x": True, "fit_gaussian": True, "domain": (1, None)},
}


def gaussian_func(x, a, offset, sigma, c):
    return a * np.exp(-np.square(x - offset) / (2 * sigma**2)) + c


def fit_gaussian_peak(
    data: np.ndarray | None = None,
    hist: np.ndarray | None = None,
    center: np.ndarray | None = None,
    peak_i: int | None = None,
    sigma0: float = 0.1,
) -> tuple[tuple[float, float, float, float], np.ndarray]:
    if hist is None or center is None:
        if data is None:
            raise ValueError("Data must be provided if hist/center isn't provided")

        hist, edges = np.histogram(data, bins="auto", density=True)
        center = (edges[1:] - edges[:-1]) / 2 + edges[:-1]

    if peak_i is None:
        dx = center[1] - center[0]
        prominence = (hist.max() - hist.min()) * 0.02
        peaks, _ = find_peaks(hist, width=2 * sigma0 / dx, prominence=prominence)
        if not len(peaks):
            peak_i = np.argmax(hist)
        else:
            peak_i = peaks[0]

    end = peak_i + np.argmax(hist[peak_i:] < hist[peak_i] * 0.97)

    (a, offset, sigma, c), pcov = curve_fit(
        gaussian_func,
        center[:end],
        hist[:end],
        p0=[hist[peak_i], center[peak_i], sigma0, hist.min()],
    )
    sigma = np.abs(sigma)
    perr = np.sqrt(np.diag(pcov))

    return (a, offset, sigma, c), perr


def find_bi_gaussian_3_points(
    hist: np.ndarray,
    center: np.ndarray,
    width: float,
) -> tuple[int, int, int]:
    # chop off end to remove peak from saturated end
    hist = hist[:-2]
    center = center[:-2]
    dx = center[1] - center[0]
    prominence = (hist.max() - hist.min()) * 0.02

    peaks, _ = find_peaks(hist, width=width / dx, prominence=prominence)
    for i in range(1, 5):
        if len(peaks) >= 2:
            break
        peaks, _ = find_peaks(hist, width=(width / dx) / 2**i, prominence=prominence)

    if len(peaks) > 2:
        peaks = peaks[:2]
    if len(peaks) < 2:
        raise ValueError("could not find at least two peaks")

    a, c = peaks
    peaks, _ = find_peaks(-hist[a:c], width=width / dx, prominence=prominence)
    if len(peaks) != 1:
        raise ValueError("Needed exactly one peak for the trough")
    b = a + peaks[0]

    return a, b, c


def parse_cell_filter(text: str) -> tuple[bool, str, str, Callable[[float, float], bool], str | None, float, float]:
    m = re.match(_cell_filter_pat, text)
    if m is None:
        raise ValueError(f'Unable to parse "{text}"')

    do_log, key, op, stat, num, num2 = m.groups()
    do_log = do_log is not None

    match op:
        case "<":
            op_f = operator.lt
        case "<=":
            op_f = operator.le
        case ">":
            op_f = operator.gt
        case ">=":
            op_f = operator.ge
        case "==":
            op_f = operator.eq
        case "!=":
            op_f = operator.ne
        case _:
            raise ValueError(f'Unable to parse operator "{op}" from "{text}"')

    value = 0
    # if in these 3, there's no number attached
    if stat in ("peak", "mean", "biGaussianTrough"):
        if num:
            raise ValueError(f"{stat} should not have a number attached ({num})")
    else:
        assert stat in ("gaussianS", "gaussian_p", "gaussian_mean", "p", "biGaussianS", None)
        try:
            value = float(num)
        except ValueError as e:
            raise ValueError(f'Could not parse number "{num}" from "{text}"') from e

    value2 = 0
    if stat == "gaussian_p":
        if not num2:
            raise ValueError(f"gaussian_p expects two numbers, (n sigma, p), got just {num}")
        assert num2[0] == ","
        try:
            value2 = float(num2[1:])
        except ValueError as e:
            raise ValueError(f'Could not parse number "{num2[1:]}" from "{text}"') from e
    elif num2:
        raise ValueError(f"Got unexpected second number {num2}")

    return do_log, key, op, op_f, stat, value, value2


def get_metadata_value(cells: list[Cell], key: str) -> np.ndarray:
    match key:
        case "intensity_rel":
            intensity = np.array([cell.metadata["intensity"] for cell in cells])
            min_intensity = np.array([cell.metadata["min_intensity"] for cell in cells])
            return intensity - min_intensity
        case "paor_cuboid_um3":
            extent = np.array([cell.metadata["paor_extent_um"] for cell in cells])
            ex1 = -extent[:, 0, 0] + extent[:, 0, 1]
            ex2 = -extent[:, 1, 0] + extent[:, 1, 1]
            ex3 = -extent[:, 2, 0] + extent[:, 2, 1]
            return ex1 * ex2 * ex3
        case "paor_mean_intensity":
            paor_intensity_total = np.array([cell.metadata["paor_intensity_total"] for cell in cells])
            paor_volume_um3 = np.array([cell.metadata["paor_volume_um3"] for cell in cells])
            return paor_intensity_total / np.maximum(paor_volume_um3, 1)
        case "paor_mean_intensity_rel":
            paor_intensity_total = np.array([cell.metadata["paor_intensity_total"] for cell in cells])
            paor_volume_um3 = np.array([cell.metadata["paor_volume_um3"] for cell in cells])
            min_intensity = np.array([cell.metadata["min_intensity"] for cell in cells])
            return paor_intensity_total / np.maximum(paor_volume_um3, 1) - min_intensity
        case "paor_d1_um5":
            return np.array([cell.metadata["paor_um5"][0] for cell in cells])
        case "paor_d2_um5":
            return np.array([cell.metadata["paor_um5"][1] for cell in cells])
        case "paor_d3_um5":
            return np.array([cell.metadata["paor_um5"][2] for cell in cells])
        case _:
            return np.array([cell.metadata[key] for cell in cells])


def filter_cells(cells: list[Cell], filters: list[str]) -> tuple[list[Cell], list[Cell], dict[str, list[float]]]:
    masks = []
    points = defaultdict(list)
    for do_log, key, op, op_f, stat, stat_value, stat_value2 in map(parse_cell_filter, filters):
        mask = np.ones(len(cells), dtype=bool)
        data = get_metadata_value(cells, key)
        if do_log:
            mask = data > 0
            data[mask] = np.log10(data[mask])
        sigma0 = 0.1 if do_log else 1.25

        match stat:
            case None:
                p_s = ""
            case "p":
                p_s = f", using percentile {stat_value}"
                stat_value = np.percentile(data[mask], stat_value)
            case "peak":
                p_s = ", using the peak"
                stat_value = get_hist_peak(data[mask], domain=MEASURES[key].get("domain"))
            case "mean":
                p_s = ", using the mean"
                stat_value = np.mean(data[mask])
            case "gaussianS":
                (a, offset, sigma, c), perr = fit_gaussian_peak(data=data[mask], sigma0=sigma0)
                stat_value = offset + stat_value * sigma
            case "gaussian_p":
                (a, offset, sigma, c), perr = fit_gaussian_peak(data=data[mask], sigma0=sigma0)
                cut = offset + stat_value * sigma
                stat_value = np.percentile(data[mask][data[mask] >= cut], stat_value2)
            case "gaussian_mean":
                (a, offset, sigma, c), perr = fit_gaussian_peak(data=data[mask], sigma0=sigma0)
                cut = offset + stat_value * sigma
                stat_value = np.mean(data[mask][data[mask] >= cut])
            case "biGaussianS":
                hist, edges = np.histogram(data[mask], bins="auto", density=True)
                center = (edges[1:] - edges[:-1]) / 2 + edges[:-1]
                ai, bi, ci = find_bi_gaussian_3_points(hist, center, 2 * sigma0)

                (a, offset, sigma, c), perr = fit_gaussian_peak(
                    hist=hist[bi:], center=center[bi:], peak_i=ci - bi, sigma0=sigma0
                )
                stat_value = offset + stat_value * sigma
            case "biGaussianTrough":
                hist, edges = np.histogram(data[mask], bins="auto", density=True)
                center = (edges[1:] - edges[:-1]) / 2 + edges[:-1]
                a, b, c = find_bi_gaussian_3_points(hist, center, 2 * sigma0)
                stat_value = center[b]
            case _:
                raise ValueError(f"Unknown stat {stat}")

        for i, c in enumerate(cells):
            if mask[i]:
                mask[i] = op_f(data[i], stat_value)

        removed = len(cells) - np.sum(mask)
        logging.info(f"Keeping cells where {key} {op} {stat_value}{p_s}. Removed {removed} cells")

        masks.append(mask[:, None])
        points[key].append(stat_value)

    mask = np.all(np.concatenate(masks, axis=1), axis=1)
    kept_cells = []
    removed_cells = []
    for cell, keep in zip(cells, mask.tolist(), strict=False):
        (kept_cells if keep else removed_cells).append(cell)

    return kept_cells, removed_cells, points


def get_hist_peak(data: np.ndarray, domain: tuple[float | None, float | None] | None = None) -> float:
    if domain:
        if domain[0] is not None:
            data = data[data >= domain[0]]
        if domain[1] is not None:
            data = data[data <= domain[1]]

    hist, edges = np.histogram(data, bins="auto", density=True)
    center = (edges[1:] + edges[:-1]) / 2
    filtered = median_filter(hist, 3, mode="nearest")
    assert filtered.shape == hist.shape

    peak_idx = np.argmax(filtered)
    peak_x = center[peak_idx].item()

    return peak_x
