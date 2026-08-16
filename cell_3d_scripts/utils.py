import logging
import operator
import re
from collections import defaultdict
from collections.abc import Callable

import numpy as np
from brainglobe_utils.cells.cells import Cell
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit

_cell_filter_pat = re.compile(
    r"^(log10:)?"
    r"([\w.]+)"
    r"([<>=!]{1,2})"
    r"(p|peak|mean|gaussian|gaussian_p|gaussian_mean)?"
    r"(-?[0-9]*\.?[0-9]*)"
    r"(,[0-9]*\.?[0-9]*)?$"
)

# region, outside of which, it's definitely a bad value
MEASURES = {
    "intensity": {"log_x": True, "fit_gaussian": True},
    "min_intensity": {"log_x": True},
    # denominator is at least 1
    "intensity_ratio": {"log_x": True, "fit_gaussian": True},
    # denominator is at least 1
    "paor_mean_intensity": {"log_x": True, "fit_gaussian": True},
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
    sigma0: float = 0.1,
) -> tuple[tuple[float, float, float, float], np.ndarray]:
    if hist is None or center is None:
        if data is None:
            raise ValueError("Data must be provided if hist/center isn't provided")

        hist, edges = np.histogram(data, bins="auto", density=True)
        center = (edges[1:] - edges[:-1]) / 2 + edges[:-1]

    max_hist_i = np.argmax(hist)
    end = max_hist_i + np.argmax(hist[max_hist_i:] < hist[max_hist_i] * 0.95)

    (a, offset, sigma, c), pcov = curve_fit(
        gaussian_func,
        center[:end],
        hist[:end],
        p0=[hist[max_hist_i], center[max_hist_i], sigma0, hist.min()],
    )
    sigma = np.abs(sigma)
    perr = np.sqrt(np.diag(pcov))

    return (a, offset, sigma, c), perr


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
    if stat in ("peak", "mean"):
        if num:
            raise ValueError(f"{stat} should not have a number attached ({num})")
    else:
        assert stat in ("gaussian", "gaussian_p", "gaussian_mean", "p", None)
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


def get_metadata_value(cell: Cell, key: str) -> float:
    match key:
        case "intensity_ratio":
            return cell.metadata["intensity"] / max(cell.metadata["min_intensity"], 1)
        case "paor_cuboid_um3":
            (ex11, ex12), (ex21, ex22), (ex31, ex32) = cell.metadata["paor_extent_um"]
            return (-ex11 + ex12) * (-ex21 + ex22) * (-ex31 + ex32)
        case "paor_mean_intensity":
            return cell.metadata["paor_intensity_total"] / max(cell.metadata["paor_volume_um3"], 1)
        case "paor_d1_um5":
            return cell.metadata["paor_um5"][0]
        case "paor_d2_um5":
            return cell.metadata["paor_um5"][1]
        case "paor_d3_um5":
            return cell.metadata["paor_um5"][2]
        case _:
            return cell.metadata[key]


def filter_cells(cells: list[Cell], filters: list[str]) -> tuple[list[Cell], list[Cell], dict[str, list[float]]]:
    masks = []
    points = defaultdict(list)
    for do_log, key, op, op_f, stat, stat_value, stat_value2 in map(parse_cell_filter, filters):
        mask = np.ones(len(cells), dtype=bool)
        data = np.array([get_metadata_value(c, key) for c in cells])
        if do_log:
            mask = data > 0
            data[mask] = np.log10(data[mask])

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
            case "gaussian":
                (a, offset, sigma, c), perr = fit_gaussian_peak(data=data[mask])
                stat_value = offset + stat_value * sigma
            case "gaussian_p":
                (a, offset, sigma, c), perr = fit_gaussian_peak(data=data[mask])
                cut = offset + stat_value * sigma
                stat_value = np.percentile(data[mask][data[mask] >= cut], stat_value2)
            case "gaussian_mean":
                (a, offset, sigma, c), perr = fit_gaussian_peak(data=data[mask])
                cut = offset + stat_value * sigma
                stat_value = np.mean(data[mask][data[mask] >= cut])
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
