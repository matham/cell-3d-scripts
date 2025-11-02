import logging
import operator
import re
from collections.abc import Callable

import numpy as np
import scipy.stats
from brainglobe_utils.cells.cells import Cell
from scipy.ndimage import median_filter

_cell_filter_pat = re.compile(r"^([\w.]+)([<>=!]{1,2})(p|peak|mean)?(-?[0-9]*.?[0-9]*)$")

# region, outside of which, it's definitely a bad value
MEASURES = {
    "center_intensity": None,
    "r_xy": (1e-9, 1e4),
    "r_z": (1e-9, 1e4),
    "r_xy_max_std": (1e-9, 100),
    "r_z_max_std": (1e-9, 100),
}


def parse_cell_filter(text: str) -> tuple[str, str, Callable[[float, float], bool], str | None, float]:
    m = re.match(_cell_filter_pat, text)
    if m is None:
        raise ValueError(f'Unable to parse "{text}"')

    key, op, measure, num = m.groups()

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
    if measure not in ("peak", "mean"):
        assert measure in ("p", None)
        try:
            value = float(num)
        except ValueError as e:
            raise ValueError(f'Could not parse number "{num}" from "{text}"') from e

    return key, op, op_f, measure, value


def filter_cells(cells: list[Cell], filters: list[str]) -> tuple[list[Cell], list[Cell]]:
    removed_cells = []
    for key, op, op_f, measure, value in map(parse_cell_filter, filters):
        p_s = ""
        if measure == "p":
            p_s = f", using percentile {value}"
            values = [c.metadata[key] for c in cells]
            value = np.percentile(values, value)
        elif measure == "peak":
            p_s = ", using the peak"
            value = get_hist_peak(
                np.array([c.metadata[key] for c in cells]),
                domain=MEASURES.get(key),
            )
        elif measure == "mean":
            p_s = ", using the mean"
            value = np.mean([c.metadata[key] for c in cells])

        n = len(cells)
        temp = []
        for c in cells:
            if op_f(c.metadata[key], value):
                temp.append(c)
            else:
                removed_cells.append(c)
        cells = temp

        removed = n - len(cells)
        logging.info(f"Keeping cells where {key} {op} {value}{p_s}. Removed {removed} cells")

    return cells, removed_cells


def get_hist_peak(values: np.ndarray, domain: tuple[float, float] | None = None) -> float:
    ql, qh = np.quantile(values, [0.01, 0.99])
    if domain:
        ql = max(ql, domain[0])
        qh = min(qh, domain[1])
    mask = np.logical_and(values >= ql, values <= qh)
    values = values[mask]

    bins = max(25, min(1000, len(values) // 1000))
    hist, edges = np.histogram(values, bins=bins, density=True)
    center = (edges[1:] + edges[:-1]) / 2
    filtered = median_filter(hist, 3, mode="nearest")
    assert filtered.shape == hist.shape

    peak_idx = np.argmax(filtered)
    peak_x = center[peak_idx].item()

    return peak_x


def get_invgamma_dist(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = values.astype(np.float64)
    min_val = np.min(values)
    norm_values = values - min_val
    max_val = np.max(norm_values)
    norm_values /= max_val

    a, loc, scale = scipy.stats.invgamma.fit(norm_values)

    x_norm = np.linspace(0, 1, 500)
    y = scipy.stats.invgamma.pdf(x_norm, a, loc=loc, scale=scale)
    x = x_norm * max_val + min_val

    return x, y
