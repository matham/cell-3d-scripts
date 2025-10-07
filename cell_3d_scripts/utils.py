import logging
import operator
import re
from collections.abc import Callable

import numpy as np
from brainglobe_utils.cells.cells import Cell

_cell_filter_pat = re.compile(r"([\w.]+)([<>=!]{1,2})(p)?([0-9]*.?[0-9]*)")


def parse_cell_filter(text: str) -> tuple[str, str, Callable[[float, float], bool], bool, float]:
    m = re.match(_cell_filter_pat, text)
    if m is None:
        raise ValueError(f'Unable to parse "{text}"')

    key, op, p, num = m.groups()

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

    percentiles = p is not None

    try:
        value = float(num)
    except ValueError as e:
        raise ValueError(f'Could not parse number "{num}" from "{text}"') from e

    return key, op, op_f, percentiles, value


def filter_cells(cells: list[Cell], filters: list[str]) -> list[Cell]:
    for key, op, op_f, percentiles, value in map(parse_cell_filter, filters):
        p_s = ""
        if percentiles:
            p_s = f", using percentile {value}"
            values = [c.metadata[key] for c in cells]
            value = np.percentile(values, value)

        n = len(cells)
        cells = [c for c in cells if op_f(c.metadata[key], value)]

        removed = n - len(cells)
        logging.info(f"Keeping cells where {key} {op} {value}{p_s}. Removed {removed} cells")

    return cells
