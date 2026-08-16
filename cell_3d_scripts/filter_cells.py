import glob
import logging
from collections import namedtuple
from collections.abc import Sequence
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import click
import fancylog
import matplotlib.pyplot as plt
import numpy as np
from brainglobe_utils.cells.cells import Cell
from brainglobe_utils.IO.cells import get_cells, save_cells
from scipy.optimize import OptimizeWarning

import cell_3d_scripts
from cell_3d_scripts.utils import MEASURES, filter_cells, fit_gaussian_peak, gaussian_func, get_metadata_value


def get_log_bins(data: np.ndarray, n_bins: int) -> np.ndarray:
    data = data[data > 0]
    bins = np.logspace(np.floor(np.log10(np.min(data))), np.ceil(np.log10(np.max(data))), n_bins)
    return bins


def export_cell_metadata_plots(
    root: Path,
    cells: list[Cell],
    measure: str,
    prefix: str = "",
    domain: tuple[float | None, float | None] | None = None,
    log_x: bool = False,
    fit_gaussian: bool = False,
    cutoff_points: Sequence[str] = (),
) -> None:
    if len(cells) < 10:
        logging.warning("Not enough cells for generating plots")
        return

    root.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(nrows=2, ncols=3)
    ax_lin = ax[0, :]  # y is linear
    ax_log = ax[1, :]  # y is log

    data = np.array([get_metadata_value(c, measure) for c in cells])
    if domain:
        if domain[0] is not None:
            data = data[domain[0] <= data]
        if domain[1] is not None:
            data = data[data <= domain[1]]
    if log_x:
        data = np.log10(data[data > 0])

    hist, edges = np.histogram(data, bins="auto", density=True)
    center = (edges[1:] - edges[:-1]) / 2 + edges[:-1]

    # plot with x on linear scale, over full range
    ax_lin[0].plot(center, hist, label="Data")
    ax_log[0].plot(center, hist)
    ax_log[0].set_yscale("log")
    ax_log[0].set_xlabel(f"{'log10 ' if log_x else ''}{measure}")
    ax_lin[0].set_title(f"({np.min(data):3g}, {np.max(data):3g})")

    if fit_gaussian:
        try:
            (a, offset, sigma, c), perr = fit_gaussian_peak(
                hist=hist,
                center=center,
            )
        except (RuntimeError, OptimizeWarning):
            pass
        else:
            pred = gaussian_func(center, a, offset, sigma, c)
            ax_lin[0].plot(center, pred, label="Gaussian")
            ax_lin[0].plot(center, hist - pred, label="Residual")
            ax_lin[0].legend()

    for p in cutoff_points:
        ax_lin[0].plot([p, p], ax_lin[0].get_ylim())

    for i, ql, qh in [(1, 0.01, 0.99), (2, 0.05, 0.95)]:
        ql_val, qh_val = np.quantile(data, [ql, qh])

        def val_to_percentile(x, min_val=ql_val, max_val=qh_val, qh=qh, ql=ql):
            return (x - min_val) / (max_val - min_val) * (qh - ql) + ql

        def percentile_to_val(x, min_val=ql_val, max_val=qh_val, qh=qh, ql=ql):
            return (x - ql) / (qh - ql) * (max_val - min_val) + min_val

        mask = np.logical_and(data >= ql_val, data <= qh_val)
        hist, edges = np.histogram(data[mask], bins="auto", density=True)
        center = (edges[1:] - edges[:-1]) / 2 + edges[:-1]
        ax_lin[i].plot(center, hist)
        ax_log[i].plot(center, hist)
        ax_log[i].set_yscale("log")
        ax_log[i].set_xlim(ql_val, qh_val)
        ax_log[i].set_xlabel(f"{'log10 ' if log_x else ''}{measure}")
        sec_ax = ax_lin[i].secondary_xaxis("top", functions=(val_to_percentile, percentile_to_val))
        sec_ax.set_xlabel("Quantile")
        ax_lin[i].set_title(f"{ql} - {qh} quantiles")

    ax_lin[0].set_ylabel("Density")
    ax_log[0].set_ylabel("Density")

    fig.set_size_inches(18, 6)
    fig.tight_layout()
    fig.savefig(root / f"{prefix}{measure}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main(
    *,
    cells: list[Cell],
    root_path: Path,
    subdir: str,
    cell_filters: list[str] | None = None,
    output_cells_name: str = "",
    output_removed_cells_name: str = "",
    output_plots_name: str = "",
) -> list[Cell]:
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.set_loglevel(level="warning")

    ts = datetime.now()
    logging.info(f"cell_3d_scripts.filter_cells: {subdir} - Starting cell filtering for {len(cells)} cells")

    if output_cells_name or output_removed_cells_name or output_plots_name:
        (root_path / subdir).mkdir(parents=True, exist_ok=True)

    cutoff_points = {}
    kept_cells = cells
    removed_cells = []
    if cell_filters:
        kept_cells, removed_cells, cutoff_points = filter_cells(cells, cell_filters)

    if output_plots_name:
        for measure, options in MEASURES.items():
            export_cell_metadata_plots(
                root_path / subdir / output_plots_name,
                cells,
                measure,
                "pre_",
                **options,
                cutoff_points=cutoff_points.get(measure, ()),
            )

    if output_cells_name:
        save_cells(kept_cells, root_path / subdir / output_cells_name)
    if output_removed_cells_name:
        save_cells(removed_cells, root_path / subdir / output_removed_cells_name)

    if output_plots_name and cell_filters and len(kept_cells) >= 10:
        for measure, options in MEASURES.items():
            options = options.copy()
            options["fit_gaussian"] = False
            export_cell_metadata_plots(
                root_path / subdir / output_plots_name,
                kept_cells,
                measure,
                "post_",
                **options,
            )

    logging.info(f"cell_3d_scripts.filter_cells: Analysis took {datetime.now() - ts}")

    return cells


@click.group(chain=True)
@click.option(
    "--cells-path",
    "-c",
    type=Path,
    required=True,
)
@click.option(
    "--cells-path-glob",
    is_flag=True,
)
@click.option(
    "--root-path",
    "-rp",
    type=Path,
    required=True,
)
@click.option(
    "--output-cells-name",
    type=str,
    required=False,
)
@click.option(
    "--output-removed-cells-name",
    type=str,
    required=False,
)
@click.option(
    "--output-plots-name",
    type=str,
    required=False,
)
@click.pass_context
def run_main(
    ctx,
    cells_path: Path,
    root_path: Path,
    cells_path_glob: bool = False,
    output_cells_name: str = "",
    output_removed_cells_name: str = "",
    output_plots_name: str = "",
):
    root_path.mkdir(parents=True, exist_ok=True)
    Args = namedtuple("Args", ctx.params.keys())

    fancylog.start_logging(
        root_path,
        cell_3d_scripts,
        variables=Args(**ctx.params),
        log_header="Cell3DScripts::FilterCells Log",
        multiprocessing_aware=False,
        filename="cell_3d_scripts-filter_cells",
        timestamp=True,
    )

    logging.debug(f"Loading cells from {cells_path}")
    if cells_path_glob:
        cells = []
        for f in glob.glob(str(cells_path)):
            logging.debug(f"Loading cells from {f}")
            cells.extend(get_cells(f, cells_only=True))
    else:
        cells = get_cells(cells_path, cells_only=True)

    ctx.ensure_object(dict)
    ctx.obj["cells"] = cells
    ctx.obj["root_path"] = root_path
    ctx.obj["output_cells_name"] = output_cells_name
    ctx.obj["output_removed_cells_name"] = output_removed_cells_name
    ctx.obj["output_plots_name"] = output_plots_name


@run_main.command()
@click.option(
    "--subdir",
    type=str,
    required=True,
)
@click.option(
    "--cell-filters",
    "-cf",
    type=str,
    required=False,
    multiple=True,
)
@click.pass_context
def apply_filter(
    ctx,
    subdir: str,
    cell_filters: list[str] | None = None,
):
    main(
        cells=deepcopy(ctx.obj["cells"]),
        cell_filters=cell_filters,
        root_path=ctx.obj["root_path"],
        output_cells_name=ctx.obj["output_cells_name"],
        output_removed_cells_name=ctx.obj["output_removed_cells_name"],
        subdir=subdir,
        output_plots_name=ctx.obj["output_plots_name"],
    )


if __name__ == "__main__":
    run_main()
