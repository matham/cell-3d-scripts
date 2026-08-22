import glob
import logging
from collections import namedtuple
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from functools import partial
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
    cells: list[Cell] | list[dict],
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

    data = get_metadata_value(cells, measure)
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
    cells: list[Cell] | list[dict],
    root_path: Path,
    subdir: str,
    cell_filters: list[str] | None = None,
    output_cells_name: str = "",
    output_removed_cells_name: str = "",
    output_plots_name: str = "",
    conserve_memory: bool = False,
) -> list[Cell] | list[dict]:
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.set_loglevel(level="warning")

    ts = datetime.now()
    logging.info(f"cell_3d_scripts.filter_cells: {subdir} - Starting cell filtering for {len(cells)} cells")

    if conserve_memory and (output_cells_name or output_removed_cells_name):
        raise ValueError("Can't save output cells if conserving memory")

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


def worker_get_cells(f: Path | str, conserve_memory: bool) -> list[Cell] | list[dict]:
    cells = get_cells(f, cells_only=True)

    if conserve_memory:
        used = set(MEASURES)
        used.update(
            {
                "intensity",
                "min_intensity",
                "paor_extent_um",
                "paor_intensity_total",
                "paor_volume_um3",
                "paor_um5",
                "intensity_factor",
            }
        )
        return [{k: v for k, v in c.metadata.items() if k in used} for c in cells]
    return cells


@click.group(chain=True)
@click.option(
    "--cells-path",
    "-c",
    type=Path,
    required=True,
    multiple=True,
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
@click.option(
    "--output-raw-name",
    type=str,
    required=False,
)
@click.option(
    "--conserve-memory",
    is_flag=True,
)
@click.pass_context
def cli_main(
    ctx,
    cells_path: list[Path],
    root_path: Path,
    cells_path_glob: bool = False,
    output_cells_name: str = "",
    output_removed_cells_name: str = "",
    output_plots_name: str = "",
    output_raw_name: str = "",
    conserve_memory: bool = False,
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

    cells = []
    if cells_path_glob:
        for p in cells_path:
            logging.debug(f"Loading cells from {p}")
            for f in glob.glob(str(p), recursive=True):
                logging.debug(f"Loading cells from {f}")
                with ProcessPoolExecutor(max_workers=1) as executor:
                    future = executor.submit(worker_get_cells, f, conserve_memory)
                    cells.extend(future.result())
    else:
        for p in cells_path:
            logging.debug(f"Loading cells from {p}")
            with ProcessPoolExecutor(max_workers=1) as executor:
                future = executor.submit(worker_get_cells, p, conserve_memory)
                cells.extend(future.result())

    ctx.ensure_object(dict)
    ctx.obj["cells"] = cells
    ctx.obj["root_path"] = root_path
    ctx.obj["output_cells_name"] = output_cells_name
    ctx.obj["output_removed_cells_name"] = output_removed_cells_name
    ctx.obj["output_plots_name"] = output_plots_name
    ctx.obj["conserve_memory"] = conserve_memory

    if output_raw_name:
        data = {}
        for measure in MEASURES:
            data[measure] = get_metadata_value(cells, measure)

        logging.debug(f"Saving raw data to {root_path / output_raw_name}")
        np.savez(root_path / output_raw_name, **data)


@cli_main.command()
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
        cells=ctx.obj["cells"],
        cell_filters=cell_filters,
        root_path=ctx.obj["root_path"],
        output_cells_name=ctx.obj["output_cells_name"],
        output_removed_cells_name=ctx.obj["output_removed_cells_name"],
        subdir=subdir,
        output_plots_name=ctx.obj["output_plots_name"],
        conserve_memory=ctx.obj["conserve_memory"],
    )


run_main = partial(cli_main, windows_expand_args=False)

if __name__ == "__main__":
    run_main()
