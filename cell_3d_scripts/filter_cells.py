import logging
from argparse import (
    ArgumentDefaultsHelpFormatter,
    ArgumentParser,
)
from datetime import datetime
from pathlib import Path

import fancylog
import matplotlib.pyplot as plt
import numpy as np
from brainglobe_utils.cells.cells import Cell
from brainglobe_utils.IO.cells import get_cells, save_cells

import cell_3d_scripts
from cell_3d_scripts import __version__
from cell_3d_scripts.utils import MEASURES, filter_cells, get_hist_peak


def arg_parser() -> ArgumentParser:
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-c",
        "--cells-path",
        dest="cells_path",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-cf",
        "--cell-filter",
        dest="cell_filter",
        type=str,
        required=False,
        action="append",
        default=[],
    )
    parser.add_argument(
        "-o",
        "--output-cells-path",
        dest="output_cells_path",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "--output-removed-cells-path",
        dest="output_removed_cells_path",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-plots",
        "--output-plots-path",
        dest="output_plots_path",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def get_log_bins(data: np.ndarray, n_bins: int) -> np.ndarray:
    data = data[data > 0]
    bins = np.logspace(np.floor(np.log10(np.min(data))), np.ceil(np.log10(np.max(data))), n_bins)
    return bins


def export_cell_metadata_plots(
    root: Path,
    cells: list[Cell],
    measure: str,
    prefix: str = "",
    domain: tuple[float, float] | None = None,
) -> None:
    root.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(nrows=2, ncols=4, sharex="col")
    ax_lin = ax[0, :]
    ax_log = ax[1, :]

    values = np.array([c.metadata[measure] for c in cells])
    n_bins = max(25, min(500, len(values) // 1000))

    peak_x = get_hist_peak(values, domain)

    # plot with x on log scale, over full range
    bins = get_log_bins(values, n_bins)
    ax_lin[0].hist(values, bins=bins, density=True)
    ax_log[0].hist(values, bins=bins, density=True, log=True)
    ax_log[0].set_xscale("log")
    ax_log[0].set_xlabel(measure)
    ax_lin[0].set_title("Log x (>0)")

    # plot with x on linear scale, over full range
    ax_lin[1].hist(values, bins=n_bins, density=True)
    ax_lin[1].plot([peak_x], [ax_lin[1].get_ylim()[1]], "*r")
    ax_log[1].hist(values, bins=n_bins, density=True, log=True)
    ax_log[1].set_xlabel(measure)
    ax_lin[1].set_title(f"({np.min(values):3g}, {peak_x:3g}*, {np.max(values):3g})")

    for i, ql, qh in [(2, 0.01, 0.99), (3, 0.05, 0.95)]:
        ql_val, qh_val = np.quantile(values, [ql, qh])

        def val_to_percentile(x, min_val=ql_val, max_val=qh_val, qh=qh, ql=ql):
            return (x - min_val) / (max_val - min_val) * (qh - ql) + ql

        def percentile_to_val(x, min_val=ql_val, max_val=qh_val, qh=qh, ql=ql):
            return (x - ql) / (qh - ql) * (max_val - min_val) + min_val

        mask = np.logical_and(values >= ql_val, values <= qh_val)
        ax_lin[i].hist(values[mask], bins=n_bins, density=True)
        ax_log[i].hist(values[mask], bins=n_bins, density=True, log=True)
        ax_log[i].set_xlim(ql_val, qh_val)
        ax_log[i].set_xlabel(measure)
        ax_lin[i].plot([peak_x], [ax_lin[i].get_ylim()[1]], "*r")
        sec_ax = ax_lin[i].secondary_xaxis("top", functions=(val_to_percentile, percentile_to_val))
        sec_ax.set_xlabel("Quantile")
        ax_lin[i].set_title(f"{ql} - {qh} quantiles")

    ax_lin[0].set_ylabel("Density")
    ax_log[0].set_ylabel("Density")

    fig.set_size_inches(18, 6)
    fig.tight_layout()
    fig.savefig(root / f"{prefix}{measure}.png", dpi=600, bbox_inches="tight")


def main(
    *,
    cells: list[Cell],
    output_cells_path: Path | None = None,
    cell_filters: list[str] | None = None,
    output_removed_cells_path: Path | None = None,
    output_plots_path: Path | None = None,
) -> list[Cell]:
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.set_loglevel(level="warning")

    ts = datetime.now()
    logging.info(f"cell_3d_scripts.filter_cells: Starting cell filtering for {len(cells)} cells")

    if output_plots_path:
        for measure, domain in MEASURES.items():
            export_cell_metadata_plots(output_plots_path, cells, measure, "pre_", domain)

    removed_cells = []
    if cell_filters:
        cells, removed_cells = filter_cells(cells, cell_filters)

    if output_cells_path:
        output_cells_path.parent.mkdir(parents=True, exist_ok=True)
        save_cells(cells, output_cells_path)
    if output_removed_cells_path:
        output_removed_cells_path.parent.mkdir(parents=True, exist_ok=True)
        save_cells(removed_cells, output_removed_cells_path)

    if output_plots_path:
        for measure, domain in MEASURES.items():
            export_cell_metadata_plots(output_plots_path, cells, measure, "post_", domain)

    logging.info(f"cell_3d_scripts.filter_cells: Analysis took {datetime.now() - ts}")

    return cells


def run_main():
    args = arg_parser().parse_args()

    if args.output_cells_path:
        output_root = args.output_cells_path.parent
    elif args.output_removed_cells_path:
        output_root = args.output_removed_cells_path.parent
    elif args.output_plots_path:
        output_root = args.output_plots_path
    else:
        output_root = args.cells_path.parent

    fancylog.start_logging(
        output_root,
        cell_3d_scripts,
        variables=[
            args,
        ],
        log_header="Cell3DScripts::FilterCells Log",
        multiprocessing_aware=True,
    )

    logging.debug(f"Loading cells from {args.cells_path}")
    cells = get_cells(args.cells_path, cells_only=True)

    main(
        cells=cells,
        cell_filters=args.cell_filter,
        output_cells_path=args.output_cells_path,
        output_removed_cells_path=args.output_removed_cells_path,
        output_plots_path=args.output_plots_path,
    )


if __name__ == "__main__":
    run_main()
