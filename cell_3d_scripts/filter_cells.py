import logging
from argparse import (
    ArgumentDefaultsHelpFormatter,
    ArgumentParser,
)
from datetime import datetime
from pathlib import Path

import fancylog
from brainglobe_utils.cells.cells import Cell
from brainglobe_utils.IO.cells import get_cells, save_cells

import cell_3d_scripts
from cell_3d_scripts import __version__
from cell_3d_scripts.utils import filter_cells


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
        required=True,
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


def export_cell_metadata_plots(root: Path, cells: list[Cell]) -> None:
    pass


def main(
    *,
    cells: list[Cell],
    output_cells_path: Path | None = None,
    cell_filters: list[str] | None = None,
    output_plots_path: Path | None = None,
) -> list[Cell]:
    ts = datetime.now()
    logging.info(f"cell_3d_scripts.filter_cells: Starting cell filtering for {len(cells)} cells")

    if output_plots_path:
        export_cell_metadata_plots(output_plots_path / "pre_filter", cells)

    if cell_filters:
        cells = filter_cells(cells, cell_filters)

    if output_cells_path:
        output_cells_path.parent.mkdir(parents=True, exist_ok=True)
        save_cells(cells, output_cells_path)

    if output_plots_path:
        export_cell_metadata_plots(output_plots_path / "post_filter", cells)

    logging.info(f"cell_3d_scripts.filter_cells: Analysis took {datetime.now() - ts}")

    return cells


def run_main():
    args = arg_parser().parse_args()

    if args.output_cells_path:
        output_root = args.output_cells_path.parent
    elif args.output_plots_path:
        output_root = args.output_plots_path.parent
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
        output_plots_path=args.output_plots_path,
    )


if __name__ == "__main__":
    run_main()
