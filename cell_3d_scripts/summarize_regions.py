import logging
from argparse import (
    ArgumentDefaultsHelpFormatter,
    ArgumentParser,
)
from datetime import datetime
from pathlib import Path

import fancylog
import numpy as np
from brainglobe_utils.cells.cells import Cell
from brainglobe_utils.IO.cells import get_cells

import cell_3d_scripts
from cell_3d_scripts import __version__
from cell_3d_scripts.atlas import AtlasNode, AtlasTree
from cell_3d_scripts.utils import parse_cell_filter


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
        "--vaa3d-atlas-path",
        dest="vaa3d_atlas_path",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--merged-atlas-path",
        dest="merged_atlas_path",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "--atlas-name",
        dest="atlas_name",
        type=str,
        required=False,
        default="",
    )
    parser.add_argument(
        "-cf",
        "--cell-filter",
        dest="cell_filter",
        type=parse_cell_filter,
        required=False,
        action="append",
        default=[],
    )
    parser.add_argument(
        "-o",
        "--output-path",
        dest="output_path",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def main(
    *,
    cells: list[Cell],
    atlas_tree: AtlasTree,
    output_path: Path,
    atlas_name: str = "",
    merged_atlas_path: Path | None = None,
) -> list[AtlasNode]:
    ts = datetime.now()
    logging.info(
        f"cell_3d_scripts.summarize_regions: Starting regions summary with atlas {atlas_name} for {len(cells)} cells"
    )

    metadata_key = f"region_id_{atlas_name}" if atlas_name else "region_id"
    n = atlas_tree.count_cells(cells, metadata_key)
    logging.info(f"cell_3d_scripts.summarize_regions: found {n} cells outside any region")

    if merged_atlas_path:
        logging.debug(f"Loading merged atlas from {merged_atlas_path}")
        nodes = atlas_tree.read_merge_tree(merged_atlas_path)
    else:
        nodes = list(atlas_tree.root.pre_order())

    atlas_tree.export_nodes_csv(output_path, nodes)

    logging.info(f"cell_3d_scripts.summarize_regions: Analysis took {datetime.now() - ts}")

    return nodes


def run_main():
    args = arg_parser().parse_args()

    output_root = args.output_path.parent
    output_root.mkdir(parents=True, exist_ok=True)

    fancylog.start_logging(
        output_root,
        cell_3d_scripts,
        variables=[
            args,
        ],
        log_header="Cell3DScripts::SummarizeRegions Log",
        multiprocessing_aware=True,
    )

    logging.debug(f"Loading vaa3d format atlas from {args.vaa3d_atlas_path}")
    atlas_tree = AtlasTree.parse_vaa3d(args.vaa3d_atlas_path)
    logging.debug(f"Loading cells from {args.cells_path}")
    cells = get_cells(args.cells_path, cells_only=True)

    for key, op, op_f, percentiles, value in args.cell_filter:
        p_s = ""
        if percentiles:
            p_s = f", using percentile {value}"
            values = [c.metadata[key] for c in cells]
            value = np.percentile(values, value)

        n = len(cells)
        cells = [c for c in cells if op_f(c.metadata[key], value)]

        removed = n - len(cells)
        logging.info(f"Keeping cells where {key} {op} {value}{p_s}. Removed {removed} cells")

    main(
        cells=cells,
        atlas_tree=atlas_tree,
        output_path=args.output_path,
        atlas_name=args.atlas_name,
        merged_atlas_path=args.merged_atlas_path,
    )


if __name__ == "__main__":
    run_main()
