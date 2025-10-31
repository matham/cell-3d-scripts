import logging
import multiprocessing as mp
import tempfile
import zipfile
from argparse import (
    ArgumentDefaultsHelpFormatter,
    ArgumentParser,
)
from copy import deepcopy
from datetime import datetime
from functools import partial
from pathlib import Path

import fancylog
import numpy as np
import tqdm
import zarr
from brainglobe_utils.cells.cells import Cell
from brainglobe_utils.IO.cells import get_cells, save_cells
from brainglobe_utils.IO.image import read_z_stack

import cell_3d_scripts
from cell_3d_scripts import __version__
from cell_3d_scripts.atlas import OUTSIDE_REGION_ID, AtlasTree

_region_cache: tuple[int, np.ndarray] | None = None


def arg_parser() -> ArgumentParser:
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-c",
        "--cells-path",
        dest="cells_path",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-r",
        "--region-id-path",
        dest="region_id_path",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--vaa3d-atlas-path",
        dest="vaa3d_atlas_path",
        type=str,
        required=False,
        default="",
    )
    parser.add_argument(
        "--atlas-name",
        dest="atlas_name",
        type=str,
        required=False,
        default="",
    )
    parser.add_argument(
        "-o",
        "--output-cells-path",
        dest="output_cells_path",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--workers",
        dest="num_workers",
        type=int,
        required=False,
        default=6,
    )

    return parser


def _get_region_from_pos(z: float, y: float, x: float, region_ids: np.ndarray, tree: AtlasTree | None) -> int:
    global _region_cache
    z = int(round(z))
    y = int(round(y))
    x = int(round(x))

    shape = region_ids.shape
    if z >= shape[0] or y >= shape[1] or x >= shape[2]:
        region_id = OUTSIDE_REGION_ID
    else:
        if _region_cache is None or _region_cache[0] != z:
            _region_cache = z, np.asarray(region_ids[z, ...])
        region_id = _region_cache[1][y, x].item()

    if tree is not None and region_id != OUTSIDE_REGION_ID:
        # if we have a region that is a potential parent, get its leaf region ID
        node = tree.get_leaf_node_from_canonical_id(region_id)
        region_id = node.region_id

    return region_id


def _look_up_cell(arg, region_ids: np.ndarray, tree: AtlasTree | None) -> tuple[int, int]:
    cell_i, pos = arg
    region_id = _get_region_from_pos(*pos, region_ids=region_ids, tree=tree)

    return cell_i, region_id


def main(
    *,
    cells: list[Cell],
    region_ids: np.ndarray,
    atlas_name: str = "",
    atlas_tree: AtlasTree | None = None,
    output_cells_path: Path | None = None,
    num_workers: int = 6,
) -> list[Cell]:
    ts = datetime.now()
    regions_voxels = region_ids.shape
    logging.info(
        f"cell_3d_scripts.locate_cells: Starting region ID lookup for atlas {atlas_name} for {len(cells)} cells"
    )
    logging.debug(f"Tiff containing the region IDs' shape is {regions_voxels}")

    cells = sorted(cells, key=lambda c: (c.z, c.x, c.y))
    output_cells = deepcopy(cells)
    metadata_key = f"region_id_{atlas_name}" if atlas_name else "region_id"

    if num_workers:
        cell_pos = [(cell.z, cell.y, cell.x) for cell in output_cells]

        progress_bar = tqdm.tqdm(total=len(output_cells), unit="cells")
        f = partial(_look_up_cell, region_ids=region_ids, tree=atlas_tree)
        ctx = mp.get_context("spawn")
        # we can't use the context manager because of coverage issues:
        pool = ctx.Pool(processes=num_workers)
        try:
            for cell_i, region_id in pool.imap(f, list(enumerate(cell_pos)), chunksize=1000):
                output_cells[cell_i].metadata[metadata_key] = region_id
                progress_bar.update()
        finally:
            pool.close()
            pool.join()
        progress_bar.close()
    else:
        for cell in tqdm.tqdm(output_cells, unit="cells"):
            cell.metadata[metadata_key] = _get_region_from_pos(cell.z, cell.y, cell.x, region_ids, atlas_tree)

    if output_cells_path:
        save_cells(output_cells, str(output_cells_path))
    logging.info(f"cell_3d_scripts.locate_cells: Analysis took {datetime.now() - ts}")

    return output_cells


def run_main():
    args = arg_parser().parse_args()

    output_cells_path = Path(args.cells_path)
    output_root = Path(args.cells_path).parent
    if args.output_cells_path:
        output_cells_path = Path(args.output_cells_path)
        output_cells_path.parent.mkdir(parents=True, exist_ok=True)
        output_root = output_cells_path.parent

    atlas_tree = None
    if args.vaa3d_atlas_path:
        atlas_tree = AtlasTree.parse_vaa3d(args.vaa3d_atlas_path)

    fancylog.start_logging(
        output_root,
        cell_3d_scripts,
        variables=[
            args,
        ],
        log_header="Cell3DScripts::LocalizeCells Log",
        multiprocessing_aware=True,
        filename="cell_3d_scripts-locate_cells",
        timestamp=True,
    )

    logging.debug(f"Loading cells from {args.cells_path}")
    cells = get_cells(args.cells_path, cells_only=True)

    region_id_path = Path(args.region_id_path)
    atlas_name = args.atlas_name

    logging.debug(f"Using region IDs from {region_id_path}")
    tempdir = None
    if region_id_path.suffix == ".zarr":
        group = zarr.open(str(region_id_path), mode="r")
        arr = group.get("data")
    elif region_id_path.name.endswith(".zarr.zip"):
        tempdir = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        with zipfile.ZipFile(region_id_path, "r") as zf:
            zf.extractall(tempdir.name)
        group = zarr.open(tempdir.name, mode="r")
        arr = group.get("data")
    else:
        arr = read_z_stack(str(region_id_path))

    try:
        main(
            cells=cells,
            region_ids=arr,
            atlas_name=atlas_name,
            atlas_tree=atlas_tree,
            output_cells_path=output_cells_path,
            num_workers=args.num_workers,
        )
    finally:
        if tempdir is not None:
            tempdir.cleanup()


if __name__ == "__main__":
    run_main()
