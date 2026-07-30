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

_region_cache: dict[int, tuple[int, np.ndarray]] = {}


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
        "--region-id-img-path",
        dest="region_id_img_path",
        nargs="+",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--vaa3d-atlas-path",
        dest="vaa3d_atlas_path",
        nargs="*",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--atlas-name",
        dest="atlas_name",
        nargs="+",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-cells-path",
        dest="output_cells_path",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--exclude-cells",
        dest="exclude_cells",
        action="store_true",
    )
    parser.add_argument(
        "--exclude-non-cells",
        dest="exclude_non_cells",
        action="store_true",
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


def _get_region_from_pos(
    z: float,
    y: float,
    x: float,
    atlases_region_ids: list[np.ndarray],
    atlases_tree: list[AtlasTree | None],
) -> list[int]:
    z = int(round(z))
    y = int(round(y))
    x = int(round(x))

    region_ids = []
    for i, (atlas_region_ids, atlas_tree) in enumerate(zip(atlases_region_ids, atlases_tree, strict=False)):
        shape = atlas_region_ids.shape
        if z >= shape[0] or y >= shape[1] or x >= shape[2]:
            region_id = OUTSIDE_REGION_ID
        else:
            if i not in _region_cache or _region_cache[i][0] != z:
                _region_cache[i] = z, np.asarray(atlas_region_ids[z, ...])
            region_id = _region_cache[i][1][y, x].item()

        if atlas_tree is not None and region_id != OUTSIDE_REGION_ID:
            # if we have a region that is a potential parent, get its leaf region ID
            node = atlas_tree.get_leaf_node_from_canonical_id(region_id)
            region_id = node.region_id

        region_ids.append(region_id)

    return region_ids


def _look_up_cell(
    arg,
    atlases_region_ids: list[np.ndarray],
    atlases_tree: list[AtlasTree | None],
) -> tuple[int, list[int]]:
    cell_i, pos = arg
    region_ids = _get_region_from_pos(*pos, atlases_region_ids=atlases_region_ids, atlases_tree=atlases_tree)

    return cell_i, region_ids


def main(
    *,
    cells: list[Cell],
    atlases_region_ids: list[np.ndarray],
    atlases_name: list[str],
    atlases_tree: list[AtlasTree | None],
    output_cells_path: Path | None = None,
    num_workers: int = 6,
) -> list[Cell]:
    global _region_cache
    _region_cache = {}

    ts = datetime.now()
    regions_voxels = [r.shape for r in atlases_region_ids]
    logging.info(f"cell_3d_scripts.locate_cells: Starting region ID lookup for {len(cells)} cells")
    logging.debug(f"Tiffs containing the region IDs' shape for atlases {atlases_name} is {regions_voxels}")

    cells = sorted(cells, key=lambda c: (c.z, c.x, c.y))
    output_cells = deepcopy(cells)
    metadata_keys = [f"region_id_{name}" for name in atlases_name]

    if num_workers:
        cell_pos = [(cell.z, cell.y, cell.x) for cell in output_cells]

        progress_bar = tqdm.tqdm(total=len(output_cells), unit="cells")
        f = partial(_look_up_cell, atlases_region_ids=atlases_region_ids, atlases_tree=atlases_tree)
        ctx = mp.get_context("spawn")
        # we can't use the context manager because of coverage issues:
        pool = ctx.Pool(processes=num_workers)
        try:
            for cell_i, ids in pool.imap(f, list(enumerate(cell_pos)), chunksize=1000):
                cell = output_cells[cell_i]
                assert len(ids) == len(metadata_keys)
                for key, id_ in zip(metadata_keys, ids, strict=False):
                    cell.metadata[key] = id_
                progress_bar.update()
        finally:
            pool.close()
            pool.join()
        progress_bar.close()
    else:
        for cell in tqdm.tqdm(output_cells, unit="cells"):
            ids = _get_region_from_pos(cell.z, cell.y, cell.x, atlases_region_ids, atlases_tree)
            assert len(ids) == len(metadata_keys)
            for key, id_ in zip(metadata_keys, ids, strict=False):
                cell.metadata[key] = id_

    if output_cells_path:
        save_cells(output_cells, str(output_cells_path))
    logging.info(f"cell_3d_scripts.locate_cells: Analysis took {datetime.now() - ts}")

    _region_cache = {}

    return output_cells


def run_main():
    args = arg_parser().parse_args()

    output_cells_path = Path(args.cells_path)
    output_root = Path(args.cells_path).parent
    if args.output_cells_path:
        output_cells_path = Path(args.output_cells_path)
        output_cells_path.parent.mkdir(parents=True, exist_ok=True)
        output_root = output_cells_path.parent

    atlas_tree = [AtlasTree.parse_vaa3d(p) for p in args.vaa3d_atlas_path]

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
    if args.exclude_cells and args.exclude_non_cells:
        raise ValueError("Excluding both cells and non-cells")
    if args.exclude_non_cells:
        cells = get_cells(args.cells_path, cells_only=True)
    elif args.exclude_cells:
        cells = get_cells(args.cells_path, cells_only=False)
        cells = [c for c in cells if not c.is_cell()]
    else:
        cells = get_cells(args.cells_path, cells_only=False)

    region_id_img_path = [Path(p) for p in args.region_id_img_path]
    atlas_name = args.atlas_name
    if len(region_id_img_path) != len(atlas_name):
        raise ValueError("Must provide same number of atlas' name as region id image paths")
    if atlas_tree and len(region_id_img_path) != len(atlas_tree):
        raise ValueError("Must provide same number of atlases names as atlases")

    region_ids = []
    tempdirs = []
    for atlas_region_id_path in region_id_img_path:
        logging.debug(f"Using region IDs from {atlas_region_id_path}")
        if atlas_region_id_path.suffix == ".zarr":
            group = zarr.open(str(atlas_region_id_path), mode="r")
            arr = group.get("data")
        elif atlas_region_id_path.name.endswith(".zarr.zip"):
            tempdir = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
            tempdirs.append(tempdir)
            with zipfile.ZipFile(atlas_region_id_path, "r") as zf:
                zf.extractall(tempdir.name)
            group = zarr.open(tempdir.name, mode="r")
            arr = group.get("data")
        else:
            arr = read_z_stack(str(atlas_region_id_path))

        region_ids.append(arr)

    if not atlas_tree:
        atlas_tree = [
            None,
        ] * len(region_ids)
    try:
        main(
            cells=cells,
            atlases_region_ids=region_ids,
            atlases_name=atlas_name,
            atlases_tree=atlas_tree,
            output_cells_path=output_cells_path,
            num_workers=args.num_workers,
        )
    finally:
        for tempdir in tempdirs:
            tempdir.cleanup()


if __name__ == "__main__":
    run_main()
