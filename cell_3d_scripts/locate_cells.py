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


def _scale_point(
    point: tuple[float, float, float], input_size: tuple[int, int, int], output_size: tuple[int, int, int]
) -> tuple[int, int, int]:
    output = []
    for p, i, o in zip(point, input_size, output_size, strict=False):
        # last value is i - 1
        v = p / (i - 1) * (o - 1)
        v = int(round(v))
        output.append(v)

    return tuple(output)


def arg_parser() -> ArgumentParser:
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-s",
        "--signal-planes-path",
        dest="signal_planes_path",
        type=str,
        default="",
        required=False,
        help="Path to the directory of the signal files. Can also be a text"
        "file pointing to the files. For a 3d tiff, data is in z, y, x order",
    )
    parser.add_argument(
        "--data-size",
        dest="signal_data_size",
        type=int,
        nargs=3,
        default=(),
        required=False,
    )
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


def _look_up_cell(arg, region_ids: np.ndarray, tree: AtlasTree | None) -> tuple[int, int]:
    cell_i, pos = arg
    region_id = region_ids[pos].item()

    if tree is not None and region_id != OUTSIDE_REGION_ID:
        node = tree.get_node_from_canonical_id(region_id)
        region_id = node.region_id

    return cell_i, region_id


def main(
    *,
    data_size_voxels: tuple[int, int, int],
    cells: list[Cell],
    region_ids: np.ndarray,
    atlas_name: str = "",
    atlas_tree: AtlasTree | None = None,
    output_cells_path: Path | None = None,
    num_workers: int = 6,
) -> list[Cell]:
    ts = datetime.now()
    regions_voxels = region_ids.shape
    logging.info(f"cell_meta_3d: Starting region ID lookup for atlas {atlas_name} for {len(cells)} cells")
    logging.debug(f"Region IDs shape is {regions_voxels}. Data shape is {data_size_voxels}")

    cells = sorted(cells, key=lambda c: (c.z, c.x, c.y))
    output_cells = deepcopy(cells)
    metadata_key = f"region_id_{atlas_name}" if atlas_name else "region_id"

    if num_workers:
        cell_pos = [_scale_point((cell.z, cell.y, cell.x), data_size_voxels, regions_voxels) for cell in output_cells]

        progress_bar = tqdm.tqdm(total=len(output_cells), unit="cells")
        f = partial(_look_up_cell, region_ids=region_ids, tree=atlas_tree)
        ctx = mp.get_context("spawn")
        # we can't use the context manager because of coverage issues:
        pool = ctx.Pool(processes=num_workers)
        try:
            for cell_i, region_id in pool.imap_unordered(f, list(enumerate(cell_pos))):
                output_cells[cell_i].metadata[metadata_key] = region_id
                progress_bar.update()
        finally:
            pool.close()
            pool.join()
        progress_bar.close()
    else:
        for cell in tqdm.tqdm(output_cells, unit="cells"):
            pos = _scale_point((cell.z, cell.y, cell.x), data_size_voxels, regions_voxels)
            cell.metadata[metadata_key] = region_ids[pos].item()

    if output_cells_path:
        save_cells(output_cells, str(output_cells_path))
    logging.info(f"cell_meta_3d: Analysis took {datetime.now() - ts}")

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
        log_header="Cell3DScripts:LocateCells Log",
        multiprocessing_aware=True,
    )

    if args.signal_data_size:
        data_size_voxels = args.signal_data_size
    elif args.signal_planes_path:
        data_size_voxels = read_z_stack(args.signal_planes_path).shape
    else:
        raise ValueError("Either --signal-planes-path or --data-size must be provided")

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
            data_size_voxels=data_size_voxels,
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
