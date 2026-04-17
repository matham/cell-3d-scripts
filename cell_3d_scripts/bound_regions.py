import logging
import multiprocessing as mp
import tempfile
import zipfile
from argparse import (
    ArgumentDefaultsHelpFormatter,
    ArgumentParser,
)
from datetime import datetime
from functools import partial
from pathlib import Path

import fancylog
import numba.typed
import numpy as np
import tqdm
import yaml
import zarr
from brainglobe_utils.IO.image import read_z_stack
from numba import njit
from numba.core import types

import cell_3d_scripts
from cell_3d_scripts import __version__

numba_logger = logging.getLogger("numba")
numba_logger.setLevel(logging.WARNING)

_region_cache: tuple[int, np.ndarray] | None = None


def arg_parser() -> ArgumentParser:
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-r",
        "--region-id-path",
        dest="region_id_path",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-b",
        "--boundaries-path",
        dest="boundaries_path",
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


list_of_dims_bound = types.ListType(types.int64)
region_id_type = types.int64


@njit()
def _bound_plane(plane: np.ndarray) -> dict[int, list[int]]:
    boundaries = numba.typed.Dict.empty(key_type=region_id_type, value_type=list_of_dims_bound)

    for i in range(plane.shape[0]):
        for j in range(plane.shape[1]):
            region = plane[i, j]
            if region not in boundaries:
                region_boundaries = numba.typed.List.empty_list(types.int64, 4)
                region_boundaries.extend([i, i + 1, j, j + 1])
                boundaries[region] = region_boundaries
            else:
                region_boundaries = boundaries[region]

                region_boundaries[0] = min(i, region_boundaries[0])
                region_boundaries[1] = max(i + 1, region_boundaries[1])

                region_boundaries[2] = min(j, region_boundaries[2])
                region_boundaries[3] = max(j + 1, region_boundaries[3])

    return boundaries


def bound_plane(plane_i: int, ids: np.ndarray) -> tuple[int, dict[int, tuple[list[int], list[int]]]]:
    boundaries = _bound_plane(np.asarray(ids[plane_i, :, :]))
    boundaries = {rid: ([a, b], [c, d]) for rid, (a, b, c, d) in boundaries.items()}
    return plane_i, boundaries


def merge_plane_into_volume_boundaries(
    global_volume_boundaries: dict[int, tuple[list[int], list[int], list[int]]],
    plane_boundaries: dict[int, tuple[list[int], list[int]]],
    plane_i: int,
) -> None:
    for region, (plane_d1_boundaries, plane_d2_boundaries) in plane_boundaries.items():
        if region in global_volume_boundaries:
            vol_d1_boundaries, vol_d2_boundaries, vol_d3_boundaries = global_volume_boundaries[region]

            vol_d1_boundaries[0] = min(plane_i, vol_d1_boundaries[0])
            vol_d1_boundaries[1] = max(plane_i + 1, vol_d1_boundaries[1])

            vol_d2_boundaries[0] = min(plane_d1_boundaries[0], vol_d2_boundaries[0])
            vol_d2_boundaries[1] = max(plane_d1_boundaries[1], vol_d2_boundaries[1])

            vol_d3_boundaries[0] = min(plane_d2_boundaries[0], vol_d3_boundaries[0])
            vol_d3_boundaries[1] = max(plane_d2_boundaries[1], vol_d3_boundaries[1])
        else:
            global_volume_boundaries[region] = [plane_i, plane_i + 1], plane_d1_boundaries, plane_d2_boundaries


def main(
    *,
    region_ids: np.ndarray,
    boundaries_path: Path | None = None,
    num_workers: int = 6,
) -> dict[int, list[list[int]]]:
    ts = datetime.now()
    regions_voxels = region_ids.shape
    logging.info("cell_3d_scripts.bound_regions: Starting region boundary detection")
    logging.debug(f"Tiff containing the region IDs' shape is {regions_voxels}")

    n_planes = region_ids.shape[0]
    boundaries = {}

    if num_workers:
        progress_bar = tqdm.tqdm(total=n_planes, unit="planes")
        f = partial(bound_plane, ids=region_ids)
        ctx = mp.get_context("spawn")
        pool = ctx.Pool(processes=num_workers)
        try:
            for plane_i, plane_boundaries in pool.imap_unordered(f, list(range(n_planes))):
                merge_plane_into_volume_boundaries(boundaries, plane_boundaries, plane_i)
                progress_bar.update()
        finally:
            pool.close()
            pool.join()
        progress_bar.close()
    else:
        for plane_i in tqdm.tqdm(range(n_planes), unit="planes", total=n_planes):
            _, plane_boundaries = bound_plane(plane_i, region_ids)
            merge_plane_into_volume_boundaries(boundaries, plane_boundaries, plane_i)

    boundaries = {
        rid: [[d11, d12], [d21, d22], [d31, d32]] for rid, ((d11, d12), (d21, d22), (d31, d32)) in boundaries.items()
    }

    if boundaries_path is not None:
        with open(boundaries_path, "w") as fh:
            yaml.dump(boundaries, fh)
    logging.info(f"cell_3d_scripts.bound_regions: Analysis took {datetime.now() - ts}")

    return boundaries


def run_main():
    args = arg_parser().parse_args()

    region_id_path = Path(args.region_id_path)
    boundaries_path = Path(args.boundaries_path)
    output_root = boundaries_path.parent
    output_root.mkdir(parents=True, exist_ok=True)

    fancylog.start_logging(
        output_root,
        cell_3d_scripts,
        variables=[
            args,
        ],
        log_header="Cell3DScripts::BoundRegions Log",
        multiprocessing_aware=True,
        filename="cell_3d_scripts-bound_regions",
        timestamp=True,
    )

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
            region_ids=arr,
            boundaries_path=boundaries_path,
            num_workers=args.num_workers,
        )
    finally:
        if tempdir is not None:
            tempdir.cleanup()


if __name__ == "__main__":
    run_main()
