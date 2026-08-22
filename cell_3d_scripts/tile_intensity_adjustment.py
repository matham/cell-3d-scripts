import logging
from argparse import (
    ArgumentDefaultsHelpFormatter,
    ArgumentParser,
)
from datetime import datetime
from pathlib import Path

import fancylog
import numpy as np
import tifffile
from brainglobe_utils.cells.cells import Cell
from brainglobe_utils.IO.cells import get_cells, save_cells
from brainglobe_utils.IO.image import read_z_stack

import cell_3d_scripts
from cell_3d_scripts import __version__
from cell_3d_scripts.utils import get_tiled_image_intensity_factor


def arg_parser() -> ArgumentParser:
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("-ic", "--input-cells-path", dest="input_cells_path", type=Path, required=False, default=None)
    parser.add_argument("-oc", "--output-cells-path", dest="output_cells_path", type=Path, required=False, default=None)
    parser.add_argument(
        "-i",
        "--image-path",
        dest="image_path",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-p",
        "--projection-path",
        dest="projection_path",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-a",
        "--adjustment-path",
        dest="adjustment_path",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-n",
        "--num-tiles",
        dest="num_tiles",
        nargs=2,
        type=int,
        required=True,
        help="Axes is y, x, starting from a 3d z, y, x volume.",
    )
    parser.add_argument(
        "--sampling-step",
        dest="sampling_step",
        type=int,
        required=False,
        default=1,
    )
    parser.add_argument(
        "--n-patch-pixels",
        dest="n_patch_pixels",
        type=int,
        required=False,
        default=200,
    )
    parser.add_argument(
        "--skip-processed-cells",
        dest="skip_processed_cells",
        action="store_true",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def main(
    *,
    num_tiles: tuple[int, int],
    image_path: Path | None = None,
    cells: list[Cell] | None = None,
    projection_path: Path | None = None,
    adjustment_path: Path | None = None,
    sampling_step: int = 1,
    n_patch_pixels: int = 200,
) -> tuple[np.ndarray, np.ndarray, list[Cell] | None]:
    ts = datetime.now()
    logging.info(f"cell_3d_scripts.tile_intensity_adjustment: Starting tile intensity estimation for {image_path}")
    if image_path:
        projection, adjustment = get_tiled_image_intensity_factor(
            read_z_stack(str(image_path)),
            n_patch_pixels=n_patch_pixels,
            num_tiles=num_tiles,
            sampling_step=sampling_step,
        )

        if projection_path:
            tifffile.imwrite(projection_path, projection)
        if adjustment_path:
            tifffile.imwrite(adjustment_path, adjustment)
    elif projection_path and projection_path.exists() and adjustment_path and adjustment_path.exists():
        projection = tifffile.imread(projection_path)
        adjustment = tifffile.imread(adjustment_path)
    else:
        raise ValueError("If image path is not provided, the projection and adjustment images must be provided")

    if cells is not None:
        logging.info(f"Processing intensity adjustment for {len(cells)} cells")
        for cell in cells:
            cell.metadata["intensity_factor"] = adjustment[int(cell.y), int(cell.x)]

    logging.info(f"cell_3d_scripts.tile_intensity_adjustment: Analysis took {datetime.now() - ts}")

    return projection, adjustment, cells


def run_main():
    args = arg_parser().parse_args()

    image_path: Path | None = args.image_path
    input_cells_path: Path | None = args.input_cells_path
    output_cells_path: Path | None = args.output_cells_path
    projection_path: Path | None = args.projection_path
    adjustment_path: Path | None = args.adjustment_path
    skip_processed_cells = args.skip_processed_cells

    cells = None
    if input_cells_path:
        cells = get_cells(input_cells_path, cells_only=False)
        if skip_processed_cells and cells and "intensity_factor" in cells[0].metadata:
            print("Already processed cells")
            return
    elif output_cells_path:
        raise ValueError("output cells path provided, but no input cells")

    if projection_path:
        output_root = projection_path.parent
    elif adjustment_path:
        output_root = adjustment_path.parent
    elif output_cells_path:
        output_root = output_cells_path.parent
    elif image_path:
        output_root = image_path.parent
    else:
        output_root = Path.home()

    fancylog.start_logging(
        output_root,
        cell_3d_scripts,
        variables=[
            args,
        ],
        log_header="Cell3DScripts::TileIntensityAdjustment Log",
        multiprocessing_aware=True,
        filename="cell_3d_scripts-tile_intensity_adjustment",
        timestamp=True,
    )

    _, _, cells = main(
        image_path=image_path,
        num_tiles=args.num_tiles,
        cells=cells,
        projection_path=projection_path,
        adjustment_path=adjustment_path,
        sampling_step=args.sampling_step,
        n_patch_pixels=args.n_patch_pixels,
    )

    if output_cells_path:
        save_cells(cells, output_cells_path)


if __name__ == "__main__":
    run_main()
