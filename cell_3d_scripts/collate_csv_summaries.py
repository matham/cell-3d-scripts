import csv
import logging
from argparse import (
    ArgumentDefaultsHelpFormatter,
    ArgumentParser,
)
from datetime import datetime
from pathlib import Path

import fancylog

import cell_3d_scripts
from cell_3d_scripts import __version__


def arg_parser() -> ArgumentParser:
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-m",
        "--metadata-csv-path",
        dest="metadata_csv_path",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-r",
        "--summaries-root-path",
        dest="summaries_root_path",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-c",
        "--channel",
        dest="channels",
        type=int,
        required=True,
        action="append",
        default=[],
    )
    parser.add_argument(
        "-k",
        "--lookup-column-key",
        dest="lookup_column_key",
        required=True,
    )
    parser.add_argument(
        "-cols",
        "--first-n-cols",
        dest="first_n_cols",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-g",
        "--glob-search-format",
        dest="glob_search_format",
        required=True,
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
    metadata_csv_path: Path,
    summaries_root_path: Path,
    channels: list[int],
    lookup_column_key: str,
    first_n_cols: int,
    glob_search_format: str,
    output_path: Path,
) -> None:
    ts = datetime.now()
    logging.info("cell_3d_scripts.collate_csv_summaries: Starting CSV summary collation")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(metadata_csv_path) as fh:
        reader = csv.reader(fh, delimiter=",")
        dataset_lines = [[c.strip() for c in row] for row in reader]

    lookup_col = dataset_lines[0].index(lookup_column_key)
    dataset_header = dataset_lines[0][:first_n_cols]

    summary_data = []
    summary_header = None
    for dataset_line in dataset_lines[1:]:
        row = dataset_line[:first_n_cols]

        name = dataset_line[lookup_col]
        if not name:
            continue

        for channel in channels:
            glob_pat = glob_search_format.format(name=name, channel=channel)
            files = list(summaries_root_path.glob(glob_pat))

            if not files:
                logging.info(f"cell_3d_scripts.collate_csv_summaries: Didn't find matching file for {glob_pat}")
                continue
            if len(files) > 1:
                raise ValueError(f"Matched multiple files for pattern '{glob_pat}'. Found {files}")
            summary_filename = files[0]
            logging.info(f"cell_3d_scripts.collate_csv_summaries: Reading {summary_filename}")

            with open(summary_filename) as fh:
                summary_lines = list(csv.reader(fh, delimiter=","))

            if summary_header:
                if summary_lines[0] != summary_header:
                    raise ValueError(
                        f"Got a different header from {summary_header} to {summary_lines[0]} for {summary_filename}"
                    )
            summary_header = summary_lines[0]

            summary_data.append((row, channel, summary_filename.name, summary_lines[1:]))

    logging.info(f"cell_3d_scripts.collate_csv_summaries: Read {len(summary_data)} files")
    if summary_data:
        header = [*dataset_header, "Channel", "Summary filename", *summary_header]
        lines = [header]

        for dataset_row, channel, name, summary_lines in summary_data:
            for summary_line in summary_lines:
                line = [*dataset_row, str(channel), name, *summary_line]
                lines.append(line)

        with open(output_path, "w", newline="\n") as fh:
            writer = csv.writer(fh, delimiter=",")
            writer.writerows(lines)

    logging.info(f"cell_3d_scripts.collate_csv_summaries: Analysis took {datetime.now() - ts}")


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
        log_header="Cell3DScripts::CollateCSVSummaries Log",
        multiprocessing_aware=False,
    )

    main(
        metadata_csv_path=args.metadata_csv_path,
        summaries_root_path=args.summaries_root_path,
        channels=args.channels,
        lookup_column_key=args.lookup_column_key,
        first_n_cols=args.first_n_cols,
        glob_search_format=args.glob_search_format,
        output_path=args.output_path,
    )


if __name__ == "__main__":
    run_main()
