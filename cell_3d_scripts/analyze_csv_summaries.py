import csv
import logging
from argparse import (
    ArgumentDefaultsHelpFormatter,
    ArgumentParser,
)
from collections import defaultdict
from copy import deepcopy
from datetime import datetime
from itertools import combinations
from pathlib import Path

import fancylog
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import stats

import cell_3d_scripts
from cell_3d_scripts import __version__
from cell_3d_scripts.atlas import AtlasTree

SUMMARY_DATA_TYPE = dict[tuple[str, ...], list[tuple[str, list[list[str]]]]]


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
        "-sk",
        "--lookup-sample-column-key",
        dest="lookup_sample_column_key",
        required=True,
    )
    parser.add_argument(
        "-gk",
        "--group-metadata-column-key",
        dest="groups_metadata_column_keys",
        type=str,
        required=True,
        action="append",
        default=[],
    )
    parser.add_argument(
        "-eg",
        "--exclude-group",
        dest="excluded_groups",
        type=str,
        nargs="+",
        required=False,
        action="append",
        default=[],
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
        type=str,
        required=True,
        action="append",
        default=[],
    )
    parser.add_argument(
        "-g",
        "--glob-search-format",
        dest="glob_search_format",
        required=True,
    )
    parser.add_argument(
        "--vaa3d-atlas-path",
        dest="vaa3d_atlas_path",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-rk",
        "--region_id-column-key",
        dest="region_id_column_key",
        required=False,
        default="ID",
    )
    parser.add_argument(
        "-ir",
        "--include-region",
        dest="include_regions",
        type=str,
        required=False,
        action="append",
        default=[],
    )
    parser.add_argument(
        "-er",
        "--exclude-region",
        dest="exclude_regions",
        type=str,
        required=False,
        action="append",
        default=[],
    )
    parser.add_argument(
        "-dk",
        "--stat-data-column-key",
        dest="stat_data_column_key",
        required=True,
    )
    parser.add_argument(
        "-ak",
        "--average-data-column-key",
        dest="average_data_columns_keys",
        type=str,
        required=True,
        action="append",
        default=[],
    )
    parser.add_argument(
        "-og",
        "--output-path-group-mean",
        dest="output_path_group_means",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-or",
        "--output-path-raw",
        dest="output_path_raw",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-oa",
        "--output-path-analyzed",
        dest="output_path_analyzed",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-oas",
        "--output-path-all-stat",
        dest="output_path_all_stat",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-oam",
        "--output-path-all-metadata",
        dest="output_path_all_metadata",
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        "--t-test",
        dest="t_test",
        action="store_true",
    )
    parser.add_argument(
        "--t-test-false-discovery",
        dest="t_test_false_discovery",
        action="store_true",
    )
    parser.add_argument(
        "--deseq2",
        dest="deseq2",
        action="store_true",
    )
    parser.add_argument(
        "--only-leafs",
        dest="only_leafs",
        action="store_true",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def _read_data(
    sample_name_col: int,
    metadata_lines: list[list[str]],
    metadata_header: list[str],
    group_cols: dict[str, int],
    excluded_groups: list[dict[str, str]],
    channels: list[str],
    glob_search_format: str,
    summaries_root_path: Path,
    group_names: list[str],
) -> tuple[SUMMARY_DATA_TYPE, list[str], int]:
    summary_data: SUMMARY_DATA_TYPE = defaultdict(list)
    summary_header = None
    count = 0

    for metadata_line in metadata_lines[1:]:
        name = metadata_line[sample_name_col]
        if not name:
            continue

        matched_excluded = False
        for excluded_group in excluded_groups:
            our_group = {k: metadata_line[metadata_header.index(k)].lower() for k, v in excluded_group.items()}
            matched_excluded = our_group == excluded_group
            if matched_excluded:
                break
        if matched_excluded:
            continue

        groups = {k: metadata_line[i].lower() for k, i in group_cols.items()}
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
                summary_lines: list[list[str]] = list(csv.reader(fh, delimiter=","))

            if summary_header:
                if summary_lines[0] != summary_header:
                    raise ValueError(
                        f"Got a different header from {summary_header} to {summary_lines[0]} for {summary_filename}"
                    )
            summary_header = summary_lines[0]

            key = [groups[g] for g in group_names]
            summary_data[(*key, channel)].append((name, summary_lines[1:]))
            count += 1

    return summary_data, summary_header, count


def _remove_excluded_regions(
    vaa3d_atlas_path: Path,
    include_regions: list[str] | None,
    exclude_regions: list[str] | None,
    region_id_column_key: str,
    summary_header: list[str],
    summary_data: SUMMARY_DATA_TYPE,
    only_leafs: bool,
) -> tuple[list[bool], list[int]]:
    tree = AtlasTree.parse_vaa3d(vaa3d_atlas_path)
    exclude_regions = {int(r) for r in (exclude_regions or [])}
    include_regions = {int(r) for r in (include_regions or [])}
    regions_col = summary_header.index(region_id_column_key)
    line_leafiness = []
    lines_regions = []

    if not summary_data:
        raise ValueError("No data read")
    sample = list(summary_data.values())[0][0][1]
    for line in sample:
        region = int(line[regions_col])
        line_leafiness.append(not only_leafs or region not in tree.nodes_map or not tree.nodes_map[region].children)
        lines_regions.append(region)

    exclusions = set()
    for region in exclude_regions:
        node = tree.nodes_map[region]
        for elem in node.pre_order():
            exclusions.add(elem.region_id)

    for region in include_regions:
        node = tree.nodes_map[region]
        if not node.parent:
            continue

        for child in node.parent.children:
            if child.region_id in include_regions:
                continue
            for elem in child.pre_order():
                exclusions.add(elem.region_id)

    summary_len = None
    cleaned_line_leafiness = False
    if exclusions:
        for key, group_data in summary_data.items():
            for name, lines in group_data:
                excluded_lines = []
                for i, line in enumerate(lines):
                    if int(line[regions_col]) in exclusions:
                        excluded_lines.append(i)

                for i in reversed(excluded_lines):
                    del lines[i]
                    if not cleaned_line_leafiness:
                        del line_leafiness[i]
                        del lines_regions[i]
                cleaned_line_leafiness = True

                if summary_len is not None and summary_len != len(lines):
                    raise ValueError("Found different file lengths after filtering for excluded regions")
                summary_len = len(lines)

    return line_leafiness, lines_regions


def _flatten_data(
    summary_data: SUMMARY_DATA_TYPE,
    summary_header: list[str],
    channel: str,
    stat_data_column_key: str,
    average_data_columns_keys: str,
) -> tuple[list[list[str]], list[list[str]], list[list[str]], dict[tuple[str, ...], list[tuple[str, np.ndarray]]]]:
    summary_data = {k: summary_data[k] for k in sorted(summary_data.keys())}

    summary_header_l = [s.lower() for s in summary_header]
    stat_col = summary_header_l.index(stat_data_column_key.lower())
    average_cols = [summary_header_l.index(k.lower()) for k in average_data_columns_keys]
    data_cols = list(set(average_cols + [stat_col]))
    data_cols_names = [summary_header[c] for c in data_cols]

    output_header = [value for i, value in enumerate(summary_header) if i not in data_cols]
    sample_lines = list(summary_data.values())[0][0][1]
    output_lines = [[value for i, value in enumerate(line) if i not in data_cols] for line in sample_lines]

    lines_bare = deepcopy(output_lines)
    lines_bare.insert(0, deepcopy(output_header))

    header_all = deepcopy(output_header)
    lines_all = deepcopy(output_lines)
    for col_name, col in zip(data_cols_names, data_cols, strict=False):
        for key, group_data in summary_data.items():
            if key[-1] != channel:
                continue

            for name, lines in group_data:
                label = f"{name}:{col_name}"
                header_all.append(label)

                for line_all, line in zip(lines_all, lines, strict=False):
                    line_all.append(line[col])

    lines_all.insert(0, header_all)

    header_grouped = deepcopy(output_header)
    lines_grouped = deepcopy(output_lines)
    for key, group_data in summary_data.items():
        if key[-1] != channel:
            continue

        label = "-".join(key[:-1]) + ":N"
        header_grouped.append(label)
        for line in lines_grouped:
            line.append(str(len(group_data)))

    for col_name, col in zip(data_cols_names, data_cols, strict=False):
        for key, group_data in summary_data.items():
            if key[-1] != channel:
                continue

            data = [[] for _ in range(len(group_data[0][1]))]
            for _, lines in group_data:
                for line_data, line in zip(data, lines, strict=False):
                    line_data.append(float(line[col]))

            mean_data = np.mean(data, axis=1).tolist()
            label = "-".join(key[:-1]) + f":{col_name}"
            header_grouped.append(label)
            for line_grouped, item in zip(lines_grouped, mean_data, strict=False):
                line_grouped.append(str(item))
    lines_grouped.insert(0, header_grouped)

    stat_data = defaultdict(list)
    for key, group_data in summary_data.items():
        if key[-1] != channel:
            continue

        for name, lines in group_data:
            data = np.array([float(line[stat_col]) for line in lines])
            stat_data[key[:-1]].append((name, data))

    return lines_bare, lines_all, lines_grouped, stat_data


def _flatten_stats_data(
    stat_data: dict[tuple[str, ...], list[tuple[str, np.ndarray]]],
    line_leafiness: list[bool],
    lines_regions: list[int],
    group_names: list[str],
    convert_data_values_to_str: bool = True,
    leafs_only: bool = True,
) -> tuple[list[list[str]], list[list[str]]]:
    stat_lines = [
        [
            "RegionID",
        ]
    ]
    for val, leaf in zip(lines_regions, line_leafiness, strict=False):
        if leaf or not leafs_only:
            stat_lines.append([str(val)])

    metadata_lines = [["Sample", *group_names]]

    for key, items in stat_data.items():
        for name, arr in items:
            metadata_lines.append([name, *key])

            stat_lines[0].append(name)
            i = 1
            for val, leaf in zip(arr.tolist(), line_leafiness, strict=False):
                if leaf or not leafs_only:
                    stat_lines[i].append(str(val) if convert_data_values_to_str else val)
                    i += 1

    return stat_lines, metadata_lines


def _generate_stats(
    data: dict[tuple[str, ...], list[tuple[str, np.ndarray]]],
    lines: list[list[str]],
    t_test: bool,
    t_test_false_discovery: bool,
    deseq2: bool,
    line_leafiness: list[bool],
    lines_regions: list[int],
    group_names: list[str],
) -> None:
    line_leafiness = np.array(line_leafiness, dtype=bool)
    assert len(line_leafiness) == len(lines_regions)

    if t_test or t_test_false_discovery:
        for group1, group2 in combinations(data.keys(), 2):
            g1 = np.column_stack([d[1] for d in data[group1]])
            g2 = np.column_stack([d[1] for d in data[group2]])
            p = stats.ttest_ind(g1, g2, axis=1, equal_var=False).pvalue

            group1_label = "-".join(group1)
            group2_label = "-".join(group2)
            label = f"{group1_label} vs. {group2_label}"

            if t_test or t_test_false_discovery:
                lines[0].append(f"p: {label}")
                assert len(p) + 1 == len(lines)
                for i, value in enumerate(p):
                    lines[i + 1].append(str(value))

            if t_test_false_discovery:
                assert len(p) == len(line_leafiness)
                valid = np.logical_and(np.logical_not(np.isnan(p)), line_leafiness)
                p_corr = np.copy(p)
                p_corr[valid] = stats.false_discovery_control(p[valid])
                p_corr[np.logical_not(valid)] = -1

                lines[0].append(f"p-corr: {label}")
                assert len(p_corr) + 1 == len(lines)
                for i, value in enumerate(p_corr):
                    lines[i + 1].append(str(value))

    if deseq2:
        stat_lines, metadata_lines = _flatten_stats_data(
            data,
            line_leafiness.tolist(),
            lines_regions,
            group_names,
            False,
            False,
        )

        sample_names = stat_lines[0][1:]
        assert len(line_leafiness) + 1 == len(stat_lines)
        region_names = [line[0] for line in stat_lines[1:]]
        data_np = np.array([line[1:] for line in stat_lines[1:]]).round()

        valid = np.logical_and(line_leafiness, np.sum(data_np, axis=1) > 0)
        region_names = [r for i, r in enumerate(region_names) if valid[i]]
        data_df = pd.DataFrame(data_np[valid, :].T, index=sample_names, columns=region_names)

        assert group_names == metadata_lines[0][1:]
        assert sample_names == [line[0] for line in metadata_lines[1:]]
        metadata_values = [line[1:] for line in metadata_lines[1:]]
        metadata_df = pd.DataFrame(metadata_values, index=sample_names, columns=group_names)

        design = "~" + " + ".join(group_names)
        dds = DeseqDataSet(counts=data_df, metadata=metadata_df, design=design, n_cpus=4, refit_cooks=True)
        dds.deseq2()

        for g, group_name in enumerate(group_names):
            levels = sorted({sample[g] for sample in metadata_values})
            if len(levels) < 2:
                continue

            for name1, name2 in combinations(levels, 2):
                label_base = f"{group_name}:{name1}-vs-{name2}"
                ds = DeseqStats(
                    dds,
                    contrast=[group_name, name1, name2],
                    n_cpus=4,
                    cooks_filter=True,
                    alpha=0.05,
                    independent_filter=True,
                )
                ds.summary()

                result = ds.results_df
                for measure in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
                    lines[0].append(f"{label_base}:{measure}")
                    vals = np.full(len(valid), -1, dtype=np.float64)
                    vals[valid] = result[measure].to_numpy()
                    for i, line in enumerate(lines[1:]):
                        line.append(str(vals[i].item()))


def save_data(
    channel: str,
    output_path_raw: Path | None,
    output_path_analyzed: Path | None,
    output_path_all_metadata: Path | None,
    output_path_all_stat: Path | None,
    output_path_group_means: Path | None,
    lines_all: list[list[str]],
    lines_analyzed: list[list[str]],
    lines_grouped: list[list[str]],
    stat_data: dict[tuple[str, ...], list[tuple[str, np.ndarray]]],
    line_leafiness: list[bool],
    lines_regions: list[int],
    group_names: list[str],
) -> None:
    if output_path_raw:
        output_path_raw.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path_raw.parent / output_path_raw.name.format(channel=channel), "w", newline="\n") as fh:
            writer = csv.writer(fh, delimiter=",")
            writer.writerows(lines_all)

    if output_path_group_means:
        output_path_group_means.parent.mkdir(parents=True, exist_ok=True)
        with open(
            output_path_group_means.parent / output_path_group_means.name.format(channel=channel),
            "w",
            newline="\n",
        ) as fh:
            writer = csv.writer(fh, delimiter=",")
            writer.writerows(lines_grouped)

    stat_lines, metadata_lines = _flatten_stats_data(stat_data, line_leafiness, lines_regions, group_names)

    if output_path_all_stat:
        output_path_all_stat.parent.mkdir(parents=True, exist_ok=True)
        with open(
            output_path_all_stat.parent / output_path_all_stat.name.format(channel=channel),
            "w",
            newline="\n",
        ) as fh:
            writer = csv.writer(fh, delimiter=",")
            writer.writerows(stat_lines)

    if output_path_all_metadata:
        output_path_all_metadata.parent.mkdir(parents=True, exist_ok=True)
        with open(
            output_path_all_metadata.parent / output_path_all_metadata.name.format(channel=channel),
            "w",
            newline="\n",
        ) as fh:
            writer = csv.writer(fh, delimiter=",")
            writer.writerows(metadata_lines)

    if output_path_analyzed:
        output_path_analyzed.parent.mkdir(parents=True, exist_ok=True)
        with open(
            output_path_analyzed.parent / output_path_analyzed.name.format(channel=channel),
            "w",
            newline="\n",
        ) as fh:
            writer = csv.writer(fh, delimiter=",")
            writer.writerows(lines_analyzed)


def main(
    *,
    metadata_csv_path: Path,
    lookup_sample_column_key: str,
    groups_metadata_column_keys: list[str],
    summaries_root_path: Path,
    channels: list[str],
    glob_search_format: str,
    stat_data_column_key: str,
    average_data_columns_keys: str,
    vaa3d_atlas_path: Path,
    output_path_raw: Path | None = None,
    output_path_analyzed: Path | None = None,
    output_path_all_metadata: Path | None = None,
    output_path_all_stat: Path | None = None,
    output_path_group_means: Path | None = None,
    excluded_groups: list[list[str]] | None = None,
    exclude_regions: list[str] | None = None,
    include_regions: list[str] | None = None,
    region_id_column_key: str = "Abbreviation",
    t_test: bool = False,
    t_test_false_discovery: bool = False,
    deseq2: bool = False,
    only_leafs: bool = False,
) -> None:
    ts = datetime.now()
    logging.info("cell_3d_scripts.analyze_csv_summaries: Starting CSV summary data analysis")

    with open(metadata_csv_path) as fh:
        reader = csv.reader(fh, delimiter=",")
        metadata_lines: list[list[str]] = [[c.strip() for c in row] for row in reader]

    groups_metadata_column_keys = [k.lower() for k in groups_metadata_column_keys]
    metadata_header = [v.lower() for v in metadata_lines[0]]

    sample_name_col: int = metadata_header.index(lookup_sample_column_key.lower())
    group_cols: dict[str, int] = {k: metadata_header.index(k) for k in groups_metadata_column_keys}
    excluded_groups: list[dict[str, str]] = [dict([g.lower().split("=") for g in groups]) for groups in excluded_groups]
    if "channel" in group_cols:
        raise ValueError('"Channel" may not be a column name in the metadata csv')
    group_names = list(group_cols.keys())

    summary_data, summary_header, count = _read_data(
        sample_name_col,
        metadata_lines,
        metadata_header,
        group_cols,
        excluded_groups,
        channels,
        glob_search_format,
        summaries_root_path,
        group_names,
    )

    logging.info(f"cell_3d_scripts.analyze_csv_summaries: Read {count} files")
    if not summary_data:
        logging.info(f"cell_3d_scripts.analyze_csv_summaries: No data. Analysis took {datetime.now() - ts}")
        return

    # filter out excluded regions from the read data
    line_leafiness, lines_regions = _remove_excluded_regions(
        vaa3d_atlas_path,
        include_regions,
        exclude_regions,
        region_id_column_key,
        summary_header,
        summary_data,
        only_leafs,
    )

    for channel in channels:
        lines_bare, lines_all, lines_grouped, stat_data = _flatten_data(
            summary_data, summary_header, channel, stat_data_column_key, average_data_columns_keys
        )

        lines_analyzed = deepcopy(lines_bare)
        _generate_stats(
            stat_data,
            lines_analyzed,
            t_test,
            t_test_false_discovery,
            deseq2,
            line_leafiness,
            lines_regions,
            group_names,
        )

        save_data(
            channel,
            output_path_raw,
            output_path_analyzed,
            output_path_all_metadata,
            output_path_all_stat,
            output_path_group_means,
            lines_all,
            lines_analyzed,
            lines_grouped,
            stat_data,
            line_leafiness,
            lines_regions,
            group_names,
        )

    logging.info(f"cell_3d_scripts.analyze_csv_summaries: Analysis took {datetime.now() - ts}")


def run_main():
    args = arg_parser().parse_args()

    output_root = args.metadata_csv_path.parent
    if args.output_path_raw is not None:
        output_root = args.output_path_raw.parent
    elif args.output_path_analyzed is not None:
        output_root = args.output_path_analyzed.parent
    elif args.output_path_all_metadata is not None:
        output_root = args.output_path_all_metadata.parent
    elif args.output_path_all_stat is not None:
        output_root = args.output_path_all_stat.parent
    elif args.output_path_group_means is not None:
        output_root = args.output_path_group_means.parent
    output_root.mkdir(parents=True, exist_ok=True)

    fancylog.start_logging(
        output_root,
        cell_3d_scripts,
        variables=[
            args,
        ],
        log_header="Cell3DScripts::CollateCSVSummaries Log",
        multiprocessing_aware=False,
        filename="cell_3d_scripts-analyze_csv_summaries",
        timestamp=True,
    )

    main(
        metadata_csv_path=args.metadata_csv_path,
        groups_metadata_column_keys=args.groups_metadata_column_keys,
        excluded_groups=args.excluded_groups,
        summaries_root_path=args.summaries_root_path,
        channels=args.channels,
        lookup_sample_column_key=args.lookup_sample_column_key,
        glob_search_format=args.glob_search_format,
        exclude_regions=args.exclude_regions,
        include_regions=args.include_regions,
        stat_data_column_key=args.stat_data_column_key,
        average_data_columns_keys=args.average_data_columns_keys,
        output_path_analyzed=args.output_path_analyzed,
        output_path_raw=args.output_path_raw,
        output_path_all_metadata=args.output_path_all_metadata,
        output_path_all_stat=args.output_path_all_stat,
        output_path_group_means=args.output_path_group_means,
        region_id_column_key=args.region_id_column_key,
        vaa3d_atlas_path=args.vaa3d_atlas_path,
        t_test=args.t_test,
        t_test_false_discovery=args.t_test_false_discovery,
        deseq2=args.deseq2,
        only_leafs=args.only_leafs,
    )


if __name__ == "__main__":
    run_main()
