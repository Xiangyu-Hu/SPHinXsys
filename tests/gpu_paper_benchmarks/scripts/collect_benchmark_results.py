#!/usr/bin/env python3
"""Validate and combine GPU-paper benchmark summary files."""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence


REQUIRED_COLUMNS = (
    "case",
    "run",
    "git",
    "build",
    "precision",
    "backend",
    "device",
    "dp",
    "initial_particle_count",
    "particle_count",
    "physical_time",
    "outer_steps",
    "advection_steps",
    "acoustic_steps",
    "solid_steps",
    "wall_seconds",
    "compute_seconds",
    "io_seconds",
    "init_seconds",
    "time_per_outer_step",
    "status",
    "requested_end_time",
    "benchmark_mode",
    "output_enabled",
    "output_interval",
    "resolution",
    "verification_seconds",
    "particle_updates",
    "particle_update_definition",
    "particle_updates_per_second",
    "mpps",
    "neighbor_interactions",
    "mnips",
    "gpips",
    "peak_rss_kb",
    "peak_gpu_memory_kb",
    "sorting_seconds",
    "cll_seconds",
    "configuration_seconds",
    "advection_seconds",
    "acoustic_component_seconds",
    "pressure_seconds",
    "density_seconds",
    "viscous_seconds",
    "solid_component_seconds",
    "contact_seconds",
    "fsi_seconds",
    "diffusion_seconds",
    "reduction_seconds",
)
BASE_GROUP_COLUMNS = (
    "case",
    "git",
    "build",
    "precision",
    "backend",
    "device",
    "benchmark_mode",
    "output_enabled",
    "output_interval",
    "resolution",
    "dp",
    "requested_end_time",
)
DEFAULT_STAT_COLUMNS = (
    "wall_seconds",
    "compute_seconds",
    "io_seconds",
    "init_seconds",
    "time_per_outer_step",
    "particle_updates_per_second",
    "mpps",
)


class CollectionError(RuntimeError):
    """A benchmark input is malformed or incompatible."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recursively collect summary.csv files, validating their required "
            "schema and deterministically unioning case-specific extra columns."
        )
    )
    parser.add_argument("result_root", type=Path, help="raw result tree to scan")
    parser.add_argument("output_csv", type=Path, help="combined CSV destination")
    parser.add_argument(
        "--stats-output",
        type=Path,
        help="optional repeated-run statistics CSV destination",
    )
    parser.add_argument(
        "--allow-failures",
        action="store_true",
        help=(
            "continue when failed run.status directories have no summary.csv; "
            "each skipped failure is still reported on stderr"
        ),
    )
    parser.add_argument(
        "--stat-column",
        action="append",
        dest="stat_columns",
        metavar="COLUMN",
        help=(
            "numeric column to summarize; may be repeated "
            "(default: standard timing columns)"
        ),
    )
    return parser.parse_args(argv)


def display_path(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def validate_header(path: Path, header: list[str]) -> None:
    if not header:
        raise CollectionError(f"{path}: empty CSV header")
    empty_positions = [str(index + 1) for index, name in enumerate(header) if not name]
    if empty_positions:
        raise CollectionError(
            f"{path}: empty header name at column(s) {', '.join(empty_positions)}"
        )
    duplicates = sorted({name for name in header if header.count(name) > 1})
    if duplicates:
        raise CollectionError(
            f"{path}: duplicate header column(s): {', '.join(duplicates)}"
        )
    missing = [name for name in REQUIRED_COLUMNS if name not in header]
    if missing:
        raise CollectionError(
            f"{path}: missing required column(s): {', '.join(missing)}"
        )
    required_prefix = header[: len(REQUIRED_COLUMNS)]
    if required_prefix != list(REQUIRED_COLUMNS):
        raise CollectionError(
            f"{path}: required columns must be the first columns in this exact "
            f"order: {', '.join(REQUIRED_COLUMNS)}"
        )
    if "source_summary" in header:
        raise CollectionError(
            f"{path}: input already contains reserved column 'source_summary'"
        )


def read_run_status(path: Path) -> dict[str, str]:
    fields: dict[str, str] = {}
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if "=" not in raw_line:
            raise CollectionError(
                f"{path}:{line_number}: expected a key=value status entry"
            )
        key, value = raw_line.split("=", 1)
        if not key or key in fields:
            raise CollectionError(
                f"{path}:{line_number}: empty or duplicate status key {key!r}"
            )
        fields[key] = value
    if not fields.get("status"):
        raise CollectionError(f"{path}: missing non-empty 'status' entry")
    return fields


def scan_run_statuses(
    root: Path, allow_failures: bool
) -> set[Path]:
    noncompleted_summaries: set[Path] = set()
    missing_failed_summaries: list[Path] = []

    for status_path in sorted(root.rglob("run.status")):
        fields = read_run_status(status_path)
        status = fields["status"]
        summary_path = status_path.parent / "summary.csv"
        relative_status = display_path(status_path, root)

        if summary_path.is_file():
            if status != "completed":
                noncompleted_summaries.add(summary_path.resolve())
                print(
                    f"WARNING: {relative_status}: status={status!r}; summary "
                    "will be included in raw output but excluded from statistics",
                    file=sys.stderr,
                )
        elif status == "failed":
            missing_failed_summaries.append(status_path)
            print(
                f"WARNING: {relative_status}: failed run has no summary.csv",
                file=sys.stderr,
            )
        else:
            raise CollectionError(
                f"{status_path}: status={status!r} but summary.csv is missing"
            )

    if missing_failed_summaries and not allow_failures:
        raise CollectionError(
            f"{len(missing_failed_summaries)} failed run(s) have no summary.csv; "
            "rerun with --allow-failures to collect the remaining summaries"
        )
    return noncompleted_summaries


def read_summaries(
    root: Path, output_paths: Iterable[Path]
) -> tuple[list[str], list[dict[str, str]], list[Path]]:
    excluded = {path.resolve() for path in output_paths}
    summaries = sorted(
        path for path in root.rglob("summary.csv") if path.resolve() not in excluded
    )
    if not summaries:
        raise CollectionError(f"{root}: no summary.csv files found")

    optional_columns: set[str] = set()
    combined_rows: list[dict[str, str]] = []
    identities: dict[tuple[str, str], tuple[Path, int]] = {}

    for path in summaries:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.reader(stream)
            try:
                header = next(reader)
            except StopIteration as error:
                raise CollectionError(f"{path}: empty CSV file") from error
            validate_header(path, header)
            optional_columns.update(header[len(REQUIRED_COLUMNS) :])

            row_count = 0
            for line_number, values in enumerate(reader, start=2):
                if len(values) != len(header):
                    raise CollectionError(
                        f"{path}:{line_number}: expected {len(header)} fields, "
                        f"found {len(values)}"
                    )
                row = dict(zip(header, values))
                identity = (row["case"], row["run"])
                if identity in identities:
                    previous_path, previous_line = identities[identity]
                    raise CollectionError(
                        f"{path}:{line_number}: duplicate (case, run)={identity!r}; "
                        f"first seen at {previous_path}:{previous_line}"
                    )
                identities[identity] = (path, line_number)
                row["source_summary"] = display_path(path, root)
                combined_rows.append(row)
                row_count += 1
            if row_count == 0:
                raise CollectionError(f"{path}: header is present but no data rows exist")

    union_header = [
        *REQUIRED_COLUMNS,
        *sorted(optional_columns),
        "source_summary",
    ]
    normalized_rows = [
        {name: row.get(name, "") for name in union_header} for row in combined_rows
    ]
    return union_header, normalized_rows, summaries


def atomic_write_csv(
    destination: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, object]]
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            delete=False,
        ) as stream:
            temporary_name = stream.name
            writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="raise")
            writer.writeheader()
            writer.writerows(rows)
        Path(temporary_name).replace(destination)
    except Exception:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
        raise


def statistic_rows(
    rows: Sequence[dict[str, str]],
    stat_columns: Sequence[str],
    group_columns: Sequence[str],
) -> list[dict[str, object]]:
    if not rows:
        return []
    missing_group = [name for name in group_columns if name not in rows[0]]
    missing_stats = [name for name in stat_columns if name not in rows[0]]
    if missing_group or missing_stats:
        pieces = []
        if missing_group:
            pieces.append(f"grouping columns: {', '.join(missing_group)}")
        if missing_stats:
            pieces.append(f"statistic columns: {', '.join(missing_stats)}")
        raise CollectionError("statistics require missing " + "; ".join(pieces))

    grouped: dict[tuple[str, ...], dict[str, list[float]]] = defaultdict(
        lambda: {name: [] for name in stat_columns}
    )
    for row in rows:
        key = tuple(row[name] for name in group_columns)
        for name in stat_columns:
            try:
                grouped[key][name].append(float(row[name]))
            except ValueError as error:
                raise CollectionError(
                    f"{row['source_summary']}: non-numeric value {row[name]!r} "
                    f"in statistic column {name!r}"
                ) from error

    output: list[dict[str, object]] = []
    for key in sorted(grouped):
        group_values = grouped[key]
        for name in stat_columns:
            values = group_values[name]
            output.append(
                {
                    **dict(zip(group_columns, key)),
                    "metric": name,
                    "count": len(values),
                    "mean": statistics.fmean(values),
                    "stddev": statistics.stdev(values) if len(values) > 1 else 0.0,
                    "min": min(values),
                    "median": statistics.median(values),
                }
            )
    return output


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.result_root.resolve()
    if not root.is_dir():
        raise CollectionError(f"{root}: result root is not a directory")

    output_paths = [args.output_csv]
    if args.stats_output is not None:
        output_paths.append(args.stats_output)
    noncompleted_summaries = scan_run_statuses(root, args.allow_failures)
    header, rows, summaries = read_summaries(root, output_paths)
    completed_rows: list[dict[str, str]] = []
    reported_csv_failures: set[str] = set()
    for row in rows:
        source_path = (root / row["source_summary"]).resolve()
        if row["status"] == "completed" and source_path not in noncompleted_summaries:
            completed_rows.append(row)
        elif (
            source_path not in noncompleted_summaries
            and row["source_summary"] not in reported_csv_failures
        ):
            print(
                f"WARNING: {row['source_summary']}: non-completed summary is "
                "included in raw output but excluded from statistics",
                file=sys.stderr,
            )
            reported_csv_failures.add(row["source_summary"])

    stats: list[dict[str, object]] | None = None
    if args.stats_output is not None:
        stat_columns = tuple(args.stat_columns or DEFAULT_STAT_COLUMNS)
        stats = statistic_rows(completed_rows, stat_columns, BASE_GROUP_COLUMNS)
        stats_header = (
            *BASE_GROUP_COLUMNS,
            "metric",
            "count",
            "mean",
            "stddev",
            "min",
            "median",
        )
    atomic_write_csv(args.output_csv, header, rows)
    if args.stats_output is not None:
        assert stats is not None
        atomic_write_csv(args.stats_output, stats_header, stats)

    print("Union columns: " + ", ".join(header[:-1]))
    print(
        f"Collected {len(rows)} row(s) from {len(summaries)} summary file(s) "
        f"into {args.output_csv}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CollectionError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
