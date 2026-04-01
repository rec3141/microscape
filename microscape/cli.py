"""Command-line interface for microscape."""

import argparse
import sys

from . import __version__


def _cmd_auto_trim(args):
    """Analyze quality profiles and recommend truncation lengths."""
    from .quality import auto_trim

    result = auto_trim(
        args.input_dir,
        min_quality=args.min_quality,
        window=args.window,
        n_reads=args.n_reads,
        n_files=args.n_files,
        verbose=args.verbose,
    )

    # Output as key=value for easy parsing
    print(f"trunc_len_fwd={result['trunc_len_fwd']}")
    print(f"trunc_len_rev={result['trunc_len_rev']}")
    print(f"fwd_read_len={result['fwd_read_len']}")
    print(f"rev_read_len={result['rev_read_len']}")
    print(f"n_reads_sampled={result['n_reads_sampled']}")

    if args.output:
        with open(args.output, "w") as f:
            for key, val in result.items():
                f.write(f"{key}\t{val}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="microscape",
        description="Downstream analysis tools for amplicon sequencing",
    )
    parser.add_argument("--version", action="version",
                        version=f"microscape {__version__}")
    sub = parser.add_subparsers(dest="command")

    aq = sub.add_parser("auto-trim",
                        help="Recommend truncation lengths from quality profiles")
    aq.add_argument("input_dir", help="Directory containing paired FASTQ files")
    aq.add_argument("--min-quality", type=float, default=25,
                    help="Min rolling median Q to keep [default: 25]")
    aq.add_argument("--window", type=int, default=10,
                    help="Rolling window size [default: 10]")
    aq.add_argument("--n-reads", type=int, default=10000,
                    help="Reads to sample [default: 10000]")
    aq.add_argument("--n-files", type=int, default=20,
                    help="Files to sample from [default: 20]")
    aq.add_argument("-o", "--output", help="Write results to TSV file")
    aq.add_argument("-v", "--verbose", action="store_true")

    # -- prep-reads --
    pr = sub.add_parser("prep-reads",
                        help="Symlink FASTQs with run/plate naming for pipeline input")
    pr.add_argument("metadata", help="Metadata CSV/TSV file")
    pr.add_argument("reads_dir", help="Directory containing input FASTQ files")
    pr.add_argument("output_dir", help="Output directory for symlinks")
    pr.add_argument("--sample-col", default="Run",
                    help="Column matching sample IDs to filenames [default: Run]")
    pr.add_argument("--run-col", default="dada_run",
                    help="Column for sequencing run grouping [default: dada_run]")
    pr.add_argument("--plate-col", default="dada_plate",
                    help="Column for plate grouping [default: dada_plate]")
    pr.add_argument("--default-run", default="run1",
                    help="Default run name if column missing [default: run1]")
    pr.add_argument("--default-plate", default="plate1",
                    help="Default plate name if column missing [default: plate1]")
    pr.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "auto-trim":
        _cmd_auto_trim(args)
    elif args.command == "prep-reads":
        from .prep import prep_reads
        prep_reads(
            args.metadata, args.reads_dir, args.output_dir,
            sample_col=args.sample_col,
            run_col=args.run_col,
            plate_col=args.plate_col,
            default_run=args.default_run,
            default_plate=args.default_plate,
            verbose=args.verbose,
        )
