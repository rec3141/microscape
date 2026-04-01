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


def _cmd_serve(args):
    """Serve the viz app pointing at pipeline results."""
    import os
    import shutil
    import subprocess

    data_dir = os.path.abspath(args.data_dir)
    if not os.path.isdir(data_dir):
        print(f"ERROR: {data_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    # Find viz app source
    viz_dir = args.viz_dir
    if viz_dir is None:
        # Try common locations
        candidates = [
            os.path.join(os.path.dirname(data_dir), "..", "viz"),  # results/../viz (pipeline dir)
            os.path.join(os.getcwd(), "viz"),
        ]
        # Also check if microscape-nf is a sibling or parent
        for parent in [os.getcwd(), os.path.dirname(os.getcwd())]:
            candidates.append(os.path.join(parent, "microscape-nf", "viz"))
        for c in candidates:
            if os.path.isfile(os.path.join(c, "package.json")):
                viz_dir = c
                break

    if viz_dir is None or not os.path.isfile(os.path.join(viz_dir, "package.json")):
        print("ERROR: Cannot find viz app. Pass --viz-dir /path/to/microscape-nf/viz",
              file=sys.stderr)
        sys.exit(1)

    viz_dir = os.path.abspath(viz_dir)

    # Install deps if needed
    if not os.path.isdir(os.path.join(viz_dir, "node_modules")):
        print(f"Installing viz dependencies in {viz_dir}...")
        subprocess.run(["npm", "install"], cwd=viz_dir, check=True)

    print(f"Serving viz at http://localhost:{args.port}")
    print(f"Data: {data_dir}")
    print(f"App:  {viz_dir}")
    print("Press Ctrl+C to stop")

    env = os.environ.copy()
    env["VIZ_DATA_DIR"] = data_dir
    subprocess.run(
        ["npx", "vite", "--host", "0.0.0.0", "--port", str(args.port)],
        cwd=viz_dir, env=env,
    )


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

    # -- serve --
    sv = sub.add_parser("serve",
                        help="Serve the viz app with pipeline results")
    sv.add_argument("data_dir", help="Directory containing viz JSON files (results/viz/)")
    sv.add_argument("--port", type=int, default=5174, help="Port [default: 5174]")
    sv.add_argument("--viz-dir", default=None,
                    help="Path to viz app source [default: auto-detect from microscape-nf]")

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
    elif args.command == "serve":
        _cmd_serve(args)
