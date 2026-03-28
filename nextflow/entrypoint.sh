#!/bin/sh
# Microscape container entrypoint
#
# Validates the output directory is writable (important when running with
# docker --user) then launches Nextflow with the container config that
# disables conda (tools are already on PATH via wrapper scripts).

set -e

# Parse --outdir from arguments
OUTDIR=""
for arg in "$@"; do
    case "$prev" in
        --outdir) OUTDIR="$arg" ;;
    esac
    prev="$arg"
done

# Default output directory
if [ -z "$OUTDIR" ]; then
    OUTDIR="results"
fi

# Validate output directory is writable
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR" 2>/dev/null || true
fi

if [ -d "$OUTDIR" ] && [ ! -w "$OUTDIR" ]; then
    echo "[ERROR] Output directory is not writable: $OUTDIR" >&2
    echo "  If using --user, ensure the output directory exists and is writable." >&2
    exit 1
fi

exec nextflow -c /pipeline/docker.config "$@"
