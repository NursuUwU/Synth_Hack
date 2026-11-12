#!/bin/sh
set -eu
cd "$(dirname "$0")"
if [ -n "${GENOMIC_PYTHON:-}" ]; then
    if ! "$GENOMIC_PYTHON" -c 'import sys, sqlite3, ssl; sys.exit(0 if sys.version_info >= (3,11) else 1)'; then
        echo 'GENOMIC_PYTHON must point to Python 3.11+.' >&2
        exit 1
    fi
else
    for candidate in python3 python; do
        if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c 'import sys, sqlite3, ssl; sys.exit(0 if sys.version_info >= (3,11) else 1)' 2>/dev/null; then
            GENOMIC_PYTHON="$candidate"
            break
        fi
    done
fi
if [ -z "${GENOMIC_PYTHON:-}" ]; then
    echo 'Install Python 3.11+ from https://www.python.org/downloads/ and run sh start.sh again.' >&2
    exit 1
fi
export PYTHONUTF8=1
exec "$GENOMIC_PYTHON" -m genomic_detective serve --open "$@"
