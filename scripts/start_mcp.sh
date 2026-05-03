#!/bin/bash
# Launcher for the synthesis_helper MCP server.
#
# Why this script exists: on macOS, when the project lives under ~/Documents or
# ~/Desktop with iCloud "Optimize Mac Storage" enabled, Python files in the
# venv (and the project itself) get evicted to cloud stubs (dataless). On
# startup Python then materializes hundreds of .py files serially from iCloud,
# turning a 1-second import into 5+ minutes. We mitigate two ways:
#   1. Use a venv outside any synced folder (~/.venvs/synthesis_helper).
#   2. Warm the project source with brctl/cat before launch so the MCP
#      handshake doesn't stall while files download on demand.
set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_PY="${SYNTHESIS_HELPER_VENV:-$HOME/.venvs/synthesis_helper/bin/python}"

if [[ ! -x "$VENV_PY" ]]; then
  echo "synthesis_helper: venv python not found at $VENV_PY" >&2
  echo "Create it with: UV_PROJECT_ENVIRONMENT=\"\$HOME/.venvs/synthesis_helper\" uv sync --project \"$PROJECT_DIR\" --python 3.13" >&2
  exit 1
fi

# Warm any iCloud-evicted project files (no-op if already local or not on iCloud).
if command -v brctl >/dev/null 2>&1; then
  brctl download "$PROJECT_DIR/synthesis_helper" >/dev/null 2>&1 || true
  brctl download "$PROJECT_DIR/data" >/dev/null 2>&1 || true
fi

cd "$PROJECT_DIR"
exec "$VENV_PY" -m synthesis_helper "$@"
