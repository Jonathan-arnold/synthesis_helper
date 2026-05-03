"""Download cytoscape.min.js into synthesis_helper/mcp/assets/.

The interactive HTML viewer embeds Cytoscape.js inline so that the rendered
page works offline. This script pulls a pinned version from jsDelivr and
verifies the SHA-256 digest.

Usage:
    python scripts/fetch_cytoscape.py
    python scripts/fetch_cytoscape.py --version 3.28.1

License: Cytoscape.js is MIT-licensed
(https://github.com/cytoscape/cytoscape.js/blob/master/LICENSE).
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

PINNED_VERSION = "3.28.1"
PINNED_SHA256 = "92d752b48ea949720675865197fd2a0001c95bc5888545e990af60321712d4c6"
URL_TEMPLATE = "https://cdn.jsdelivr.net/npm/cytoscape@{version}/dist/cytoscape.min.js"


def fetch(version: str, expected_sha256: str | None) -> bytes:
    url = URL_TEMPLATE.format(version=version)
    print(f"downloading {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = resp.read()
    digest = hashlib.sha256(data).hexdigest()
    print(f"sha256 = {digest}", file=sys.stderr)
    if expected_sha256 and digest != expected_sha256:
        print(
            f"sha256 mismatch: expected {expected_sha256}, got {digest}",
            file=sys.stderr,
        )
        sys.exit(1)
    return data


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default=PINNED_VERSION)
    ap.add_argument(
        "--skip-verify",
        action="store_true",
        help="do not enforce the pinned SHA-256 (useful when bumping versions)",
    )
    args = ap.parse_args()

    expected = None if (args.skip_verify or args.version != PINNED_VERSION) else PINNED_SHA256
    data = fetch(args.version, expected)

    out_dir = Path(__file__).resolve().parents[1] / "synthesis_helper" / "mcp" / "assets"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "cytoscape.min.js"
    out_path.write_bytes(data)
    print(f"wrote {len(data):,} bytes -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
