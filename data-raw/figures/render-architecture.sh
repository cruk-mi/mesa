#!/usr/bin/env bash
#
# Regenerate man/figures/README-architecture.svg from architecture.mmd.
#
# The SVG is committed because README images must be present in the repository
# and in the source tarball; this script is how it gets rebuilt. data-raw/ is
# .Rbuildignore'd, so neither this script nor the Mermaid source ships.
#
# Requires node. Run from anywhere:
#     data-raw/figures/render-architecture.sh
#
# mermaid-cli is pinned: an unpinned `npx` picks up whatever is current, and
# the rendered SVG changes across mermaid releases. The committed SVG was
# rendered with the version below; bump it deliberately, in its own commit,
# together with a re-render.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "${here}/../.." && pwd)"

MERMAID_CLI_VERSION="11.17.0"

npx -y "@mermaid-js/mermaid-cli@${MERMAID_CLI_VERSION}" \
    --input "${here}/architecture.mmd" \
    --output "${root}/man/figures/README-architecture.svg" \
    --backgroundColor '#ffffff' \
    --width 2400

echo "wrote ${root}/man/figures/README-architecture.svg"
ls -lh "${root}/man/figures/README-architecture.svg"
