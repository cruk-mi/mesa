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

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "${here}/../.." && pwd)"

npx -y @mermaid-js/mermaid-cli \
    --input "${here}/architecture.mmd" \
    --output "${root}/man/figures/README-architecture.svg" \
    --backgroundColor '#ffffff' \
    --width 2400

echo "wrote ${root}/man/figures/README-architecture.svg"
ls -lh "${root}/man/figures/README-architecture.svg"
