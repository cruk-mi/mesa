#!/usr/bin/env bash
# ============================================================
# resolve_versions.sh
#
# Single place that maps R -> Bioconductor for the whole repo.
#
# It prints three KEY=VALUE lines (R_VERSION, BIOC_VERSION,
# BIOC_RELEASE) that can be appended to $GITHUB_OUTPUT /
# $GITHUB_ENV in CI, or read locally.
#
# Resolution order (DESCRIPTION is the single source of truth):
#   1. R_VERSION   = DESCRIPTION "R (>= X.Y.Z)"  -> "X.Y"
#                    (overridable via versions.env, rarely needed)
#   2. BIOC_VERSION= latest RELEASED Bioconductor for that R,
#                    looked up in https://bioconductor.org/config.yaml
#                    (overridable via versions.env to pin a release)
#   3. BIOC_RELEASE= RELEASE_<BIOC_VERSION with . -> _>  (Docker tag)
#
# No R or extra tooling required: runs on a bare runner with
# curl, grep, sed and POSIX awk, before any container is pulled.
# Deliberately avoids gawk-only awk extensions so it behaves the
# same on macOS (BSD awk) and Linux CI runners.
# ============================================================
set -euo pipefail

# Resolve repo paths relative to this script so it works from any CWD.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
description="${repo_root}/DESCRIPTION"
versions_env="${script_dir}/versions.env"
config_url="${BIOC_CONFIG_URL:-https://bioconductor.org/config.yaml}"

# 1. Optional overrides (R_VERSION / BIOC_VERSION). Absent by design.
if [ -f "${versions_env}" ]; then
    set -a
    # shellcheck disable=SC1090
    . "${versions_env}"
    set +a
fi

# 2. R_VERSION from DESCRIPTION's "R (>= X.Y.Z)" -> "X.Y" unless overridden.
#    R_VERSION_FULL preserves the full patch string (e.g. "4.6.0") so CI
#    workflows don't need to hardcode a ".0" suffix.
if [ -z "${R_VERSION:-}" ]; then
    R_VERSION_FULL="$(grep -oE 'R \(>=[ ]*[0-9]+\.[0-9]+(\.[0-9]+)?' "${description}" \
        | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
    R_VERSION="$(printf '%s\n' "${R_VERSION_FULL}" | grep -oE '^[0-9]+\.[0-9]+')"
fi
if [ -z "${R_VERSION:-}" ]; then
    echo "resolve_versions.sh: could not parse 'R (>= X.Y.Z)' from ${description}" >&2
    exit 1
fi
# If R_VERSION was provided via versions.env override, R_VERSION_FULL may be
# unset; fall back to the bare major.minor in that case.
: "${R_VERSION_FULL:=${R_VERSION}}"

# 3. BIOC_VERSION: latest released Bioc for this R, from config.yaml,
#    unless pinned via versions.env.
if [ -z "${BIOC_VERSION:-}" ]; then
    config="$(curl -fsSL "${config_url}")" || {
        echo "resolve_versions.sh: failed to fetch ${config_url}" >&2
        exit 1
    }

    # Newest released Bioc; anything above it is devel and excluded.
    release_version="$(printf '%s\n' "${config}" \
        | grep -E '^release_version:' | grep -oE '[0-9]+\.[0-9]+' | head -1)"

    # Bioc versions mapped to this R, from the r_ver_for_bioc_ver: block.
    # Lines look like:  "3.23": "4.6"   (key = Bioc, value = R).
    candidates="$(printf '%s\n' "${config}" \
        | sed -n '/^r_ver_for_bioc_ver:/,/^[^[:space:]]/p' \
        | sed -E 's/#.*//' \
        | grep -E "\"[0-9.]+\"[[:space:]]*:[[:space:]]*\"${R_VERSION}\"[[:space:]]*\$" \
        | sed -E 's/^[[:space:]]*"([0-9]+\.[0-9]+)".*/\1/')"

    # Pick the highest candidate that is <= release_version (released, not devel).
    BIOC_VERSION="$(printf '%s\n' "${candidates}" | awk -v rel="${release_version}" '
        function vnum(v,   a) { split(v, a, "."); return a[1] * 1000 + a[2] }
        BEGIN { rn = vnum(rel); best = ""; bn = -1 }
        NF {
            n = vnum($1)
            if (n <= rn && n > bn) { bn = n; best = $1 }
        }
        END { print best }
    ')"
fi
if [ -z "${BIOC_VERSION:-}" ]; then
    echo "resolve_versions.sh: no released Bioconductor found for R ${R_VERSION}" >&2
    exit 1
fi

BIOC_RELEASE="RELEASE_${BIOC_VERSION/./_}"

printf 'R_VERSION=%s\n'      "${R_VERSION}"
printf 'R_VERSION_FULL=%s\n' "${R_VERSION_FULL}"
printf 'BIOC_VERSION=%s\n'   "${BIOC_VERSION}"
printf 'BIOC_RELEASE=%s\n'   "${BIOC_RELEASE}"
