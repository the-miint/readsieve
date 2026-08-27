#!/usr/bin/env bash
#
# install.sh -- set up readsieve's only dependency: a DuckDB CLI that can load miint.
#
# Idempotent; anything already present is detected and skipped. It installs:
#   1. the DuckDB v1.5.4 CLI, into ./bin/duckdb
#   2. the 'miint' extension, from https://ftp.microbio.me/pub/miint
#
# Usage:
#   ./install.sh
#   DUCKDB_DEST=~/.local/bin ./install.sh   # put the duckdb CLI somewhere else
#
# readsieve finds ./bin/duckdb by itself, so on the default path there is nothing to
# configure afterwards. If you already build duckdb-miint from source you do not need
# this at all -- readsieve will find a sibling ../duckdb-miint build.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The DuckDB version must match the one miint was built against; bump both together, and
# keep this in lockstep with lib/duckdb.sh.
DUCKDB_VERSION="v1.5.4"
MIINT_REPO_URL="https://ftp.microbio.me/pub/miint"
DUCKDB_DEST="${DUCKDB_DEST:-${SCRIPT_DIR}/bin}"

c_grn=$'\033[1;32m'; c_red=$'\033[1;31m'; c_off=$'\033[0m'
log()  { printf '%s[install]%s %s\n' "$c_grn" "$c_off" "$*"; }
die()  { printf '%s[install] ERROR:%s %s\n' "$c_red" "$c_off" "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1; }

case "${1:-}" in
  "")        : ;;
  -h|--help) sed -n '2,15p' "$0" | sed 's/^# \?//'; exit 0 ;;
  *)         die "unknown option: $1" ;;
esac

need curl || die "required tool 'curl' not found on PATH"
OS="$(uname -s)"; ARCH="$(uname -m)"

# One scratch dir for downloads (project-local, never /tmp), cleaned on exit.
WORK="$(mktemp -d "${SCRIPT_DIR}/.install.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

#############################################################################
# 1. the DuckDB CLI
#############################################################################
duckdb_ok() { [ -x "$1" ] && [ "$("$1" --version 2>/dev/null | awk '{print $1}')" = "$DUCKDB_VERSION" ]; }

DUCKDB_BIN=""
for cand in "$DUCKDB_DEST/duckdb" "$(command -v duckdb 2>/dev/null || true)"; do
  if [ -n "$cand" ] && duckdb_ok "$cand"; then DUCKDB_BIN="$cand"; break; fi
done

if [ -n "$DUCKDB_BIN" ]; then
  log "DuckDB $DUCKDB_VERSION already present: $DUCKDB_BIN"
else
  case "${OS}-${ARCH}" in
    Linux-x86_64)              asset="duckdb_cli-linux-amd64.zip" ;;
    Linux-aarch64|Linux-arm64) asset="duckdb_cli-linux-arm64.zip" ;;
    # miint publishes linux_amd64, linux_arm64 and osx_arm64 -- there is no osx_amd64
    # build, so an Intel Mac is refused here rather than after a 40 MB download.
    Darwin-arm64)              asset="duckdb_cli-osx-universal.zip" ;;
    Darwin-x86_64)             die "miint publishes no osx_amd64 build; use an arm64 Mac or Linux" ;;
    *) die "no prebuilt DuckDB asset for ${OS}-${ARCH}; install ${DUCKDB_VERSION} manually and re-run" ;;
  esac
  need unzip || die "'unzip' is required to install the DuckDB CLI"
  mkdir -p "$DUCKDB_DEST"
  url="https://github.com/duckdb/duckdb/releases/download/${DUCKDB_VERSION}/${asset}"
  log "Downloading DuckDB ${DUCKDB_VERSION} ($asset)..."
  curl -fsSL "$url" -o "$WORK/duckdb.zip" || die "DuckDB download failed: $url"
  unzip -qo "$WORK/duckdb.zip" -d "$WORK"
  install -m 0755 "$WORK/duckdb" "$DUCKDB_DEST/duckdb"
  DUCKDB_BIN="$DUCKDB_DEST/duckdb"
  log "Installed DuckDB CLI -> $DUCKDB_BIN"
fi

#############################################################################
# 2. the miint extension
#############################################################################
log "Installing the 'miint' extension (unsigned, from $MIINT_REPO_URL)..."
"$DUCKDB_BIN" -unsigned \
    -c "INSTALL miint FROM '${MIINT_REPO_URL}'; LOAD miint;" \
  || die "miint install failed -- confirm ${MIINT_REPO_URL}/${DUCKDB_VERSION}/ publishes a build for ${OS}-${ARCH} and the host is reachable"

# Same check lib/duckdb.sh makes at run time: a loadable miint is not enough, it has to be
# new enough to provide what the SQL calls.
have=$("$DUCKDB_BIN" -unsigned -cmd "LOAD miint;" -noheader -list -c "
    select count(distinct function_name) from duckdb_functions()
    where function_name in ('read_fastx', 'align_minimap2', 'save_minimap2_index',
                            'cigar_sequence_identity', 'cigar_query_coverage');" 2>/dev/null) || have=""
[ "$have" = "5" ] || die "miint loaded but is missing functions readsieve needs (${have:-0} of 5); the published build for ${DUCKDB_VERSION} may be older than this checkout expects"
log "miint installed; all 5 required functions present"

#############################################################################
# done
#############################################################################
printf '\n%s[install] readsieve is ready.%s\n\n' "$c_grn" "$c_off"
printf '  DuckDB CLI : %s  (%s, miint loaded)\n\n' "$DUCKDB_BIN" "$DUCKDB_VERSION"
if [ "$DUCKDB_BIN" = "${SCRIPT_DIR}/bin/duckdb" ]; then
  printf 'readsieve finds this automatically. Try it:\n\n'
  printf '  %s/test/smoke.sh\n\n' "$SCRIPT_DIR"
else
  printf 'Point readsieve at it:\n\n'
  printf '  export READSIEVE_DUCKDB="%s"\n\n' "$DUCKDB_BIN"
fi
