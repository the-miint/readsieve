# Shared DuckDB resolution. Sourced by ./readsieve and ./test/smoke.sh, both of which set
# $ROOT to the checkout first. Defines $DUCKDB and a duckdb() wrapper -- callers go through
# the wrapper, never "$DUCKDB" directly, so every session loads miint the same way.
#
# readsieve needs a DuckDB CLI that can LOAD the miint extension. Two shapes qualify and
# are treated identically:
#   - a stock DuckDB CLI with miint installed from the repository below, which is what
#     ./install.sh sets up in ./bin/duckdb;
#   - a duckdb-miint source build, which links the extension in already.

# Keep in lockstep with install.sh.
RS_DUCKDB_VERSION="v1.5.4"
RS_MIINT_REPO_URL="https://ftp.microbio.me/pub/miint"

rs_die() { echo "readsieve: $*" >&2; exit 2; }

# An explicit READSIEVE_DUCKDB is never silently overridden: if it is wrong the run stops
# here rather than quietly mapping against a different binary than the operator named.
if [ -n "${READSIEVE_DUCKDB:-}" ]; then
    DUCKDB="$READSIEVE_DUCKDB"
    [ -x "$DUCKDB" ] || rs_die "READSIEVE_DUCKDB=$DUCKDB is not an executable"
else
    DUCKDB=""
    # type -P, not command -v: command -v would match the duckdb() wrapper defined below.
    for rs_cand in "$ROOT/bin/duckdb" "$(type -P duckdb || true)" \
                   "${ROOT%/*}/duckdb-miint/build/release/duckdb"; do
        if [ -n "$rs_cand" ] && [ -x "$rs_cand" ]; then DUCKDB="$rs_cand"; break; fi
    done
    unset rs_cand
    [ -n "$DUCKDB" ] || rs_die "no duckdb found in $ROOT/bin, on PATH, or at ${ROOT%/*}/duckdb-miint.
  Run ./install.sh to fetch DuckDB $RS_DUCKDB_VERSION and the miint extension, or set
  READSIEVE_DUCKDB to a duckdb binary that has miint."
fi

# -unsigned is required when miint was installed from a non-default repository; both flags
# are harmless on a duckdb-miint build that has the extension linked in.
duckdb() { "$DUCKDB" -unsigned -cmd "LOAD miint;" "$@"; }

# Fail once, here, instead of inside a query. A duckdb that cannot load miint -- or a miint
# too old for what the SQL calls -- otherwise surfaces as a mid-run parse error that says
# nothing about the cause, after pack/filter has already started writing.
rs_have=$(duckdb -noheader -list -c "
    select count(distinct function_name) from duckdb_functions()
    where function_name in ('read_fastx', 'align_minimap2', 'save_minimap2_index',
                            'cigar_sequence_identity', 'cigar_query_coverage');" 2>/dev/null) || rs_have=""
if [ "$rs_have" != "5" ]; then
    rs_die "$DUCKDB does not provide the miint functions readsieve needs (${rs_have:-0} of 5 present).
  Install the extension:  $DUCKDB -unsigned -c \"INSTALL miint FROM '$RS_MIINT_REPO_URL'\"
  or run ./install.sh, which does that against a pinned DuckDB $RS_DUCKDB_VERSION."
fi
unset rs_have
