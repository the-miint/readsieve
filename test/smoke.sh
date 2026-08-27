#!/usr/bin/env bash
# End-to-end test on synthetic data. No external inputs: the reference set, the reads and
# the truth are all generated here, so this runs anywhere the duckdb-miint binary does.
#
#   bash test/smoke.sh [workdir]
#
# The design mirrors the real problem. One ancestral sequence; TARGET references are close
# relatives of each other, DECOY references are a diverged sister group, and OUTGROUP
# sequences (never indexed) stand in for the ordinary background the filter must not touch.
# Reads are cut from all three.
#
# WHAT THIS CAN AND CANNOT SHOW. Target and decoy reads are cut from the INDEXED sequences,
# so recall on them is a CEILING, not an estimate -- the filter is being asked whether a
# reference recognises its own sequence. Only the outgroup arm tests anything unseen. That
# is the same limitation the smoke test in any such package has, and it is the reason the
# bounds below are loose on recall and tight on the two error rates.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Resolves the duckdb binary, defines the duckdb() wrapper, and checks that miint really
# loads -- before the workdir is cleared, so a missing dependency costs nothing. Exported
# so the readsieve runs below use this exact binary rather than resolving one of their own.
. "$ROOT/lib/duckdb.sh"
export READSIEVE_DUCKDB="$DUCKDB"
WORK="${1:-$ROOT/tmp/smoke}"
rm -rf "$WORK"; mkdir -p "$WORK"

echo ">>> generating a synthetic reference set and reads"
duckdb -box <<SQL
set threads=4;
set preserve_insertion_order=true;

-- Deterministic base from an integer; hash() is UBIGINT so the modulo is cast before use.
create or replace macro nt(x) as (['A','C','G','T'])[1 + ((x % 4)::bigint)];

create or replace macro rnd(n, seed) as (
    select string_agg(nt(hash(seed * 7919 + i)), '' order by i)
    from range(1, n + 1) t(i));

-- Mutate every Nth base to a DIFFERENT base: take the current base's 0-based index
-- and add 1..3, so the result can never equal the original.
create or replace macro mutate(s, every, seed) as (
    select string_agg(
             case when i % every = 0
                  then nt(list_position(['A','C','G','T'], substr(s, i, 1)) - 1
                          + 1 + (hash(seed * 104729 + i) % 3)::bigint)
                  else substr(s, i, 1) end, '' order by i)
    from range(1, len(s) + 1) t(i));

-- Mate 2 must be a real reverse COMPLEMENT; a plain reverse aligns to nothing and would
-- make the paired arm of this test vacuous.
create or replace macro revcomp(s) as (
    select string_agg(case substr(s, i, 1) when 'A' then 'T' when 'C' then 'G'
                                            when 'G' then 'C' else 'A' end, '' order by -i)
    from range(1, len(s) + 1) t(i));

create or replace table anc as select rnd(1500, 1) as s;

-- 12 targets ~2% from the ancestor, 12 decoys off a 12.5%-diverged sister, and 12 unseen
-- outgroups at ~25% that go into NO index.
create or replace table refs as
    select 'T' || i as ref_id, 'target' as class, 'targetgrp' as label,
           mutate((select s from anc), 50, i) as sequence
    from range(1, 13) t(i)
union all
    select 'D' || i, 'decoy', 'decoygrp',
           mutate(mutate((select s from anc), 8, 999), 50, 100 + i)
    from range(1, 13) t(i);

copy (select ref_id, class, label, sequence from refs)
  to '${WORK}/refs.parquet' (format parquet);

create or replace table src as
    select ref_id, class, sequence from refs
union all
    select 'O' || i, 'outgroup', mutate((select s from anc), 4, 5000 + i) from range(1, 13) t(i);

create or replace table reads as
select class || '_' || ref_id || '_' || p as read_id,
       substr(sequence, 100 + p * 90, 150)                    as sequence1,
       revcomp(substr(sequence, 100 + p * 90 + 200, 150))     as sequence2
from src, range(0, 6) t(p)
where len(substr(sequence, 100 + p * 90 + 200, 150)) = 150;

copy (select read_id, sequence1, sequence2,
             [40 for _ in range(len(sequence1))]::utinyint[] as qual1,
             [40 for _ in range(len(sequence2))]::utinyint[] as qual2
      from reads order by read_id)
  to '${WORK}/reads_{ORIENTATION}.fastq.gz' (format fastq, overwrite true);

select class, count(*) as reference_seqs from refs group by 1 order by 1;
select split_part(read_id,'_',1) as truth, count(*) as read_pairs from reads group by 1 order by 1;
SQL

echo
echo ">>> pack"
bash "$ROOT/readsieve" pack --source "$WORK/refs.parquet" --out "$WORK/pack" \
    --label-col label --threads 4

echo
echo ">>> filter"
bash "$ROOT/readsieve" filter --pack "$WORK/pack" \
    --r1 "$WORK/reads_R1.fastq.gz" --r2 "$WORK/reads_R2.fastq.gz" \
    --out-prefix "$WORK/out" --threads 4

echo
echo ">>> confusion matrix by truth class"
duckdb -box <<SQL
create or replace view c as
select split_part(read_id, '_', 1) as truth, is_positive::boolean as pos
from read_csv('${WORK}/out.calls.tsv.gz', delim='\t', header=true);
select truth, count(*) as pairs, count(*) filter (where pos) as called,
       round(count(*) filter (where pos)::double / count(*), 4) as rate
from c group by 1 order by 1;
SQL

echo
echo ">>> asserted bounds"
rate() { duckdb -noheader -list -c "
select coalesce(round(count(*) filter (where is_positive::boolean)::double / nullif(count(*),0), 4), -1)
from read_csv('${WORK}/out.calls.tsv.gz', delim='\t', header=true)
where split_part(read_id, '_', 1) = '$1';"; }

fail=0
check() { # label  actual  op  bound
    if awk -v a="$2" -v b="$4" "BEGIN{exit !(a $3 b)}"; then
        printf '  %-10s %-8s %s %-6s PASS\n' "$1" "$2" "$3" "$4"
    else
        printf '  %-10s %-8s %s %-6s FAIL\n' "$1" "$2" "$3" "$4"; fail=1
    fi
}
# Recall is a ceiling (positives are cut from indexed sequences), so this bound asserts
# the pipeline works end to end, not that it generalises.
check target   "$(rate target)"   ">=" 0.95
# The costliest error the design exists to prevent: a decoy read called positive.
check decoy    "$(rate decoy)"    "==" 0
# Unseen background. Nothing in either index is its relative.
check outgroup "$(rate outgroup)" "==" 0

viol=$(awk -F'\t' 'NR>1{s+=$2} END{print s+0}' "$WORK/out.verify.tsv")
check verify "$viol" "==" 0

# FASTA in, FASTA out. Qualities are the only thing that changes, and nothing downstream
# of the mapper reads them, so every verdict must survive the round trip unchanged -- and
# the outputs must come back as .fasta.gz rather than as an unwritable .fastq.gz.
echo
echo ">>> same reads as FASTA"
for o in R1 R2; do
    duckdb -c "copy (select read_id, sequence1 from read_fastx('${WORK}/reads_${o}.fastq.gz'))
               to '${WORK}/reads_${o}.fasta.gz' (format fasta, overwrite true);"
done
bash "$ROOT/readsieve" filter --pack "$WORK/pack" \
    --r1 "$WORK/reads_R1.fasta.gz" --r2 "$WORK/reads_R2.fasta.gz" \
    --out-prefix "$WORK/fa" --threads 4 | tail -6
for f in kept_R1 kept_R2 removed_R1 removed_R2; do
    [ -s "$WORK/fa.${f}.fasta.gz" ] || { echo "  MISSING: fa.${f}.fasta.gz"; fail=1; }
done
disagree=$(duckdb -noheader -list -c "
select count(*) filter (where a.is_positive is distinct from b.is_positive)
from read_csv('${WORK}/out.calls.tsv.gz', delim='\t', header=true) a
join read_csv('${WORK}/fa.calls.tsv.gz', delim='\t', header=true) b using (read_id);")
check fasta_calls "$disagree" "==" 0

[ "$fail" = "0" ] || { echo "SMOKE FAILED"; exit 1; }
echo "SMOKE PASSED"
