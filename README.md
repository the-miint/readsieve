<p align="center">
  <img src="docs/readsieve-logo.png" width="330"
       alt="readsieve — target · decoy · DNA filtering">
</p>

Competitive-mapping read filter. You give it two sets of reference sequences — **target**
and **decoy** — and it removes read pairs that clearly look more like target than decoy.

It does not know what the two sets mean. Chloroplast vs cyanobacteria, host vs microbiome,
vector vs insert, one species vs its sister group: same tool, different pack.

Bash plus DuckDB with the [miint](https://github.com/the-miint/duckdb-miint) extension. No
Python, no conda, no pipeline framework, and no SQL required to use it.

## Install

```bash
./install.sh
```

That is the whole dependency setup. It fetches the DuckDB v1.5.4 CLI into `./bin/duckdb`
and installs `miint` into it from `https://ftp.microbio.me/pub/miint`, then checks that the
functions readsieve calls are actually there. It is idempotent, so re-running it is cheap,
and `DUCKDB_DEST=~/.local/bin ./install.sh` puts the CLI somewhere else.

You do not have to use it. readsieve takes the first working duckdb it finds, in order:

| candidate | what it is |
|---|---|
| `$READSIEVE_DUCKDB` | an explicit path; never silently overridden |
| `./bin/duckdb` | where `install.sh` puts it |
| `duckdb` on `PATH` | any stock CLI that has `miint` installed |
| `../duckdb-miint/build/release/duckdb` | a sibling source build, extension linked in |

A stock CLI plus the published extension and a duckdb-miint source build are equivalent
here — the same pack and the same filter come out of both. Whichever it picks, readsieve
verifies up front that `miint` loads and provides the functions it needs, so a missing or
outdated dependency fails immediately with a message that names the fix, rather than as a
SQL error partway through a run.

Intel Macs are not supported: `miint` publishes `linux_amd64`, `linux_arm64` and
`osx_arm64` builds, and `install.sh` says so rather than downloading a CLI that cannot
load it.

## Quick start

Two FASTA files in, filtered reads out:

```bash
readsieve pack   --target targets.fasta --decoy decoys.fasta --out mypack
readsieve filter --pack mypack --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
                 --out-prefix out/sampleA
```

That is the whole common case. `pack` builds a reference pack once; `filter` runs per
sample and never touches the source files again.

Sequence ids come from the FASTA ids. If a description follows the id, its first word
becomes a **label** used for reporting — so with headers like `>NC_001320 Oryza`, the run
tells you which labels your removed reads matched, which is how you notice that everything
you stripped came from one reference.

## Other input formats

If your references already live in a table, point at it and name the columns:

```bash
readsieve pack --source refs.parquet --out mypack \
    --id-col accession --class-col kind --sequence-col seq --label-col genus \
    --target-class chloroplast --decoy-class cyanobacteria
```

`--source` takes `.parquet`, `.csv`, `.tsv` (gzip fine), or a `.duckdb` file with
`--table NAME`. Rows whose class is neither `--target-class` nor `--decoy-class` are ignored,
and `--where` can narrow the source further.

There is also `--source-sql` for arbitrary queries, if you want to join or filter something
complicated on the way in. You never need it for normal use.

**Which references deserve to be in the pack is your decision.** That is a
question about your data and your risk tolerance, so it belongs in how you prepare the
input, not in here.

A pack is self-contained — two minimap2 indices, `manifest.parquet`, and `pack.tsv`
recording the settings the indices were built with. Copy the directory and it works
anywhere.

## The call

A read pair is **positive** — and dropped — when **either mate** satisfies all of:

```
target_score > decoy_score                                (strictly; ties keep the read)
(target_score - decoy_score) / target_score >= min_margin_frac
target_identity >= min_identity
target_coverage >= min_coverage
```

Defaults: identity 0.95, coverage 0.90, margin 0.20. **These are a starting point, not a
recommendation** — the right values depend on how close your two classes are, and the sweep
that justifies a particular choice is yours to run. `P.calls.tsv.gz` keeps every raw
measurement so you can re-sweep thresholds without re-mapping anything.

What each term is for:

- **The competitive terms** are what stop the decoy class being stripped out of your data.
- **The absolute terms** stop an unrelated sequence whose best hit happens to be a target
  over a short conserved stretch. Of the two, **coverage usually does far more work than
  identity**: a short conserved window matches at very high identity, so identity alone is
  a weak filter and coverage is what rejects partial homology.
- **The margin is relative, not absolute.** minimap2's alignment score scales with read
  length, so a fixed margin silently loosens as reads get longer.

## Six things to watch out for

Each of these produced a wrong answer that looked right.

1. **Unmapped is not negative evidence.** "no hit" means *not measured*, not *distant*. A missing measurement read as a
   confident zero inverts the safe default. `mates_unmeasured_frac` is in the summary for
   this reason: a high value invalidates the negative calls rather than confirming them.

2. **A null decoy score SKIPS the competitive terms** For those
   reads the call is identity + coverage only. This is deliberate — it keeps genuinely
   divergent target reads from being silently discarded — but it is usually a large
   minority of positives, not an edge case, so `positive_no_decoy_measured` is reported.
   On one real dataset it was 53% of positives; the false-positive rate in that stratum was
   *lower* than in the competitively tested one, but the risk of stripping a genuine decoy
   was about 10× higher, because a decoy with no representation in the pack is exactly the
   one that cannot defend itself.

3. **Two passes, one per class — not one combined index.** A single pass reports only a
   bounded number of secondary alignments, so when several target references tie ahead of
   the best decoy, *no decoy row is emitted at all* — and by (1) that is indistinguishable
   from distance. Two passes cost the same total work and give a per-class best for every
   read, so the margin is always a measured number.

4. **Class comes from the reference name.** Names are
   `class|label|id` and the filter parses the class out. Deriving it from the pass would
   make everything silently wrong the moment someone pointed both indices at the same file.
   `pack` rejects a `|` in any id or label for this reason, and the filter asserts that
   every name it emits resolves in the manifest.

5. **If your unmeasured rate is high, consider minimap2 k/w parameters at packing.**
   On 183 near-identical full-length sequences compared against each other, unmeasured went 47.5% → 42.6% (k=15) →
   **16.4%** (k=11, w=5), roughly tripling the usable alignments (487 → 1,460), and those
   recovered alignments are real — median query coverage 0.998 at median identity 0.966 —
   rather than short partial hits. Raising the occupancy limit on top of fine seeding added
   nothing.

6. **k and w are baked into the `.mmi`.** An index only works with a preset whose k/w
   match, so `pack.tsv` records the preset and `filter` reads it from the pack rather than
   defaulting and hoping. Reported identity is also a **lower bound**, so treat it as a
   floor rather than a fine discriminator — and never gate on `mapq`: near-identical
   references tie constantly, so it is 0 for most true hits.

## Outputs

`pack --out PACK` writes four files, and all four are needed:

| file | what |
|---|---|
| `PACK/target.mmi` | minimap2 index over the target class |
| `PACK/decoy.mmi` | minimap2 index over the decoy class |
| `PACK/manifest.parquet` | one row per reference: `ref_name`, `class`, `ref_id`, `label`, `length`. Every reference name a run can emit has to resolve here, and `filter` asserts it |
| `PACK/pack.tsv` | the settings the indices were built with — preset, k, w, counts, median length. `filter` reads the preset back out of this, because k and w are baked into the `.mmi` and an index is only usable with a preset whose k and w match |

`filter --out-prefix P` writes five:

| file | what |
|---|---|
| `P.kept[_R1/_R2].fastq.gz` | reads that survive |
| `P.removed[_R1/_R2].fastq.gz` | reads removed (written by default — a filter you cannot audit is a filter you cannot trust); suppress with `--no-keep-removed` |
| `P.calls.tsv.gz` | one row per pair: the verdict, the reject reason, **and every raw measurement behind it** — per-mate target and decoy reference, score, identity and coverage. This is the file to re-sweep thresholds against without re-mapping |
| `P.summary.tsv` | every number the run printed, as `metric`/`value` |
| `P.verify.tsv` | self-checks, one row each with a violation count |

Reads go in as FASTQ or FASTA (gzipped or not) and come back out in the same format, so a
FASTA sample yields `P.kept.fasta.gz`. Qualities are carried through untouched; nothing
downstream of the mapper reads them, and the smoke test asserts that the same reads as
FASTA and as FASTQ produce identical verdicts.

Nothing else is left behind. Both subcommands write only the files above; `filter` uses
`tmp/` for DuckDB spill and one scratch file, and removes the scratch file even when the
run fails.

### What a run prints

The report goes to **stdout** and everything that went wrong goes to **stderr**, so
`readsieve filter … > run.txt` keeps the report and still shows failures on the terminal.
Both subcommands print the same shape: what they were given, then what came of it.

```
readsieve filter
  reads       sampleA_R1.fastq.gz
              sampleA_R2.fastq.gz
  pack        chloro_pack
  format      fastq, paired, preset sr
  thresholds  identity>=0.95  coverage>=0.90  margin>=0.20  secondary<=5

verdicts                                   pairs   share
  removed                                    720   33.3%
  kept                                     1,440   66.7%

why each pair landed there                 pairs   share
  called (removed)                           720   33.3%
  unmapped                                   720   33.3%
  no_target_hit                              720   33.3%
  tied_with_decoy                              0    0.0%
  decoy_won                                    0    0.0%
  thin_margin                                  0    0.0%
  below_threshold                              0    0.0%

what the removed reads matched             pairs   share
  Oryza                                      720  100.0%

evidence
  median read length                         150
  mate queries                             4,320
  mates with an alignment                  2,880   66.7%
  mates unmeasured                         1,440   33.3%  not evidence of a negative
  removed with no decoy measured             720  100.0%  identity and coverage only
  removed on R1 only                           0    0.0%
  removed on R2 only                           0    0.0%
  removed on both mates                      720  100.0%

self-checks
  checks run                                  12
  violations                                   0

outputs
  sampleA.kept_R1.fastq.gz                 1,440  records
  …
```

Four of those sections are worth reading closely:

**why each pair landed there.** The reject reasons partition the pairs, so the column sums
to the input: `called` · `unmapped` · `no_target_hit` · `decoy_won` · `tied_with_decoy` ·
`thin_margin` · `below_threshold`. Ties are counted separately because they are a real
cost, not a rounding detail.

**what the removed reads matched.** The labels the removed reads hit, largest first. This
is how you notice that everything you stripped came from a single reference — the top 15
are shown, and the rest are in `P.calls.tsv.gz`.

**evidence.** `mates unmeasured` is the one number here that can invalidate a run. An
unmeasured mate is not evidence of a negative, so a high share means the negative calls are
unsupported rather than confirmed; rebuild the pack with a smaller `--k` if it is high.
`removed with no decoy measured` counts positives whose competitive term was skipped
because no decoy aligned at all — for those reads the call rests on identity and coverage
alone, which is usually a large minority and is reported rather than hidden.

**self-checks.** These run in the same process as the work, and the run exits non-zero on
any violation. They are written so they can actually fail — a row-count comparison between
input and output, for instance, is deliberately absent because it is equal by construction
and could never catch anything. The record counts under `outputs` are read back off the
*written* files, which is what catches an output that was truncated or never written.

### When something is wrong

`pack` validates the reference set before it builds anything, and a refusal names both the
problem and the fix:

```
Invalid Input Error:
  sequence_has_gap_or_space (1)
    this is an ALIGNMENT, not sequences; degap it before packing
readsieve: pack: refused, see above -- nothing was written to mypack
```

That check runs ahead of the index build, so "nothing was written" is literal: no
half-built pack is left on disk, and an output directory created for the attempt is removed
again. A pack directory that exists is a pack directory that was finished.

## Test

```bash
bash test/smoke.sh
```

Self-contained: it generates the reference set, the reads and the truth, so it runs
anywhere the binary does. Target and decoy reads are cut from the *indexed* sequences, so
recall on them is a **ceiling, not an estimate** — the filter is only being asked whether a
reference recognises its own sequence. Only the outgroup arm tests anything unseen. A real
generalisation estimate needs a lineage holdout, which is a property of your data rather
than of this tool.

## Where this came from

The design, and every trap in the list above, was extracted from an in-progress 16S
reference-construction project, where this filter shape was developed to strip chloroplast
reads from amplicon data. Nothing here is tied to that project or to any published
reference release, and the sequences that filter used are not necessarily part of any
release either — readsieve is the mechanism only, with the biology and the reference
selection left to the caller.

As an equivalence check, a pack built from that project's chloroplast and cyanobacterial
reference sets reproduces its filter on 9,178 synthetic read pairs with **zero disagreeing
verdicts**, via both the FASTA and the tabular input paths. That establishes that the
extraction is faithful. It is not an independent validation of the filter's accuracy, which
depends entirely on the pack you build.
