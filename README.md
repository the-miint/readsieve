# readsieve

Competitive-mapping read filter. You give it two sets of reference sequences — **target**
and **decoy** — and it removes read pairs that clearly look more like target than decoy.

It does not know what the two sets mean. Chloroplast vs cyanobacteria, host vs microbiome,
vector vs insert, one species vs its sister group: same tool, different pack.

Bash plus a [duckdb-miint](https://github.com/wasade/duckdb-miint) build. No Python, no
conda, no pipeline framework, and no SQL required to use it.

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

If your references already live in a table, point at it and name the columns — still no SQL:

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

**Which references deserve to be in the pack is your decision, not the tool's.** That is a
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

## Six things that will bite you

Each of these produced a wrong answer that looked right.

1. **Unmapped is not negative evidence.** minimap2 emits no row for a query it cannot
   align, so "no hit" means *not measured*, not *distant*. A missing measurement read as a
   confident zero inverts the safe default. `mates_unmeasured_frac` is in the summary for
   this reason: a high value invalidates the negative calls rather than confirming them.

2. **A null decoy score SKIPS the competitive terms, it does not fail them.** For those
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

4. **Class comes from the reference name, never from which pass found it.** Names are
   `class|label|id` and the filter parses the class out. Deriving it from the pass would
   make everything silently wrong the moment someone pointed both indices at the same file.
   `pack` rejects a `|` in any id or label for this reason, and the filter asserts that
   every name it emits resolves in the manifest.

5. **If your unmeasured rate is high, the lever is seeding — not occupancy masking.**
   This one is measured, and it is the opposite of what we assumed. minimap2's
   high-occurrence minimizer mask (`-f`) is a natural suspect on a dense, near-identical
   reference set, and it is **not** the cause. Swept over `0.01 / 0.5 / 1000 / 100000 /
   '2000,10000'` on 9,178 pairs against 890 target + 757 decoy near-identical 16S sequences:
   **not one of the 45,890 per-pair verdicts changed**, and no setting reduced the unmeasured
   rate. So readsieve deliberately exposes **no** flag for it.

   Be precise about what that does and does not say, because an earlier version of this
   README overstated it. The parameter is *live*, not ignored: at the two **fractional**
   settings — which minimap2 reads as a proportion to mask rather than an occurrence count —
   the reject-reason attribution moves (197 pairs at `0.01`, 155 at `0.5`), and `0.5` nudges
   the unmeasured rate the wrong way, 0.09713 → 0.09719. The three integer settings change
   nothing at all. The claim is "masking changes no kept-or-removed decision and recovers no
   measurement", not "masking changes nothing".

   What does move the unmeasured rate is **k/w at pack time**. On 183 near-identical
   full-length sequences compared against each other, unmeasured went 47.5% → 42.6% (k=15) →
   **16.4%** (k=11, w=5), roughly tripling the usable alignments (487 → 1,460), and those
   recovered alignments are real — median query coverage 0.998 at median identity 0.966 —
   rather than short partial hits. Raising the occupancy limit on top of fine seeding added
   exactly nothing.

6. **k and w are baked into the `.mmi`.** An index only works with a preset whose k/w
   match, so `pack.tsv` records the preset and `filter` reads it from the pack rather than
   defaulting and hoping. Reported identity is also a **lower bound**, so treat it as a
   floor rather than a fine discriminator — and never gate on `mapq`: near-identical
   references tie constantly, so it is 0 for most true hits.

## Outputs

Given `--out-prefix P`:

| file | what |
|---|---|
| `P.kept[_R1/_R2].fastq.gz` | reads that survive |
| `P.removed[_R1/_R2].fastq.gz` | reads removed (written by default — a filter you cannot audit is a filter you cannot trust) |
| `P.calls.tsv.gz` | per-pair verdict **and every raw measurement** |
| `P.summary.tsv` | the numbers, also printed |
| `P.verify.tsv` | self-checks, PASS/FAIL per row |

Reject reasons partition the pairs, so the buckets sum to the input: `unmapped` ·
`no_target_hit` · `decoy_won` · `tied_with_decoy` · `thin_margin` · `below_threshold` ·
`called`. Ties are counted separately because they are a real cost, not a rounding detail.

The self-checks run in the same process as the work, and the driver exits non-zero on any
violation. They are written so they can actually fail — a row-count comparison between
input and output, for instance, is deliberately absent because it is equal by construction
and could never catch anything.

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
