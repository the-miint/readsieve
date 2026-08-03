-- The call. A read is POSITIVE (and its pair is dropped) when a target reference beats
-- every decoy by a clear relative margin AND the target alignment is good:
--
--     target_score > decoy_score                                   (strictly)
--     (target_score - decoy_score) / target_score >= rs_min_margin_frac
--     target_identity >= rs_min_identity
--     target_coverage >= rs_min_coverage
--
-- The competitive terms stop the decoy class being stripped out of user data; the
-- absolute terms stop an unrelated sequence whose best hit is a target over a conserved
-- stretch. See README.md for why each is load-bearing and what was measured.
--
-- Variables the driver sets: rs_target_index, rs_decoy_index, rs_manifest, rs_preset,
-- rs_max_secondary, rs_min_identity, rs_min_coverage, rs_min_margin_frac.
-- The driver has already built rs_input (pair_ix, read_id, comment, sequence1/2, qual1/2).

-- ---- one query row per MATE ---------------------------------------------------
-- Mates are mapped INDEPENDENTLY, as separate single-end queries, so a mate's own call
-- never depends on where its partner landed. Paired minimap2 would let one mate be
-- rescued to its partner's reference and manufacture agreement.
create or replace table rs_query as
    select pair_ix::varchar || ':1' as read_id, sequence1
    from rs_input where sequence1 is not null and len(sequence1) > 0
union all
    select pair_ix::varchar || ':2', sequence2
    from rs_input where sequence2 is not null and len(sequence2) > 0;

-- ---- competitive mapping, one pass per class ----------------------------------
-- MEASURED, so it is not a knob here: minimap2's high-occurrence minimizer mask (-f,
-- exposed as align_minimap2's occ_filter) changes no VERDICT on this workload. Swept over
-- 0.01 / 0.5 / 1000 / 100000 / '2000,10000' on 9,178 pairs against 890 target + 757 decoy
-- near-identical 16S: 0 of 45,890 per-pair verdicts changed and no setting reduced the
-- unmeasured rate. It is not inert, though -- the two FRACTIONAL settings move the
-- reject-reason attribution (197 pairs at 0.01, 155 at 0.5) and 0.5 raises the unmeasured
-- rate slightly, 0.09713 -> 0.09719. So the claim is "no kept-or-removed decision changes
-- and no measurement is recovered", not "nothing changes". If the unmeasured rate is high on
-- a dense reference set, the lever is SEEDING (k/w at pack time), not masking -- see README.md.
create or replace table rs_hit as
    select read_id as qkey, reference, tag_as as score, cigar
    from align_minimap2('rs_query', index_path := getvariable('rs_target_index'),
                        preset := getvariable('rs_preset'), eqx := true,
                        max_secondary := getvariable('rs_max_secondary'))
union all
    select read_id, reference, tag_as, cigar
    from align_minimap2('rs_query', index_path := getvariable('rs_decoy_index'),
                        preset := getvariable('rs_preset'), eqx := true,
                        max_secondary := getvariable('rs_max_secondary'));

-- Best hit per (mate, class). CLASS COMES FROM THE REFERENCE NAME, not from which pass
-- found it -- deriving it from the pass makes the result silently wrong the moment both
-- --target-index and --decoy-index point at the same file. Ordered on
-- (score, identity, reference) so the pick among ties is reproducible rather than
-- whichever row arrived first. The CIGAR is parsed once in the inner query and the alias
-- reused, because calling it again in the window ORDER BY re-parses every row.
create or replace table rs_best as
select hit_class, qkey, reference, score, ident, qcov
from (
    select split_part(reference, '|', 1)   as hit_class,
           qkey, reference, score,
           cigar_sequence_identity(cigar)  as ident,
           cigar_query_coverage(cigar)     as qcov
    from rs_hit)
qualify row_number() over (partition by qkey, hit_class
                           order by score desc, ident desc, reference) = 1;

-- rs_best holds at most one row per (qkey, hit_class), so these aggregates each collapse
-- a single row and stay deterministic.
create or replace table rs_mate as
select qkey,
       max(score) filter (where hit_class = 'target') as target_score,
       max(ident) filter (where hit_class = 'target') as target_ident,
       max(qcov)  filter (where hit_class = 'target') as target_qcov,
       any_value(reference) filter (where hit_class = 'target') as target_ref,
       max(score) filter (where hit_class = 'decoy')  as decoy_score,
       max(ident) filter (where hit_class = 'decoy')  as decoy_ident,
       max(qcov)  filter (where hit_class = 'decoy')  as decoy_qcov,
       any_value(reference) filter (where hit_class = 'decoy')  as decoy_ref
from rs_best group by qkey;

-- A NULL decoy score means no decoy alignment was REPORTED, which is not the same as
-- "distant from every decoy" -- so the competitive terms cannot be evaluated and are
-- SKIPPED, leaving the absolute terms to carry the call. This is the lenient reading and
-- it is deliberate; the summary reports the rate because it is usually a large minority
-- of positives, not an edge case. README.md carries the measurement behind that choice.
create or replace macro rs_mate_call(tscore, dscore, tident, tqcov) as
    coalesce(tscore is not null
             and (dscore is null
                  or (tscore > dscore
                      and (tscore - dscore)::double / tscore
                          >= getvariable('rs_min_margin_frac')))
             and tident >= getvariable('rs_min_identity')
             and tqcov  >= getvariable('rs_min_coverage'), false);

-- Why a mate failed, over EXACTLY the quantities rs_mate_call tests, in the order it
-- tests them. Per MATE, not per pair: a per-pair version comparing greatest(target) to
-- greatest(decoy) across both mates pits mate 2's target score against mate 1's decoy
-- score, a pairing the call never makes.
-- 'unmapped' = no alignment row in EITHER pass. 'no_target_hit' = aligned to a decoy but
-- to no target, which is a measurement.
create or replace macro rs_mate_reason(tscore, dscore, tident, tqcov) as
    case
      when tscore is null and dscore is null then 'unmapped'
      when tscore is null                    then 'no_target_hit'
      when dscore is not null and tscore < dscore then 'decoy_won'
      -- a tie means the read is equally consistent with both classes; the comparison is
      -- strict, so ties KEEP the read. They are counted, not hidden.
      when dscore is not null and tscore = dscore then 'tied_with_decoy'
      when dscore is not null
           and (tscore - dscore)::double / tscore
               < getvariable('rs_min_margin_frac')   then 'thin_margin'
      else 'below_threshold'
    end;

-- How far a mate got before failing, so a pair reports the reason of its BETTER mate. A
-- mate that does not exist (single-end) ranks -1 and never wins.
create or replace macro rs_reason_rank(reason) as
    case reason when 'unmapped' then 0 when 'no_target_hit' then 1 when 'decoy_won' then 2
                when 'tied_with_decoy' then 3 when 'thin_margin' then 4
                when 'below_threshold' then 5 else -1 end;

-- ---- per-PAIR verdict: EITHER mate positive drops the pair --------------------
-- 'called' is itself a reject_reason value so the buckets PARTITION the pairs: a positive
-- pair must not also land in a reject bucket. rs_verify asserts that.
create or replace table rs_call as
select *,
       case when is_positive then 'called'
            when rs_reason_rank(r1_reason) >= rs_reason_rank(r2_reason) then r1_reason
            else r2_reason end as reject_reason
from (
select i.pair_ix, i.read_id, i.comment, (i.sequence2 is not null) as paired,
       m1.target_ref as r1_target_ref, m1.target_score as r1_target_score,
       m1.target_ident as r1_target_ident, m1.target_qcov as r1_target_qcov,
       m1.decoy_ref as r1_decoy_ref, m1.decoy_score as r1_decoy_score,
       m1.decoy_ident as r1_decoy_ident, m1.decoy_qcov as r1_decoy_qcov,
       rs_mate_call(m1.target_score, m1.decoy_score, m1.target_ident, m1.target_qcov) as r1_positive,
       m2.target_ref as r2_target_ref, m2.target_score as r2_target_score,
       m2.target_ident as r2_target_ident, m2.target_qcov as r2_target_qcov,
       m2.decoy_ref as r2_decoy_ref, m2.decoy_score as r2_decoy_score,
       m2.decoy_ident as r2_decoy_ident, m2.decoy_qcov as r2_decoy_qcov,
       rs_mate_call(m2.target_score, m2.decoy_score, m2.target_ident, m2.target_qcov) as r2_positive,
       (rs_mate_call(m1.target_score, m1.decoy_score, m1.target_ident, m1.target_qcov)
        or rs_mate_call(m2.target_score, m2.decoy_score, m2.target_ident, m2.target_qcov)) as is_positive,
       rs_mate_reason(m1.target_score, m1.decoy_score, m1.target_ident, m1.target_qcov) as r1_reason,
       case when i.sequence2 is not null
            then rs_mate_reason(m2.target_score, m2.decoy_score, m2.target_ident, m2.target_qcov) end as r2_reason
from rs_input i
left join rs_mate m1 on m1.qkey = i.pair_ix::varchar || ':1'
left join rs_mate m2 on m2.qkey = i.pair_ix::varchar || ':2');

-- ---- what the removed reads matched -------------------------------------------
create or replace view rs_removed_labels as
select coalesce(m.label, '(unresolved)') as label, count(*) as pairs
from rs_call c
left join read_parquet(getvariable('rs_manifest')) m
       on m.ref_name = coalesce(case when c.r1_positive then c.r1_target_ref end,
                                case when c.r2_positive then c.r2_target_ref end)
where c.is_positive group by all;

-- ---- self-checks the caller must act on ----------------------------------------
-- Every check here has to be able to FAIL on real data, or it is decoration. A row-count
-- comparison between rs_call and rs_input, for instance, cannot fail: rs_call is
-- rs_input left-joined to a table unique on the join key, so it is equal by construction.
create or replace table rs_verify as
with names as (
    select unnest([r1_target_ref, r1_decoy_ref, r2_target_ref, r2_decoy_ref]) as rn
    from rs_call),
tests(test, violations) as (
    select 'pair_ix_unique',
        (select count(*) - count(distinct pair_ix) from rs_input)

    union all
    -- Every mate result must belong to a real input pair. Written as an EQUALITY on the
    -- parsed key: the obvious `m.qkey in (i.pair_ix||':1', i.pair_ix||':2')` puts input
    -- columns in a computed IN-list, which no planner can hash-join, and it degrades to a
    -- nested loop. On a 130k x 317k run that ONE check took 12m50s of a 12m55s job while
    -- the mapping took ~1s. pair_ix::varchar never contains ':', so a well-formed qkey
    -- splits into exactly two parts and this is equivalent, not merely similar.
    select 'mate_keys_resolve',
        (select count(*) from rs_mate m
         where len(str_split(m.qkey, ':')) <> 2
            or split_part(m.qkey, ':', 2) not in ('1', '2')
            or not exists (select 1 from rs_input i
                           where i.pair_ix::varchar = split_part(m.qkey, ':', 1)))

    union all
    select 'reject_reason_partitions_calls',
        (select count(*) from rs_call
         where (is_positive and reject_reason <> 'called')
            or (not is_positive and reject_reason = 'called')
            or reject_reason is null)

    union all
    select 'reject_reason_agrees_with_mates',
        (select count(*) from rs_call
         where not is_positive and (r1_positive or coalesce(r2_positive, false)))

    union all
    select 'verdict_not_null',
        (select count(*) from rs_call where is_positive is null)

    union all
    -- A name that does not resolve means the class read off it is not a class the pack
    -- shipped. NOT EXISTS rather than NOT IN: one NULL ref_name in the manifest makes
    -- NOT IN unknown for every row, so the check would silently report zero.
    select 'ref_names_resolve',
        (select count(*) from names n
         where n.rn is not null
           and not exists (select 1 from read_parquet(getvariable('rs_manifest')) m
                           where m.ref_name = n.rn))

    union all
    -- A pack with one class is the dangerous one: the competitive term is dead and every
    -- read that hits anything is called. Checked against the MANIFEST, not this sample's
    -- hits, because a sample may legitimately produce no decoy alignment at all.
    select 'manifest_has_both_classes',
        (select ((count(*) filter (where class = 'target') = 0)::int
               + (count(*) filter (where class = 'decoy')  = 0)::int)
         from read_parquet(getvariable('rs_manifest')))

    union all
    select 'manifest_ref_name_not_null',
        (select count(*) from read_parquet(getvariable('rs_manifest')) where ref_name is null)

    union all
    select 'positive_has_target_hit',
        (select count(*) from rs_call
         where is_positive and coalesce(r1_target_ref, r2_target_ref) is null)

    union all
    -- The thresholds are resolved by the driver and consumed by the macro. If those ever
    -- read different values -- a variable set after expansion, an env var that did not
    -- survive -- every call would be made against one threshold and reported against
    -- another with nothing visibly wrong. So assert against the floors the summary prints.
    select 'positives_meet_identity_floor',
        (select count(*) from rs_call
         where (r1_positive and r1_target_ident < getvariable('rs_min_identity'))
            or (coalesce(r2_positive, false) and r2_target_ident < getvariable('rs_min_identity')))

    union all
    select 'positives_meet_coverage_floor',
        (select count(*) from rs_call
         where (r1_positive and r1_target_qcov < getvariable('rs_min_coverage'))
            or (coalesce(r2_positive, false) and r2_target_qcov < getvariable('rs_min_coverage')))

    union all
    -- A malformed threshold reaching SQL as NULL would make every comparison unknown and
    -- every call false, silently.
    select 'thresholds_sane',
        (select (getvariable('rs_min_identity') is null or getvariable('rs_min_identity') <= 0
                 or getvariable('rs_min_identity') > 1
                 or getvariable('rs_min_coverage') is null or getvariable('rs_min_coverage') <= 0
                 or getvariable('rs_min_coverage') > 1
                 or getvariable('rs_min_margin_frac') is null
                 or getvariable('rs_min_margin_frac') < 0
                 or getvariable('rs_min_margin_frac') >= 1)::int)
)
select test, violations, case when violations = 0 then 'PASS' else 'FAIL' end as status
from tests order by test;
