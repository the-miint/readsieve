-- Build a reference pack: two minimap2 indices plus a manifest.
--
-- Input is one table, rs_src, that the driver has already defined, yielding ref_id,
-- class ('target'|'decoy'), sequence and label.
--
-- Variables the driver sets: rs_out, rs_preset, rs_k, rs_w.

-- ---- normalise, and reject what would silently corrupt a run -----------------
-- ref_name is 'class|label|ref_id' and the CLASS IS READ OFF THIS NAME at filter time,
-- never off which index a hit came from. That is what stops a run being silently wrong
-- when someone points both --target-index and --decoy-index at the same file.
--
-- Only the FIRST field therefore has to be free of the separator, and it always is:
-- class_not_target_or_decoy below leaves nothing but 'target' and 'decoy'. A '|' inside
-- label or ref_id is harmless, because nothing splits those back out -- the filter joins
-- the manifest on the whole ref_name and reads label and ref_id from their own columns.
-- The one way a stray separator could do damage is by letting two references collide
-- ('a|b' + 'c' and 'a' + 'b|c' concatenate alike), and ref_name_not_unique catches that.
-- Banning '|' outright instead cost nothing on invented data and rejected every real
-- record from sources that use it inside accessions, NCBI-style.
create or replace table rs_ref as
select class || '|' || label || '|' || ref_id as ref_name,
       ref_id, class, label, upper(sequence) as sequence, len(sequence) as length
from rs_src;

-- Each row carries what to DO about it: 'sequence_has_gap_or_space' names a symptom, and
-- the operator staring at it still has to be told that their input is an alignment.
create or replace table rs_reject as
select 'class_not_target_or_decoy' as problem, count(*) as n,
       'class is neither --target-class nor --decoy-class; check --class-col' as remedy
    from rs_ref where class not in ('target', 'decoy')
union all select 'ref_id_null_or_empty', count(*),
       'every reference needs an id; check --id-col, or the FASTA headers'
    from rs_ref where ref_id is null or ref_id = ''
union all select 'sequence_null_or_empty', count(*),
       'empty sequences; drop them with --where, or fix the source'
    from rs_ref where sequence is null or length = 0
union all select 'sequence_has_gap_or_space', count(*),
       'this is an ALIGNMENT, not sequences; degap it before packing'
    from rs_ref where regexp_matches(sequence, '[-. \t]')
union all select 'ref_name_not_unique',
       (select count(*) - count(distinct ref_name) from rs_ref),
       'two references share class+label+ref_id; ids must be unique within a class'
union all select 'no_target_sequences',
       (select (count(*) filter (where class='target') = 0)::int from rs_ref),
       'nothing matched --target-class, so there would be nothing to detect'
-- A pack with no decoys is the dangerous case, not an empty one: every competitive term
-- would be skipped and every read that hits anything would be called positive.
union all select 'no_decoy_sequences',
       (select (count(*) filter (where class='decoy') = 0)::int from rs_ref),
       'nothing matched --decoy-class; with no decoy every hit is called positive';

-- Nothing below this line runs when a reference would corrupt the pack, and nothing above
-- it has written a file -- the indices are built by the next statement, and the driver's
-- COPYs come later still. That ordering is what makes "nothing was written" true. It used
-- to be checked by the driver AFTER the whole script had finished, which left a complete,
-- usable-looking pack directory sitting next to the refusal to have built one.
-- (`.bail on`, set by the driver, is what makes this stop the script rather than just
-- print; the CLI otherwise carries on to the next statement.)
-- The WHERE, not a CASE, is what keeps this silent: with nothing rejected it yields no
-- rows, so error() is never reached and the statement prints nothing at all.
select error(chr(10) || (select string_agg('  ' || problem || ' (' || n || ')' ||
                                           chr(10) || '    ' || remedy,
                                           chr(10) order by problem)
                         from rs_reject where n > 0))
where (select coalesce(sum(n), 0) from rs_reject) > 0;

-- ---- the indices -------------------------------------------------------------
-- One per class, and the filter maps against them separately. A single combined index
-- would decide correctly but reports only a bounded number of secondary alignments, so
-- whenever several target references tie ahead of the best decoy NO decoy row is emitted
-- -- and a missing row is not evidence of distance. Two passes cost the same total work
-- and yield a per-class best for every read, so the margin is always a measured number.
--
-- k and w are BAKED INTO the .mmi. An index is only usable with a preset whose k/w match,
-- which is why the preset is recorded in pack.tsv and asserted at filter time.
create or replace table rs_target as select ref_name as read_id, sequence as sequence1 from rs_ref where class = 'target';
create or replace table rs_decoy  as select ref_name as read_id, sequence as sequence1 from rs_ref where class = 'decoy';

create or replace table rs_index as
    select 'target' as role, success, num_subjects
    from save_minimap2_index('rs_target', getvariable('rs_out') || '/target.mmi',
                             preset := getvariable('rs_preset'),
                             k := getvariable('rs_k'), w := getvariable('rs_w'))
union all
    select 'decoy', success, num_subjects
    from save_minimap2_index('rs_decoy', getvariable('rs_out') || '/decoy.mmi',
                             preset := getvariable('rs_preset'),
                             k := getvariable('rs_k'), w := getvariable('rs_w'));

-- ---- manifest and pack metadata ------------------------------------------------
-- Built as tables here and WRITTEN by the driver: COPY ... TO takes a string literal, not
-- an expression, so getvariable() cannot name the target. (save_minimap2_index above is a
-- table function, whose arguments are ordinary expressions -- hence the asymmetry.)
--
-- Every name the filter can emit must resolve in the manifest, and the filter asserts it.
-- Extra attribution columns can be carried through by selecting them in --source-sql.
create or replace table rs_manifest as
select ref_name, class, ref_id, label, length from rs_ref order by ref_name;

create or replace table rs_packmeta as
    select 'preset' as key, getvariable('rs_preset')::varchar as value
    union all select 'k', coalesce(getvariable('rs_k')::varchar, 'preset_default')
    union all select 'w', coalesce(getvariable('rs_w')::varchar, 'preset_default')
    union all select 'n_target', (select count(*) from rs_ref where class='target')::varchar
    union all select 'n_decoy',  (select count(*) from rs_ref where class='decoy')::varchar
    union all select 'median_ref_length', (select median(length) from rs_ref)::varchar
    union all select 'target_index_subjects', (select num_subjects from rs_index where role='target')::varchar
    union all select 'decoy_index_subjects',  (select num_subjects from rs_index where role='decoy')::varchar;
