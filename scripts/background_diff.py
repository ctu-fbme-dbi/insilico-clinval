"""Difference of two background genomes around the planted variants of existing suites.

Compares the background of an old suite (e.g. GIAB HG002 v4.2.1) with that of a new one (e.g.
v5.0q, with its record filters) without running any pipeline. Each background is read exactly as
``build-truth`` reads it (``BackgroundMerger``, no planted variants). Writes to ``--out``:

- ``bed_sizes.tsv``: per chromosome, bp of the old BED, the new BED, both, old only, new only;
- ``records.tsv``: background records kept per background, by chromosome, variant type and
  stratum of the new suite's ``sources.strata`` (``*`` = all);
- ``planted_context.tsv``: per planted variant of the given runs (``variants.selected.json``):
  inside each BED, and the number of background records of each version within
  ± ``min_distance_bp`` and ± ``--window`` bp;
- ``skipped.tsv``: records skipped per background and reason.

Genotype concordance of the two backgrounds in the BED intersection is a separate hap.py run, e.g.
``hap.py v5.0q.vcf.gz v4.2.1.vcf.gz -f <both.bed> -r <fasta> --engine vcfeval --gender male``.

    python scripts/background_diff.py --old configs/suites/acmg_sf32_full1.yaml \\
        --new configs/suites/difficult_v5q.yaml --run runs/acmg_sf32_short1 ... --out DIR
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
from collections import Counter
from collections.abc import Iterable, Sequence
from pathlib import Path

from insilico_clinval.config import LoadedSuite, load_suite
from insilico_clinval.genome.alleles import classify, ref_span0
from insilico_clinval.genome.bed import IntervalSet
from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.selection import SelectionResult
from insilico_clinval.schemas.suite import BackgroundConfig
from insilico_clinval.truth.background import BackgroundMerger
from insilico_clinval.variants.report import SELECTED_JSON

Intervals = dict[str, list[tuple[int, int]]]


def read_bed(path: Path) -> Intervals:
    """Merged, sorted half-open intervals per chromosome."""
    raw: Intervals = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    merged: Intervals = {}
    for chrom, ivs in raw.items():
        out: list[tuple[int, int]] = []
        for s, e in sorted(ivs):
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
    return merged


def size(ivs: Iterable[tuple[int, int]]) -> int:
    return sum(e - s for s, e in ivs)


def intersect(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    i = j = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def bed_sizes(old: Intervals, new: Intervals, contigs: Sequence[str]) -> list[list[str]]:
    rows = []
    for chrom in contigs:
        a, b = old.get(chrom, []), new.get(chrom, [])
        both = size(intersect(a, b))
        rows.append([chrom, *map(str, (size(a), size(b), both, size(a) - both, size(b) - both))])
    return rows


def _merger(suite: LoadedSuite, reference: Reference) -> BackgroundMerger:
    bg = suite.config.background
    if bg is None:
        raise SystemExit(f"{suite.path}: no background configured")
    return BackgroundMerger(
        config=bg,
        vcf_path=suite.resolve(bg.vcf),
        bed_path=suite.resolve(bg.bed),
        reference=reference,
        seed=suite.config.suite.seed,
        sample_id="background_diff",
        planted=[],
        min_distance_bp=0,
    )


def _label(bg: BackgroundConfig) -> str:
    return f"{bg.name} {bg.release}"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--old", type=Path, required=True, help="suite with the old background")
    parser.add_argument("--new", type=Path, required=True, help="suite with the new background")
    parser.add_argument("--run", type=Path, action="append", default=[], help="run directory")
    parser.add_argument("--window", type=int, default=150, help="second window (read length)")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    old, new = load_suite(args.old), load_suite(args.new)
    assert old.config.background is not None and new.config.background is not None
    reference = Reference(new.resolve(new.config.reference.fasta))
    contigs = [c for c, _ in reference.contigs if c in reference.lengths and "_" not in c]
    strata = {k: IntervalSet(new.resolve(v)) for k, v in sorted(new.config.sources.strata.items())}
    args.out.mkdir(parents=True, exist_ok=True)

    beds = {
        "old": read_bed(old.resolve(old.config.background.bed)),
        "new": read_bed(new.resolve(new.config.background.bed)),
    }
    _write(args.out / "bed_sizes.tsv",
           ["chrom", "old_bp", "new_bp", "both_bp", "old_only_bp", "new_only_bp"],
           bed_sizes(beds["old"], beds["new"], contigs))  # fmt: skip

    positions: dict[str, dict[str, list[int]]] = {}
    counts: Counter[tuple[str, str, str, str]] = Counter()
    skipped: list[list[str]] = []
    for key, suite in (("old", old), ("new", new)):
        merger = _merger(suite, reference)
        per_chrom: dict[str, list[int]] = {}
        for chrom in contigs:
            pos: list[int] = []
            for rec in merger.records(chrom):
                pos.append(rec.pos)
                vtype = "|".join(sorted({classify(rec.ref, a).value for a in rec.alts}))
                start, end = ref_span0(rec.pos, rec.ref)
                names = ["*", *(n for n, s in strata.items() if s.overlaps(chrom, start, end))]
                for name in names:
                    counts[(key, chrom, vtype, name)] += 1
            per_chrom[chrom] = pos
        positions[key] = per_chrom
        bg = suite.config.background
        assert bg is not None
        label = _label(bg)
        skipped += [[key, label, r, str(n)] for r, n in sorted(merger.skipped.items())]
        skipped.append(
            [key, label, "spanning_deletion_to_ref", str(merger.spanning_deletion_to_ref)]
        )
    _write(args.out / "records.tsv", ["background", "chrom", "variant_type", "stratum", "records"],
           [[*k, str(n)] for k, n in sorted(counts.items())])  # fmt: skip
    _write(args.out / "skipped.tsv", ["background", "label", "reason", "records"], skipped)

    min_dist = old.config.constraints.min_distance_bp
    rows: list[list[str]] = []
    for run in args.run:
        result = SelectionResult.model_validate_json((run / SELECTED_JSON).read_text())
        for v in result.variants:
            start, end = ref_span0(v.pos, v.ref)
            row = [run.name, v.variant_id, v.chrom, str(v.pos), v.gene or ""]
            for key in ("old", "new"):
                inside = any(s <= start and end <= e for s, e in beds[key].get(v.chrom, []))
                pos = positions[key].get(v.chrom, [])
                near = [_count(pos, start - w, end + w) for w in (min_dist, args.window)]
                row += ["yes" if inside else "no", *map(str, near)]
            rows.append(row)
    w = args.window
    header = ["run", "variant_id", "chrom", "pos", "gene"]
    for key in ("old", "new"):
        header += [f"{key}_in_bed", f"{key}_bg_within_{min_dist}bp", f"{key}_bg_within_{w}bp"]
    _write(args.out / "planted_context.tsv", header, rows)
    reference.close()


def _count(sorted_pos: list[int], start0: int, end0: int) -> int:
    """Records whose POS (1-based) lies in the half-open 0-based window [start0, end0)."""
    return bisect.bisect_right(sorted_pos, end0) - bisect.bisect_right(sorted_pos, start0)


def _write(path: Path, header: list[str], rows: Iterable[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
