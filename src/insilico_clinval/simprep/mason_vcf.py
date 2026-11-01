"""Mason input VCF: writer and validator.

Constraints verified against ``mason_simulator`` (seqan-v2.4.0, ``vcf_materialization.cpp``)
and test runs:

- one sample column (only the first is read);
- FORMAT ``GT`` only (indel genotype parsing does not stop at ``:``; ``GT:PS`` crashes);
- one ploidy for the whole file (taken from the first record);
- header contigs identical to the FASTA, in the same order (matched by index);
- biallelic SNV, insertion or deletion records only (MNV/delins are skipped silently);
- records sorted; no conflicting parts on one haplotype.
"""

from __future__ import annotations

from pathlib import Path

import pysam

from insilico_clinval.schemas.simulation import SegmentRecord
from insilico_clinval.simprep import SimPrepError
from insilico_clinval.simprep.events import HaplotypeEvents, Part


class MasonVcfWriter:
    def __init__(self, path: Path, sample_id: str, segments: list[SegmentRecord]):
        header = pysam.VariantHeader()
        for seg in segments:
            header.contigs.add(seg.name, length=seg.end - seg.start + 1)
        header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        header.add_sample(sample_id)
        self.sample_id = sample_id
        self.path = path
        self._out = pysam.VariantFile(str(path), "w", header=header)

    def write(
        self, segment: str, pos: int, ref: str, alt: str, gt: tuple[int, ...], vid: str
    ) -> None:
        rec = self._out.new_record(
            contig=segment, start=pos - 1, alleles=(ref, alt), id=None if vid == "." else vid
        )
        rec.filter.add("PASS")
        sample = rec.samples[self.sample_id]
        sample["GT"] = gt
        sample.phased = len(gt) == 2
        self._out.write(rec)

    def close(self) -> None:
        self._out.close()


def validate_mason_vcf(vcf_path: Path, fasta_path: Path) -> int:
    """Assert all Mason input constraints; returns the number of records."""
    errors: list[str] = []
    fai = [line.split("\t")[:2] for line in Path(f"{fasta_path}.fai").read_text().splitlines()]
    fasta_contigs = [(name, int(length)) for name, length in fai]
    events: dict[tuple[str, int], HaplotypeEvents] = {}
    n = 0
    with pysam.VariantFile(str(vcf_path)) as vcf:
        header_contigs = [
            (str(c), int(vcf.header.contigs[c].length or 0)) for c in vcf.header.contigs
        ]
        if header_contigs != fasta_contigs:
            errors.append("header contigs differ from the FASTA index (names, lengths or order)")
        if len(vcf.header.samples) != 1:
            errors.append(f"{len(vcf.header.samples)} sample columns (Mason reads only the first)")
        formats = list(vcf.header.formats)
        if formats != ["GT"]:
            errors.append(f"FORMAT fields {formats} (only GT is allowed)")
        sample = next(iter(vcf.header.samples), None)
        order = {name: i for i, (name, _) in enumerate(fasta_contigs)}
        previous = (-1, 0)
        ploidy: int | None = None
        for r in vcf:
            n += 1
            where = f"{r.chrom}:{r.pos}"
            key = (order.get(r.chrom, -1), r.pos)
            if key < previous:
                errors.append(f"{where}: not sorted")
            previous = key
            alts = r.alts or ()
            if len(alts) != 1 or r.ref is None:
                errors.append(f"{where}: not biallelic")
                continue
            ref, alt = r.ref, alts[0]
            if len(ref) > 1 and len(alt) > 1:
                errors.append(f"{where}: MNV/delins record (Mason skips it)")
            if list(r.format) != ["GT"]:
                errors.append(f"{where}: FORMAT {list(r.format)}")
            gt = r.samples[sample]["GT"] if sample is not None else None
            if gt is None or any(a is None for a in gt):
                errors.append(f"{where}: missing GT")
                continue
            if ploidy is None:
                ploidy = len(gt)
            elif len(gt) != ploidy:
                errors.append(f"{where}: ploidy {len(gt)} differs from the first record ({ploidy})")
            part = Part(r.pos, ref.upper(), alt.upper())
            for hap, allele in enumerate(gt):
                if allele == 0:
                    continue
                hap_events = events.setdefault((r.chrom, hap), HaplotypeEvents(planted=[]))
                if hap_events.background_conflict(part):
                    errors.append(
                        f"{where}: conflicts with an earlier record (haplotype {hap + 1})"
                    )
                hap_events.accept_background(part)
    if errors:
        raise SimPrepError(f"{vcf_path.name} is not a valid Mason input:\n  " + "\n  ".join(errors))
    return n
