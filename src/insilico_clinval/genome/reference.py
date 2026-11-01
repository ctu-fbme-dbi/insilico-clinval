"""Reference FASTA access and contig-name resolution."""

from __future__ import annotations

from pathlib import Path

import pysam


class Reference:
    """Contig order and lengths from the FASTA index; sequence access through pysam."""

    def __init__(self, fasta: Path):
        fai = Path(f"{fasta}.fai")
        if not fai.is_file():
            raise FileNotFoundError(f"FASTA index not found: {fai} (run `samtools faidx`)")
        self.fasta = fasta
        self.contigs: list[tuple[str, int]] = []
        with fai.open() as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                self.contigs.append((fields[0], int(fields[1])))
        self.lengths = dict(self.contigs)
        self.order = {name: i for i, (name, _) in enumerate(self.contigs)}
        self._fasta: pysam.FastaFile | None = None

    def resolve_contig(self, name: str) -> str | None:
        """Map a source contig name (``1``, ``X``, ``MT``) to the reference name (``chr1`` …)."""
        if name in self.lengths:
            return name
        candidates = [f"chr{name}"]
        if name in ("MT", "M"):
            candidates = ["chrM", "chrMT"]
        elif name.startswith("chr"):
            candidates = [name[3:]]
        for candidate in candidates:
            if candidate in self.lengths:
                return candidate
        return None

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        if self._fasta is None:
            self._fasta = pysam.FastaFile(str(self.fasta))
        return str(self._fasta.fetch(chrom, start0, end0)).upper()

    def sort_key(self, chrom: str, pos: int) -> tuple[int, int]:
        return (self.order[chrom], pos)

    def close(self) -> None:
        if self._fasta is not None:
            self._fasta.close()
            self._fasta = None
