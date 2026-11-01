"""Step 3a: turn the truth VCF into Mason inputs (segments, segment FASTA, Mason VCF)."""


class SimPrepError(RuntimeError):
    """The truth set cannot be represented for Mason (e.g. overlapping planted variants)."""
