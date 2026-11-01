"""Step 2: truth encoding (REF check, genotypes, truth VCF, expected outcomes, sample spec)."""


class TruthError(RuntimeError):
    """The truth set cannot be encoded as configured (REF mismatch, ploidy, phasing …)."""
