"""mason_simulator command line, fragment count and seed.

``mason_simulator`` distributes the requested fragments over the contigs of its FASTA in proportion
to their length and uniformly over the haplotypes (seqan-v2.4.0, ``mason_simulator.cpp``), so a
paired-end run of ``n`` fragments gives a mean depth of ``2 · n · read_length / L`` over all
haplotypes together. The profile coverage is the depth of a diploid locus; a haploid segment
(male chrX/chrY) gets half of it:

    n = ceil(coverage / 2 · ploidy · L / (2 · read_length))

``L`` includes N bases (Mason places no fragments in gaps), so contigs with long gaps are
simulated at slightly higher depth in their non-N part. The output is a deterministic function of
``--seed``, ``--num-threads`` and ``--chunk-size`` (one RNG per thread).
"""

from __future__ import annotations

import math

from insilico_clinval.rng import derive_int
from insilico_clinval.schemas.simulation import SimulationProfile

MASON_SEED_MODULUS = 2**31 - 1


def fragment_count(profile: SimulationProfile, length: int, ploidy: int) -> int:
    return max(1, math.ceil(profile.coverage / 2 * ploidy * length / (2 * profile.read_length)))


def mason_seed(seed: int, sample_id: str, group: str) -> int:
    return derive_int(seed, "mason", f"{sample_id}|{group}") % MASON_SEED_MODULUS


_ILLUMINA_OPTIONS: dict[str, str] = {
    "prob_mismatch": "--illumina-prob-mismatch",
    "prob_mismatch_scale": "--illumina-prob-mismatch-scale",
    "prob_mismatch_begin": "--illumina-prob-mismatch-begin",
    "prob_mismatch_end": "--illumina-prob-mismatch-end",
    "prob_insert": "--illumina-prob-insert",
    "prob_deletion": "--illumina-prob-deletion",
}


def mason_args(
    profile: SimulationProfile,
    *,
    fasta: str,
    vcf: str | None,
    fragments: int,
    seed: int,
    read_name_prefix: str,
    out_r1: str,
    out_r2: str,
    out_alignment: str,
) -> list[str]:
    """Arguments of one mason_simulator run (paths as given, normally relative to its cwd)."""
    frag = profile.fragment
    args = ["-ir", fasta]
    if vcf is not None:
        args += ["-iv", vcf]
    args += [
        "-n",
        str(fragments),
        "--seed",
        str(seed),
        "--num-threads",
        str(profile.threads),
        "--chunk-size",
        str(profile.chunk_size),
        "--seq-technology",
        "illumina",
        "--illumina-read-length",
        str(profile.read_length),
        "--fragment-size-model",
        frag.model,
        "--fragment-mean-size",
        str(frag.mean),
        "--fragment-size-std-dev",
        str(frag.sd),
    ]
    if frag.min is not None:
        args += ["--fragment-min-size", str(frag.min)]
    if frag.max is not None:
        args += ["--fragment-max-size", str(frag.max)]
    errors = profile.illumina.model_dump()
    for key, option in _ILLUMINA_OPTIONS.items():
        if errors[key] is not None:
            args += [option, repr(float(errors[key]))]
    args += [
        "--read-name-prefix",
        read_name_prefix,
        "-o",
        out_r1,
        "-or",
        out_r2,
        "-oa",
        out_alignment,
    ]
    return args
