"""Steps of the fidelity experiment and assembly of `fidelity.json`.

Each step writes the raw output of its tool under the output directory; :func:`assemble` reads
whatever raw outputs are present, so steps can run separately (or on another machine).

Layout::

    planting.tsv
    qc/regions.samtools.tsv            1-based regions of `samtools stats -t`
    qc/<id>.samtools_stats.txt
    qc/<id>.bcftools_stats.txt
    strata.tsv                         hap.py --stratification
    happy/<id>/<id>.extended.csv (+ other hap.py outputs)
    fidelity.json, fidelity_calling.tsv, fidelity_qc.tsv
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from insilico_clinval._runner import Toolbox
from insilico_clinval.fidelity import LoadedFidelity
from insilico_clinval.fidelity.planting import (
    collect_planting,
    read_planting_tsv,
    summarise_planting,
    write_planting_tsv,
)
from insilico_clinval.fidelity.stats import (
    parse_bcftools_stats,
    parse_happy_extended,
    parse_samtools_stats,
    region_length,
    write_samtools_regions,
)
from insilico_clinval.jsonio import write_json
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.schemas.fidelity import (
    CallingMetric,
    CallsetQc,
    FidelityDataset,
    FidelityResult,
    ReadQc,
    RecallGap,
)

log = logging.getLogger(__name__)

STEPS = ("planting", "qc", "calling")
SCHEMA_VERSION = 1


def samtools_stats_path(out_dir: Path, dataset_id: str) -> Path:
    return out_dir / "qc" / f"{dataset_id}.samtools_stats.txt"


def bcftools_stats_path(out_dir: Path, dataset_id: str) -> Path:
    return out_dir / "qc" / f"{dataset_id}.bcftools_stats.txt"


def happy_prefix(out_dir: Path, dataset_id: str) -> Path:
    return out_dir / "happy" / dataset_id / dataset_id


def regions_path(out_dir: Path) -> Path:
    return out_dir / "qc" / "regions.samtools.tsv"


# -- steps -----------------------------------------------------------------------------------


def run_planting(fid: LoadedFidelity, out_dir: Path) -> Path:
    rows = collect_planting([fid.resolve(p) for p in fid.config.planting])
    return write_planting_tsv(rows, out_dir / "planting.tsv")


def samtools_stats_args(fid: LoadedFidelity, d: FidelityDataset, regions: Path) -> list[str]:
    return [
        "stats", "--reference", str(fid.reference), "-@", str(fid.config.threads),
        "-t", str(regions), str(fid.resolve(d.alignment)),
    ]  # fmt: skip


def bcftools_stats_args(fid: LoadedFidelity, d: FidelityDataset) -> list[str]:
    args = ["stats", "--threads", str(fid.config.threads), "-s", "-", "-T", str(fid.qc_regions)]
    if fid.config.happy.pass_only:
        args += ["-f", "PASS"]
    return [*args, str(fid.resolve(d.vcf))]


def run_qc(
    fid: LoadedFidelity, out_dir: Path, toolbox: Toolbox, manifest: ManifestBuilder
) -> list[Path]:
    (out_dir / "qc").mkdir(parents=True, exist_ok=True)
    regions = regions_path(out_dir)
    write_samtools_regions(fid.qc_regions, regions)
    written: list[Path] = []
    for d in fid.config.datasets:
        mounts = [fid.reference.parent, fid.resolve(d.alignment).parent, regions.parent]
        log.info("samtools stats %s", d.id)
        result = toolbox.run(
            "samtools",
            samtools_stats_args(fid, d, regions),
            step=f"fidelity samtools stats ({d.id})",
            manifest=manifest,
            mounts=mounts,
        )
        path = samtools_stats_path(out_dir, d.id)
        path.write_text(result.stdout, encoding="utf-8")
        written.append(path)
        log.info("bcftools stats %s", d.id)
        result = toolbox.run(
            "bcftools",
            bcftools_stats_args(fid, d),
            step=f"fidelity bcftools stats ({d.id})",
            manifest=manifest,
            mounts=[fid.qc_regions.parent, fid.resolve(d.vcf).parent],
        )
        path = bcftools_stats_path(out_dir, d.id)
        # The first lines repeat the command line (absolute paths); keep only the data lines.
        path.write_text(
            "".join(line + "\n" for line in result.stdout.splitlines() if not line.startswith("#")),
            encoding="utf-8",
        )
        written.append(path)
    return written


def write_strata_tsv(fid: LoadedFidelity, out_dir: Path) -> Path:
    path = out_dir / "strata.tsv"
    path.write_text("".join(f"{n}\t{p}\n" for n, p in fid.strata()), encoding="utf-8")
    return path


def happy_args(
    fid: LoadedFidelity, d: FidelityDataset, prefix: Path, strata: Path | None
) -> list[str]:
    h = fid.config.happy
    args = [
        str(fid.resolve(d.truth_vcf)), str(fid.resolve(d.vcf)),
        "-r", str(fid.reference), "-f", str(fid.confident_bed), "-o", str(prefix),
        "--engine", h.engine, "--threads", str(h.threads),
    ]  # fmt: skip
    if h.engine == "vcfeval":
        args += ["--engine-vcfeval-path", "rtg"]
        if h.vcfeval_template is not None:
            args += ["--engine-vcfeval-template", str(fid.resolve(h.vcfeval_template))]
    if h.pass_only:
        args.append("--pass-only")
    if strata is not None:
        args += ["--stratification", str(strata)]
    if h.gender is not None:
        args += ["--gender", h.gender.value]
    return args


def run_calling(
    fid: LoadedFidelity, out_dir: Path, toolbox: Toolbox, manifest: ManifestBuilder
) -> list[Path]:
    strata = write_strata_tsv(fid, out_dir) if fid.config.stratifications else None
    common = [fid.reference.parent, fid.confident_bed.parent, out_dir]
    common += [p.parent for _, p in fid.strata()]
    if fid.config.happy.vcfeval_template is not None:
        common.append(fid.resolve(fid.config.happy.vcfeval_template).parent)
    written: list[Path] = []
    for d in fid.config.datasets:
        prefix = happy_prefix(out_dir, d.id)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        log.info("hap.py %s", d.id)
        toolbox.run(
            "hap.py",
            happy_args(fid, d, prefix, strata),
            step=f"fidelity hap.py ({d.id})",
            manifest=manifest,
            mounts=[*common, fid.resolve(d.vcf).parent, fid.resolve(d.truth_vcf).parent],
            cwd=prefix.parent,
            stderr_path=prefix.parent / "happy.stderr.log",
        )
        written.append(Path(f"{prefix}.extended.csv"))
    return written


# -- assembly --------------------------------------------------------------------------------


def recall_gaps(
    calling: list[CallingMetric], baseline: str, threshold: float | None
) -> list[RecallGap]:
    """Recall of every non-baseline dataset minus the baseline's, per type and stratum."""
    base = {(m.type, m.stratum): m for m in calling if m.dataset == baseline and m.subtype == "*"}
    gaps: list[RecallGap] = []
    for m in calling:
        if m.dataset == baseline or m.subtype != "*":
            continue
        b = base.get((m.type, m.stratum))
        r, br = m.recall.value, b.recall.value if b is not None else None
        gap = round(r - br, 6) if r is not None and br is not None else None
        flagged = gap is not None and threshold is not None and abs(gap) > threshold
        gaps.append(
            RecallGap(
                dataset=m.dataset,
                type=m.type,
                stratum=m.stratum,
                recall=r,
                baseline_recall=br,
                gap=gap,
                flagged=flagged,
            )
        )
    gaps.sort(key=lambda g: (g.dataset, g.type, g.stratum))
    return gaps


def assemble(fid: LoadedFidelity, out_dir: Path) -> FidelityResult:
    cfg = fid.config
    missing: list[str] = []
    planting = None
    planting_tsv = out_dir / "planting.tsv"
    if planting_tsv.is_file():
        planting = summarise_planting(read_planting_tsv(planting_tsv))
    elif cfg.planting:
        missing.append("planting.tsv")
    read_qc: list[ReadQc] = []
    callset_qc: list[CallsetQc] = []
    calling: list[CallingMetric] = []
    regions = regions_path(out_dir)
    region_bp = region_length(regions) if regions.is_file() else 0
    for d in cfg.datasets:
        sam = samtools_stats_path(out_dir, d.id)
        if sam.is_file() and region_bp:
            text = sam.read_text(encoding="utf-8")
            read_qc.append(parse_samtools_stats(text, dataset=d.id, region_bp=region_bp))
        else:
            missing.append(str(sam.relative_to(out_dir)))
        bcf = bcftools_stats_path(out_dir, d.id)
        if bcf.is_file():
            callset_qc.append(parse_bcftools_stats(bcf.read_text(encoding="utf-8"), dataset=d.id))
        else:
            missing.append(str(bcf.relative_to(out_dir)))
        ext = Path(f"{happy_prefix(out_dir, d.id)}.extended.csv")
        if ext.is_file():
            calling += parse_happy_extended(
                ext.read_text(encoding="utf-8"), dataset=d.id, pass_only=cfg.happy.pass_only
            )
        else:
            missing.append(str(ext.relative_to(out_dir)))
    return FidelityResult(
        schema_version=SCHEMA_VERSION,
        name=cfg.name,
        baseline=cfg.baseline,
        datasets={d.id: d.label for d in cfg.datasets},
        pass_only=cfg.happy.pass_only,
        recall_gap=cfg.recall_gap,
        planting=planting,
        read_qc=read_qc,
        callset_qc=callset_qc,
        calling=calling,
        gaps=recall_gaps(calling, cfg.baseline, cfg.recall_gap),
        missing=missing,
    )


def _fmt(value: float | int | None) -> str:
    return "" if value is None else str(value)


def write_result(result: FidelityResult, out_dir: Path) -> dict[str, Path]:
    paths = {"fidelity": write_json(out_dir / "fidelity.json", result)}
    calling = out_dir / "fidelity_calling.tsv"
    with calling.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            ["dataset", "type", "stratum", "subtype", "truth_tp", "truth_tp_fn", "recall",
             "recall_lower", "recall_upper", "query_tp", "query_tp_fp", "precision",
             "precision_lower", "precision_upper", "f1"]
        )  # fmt: skip
        for m in result.calling:
            r, p = m.recall, m.precision
            w.writerow(
                [m.dataset, m.type, m.stratum, m.subtype, r.k, r.n, _fmt(r.value), _fmt(r.lower),
                 _fmt(r.upper), p.k, p.n, _fmt(p.value), _fmt(p.lower), _fmt(p.upper), _fmt(m.f1)]
            )  # fmt: skip
    paths["fidelity_calling"] = calling
    qc = out_dir / "fidelity_qc.tsv"
    with qc.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["dataset", "metric", "value"])
        items: list[ReadQc | CallsetQc] = [*result.read_qc, *result.callset_qc]
        for item in items:
            for key, value in item.model_dump().items():
                if key != "dataset":
                    w.writerow([item.dataset, key, _fmt(value)])
    paths["fidelity_qc"] = qc
    return paths
