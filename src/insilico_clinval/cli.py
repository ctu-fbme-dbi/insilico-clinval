"""Command-line interface (entry point ``clinval``); one sub-command per step.

Exit codes: 0 success; 1 invalid suite, missing input or tool failure; 2 selection
unsatisfied (required rule or empty sample); 3 truth encoding failed (e.g. REF mismatch, including
reference bases stated in eRepo HGVS); 4 simulation failed (a planted variant cannot be simulated or
is missing from the simulated reads); 5 the system under test failed (or broke its output contract).
"""

from __future__ import annotations

import json
import logging
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import typer

from insilico_clinval import __version__
from insilico_clinval import fidelity as fidelity_config
from insilico_clinval._exec import ToolError
from insilico_clinval._runner import RuntimeConfigError, Toolbox, load_runtime
from insilico_clinval.adapters import (
    LoadedPipeline,
    PipelineConfigError,
    SampleInput,
    SutFailure,
    load_pipeline,
    missing_files,
)
from insilico_clinval.adapters.run import load_secondary, run_secondary, run_tertiary
from insilico_clinval.compare.run import (
    CompareError,
    CompareParams,
    compare_dir,
    compare_sample,
    write_outcomes_tsv,
)
from insilico_clinval.config import (
    LoadedSuite,
    SuiteError,
    check_inputs_exist,
    load_suite,
    required_inputs,
)
from insilico_clinval.diff.run import diff as diff_runs
from insilico_clinval.diff.run import write_diff_tsv
from insilico_clinval.fidelity.planting import PlantingError
from insilico_clinval.fidelity.run import STEPS as FIDELITY_STEPS
from insilico_clinval.fidelity.run import assemble, run_calling, run_planting, run_qc, write_result
from insilico_clinval.fidelity.stats import StatsParseError
from insilico_clinval.genome import refdiff
from insilico_clinval.genome.reference import Reference
from insilico_clinval.jsonio import write_json
from insilico_clinval.metrics.acceptance import AcceptanceError, evaluate, load_acceptance
from insilico_clinval.metrics.compute import compute
from insilico_clinval.metrics.load import MetricsInputError, load_run, read_par
from insilico_clinval.provenance import ManifestBuilder, ManifestError, sha256_file
from insilico_clinval.report.markdown import diff_markdown, metrics_markdown, write_metrics_tsv
from insilico_clinval.schemas.fidelity import FidelityConfig, FidelityResult
from insilico_clinval.schemas.finding import FindingSet, SecondaryResult
from insilico_clinval.schemas.manifest import RunManifest
from insilico_clinval.schemas.metrics import AcceptanceConfig, DiffReport, MetricsReport
from insilico_clinval.schemas.outcome import ExpectedOutcome, ExpectedOutcomeSet, OutcomeSet
from insilico_clinval.schemas.pipeline import PipelineConfig
from insilico_clinval.schemas.runtime import RuntimeConfig
from insilico_clinval.schemas.sample import SampleSpec
from insilico_clinval.schemas.selection import SelectionResult
from insilico_clinval.schemas.simulation import PlantingCheck, SimulationProfile, SimulationSpec
from insilico_clinval.schemas.suite import SimulationConfig, SuiteConfig
from insilico_clinval.schemas.variant import TestVariant
from insilico_clinval.simprep import SimPrepError
from insilico_clinval.simprep.prepare import SimPrep, prepare_sample
from insilico_clinval.simulate import SimulationError, load_profile, load_truth_outputs
from insilico_clinval.simulate.run import recheck_planting, sim_dir_name, simulate_sample
from insilico_clinval.truth import TruthError
from insilico_clinval.truth.build import build_truth
from insilico_clinval.truth.genotype import par_regions
from insilico_clinval.variants import held_out as held_out_mod
from insilico_clinval.variants.erepo import (
    ErepoError,
    ErepoMismatchError,
    approval_dates,
    load_held_out_reference,
)
from insilico_clinval.variants.report import SELECTED_JSON, write_selection_outputs
from insilico_clinval.variants.select import SelectionError, select_variants, selection_failures

app = typer.Typer(
    no_args_is_help=True,
    help="End-to-end in silico validation of clinical genome sequencing pipelines.",
)

EXIT_INPUT = 1
EXIT_SELECTION = 2
EXIT_TRUTH = 3
EXIT_SIMULATION = 4
EXIT_SUT = 5
EXIT_COMPARE = 6
EXIT_REJECTED = 7
LOG_NAME = "clinval.log"

SCHEMA_MODELS: dict[str, type] = {
    "SuiteConfig": SuiteConfig,
    "TestVariant": TestVariant,
    "SelectionResult": SelectionResult,
    "ExpectedOutcome": ExpectedOutcome,
    "ExpectedOutcomeSet": ExpectedOutcomeSet,
    "OutcomeSet": OutcomeSet,
    "AcceptanceConfig": AcceptanceConfig,
    "MetricsReport": MetricsReport,
    "DiffReport": DiffReport,
    "SampleSpec": SampleSpec,
    "RunManifest": RunManifest,
    "RuntimeConfig": RuntimeConfig,
    "SimulationProfile": SimulationProfile,
    "SimulationSpec": SimulationSpec,
    "PlantingCheck": PlantingCheck,
    "PipelineConfig": PipelineConfig,
    "SecondaryResult": SecondaryResult,
    "FindingSet": FindingSet,
    "FidelityConfig": FidelityConfig,
    "FidelityResult": FidelityResult,
}

SuiteOption = typer.Option(..., "--suite", help="Suite YAML file.", dir_okay=False)
OutOption = typer.Option(None, "--out", help="Run directory (default: output.directory).")
SchemaOutOption = typer.Option(Path("docs/json_schema"), "--out", help="Target directory.")
SuiteArgument = typer.Argument(..., help="Suite YAML file.", dir_okay=False)
RuntimeOption = typer.Option(
    None, "--runtime", help="Runtime YAML (runner per tool); default: all local.", dir_okay=False
)
SampleOption = typer.Option(None, "--sample", help="Only these samples (repeatable).")
KeepWorkOption = typer.Option(False, "--keep-work", help="Keep Mason's raw outputs.")
PipelineOption = typer.Option(
    ..., "--pipeline", help="Pipeline config (system under test).", dir_okay=False
)
PipelineArgument = typer.Argument(..., help="Pipeline config YAML.", dir_okay=False)
StageOption = typer.Option("all", "--stage", help="secondary, tertiary or all.")
SarekResultsOption = typer.Option(
    None,
    "--sarek-results",
    help="Adopt a finished sarek run (its outdir; params.yaml and samplesheet.csv in the parent) "
    "instead of running sarek; it must match the pipeline config.",
)
PipelineReferenceOption = typer.Option(
    None,
    "--pipeline-reference",
    help="FASTA the pipeline aligned to (decodes CRAM, used by hap.py); default: the suite's.",
    dir_okay=False,
)
ConfidentBedOption = typer.Option(
    None,
    "--confident-bed",
    help="BED reported as stratum giab_confident (default: the suite's background BED).",
    dir_okay=False,
)
WindowPadOption = typer.Option(50, "--window-pad", min=0, help="Padding of planted windows (bp).")
MismappedOption = typer.Option(
    0.5, "--mismapped-fraction", min=0.0, max=1.0,
    help="mismapped if the alignment has fewer ALT reads than this fraction of the true ones.",
)  # fmt: skip
EdgeOption = typer.Option(500, "--edge-bp", min=0, help="Trim of simulated segment ends (bp).")
ThreadsOption = typer.Option(4, "--threads", min=1, help="hap.py threads.")
SimReferenceOption = typer.Option(
    ...,
    "--reference",
    help="FASTA (with .fai) to simulate from, e.g. GIAB-masked GRCh38.",
    exists=True,
    dir_okay=False,
)
AgainstReferenceOption = typer.Option(
    ..., "--against", help="FASTA (with .fai) the pipeline aligns to.", exists=True, dir_okay=False
)
MaskedBedOption = typer.Option(
    None, "--masked-bed", help="Write the bases that are N only in --reference (BED).",
    dir_okay=False,
)  # fmt: skip
HappyGenderOption = typer.Option(
    False,
    "--happy-gender/--no-happy-gender",
    help="Pass the sample's sex to hap.py (--gender); needed for haploid chrX/chrY truth.",
)
HeldOutReferenceOption = typer.Option(
    ...,
    "--reference",
    help="Held-out reference: GRCh38 VCF of its alleles, or an eRepo export.",
    dir_okay=False,
)
HeldOutReleaseOption = typer.Option(..., "--release", help="Release of the reference (v2.2.0).")
ApprovedAfterOption = typer.Option(
    None,
    "--approved-after",
    help="Snapshot date of the reference (YYYY-MM-DD) for the approval-date check.",
)
RunsOption = typer.Option(..., "--run", help="Run directory (repeatable).", file_okay=False)
MetricsOutOption = typer.Option(..., "--out", help="Output directory.", file_okay=False)
AcceptanceOption = typer.Option(
    None, "--acceptance", help="Acceptance criteria YAML (configs/acceptance/).", dir_okay=False
)
ParBedOption = typer.Option(
    None, "--par-bed", help="PAR BED; splits chrX/chrY into PAR and non-PAR groups.", dir_okay=False
)
FailOnRejectOption = typer.Option(
    False, "--fail-on-reject", help="Exit with code 7 when a criterion fails."
)
BaselineOption = typer.Option(..., "--baseline", help="Baseline pipeline config.", dir_okay=False)
CandidateOption = typer.Option(
    ..., "--candidate", help="Candidate pipeline config.", dir_okay=False
)
NoiseOption = typer.Option(
    None, "--noise", help="Second run of the baseline (labels noise).", dir_okay=False
)
FailOnRegressionOption = typer.Option(
    False, "--fail-on-regression", help="Exit with code 7 on a regression."
)
MetricsDirOption = typer.Option(
    ..., "--metrics", help="Output of clinval metrics.", file_okay=False
)
DiffDirOption = typer.Option(None, "--diff", help="Output of clinval diff.", file_okay=False)
ReportOutOption = typer.Option(..., "--out", help="Markdown file.", dir_okay=False)
FidelityConfigOption = typer.Option(
    ..., "--config", help="Fidelity config YAML (configs/fidelity/).", dir_okay=False
)
FidelityOutOption = typer.Option(..., "--out", help="Output directory.", file_okay=False)
FidelityStepOption = typer.Option(
    None, "--step", help="planting, qc or calling (repeatable; default: all)."
)


def _fail(message: str, code: int) -> typer.Exit:
    typer.echo(f"error: {message}", err=True)
    return typer.Exit(code=code)


# Loggers whose messages (phases and progress of the long-running steps) are also shown on stderr.
PROGRESS_LOGGERS = (
    "insilico_clinval.variants.select.progress",
    "insilico_clinval.truth",
    "insilico_clinval.simprep",
    "insilico_clinval.simulate",
)


class _LoggerFilter(logging.Filter):
    def __init__(self, names: tuple[str, ...]):
        super().__init__()
        self._filters = [logging.Filter(name) for name in names]

    def filter(self, record: logging.LogRecord) -> bool:
        return any(f.filter(record) for f in self._filters)


@contextmanager
def _run_logging(out_dir: Path, *, console: tuple[str, ...] = ()) -> Iterator[None]:
    """Log to ``<out_dir>/clinval.log``; records of the ``console`` loggers also go to stderr."""
    out_dir.mkdir(parents=True, exist_ok=True)
    handlers: list[logging.Handler] = []
    file_handler = logging.FileHandler(out_dir / LOG_NAME, encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    handlers.append(file_handler)
    if console:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(logging.Formatter("%(asctime)s %(message)s", "%H:%M:%S"))
        stream_handler.addFilter(_LoggerFilter(console))
        handlers.append(stream_handler)
    logger = logging.getLogger("insilico_clinval")
    logger.setLevel(logging.INFO)
    for handler in handlers:
        logger.addHandler(handler)
    try:
        yield
    finally:
        for handler in handlers:
            logger.removeHandler(handler)
            handler.close()


def _prepare(suite_path: Path, out: Path | None) -> tuple[LoadedSuite, Path]:
    try:
        suite = load_suite(suite_path)
        check_inputs_exist(suite)
    except SuiteError as exc:
        raise _fail(str(exc), EXIT_INPUT) from exc
    return suite, suite.output_dir(out)


def _manifest(suite: LoadedSuite, out_dir: Path) -> ManifestBuilder:
    try:
        manifest = ManifestBuilder(
            out_dir,
            config_path=suite.path,
            config_sha256=suite.sha256,
            seed=suite.config.suite.seed,
        )
    except ManifestError as exc:
        raise _fail(str(exc), EXIT_INPUT) from exc
    for role, path in required_inputs(suite):
        manifest.add_input(role, path)
    sources = suite.config.sources
    hashes = {ref.role: ref.sha256 for ref in manifest.manifest.inputs}
    erepo = sources.clingen_erepo
    if erepo is not None and "clingen_erepo" in hashes:
        manifest.add_database("ClinGen eRepo", erepo.release, hashes["clingen_erepo"])
        if erepo.held_out_reference is not None:
            manifest.add_database(
                "ClinGen eRepo (held-out reference)",
                erepo.held_out_reference.release,
                hashes["clingen_erepo_held_out"],
            )
    bg = suite.config.background
    if bg is not None and "background_vcf" in hashes:
        manifest.add_database(bg.name, bg.release, hashes["background_vcf"])
    if sources.clinvar is not None and "clinvar_vcf" in hashes:
        manifest.add_database("ClinVar", sources.clinvar.release, hashes["clinvar_vcf"])
    if sources.giab is not None and "giab_truth_vcf" in hashes:
        manifest.add_database(sources.giab.name, sources.giab.release, hashes["giab_truth_vcf"])
    return manifest


def _toolbox(runtime: Path | None, manifest: ManifestBuilder | None) -> Toolbox:
    try:
        config = load_runtime(runtime)
    except RuntimeConfigError as exc:
        raise _fail(str(exc), EXIT_INPUT) from exc
    if runtime is not None and manifest is not None:
        manifest.add_input("runtime_config", runtime)
    return Toolbox(config)


def _select(
    suite: LoadedSuite, out_dir: Path, manifest: ManifestBuilder, reference: Reference
) -> SelectionResult:
    try:
        with manifest.step("select"):
            result = select_variants(suite, reference)
    except SelectionError as exc:
        manifest.write()
        raise _fail(str(exc), EXIT_INPUT) from exc
    except ErepoMismatchError as exc:
        manifest.write()
        raise _fail(str(exc), EXIT_TRUTH) from exc
    for role, path in write_selection_outputs(suite, result, out_dir).items():
        manifest.add_output(role, path)
    for w in result.warnings:
        manifest.warn(w)
    failures = selection_failures(result, [s.sample_id for s in suite.config.samples])
    if failures:
        manifest.write()
        raise _fail(
            "selection failed (see selection_report.md):\n  " + "\n  ".join(failures),
            EXIT_SELECTION,
        )
    return result


@app.callback()
def main() -> None:
    """insilico-clinval command-line interface."""


@app.command()
def version() -> None:
    """Print the framework version."""
    typer.echo(__version__)


@app.command()
def select(suite: Path = SuiteOption, out: Path | None = OutOption) -> None:
    """Step 1: select test variants and write variants.selected.{json,tsv} + selection_report.md."""
    loaded, out_dir = _prepare(suite, out)
    with _run_logging(out_dir, console=PROGRESS_LOGGERS):
        manifest = _manifest(loaded, out_dir)
        reference = Reference(loaded.resolve(loaded.config.reference.fasta))
        try:
            result = _select(loaded, out_dir, manifest, reference)
        finally:
            reference.close()
        manifest.write()
    typer.echo(f"selected {len(result.variants)} variants -> {out_dir}")


@app.command("build-truth")
def build_truth_cmd(
    suite: Path = SuiteOption, out: Path | None = OutOption, runtime: Path | None = RuntimeOption
) -> None:
    """Steps 1+2: select variants and write per-sample truth VCF, expected outcomes, sample spec."""
    loaded, out_dir = _prepare(suite, out)
    with _run_logging(out_dir, console=PROGRESS_LOGGERS):
        manifest = _manifest(loaded, out_dir)
        toolbox = _toolbox(runtime, manifest)
        reference = Reference(loaded.resolve(loaded.config.reference.fasta))
        try:
            result = _select(loaded, out_dir, manifest, reference)
            try:
                with manifest.step("build-truth"):
                    outputs = build_truth(
                        loaded, result, reference, out_dir, manifest, toolbox=toolbox
                    )
            except TruthError as exc:
                manifest.write()
                raise _fail(str(exc), EXIT_TRUTH) from exc
            except ToolError as exc:
                manifest.write()
                raise _fail(str(exc), EXIT_INPUT) from exc
        finally:
            reference.close()
        for o in outputs:
            manifest.add_output(f"truth_vcf:{o.sample_id}", o.truth_vcf)
            manifest.add_output(f"truth_vcf_index:{o.sample_id}", o.truth_tbi)
            manifest.add_output(f"expected_outcomes:{o.sample_id}", o.expected_outcomes)
            manifest.add_output(f"sample_spec:{o.sample_id}", o.sample_spec)
            if o.background_removed is not None:
                manifest.add_output(f"background_removed:{o.sample_id}", o.background_removed)
        manifest.write()
    typer.echo(f"truth for {len(outputs)} samples -> {out_dir}")


def _simulation_inputs(
    suite: LoadedSuite, manifest: ManifestBuilder
) -> tuple[SimulationConfig, SimulationProfile, str, Path | None]:
    sim = suite.config.simulation
    if sim is None:
        raise _fail("the suite has no `simulation` section", EXIT_INPUT)
    profile_path = suite.resolve(sim.profile)
    try:
        profile, profile_sha = load_profile(profile_path)
    except SimulationError as exc:
        raise _fail(str(exc), EXIT_INPUT) from exc
    manifest.add_input("simulation_profile", profile_path)
    extra = suite.resolve(sim.extra_regions_bed) if sim.extra_regions_bed else None
    if extra is not None:
        if not extra.is_file():
            raise _fail(f"simulation.extra_regions_bed not found: {extra}", EXIT_INPUT)
        manifest.add_input("simulation_extra_regions", extra)
    return sim, profile, profile_sha, extra


def _sample_ids(suite: LoadedSuite, requested: list[str] | None) -> list[str]:
    known = [s.sample_id for s in suite.config.samples]
    if not requested:
        return known
    unknown = [s for s in requested if s not in known]
    if unknown:
        raise _fail(f"unknown samples: {', '.join(unknown)}", EXIT_INPUT)
    return [s for s in known if s in requested]


def _prepare_samples(
    suite: LoadedSuite,
    out_dir: Path,
    manifest: ManifestBuilder,
    reference: Reference,
    sample_ids: list[str],
) -> list[tuple[SimPrep, SampleSpec]]:
    sim, profile, _, extra = _simulation_inputs(suite, manifest)
    pars = par_regions(suite.config)
    sexes = {s.sample_id: s.sex for s in suite.config.samples}
    preps: list[tuple[SimPrep, SampleSpec]] = []
    for sid in sample_ids:
        try:
            spec, truth = load_truth_outputs(out_dir, sid, manifest)
            with manifest.step(f"prepare-sim ({sid})"):
                prep = prepare_sample(
                    sample_id=sid,
                    sex=sexes[sid],
                    truth_vcf=truth,
                    planted_ids=[v.variant_id for v in spec.variants],
                    reference=reference,
                    pars=pars,
                    sim=sim,
                    extra_bed=extra,
                    min_segment_bp=profile.fragment.longest,
                    sim_dir=out_dir / sim_dir_name(sid),
                )
        except SimulationError as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_INPUT) from exc
        except (SimPrepError, ValueError) as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_SIMULATION) from exc
        for group in prep.groups:
            if group.vcf is not None:
                manifest.add_output(f"mason_vcf:{sid}:{group.group}", group.vcf)
        manifest.add_output(f"segments:{sid}", prep.segments_tsv)
        if prep.dropped:
            manifest.warn(
                f"sample {sid}: {len(prep.dropped)} background records cannot be simulated "
                f"and were dropped (see {prep.dropped_log})"
            )
        if prep.decomposed:
            manifest.warn(
                f"sample {sid}: {len(prep.decomposed)} planted MNV/delins split into Mason records"
            )
        preps.append((prep, spec))
    return preps


@app.command("prepare-sim")
def prepare_sim(
    suite: Path = SuiteOption,
    out: Path | None = OutOption,
    sample: list[str] | None = SampleOption,
) -> None:
    """Step 3a: write Mason inputs (segments, segment FASTA, Mason VCF) from the truth outputs."""
    loaded, out_dir = _prepare(suite, out)
    with _run_logging(out_dir, console=PROGRESS_LOGGERS):
        manifest = _manifest(loaded, out_dir)
        reference = Reference(loaded.resolve(loaded.config.reference.fasta))
        try:
            preps = _prepare_samples(
                loaded, out_dir, manifest, reference, _sample_ids(loaded, sample)
            )
        finally:
            reference.close()
        manifest.write()
    typer.echo(f"Mason inputs for {len(preps)} samples -> {out_dir}")


@app.command()
def simulate(
    suite: Path = SuiteOption,
    out: Path | None = OutOption,
    runtime: Path | None = RuntimeOption,
    sample: list[str] | None = SampleOption,
    keep_work: bool = KeepWorkOption,
) -> None:
    """Steps 3a+3b: Mason inputs, mason_simulator, true alignments and planting check."""
    loaded, out_dir = _prepare(suite, out)
    incomplete: list[str] = []
    with _run_logging(out_dir, console=PROGRESS_LOGGERS):
        manifest = _manifest(loaded, out_dir)
        toolbox = _toolbox(runtime, manifest)
        sim, profile, profile_sha, _ = _simulation_inputs(loaded, manifest)
        reference = Reference(loaded.resolve(loaded.config.reference.fasta))
        try:
            preps = _prepare_samples(
                loaded, out_dir, manifest, reference, _sample_ids(loaded, sample)
            )
            for prep, spec in preps:
                try:
                    with manifest.step(f"simulate ({prep.sample_id})"):
                        result = simulate_sample(
                            prep=prep,
                            spec=spec,
                            sim=sim,
                            profile=profile,
                            profile_sha256=profile_sha,
                            seed=loaded.config.suite.seed,
                            reference=reference,
                            out_dir=out_dir,
                            toolbox=toolbox,
                            manifest=manifest,
                            keep_work=keep_work,
                        )
                except ToolError as exc:
                    manifest.write()
                    raise _fail(str(exc), EXIT_INPUT) from exc
                sid = result.sample_id
                manifest.add_output(f"fastq_r1:{sid}", result.fastq_r1)
                manifest.add_output(f"fastq_r2:{sid}", result.fastq_r2)
                manifest.add_output(f"true_bam:{sid}", result.true_bam)
                manifest.add_output(f"true_bam_index:{sid}", result.true_bai)
                manifest.add_output(f"simulation_spec:{sid}", result.spec)
                manifest.add_output(f"planting_check:{sid}", result.check_json)
                manifest.add_output(f"planting_check_tsv:{sid}", result.check_tsv)
                if not result.complete:
                    incomplete.append(sid)
        finally:
            reference.close()
        if incomplete:
            manifest.warn(f"planting check incomplete for samples: {', '.join(incomplete)}")
        manifest.write()
    if incomplete:
        raise _fail(
            "planted variants missing from the simulated reads (see <sample>.planting_check.tsv): "
            + ", ".join(incomplete),
            EXIT_SIMULATION,
        )
    typer.echo(f"simulated {len(preps)} samples -> {out_dir}")


def _load_pipeline(path: Path, *, check_files: bool) -> LoadedPipeline:
    try:
        pipeline = load_pipeline(path)
    except PipelineConfigError as exc:
        raise _fail(str(exc), EXIT_INPUT) from exc
    if check_files:
        missing = missing_files(pipeline)
        if missing:
            raise _fail(
                f"{path}: missing files referenced by the pipeline:\n  " + "\n  ".join(missing),
                EXIT_INPUT,
            )
    return pipeline


def _sample_inputs(
    suite: LoadedSuite, out_dir: Path, manifest: ManifestBuilder, sample_ids: list[str]
) -> list[SampleInput]:
    """Simulated FASTQ of each sample; hashes must match the outputs recorded by `simulate`."""
    recorded = {o.role: o.sha256 for o in manifest.manifest.outputs}
    sexes = {s.sample_id: s.sex for s in suite.config.samples}
    inputs: list[SampleInput] = []
    for sid in sample_ids:
        paths = []
        for mate in ("r1", "r2"):
            path = out_dir / f"{sid}.{mate.upper()}.fastq.gz"
            sha = recorded.get(f"fastq_{mate}:{sid}")
            if not path.is_file() or sha is None:
                raise _fail(
                    f"{path.name} not found or not recorded (run `clinval simulate`)", EXIT_INPUT
                )
            if sha256_file(path) != sha:
                raise _fail(
                    f"{path.name} differs from the file recorded in the manifest", EXIT_INPUT
                )
            paths.append(path)
        inputs.append(SampleInput(sid, sexes[sid], paths[0], paths[1]))
    return inputs


@app.command("recheck-planting")
def recheck_planting_cmd(
    suite: Path = SuiteOption,
    out: Path | None = OutOption,
    sample: list[str] | None = SampleOption,
) -> None:
    """Repeat the planting check of finished simulations from their true BAMs (no Mason run)."""
    loaded, out_dir = _prepare(suite, out)
    sim = loaded.config.simulation
    if sim is None:
        raise _fail("the suite has no simulation section", EXIT_INPUT)
    incomplete: list[str] = []
    sample_ids = _sample_ids(loaded, sample)
    with _run_logging(out_dir, console=PROGRESS_LOGGERS):
        manifest = _manifest(loaded, out_dir)
        reference = Reference(loaded.resolve(loaded.config.reference.fasta))
        try:
            for sid in sample_ids:
                with manifest.step(f"recheck-planting ({sid})"):
                    spec, truth = load_truth_outputs(out_dir, sid, manifest)
                    result = recheck_planting(
                        out_dir, spec, truth, reference, min_alt_reads=sim.min_alt_reads
                    )
                manifest.add_output(f"planting_check:{sid}", result.check_json)
                manifest.add_output(f"planting_check_tsv:{sid}", result.check_tsv)
                manifest.add_output(f"simulation_spec:{sid}", result.spec)
                if not result.complete:
                    incomplete.append(sid)
        except (SimulationError, SimPrepError) as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_INPUT) from exc
        finally:
            reference.close()
        if incomplete:
            manifest.warn(f"planting check incomplete for samples: {', '.join(incomplete)}")
        manifest.write()
    if incomplete:
        raise _fail(
            "planted variants missing from the simulated reads (see <sample>.planting_check.tsv): "
            + ", ".join(incomplete),
            EXIT_SIMULATION,
        )
    typer.echo(f"planting check repeated for {len(sample_ids)} samples -> {out_dir}")


@app.command("validate-pipeline")
def validate_pipeline(pipeline: Path = PipelineArgument) -> None:
    """Validate a pipeline config (system under test) and the files it references."""
    loaded = _load_pipeline(pipeline, check_files=True)
    cfg = loaded.config
    typer.echo(
        f"OK: {cfg.name} — secondary {cfg.secondary.adapter}, tertiary {cfg.tertiary.adapter}"
    )


@app.command("run-sut")
def run_sut(
    suite: Path = SuiteOption,
    pipeline: Path = PipelineOption,
    out: Path | None = OutOption,
    runtime: Path | None = RuntimeOption,
    sample: list[str] | None = SampleOption,
    stage: str = StageOption,
    sarek_results: Path | None = SarekResultsOption,
) -> None:
    """Steps 4+5: run the system under test (secondary and/or tertiary) on the simulated reads."""
    if stage not in ("secondary", "tertiary", "all"):
        raise _fail(f"--stage must be secondary, tertiary or all (got {stage})", EXIT_INPUT)
    loaded, out_dir = _prepare(suite, out)
    pipe = _load_pipeline(pipeline, check_files=True)
    name = pipe.config.name
    with _run_logging(out_dir):
        manifest = _manifest(loaded, out_dir)
        toolbox = _toolbox(runtime, manifest)
        manifest.add_input(f"pipeline:{name}", pipe.path)
        sample_ids = _sample_ids(loaded, sample)
        fasta = loaded.resolve(loaded.config.reference.fasta).resolve()
        try:
            if stage in ("secondary", "all"):
                with manifest.step(f"run-sut secondary ({name})"):
                    secondary = run_secondary(
                        pipe,
                        _sample_inputs(loaded, out_dir, manifest, sample_ids),
                        run_dir=out_dir,
                        fasta=fasta,
                        toolbox=toolbox,
                        manifest=manifest,
                        sarek_results=sarek_results,
                    )
                for result, path in secondary:
                    manifest.add_output(f"secondary:{name}:{result.sample_id}", path)
                sec_results = [r for r, _ in secondary]
            else:
                sec_results = [load_secondary(out_dir, pipe, sid) for sid in sample_ids]
            if stage in ("tertiary", "all"):
                with manifest.step(f"run-sut tertiary ({name})"):
                    tertiary = run_tertiary(
                        pipe,
                        sec_results,
                        run_dir=out_dir,
                        fasta=fasta,
                        toolbox=toolbox,
                        manifest=manifest,
                    )
                for finding_set, path in tertiary:
                    manifest.add_output(f"findings:{name}:{finding_set.sample_id}", path)
        except (SutFailure, ToolError) as exc:
            manifest.warn(f"system under test {name} failed: {str(exc).splitlines()[0]}")
            manifest.write()
            raise _fail(str(exc), EXIT_SUT) from exc
        except PipelineConfigError as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_INPUT) from exc
        manifest.write()
    typer.echo(f"{name}: {len(sample_ids)} samples ({stage}) -> {out_dir / 'sut' / name}")


@app.command("held-out")
def held_out(
    suite: Path = SuiteOption,
    out: Path | None = OutOption,
    reference_file: Path = HeldOutReferenceOption,
    release: str = HeldOutReleaseOption,
    approved_after: str | None = ApprovedAfterOption,
) -> None:
    """Recompute held_out of an existing run's eRepo variants (writes held_out.tsv)."""
    loaded, out_dir = _prepare(suite, out)
    selected = out_dir / SELECTED_JSON
    if not selected.is_file():
        raise _fail(f"{selected} not found (run clinval select first)", EXIT_INPUT)
    if not reference_file.is_file():
        raise _fail(f"held-out reference not found: {reference_file}", EXIT_INPUT)
    with _run_logging(out_dir):
        manifest = _manifest(loaded, out_dir)
        manifest.add_input("held_out_reference", reference_file)
        manifest.add_database(
            "ClinGen eRepo (held-out reference)", release, sha256_file(reference_file)
        )
        try:
            with manifest.step("held-out"):
                reference = Reference(loaded.resolve(loaded.config.reference.fasta))
                ref_set = load_held_out_reference(reference_file, release, reference)
                erepo = loaded.config.sources.clingen_erepo
                approvals = approval_dates(loaded.resolve(erepo.file)) if erepo else None
                variants = SelectionResult.model_validate_json(selected.read_text()).variants
                rows = held_out_mod.annotate(
                    variants, ref_set, reference, approvals=approvals, cutoff=approved_after
                )
        except ErepoError as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_INPUT) from exc
        path = held_out_mod.write_tsv(rows, out_dir / held_out_mod.HELD_OUT_TSV)
        manifest.add_output("held_out", path)
        n_out = sum(r.held_out for r in rows)
        checked = [r for r in rows if r.agrees_with_date is not None]
        disagree = [r.variant_id for r in checked if not r.agrees_with_date]
        if disagree:
            manifest.warn(
                f"held-out: {len(disagree)} of {len(checked)} eRepo variants disagree with the "
                f"approval-date check (cutoff {approved_after}): {', '.join(disagree)}"
            )
        manifest.write()
    typer.echo(
        f"{n_out} of {len(rows)} eRepo variants held out (absent from {release}); approval-date "
        f"check: {len(checked) - len(disagree)} of {len(checked)} agree -> {path}"
    )


@app.command("compare")
def compare(
    suite: Path = SuiteOption,
    pipeline: Path = PipelineOption,
    out: Path | None = OutOption,
    runtime: Path | None = RuntimeOption,
    sample: list[str] | None = SampleOption,
    pipeline_reference: Path | None = PipelineReferenceOption,
    confident_bed: Path | None = ConfidentBedOption,
    window_pad: int = WindowPadOption,
    mismapped_fraction: float = MismappedOption,
    edge_bp: int = EdgeOption,
    threads: int = ThreadsOption,
    happy_gender: bool = HappyGenderOption,
) -> None:
    """Step 6: compare the pipeline's outputs with the truth and attribute each failure."""
    loaded, out_dir = _prepare(suite, out)
    pipe = _load_pipeline(pipeline, check_files=False)
    name = pipe.config.name
    cfg = loaded.config
    fasta = pipeline_reference or loaded.resolve(cfg.reference.fasta)
    if confident_bed is None and cfg.background is not None:
        confident_bed = loaded.resolve(cfg.background.bed)
    strata = [(k, loaded.resolve(v)) for k, v in sorted(cfg.sources.strata.items())]
    if confident_bed is not None:
        strata.insert(0, ("giab_confident", confident_bed))
    for role, path in strata:
        if not path.is_file():
            raise _fail(f"stratum {role}: {path} not found", EXIT_INPUT)
    params = CompareParams(
        window_pad=window_pad, mismapped_fraction=mismapped_fraction, edge_bp=edge_bp,
        threads=threads, happy_gender=happy_gender,
    )  # fmt: skip
    sexes = {s.sample_id: s.sex for s in cfg.samples}
    with _run_logging(out_dir, console=("insilico_clinval.compare",)):
        manifest = _manifest(loaded, out_dir)
        toolbox = _toolbox(runtime, manifest)
        manifest.add_input(f"pipeline:{name}", pipe.path)
        for role, path in strata:
            manifest.add_input(f"compare_stratum:{role}", path)
        if pipeline_reference is not None:
            manifest.add_input("pipeline_reference", pipeline_reference)
        try:
            reference = Reference(fasta)
            for sid in _sample_ids(loaded, sample):
                with manifest.step(f"compare ({name}, {sid})"):
                    compare_sample(
                        out_dir, pipe, sid, reference=reference,
                        strata=strata, toolbox=toolbox, manifest=manifest, params=params,
                        sex=sexes.get(sid),
                    )  # fmt: skip
            manifest.add_output(
                f"outcomes_tsv:{name}", write_outcomes_tsv(compare_dir(out_dir, pipe))
            )
        except (CompareError, SutFailure, FileNotFoundError) as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_INPUT) from exc
        except (ToolError, StatsParseError) as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_COMPARE) from exc
        manifest.write()
    typer.echo(f"{name}: compared -> {compare_dir(out_dir, pipe)}")


def _output_manifest(out: Path, config: Path) -> ManifestBuilder:
    """Manifest of an output directory that is not a run directory (seed 0: nothing random)."""
    try:
        return ManifestBuilder(out, config_path=config, config_sha256=sha256_file(config), seed=0)
    except ManifestError as exc:
        raise _fail(str(exc), EXIT_INPUT) from exc


def _metrics_for(
    pipe_name: str, runs: list[Path], par: Path | None, acceptance: Path | None
) -> MetricsReport:
    report = compute(pipe_name, [load_run(r, pipe_name) for r in runs], read_par(par))
    if acceptance is not None:
        cfg = load_acceptance(acceptance)
        report = report.model_copy(
            update={"acceptance": cfg.name, "criteria": evaluate(cfg, report)}
        )
    return report


@app.command("metrics")
def metrics(
    pipeline: Path = PipelineOption,
    run: list[Path] = RunsOption,
    out: Path = MetricsOutOption,
    acceptance: Path | None = AcceptanceOption,
    par_bed: Path | None = ParBedOption,
    fail_on_reject: bool = FailOnRejectOption,
) -> None:
    """Metrics with 95 % Wilson intervals over one or more runs; optional acceptance criteria."""
    pipe = _load_pipeline(pipeline, check_files=False)
    out.mkdir(parents=True, exist_ok=True)
    with _run_logging(out):
        manifest = _output_manifest(out, pipe.path)
        try:
            with manifest.step(f"metrics ({pipe.config.name})"):
                report = _metrics_for(pipe.config.name, run, par_bed, acceptance)
        except (MetricsInputError, AcceptanceError) as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_INPUT) from exc
        for r in run:
            for f in sorted((r / "sut" / pipe.config.name / "compare").glob("*.outcomes.json")):
                manifest.add_input(f"outcomes:{r.name}:{f.name}", f)
        if acceptance is not None:
            manifest.add_input("acceptance", acceptance)
        manifest.add_output("metrics_json", write_json(out / "metrics.json", report))
        manifest.add_output("metrics_tsv", write_metrics_tsv(report, out / "metrics.tsv"))
        manifest.write()
    failed = [c for c in report.criteria if c.status == "fail"]
    typer.echo(
        f"{pipe.config.name}: metrics -> {out}"
        + (f"; {len(failed)} criteria failed" if report.criteria else "")
    )
    if failed and fail_on_reject:
        raise typer.Exit(code=EXIT_REJECTED)


@app.command("diff")
def diff(
    baseline: Path = BaselineOption,
    candidate: Path = CandidateOption,
    run: list[Path] = RunsOption,
    out: Path = MetricsOutOption,
    noise: Path | None = NoiseOption,
    acceptance: Path | None = AcceptanceOption,
    par_bed: Path | None = ParBedOption,
    fail_on_regression: bool = FailOnRegressionOption,
) -> None:
    """Regression diff of a candidate pipeline against a baseline on the same samples."""
    base = _load_pipeline(baseline, check_files=False)
    cand = _load_pipeline(candidate, check_files=False)
    noise_name = _load_pipeline(noise, check_files=False).config.name if noise else None
    out.mkdir(parents=True, exist_ok=True)
    with _run_logging(out):
        manifest = _output_manifest(out, cand.path)
        manifest.add_input("baseline_pipeline", base.path)
        try:
            with manifest.step(f"diff ({cand.config.name} vs {base.config.name})"):
                report = diff_runs(
                    run, base.config.name, cand.config.name, noise=noise_name,
                    par=read_par(par_bed),
                )  # fmt: skip
                if acceptance is not None:
                    cfg = load_acceptance(acceptance)
                    metrics_report = _metrics_for(cand.config.name, run, par_bed, None)
                    report = report.model_copy(
                        update={"acceptance": cfg.name, "criteria": evaluate(cfg, metrics_report)}
                    )
        except (MetricsInputError, AcceptanceError) as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_INPUT) from exc
        manifest.add_output("diff_json", write_json(out / "diff.json", report))
        manifest.add_output("diff_tsv", write_diff_tsv(report, out / "diff.tsv"))
        md = out / "diff.md"
        md.write_text(diff_markdown(report), encoding="utf-8")
        manifest.add_output("diff_md", md)
        manifest.write()
    typer.echo(f"{cand.config.name} vs {base.config.name}: {report.verdict} -> {out}")
    if report.verdict == "regression" and fail_on_regression:
        raise typer.Exit(code=EXIT_REJECTED)


@app.command("report")
def report(
    metrics_dir: Path = MetricsDirOption,
    diff_dir: Path | None = DiffDirOption,
    out: Path = ReportOutOption,
) -> None:
    """Markdown summary of clinval metrics (and optionally clinval diff) outputs."""
    path = metrics_dir / "metrics.json"
    if not path.is_file():
        raise _fail(f"{path} not found (run clinval metrics first)", EXIT_INPUT)
    text = metrics_markdown(MetricsReport.model_validate_json(path.read_text()))
    if diff_dir is not None:
        dpath = diff_dir / "diff.json"
        if not dpath.is_file():
            raise _fail(f"{dpath} not found (run clinval diff first)", EXIT_INPUT)
        text += "\n" + diff_markdown(DiffReport.model_validate_json(dpath.read_text()))
    out.write_text(text, encoding="utf-8")
    typer.echo(f"report -> {out}")


# Input roles (prefixes) recorded in the manifest per fidelity step: alignments are hashed only
# when the step reads them.
FIDELITY_INPUT_ROLES = {
    "planting": (),  # the planting checks themselves are recorded by the step
    "qc": ("reference", "confident_bed", "qc_regions", "alignment:", "vcf:"),
    "calling": ("reference", "confident_bed", "stratum:", "vcf:", "truth_vcf:", "vcfeval_template"),
}


def _load_fidelity(path: Path) -> fidelity_config.LoadedFidelity:
    try:
        fid = fidelity_config.load_fidelity(path)
    except fidelity_config.FidelityConfigError as exc:
        raise _fail(str(exc), EXIT_INPUT) from exc
    missing = fidelity_config.missing_files(fid)
    if missing:
        raise _fail("missing inputs:\n  " + "\n  ".join(missing), EXIT_INPUT)
    return fid


@app.command()
def fidelity(
    config: Path = FidelityConfigOption,
    out: Path = FidelityOutOption,
    runtime: Path | None = RuntimeOption,
    step: list[str] | None = FidelityStepOption,
) -> None:
    """Fidelity of simulated vs real data: planting summary, read QC, stratified hap.py."""
    steps = step or list(FIDELITY_STEPS)
    unknown = sorted(set(steps) - set(FIDELITY_STEPS))
    if unknown:
        raise _fail(f"unknown --step {unknown} (known: {', '.join(FIDELITY_STEPS)})", EXIT_INPUT)
    fid = _load_fidelity(config)
    with _run_logging(out, console=("insilico_clinval.fidelity",)):
        try:
            manifest = ManifestBuilder(out, config_path=fid.path, config_sha256=fid.sha256, seed=0)
        except ManifestError as exc:
            raise _fail(str(exc), EXIT_INPUT) from exc
        toolbox = _toolbox(runtime, manifest)
        for role, path in fidelity_config.referenced_files(fid):
            prefixes = [r for s in steps for r in FIDELITY_INPUT_ROLES[s]]
            if path.is_file() and role.startswith(tuple(prefixes)):
                manifest.add_input(role, path)
        try:
            if "planting" in steps and fid.config.planting:
                for run_dir in (fid.resolve(p) for p in fid.config.planting):
                    for check in sorted(run_dir.glob("*.planting_check.json")):
                        manifest.add_input(f"planting_check:{run_dir.name}", check)
                with manifest.step("fidelity planting"):
                    manifest.add_output("planting", run_planting(fid, out))
            if "qc" in steps:
                for path in run_qc(fid, out, toolbox, manifest):
                    manifest.add_output(f"fidelity_qc:{path.name}", path)
            if "calling" in steps:
                for path in run_calling(fid, out, toolbox, manifest):
                    manifest.add_output(f"happy:{path.name}", path)
            with manifest.step("fidelity assemble"):
                result = assemble(fid, out)
            for role, path in write_result(result, out).items():
                manifest.add_output(role, path)
        except (ToolError, PlantingError, StatsParseError) as exc:
            manifest.write()
            raise _fail(str(exc), EXIT_INPUT) from exc
        for item in result.missing:
            manifest.warn(f"fidelity: raw output not found: {item}")
        manifest.write()
    typer.echo(f"{fid.config.name}: steps {', '.join(steps)} -> {out / 'fidelity.json'}")
    if result.missing:
        typer.echo("not yet available: " + ", ".join(result.missing), err=True)


@app.command("check-references")
def check_references(
    reference: Path = SimReferenceOption,
    against: Path = AgainstReferenceOption,
    masked_bed: Path | None = MaskedBedOption,
) -> None:
    """Check that two reference builds share the primary contigs (name, order, length)."""
    a, b = Reference(reference), Reference(against)
    try:
        problems = refdiff.primary_problems(a, b)
        if problems:
            raise _fail(f"{reference.name} vs {against.name}: " + "; ".join(problems), EXIT_INPUT)
        typer.echo(f"{reference.name} and {against.name}: primary contigs identical")
        if masked_bed is not None:
            n = refdiff.write_masked_bed(a, b, masked_bed)
            typer.echo(f"{n:,} bp masked only in {reference.name} -> {masked_bed}")
    finally:
        a.close()
        b.close()


@app.command("export-schemas")
def export_schemas(
    out: Path = SchemaOutOption,
) -> None:
    """Write the JSON Schema of every persisted contract."""
    out.mkdir(parents=True, exist_ok=True)
    for name, model in SCHEMA_MODELS.items():
        schema = model.model_json_schema()  # type: ignore[attr-defined]
        text = json.dumps(schema, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
        (out / f"{name}.schema.json").write_text(text, encoding="utf-8")
    typer.echo(f"wrote {len(SCHEMA_MODELS)} schemas -> {out}")


@app.command("validate-suite")
def validate_suite(suite: Path = SuiteArgument) -> None:
    """Validate a suite file: schema, cross-references and presence of all inputs."""
    loaded, _ = _prepare(suite, None)
    cfg = loaded.config
    typer.echo(
        f"OK: {cfg.suite.name} {cfg.suite.version} — {len(cfg.selection)} rules, "
        f"{len(cfg.samples)} samples"
    )


if __name__ == "__main__":
    app()
