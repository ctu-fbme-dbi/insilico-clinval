"""RunManifest: provenance of one run directory."""

from pydantic import BaseModel, ConfigDict, Field

from insilico_clinval.schemas.base import Contract


class FrameworkInfo(Contract):
    version: str
    git_commit: str | None
    dirty: bool | None
    python_packages: dict[str, str] = Field(
        default_factory=dict,
        description="versions of Python libraries that shape output bytes (pysam, isal)",
    )


class ConfigRef(Contract):
    path: str
    sha256: str


class FileRef(Contract):
    role: str
    path: str
    sha256: str


class ToolRef(Contract):
    name: str
    version: str
    container_digest: str | None = None


class DatabaseRef(Contract):
    name: str
    release: str
    sha256: str | None


class StepRecord(Contract):
    name: str
    started_at: str
    wall_s: float
    peak_rss_mb: float | None
    exit_code: int


class RunManifest(BaseModel):
    """Mutable while a command runs; extended by every command writing into the run directory."""

    model_config = ConfigDict(extra="forbid")

    schema_version: int
    run_id: str
    created_at: str
    framework: FrameworkInfo
    config: ConfigRef
    seed: int
    inputs: list[FileRef]
    outputs: list[FileRef]
    tools: list[ToolRef]
    databases: list[DatabaseRef]
    steps: list[StepRecord]
    warnings: list[str]
