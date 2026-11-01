"""RuntimeConfig: how external tools are executed (independent of the suite)."""

from typing import Literal

from pydantic import Field, model_validator

from insilico_clinval.schemas.base import Contract

RunnerName = Literal["local", "docker"]

# Closed set of external tools the framework may call (typos are configuration errors).
KNOWN_TOOLS: tuple[str, ...] = (
    "bcftools",
    "samtools",
    "mason_simulator",
    "hap.py",
    "rtg",
    "nextflow",
    "nirvana",
    "vep",
    "bias_2015",
    "bias_2015_paths",
)


class DockerSettings(Contract):
    binary: str = "docker"
    extra_args: list[str] = Field(
        default_factory=list, description="added to every `docker run`, e.g. --platform"
    )


class ToolRuntime(Contract):
    runner: RunnerName | None = Field(default=None, description="null = default_runner")
    image: str | None = Field(default=None, description="docker image, name@sha256:… preferred")
    executable: str | None = Field(default=None, description="null = the tool name")
    version_args: list[str] | None = Field(
        default=None, description="arguments that print the version; null = the tool's default"
    )


class RuntimeConfig(Contract):
    schema_version: Literal[1] = 1
    default_runner: RunnerName = "local"
    docker: DockerSettings = Field(default_factory=DockerSettings)
    tools: dict[str, ToolRuntime] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _tools(self) -> "RuntimeConfig":
        errors: list[str] = []
        for name, tool in sorted(self.tools.items()):
            if name not in KNOWN_TOOLS:
                errors.append(f"unknown tool {name!r} (known: {', '.join(KNOWN_TOOLS)})")
            if self.runner_of(name) == "docker" and tool.image is None:
                errors.append(f"tool {name}: the docker runner needs an image")
        if errors:
            raise ValueError("; ".join(errors))
        return self

    def runner_of(self, tool: str) -> RunnerName:
        spec = self.tools.get(tool)
        if spec is not None and spec.runner is not None:
            return spec.runner
        return self.default_runner
