"""Common pydantic base class for all persisted contracts."""

from pydantic import BaseModel, ConfigDict


class Contract(BaseModel):
    """Immutable model that rejects unknown keys."""

    model_config = ConfigDict(extra="forbid", frozen=True)
