"""Deterministic JSON serialisation of contracts."""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from pydantic import BaseModel


def dumps(data: BaseModel | Sequence[BaseModel] | dict[str, Any]) -> str:
    if isinstance(data, BaseModel):
        payload: Any = data.model_dump(mode="json")
    elif isinstance(data, (list, tuple)):
        payload = [item.model_dump(mode="json") for item in data]
    else:
        payload = data
    return json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def write_json(path: Path, data: BaseModel | Sequence[BaseModel] | dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dumps(data), encoding="utf-8")
    return path
