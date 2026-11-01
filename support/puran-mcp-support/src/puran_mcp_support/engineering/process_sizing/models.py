from __future__ import annotations
import hashlib
import json
from typing import Any, Literal
from pydantic import BaseModel, ConfigDict, Field, model_serializer, model_validator

def content_hash(value: Any) -> str:
    if isinstance(value, BaseModel):
        value = value.model_dump(mode='json')

    def encode(item):
        if isinstance(item, BaseModel):
            return item.model_dump(mode='json')
        raise TypeError(f'not canonical JSON: {type(item).__name__}')
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False, default=encode).encode()).hexdigest()
