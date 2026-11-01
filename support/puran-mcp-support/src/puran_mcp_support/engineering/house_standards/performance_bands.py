from __future__ import annotations
from typing import Literal
from pydantic import BaseModel, ConfigDict, Field, model_validator
Text = str

class BandSource(BaseModel):
    """Where a band came from, in enough detail that a reader can go and check it."""
    model_config = ConfigDict(extra='forbid', frozen=True)
    kind: Literal['handbook', 'peer_reviewed', 'primary_study', 'plant_data', 'derived', 'owner_ruling']
    citation: Text = Field(min_length=8)
    locator: Text | None = None
    doi: Text | None = None

    @model_validator(mode='after')
    def _citation_is_checkable(self) -> 'BandSource':
        if self.kind == 'handbook' and (not self.locator):
            raise ValueError(f'{self.citation}: a handbook band names its edition and table or section; a bare title cannot be checked and so cannot establish a band')
        if self.kind == 'derived' and (not self.locator):
            raise ValueError(f'{self.citation}: a derived value states its derivation')
        if self.kind in {'peer_reviewed', 'primary_study'} and (not (self.doi or self.locator)):
            raise ValueError(f'{self.citation}: a study needs its DOI or a resolvable handle')
        if self.kind == 'owner_ruling' and (not self.locator):
            raise ValueError(f'{self.citation}: an owner ruling names who ruled it and when; the house owns its criteria, but a criterion nobody can trace is not a house one')
        return self
