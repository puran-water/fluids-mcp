from __future__ import annotations
from datetime import date as _Date
from typing import Annotated, Any, Literal, Mapping
from pydantic import BaseModel, ConfigDict, Field, model_validator

class DeviationAuthority(BaseModel):
    """Who ruled the deviation, or an explicit statement that nobody has yet.

    A design run can propose a deviation; it cannot grant one.  The unruled
    state is therefore representable and is the only shape a run may author
    for itself, because a required ``person`` with no alternative makes
    typing somebody's name the path of least resistance -- which is how a
    ruling nobody gave reaches a workbook, a schedule and a client.  Naming a
    person is an affirmative act and must cite where the ruling is recorded.
    """
    model_config = ConfigDict(extra='forbid', frozen=True, str_strip_whitespace=True)
    status: Literal['proposed_pending_ruling', 'ruled']
    person: Annotated[str, Field(min_length=1)] | None = None
    date: _Date | None = None
    ruling_ref: Annotated[str, Field(min_length=1)] | None = None
    proposed_by: Annotated[str, Field(min_length=1)] | None = None

    @model_validator(mode='after')
    def _authority_is_earned(self):
        if self.status == 'ruled':
            missing = [f for f in ('person', 'date', 'ruling_ref') if getattr(self, f) is None]
            if missing:
                raise ValueError('a ruled guideline deviation requires ' + ', '.join(missing) + '; ruling_ref must cite where the ruling is recorded (artifact id, document path or dated instruction) -- a name alone is not authority')
        else:
            named = [f for f in ('person', 'date', 'ruling_ref') if getattr(self, f) is not None]
            if named:
                raise ValueError('a deviation pending a ruling must not name ' + ', '.join(named) + "; set status='ruled' with a ruling_ref once the owner has ruled")
            if self.proposed_by is None:
                raise ValueError('a proposed deviation must name proposed_by (the run that raised it)')
        return self
