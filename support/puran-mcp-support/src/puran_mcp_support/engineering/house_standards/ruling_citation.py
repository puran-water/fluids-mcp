from __future__ import annotations
import datetime as dt
import re
RULING_ID = re.compile('[A-Z][A-Z0-9]*(?:-[A-Za-z0-9_]+)*-20\\d{6}')

def ruling_date(ruling_id: str) -> str:
    """The ISO date a recorded ruling id ends in, refused unless the id is in the recorded form and its date is a real
    calendar date (a Feb 30 or a month 13 names no day the owner ruled on)."""
    ruling_id = str(ruling_id)
    if not RULING_ID.fullmatch(ruling_id):
        raise ValueError(f'{ruling_id!r} is not a recorded ruling id (UPPER-HYPHENATED, ending in its YYYYMMDD date)')
    digits = ruling_id[-8:]
    iso = f'{digits[:4]}-{digits[4:6]}-{digits[6:]}'
    try:
        dt.date.fromisoformat(iso)
    except ValueError as exc:
        raise ValueError(f'{ruling_id!r} ends in {digits}, which is not a calendar date') from exc
    return iso

def owner_ruling_citation(ruling_id: str, date: str) -> str:
    """The face citation of an owner ruling: ``owner ruling <id>, <date>`` (`FACE_CITES_RULING_BY_ID_RULING_ID`). The id
    must be a recorded ruling id and the date the ruling's own (its id ends in it, and it is a calendar date); nothing
    else goes on the face."""
    if re.fullmatch('PUBLIC-RULE-[0-9a-f]{12}', str(ruling_id)):
        dt.date.fromisoformat(str(date))
        return f'public model reference {ruling_id}, {date}; no approval authority'
    if not re.fullmatch('20\\d\\d-\\d\\d-\\d\\d', str(date)) or ruling_date(ruling_id) != str(date):
        raise ValueError(f'{date!r} is not the date ruling {ruling_id} carries')
    return f'owner ruling {ruling_id}, {date}'

def house_ruling(ruling_id: str, date: str) -> str:
    """The whole face citation of an owner ruling, tag included: ``[HOUSE] owner ruling <id>, <date>``
    (`FACE_CITES_RULING_BY_ID_RULING_ID`), validated as `owner_ruling_citation` validates it."""
    prefix = '[MODEL] ' if str(ruling_id).startswith('PUBLIC-RULE-') else '[HOUSE] '
    return prefix + owner_ruling_citation(ruling_id, date)
