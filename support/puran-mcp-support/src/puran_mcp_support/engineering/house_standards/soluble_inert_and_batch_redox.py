'Public engineering calculation interface.'
from __future__ import annotations
from puran_mcp_support.engineering.house_standards.ruling_citation import house_ruling
RULING_ID = 'PUBLIC-RULE-75274692f42e'
RULING_DATE = '2026-10-15'
SOLUBLE_INERT_N_CONTENT = 0.01
MODIFIED_ASM2D_I_NSI = 0.06003
SOLUBLE_INERT_N_SOURCE = f"{house_ruling(RULING_ID, RULING_DATE)}: the N content of soluble inert COD is ASM2d's i_N,SI = 0.01 g N/g COD (Henze, Gujer, Mino, Matsuo, Wentzel, Marais and van Loosdrecht 1999, Activated Sludge Model No.2d, ASM2d, Water Sci. Technol. 39(1):165-182, doi:10.2166/wst.1999.0036, as WaterTAP's asm2d_properties.py and asm2d_reactions.py implement it), the nominal the organic-N fit scales; not ModifiedASM2d's 0.06003, the particulate inert's content the BSM2-P harmonisation gave both inerts (Flores-Alsina et al. 2016; Solon et al. 2017)"
BATCH_REDOX_STATEMENT = f"ammonium and sulfide decoupled from redox in the equilibrium (phreeqc.dat's Amm and a sulfide element derived the same way; no redox state moves: {house_ruling(RULING_ID, RULING_DATE)})"
RULING_PROVENANCE: dict[str, str] = {RULING_ID: 'Published soluble-inert and held-redox model conventions; this metadata grants no approval.'}
__all__ = ['BATCH_REDOX_STATEMENT', 'MODIFIED_ASM2D_I_NSI', 'RULING_DATE', 'RULING_ID', 'RULING_PROVENANCE', 'SOLUBLE_INERT_N_CONTENT', 'SOLUBLE_INERT_N_SOURCE']
