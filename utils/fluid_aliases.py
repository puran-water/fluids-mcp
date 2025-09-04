"""
Fluid name aliases and mapping for common fluid names to CoolProp-compatible names.
"""

import logging
from typing import Optional

logger = logging.getLogger("fluids-mcp.fluid_aliases")

# Mapping of common aliases to CoolProp names
FLUID_NAME_MAP = {
    # Ammonia variants
    'ammonia': 'Ammonia',
    'nh3': 'Ammonia',
    'r717': 'Ammonia',
    'r-717': 'Ammonia',
    
    # Natural gas / methane
    'natural gas': 'Methane',
    'naturalgas': 'Methane',
    'natural_gas': 'Methane',
    'ng': 'Methane',
    'ch4': 'Methane',
    
    # Common refrigerants
    'r134a': 'R134a',
    'r-134a': 'R134a',
    'r404a': 'R404A',
    'r-404a': 'R404A',
    'r410a': 'R410A',
    'r-410a': 'R410A',
    
    # Water variants
    'h2o': 'Water',
    'steam': 'Water',
    
    # Air
    'compressed air': 'Air',
    'compressed_air': 'Air',
    
    # CO2 variants
    'carbon dioxide': 'CarbonDioxide',
    'carbondioxide': 'CarbonDioxide',
    'co2': 'CarbonDioxide',
    'r744': 'CarbonDioxide',
    'r-744': 'CarbonDioxide',
    
    # Nitrogen
    'n2': 'Nitrogen',
    
    # Oxygen
    'o2': 'Oxygen',
    
    # Hydrogen
    'h2': 'Hydrogen',
    
    # Sulfur compounds
    'sulfur dioxide': 'SulfurDioxide',
    'sulphur dioxide': 'SulfurDioxide',
    'so2': 'SulfurDioxide',
    
    # Ethanol
    'ethyl alcohol': 'Ethanol',
    'alcohol': 'Ethanol',
    'etoh': 'Ethanol',
    
    # Methanol
    'methyl alcohol': 'Methanol',
    'meoh': 'Methanol',
    
    # Propane
    'lpg': 'Propane',  # Simplified - LPG is actually a mixture
    'c3h8': 'Propane',
    'r290': 'Propane',
    
    # Glycols (commented out - FluidProperties doesn't support INCOMP:: fluids)
    # 'glycol': 'INCOMP::MEG[0.3]',  # Default to 30% monoethylene glycol
    # 'ethylene glycol': 'INCOMP::MEG[0.3]',
    # 'ethyleneglycol': 'INCOMP::MEG[0.3]',
    # 'eg': 'INCOMP::MEG[0.3]',
    # 'meg': 'INCOMP::MEG[0.3]',
    # 'propylene glycol': 'INCOMP::MPG[0.3]',
    # 'propyleneglycol': 'INCOMP::MPG[0.3]',
    # 'pg': 'INCOMP::MPG[0.3]',
    # 'mpg': 'INCOMP::MPG[0.3]',
}


def normalize_fluid_name(name: str) -> str:
    """
    Normalize a fluid name by removing spaces, converting to lowercase, etc.
    
    Args:
        name: Raw fluid name input
        
    Returns:
        Normalized name for lookup
    """
    if not name:
        return name
    
    # Convert to lowercase and strip whitespace
    normalized = name.lower().strip()
    
    # Remove hyphens in refrigerant names
    if normalized.startswith('r-'):
        normalized = 'r' + normalized[2:]
    
    return normalized


def map_fluid_name(name: str, warn_on_fallback: bool = True) -> str:
    """
    Map a common fluid name or alias to a CoolProp-compatible name.
    
    Args:
        name: Input fluid name (can be an alias)
        warn_on_fallback: Whether to log warnings for approximations
        
    Returns:
        CoolProp-compatible fluid name
    """
    if not name:
        return name
    
    # First check if it's already a valid CoolProp name (case-sensitive)
    # This includes INCOMP:: and HEOS:: prefixes
    if name.startswith('INCOMP::') or name.startswith('HEOS::'):
        return name
    
    # Normalize the input
    normalized = normalize_fluid_name(name)
    
    # Check the alias map
    if normalized in FLUID_NAME_MAP:
        mapped_name = FLUID_NAME_MAP[normalized]
        
        # Log warnings for approximations
        if warn_on_fallback:
            if normalized in ['natural gas', 'naturalgas', 'natural_gas', 'ng']:
                logger.warning(
                    f"'{name}' mapped to 'Methane'. This is an approximation. "
                    "For accurate natural gas calculations, provide gas composition."
                )
            elif normalized in ['glycol', 'ethylene glycol', 'ethyleneglycol', 'eg', 'meg',
                               'propylene glycol', 'propyleneglycol', 'pg', 'mpg']:
                logger.warning(
                    f"'{name}' mapped to '{mapped_name}'. "
                    "This assumes 30% concentration. Adjust if needed."
                )
            elif normalized == 'lpg':
                logger.warning(
                    f"'{name}' mapped to 'Propane'. "
                    "LPG is typically a propane/butane mixture. This is an approximation."
                )
        
        return mapped_name
    
    # If no alias found, return the original name (CoolProp will validate)
    # Try with proper capitalization for common patterns
    if normalized == name.lower():
        # Try capitalizing first letter
        capitalized = name[0].upper() + name[1:] if len(name) > 1 else name.upper()
        return capitalized
    
    return name


def get_glycol_fluid_string(glycol_type: str = 'MEG', concentration_percent: float = 30.0) -> str:
    """
    Create a CoolProp incompressible fluid string for glycol solutions.
    
    Args:
        glycol_type: 'MEG' for monoethylene glycol or 'MPG' for monopropylene glycol
        concentration_percent: Mass concentration in percent (0-100)
        
    Returns:
        CoolProp incompressible fluid string
    """
    if glycol_type.upper() not in ['MEG', 'MPG']:
        raise ValueError(f"Invalid glycol type: {glycol_type}. Use 'MEG' or 'MPG'.")
    
    if not 0 <= concentration_percent <= 100:
        raise ValueError(f"Concentration must be between 0 and 100%, got {concentration_percent}")
    
    # Convert percentage to fraction (0-1)
    fraction = concentration_percent / 100.0
    return f"INCOMP::{glycol_type.upper()}[{fraction:.2f}]"


def get_mixture_string(components: dict) -> str:
    """
    Create a CoolProp HEOS mixture string from component fractions.
    
    Args:
        components: Dictionary of component names to mole fractions
                   e.g., {'Methane': 0.9, 'Ethane': 0.05, 'Propane': 0.05}
        
    Returns:
        CoolProp HEOS mixture string
    """
    if not components:
        raise ValueError("No components provided for mixture")
    
    # Normalize fractions
    total = sum(components.values())
    if abs(total - 1.0) > 0.01:
        logger.warning(f"Component fractions sum to {total}, normalizing to 1.0")
        components = {k: v/total for k, v in components.items()}
    
    # Build mixture string
    mixture_parts = []
    for component, fraction in components.items():
        # Map component names if needed
        mapped_component = map_fluid_name(component, warn_on_fallback=False)
        mixture_parts.append(f"{mapped_component}[{fraction:.6f}]")
    
    return "HEOS::" + "&".join(mixture_parts)