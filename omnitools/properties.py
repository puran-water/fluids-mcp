"""Unified property lookup for fluids and pipes."""

from typing import Optional, Literal
from tools.fluid_properties import get_fluid_properties, list_available_fluids
from tools.pipe_properties import get_pipe_properties


def properties(
    lookup_type: Literal["fluid", "pipe", "list_fluids"] = "fluid",
    
    # Fluid properties parameters
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = 1.0,
    
    # Pipe properties parameters
    nominal_size: Optional[float] = None,
    schedule: str = "40",
    inner_diameter: Optional[float] = None,
    outer_diameter: Optional[float] = None,
    material: Optional[str] = None,
) -> str:
    """Unified property lookup for fluids and pipes.
    
    Provides access to thermodynamic properties and dimensions:
    - lookup_type='fluid': Get fluid properties at specified conditions
    - lookup_type='pipe': Get pipe dimensions and roughness
    - lookup_type='list_fluids': List all available fluids
    
    Args:
        lookup_type: Type of property lookup
        
        Fluid parameters:
            fluid_name: Name of the fluid (e.g., "Water", "Air", "Methane")
            temperature_c: Temperature in Celsius
            pressure_bar: Pressure in bar (default: 1.0)
        
        Pipe parameters:
            nominal_size: Nominal pipe size in inches
            schedule: Pipe schedule designation
            inner_diameter: Inner diameter in m (to find closest)
            outer_diameter: Outer diameter in m (to find closest)
            material: Pipe material name
    
    Returns:
        JSON string with property data
    
    Examples:
        Fluid properties:
        >>> properties(lookup_type="fluid", fluid_name="Water", 
        ...           temperature_c=25, pressure_bar=1)
        
        Pipe properties:
        >>> properties(lookup_type="pipe", nominal_size=2, schedule="40")
        
        List fluids:
        >>> properties(lookup_type="list_fluids")
    """
    if lookup_type == "fluid":
        if not fluid_name or temperature_c is None:
            return '{"error": "fluid_name and temperature_c required for fluid lookup"}'
        return get_fluid_properties(
            fluid_name=fluid_name,
            temperature_c=temperature_c,
            pressure_bar=pressure_bar
        )
    
    elif lookup_type == "pipe":
        return get_pipe_properties(
            nominal_size=nominal_size,
            schedule=schedule,
            inner_diameter=inner_diameter,
            outer_diameter=outer_diameter,
            material=material
        )
    
    elif lookup_type == "list_fluids":
        return list_available_fluids()
    
    else:
        return f'{{"error": "Invalid lookup_type: {lookup_type}"}}'