import json
import logging
import math
from typing import Optional

# Import the fluids library - this typically makes submodules available
import fluids

# Configure logging
logger = logging.getLogger("fluids-mcp.get_pipe_properties")

def get_pipe_properties(
    nominal_size: Optional[float] = None,  # Nominal pipe size in inches
    schedule: str = "40",        # Pipe schedule designation
    inner_diameter: Optional[float] = None, # Inner diameter in m
    outer_diameter: Optional[float] = None, # Outer diameter in m
    material: Optional[str] = None  # Pipe material name
) -> str:
    """Look up standard pipe properties.
    
    Args:
        nominal_size: Nominal pipe size in inches
        schedule: Pipe schedule designation (40, 80, etc.)
        inner_diameter: Inner diameter in m (to find closest pipe)
        outer_diameter: Outer diameter in m (to find closest pipe)
        material: Pipe material name (e.g., "Steel", "Cast iron", "Concrete", etc.)
    
    Returns:
        Standard pipe dimensions and properties
    """
    try:
        result = {}
        
        # Get pipe dimensions
        if nominal_size is not None or inner_diameter is not None or outer_diameter is not None:
            if nominal_size is not None:
                NPS, Di, Do, t = fluids.piping.nearest_pipe(NPS=nominal_size, schedule=schedule)
            elif inner_diameter is not None:
                NPS, Di, Do, t = fluids.piping.nearest_pipe(Di=inner_diameter, schedule=schedule)
            elif outer_diameter is not None:
                NPS, Di, Do, t = fluids.piping.nearest_pipe(Do=outer_diameter, schedule=schedule)
            
            # Calculate other properties
            area = math.pi * (Di/2)**2
            circumference = math.pi * Di
            
            # Add pipe dimensions to result
            result.update({
                "nominal_pipe_size_inches": NPS,
                "inner_diameter_m": round(Di, 6),
                "outer_diameter_m": round(Do, 6),
                "wall_thickness_m": round(t, 6),
                "cross_sectional_area_m2": round(area, 6),
                "inner_circumference_m": round(circumference, 6),
                "schedule": schedule
            })
        
        # Get roughness information if material is provided
        if material is not None:
            # Get the roughness dictionary from fluids library
            from fluids.friction import _roughness
            
            # Create a case-insensitive lookup
            materials_lower = {k.lower(): k for k in _roughness.keys()}
            
            if material.lower() in materials_lower:
                # Get the actual case-sensitive key
                actual_key = materials_lower[material.lower()]
                roughness_value = _roughness[actual_key]
                
                # Add roughness info to result
                result.update({
                    "material": actual_key,
                    "absolute_roughness_m": roughness_value
                })
                
                # Calculate relative roughness if we have pipe dimensions
                if "inner_diameter_m" in result:
                    relative_roughness = roughness_value / result["inner_diameter_m"]
                    result["relative_roughness"] = round(relative_roughness, 6)
            else:
                # Material not found, provide list of available materials
                all_materials = list(_roughness.keys())
                all_materials.sort()
                
                result.update({
                    "error_material": f"Material '{material}' not found in database",
                    "available_materials": all_materials
                })
        
        # Check if we have any results
        if not result:
            return json.dumps({
                "error": "Must provide either pipe dimensions (nominal_size, inner_diameter, or outer_diameter) or material name"
            })
        
        result["successful"] = True
        return json.dumps(result)
        
    except Exception as e:
        logger.error(f"Error in get_pipe_properties: {e}", exc_info=True)
        return json.dumps({
            "error": f"Calculation error: {str(e)}",
            "successful": False
        })
