"""
Constants used across the Fluids MCP server.

This module defines unit conversion factors and other constants used throughout the server.
"""

# Conversion factors for unit flexibility
GPM_to_M3S = 0.0000630902   # US GPM to m³/s
INCH_to_M = 0.0254           # inch to meter
PSI_to_PA = 6894.76          # psi to Pascal
CENTIPOISE_to_PAS = 0.001    # centipoise to Pa·s
FT_to_M = 0.3048             # foot to meter
LBFT3_to_KGM3 = 16.0185      # lb/ft³ to kg/m³
DEG_C_to_K = 273.15          # Celsius to Kelvin (offset)

# Standard conditions
P_NORM = 101325.0            # Normal pressure, Pa (0°C, 1 atm)
T_NORM = 273.15              # Normal temperature, K (0°C)
P_STD = 101325.0             # Standard pressure, Pa (15°C, 1 atm)
T_STD = 288.15               # Standard temperature, K (15°C)
R_UNIV = 8314.462            # Universal gas constant, J/(kmol·K)

# Physical constants
G_GRAVITY = 9.80665          # Standard gravity acceleration, m/s² (NIST value)

# Default values
DEFAULT_ROUGHNESS = 1.5e-5   # Default pipe roughness for general calculations, m
DEFAULT_GAS_ROUGHNESS = 4.5e-5  # Default pipe roughness for gas calculations, m
DEFAULT_ATMOSPHERIC_PRESSURE = 101325.0  # Default atmospheric pressure, Pa
