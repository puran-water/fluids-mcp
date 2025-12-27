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

# Default values for pipe calculations
DEFAULT_ROUGHNESS = 1.5e-5   # Default pipe roughness for liquid calculations, m (smooth commercial steel)
DEFAULT_GAS_ROUGHNESS = 4.5e-5  # Default pipe roughness for gas calculations, m (typical commercial steel)
DEFAULT_ATMOSPHERIC_PRESSURE = 101325.0  # Default atmospheric pressure, Pa

# Default gas properties (for Air at standard conditions)
DEFAULT_GAS_MW = 28.96       # Molecular weight of air, kg/kmol
DEFAULT_GAS_GAMMA = 1.4      # Specific heat ratio Cp/Cv for diatomic gases (air, N2, O2)
DEFAULT_GAS_Z_FACTOR = 1.0   # Compressibility factor for ideal gas

# Special gamma values for common gases
GAMMA_METHANE = 1.3          # Cp/Cv for methane and biogas
GAMMA_CO2 = 1.3              # Cp/Cv for carbon dioxide
GAMMA_MONATOMIC = 1.67       # Cp/Cv for monatomic gases (He, Ar)

# Default liquid properties
DEFAULT_WATER_VAPOR_PRESSURE = 3170.0  # Saturation pressure of water at 25°C, Pa
DEFAULT_WATER_DENSITY = 1000.0  # Density of water at 25°C, kg/m³
DEFAULT_WATER_VISCOSITY = 0.001  # Dynamic viscosity of water at 20°C, Pa·s

# Friction factor defaults (when correlation fails)
DEFAULT_FRICTION_FACTOR = 0.02  # Typical fully turbulent friction factor

# Solver bounds (for numerical calculations)
SOLVER_VELOCITY_MIN = 0.01   # Minimum velocity for flow calculations, m/s
SOLVER_VELOCITY_MAX = 30.0   # Maximum velocity for flow calculations, m/s
SOLVER_FLOW_MIN = 1e-8       # Minimum flow rate, m³/s
SOLVER_FLOW_MAX = 100.0      # Maximum flow rate, m³/s
