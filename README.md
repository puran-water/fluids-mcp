# Fluids MCP Server

[![MCP](https://img.shields.io/badge/MCP-1.0-blue)](https://modelcontextprotocol.io)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org)
[![CoolProp](https://img.shields.io/badge/CoolProp-6.6.0-green)](http://www.coolprop.org)
[![Fluids](https://img.shields.io/badge/Fluids-1.0.26-green)](https://github.com/CalebBell/fluids)
[![Thermo](https://img.shields.io/badge/Thermo-0.3.0-green)](https://github.com/CalebBell/thermo)
[![License](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Production-success)](https://github.com/puran-water/fluids-mcp)

Model Context Protocol server providing comprehensive fluid mechanics, thermodynamics, and hydraulic calculations for industrial process engineering.

## Overview

This MCP server provides 6 consolidated omnitools that encapsulate 16+ specialized calculation methods:

1. **pipe_flow**: Unified liquid and gas pipe flow calculations with pressure drop analysis
2. **control_valve**: Valve sizing for both liquid and gas service per IEC 60534
3. **pipe_sizing**: Optimal pipe diameter selection based on constraints
4. **parameter_sweep**: Performance optimization through systematic parameter variation
5. **properties**: Thermodynamic and physical property lookup for 120+ fluids and pipe dimensions
6. **machine_requirements**: Pump, compressor, and hydraulic system design calculations

## Technical Capabilities

### Thermodynamic Properties
- **CoolProp Integration**: NIST-validated equations of state for 120+ pure fluids
- **Mixture Support**: Composition-based property calculation via thermo.Mixture
- **Fallback Strategy**: Automatic failover to thermo.Chemical for incomplete datasets
- **Temperature/Pressure Dependent**: Real-time property calculation at operating conditions
- **Unit Safety**: Explicit unit suffixes (no suffix = SI, suffix = alternative units)

### Pipe Flow Analysis
- **Liquid Phase**: Darcy-Weisbach with Colebrook-White friction factor
- **Gas Phase**: Multiple correlations (Weymouth, Panhandle A/B, IGT, Oliphant, Spitzglass, isothermal)
- **Solve-for Capabilities**: Automatic solution for any unknown (P1, P2, L, Q, D)
- **Fitting Losses**: K-factor and equivalent length methods with extensive fitting library
- **Mixture Support**: Gas composition handling with proper mixing rules

### Equipment Sizing
- **Pumps**: TDH calculation with velocity head corrections for different nozzle sizes
- **Compressors**: Power requirements with isentropic/polytropic efficiency models
- **Control Valves**: Cv/Kv calculation per IEC 60534-2-1 standards
- **NPSH Analysis**: Available NPSH calculation with vapor pressure corrections
- **Reynolds Number**: Flow regime determination for any fluid
- **Open Channel Flow**: Manning equation and normal depth calculations

### Optimization Tools
- **Parameter Sweeps**: Vectorized calculations for 10,000+ data points
- **Design Space Exploration**: Systematic variation of operating parameters
- **Constraint-Based Sizing**: Automatic pipe selection within pressure drop limits
- **Multi-Variable Analysis**: Sweep any parameter while holding others constant

## Installation

### Prerequisites
- Python 3.10+ (3.12 tested)
- Virtual environment
- MCP client (Claude Desktop or compatible)

### Setup
```bash
git clone https://github.com/puran-water/fluids-mcp.git
cd fluids-mcp
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## MCP Configuration

### Claude Desktop Integration

Add to your configuration file:

**Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
**Linux**: `~/.config/claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "fluids-mcp": {
      "command": "python",
      "args": ["/absolute/path/to/fluids-mcp/server.py"],
      "env": {}
    }
  }
}
```

## API Reference

### Registered Omnitools

#### pipe_flow
Unified pipe pressure drop calculations for liquid and gas phases.

**Parameters:**
- `phase`: "liquid" or "gas" - Determines calculation method
- **Common**: `pipe_length` (m), `pipe_diameter` (m), `nominal_size_in` (inches), `schedule`, `material`, `pipe_roughness`, `fittings`, `temperature_c`
- **Liquid**: `flow_rate` (m³/s), `flow_rate_gpm` (GPM), `fluid_density`, `fluid_viscosity`, `pressure_drop`
- **Gas**: `inlet_pressure`, `outlet_pressure` (Pa absolute), `flow_rate_kg_s`, `flow_rate_norm_m3_hr`, `gas_mw`, `gas_gamma`, `gas_z_factor`, `gas_composition_mol`, `method`

**Returns:** JSON with pressure drop, velocity, Reynolds number, fitting losses

#### control_valve
Valve sizing calculations per IEC 60534 standards.

**Parameters:**
- `phase`: "liquid" or "gas" - Service type
- **Common**: `inlet_pressure`, `outlet_pressure`, `temperature_c`, `valve_type`, `size_units`
- **Liquid**: `flow_rate` (m³/s), `fluid_density`, `fluid_viscosity`, `fluid_saturation_pressure`
- **Gas**: `flow_rate_kg_s`, `gas_mw`, `gas_gamma`, `gas_z_factor`

**Returns:** JSON with Cv, Kv, recommended valve size, flow characteristics

#### pipe_sizing
Optimal pipe diameter selection based on constraints.

**Parameters:**
- `phase`: "liquid" or "gas" - Fluid type
- `flow_rate` or variant: Flow rate in appropriate units
- `pipe_length`: Length in meters
- `allowable_dp_pa`: Maximum pressure drop in Pa
- `velocity_min_m_s`, `velocity_max_m_s`: Velocity constraints
- `nps_candidates`: List of nominal sizes to evaluate
- `return_all_candidates`: Return all evaluated sizes

**Returns:** JSON with recommended NPS, actual pressure drop, velocity

#### parameter_sweep
Systematic parameter variation for optimization studies.

**Parameters:**
- `calculation`: "pipe_liquid", "pipe_gas", or "blower"
- `variable`: Parameter to vary (e.g., "flow_rate", "pipe_diameter")
- `start`, `stop`: Range for variation
- `n`: Number of points
- **Base parameters**: Fixed parameters for the calculation

**Returns:** JSON with arrays of sweep results

#### properties
Thermodynamic and physical property lookup.

**Parameters:**
- `lookup_type`: "fluid", "pipe", or "list_fluids"
- **Fluid**: `fluid_name`, `temperature_c`, `pressure_bar`
- **Pipe**: `nominal_size`, `schedule`, `material`

**Returns:** 
- Fluid: Comprehensive properties (density, viscosity, Cp, Cv, thermal conductivity, etc.)
- Pipe: Dimensions (OD, ID, thickness, roughness)
- List: Available fluid names

#### machine_requirements
Pump, compressor, and hydraulic calculations.

**Parameters:**
- `machine_type`: "pump", "compressor", "reynolds", or "open_channel"
- **Pump**: Flow rate, static heads, pipe configurations, nozzle sizes
- **Compressor**: Flow rate, pressures, efficiency, temperature
- **Reynolds**: Velocity, characteristic length, fluid properties
- **Open Channel**: Channel geometry, slope, Manning coefficient

**Returns:** 
- Pump: TDH, NPSHA, power requirements
- Compressor: Discharge temperature, power, efficiency
- Reynolds: Reynolds number, flow regime
- Open Channel: Flow rate, velocity, normal depth

## Technical Implementation

### Architecture
```
fluids-mcp/
├── server.py                # MCP server with 6 omnitool registrations
├── omnitools/              # Consolidated wrapper tools
│   ├── pipe_flow.py        # Wraps liquid/gas pipe calculations
│   ├── control_valve.py    # Wraps liquid/gas valve sizing
│   ├── pipe_sizing.py      # Wraps optimal sizing algorithms
│   ├── parameter_sweep.py  # Wraps all sweep functions
│   ├── properties.py       # Wraps fluid/pipe properties
│   └── machine_requirements.py  # Wraps pump/compressor/hydraulics
├── tools/                  # 16+ core calculation engines
│   ├── pipe_pressure_drop.py       # Liquid pipe calculations
│   ├── gas_pipe_pressure_drop.py   # Gas pipe with mixtures
│   ├── pump_requirements.py        # Pump TDH and NPSH
│   ├── blower_compressor.py        # Compressor power
│   ├── liquid_control_valve.py     # Liquid Cv calculation
│   ├── gas_control_valve.py        # Gas Cv calculation
│   ├── fluid_properties.py         # CoolProp interface
│   ├── pipe_properties.py          # ASME dimensions
│   ├── pipe_sizing.py              # Sizing algorithms
│   ├── reynolds_number.py          # Flow regime
│   ├── open_channel_flow_new.py    # Channel hydraulics
│   └── [sweep functions]            # Optimization tools
├── utils/                  # Shared utilities
│   ├── constants.py        # Unit conversions (GPM_to_M3S, PSI_to_PA, etc.)
│   ├── helpers.py          # Common functions (get_fitting_K)
│   ├── import_helpers.py   # Dependency management
│   ├── property_cache.py   # Performance optimization
│   └── fluid_aliases.py    # Name mapping ("natural gas" → "Methane")
├── pydraulics/            # Open channel hydraulics engine
├── tests/                 # Comprehensive test suite
└── UNITS.md              # Unit convention documentation
```

### Key Design Patterns

#### Omnitool Pattern
Each omnitool acts as a discriminator-based wrapper:
```python
def pipe_flow(phase: Literal["liquid", "gas"], **kwargs):
    if phase == "liquid":
        return calculate_pipe_pressure_drop(**filtered_kwargs)
    else:
        return calculate_gas_pipe_pressure_drop(**filtered_kwargs)
```

#### Unit Convention
- **Base SI Units**: Parameters without suffix use SI units
- **Alternative Units**: Suffixed parameters indicate unit type
- **Examples**:
  - `flow_rate` → m³/s
  - `flow_rate_gpm` → US gallons per minute
  - `pressure` → Pa (absolute)
  - `pressure_psi` → pounds per square inch
  - `temperature_c` → Celsius (converted to K internally)

#### Property Resolution Hierarchy
1. User-specified values (highest priority)
2. Mixture calculations via thermo.Mixture (for gas_composition_mol)
3. FluidProp/CoolProp lookup (primary source)
4. thermo.Chemical fallback (for NaN or missing properties)
5. Engineering defaults with warnings

#### Fluid Name Aliasing
Common names are automatically mapped:
- "natural gas" → "Methane"
- "ammonia" → "Ammonia"
- "nh3" → "Ammonia"
- "co2" → "CarbonDioxide"
- "biogas" → Use gas_composition_mol instead

### Performance Optimizations
- **Vectorized Calculations**: NumPy arrays for parameter sweeps
- **Property Caching**: Avoid redundant thermodynamic calculations
- **Lazy Imports**: Optional dependencies loaded on demand
- **Parameter Filtering**: inspect.signature() ensures valid parameters

## Calculation Methods

### Fluid Flow
- **Friction Factor**: Colebrook-White equation with Swamee-Jain approximation
- **Gas Flow**: Weymouth, Panhandle A/B, IGT, Oliphant, Spitzglass methods
- **Two-Phase Flow**: Homogeneous equilibrium model
- **Fitting Losses**: Crane TP-410 K-factors

### Thermodynamics
- **Pure Fluids**: CoolProp with REFPROP-quality equations
- **Mixtures**: thermo.Mixture with EOS calculations
- **Mixing Rules**: Automated selection (Kay's rule, BROKAW for viscosity)
- **Fallback**: thermo.Chemical for missing CoolProp data

### Standards Compliance
- **IEC 60534-2-1**: Control valve flow coefficients
- **ASME B36.10M**: Welded and seamless pipe dimensions
- **ISA-75.01.01**: Control valve sizing equations
- **Crane TP-410**: Flow of fluids through valves and fittings

## Usage Examples

### Example 1: Biogas Pipeline Design
```python
# Calculate pressure drop for biogas (65% CH4, 35% CO2)
pipe_flow(
    phase="gas",
    flow_rate_norm_m3_hr=150,
    pipe_length=25,
    nominal_size_in=3,
    temperature_c=25,
    inlet_pressure=101325,
    gas_composition_mol={"CH4": 0.65, "CO2": 0.35},
    fittings=[
        {"type": "elbow_90_long_radius", "quantity": 10},
        {"type": "gate_valve", "quantity": 1}
    ]
)
```

### Example 2: Pump Station Design
```python
# Size pump for water transfer system
machine_requirements(
    machine_type="pump",
    flow_rate_gpm=500,
    static_suction_head_ft=-10,
    static_discharge_head_ft=100,
    suction_pipe_length_ft=50,
    discharge_pipe_length_ft=1000,
    nominal_size_in=6,
    fluid_name="Water",
    temperature_c=20
)
```

### Example 3: Control Valve Selection
```python
# Size control valve for natural gas service
control_valve(
    phase="gas",
    flow_rate_norm_m3_hr=5000,
    inlet_pressure_psi=100,
    outlet_pressure_psi=50,
    fluid_name="natural gas",
    temperature_c=15,
    valve_type="globe"
)
```

## Testing

### Unit Tests
```bash
pytest tests/ -v
```

### Coverage Analysis
```bash
pytest tests/ --cov=tools --cov=utils --cov=omnitools --cov-report=term-missing
```

### Integration Tests
```bash
python tests/integration/test_mcp_server.py
```

## Version History

### v3.0.0 (Current)
- **Architecture**: Refactored 16 tools into 6 discriminator-based omnitools
- **Mixture Support**: Full gas composition handling via thermo.Mixture
- **Fallback Strategy**: Automatic CoolProp → thermo.Chemical failover
- **Fluid Aliasing**: Common name mapping system
- **Documentation**: Comprehensive UNITS.md reference

### v2.1.0
- **Solve-for**: Automatic unknown variable solution
- **Parameter Sweeps**: Optimization tools for all calculations
- **EOS Integration**: Enhanced thermodynamic accuracy
- **Velocity Corrections**: Pump calculations with nozzle effects

### v2.0.0
- **MCP Implementation**: Model Context Protocol support
- **CoolProp Integration**: 120+ fluid properties
- **Comprehensive Testing**: Full test coverage

## Dependencies

### Core Libraries
- **[CoolProp](http://www.coolprop.org)** - Thermodynamic and transport properties (I. Bell et al.)
- **[Fluids](https://github.com/CalebBell/fluids)** - Fluid dynamics calculations (C. Bell)
- **[Thermo](https://github.com/CalebBell/thermo)** - Chemical engineering thermodynamics (C. Bell)
- **[Chemicals](https://github.com/CalebBell/chemicals)** - Chemical property correlations (C. Bell)

### Infrastructure
- **[MCP SDK](https://github.com/modelcontextprotocol/python-sdk)** - Protocol implementation
- **[NumPy](https://numpy.org)** - Numerical computing
- **[SciPy](https://scipy.org)** - Scientific algorithms

## Contributing

Contributions are welcome. Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/enhancement`)
3. Write tests for new functionality
4. Ensure all tests pass (`pytest tests/`)
5. Update documentation as needed
6. Submit a pull request with clear description

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

This project builds upon the exceptional work of:

- **Caleb Bell** - Creator of [Fluids](https://github.com/CalebBell/fluids), [Thermo](https://github.com/CalebBell/thermo), and [Chemicals](https://github.com/CalebBell/chemicals) libraries providing industrial-strength calculations
- **Ian Bell and CoolProp Team** - Maintainers of the [CoolProp](http://www.coolprop.org) thermodynamic property library with NIST-quality correlations
- **Anthropic** - Developers of the [Model Context Protocol](https://modelcontextprotocol.io) specification enabling tool integration
- **NIST** - Reference data and validation for thermodynamic properties

Special recognition to the open-source community for maintaining these critical engineering tools.

## Support

For technical issues, feature requests, or contributions:
- GitHub Issues: [fluids-mcp/issues](https://github.com/puran-water/fluids-mcp/issues)
- Documentation: [UNITS.md](UNITS.md) for unit conventions
- Examples: See [tests/](tests/) directory for usage patterns