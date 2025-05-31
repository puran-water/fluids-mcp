# Fluids MCP Server

A Model Context Protocol (MCP) server providing comprehensive fluid mechanics and hydraulics calculations for industrial applications, with a focus on wastewater treatment systems.

## Features

### Gas Flow Calculations
- **Gas Pipe Pressure Drop**: Calculate pressure losses in gas piping systems using industry-standard methods (Weymouth, Panhandle A/B)
- **Blower/Compressor Sizing**: Determine power requirements and discharge conditions for gas compression
- **Gas Control Valve Sizing**: Size control valves for gas service per IEC 60534 standards

### Liquid Flow Calculations  
- **Pipe Pressure Drop**: Calculate friction losses and head requirements for liquid piping systems
- **Pump Requirements**: Determine Total Dynamic Head (TDH) and Net Positive Suction Head Available (NPSHa)
- **Liquid Control Valve Sizing**: Size control valves for liquid service applications

### Fluid Properties
- **Property Lookup**: Access thermodynamic properties via CoolProp integration
- **Common Fluids**: Support for water, air, methane, and 120+ other fluids
- **Temperature/Pressure Dependent**: Properties calculated at actual operating conditions

### Additional Tools
- **Reynolds Number**: Flow regime determination
- **Pipe Properties**: Standard pipe dimensions lookup (ASME schedules)
- **Open Channel Flow**: Hydraulic calculations for open channels

## Installation

### Prerequisites
- Python 3.10+
- MCP client (Claude Desktop or compatible)

### Setup

1. Clone the repository:
```bash
git clone https://github.com/puran-water/fluids-mcp.git
cd fluids-mcp
```

2. Create a virtual environment:
```bash
python -m venv venv
# Windows
venv\Scripts\activate
# Linux/Mac
source venv/bin/activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Configuration

### Claude Desktop Configuration

Add to your Claude Desktop configuration file:

**Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
**Linux**: `~/.config/claude/claude_desktop_config.json`

**Windows Example:**
```json
{
  "mcpServers": {
    "fluids-mcp": {
      "command": "python",
      "args": [
        "C:\\Users\\YOUR_USERNAME\\mcp-servers\\fluids-mcp\\server.py"
      ],
      "env": {}
    }
  }
}
```

**macOS/Linux Example:**
```json
{
  "mcpServers": {
    "fluids-mcp": {
      "command": "python",
      "args": [
        "/path/to/your/mcp-servers/fluids-mcp/server.py"
      ],
      "env": {}
    }
  }
}
```

## Usage Examples

### Example 1: Pump Sizing for Water Transfer
```
Calculate pump requirements for:
- Flow rate: 500 GPM
- Suction pipe: 50 ft of 6" steel pipe  
- Discharge pipe: 1000 ft of 6" steel pipe
- Static lift: 110 ft (10 ft suction lift + 100 ft discharge)
- Fluid: Water at 20°C
```

### Example 2: Aeration Blower Sizing
```
Size an aeration blower for:
- Air flow: 5000 m³/hr at normal conditions
- Inlet: Atmospheric pressure (14.7 psi)
- Discharge: 22 psi (for 5m water depth + losses)
- Ambient temperature: 20°C
- Required: Power consumption and discharge temperature
```

### Example 3: Biogas Pipe Sizing
```
Calculate pressure drop for biogas piping:
- Flow rate: 500 m³/hr (normal conditions)
- Pipe: 200m of 6" Schedule 40 steel
- Fluid: Methane at 35°C (biogas approximation)
- Inlet pressure: 103 kPa (1.5 psig)
```

## API Reference

### Tools Available

| Tool | Description | Key Parameters |
|------|-------------|----------------|
| `calculate_pipe_pressure_drop` | Liquid pipe friction losses | flow_rate, pipe_diameter, length, fluid |
| `calculate_pump_requirements` | Pump TDH and NPSHa | flow_rate, pipe_config, static_heads |
| `calculate_gas_pipe_pressure_drop` | Gas pipe pressure drop | flow_rate, pressure, temperature, method |
| `calculate_blower_compressor_requirements` | Compressor power | flow_rate, pressures, efficiency |
| `calculate_gas_control_valve` | Gas valve Cv | flow_rate, pressures, valve_type |
| `calculate_liquid_control_valve` | Liquid valve Cv | flow_rate, pressure_drop, fluid |
| `get_fluid_properties` | Thermodynamic properties | fluid_name, temperature, pressure |
| `calculate_reynolds_number` | Flow regime | velocity, diameter, fluid_properties |
| `get_pipe_properties` | Pipe dimensions | nominal_size, schedule, material |
| `calculate_open_channel_flow` | Open channel hydraulics | channel_type, dimensions, slope |

## Technical Details

### Calculation Methods
- **Pipe Friction**: Darcy-Weisbach equation with Colebrook-White friction factor
- **Gas Flow**: Weymouth, Panhandle A, and Panhandle B equations
- **Valve Sizing**: IEC 60534 / ISA-75.01.01 standards
- **Fluid Properties**: CoolProp library with NIST-validated equations of state

### Unit Systems
- **Flexible Units**: Accepts common engineering units (GPM, psi, ft, m³/hr, etc.)
- **Internal Calculations**: SI units for accuracy
- **Output**: User-preferred units with clear labeling

### Latest Version (v1.2.0)
- Professional-grade testing with comprehensive test suite
- Robust property lookup for all common fluids including methane, air, water
- Accurate unit conversions and valve sizing calculations
- Enhanced package structure for reliable deployment
- Cross-platform compatibility

## Development

### Project Structure
```
fluids-mcp/
├── server.py              # Main MCP server implementation
├── tools/                 # Calculation tools
│   ├── pipe_pressure_drop.py
│   ├── pump_requirements.py
│   ├── gas_pipe_pressure_drop.py
│   ├── blower_compressor.py
│   ├── gas_control_valve.py
│   ├── liquid_control_valve.py
│   ├── fluid_properties.py
│   ├── reynolds_number.py
│   ├── pipe_properties.py
│   └── open_channel_flow.py
├── utils/                 # Shared utilities
│   ├── constants.py      # Unit conversions and constants
│   ├── helpers.py        # Common functions
│   └── import_helpers.py # Optional dependency handling
├── pydraulics/           # Open channel hydraulics
└── tests/                # Test suite
```

### Testing

Run the test suite:
```bash
pytest tests/ -v
```

For coverage reporting:
```bash
pytest tests/ --cov=tools --cov=utils --cov=pydraulics --cov-report=term-missing
```

### Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Built on the [Fluids](https://github.com/CalebBell/fluids) library by Caleb Bell
- Thermodynamic properties provided by [CoolProp](http://www.coolprop.org/)
- Implements the [Model Context Protocol](https://modelcontextprotocol.io) specification

## Support

For issues, questions, or contributions, please visit our [GitHub repository](https://github.com/puran-water/fluids-mcp).

## Roadmap

### Planned Features
- MCP Resources for fluid property tables and pipe standards
- MCP Prompts for complete system design workflows
- Additional calculation methods (Hazen-Williams, Manning equation)
- Equipment selection databases
- Integration with other Puran Water MCP servers

### Future Enhancements
- Web-based testing interface
- Calculation report generation
- Multi-fluid mixture properties
- Economic optimization tools