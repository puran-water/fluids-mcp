# Installation Guide - Fluids MCP Server v2.1

## Quick Install

```bash
pip install -r requirements.txt
```

## Enhanced Features (v2.1)

The enhanced blower compressor tool with expert reviewer recommendations requires:

### Core Dependencies
- `thermo>=0.4.0` - For equation of state calculations with Chemical objects
- `pandas>=2.0.0` - For sweep() function DataFrame output  
- `scipy>=1.10.0` - For numerical integration (quad) and optimization
- `numpy>=1.24.0` - For sweep() function array operations

### Verification Test

After installation, verify the enhanced features work:

```python
from tools.blower_compressor import calculate_blower_compressor_requirements, blower_sweep
import json

# Test EOS calculations
result = calculate_blower_compressor_requirements(
    flow_rate_norm_m3_hr=1000.0,
    inlet_pressure=101325.0,
    outlet_pressure=202650.0,
    inlet_temperature_c=20.0,
    fluid_name='Air',
    use_eos_calculations=True
)
print("EOS test successful")

# Test sweep function (returns JSON now, not DataFrame)
result_json = blower_sweep(
    variable='pressure_ratio',
    start=1.5,
    stop=3.0,
    n=5,
    flow_rate_norm_m3_hr=1000.0,
    inlet_pressure=101325.0,
    inlet_temperature_c=20.0,
    fluid_name='Air'
)
result_data = json.loads(result_json)
print(f"Sweep test successful: {result_data['summary']['successful_points']} points calculated")
```

## Troubleshooting

If you get import errors:
1. Ensure you're using Python 3.8+
2. Install thermo separately: `pip install thermo>=0.4.0`
3. For Windows, may need: `pip install --no-deps thermo` then `pip install -r requirements.txt`

## Enhanced Features in v2.1

The latest implementation includes:
- ✅ Integral-averaged Z-factor and Cp calculations using EOS
- ✅ solve_for parameter framework for all major calculations
- ✅ sweep functions (pipe_pressure_drop_sweep, gas_pipe_sweep, blower_sweep) returning JSON
- ✅ Real gas properties via Peng-Robinson EOS
- ✅ Velocity head corrections for pump calculations with different nozzle sizes
- ✅ Support for solving gas pipe diameter and length (not just pressure/flow)
- ✅ MCP-compatible parameter structure (no **kwargs)