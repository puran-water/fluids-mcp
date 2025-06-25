# Changelog

All notable changes to the Fluids MCP Server will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.0] - 2025-06-18

### Added
- Parameter sweep functions for all major calculations (pipe_pressure_drop_sweep, gas_pipe_sweep, blower_sweep)
- Solve_for capabilities - automatically solve for any unknown variable in calculations
- Enhanced pump calculations with velocity head corrections for different nozzle sizes
- EOS-based thermodynamic calculations using thermo library for improved accuracy
- Support for solving gas pipe diameter and length (not just pressure/flow)
- Integral-averaged Z-factor and Cp calculations for compressor calculations
- Performance optimizations with property caching and vectorized calculations
- Input resolver for flexible parameter handling

### Fixed
- MCP client parameter validation errors - replaced **kwargs with explicit parameters
- Sweep functions now return JSON instead of pandas DataFrames for MCP compatibility
- Import error with blower_sweep function alias

### Changed
- All sweep functions updated to use explicit parameters instead of **base_kwargs
- Pump TDH calculations now properly account for velocity head differences between nozzles
- NPSHA calculations include velocity head and closed tank pressure corrections
- Gas pipe calculations support 5 variables (P1, P2, L, Q, D) with solve_for any 1
- Enhanced error messages with clear guidance on missing parameters

### Technical Details
- Sweep functions: Converted from pandas DataFrame returns to JSON for MCP compatibility
- Parameter structure: Eliminated **kwargs pattern for proper MCP parameter validation
- Pump velocity head: V²/2g corrections applied based on actual nozzle diameters
- Gas pipe solver: Extended from 4 to 5 variables with diameter solving capability

## [1.2.0] - 2025-01-31

### Added
- Proper pytest test suite with assertions and coverage reporting
- Cross-platform configuration examples for Claude Desktop
- Code quality improvements based on expert review
- Future improvements documentation for advanced features

### Fixed
- Removed sys.path manipulation in pydraulics imports (now properly packaged)
- Removed redundant open_channel_flow.py (replaced by open_channel_flow_new.py)
- Added sympy dependency for pydraulics symbolic math support
- Updated CI/CD to use proper pytest instead of print-based tests

### Changed
- Improved package structure for better deployment reliability
- Enhanced testing with code coverage reporting
- Updated documentation with cross-platform examples

## [1.1.0] - 2025-01-31

### Added
- Comprehensive test suite for wastewater treatment scenarios
- Detailed documentation and API reference
- Windows batch script for easy testing
- Debug and diagnostic utilities

### Fixed
- FluidProperties gamma attribute access error - now safely checks for attribute existence
- Molecular weight unit conversion issues causing 1000x errors in density calculations
- Incorrect Kv to Cv conversion factor (was ÷0.865, now ×1.156)
- Wrong units passed to fluids library for valve sizing calculations
- Property lookup failures for common fluids like Methane - now tries direct lookup first

### Changed
- Improved fluid property lookup logic to avoid unnecessary fallbacks
- Enhanced error messages and logging for better debugging
- Updated all gas tools to use consistent unit handling

### Technical Details
- MW conversion: Added ×1000 factor to convert kg/mol to kg/kmol for ideal gas calculations
- Flow rate: Added ÷3600 conversion for m³/hr to m³/s in valve sizing
- Hardcoded MW values: Updated from g/mol to kg/kmol for consistency
- Property lookup: Changed to try FluidProperties directly before FLUID_SELECTION validation

## [1.0.0] - 2025-01-30

### Initial Release
- Gas flow calculations: pipe pressure drop, blower sizing, control valve sizing
- Liquid flow calculations: pipe pressure drop, pump requirements, control valve sizing
- Fluid property lookup with CoolProp integration
- Reynolds number and flow regime determination
- Pipe property database (ASME schedules)
- Open channel flow calculations
- Support for common engineering units
- Flexible unit conversion system