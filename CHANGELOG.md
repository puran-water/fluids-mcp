# Changelog

All notable changes to the Fluids MCP Server will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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