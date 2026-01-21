# Changelog

All notable changes to the Fluids MCP Server will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.2.2] - 2025-01-21

### Fixed
- **Critical**: Gas orifice analysis functions now correctly derive pressure for property lookup
  - `calculate_orifice_analysis_dp` and `calculate_orifice_analysis_flow` now derive `pressure_bar` from `inlet_pressure` for gas phase
  - Previously used default 1 bar, causing ~10x density errors at elevated pressures
- **High**: `CachedFluidProperties` now queries real gas properties from CoolProp
  - Z-factor queried via `PropsSI('Z', ...)` instead of hard-coded 1.0
  - Cv queried via `PropsSI('O', ...)` instead of ideal gas approximation (Cp - R)
  - Gamma now calculated from actual Cp/Cv values for real gas accuracy
- **High**: ISO 5167 gas expansibility validity check added to all orifice modes
  - Warning issued when P2/P1 < 0.80 (outside ISO 5167 correlation validity)
  - Applied to sizing, analysis_dp, analysis_flow, and plate_kit modes
- **Medium**: `isothermal_darcy` method documentation corrected
  - Fixed misleading comments about "internal iteration" - upstream does NOT iterate
  - Added warning when ΔP/P1 > 20% (reduced accuracy per upstream docs)
  - Fixed warning propagation bug ensuring late warnings are included in results

### Added
- 8 regression tests for bug fixes with proper assertions:
  - Gas pressure ratio warnings (P2/P1 < 0.80)
  - Gas orifice pressure derivation verification
  - CachedFluidProperties Z-factor and gamma validation

## [2.2.1] - 2025-01-20

### Fixed
- Gas pipe flow `pipe_length_ft` parameter now correctly forwarded from omnitool
  - Added `pipe_length_ft` parameter to `calculate_gas_pipe_pressure_drop()` signature
  - Added ft-to-m conversion logic using `FT_to_M` constant
  - Previously, `pipe_flow(phase="gas", pipe_length_ft=...)` would fail with "Incorrect number of primary variables" error because `pipe_length_ft` was filtered out by signature inspection

## [2.2.0] - 2025-12-22

### Added
- **orifice_restrictor**: New omnitool for fixed orifice (restriction plate) sizing per ISO 5167
  - Sizing mode: Calculate orifice diameter from flow rate and target pressure drop
  - Analysis ΔP mode: Calculate pressure drop for given orifice and flow
  - Analysis flow mode: Calculate flow rate for given orifice and pressure drop
  - Plate kit mode: Generate sensitivity table for commissioning flexibility
- Reader-Harris-Gallagher discharge coefficient correlation (ISO 5167-2)
- ISO 5167 validity checking with warnings for out-of-range conditions
- Support for both liquid (incompressible) and gas (compressible) phases
- Comprehensive test suite with 20 tests covering all modes and edge cases

### Technical Details
- Uses `differential_pressure_meter_solver()` from fluids library
- Beta ratio limits: 0.1 ≤ β ≤ 0.75 per ISO 5167
- Automatic fluid property lookup via FluidProp/CoolProp
- Unit flexibility: mm, inch, or m for orifice diameter output

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