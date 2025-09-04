Units and Parameter Conventions

This project uses a simple, explicit convention for units:

- No suffix = SI units
- With suffix = Non‑SI variant or alternative basis

Unless noted, numeric inputs and outputs are scalar floats in base SI.

Core Conventions

- Pressure: `Pa` (absolute). Alternatives include `pressure_drop_psi` (psi), `pressure_bar` (bar).
- Temperature: `temperature_c` (°C). Internally converted to Kelvin as needed.
- Length: `m`. Alternatives include `pipe_length_ft` (ft) and `pipe_diameter_in`/`nominal_size_in` (in).
- Diameter: `m` when `pipe_diameter`; inches when `pipe_diameter_in` or `nominal_size_in`.
- Roughness: `m` (absolute roughness).
- Mass flow: `kg/s` (`flow_rate_kg_s`).
- Volumetric flow (gas):
  - `flow_rate_norm_m3_hr`: m³/hr at Normal conditions (0°C, 1 atm).
  - `flow_rate_std_m3_hr`: m³/hr at Standard conditions (15°C, 1 atm).
  - Internally, some gas methods use m³/s at base conditions; conversions are handled inside tools.
- Volumetric flow (liquid): `flow_rate` in m³/s; alternative `flow_rate_gpm` in US GPM.
- Density: `kg/m³` (`fluid_density`, results). Alternative input `fluid_density_lbft3` (lb/ft³).
- Viscosity: dynamic viscosity `Pa·s` (`fluid_viscosity`, results). Alternative input `fluid_viscosity_cp` (centipoise).
- Molecular weight: `kg/kmol` (numerically equal to g/mol). Parameter `gas_mw` accepts either interpretation.
- Specific heat ratio (gamma): dimensionless (`gas_gamma`).
- Compressibility factor (Z): dimensionless (`gas_z_factor`).
- Pipe material: string name; roughness is looked up when possible.
- Gas composition: `gas_composition_mol` as `{component: mole_fraction}` with fraction sum ≈ 1.0.

Common Parameters by Omnitool

Pipe Flow (`omnitools/pipe_flow.py` → liquid or gas):

- Common:
  - `pipe_length` (m), `pipe_length_ft` (ft)
  - `pipe_diameter` (m), `pipe_diameter_in` (in), `nominal_size_in` (in), `schedule` (string)
  - `material` (string), `pipe_roughness` (m), `fittings` (list)
  - `temperature_c` (°C), `fluid_name` (string)
- Liquid phase:
  - `flow_rate` (m³/s), `flow_rate_gpm` (US GPM)
  - `fluid_density` (kg/m³), `fluid_density_lbft3` (lb/ft³)
  - `fluid_viscosity` (Pa·s), `fluid_viscosity_cp` (cP)
  - `pressure_drop` (Pa), `pressure_drop_psi` (psi), `pressure_bar` (bar)
- Gas phase:
  - `inlet_pressure`, `outlet_pressure` (Pa, absolute)
  - `flow_rate_kg_s` (kg/s), `flow_rate_norm_m3_hr` (m³/hr @ 0°C, 1 atm), `flow_rate_std_m3_hr` (m³/hr @ 15°C, 1 atm)
  - `gas_mw` (kg/kmol = g/mol), `gas_gamma` (-), `gas_z_factor` (-), `gas_viscosity` (Pa·s)
  - `gas_composition_mol` ({component: mole_fraction})
  - `method` (Weymouth, Panhandle_A, Panhandle_B, IGT, Oliphant, Spitzglass_low, Spitzglass_high, isothermal_darcy)

Control Valve (`omnitools/control_valve.py` → liquid or gas):

- Liquid: uses `flow_rate` (m³/s) or `flow_rate_gpm`, `pressure_drop` (Pa) or `pressure_drop_psi`, `fluid_density` (kg/m³) or lookup via `fluid_name`+`temperature_c`+`pressure_bar`.
- Gas: uses `flow_rate_norm_m3_hr`/`flow_rate_std_m3_hr`/`flow_rate_kg_s`, `inlet_pressure`/`outlet_pressure` (Pa), and gas properties (`gas_mw`, `gas_gamma`, `gas_z_factor`, `gas_viscosity`) or `gas_composition_mol`.
- Result `Cv` is typically dimensionless in US customary context; internal calculations respect SI with appropriate conversions.

Pipe Sizing (`omnitools/pipe_sizing.py`):

- Inputs mirror Pipe Flow; optimizes for diameter (`pipe_diameter` in m) under constraints. Gas/liquid options share unit conventions as above.

Parameter Sweep (`omnitools/parameter_sweep.py`):

- Sweeps one parameter in the same units defined above. Example: `inlet_pressure` in Pa, `pipe_diameter` in m, `flow_rate_norm_m3_hr` in m³/hr.

Properties (`omnitools/properties.py`):

- Fluid lookup: `fluid_name` (string), `temperature_c` (°C), `pressure_bar` (bar). Returns density (kg/m³), viscosity (Pa·s), Cp/Cv (J/kg/K), etc.
- Pipe lookup: `nominal_size` (in), `schedule` (string), `inner_diameter` (m), `outer_diameter` (m), `material` (string).

Machine Requirements (`omnitools/machine_requirements.py`):

- Pump/compressor/hydraulics; follows the same unit conventions for flow, pressure, temperature, and properties.

Standard Conditions and Constants

- Normal: `P_NORM = 101325 Pa`, `T_NORM = 273.15 K` (0°C)
- Standard: `P_STD = 101325 Pa`, `T_STD = 288.15 K` (15°C)
- Gas constant: `R_UNIV = 8314.462 J/(kmol·K)`

Conversions Used (utils/constants.py)

- `GPM_to_M3S = 0.0000630902` (US GPM → m³/s)
- `INCH_to_M = 0.0254` (in → m)
- `FT_to_M = 0.3048` (ft → m)
- `PSI_to_PA = 6894.76` (psi → Pa)
- `LBFT3_to_KGM3 = 16.0185` (lb/ft³ → kg/m³)
- `CENTIPOISE_to_PAS = 0.001` (cP → Pa·s)
- `DEG_C_to_K = 273.15` (°C → K offset)

Notes and Best Practices

- Provide exactly the required number of independent variables for solvers. For gas pipe flow, specify exactly four of {`P1`, `P2`, `L`, `Q`, `D`} depending on method.
- For gas mixtures, prefer `gas_composition_mol` over a single `fluid_name` when available; MW, gamma, Z, and viscosity are computed via robust mixture rules.
- If you must input non‑SI values, use the explicit “_in”, “_ft”, “_psi”, “_gpm”, “_std_m3_hr”, “_norm_m3_hr” parameter variants. The base parameter names always assume SI.
