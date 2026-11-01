"""Pydantic contracts for the eight fluids-mcp omnitools."""

from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional

from pydantic import ConfigDict, Field, field_validator, model_validator

from fluids_mcp.resource_limits import MAX_NPS_CANDIDATES, MAX_PLATE_KIT_SIZES

from fluids_mcp.errors import redact_public_payload, redact_public_text
from puran_mcp_support.pydantic_base import StrictFiniteBaseModel

MAX_SWEEP_POINTS = 100


class StrictModel(StrictFiniteBaseModel):
    """Base model with stable public-contract defaults."""

    model_config = ConfigDict(
        allow_inf_nan=False,
        extra="forbid",
        populate_by_name=True,
        validate_assignment=True,
    )


class CalculationMessage(StrictModel):
    code: str
    message: str
    severity: Literal["info", "warning", "error"] = "warning"
    details: Dict[str, Any] = Field(default_factory=dict)
    correlation_id: Optional[str] = None

    @field_validator("message", mode="before")
    @classmethod
    def redact_message(cls, value: Any) -> Any:
        return redact_public_text(value)

    @field_validator("details", mode="before")
    @classmethod
    def redact_details(cls, value: Any) -> Any:
        return redact_public_payload(value or {})


class OutputBase(StrictModel):
    ok: bool = True
    errors: List[CalculationMessage] = Field(default_factory=list)
    warnings: List[CalculationMessage] = Field(default_factory=list)
    inputs_resolved: List[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def align_ok_with_errors(self) -> "OutputBase":
        if self.errors:
            object.__setattr__(self, "ok", False)
        return self


class FluidPropertyFields(StrictModel):
    fluid_name: Optional[str] = None
    temperature_c: Optional[float] = None
    pressure_bar: Optional[float] = None
    fluid_density: Optional[float] = None
    fluid_density_lbft3: Optional[float] = None
    fluid_viscosity: Optional[float] = None
    fluid_viscosity_cp: Optional[float] = None


class GasPropertyFields(StrictModel):
    gas_mw: Optional[float] = None
    gas_gamma: Optional[float] = None
    gas_z_factor: Optional[float] = None
    gas_viscosity: Optional[float] = None
    gas_composition_mol: Optional[Dict[str, float]] = None
    allow_property_defaults: bool = False


class PipeGeometryFields(StrictModel):
    pipe_diameter: Optional[float] = None
    pipe_diameter_in: Optional[float] = None
    nominal_size_in: Optional[float] = None
    schedule: str = "40"
    material: Optional[str] = None
    pipe_roughness: Optional[float] = None
    fittings: Optional[List[Dict[str, Any]]] = None


class PipeFlowInput(FluidPropertyFields, GasPropertyFields, PipeGeometryFields):
    phase: Literal["liquid", "gas"] = "liquid"
    pipe_length: Optional[float] = None
    pipe_length_ft: Optional[float] = None
    flow_rate: Optional[float] = None
    flow_rate_gpm: Optional[float] = None
    pressure_drop: Optional[float] = None
    pressure_drop_psi: Optional[float] = None
    inlet_pressure: Optional[float] = None
    outlet_pressure: Optional[float] = None
    flow_rate_kg_s: Optional[float] = None
    flow_rate_norm_m3_hr: Optional[float] = None
    flow_rate_std_m3_hr: Optional[float] = None
    method: Literal[
        "Weymouth",
        "Panhandle_A",
        "Panhandle_B",
        "IGT",
        "Oliphant",
        "Spitzglass_low",
        "Spitzglass_high",
        "isothermal_darcy",
    ] = "Weymouth"


class PipeFlowOutput(OutputBase):
    phase: Optional[Literal["liquid", "gas"]] = Field(...)
    solved_variable: Optional[str] = None
    flow_rate_m3s: Optional[float] = None
    flow_rate_gpm: Optional[float] = None
    flow_rate_kg_s: Optional[float] = None
    flow_rate_norm_m3_hr: Optional[float] = None
    inlet_pressure_pa: Optional[float] = None
    outlet_pressure_pa: Optional[float] = None
    pressure_drop_pa: Optional[float] = None
    pressure_drop_total_pa: Optional[float] = None
    pipe_diameter_m: Optional[float] = None
    pipe_length_m: Optional[float] = None
    reynolds_number: Optional[float] = None
    friction_factor: Optional[float] = None
    flow_velocity_m_s: Optional[float] = None
    average_velocity_m_s: Optional[float] = None
    average_mach_number: Optional[float] = None
    head_loss_m: Optional[float] = None
    head_loss_ft: Optional[float] = None
    choked: Optional[bool] = None
    pressure_ratio: Optional[float] = None
    critical_pressure_ratio: Optional[float] = None
    expansion_factor: Optional[float] = None
    details: Dict[str, Any] = Field(default_factory=dict)


class ControlValveInput(FluidPropertyFields, GasPropertyFields):
    phase: Literal["liquid", "gas"] = "liquid"
    flow_rate: Optional[float] = None
    flow_rate_gpm: Optional[float] = None
    flow_rate_kg_s: Optional[float] = None
    flow_rate_norm_m3_hr: Optional[float] = None
    flow_rate_std_m3_hr: Optional[float] = None
    inlet_pressure: Optional[float] = None
    outlet_pressure: Optional[float] = None
    inlet_pressure_psi: Optional[float] = None
    outlet_pressure_psi: Optional[float] = None
    fluid_saturation_pressure_pa: Optional[float] = None
    fluid_saturation_pressure_psi: Optional[float] = None
    fluid_critical_pressure_pa: Optional[float] = None
    fluid_critical_pressure_psi: Optional[float] = None
    valve_fl: Optional[float] = None
    valve_fd: Optional[float] = None
    inlet_temperature_c: Optional[float] = None
    valve_xt: Optional[float] = None
    valve_type: Literal["globe", "butterfly", "ball", "gate", "other"] = "globe"


class ControlValveOutput(OutputBase):
    phase: Optional[Literal["liquid", "gas"]] = Field(...)
    flow_rate_kg_s: Optional[float] = None
    inlet_pressure_pa: Optional[float] = None
    outlet_pressure_pa: Optional[float] = None
    pressure_drop_pa: Optional[float] = None
    recommended_nps_inches: Optional[float] = None
    required_cv: Optional[float] = None
    rated_cv: Optional[float] = None
    cv_utilization: Optional[float] = None
    cavitation: Optional[bool] = None
    flashing: Optional[bool] = None
    choked: Optional[bool] = None
    expansion_factor: Optional[float] = None
    pressure_ratio: Optional[float] = None
    critical_pressure_ratio: Optional[float] = None
    mach_outlet: Optional[float] = None
    details: Dict[str, Any] = Field(default_factory=dict)


class PipeSizingInput(FluidPropertyFields, GasPropertyFields):
    phase: Literal["liquid", "gas"] = "liquid"
    flow_rate: Optional[float] = None
    flow_rate_gpm: Optional[float] = None
    flow_rate_kg_s: Optional[float] = None
    flow_rate_norm_m3_hr: Optional[float] = None
    flow_rate_std_m3_hr: Optional[float] = None
    allowable_dp_pa: Optional[float] = None
    allowable_dp_psi: Optional[float] = None
    inlet_pressure: Optional[float] = None
    inlet_pressure_psi: Optional[float] = None
    pipe_length: Optional[float] = None
    pipe_length_ft: Optional[float] = None
    schedule: str = "40"
    nps_candidates: Optional[List[float]] = Field(default=None, min_length=1, max_length=MAX_NPS_CANDIDATES)
    material: Optional[str] = None
    pipe_roughness: Optional[float] = None
    fittings: Optional[List[Dict[str, Any]]] = None
    method: str = "isothermal_darcy"
    velocity_min_m_s: Optional[float] = None
    velocity_max_m_s: Optional[float] = None
    return_all_candidates: bool = False


class PipeSizingOutput(OutputBase):
    phase: Optional[Literal["liquid", "gas"]] = Field(...)
    recommended_nps_in: Optional[float] = None
    recommended_details: Dict[str, Any] = Field(default_factory=dict)
    feasible_candidates: List[Dict[str, Any]] = Field(default_factory=list)
    all_candidates: List[Dict[str, Any]] = Field(default_factory=list)
    details: Dict[str, Any] = Field(default_factory=dict)


class SweepRow(StrictModel):
    ok: bool
    input_value: Any
    result: Dict[str, Any] = Field(default_factory=dict)
    error_code: Optional[str] = None
    message: Optional[str] = None
    error: Optional[CalculationMessage] = None


class ParameterSweepInput(FluidPropertyFields, GasPropertyFields):
    calculation: Literal["pipe_liquid", "pipe_gas", "blower"]
    variable: str
    start: float
    stop: float
    n: int = Field(ge=1, le=MAX_SWEEP_POINTS, description="Number of sweep points, from 1 to 100.")
    flow_rate: Optional[float] = None
    flow_rate_gpm: Optional[float] = None
    flow_rate_kg_s: Optional[float] = None
    flow_rate_norm_m3_hr: Optional[float] = None
    flow_rate_std_m3_hr: Optional[float] = None
    pipe_diameter: Optional[float] = None
    pipe_diameter_in: Optional[float] = None
    pipe_length: Optional[float] = None
    pipe_length_ft: Optional[float] = None
    nominal_size_in: Optional[float] = None
    schedule: str = "40"
    material: Optional[str] = None
    pipe_roughness: Optional[float] = None
    inlet_pressure: Optional[float] = None
    outlet_pressure: Optional[float] = None
    method: str = "Weymouth"
    inlet_pressure_psi: Optional[float] = None
    outlet_pressure_psi: Optional[float] = None
    inlet_temperature_c: Optional[float] = None
    efficiency: Optional[float] = None
    efficiency_type: str = "isentropic"
    use_eos_calculations: Optional[bool] = None
    solve_for: Optional[str] = None

    @field_validator("n")
    @classmethod
    def positive_n(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("n must be positive")
        return value


class ParameterSweepOutput(OutputBase):
    sweep_variable: str
    sweep_range: Dict[str, Any]
    rows: List[SweepRow] = Field(default_factory=list)
    summary: Dict[str, int] = Field(default_factory=dict)


class PropertiesInput(StrictModel):
    lookup_type: Literal["fluid", "pipe", "list_fluids"] = "fluid"
    fluid_name: Optional[str] = None
    temperature_c: Optional[float] = None
    pressure_bar: Optional[float] = 1.0
    nominal_size: Optional[float] = None
    schedule: str = "40"
    inner_diameter: Optional[float] = None
    outer_diameter: Optional[float] = None
    material: Optional[str] = None


class PropertiesOutput(OutputBase):
    lookup_type: Optional[Literal["fluid", "pipe", "list_fluids"]] = Field(...)
    properties: Dict[str, Any] = Field(default_factory=dict)
    resources: Dict[str, Any] = Field(default_factory=dict)


class MachineRequirementsInput(FluidPropertyFields, GasPropertyFields):
    machine_type: Literal["pump", "compressor", "reynolds", "open_channel"] = "pump"
    pump_mode: Literal["flow_rate_given", "target_tdh_solve_flow", "target_tdh_solve_diameter"] = "flow_rate_given"
    endpoint_basis: Literal["tank_surface", "nozzle"] = "tank_surface"
    compressor_mode: Literal["specified_outlet_pressure", "aeration_system_pressure"] = "specified_outlet_pressure"
    compressor_type: Literal["turbo", "roots_pd"] = "turbo"
    open_channel_mode: Literal["normal_flow", "gravity_leg"] = "normal_flow"
    flow_rate: Optional[float] = None
    flow_rate_gpm: Optional[float] = None
    flow_rate_kg_s: Optional[float] = None
    flow_rate_norm_m3_hr: Optional[float] = None
    flow_rate_std_m3_hr: Optional[float] = None
    flow_rate_m3s: Optional[float] = None
    pipe_diameter: Optional[float] = None
    pipe_diameter_in: Optional[float] = None
    suction_pipe_diameter: Optional[float] = None
    suction_pipe_diameter_in: Optional[float] = None
    discharge_pipe_diameter: Optional[float] = None
    discharge_pipe_diameter_in: Optional[float] = None
    pipe_length: Optional[float] = None
    pipe_length_ft: Optional[float] = None
    nominal_size_in: Optional[float] = None
    schedule: str = "40"
    material: Optional[str] = None
    pipe_roughness: Optional[float] = None
    suction_material: Optional[str] = None
    discharge_material: Optional[str] = None
    suction_pipe_roughness: Optional[float] = None
    discharge_pipe_roughness: Optional[float] = None
    fittings: Optional[List[Dict[str, Any]]] = None
    suction_pipe_length: Optional[float] = None
    suction_pipe_length_ft: Optional[float] = None
    discharge_pipe_length: Optional[float] = None
    discharge_pipe_length_ft: Optional[float] = None
    static_suction_head: Optional[float] = None
    static_suction_head_ft: Optional[float] = None
    static_discharge_head: Optional[float] = None
    static_discharge_head_ft: Optional[float] = None
    suction_fittings: Optional[List[Dict[str, Any]]] = None
    discharge_fittings: Optional[List[Dict[str, Any]]] = None
    suction_nozzle_diameter: Optional[float] = None
    discharge_nozzle_diameter: Optional[float] = None
    suction_nozzle_diameter_in: Optional[float] = None
    discharge_nozzle_diameter_in: Optional[float] = None
    atmospheric_pressure: Optional[float] = None
    atmospheric_pressure_psi: Optional[float] = None
    suction_tank_pressure: Optional[float] = None
    suction_tank_pressure_psi: Optional[float] = None
    discharge_tank_pressure: Optional[float] = None
    discharge_tank_pressure_psi: Optional[float] = None
    suction_nozzle_pressure: Optional[float] = None
    suction_nozzle_pressure_psi: Optional[float] = None
    discharge_nozzle_pressure: Optional[float] = None
    discharge_nozzle_pressure_psi: Optional[float] = None
    fluid_vapor_pressure: Optional[float] = None
    fluid_vapor_pressure_psi: Optional[float] = None
    target_tdh: Optional[float] = None
    target_tdh_ft: Optional[float] = None
    npshr: Optional[float] = None
    npshr_ft: Optional[float] = None
    required_npsh_margin: Optional[float] = None
    required_npsh_margin_ft: Optional[float] = None
    inlet_pressure: Optional[float] = None
    outlet_pressure: Optional[float] = None
    inlet_pressure_psi: Optional[float] = None
    outlet_pressure_psi: Optional[float] = None
    inlet_temperature_c: Optional[float] = None
    efficiency: Optional[float] = None
    efficiency_type: str = "isentropic"
    polytropic_n: Optional[float] = None
    submergence_m: Optional[float] = None
    submergence_liquid_density_kg_m3: float = 998.2
    diffuser_dp_clean_mbar: Optional[float] = None
    diffuser_dp_fouled_mbar: Optional[float] = None
    use_eos_calculations: Optional[bool] = None
    solve_for: Optional[str] = None
    velocity: Optional[float] = None
    velocity_ft_s: Optional[float] = None
    characteristic_length: Optional[float] = None
    characteristic_length_in: Optional[float] = None
    fluid_kinematic_viscosity: Optional[float] = None
    channel_type: Optional[str] = None
    width_m: Optional[float] = None
    diameter_m: Optional[float] = None
    side_slope_z: Optional[float] = None
    slope: Optional[float] = None
    manning_n: Optional[float] = None
    normal_depth_m: Optional[float] = None
    circular_depth_branch: Literal["shallow", "near_full"] = "shallow"
    invert_elevation_m: Optional[float] = None
    upstream_head_m: Optional[float] = None
    downstream_head_m: Optional[float] = None


class MachineRequirementsOutput(OutputBase):
    machine_type: Optional[Literal["pump", "compressor", "reynolds", "open_channel"]] = Field(...)
    pump_mode: Optional[Literal["flow_rate_given", "target_tdh_solve_flow", "target_tdh_solve_diameter"]] = None
    compressor_mode: Optional[Literal["specified_outlet_pressure", "aeration_system_pressure"]] = None
    compressor_type: Optional[Literal["turbo", "roots_pd"]] = None
    open_channel_mode: Optional[Literal["normal_flow", "gravity_leg"]] = None
    flow_rate_m3s: Optional[float] = None
    flow_rate_gpm: Optional[float] = None
    pipe_diameter_m: Optional[float] = None
    pipe_velocity_m_s: Optional[float] = None
    static_head_m: Optional[float] = None
    static_head_ft: Optional[float] = None
    target_tdh_m: Optional[float] = None
    target_tdh_ft: Optional[float] = None
    total_dynamic_head_m: Optional[float] = None
    total_dynamic_head_ft: Optional[float] = None
    npsha_m: Optional[float] = None
    npsha_ft: Optional[float] = None
    npshr_m: Optional[float] = None
    npshr_ft: Optional[float] = None
    npshr_status: Optional[str] = None
    npsh_margin_m: Optional[float] = None
    npsh_margin_ft: Optional[float] = None
    required_npsh_margin_m: Optional[float] = None
    npsh_margin_pass: Optional[bool] = None
    estimated_power_kw: Optional[float] = None
    estimated_power_hp: Optional[float] = None
    required_power_kw: Optional[float] = None
    shaft_power_kw: Optional[float] = None
    required_discharge_pressure_gauge_pa: Optional[float] = None
    required_discharge_pressure_absolute_pa: Optional[float] = None
    pressure_component_breakdown: Optional[Dict[str, Any]] = None
    blower_scenarios: Optional[Dict[str, Any]] = None
    actual_discharge_temp_c: Optional[float] = None
    reynolds_number: Optional[float] = None
    normal_depth_m: Optional[float] = None
    critical_depth_m: Optional[float] = None
    specific_energy_m: Optional[float] = None
    critical_specific_energy_m: Optional[float] = None
    froude_number: Optional[float] = None
    flow_regime: Optional[str] = None
    flow_state: Optional[str] = None
    hgl_elevation_m: Optional[float] = None
    hydraulic_grade_above_invert_m: Optional[float] = None
    energy_grade_above_invert_m: Optional[float] = None
    energy_grade_elevation_m: Optional[float] = None
    available_fall_m: Optional[float] = None
    required_fall_m: Optional[float] = None
    residual_fall_m: Optional[float] = None
    gravity_capacity_pass: Optional[bool] = None
    details: Dict[str, Any] = Field(default_factory=dict)


class OrificeRestrictorInput(FluidPropertyFields):
    mode: Literal["sizing", "analysis_dp", "analysis_flow", "plate_kit"] = "sizing"
    phase: Literal["liquid", "gas"] = "liquid"
    flow_rate: Optional[float] = None
    flow_rate_gpm: Optional[float] = None
    flow_rate_m3_hr: Optional[float] = None
    pressure_drop: Optional[float] = None
    pressure_drop_psi: Optional[float] = None
    pressure_drop_bar: Optional[float] = None
    inlet_pressure: Optional[float] = None
    inlet_pressure_psi: Optional[float] = None
    inlet_pressure_bar: Optional[float] = None
    orifice_diameter: Optional[float] = None
    orifice_diameter_in: Optional[float] = None
    orifice_diameter_mm: Optional[float] = None
    pipe_diameter: Optional[float] = None
    pipe_diameter_in: Optional[float] = None
    nominal_size_in: Optional[float] = None
    schedule: str = "40"
    gas_gamma: Optional[float] = None
    gas_mw: Optional[float] = None
    gas_z_factor: Optional[float] = None
    cd: Optional[float] = None
    meter_type: Optional[str] = None
    taps: str = "D"
    # A restriction plate is sized against the loss the system keeps, not the
    # tap-to-tap differential a meter reads; dp_basis picks which one the supplied
    # pressure_drop_* means. duty relaxes nothing - it annotates the ISO 5167
    # metering-envelope warnings for a device legitimately outside them.
    dp_basis: Literal["tap", "permanent"] = "tap"
    duty: Literal["metering", "restriction"] = "metering"
    # Kept separate on purpose: a gas-saturated liquid breaks out at its bubble
    # point, which can be an order of magnitude above the water vapour pressure.
    vapor_pressure_pa: Optional[float] = None
    bubble_point_pa: Optional[float] = None
    diameter_variability_percent: float = 20.0
    kit_sizes: int = Field(default=5, strict=True, ge=2, le=MAX_PLATE_KIT_SIZES)
    size_units: Literal["m", "inch", "mm"] = "mm"


class OrificeRestrictorOutput(OutputBase):
    mode: Optional[Literal["sizing", "analysis_dp", "analysis_flow", "plate_kit"]] = Field(...)
    phase: Optional[Literal["liquid", "gas"]] = Field(...)
    orifice_diameter_m: Optional[float] = None
    orifice_diameter_mm: Optional[float] = None
    orifice_diameter_in: Optional[float] = None
    beta_ratio: Optional[float] = None
    pressure_drop_pa: Optional[float] = None
    permanent_pressure_loss_pa: Optional[float] = None
    permanent_pressure_loss_psi: Optional[float] = None
    permanent_pressure_loss_bar: Optional[float] = None
    pressure_recovery_fraction: Optional[float] = None
    discharge_coefficient: Optional[float] = None
    cavitation: Optional[Dict[str, Any]] = None
    flow_rate_m3_hr: Optional[float] = None
    mass_flow_kg_s: Optional[float] = None
    choked: Optional[bool] = None
    plates: List[Dict[str, Any]] = Field(default_factory=list)
    details: Dict[str, Any] = Field(default_factory=dict)


class HelpResourcesInput(StrictModel):
    resource_type: Literal["fittings", "fluids", "materials", "methods", "all"] = "all"
    category: Optional[str] = None


class HelpResourcesOutput(OutputBase):
    resource_type: Literal["fittings", "fluids", "materials", "methods", "all"] = Field(...)
    resources: Dict[str, Any] = Field(default_factory=dict)
