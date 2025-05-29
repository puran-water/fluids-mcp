"""
Enhanced MCP Server for fluid dynamics calculations.

This server provides hydraulics calculation tools for a variety of fluid mechanics 
problems including pipe flow, pump sizing, open channel flow, gas flow, and more.

This is a modular version where each tool is in its own module for easier maintenance
and future development.

Author: Claude AI
"""

import logging
from mcp.server.fastmcp import FastMCP

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("fluids-mcp")

# Initialize the MCP server
mcp = FastMCP("fluids-calculator")

# Import all tools
from tools.pipe_pressure_drop import calculate_pipe_pressure_drop
from tools.fluid_properties import get_fluid_properties, list_available_fluids
from tools.reynolds_number import calculate_reynolds_number
from tools.pump_requirements import calculate_pump_requirements
from tools.pipe_properties import get_pipe_properties
from tools.liquid_control_valve import calculate_liquid_control_valve
from tools.gas_pipe_pressure_drop import calculate_gas_pipe_pressure_drop
from tools.blower_compressor import calculate_blower_compressor_requirements
from tools.gas_control_valve import calculate_gas_control_valve
from tools.open_channel_flow_new import calculate_open_channel_flow

# Register all tools with MCP
mcp.tool()(calculate_pipe_pressure_drop)
mcp.tool()(get_fluid_properties)
mcp.tool()(list_available_fluids)
mcp.tool()(calculate_reynolds_number)
mcp.tool()(calculate_pump_requirements)
mcp.tool()(get_pipe_properties)
mcp.tool()(calculate_liquid_control_valve)
mcp.tool()(calculate_gas_pipe_pressure_drop)
mcp.tool()(calculate_blower_compressor_requirements)
mcp.tool()(calculate_gas_control_valve)
mcp.tool()(calculate_open_channel_flow)

# Log information about available dependencies
from utils.import_helpers import FLUIDPROP_AVAILABLE, COOLPROP_AVAILABLE

if __name__ == "__main__":
    logger.info("Starting Fluids MCP server...")
    logger.info(f"FluidProp available: {FLUIDPROP_AVAILABLE}")
    logger.info(f"CoolProp available: {COOLPROP_AVAILABLE}")
    
    if COOLPROP_AVAILABLE:
        from utils.import_helpers import get_coolprop_fluids_list
        fluids = get_coolprop_fluids_list()
        logger.info(f"CoolProp fluids count: {len(fluids)}")
    
    # Log which tools are registered
    logger.info("Registered tools:")
    logger.info("  - calculate_pipe_pressure_drop")
    logger.info("  - get_fluid_properties")
    logger.info("  - list_available_fluids")
    logger.info("  - calculate_reynolds_number")
    logger.info("  - calculate_pump_requirements")
    logger.info("  - get_pipe_properties")
    logger.info("  - calculate_liquid_control_valve")
    logger.info("  - calculate_gas_pipe_pressure_drop")
    logger.info("  - calculate_blower_compressor_requirements")
    logger.info("  - calculate_gas_control_valve")
    logger.info("  - calculate_open_channel_flow")
    
    # Start the server
    mcp.run()
