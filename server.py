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

# Import omnitools (consolidated tools)
from omnitools.pipe_flow import pipe_flow
from omnitools.control_valve import control_valve
from omnitools.pipe_sizing import pipe_sizing
from omnitools.parameter_sweep import parameter_sweep
from omnitools.properties import properties
from omnitools.machine_requirements import machine_requirements
from omnitools.orifice_restrictor import orifice_restrictor

# Register omnitools with MCP
mcp.tool()(pipe_flow)
mcp.tool()(control_valve)
mcp.tool()(pipe_sizing)
mcp.tool()(parameter_sweep)
mcp.tool()(properties)
mcp.tool()(machine_requirements)
mcp.tool()(orifice_restrictor)

# Log information about available dependencies
from utils.import_helpers import FLUIDPROP_AVAILABLE, COOLPROP_AVAILABLE

if __name__ == "__main__":
    logger.info("Starting Fluids MCP server...")
    logger.info("FluidProp available: %s", FLUIDPROP_AVAILABLE)
    logger.info("CoolProp available: %s", COOLPROP_AVAILABLE)
    
    if COOLPROP_AVAILABLE:
        from utils.import_helpers import get_coolprop_fluids_list
        fluids = get_coolprop_fluids_list()
        logger.info("CoolProp fluids count: %d", len(fluids))
    
    # Log which omnitools are registered
    logger.info("Registered omnitools (consolidated from 17 tools):")
    logger.info("  - pipe_flow: Unified liquid/gas pipe pressure drop calculations")
    logger.info("  - control_valve: Unified liquid/gas control valve sizing")
    logger.info("  - pipe_sizing: Unified liquid/gas pipe sizing optimization")
    logger.info("  - parameter_sweep: Unified parameter sweeps for all calculations")
    logger.info("  - properties: Unified fluid/pipe property lookups")
    logger.info("  - machine_requirements: Unified pump/compressor/hydraulic calculations")
    logger.info("  - orifice_restrictor: Fixed orifice/restriction plate sizing")
    
    # Start the server
    mcp.run()
