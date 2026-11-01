"""
Enhanced MCP Server for fluid dynamics calculations.

This server provides hydraulics calculation tools for a variety of fluid mechanics
problems including pipe flow, pump sizing, open channel flow, gas flow, and more.

This is a modular version where each tool is in its own module for easier maintenance
and future development.

Author: Claude AI
"""
import logging
from fastmcp import FastMCP
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('fluids-mcp')
mcp = FastMCP('fluids-calculator')
from fluids_mcp.omnitools.pipe_flow import pipe_flow
from fluids_mcp.omnitools.control_valve import control_valve
from fluids_mcp.omnitools.pipe_sizing import pipe_sizing
from fluids_mcp.omnitools.parameter_sweep import parameter_sweep
from fluids_mcp.omnitools.properties import properties
from fluids_mcp.omnitools.machine_requirements import machine_requirements
from fluids_mcp.omnitools.orifice_restrictor import orifice_restrictor
from fluids_mcp.omnitools.help_resources import help_resources
mcp.tool()(pipe_flow)
mcp.tool()(control_valve)
mcp.tool()(pipe_sizing)
mcp.tool()(parameter_sweep)
mcp.tool()(properties)
mcp.tool()(machine_requirements)
mcp.tool()(orifice_restrictor)
mcp.tool()(help_resources)
from fluids_mcp.utils.import_helpers import FLUIDPROP_AVAILABLE, COOLPROP_AVAILABLE

def main():
    mcp.run(transport='stdio')
if __name__ == '__main__':
    main()
