"""Resolve installed public companion modules without filesystem fallbacks."""
import importlib

def ensure_sibling_importable(package_name, *args, **kwargs):
    aliases = {"process_model_schema": "puran_mcp_support.schema", "process_chemistry_reference": "puran_mcp_support.chemistry"}
    importlib.import_module(aliases.get(package_name, package_name))
