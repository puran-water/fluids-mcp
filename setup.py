#!/usr/bin/env python
"""Setup script for fluids-mcp package."""

from setuptools import setup, find_packages

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="fluids-mcp",
    version="1.1.0",
    author="Puran Water LLC",
    author_email="engineering@puranwater.com",
    description="MCP server for fluid mechanics and hydraulics calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/puran-water/fluids-mcp",
    packages=find_packages(exclude=["tests", "tests.*"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    install_requires=[
        "mcp>=1.0.0",
        "httpx>=0.25.0", 
        "pydantic>=2.0.0",
        "fluids>=1.0.26",
        "scipy>=1.10.0",
        "numpy>=1.24.0",
        "CoolProp>=6.4.1",
        "fluidprop>=1.0.0",
    ],
    entry_points={
        "console_scripts": [
            "fluids-mcp=server:main",
        ],
    },
)