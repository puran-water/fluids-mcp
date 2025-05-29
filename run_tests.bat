@echo off
echo Running Fluids MCP Server Tests
echo ================================

REM Activate virtual environment if it exists
if exist "..\..\venv\Scripts\activate.bat" (
    call "..\..\venv\Scripts\activate.bat"
) else if exist "..\venv\Scripts\activate.bat" (
    call "..\venv\Scripts\activate.bat"
) else if exist "venv\Scripts\activate.bat" (
    call "venv\Scripts\activate.bat"
)

REM Run tests
echo.
echo Running unit tests...
python -m pytest tests/ -v

echo.
echo Running integration tests...
python tests/test_fluids_mcp.py

echo.
echo Running fix validation tests...
python tests/test_fixes.py

echo.
echo Tests complete!
pause