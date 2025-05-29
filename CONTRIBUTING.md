# Contributing to Fluids MCP Server

Thank you for your interest in contributing to the Fluids MCP Server! This document provides guidelines and instructions for contributing.

## Code of Conduct

By participating in this project, you agree to maintain a respectful and inclusive environment for all contributors.

## How to Contribute

### Reporting Issues

1. Check if the issue already exists in the [issue tracker](https://github.com/puran-water/fluids-mcp/issues)
2. If not, create a new issue with:
   - Clear, descriptive title
   - Detailed description of the problem
   - Steps to reproduce
   - Expected vs actual behavior
   - System information (OS, Python version, etc.)

### Suggesting Features

1. Open an issue with the "enhancement" label
2. Describe the feature and its use case
3. Provide examples of how it would work
4. Explain why this would be valuable

### Code Contributions

#### Setup Development Environment

1. Fork the repository
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR-USERNAME/fluids-mcp.git
   cd fluids-mcp
   ```

3. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # or venv\Scripts\activate on Windows
   ```

4. Install in development mode:
   ```bash
   pip install -r requirements.txt
   pip install -e .
   ```

#### Development Workflow

1. Create a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes following our coding standards

3. Add tests for new functionality

4. Run the test suite:
   ```bash
   python test_fluids_mcp.py
   python test_fixes.py
   ```

5. Commit your changes:
   ```bash
   git add .
   git commit -m "Add clear description of changes"
   ```

6. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

7. Create a Pull Request

## Coding Standards

### Python Style

- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/)
- Use meaningful variable and function names
- Add type hints where appropriate
- Maximum line length: 100 characters

### Documentation

- Add docstrings to all functions and classes
- Use Google-style docstrings:
  ```python
  def calculate_pressure_drop(flow_rate: float, pipe_diameter: float) -> float:
      """Calculate pressure drop in a pipe.
      
      Args:
          flow_rate: Volumetric flow rate in m³/s
          pipe_diameter: Internal pipe diameter in meters
          
      Returns:
          Pressure drop in Pascals
          
      Raises:
          ValueError: If flow_rate or pipe_diameter is negative
      """
  ```

### Testing

- Write tests for all new functionality
- Maintain or improve code coverage
- Test edge cases and error conditions
- Use descriptive test names

### Commit Messages

- Use clear, descriptive commit messages
- Start with a verb in present tense
- Keep the first line under 50 characters
- Add detailed description if needed

Examples:
```
Add support for Hazen-Williams equation
Fix unit conversion in pump power calculation
Update documentation for valve sizing methods
```

## Code Structure

### Adding New Tools

1. Create a new file in the `tools/` directory
2. Follow the existing pattern:
   ```python
   import json
   from typing import Optional
   
   def calculate_new_feature(
       param1: Optional[float] = None,
       param2: Optional[float] = None,
       # ... parameters
   ) -> str:
       """Tool description and parameter documentation."""
       
       # Input validation
       # Unit conversion
       # Calculations
       # Return JSON result
   ```

3. Add the tool to `server.py`
4. Update documentation

### Adding Unit Conversions

Add new conversion factors to `utils/constants.py`:
```python
NEW_UNIT_to_SI = 1.234  # Conversion factor with comment
```

### Error Handling

- Validate all inputs
- Provide clear error messages
- Handle edge cases gracefully
- Log warnings for unusual conditions

## Pull Request Process

1. **Title**: Clear description of changes
2. **Description**: 
   - What changes were made
   - Why they were needed
   - How they were tested
3. **Checklist**:
   - [ ] Tests pass
   - [ ] Documentation updated
   - [ ] Follows coding standards
   - [ ] No breaking changes (or documented)
4. **Review**: Address feedback promptly

## Release Process

1. Update version in relevant files
2. Update CHANGELOG.md
3. Create release notes
4. Tag the release

## Questions?

Feel free to open an issue for any questions about contributing!