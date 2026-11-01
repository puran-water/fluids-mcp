# puran-mcp-support

Shared validation, numerical contracts and engineering helpers for Puran Water MCP servers. It runs no agent and supplies no private data, credentials, database connector or approval authority.

Install from the public release:

```sh
uv pip install 'puran-mcp-support @ git+https://github.com/puran-water/fluids-mcp.git@v0.1.0#subdirectory=support/puran-mcp-support'
```

The base package supplies validation. `thermal`, `chemistry` and `simulation` extras supply the optional calculation dependencies. There is no MCP entrypoint: install the appropriate server for agent access. Models and correlations require engineering review; validation is not a performance guarantee. See the repository's known issues and qualification report.

Built by [Puran Water](https://puranwater.com).
