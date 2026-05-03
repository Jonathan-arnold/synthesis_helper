"""MCP server exposing the synthesis_helper API to Claude Code.

This package name (`synthesis_helper.mcp`) does not shadow the top-level
`mcp` SDK package: within server.py, `from mcp.server.fastmcp import FastMCP`
still resolves to the installed SDK via absolute imports.
"""

from synthesis_helper.mcp.server import mcp

__all__ = ["mcp"]
