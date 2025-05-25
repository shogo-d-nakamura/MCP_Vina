"""
Entry point for the molDocking package
"""
import sys
import os
from pathlib import Path

# Add the parent directory to the system path
parent_dir = str(Path(__file__).resolve().parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

if __name__ == "__main__":
    try:
        # Import and run the server
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from server import mcp
        print("Starting MCP server for molecular docking simulation...", file=sys.stderr)
        mcp.run()
    except ImportError as e:
        print(f"Import Error: {e}", file=sys.stderr)
        print("Please install the required packages: uv pip install -e .", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Server Error: {e}", file=sys.stderr)
        sys.exit(1)
