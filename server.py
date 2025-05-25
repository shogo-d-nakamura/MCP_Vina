#!/usr/bin/env python
"""
MCP server for molecular docking simulation
"""
import logging
import sys
from typing import Any, Dict, List, Optional

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Check for required module imports
try:
    from mcp.server.fastmcp import FastMCP
    
    # Import molDocking library
    from molDocking.docking import dock_molecule, get_available_targets
    
    required_modules_available = True
except ImportError as e:
    print(f"Failed to import required modules: {str(e)}", file=sys.stderr)
    print("Please install required packages: uv pip install -e .", file=sys.stderr)
    # Implement minimal MCP server
    if 'mcp' in str(e):
        print("MCP module is not installed. Starting minimal server.", file=sys.stderr)
        try:
            # Minimal implementation using standard library only
            print("Starting minimal MCP server...", file=sys.stderr)
            import http.server
            import socketserver
            import json
            
            PORT = 8080
            
            class MinimalHandler(http.server.SimpleHTTPRequestHandler):
                def do_GET(self):
                    self.send_response(200)
                    self.send_header('Content-type', 'application/json')
                    self.end_headers()
                    response = {
                        "status": "error",
                        "message": "Required packages are not installed. Please run `uv pip install -e .`."
                    }
                    self.wfile.write(json.dumps(response).encode())
                
                def do_POST(self):
                    self.send_response(200)
                    self.send_header('Content-type', 'application/json')
                    self.end_headers()
                    response = {
                        "error": "Required packages are not installed. Please run `uv pip install -e .`."
                    }
                    self.wfile.write(json.dumps(response).encode())
            
            print(f"Starting minimal server on port {PORT}...", file=sys.stderr)
            with socketserver.TCPServer(("", PORT), MinimalHandler) as httpd:
                print(f"Server is running on port {PORT}", file=sys.stderr)
                httpd.serve_forever()
        except Exception as server_error:
            print(f"Failed to start minimal server: {str(server_error)}", file=sys.stderr)
            sys.exit(1)
    sys.exit(1)

# Initialize MCP server
mcp = FastMCP("Molecular Docking Server")


@mcp.tool()
def run_molecular_docking(
    smiles: str, 
    target: str = "akt1"
) -> Dict[str, Any]:
    """
    Dock a molecule against a target protein and calculate binding affinity
    
    Args:
        smiles: SMILES string of the molecule to dock
        target: Name of the target protein (default: 'akt1')
        
    Returns:
        Dictionary containing docking results including binding scores and energies
    """
    try:
        # Check if SMILES is empty
        if not smiles:
            return {"error": "No SMILES string provided"}
            
        logger.info(f"Running docking simulation for molecule {smiles} against target {target}")
        
        # Run docking using the core function
        results = dock_molecule(smiles, target)
        
        # Check for error
        if "error" in results:
            logger.error(f"Error in docking: {results['error']}")
            return results
            
        logger.info(f"Docking completed successfully. Best score: {results.get('best_score')}")
        
        return {
            "smiles": smiles,
            "target": results.get("target"),
            "binding_affinity": results.get("energy"),
            "all_scores": results.get("scores"),
            "best_score": results.get("best_score"),
            "num_poses": len(results.get("scores", [])),
            "message": f"Successfully docked molecule {smiles} against {target}. Best binding affinity: {results.get('best_score')} kcal/mol"
        }
            
    except Exception as e:
        logger.exception("Error occurred during docking simulation")
        return {"error": f"An error occurred: {str(e)}"}


@mcp.tool()
def list_available_targets() -> Dict[str, Any]:
    """
    Get list of all available target proteins for docking simulation
    
    Returns:
        Dictionary containing list of available targets with details
    """
    try:
        # Get list of available targets
        targets = get_available_targets()
        
        # Add target details if available
        target_details = {}
        for target in targets:
            target_details[target] = {
                "name": target,
                "description": f"Target protein {target} available for docking"
            }
        
        return {
            "targets": targets,
            "target_count": len(targets),
            "target_details": target_details,
            "message": f"Found {len(targets)} available target proteins for docking simulation"
        }
        
    except Exception as e:
        logger.exception("Error occurred while retrieving available targets")
        return {"error": f"An error occurred: {str(e)}"}


@mcp.tool()
def get_target_info(target: str = "akt1") -> Dict[str, Any]:
    """
    Get detailed information about a specific target protein
    
    Args:
        target: Name of the target protein
        
    Returns:
        Dictionary containing target configuration and details
    """
    try:
        from molDocking.docking import config, DOCKING_DIR
        
        if target not in config:
            available_targets = get_available_targets()
            return {
                "error": f"Target '{target}' not found. Available targets: {available_targets}"
            }
        
        target_config = config[target]
        target_dir = DOCKING_DIR / target
        
        return {
            "target": target,
            "receptor_file": target_config.get("receptor"),
            "center_coordinates": {
                "x": target_config.get("center_x"),
                "y": target_config.get("center_y"), 
                "z": target_config.get("center_z")
            },
            "search_space": {
                "size_x": target_config.get("size_x"),
                "size_y": target_config.get("size_y"),
                "size_z": target_config.get("size_z")
            },
            "computational_settings": {
                "cpu": target_config.get("cpu"),
                "exhaustiveness": target_config.get("exhaustiveness")
            },
            "target_directory": str(target_dir),
            "message": f"Configuration details for target {target}"
        }
        
    except Exception as e:
        logger.exception("Error occurred while retrieving target information")
        return {"error": f"An error occurred: {str(e)}"}


if __name__ == "__main__":
    try:
        print("Starting MCP server for molecular docking simulation...", file=sys.stderr)
        mcp.run()
    except Exception as e:
        print(f"Server startup error: {str(e)}", file=sys.stderr)
        sys.exit(1)
