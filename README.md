# Molecular Docking MCP Server

A Model Context Protocol (MCP) server for molecular docking simulation using AutoDock Vina. This server enables docking small molecules against protein targets by providing SMILES strings and target names.

## Features

- **Molecular Docking**: Dock small molecules against protein targets using AutoDock Vina
- **SMILES Input**: Accept molecular structures in SMILES notation
- **Multiple Targets**: Support for multiple protein targets (currently includes AKT1)
- **MCP Protocol**: Standardized API through Model Context Protocol
- **Comprehensive Logging**: Detailed logging for debugging and monitoring
- **Error Handling**: Robust error handling with informative messages

## Installation

### Prerequisites

1. **AutoDock Vina**: Install AutoDock Vina for molecular docking from https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.2/vina_1.2.2_macos_arm64
   Please make sure that the vina command is in your PATH.
   ```bash
    which vina
   ```

2. **Python Dependencies**: Install using uv or pip
   ```bash
   # Using uv (recommended)
   uv pip install -e .
   ```

### Required Dependencies

- `rdkit>=2024.9.6` - Chemical informatics toolkit
- `mcp>=1.2.0` - Model Context Protocol
- `meeko>=0.5.0` - PDBQT file preparation
- `pathlib` - Path handling utilities

## Usage

### Starting the MCP Server

```bash
# From the project directory
python server.py
```

### Available Tools

The MCP server provides three main tools:

#### 1. `run_molecular_docking`
Dock a molecule against a target protein.

**Parameters:**
- `smiles` (str): SMILES string of the molecule to dock
- `target` (str, optional): Target protein name (default: "akt1")

**Returns:**
- `smiles`: Input SMILES string
- `target`: Target protein name
- `binding_affinity`: Best binding affinity (kcal/mol)
- `all_scores`: List of all binding scores
- `best_score`: Best (lowest) binding score
- `num_poses`: Number of generated poses

#### 2. `list_available_targets`
Get list of all available target proteins.

**Returns:**
- `targets`: List of target names
- `target_count`: Number of available targets
- `target_details`: Detailed information for each target

#### 3. `get_target_info`
Get detailed configuration for a specific target.

**Parameters:**
- `target` (str, optional): Target protein name (default: "akt1")

**Returns:**
- `target`: Target name
- `receptor_file`: Path to receptor PDBQT file
- `center_coordinates`: Binding site center coordinates
- `search_space`: Search space dimensions
- `computational_settings`: CPU and exhaustiveness settings

### Example Usage

```python
# Example docking simulation
from molDocking import dock_molecule, get_available_targets

# List available targets
targets = get_available_targets()
print(f"Available targets: {targets}")

# Dock aspirin against AKT1
result = dock_molecule("CC(=O)OC1=CC=CC=C1C(=O)O", "akt1")
print(f"Best binding affinity: {result['best_score']} kcal/mol")
```

## Target Configuration

### Current Targets

- **AKT1**: Protein kinase involved in cell survival and metabolism

### Adding New Targets

To add a new target:

1. Create a directory in `../docking/` with the target name
2. Add the receptor PDBQT file as `receptor.pdbqt`
3. Update the configuration in `molDocking/docking.py`:

```python
config = {
    'your_target': {
        'receptor': f'{DOCKING_DIR}/your_target/receptor.pdbqt',
        'center_x': 0.0,      # Binding site center X
        'center_y': 0.0,      # Binding site center Y  
        'center_z': 0.0,      # Binding site center Z
        'size_x': 25,         # Search space size X
        'size_y': 25,         # Search space size Y
        'size_z': 25,         # Search space size Z
        'cpu': 8,             # Number of CPU cores
        'exhaustiveness': 8   # Search exhaustiveness
    }
}
```

## MCP Configuration

To use this server with an MCP client, use the configuration in [config file](./calude-desktop-config.json)


## File Structure

```
molDocking/
├── README.md
├── pyproject.toml
├── server.py              # MCP server entry point
├── docking/               # Target protein files
│   └── akt1/
│       ├── config.txt
│       └── receptor.pdbqt
└── molDocking/
    ├── __init__.py        # Package initialization
    ├── __main__.py        # Module entry point
    └── docking.py         # Core docking functionality

```

## Troubleshooting

### Common Issues

1. **Vina not found**: Ensure AutoDock Vina is installed and in PATH
2. **Missing dependencies**: Install all required Python packages
3. **Invalid SMILES**: Check SMILES string format
4. **Target not found**: Verify target exists in configuration

```

## Development

### Running Tests

```bash
# Install development dependencies
uv pip install -e .[dev]

# Run tests
pytest
```

### Code Formatting

```bash
# Format code
black molDocking/
isort molDocking/
```

## License

MIT License

## Dependencies Credits

- **RDKit**: Chemical informatics and machine learning toolkit
- **AutoDock Vina**: Molecular docking program
- **Meeko**: PDBQT file preparation toolkit
- **MCP**: Model Context Protocol for AI integration

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Support

For issues and questions:
1. Check the troubleshooting section
2. Review the logs for error details
3. Create an issue on the project repository