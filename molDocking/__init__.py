"""
molDocking - Library for molecular docking simulation using AutoDock Vina
"""

from .docking import dock_molecule, get_available_targets, run_docking, prepare_ligand_from_smiles

__version__ = "0.1.0"