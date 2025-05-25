import logging
import os
import subprocess
import tempfile
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Import RDKit related modules
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    rdkit_available = True
except ImportError:
    logger.warning("RDKit is not installed. Some features may not be available.")
    rdkit_available = False

# Import Meeko for PDBQT conversion
try:
    from meeko import MoleculePreparation
    meeko_available = True
except ImportError:
    logger.warning("Meeko is not installed. Using fallback method for PDBQT conversion.")
    meeko_available = False

# Constants
PROJECT_ROOT = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DOCKING_DIR = PROJECT_ROOT.parent / "docking"

# Load configuration
config = {
    'akt1': {
        'receptor': f'{DOCKING_DIR}/akt1/receptor.pdbqt',
        'ligand': 'ligand.pdbqt',
        'out': 'output.pdbqt',
        'center_x': 5.9902,
        'center_y': 3.0140,
        'center_z': 17.3449,
        'size_x': 25,
        'size_y': 25,
        'size_z': 25,
        'cpu': 8,
        'exhaustiveness': 8
    },
}


def prepare_ligand_from_smiles(smiles: str) -> str:
    """
    Prepare ligand PDBQT file from SMILES string
    
    Args:
        smiles: SMILES string of the ligand
    
    Returns:
        Path to the prepared PDBQT file
        
    Raises:
        ValueError: If SMILES is invalid or conversion fails
    """
    if not smiles or not smiles.strip():
        raise ValueError("SMILES string cannot be empty")
        
    logger.info(f"Preparing ligand from SMILES: {smiles}")
    # 一時ディレクトリを使用
    pdbqt_path = tempfile.mkdtemp()
    output_pdbqt = os.path.join(pdbqt_path, f'{smiles}.pdbqt')
    
    # Meekeが利用可能ならそれを使用
    if meeko_available and rdkit_available:
        try:
            # RDKitで分子を生成
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
                
            # 3D構造を生成
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Meekeで分子を準備
            preparator = MoleculePreparation(hydrate=True)
            preparator.prepare(mol)
            pdbqt_string = preparator.write_pdbqt_string()
            
            # PDBQTファイルに書き込み
            with open(output_pdbqt, 'w') as f:
                f.write(pdbqt_string)
                
            return output_pdbqt
            
        except Exception as e:
            logger.warning(f"Meeko conversion failed: {str(e)}. Falling back to scrub.py")
    
    # フォールバック：scrub.pyを使用
    try:
        subprocess.run(
            ['scrub.py', f'{smiles}', '-o', f'{smiles}.sdf', '--skip_tautomers', '--ph_low', '5', '--ph_high', '9'],
            capture_output=True, text=True, cwd=pdbqt_path, check=True
        )
        
        # mk_prepare_ligand.pyで分子を準備
        subprocess.run(
            ['mk_prepare_ligand.py', '-i', f'{smiles}.sdf', '-o', output_pdbqt],
            capture_output=True, text=True, cwd=pdbqt_path, check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"Fallback ligand preparation failed: {e}")
        raise ValueError(f"Failed to prepare ligand from SMILES {smiles}: {e}")
    
    # SDFファイルを削除
    sdf_path = os.path.join(pdbqt_path, f'{smiles}.sdf')
    if os.path.exists(sdf_path):
        os.unlink(sdf_path)
    
    return output_pdbqt


def run_docking(target_name: str, ligand_file: str) -> Dict[str, Any]:
    """
    Run molecular docking using Vina as a command line tool
    
    Args:
        target_name: Name of the target protein (e.g., 'akt1')
        ligand_file: Path to the ligand PDBQT file
    
    Returns:
        Dictionary containing docking results
        
    Raises:
        ValueError: If target not found or docking fails
    """
    logger.info(f"Running docking for target: {target_name}")
    target_dir = DOCKING_DIR / target_name
    
    if not target_dir.exists() or not (target_dir / "receptor.pdbqt").exists():
        raise ValueError(f"Target '{target_name}' not found or receptor file missing")
    
    # 一時ディレクトリに出力ファイルを配置
    output_dir = tempfile.mkdtemp()
    output_file = os.path.join(output_dir, "output.pdbqt")
    
    # 設定ファイルを作成
    config_path = os.path.join(output_dir, "config.txt")
    config_target = config[target_name]
    
    with open(config_path, 'w') as f:
        f.write(f'receptor = {config_target["receptor"]}\n')
        f.write(f'ligand = {ligand_file}\n')
        f.write(f'out = {output_file}\n')
        f.write(f'center_x = {config_target["center_x"]}\n')
        f.write(f'center_y = {config_target["center_y"]}\n')
        f.write(f'center_z = {config_target["center_z"]}\n')
        f.write(f'size_x = {config_target["size_x"]}\n')
        f.write(f'size_y = {config_target["size_y"]}\n')
        f.write(f'size_z = {config_target["size_z"]}\n')
        f.write(f'cpu = {config_target["cpu"]}\n')
        f.write(f'exhaustiveness = {config_target["exhaustiveness"]}\n')
    
    # Run Vina as command line process with config file
    vina_cmd = ["vina", "--config", config_path]
    try:
        result = subprocess.run(
            vina_cmd, 
            capture_output=True, 
            text=True, 
            check=True
        )
        logger.info("Vina docking completed successfully")
    except subprocess.CalledProcessError as e:
        logger.error(f"Vina docking failed: {e}")
        raise ValueError(f"Docking failed for target {target_name}: {e}")
    except FileNotFoundError:
        logger.error("Vina executable not found")
        raise ValueError("AutoDock Vina is not installed or not in PATH")
    
    # Parse output to get scores
    scores = []
    energy = None
    
    # Look for the results table in stdout
    lines = result.stdout.split("\n")
    in_results_section = False
    
    for line in lines:
        # Check if we've reached the results section
        if "mode |   affinity | dist from best mode" in line:
            in_results_section = True
            continue
        
        # Skip header lines
        if in_results_section and ("-----" in line or "rmsd" in line or not line.strip()):
            continue
            
        # Parse score lines like "   1       -5.842          0          0"
        if in_results_section and line.strip():
            try:
                parts = line.split()
                if len(parts) >= 2:
                    # Second column contains the affinity score
                    score = float(parts[1])
                    scores.append(score)
                    # First score is the best score
                    if energy is None:
                        energy = score
            except (ValueError, IndexError):
                # If we can't parse a line, we might be past the results section
                break
    
    # Generate results
    results = {
        "target": target_name,
        "energy": energy,
        "scores": scores,
        "best_score": min(scores) if scores else None,
        "output_file": output_file  # Include path to output file with docked poses
    }
    
    return results


def get_available_targets() -> List[str]:
    """
    List all available target proteins
    
    Returns:
        List of target names
    """
    logger.info("Retrieving available targets")
    targets = []
    try:
        if not DOCKING_DIR.exists():
            logger.warning(f"Docking directory not found: {DOCKING_DIR}")
            return targets
            
        for d in DOCKING_DIR.iterdir():
            if d.is_dir() and (d / "receptor.pdbqt").exists():
                targets.append(d.name)
                
        logger.info(f"Found {len(targets)} available targets: {targets}")
        return targets
    except Exception as e:
        logger.error(f"Error retrieving targets: {e}")
        return []


def dock_molecule(smiles: str, target: str = "akt1") -> Dict[str, Any]:
    """
    Dock a molecule against a target protein
    
    Args:
        smiles: SMILES string of the molecule
        target: Name of the target protein (default: 'akt1')
        
    Returns:
        Dictionary containing docking results with error handling
    """
    logger.info(f"Starting docking simulation: SMILES={smiles}, target={target}")
    try:
        # Validate inputs
        if not smiles or not smiles.strip():
            raise ValueError("SMILES string cannot be empty")
            
        if target not in config:
            available_targets = get_available_targets()
            if not available_targets:
                raise ValueError("No targets available")
            if target not in available_targets:
                raise ValueError(f"Target '{target}' not available. Available targets: {available_targets}")
        
        # Prepare ligand
        logger.info("Preparing ligand...")
        ligand_file = prepare_ligand_from_smiles(smiles)
        
        # Run docking
        logger.info("Running docking simulation...")
        results = run_docking(target, ligand_file)
        
        # Clean up - remove ligand file and its parent directory
        try:
            ligand_dir = os.path.dirname(ligand_file)
            if os.path.exists(ligand_file):
                os.unlink(ligand_file)
            if os.path.exists(ligand_dir):
                os.rmdir(ligand_dir)
            
            # Clean up output file and its directory
            output_file = results.get("output_file")
            if output_file:
                output_dir = os.path.dirname(output_file)
                if os.path.exists(output_file):
                    os.unlink(output_file)
                    
                # Remove config.txt file in the output directory
                config_file = os.path.join(output_dir, "config.txt")
                if os.path.exists(config_file):
                    os.unlink(config_file)
                    
                # Remove the empty directory
                if os.path.exists(output_dir):
                    os.rmdir(output_dir)
        except Exception as cleanup_error:
            logger.warning(f"Cleanup failed: {cleanup_error}")
        
        # Remove output_file from results as it's no longer valid
        if "output_file" in results:
            del results["output_file"]
            
        logger.info(f"Docking completed successfully. Best score: {results.get('best_score')}")
        return results
    
    except Exception as e:
        logger.error(f"Docking failed: {e}")
        return {"error": str(e)}
