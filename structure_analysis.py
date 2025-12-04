import numpy as np
from Bio.PDB import PDBParser, Superimposer, PDBIO
import warnings
from pathlib import Path

############################################################
# Main functions for superimposing and RMSD/TM calculations
# Peneeta Wojcik
# Gen AI Biomed - Fall 2025
############################################################

def GetAtoms(pdb_structure, atom_type='C3\''):
    """
    Extract specific atoms from PDB structure.
    
    Args:
        pdb_structure: BioPython Structure object
        atom_type: Atom name to extract (default C3' for RNA backbone)
    
    Returns:
        List of atom objects
    """
    atoms = []
    for model in pdb_structure:
        for chain in model:
            for residue in chain:
                atom_names = [atom_type, atom_type.replace("'", "*")]
                for name in atom_names:
                    if name in residue:
                        atoms.append(residue[name])
                        break
    return atoms

def CalculateRMSDScore(pdb1, pdb2, atom_type='C3\''):
    """
    Calculate RMSD between two RNA structures using BioPython's Superimposer.
    
    Args:
        pdb1, pdb2: BioPython Structure objects or file paths
        atom_type: Atom type to use (default C3' for RNA backbone)
    
    Returns:
        RMSD value in Angstroms
    """
    # Parse PDB files if paths are provided
    parser = PDBParser(QUIET=True)
    
    if isinstance(pdb1, str):
        struct1 = parser.get_structure('struct1', pdb1)
    else:
        struct1 = pdb1
        
    if isinstance(pdb2, str):
        struct2 = parser.get_structure('struct2', pdb2)
    else:
        struct2 = pdb2
    
    # Extract atoms
    atoms1 = GetAtoms(struct1, atom_type)
    atoms2 = GetAtoms(struct2, atom_type)
    
    if len(atoms1) == 0 or len(atoms2) == 0:
        raise ValueError(f"No {atom_type} atoms found for alignment")
    
    # Handle different lengths
    if len(atoms1) != len(atoms2):
        warnings.warn(f"Structures have different lengths: {len(atoms1)} vs {len(atoms2)}. Using minimum length.")
        min_len = min(len(atoms1), len(atoms2))
        atoms1 = atoms1[:min_len]
        atoms2 = atoms2[:min_len]
    
    # Use Superimposer
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    
    return super_imposer.rms

def CalculateTMScore(pdb1, pdb2, atom_type='C3\''):
    """
    Calculate TM-score between two RNA structures.
    Uses BioPython's Superimposer for alignment.
    
    Args:
        pdb1, pdb2: BioPython Structure objects or file paths
        atom_type: Atom type to use (default C3' for RNA backbone)
    
    Returns:
        TM-score value (0 to 1)
    """
    # Parse PDB files if paths are provided
    parser = PDBParser(QUIET=True)
    
    if isinstance(pdb1, str):
        struct1 = parser.get_structure('struct1', pdb1)
    else:
        struct1 = pdb1
        
    if isinstance(pdb2, str):
        struct2 = parser.get_structure('struct2', pdb2)
    else:
        struct2 = pdb2
    
    # Extract atoms
    atoms1 = GetAtoms(struct1, atom_type)
    atoms2 = GetAtoms(struct2, atom_type)
    
    if len(atoms1) == 0 or len(atoms2) == 0:
        raise ValueError(f"No {atom_type} atoms found for alignment")
    
    # Handle different lengths
    if len(atoms1) != len(atoms2):
        warnings.warn(f"Structures have different lengths: {len(atoms1)} vs {len(atoms2)}. Using minimum length.")
        min_len = min(len(atoms1), len(atoms2))
        atoms1 = atoms1[:min_len]
        atoms2 = atoms2[:min_len]
    
    # Use Superimposer
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    
    # Get coordinates after superposition
    coords1 = np.array([atom.get_coord() for atom in atoms1])
    coords2 = np.array([atom.get_coord() for atom in atoms2])
    
    # Apply the rotation and translation to coords2
    super_imposer.apply(atoms2)
    coords2_aligned = np.array([atom.get_coord() for atom in atoms2])
    
    # Calculate distances
    distances = np.sqrt(np.sum((coords1 - coords2_aligned)**2, axis=1))
    
    # Calculate TM-score
    L_target = len(atoms1)
    if L_target > 21:
        d0 = 1.24 * ((L_target - 15) ** (1/3)) - 1.8
    else:
        d0 = 0.5
    
    tm_score = np.mean(1.0 / (1.0 + (distances / d0) ** 2))
    
    return tm_score

def SaveSuperimposedStructures(pdb1, pdb2, output_prefix='aligned', atom_type='C3\''):
    """
    Superimpose two structures using BioPython's Superimposer and save as PDB files.
    
    Args:
        pdb1, pdb2: BioPython Structure objects or file paths
        output_prefix: Prefix for output files (default 'aligned')
        atom_type: Atom type to use for alignment (default C3')
    
    Returns:
        Dictionary with alignment info including filenames and scores
    """
    # Parse PDB files if paths are provided
    parser = PDBParser(QUIET=True)
    
    if isinstance(pdb1, str):
        struct1 = parser.get_structure('struct1', pdb1)
    else:
        struct1 = pdb1
        
    if isinstance(pdb2, str):
        struct2 = parser.get_structure('struct2', pdb2)
    else:
        struct2 = pdb2
    
    # Extract atoms for alignment
    atoms1 = GetAtoms(struct1, atom_type)
    atoms2 = GetAtoms(struct2, atom_type)
    
    # Handle different lengths
    if len(atoms1) != len(atoms2):
        warnings.warn(f"Structures have different lengths: {len(atoms1)} vs {len(atoms2)}. Using minimum length.")
        min_len = min(len(atoms1), len(atoms2))
        atoms1 = atoms1[:min_len]
        atoms2 = atoms2[:min_len]
    
    # Use Superimposer
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    rmsd = super_imposer.rms
    
    # Apply transformation to ALL atoms in structure 2
    # Get all atoms from structure 2
    all_atoms2 = []
    for model in struct2:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    all_atoms2.append(atom)
    
    # Apply the superposition transformation
    super_imposer.apply(all_atoms2)
    
    # Calculate TM-score
    coords1 = np.array([atom.get_coord() for atom in atoms1])
    coords2 = np.array([atom.get_coord() for atom in atoms2])
    distances = np.sqrt(np.sum((coords1 - coords2)**2, axis=1))
    
    L_target = len(atoms1)
    if L_target > 21:
        d0 = 1.24 * ((L_target - 15) ** (1/3)) - 1.8
    else:
        d0 = 0.5
    
    tm_score = np.mean(1.0 / (1.0 + (distances / d0) ** 2))
    
    # Save structures
    io = PDBIO()
    
    output1 = f"{output_prefix}_struct1.pdb"
    io.set_structure(struct1)
    io.save(output1)
    
    output2 = f"{output_prefix}_struct2_superimposed.pdb"
    io.set_structure(struct2)
    io.save(output2)
    
    print(f"Saved superimposed structures:")
    print(f"  Structure 1 (reference): {output1}")
    print(f"  Structure 2 (superimposed): {output2}")
    print(f"\nAlignment statistics:")
    print(f"Aligned atoms: {len(atoms1)}")
    print(f"RMSD: {rmsd:.3f} Å")
    print(f"TM-score: {tm_score:.4f}")
    print("\n")
    
    return {
        'file1': output1,
        'file2': output2,
        'rmsd': rmsd,
        'tm_score': tm_score,
        'aligned_atoms': len(atoms1),
        'total_atoms': (len(GetAtoms(struct1, atom_type)), 
                       len(GetAtoms(struct2, atom_type)))
    }

def CalculateRMSD_TM_Multiple(ref_aptamer, base_dir):
    base_dir = Path(base_dir)
    
    rmsd_arr = []
    tm_arr = []

    for subdir in base_dir.iterdir():
        if subdir.is_dir():
            # Construct path to the PDB file
            pdb_path = subdir / "relaxed_1000_model.pdb"
            
            # Calculate RMSD and TM scores
            rmsd = CalculateRMSDScore(ref_aptamer, str(pdb_path))
            tm = CalculateTMScore(ref_aptamer, str(pdb_path))
            
            # Append to arrays
            rmsd_arr.append(rmsd)
            tm_arr.append(tm)
            
            print(f"Processed {subdir.name}: RMSD={rmsd:.3f}, TM={tm:.3f}")
                    

    print(f"\nTotal processed: {len(rmsd_arr)} structures")
    print(f"Mean RMSD: {np.mean(rmsd_arr):.3f} ± {np.std(rmsd_arr):.3f}")
    print(f"Mean TM-score: {np.mean(tm_arr):.3f} ± {np.std(tm_arr):.3f}")
    
    return rmsd_arr, tm_arr