import numpy as np
from Bio.PDB import PDBParser, Superimposer
import warnings

def GetCoords(pdb_structure, atom_type='P'):
    """
    Extract coordinates from PDB structure for specified atom type.
    For RNA, phosphorus (P) atoms are commonly used as they define the backbone.
    
    Args:
        pdb_structure: BioPython Structure object
        atom_type: Atom name to extract (default 'P' for phosphorus backbone)
    
    Returns:
        numpy array of coordinates (N x 3)
    """
    coords = []
    for model in pdb_structure:
        for chain in model:
            for residue in chain:
                if atom_type in residue:
                    atom = residue[atom_type]
                    coords.append(atom.get_coord())
    return np.array(coords)

def kabsch_alignment(coords1, coords2):
    """
    Perform Kabsch algorithm to find optimal rotation matrix.
    
    Args:
        coords1, coords2: numpy arrays of coordinates (N x 3)
    
    Returns:
        Rotation matrix and translation vectors
    """
    # Center the coordinates
    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)
    
    centered1 = coords1 - centroid1
    centered2 = coords2 - centroid2
    
    # Compute covariance matrix
    H = centered1.T @ centered2
    
    # SVD
    U, S, Vt = np.linalg.svd(H)
    
    # Compute rotation matrix
    R = Vt.T @ U.T
    
    # Ensure proper rotation (det(R) = 1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    return R, centroid1, centroid2

def CalculateRMSDScore(pdb1, pdb2, atom_type='P'):
    """
    Calculate RMSD (Root Mean Square Deviation) between two RNA structures.
    
    Args:
        pdb1, pdb2: BioPython Structure objects or file paths
        atom_type: Atom type to use for alignment (default 'P' for RNA backbone)
    
    Returns:
        RMSD value in Angstroms
    """
    # Parse PDB files if paths are provided
    parser = PDBParser(QUIET=True)
    
    if isinstance(pdb1, str):
        pdb1 = parser.get_structure('struct1', pdb1)
    if isinstance(pdb2, str):
        pdb2 = parser.get_structure('struct2', pdb2)
    
    # Extract coordinates
    coords1 = GetCoords(pdb1, atom_type)
    coords2 = GetCoords(pdb2, atom_type)
    
    # Check if structures have same number of atoms
    if len(coords1) != len(coords2):
        warnings.warn(f"Structures have different lengths: {len(coords1)} vs {len(coords2)}. Using minimum length.")
        min_len = min(len(coords1), len(coords2))
        coords1 = coords1[:min_len]
        coords2 = coords2[:min_len]
    
    if len(coords1) == 0:
        raise ValueError("No atoms found for alignment")
    
    # Perform Kabsch alignment
    R, centroid1, centroid2 = kabsch_alignment(coords1, coords2)
    
    # Apply transformation
    centered1 = coords1 - centroid1
    centered2 = coords2 - centroid2
    aligned2 = (R @ centered2.T).T
    
    # Calculate RMSD
    diff = centered1 - aligned2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd

def CalculateTMScore(pdb1, pdb2, atom_type='P'):
    """
    Calculate TM-score (Template Modeling score) between two RNA structures.
    TM-score ranges from 0 to 1, where 1 indicates identical structures.
    TM-score > 0.5 generally indicates similar fold.
    
    Args:
        pdb1, pdb2: BioPython Structure objects or file paths
        atom_type: Atom type to use for alignment (default 'P' for RNA backbone)
    
    Returns:
        TM-score value (0 to 1)
    """
    # Parse PDB files if paths are provided
    parser = PDBParser(QUIET=True)
    
    if isinstance(pdb1, str):
        pdb1 = parser.get_structure('struct1', pdb1)
    if isinstance(pdb2, str):
        pdb2 = parser.get_structure('struct2', pdb2)
    
    # Extract coordinates
    coords1 = GetCoords(pdb1, atom_type)
    coords2 = GetCoords(pdb2, atom_type)
    
    # Check if structures have same number of atoms
    if len(coords1) != len(coords2):
        warnings.warn(f"Structures have different lengths: {len(coords1)} vs {len(coords2)}. Using minimum length.")
        min_len = min(len(coords1), len(coords2))
        coords1 = coords1[:min_len]
        coords2 = coords2[:min_len]
    
    if len(coords1) == 0:
        raise ValueError("No atoms found for alignment")
    
    # Length for normalization (typically use target length)
    L_target = len(coords1)
    
    # Calculate d0 (distance scale)
    # d0 = 1.24 * (L_target - 15)^(1/3) - 1.8 for L > 21
    # This is the standard TM-score normalization
    if L_target > 21:
        d0 = 1.24 * ((L_target - 15) ** (1/3)) - 1.8
    else:
        d0 = 0.5  # For short sequences
    
    # Perform Kabsch alignment
    R, centroid1, centroid2 = kabsch_alignment(coords1, coords2)
    
    # Apply transformation
    centered1 = coords1 - centroid1
    centered2 = coords2 - centroid2
    aligned2 = (R @ centered2.T).T
    
    # Calculate distances
    diff = centered1 - aligned2
    distances = np.sqrt(np.sum(diff**2, axis=1))
    
    # Calculate TM-score
    # TM-score = (1/L_target) * sum(1 / (1 + (di/d0)^2))
    tm_score = np.mean(1.0 / (1.0 + (distances / d0) ** 2))
    
    return tm_score


# Example usage
if __name__ == "__main__":
    # Example with file paths
    # rmsd = CalculateRMSDScore('structure1.pdb', 'structure2.pdb')
    # tm_score = CalculateTMScore('structure1.pdb', 'structure2.pdb')
    
    # Example with BioPython structures
    parser = PDBParser(QUIET=True)
    # struct1 = parser.get_structure('s1', 'structure1.pdb')
    # struct2 = parser.get_structure('s2', 'structure2.pdb')
    # 
    # rmsd = CalculateRMSDScore(struct1, struct2)
    # tm_score = CalculateTMScore(struct1, struct2)
    # 
    # print(f"RMSD: {rmsd:.3f} Å")
    # print(f"TM-score: {tm_score:.4f}")
    
    # You can also use different atom types for alignment:
    # rmsd_c4 = CalculateRMSDScore(struct1, struct2, atom_type="C4'")  # Sugar carbon
    # rmsd_all = CalculateRMSDScore(struct1, struct2, atom_type='CA')  # All heavy atoms
    
    print("Functions ready to use!")
    print("Usage:")
    print("  rmsd = CalculateRMSDScore(pdb1, pdb2)")
    print("  tm_score = CalculateTMScore(pdb1, pdb2)")