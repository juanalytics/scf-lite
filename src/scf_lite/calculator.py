"""
Main SCF calculation logic.
"""

from pyscf import gto, scf


def calculate_h2o_energy():
    """
    Minimal working example: Calculate energy for H₂O molecule.
    Hardcoded geometry and parameters.
    
    Returns:
        float: Electronic energy in Hartree
    """
    # Build H₂O molecule (hardcoded)
    # Geometry: O at origin, two H atoms
    # O-H distance ~0.957 Å, H-O-H angle ~104.5°
    mol = gto.Mole()
    mol.atom = '''
    O  0.000000  0.000000  0.000000
    H  0.957000  0.000000  0.000000
    H -0.240000  0.927000  0.000000
    '''
    mol.basis = '6-31g'
    mol.charge = 0
    mol.spin = 0  # singlet (multiplicity = 1)
    mol.build()
    
    # Run Hartree-Fock calculation
    mf = scf.RHF(mol)
    energy = mf.kernel()
    
    return energy


if __name__ == "__main__":
    # Test the minimal example
    print("Running minimal H₂O calculation...")
    energy = calculate_h2o_energy()
    print(f"Energy: {energy:.6f} Hartree")
    print(f"Energy: {energy * 27.2114:.4f} eV")
