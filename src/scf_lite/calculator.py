"""
Módulo principal para ejecutar cálculos SCF usando PySCF
"""

from typing import List, Dict, Any
from pyscf import gto, scf


def calculate_scf(
    symbols: List[str],
    coordinates: List[List[float]],
    charge: int = 0,
    spin: int = 0,
    basis: str = "sto-3g"
) -> Dict[str, Any]:
    """
    Ejecuta un cálculo SCF (RHF o UHF según corresponda).
    
    Args:
        symbols: Lista de símbolos químicos
        coordinates: Lista de coordenadas 3D en Angstrom
        charge: Carga total
        spin: Spin total (0 = singlete/RHF, != 0 = UHF)
        basis: Base a usar
    
    Returns:
        Dict con los resultados del cálculo
    """
    # Construir la cadena de átomos para PySCF
    atom_string = []
    for symbol, coord in zip(symbols, coordinates):
        atom_string.append(f"{symbol} {coord[0]} {coord[1]} {coord[2]}")
    
    # Crear objeto molecular
    mol = gto.M(
        atom=atom_string,
        basis=basis,
        charge=charge,
        spin=spin
    )
    
    # Elegir método SCF según el spin
    if spin == 0:
        mf = scf.RHF(mol)
    else:
        mf = scf.UHF(mol)
    
    # Ejecutar cálculo
    energia = mf.kernel()
    
    # Recopilar resultados
    resultados = {
        "energia": float(energia),
        "convergio": mf.converged,
        "iteraciones": mf.scf_summary.get("niter", 0),
        "metodo": "RHF" if spin == 0 else "UHF",
        "spin": spin,
        "charge": charge,
        "basis": basis,
        "natom": mol.natm,
        "nelec": mol.nelec
    }
    
    return resultados
