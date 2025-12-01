"""
Módulo para validar y cargar inputs del usuario
"""

import json
from typing import List, Dict, Any, Tuple


def validate_input(
    symbols: List[str],
    coordinates: List[List[float]],
    charge: int = None,
    spin: int = None,
    basis: str = "sto-3g"
) -> Tuple[bool, str]:
    """
    Valida que el input sea correcto para PySCF.
    
    Args:
        symbols: Lista de símbolos químicos
        coordinates: Lista de coordenadas 3D en Angstrom
        charge: Carga total (opcional)
        spin: Spin total (opcional)
        basis: Base a usar (default: sto-3g)
    
    Returns:
        Tuple[bool, str]: (es_válido, mensaje_error)
    """
    # Validar que haya el mismo número de símbolos y coordenadas
    if len(symbols) != len(coordinates):
        return False, f"Número de símbolos ({len(symbols)}) no coincide con número de coordenadas ({len(coordinates)})"
    
    # Validar símbolos químicos básicos
    simbolos_validos = {
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca"
    }
    
    for symbol in symbols:
        if symbol not in simbolos_validos:
            return False, f"Símbolo químico no soportado: {symbol}"
    
    # Validar coordenadas
    for i, coord in enumerate(coordinates):
        if len(coord) != 3:
            return False, f"Coordenada {i} debe tener 3 componentes, tiene {len(coord)}"
        
        for j, val in enumerate(coord):
            if not isinstance(val, (int, float)):
                return False, f"Coordenada {i}[{j}] debe ser numérica, es {type(val)}"
    
    # Validar base
    bases_validas = ["sto-3g", "6-31g", "cc-pvdz", "def2-svp"]
    if basis not in bases_validas:
        return False, f"Base no soportada: {basis}. Bases válidas: {', '.join(bases_validas)}"
    
    return True, ""


def load_input_file(filepath: str) -> Dict[str, Any]:
    """
    Carga un archivo JSON con la configuración de la molécula.
    
    Formato esperado:
    {
        "symbols": ["O", "H", "H"],
        "coordinates": [[0.0, 0.0, 0.0], [0.0, -0.757, 0.587], [0.0, 0.757, 0.587]],
        "charge": 0,
        "spin": 0,
        "basis": "sto-3g"
    }
    
    Args:
        filepath: Ruta al archivo JSON
    
    Returns:
        Dict con la configuración
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    return data
