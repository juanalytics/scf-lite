"""
Módulo para validar y cargar inputs del usuario
"""

import json
from pathlib import Path
from typing import List, Dict, Any, Tuple


def validate_input(
    symbols: List[str],
    coordinates: List[List[float]],
    charge: int | None = None,
    spin: int | None = None,
    basis: str = "sto-3g",
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
        return (
            False,
            f"Número de símbolos ({len(symbols)}) no coincide con número de coordenadas ({len(coordinates)})",
        )

    # Validar símbolos químicos básicos (restringido para mantenerlo simple)
    simbolos_validos = {
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
    }

    for symbol in symbols:
        if symbol not in simbolos_validos:
            return False, f"Símbolo químico no soportado: {symbol}"

    # Validar coordenadas
    for i, coord in enumerate(coordinates):
        if len(coord) != 3:
            return (
                False,
                f"Coordenada {i} debe tener 3 componentes, tiene {len(coord)}",
            )

        for j, val in enumerate(coord):
            if not isinstance(val, (int, float)):
                return (
                    False,
                    f"Coordenada {i}[{j}] debe ser numérica, es {type(val)}",
                )

    # Validar base
    bases_validas = ["sto-3g", "6-31g", "cc-pvdz", "def2-svp"]
    if basis not in bases_validas:
        return (
            False,
            f"Base no soportada: {basis}. Bases válidas: {', '.join(bases_validas)}",
        )

    return True, ""


def _load_json_input(filepath: str) -> Dict[str, Any]:
    """
    Carga un archivo JSON con la configuración de la molécula.

    Formato esperado (estándar):
    {
        "name": "agua",        # opcional
        "symbols": ["O", "H", "H"],
        "coordinates": [[0.0, 0.0, 0.0], [0.0, -0.757, 0.587], [0.0, 0.757, 0.587]],
        "charge": 0,
        "spin": 0,
        "basis": "sto-3g"
    }
    """
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)

    return data


def _load_xyz_input(filepath: str) -> Dict[str, Any]:
    """
    Carga un archivo XYZ simple y lo convierte al formato interno.

    Formato esperado:
        N
        comentario opcional
        Sym x y z
        Sym x y z
        ...

    La carga, el spin y la base se toman por defecto como:
        charge = 0, spin = 0, basis = "sto-3g"
    y pueden sobreescribirse luego desde la CLI con --charge/--spin/--basis.
    """
    symbols: List[str] = []
    coordinates: List[List[float]] = []

    with open(filepath, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    if len(lines) < 3:
        raise ValueError("Archivo XYZ demasiado corto.")

    try:
        n_atoms = int(lines[0])
    except ValueError as exc:
        raise ValueError("Primera línea del XYZ debe ser el número de átomos.") from exc

    atom_lines = lines[2 : 2 + n_atoms]
    if len(atom_lines) != n_atoms:
        raise ValueError(
            f"El archivo XYZ indica {n_atoms} átomos pero se encontraron {len(atom_lines)} líneas de coordenadas."
        )

    for line in atom_lines:
        parts = line.split()
        if len(parts) != 4:
            raise ValueError(
                f"Línea XYZ inválida: '{line}'. Se espera: 'Simbolo x y z'."
            )
        sym = parts[0]
        x, y, z = map(float, parts[1:4])
        symbols.append(sym)
        coordinates.append([x, y, z])

    return {
        "name": Path(filepath).stem,
        "symbols": symbols,
        "coordinates": coordinates,
        "charge": 0,
        "spin": 0,
        "basis": "sto-3g",
    }


def load_input_file(filepath: str) -> Dict[str, Any]:
    """
    Carga un archivo de entrada de usuario y lo normaliza a un diccionario estándar.

    Soporta:
        - JSON con campos: symbols, coordinates, charge, spin, basis (name opcional)
        - XYZ: se convierte a la misma estructura, con charge=0, spin=0, basis="sto-3g"
    """
    suffix = Path(filepath).suffix.lower()

    if suffix == ".json":
        return _load_json_input(filepath)
    if suffix == ".xyz":
        return _load_xyz_input(filepath)

    raise ValueError(
        f"Formato de archivo no soportado: {suffix}. Usa .json o .xyz para moléculas."
    )




