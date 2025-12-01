"""
Interfaz de línea de comandos para SCF Lite
"""

import argparse
import json
import sys
from pathlib import Path

from .input_validator import validate_input, load_input_file
from .calculator import calculate_scf
from .output_formatter import format_results


def run_h2_scan() -> None:
    """
    Escaneo rápido de la molécula H2 (energía vs distancia)
    y gráfico interactivo usando matplotlib.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Rango de distancias H-H en Angstrom
    distances = np.linspace(0.4, 2.5, 30)
    energies = []

    symbols = ["H", "H"]

    for r in distances:
        coordinates = [[0.0, 0.0, 0.0], [float(r), 0.0, 0.0]]
        resultados = calculate_scf(
            symbols,
            coordinates,
            charge=0,
            spin=0,
            basis="sto-3g",
        )
        energies.append(resultados["energia"])

    energies_arr = np.array(energies)
    distances_arr = np.array(distances)

    # Encontrar mínimo
    idx_min = int(np.argmin(energies_arr))
    r_min = float(distances_arr[idx_min])
    e_min = float(energies_arr[idx_min])

    # Conversión Hartree -> kcal/mol (para referencia)
    HARTREE_TO_KCAL = 627.509

    # Imprimir tabla en terminal
    print("Escaneo H2 (RHF / STO-3G)")
    print("R_HH (Å)   E (Hartree)     ΔE (kcal/mol)   nota")
    for r, e in zip(distances_arr, energies_arr):
        delta_e = (e - e_min) * HARTREE_TO_KCAL
        mark = "<-- mínimo" if abs(r - r_min) < 1e-8 else ""
        print(f"{r:7.3f}   {e: .8f}   {delta_e:10.3f}   {mark}")

    # Intentar ajuste tipo Morse: E(R) = E0 + De (1 - exp(-a (R-Re)))^2
    morse_params = None

    try:
        from scipy.optimize import curve_fit  # type: ignore[import]

        def morse(R, De, a, Re, E0):
            return E0 + De * (1.0 - np.exp(-a * (R - Re))) ** 2

        # Estimaciones iniciales razonables
        De0 = e_min - float(energies_arr[-1])
        a0 = 1.5
        Re0 = r_min
        E00 = float(energies_arr[-1])

        popt, _ = curve_fit(
            morse,
            distances_arr,
            energies_arr,
            p0=[De0, a0, Re0, E00],
            maxfev=10000,
        )
        morse_params = popt  # type: ignore[assignment]
    except Exception:
        morse_params = None

    # Gráfico
    plt.figure()
    plt.plot(distances_arr, energies_arr, marker="o", label="E(R) puntos SCF")
    plt.scatter([r_min], [e_min], color="red", zorder=5, label="Mínimo")
    plt.axvline(r_min, color="red", linestyle="--", alpha=0.5)

    if morse_params is not None:
        De_fit, a_fit, Re_fit, E0_fit = morse_params

        R_fine = np.linspace(distances_arr.min(), distances_arr.max(), 400)
        E_fine = E0_fit + De_fit * (1.0 - np.exp(-a_fit * (R_fine - Re_fit))) ** 2
        plt.plot(R_fine, E_fine, color="orange", label="Ajuste tipo Morse")

        # Texto formateado con sintaxis tipo LaTeX para que se vea limpio
        eq_text = (
            r"$E(R) \approx E_0 + D_e \left(1 - e^{-a (R-R_e)}\right)^2$" "\n"
            rf"$D_e = {De_fit:.3f}\ \mathrm{{Ha}}$" "\n"
            rf"$R_e = {Re_fit:.3f}\ \mathrm{{\AA}},\ a = {a_fit:.3f}\ \mathrm{{\AA}}^{{-1}}$" "\n"
            rf"$E_0 = {E0_fit:.3f}\ \mathrm{{Ha}}$"
        )
        plt.annotate(
            eq_text,
            xy=(0.05, 0.95),
            xycoords="axes fraction",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round", fc="white", alpha=0.85),
        )

        print("\nAjuste tipo Morse (aprox.):")
        print("E(R) = E0 + De (1 - exp(-a (R-Re)))^2")
        print(f"De  = {De_fit:.6f} Ha")
        print(f"Re  = {Re_fit:.6f} Å")
        print(f"a   = {a_fit:.6f} Å^-1")
        print(f"E0  = {E0_fit:.6f} Ha")

    plt.xlabel("Distancia H-H (Å)")
    plt.ylabel("Energía (Hartree)")
    plt.title("Curva de energía para H₂ (RHF / STO-3G)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def run_bond_scan(
    symbols: list[str],
    coordinates: list[list[float]],
    i: int,
    j: int,
    charge: int,
    spin: int,
    basis: str,
) -> None:
    """
    Escaneo genérico de un enlace entre los átomos i y j
    usando la geometría proporcionada por el usuario.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    natom = len(symbols)
    if not (0 <= i < natom and 0 <= j < natom):
        raise ValueError(f"Índices de enlace fuera de rango: i={i}, j={j}, natom={natom}")

    # Clonar coordenadas y pasar a numpy
    coords = np.array(coordinates, dtype=float)
    ri = coords[i]
    rj = coords[j]

    r0_vec = rj - ri
    r0 = float(np.linalg.norm(r0_vec))
    if r0 == 0.0:
        raise ValueError("Los átomos i y j tienen la misma posición; no se puede escanear el enlace.")

    direction = r0_vec / r0

    # Escanear de 70% a 130% de la distancia original
    distances = np.linspace(0.7 * r0, 1.3 * r0, 30)
    energies: list[float] = []

    sym_i = symbols[i]
    sym_j = symbols[j]

    for r in distances:
        coords_scan = coords.copy()
        # Fijamos átomo i y movemos j a lo largo de la dirección original
        coords_scan[j] = ri + direction * float(r)
        resultados = calculate_scf(
            symbols,
            coords_scan.tolist(),
            charge=charge,
            spin=spin,
            basis=basis,
        )
        energies.append(resultados["energia"])

    energies_arr = np.array(energies)
    distances_arr = np.array(distances)

    # Encontrar mínimo
    idx_min = int(np.argmin(energies_arr))
    r_min = float(distances_arr[idx_min])
    e_min = float(energies_arr[idx_min])

    HARTREE_TO_KCAL = 627.509

    print(f"Escaneo enlace {sym_i}{i+1}–{sym_j}{j+1} (RHF/UHF, base {basis})")
    print("R (Å)      E (Hartree)     ΔE (kcal/mol)   nota")
    for r, e in zip(distances_arr, energies_arr):
        delta_e = (e - e_min) * HARTREE_TO_KCAL
        mark = "<-- mínimo" if abs(r - r_min) < 1e-8 else ""
        print(f"{r:7.3f}   {e: .8f}   {delta_e:10.3f}   {mark}")

    # Ajuste tipo Morse
    morse_params = None
    try:
        from scipy.optimize import curve_fit  # type: ignore[import]

        def morse(R, De, a, Re, E0):
            return E0 + De * (1.0 - np.exp(-a * (R - Re))) ** 2

        De0 = e_min - float(energies_arr[-1])
        a0 = 1.5
        Re0 = r_min
        E00 = float(energies_arr[-1])

        popt, _ = curve_fit(
            morse,
            distances_arr,
            energies_arr,
            p0=[De0, a0, Re0, E00],
            maxfev=10000,
        )
        morse_params = popt  # type: ignore[assignment]
    except Exception:
        morse_params = None

    plt.figure()
    plt.plot(distances_arr, energies_arr, marker="o", label="E(R) puntos SCF")
    plt.scatter([r_min], [e_min], color="red", zorder=5, label="Mínimo")
    plt.axvline(r_min, color="red", linestyle="--", alpha=0.5)

    if morse_params is not None:
        De_fit, a_fit, Re_fit, E0_fit = morse_params

        R_fine = np.linspace(distances_arr.min(), distances_arr.max(), 400)
        E_fine = E0_fit + De_fit * (1.0 - np.exp(-a_fit * (R_fine - Re_fit))) ** 2
        plt.plot(R_fine, E_fine, color="orange", label="Ajuste tipo Morse")

        eq_text = (
            r"$E(R) \approx E_0 + D_e \left(1 - e^{-a (R-R_e)}\right)^2$" "\n"
            rf"$D_e = {De_fit:.3f}\ \mathrm{{Ha}}$" "\n"
            rf"$R_e = {Re_fit:.3f}\ \mathrm{{\AA}},\ a = {a_fit:.3f}\ \mathrm{{\AA}}^{{-1}}$" "\n"
            rf"$E_0 = {E0_fit:.3f}\ \mathrm{{Ha}}$"
        )
        plt.annotate(
            eq_text,
            xy=(0.05, 0.95),
            xycoords="axes fraction",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round", fc="white", alpha=0.85),
        )

        print("\nAjuste tipo Morse (enlace, aprox.):")
        print("E(R) = E0 + De (1 - exp(-a (R-Re)))^2")
        print(f"De  = {De_fit:.6f} Ha")
        print(f"Re  = {Re_fit:.6f} Å")
        print(f"a   = {a_fit:.6f} Å^-1")
        print(f"E0  = {E0_fit:.6f} Ha")

    plt.xlabel(f"Distancia {sym_i}{i+1}–{sym_j}{j+1} (Å)")
    plt.ylabel("Energía (Hartree)")
    plt.title(f"Curva de energía para enlace {sym_i}{i+1}–{sym_j}{j+1}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def run_oh_scan() -> None:
    """
    Escaneo del enlace O-H en la molécula de agua (H2O)
    manteniendo fija la geometría del resto.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Geometría base de agua (la misma que en examples/water.json)
    O = np.array([0.0, 0.0, 0.0])
    H1 = np.array([0.0, -0.757, 0.587])
    H2 = np.array([0.0, 0.757, 0.587])

    r0 = np.linalg.norm(H1 - O)
    direction = (H1 - O) / r0

    # Escaneamos alrededor de la distancia de equilibrio
    distances = np.linspace(0.7 * r0, 1.3 * r0, 30)
    energies = []

    symbols = ["O", "H", "H"]

    for r in distances:
        H1_new = O + direction * float(r)
        coordinates = [
            O.tolist(),
            H1_new.tolist(),
            H2.tolist(),
        ]
        resultados = calculate_scf(
            symbols,
            coordinates,
            charge=0,
            spin=0,
            basis="sto-3g",
        )
        energies.append(resultados["energia"])

    energies_arr = np.array(energies)
    distances_arr = np.array(distances)

    # Encontrar mínimo
    idx_min = int(np.argmin(energies_arr))
    r_min = float(distances_arr[idx_min])
    e_min = float(energies_arr[idx_min])

    HARTREE_TO_KCAL = 627.509

    print("Escaneo O-H en H2O (RHF / STO-3G)")
    print("R_OH (Å)   E (Hartree)     ΔE (kcal/mol)   nota")
    for r, e in zip(distances_arr, energies_arr):
        delta_e = (e - e_min) * HARTREE_TO_KCAL
        mark = "<-- mínimo" if abs(r - r_min) < 1e-8 else ""
        print(f"{r:7.3f}   {e: .8f}   {delta_e:10.3f}   {mark}")

    # Ajuste tipo Morse para el enlace O-H
    morse_params = None

    try:
        from scipy.optimize import curve_fit  # type: ignore[import]

        def morse(R, De, a, Re, E0):
            return E0 + De * (1.0 - np.exp(-a * (R - Re))) ** 2

        De0 = e_min - float(energies_arr[-1])
        a0 = 2.0
        Re0 = r_min
        E00 = float(energies_arr[-1])

        popt, _ = curve_fit(
            morse,
            distances_arr,
            energies_arr,
            p0=[De0, a0, Re0, E00],
            maxfev=10000,
        )
        morse_params = popt  # type: ignore[assignment]
    except Exception:
        morse_params = None

    plt.figure()
    plt.plot(distances_arr, energies_arr, marker="o", label="E(R) puntos SCF")
    plt.scatter([r_min], [e_min], color="red", zorder=5, label="Mínimo")
    plt.axvline(r_min, color="red", linestyle="--", alpha=0.5)

    if morse_params is not None:
        De_fit, a_fit, Re_fit, E0_fit = morse_params

        R_fine = np.linspace(distances_arr.min(), distances_arr.max(), 400)
        E_fine = E0_fit + De_fit * (1.0 - np.exp(-a_fit * (R_fine - Re_fit))) ** 2
        plt.plot(R_fine, E_fine, color="orange", label="Ajuste tipo Morse")

        eq_text = (
            r"$E(R) \approx E_0 + D_e \left(1 - e^{-a (R-R_e)}\right)^2$" "\n"
            rf"$D_e = {De_fit:.3f}\ \mathrm{{Ha}}$" "\n"
            rf"$R_e = {Re_fit:.3f}\ \mathrm{{\AA}},\ a = {a_fit:.3f}\ \mathrm{{\AA}}^{{-1}}$" "\n"
            rf"$E_0 = {E0_fit:.3f}\ \mathrm{{Ha}}$"
        )
        plt.annotate(
            eq_text,
            xy=(0.05, 0.95),
            xycoords="axes fraction",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round", fc="white", alpha=0.85),
        )

        print("\nAjuste tipo Morse (O-H, aprox.):")
        print("E(R) = E0 + De (1 - exp(-a (R-Re)))^2")
        print(f"De  = {De_fit:.6f} Ha")
        print(f"Re  = {Re_fit:.6f} Å")
        print(f"a   = {a_fit:.6f} Å^-1")
        print(f"E0  = {E0_fit:.6f} Ha")

    plt.xlabel("Distancia O-H (Å)")
    plt.ylabel("Energía (Hartree)")
    plt.title("Curva de energía para enlace O–H en H₂O (RHF / STO-3G)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    """Función principal de la CLI"""
    parser = argparse.ArgumentParser(
        description="SCF Lite - Sistema minimalista para cálculos SCF usando PySCF"
    )

    # Modos especiales de escaneo con gráfico
    parser.add_argument(
        "--scan-h2",
        action="store_true",
        help="Escaneo energía vs distancia H-H (H2) con gráfico",
    )
    parser.add_argument(
        "--scan-oh",
        action="store_true",
        help="Escaneo energía vs distancia O-H en H2O con gráfico",
    )
    parser.add_argument(
        "--scan-bond",
        nargs=2,
        type=int,
        metavar=("I", "J"),
        help="Escaneo genérico energía vs distancia entre átomos I y J (índices 1-based)",
    )

    # Opción para usar archivo de input
    parser.add_argument(
        "-f", "--file",
        type=str,
        help="Archivo JSON con la configuración de la molécula"
    )
    
    # Opciones para entrada directa (alternativa al archivo)
    parser.add_argument(
        "-s", "--symbols",
        nargs="+",
        help="Símbolos químicos (ej: O H H)"
    )
    
    parser.add_argument(
        "-c", "--coordinates",
        nargs="+",
        help="Coordenadas como lista plana (ej: 0.0 0.0 0.0 0.0 -0.757 0.587 0.0 0.757 0.587)"
    )
    
    parser.add_argument(
        "--charge",
        type=int,
        default=0,
        help="Carga total (default: 0)"
    )
    
    parser.add_argument(
        "--spin",
        type=int,
        default=0,
        help="Spin total (default: 0, singlete/RHF)"
    )
    
    parser.add_argument(
        "--basis",
        type=str,
        default="sto-3g",
        help="Base a usar (default: sto-3g)"
    )
    
    # Opciones de salida
    parser.add_argument(
        "-o", "--output",
        type=str,
        help="Archivo de salida (JSON)"
    )
    
    parser.add_argument(
        "--format",
        choices=["dict", "json"],
        default="json",
        help="Formato de salida (default: json)"
    )
    
    args = parser.parse_args()

    # Cargar input
    if args.file:
        # Cargar desde archivo
        try:
            data = load_input_file(args.file)
            symbols = data["symbols"]
            coordinates = data["coordinates"]
            charge = data.get("charge", 0)
            spin = data.get("spin", 0)
            basis = data.get("basis", "sto-3g")
        except FileNotFoundError:
            print(f"Error: No se encontró el archivo {args.file}", file=sys.stderr)
            sys.exit(1)
        except KeyError as e:
            print(f"Error: Falta el campo requerido en el archivo: {e}", file=sys.stderr)
            sys.exit(1)
        except json.JSONDecodeError as e:
            print(f"Error: El archivo JSON no es válido: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.symbols and args.coordinates:
        # Cargar desde argumentos de línea de comandos
        symbols = args.symbols
        coords_flat = [float(x) for x in args.coordinates]
        
        if len(coords_flat) % 3 != 0:
            print("Error: El número de coordenadas debe ser múltiplo de 3", file=sys.stderr)
            sys.exit(1)
        
        coordinates = [
            coords_flat[i:i+3] 
            for i in range(0, len(coords_flat), 3)
        ]
        charge = args.charge
        spin = args.spin
        basis = args.basis
    else:
        print("Error: Debes proporcionar un archivo (-f) o símbolos y coordenadas (-s, -c)", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    # Modos gráficos de escaneo dedicados (no requieren input previo)
    if args.scan_h2:
        run_h2_scan()
        return
    if args.scan_oh:
        run_oh_scan()
        return

    # Validar input
    es_valido, mensaje = validate_input(symbols, coordinates, charge, spin, basis)
    if not es_valido:
        print(f"Error de validación: {mensaje}", file=sys.stderr)
        sys.exit(1)

    # Si se solicita escanear un enlace específico, usar geometría del usuario
    if args.scan_bond:
        i_cli, j_cli = args.scan_bond
        # Convertimos de 1-based (CLI) a 0-based (interno)
        i_idx = i_cli - 1
        j_idx = j_cli - 1
        try:
            run_bond_scan(symbols, coordinates, i_idx, j_idx, charge, spin, basis)
        except Exception as e:
            print(f"Error durante el escaneo de enlace: {e}", file=sys.stderr)
            sys.exit(1)
        return

    # Ejecutar cálculo simple
    try:
        resultados = calculate_scf(symbols, coordinates, charge, spin, basis)
    except Exception as e:
        print(f"Error durante el cálculo: {e}", file=sys.stderr)
        sys.exit(1)

    # Formatear y mostrar resultados
    output = format_results(resultados, formato=args.format, output_file=args.output)

    if args.format == "json":
        print(output)
    else:
        print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()




