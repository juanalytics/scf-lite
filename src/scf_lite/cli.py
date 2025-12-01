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

        eq_text = (
            "E(R) ≈ E₀ + Dₑ (1 - e^{-a (R-Rₑ)})²\n"
            f"Dₑ = {De_fit:.3f} Ha, Rₑ = {Re_fit:.3f} Å\n"
            f"a = {a_fit:.3f} Å⁻¹, E₀ = {E0_fit:.3f} Ha"
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


def main():
    """Función principal de la CLI"""
    parser = argparse.ArgumentParser(
        description="SCF Lite - Sistema minimalista para cálculos SCF usando PySCF"
    )

    # Modo especial: escaneo H2 con gráfico
    parser.add_argument(
        "--scan-h2",
        action="store_true",
        help="Realiza un escaneo de energía vs distancia para H2 y muestra un gráfico",
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

    # Modo gráfico: escaneo de H2
    if args.scan_h2:
        run_h2_scan()
        return

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
    
    # Validar input
    es_valido, mensaje = validate_input(symbols, coordinates, charge, spin, basis)
    if not es_valido:
        print(f"Error de validación: {mensaje}", file=sys.stderr)
        sys.exit(1)
    
    # Ejecutar cálculo
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




