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


def main():
    """Función principal de la CLI"""
    parser = argparse.ArgumentParser(
        description="SCF Lite - Sistema minimalista para cálculos SCF usando PySCF"
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




