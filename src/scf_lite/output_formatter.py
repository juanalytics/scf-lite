"""
Módulo para formatear y exportar resultados
"""

import json
from typing import Dict, Any


def format_results(
    resultados: Dict[str, Any],
    formato: str = "dict",
    output_file: str = None
) -> Any:
    """
    Formatea los resultados según el formato solicitado.
    
    Args:
        resultados: Diccionario con los resultados del cálculo
        formato: "dict" para diccionario Python, "json" para string JSON
        output_file: Si se proporciona, guarda el resultado en este archivo
    
    Returns:
        Diccionario o string JSON según formato
    """
    # Crear diccionario con solo los campos esenciales para salida simple
    output = {
        "energia": resultados["energia"],
        "convergio": resultados["convergio"],
        "iteraciones": resultados["iteraciones"]
    }
    
    # Agregar campos opcionales si están disponibles
    if "metodo" in resultados:
        output["metodo"] = resultados["metodo"]
    
    if formato == "json":
        json_str = json.dumps(output, indent=2, ensure_ascii=False)
        
        if output_file:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(json_str)
        
        return json_str
    
    # Formato dict
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(output, f, indent=2, ensure_ascii=False)
    
    return output
