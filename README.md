# SCF Lite

Sistema minimalista para cálculos SCF (Self-Consistent Field) usando PySCF.

## Instalación

1. Crear entorno virtual (recomendado):
```bash
python -m venv venv
source venv/bin/activate  # En Windows: venv\Scripts\activate
```

2. Instalar dependencias:
```bash
pip install -r requirements.txt
```

## Uso

### Desde línea de comandos

#### Usando archivo de input:
```bash
python -m scf_lite.cli -f examples/water.json
```

#### Usando argumentos directos:
```bash
python -m scf_lite.cli -s O H H -c 0.0 0.0 0.0 0.0 -0.757 0.587 0.0 0.757 0.587
```

#### Guardar resultado en archivo:
```bash
python -m scf_lite.cli -f examples/water.json -o resultado.json
```

### Desde Python

```python
from scf_lite import calculate_scf, validate_input

# Validar input
symbols = ["O", "H", "H"]
coordinates = [[0.0, 0.0, 0.0], [0.0, -0.757, 0.587], [0.0, 0.757, 0.587]]

es_valido, mensaje = validate_input(symbols, coordinates)
if es_valido:
    resultados = calculate_scf(symbols, coordinates)
    print(resultados)
```

## Formato de archivo de input

El archivo JSON debe tener el siguiente formato:

```json
{
  "symbols": ["O", "H", "H"],
  "coordinates": [[0.0, 0.0, 0.0], [0.0, -0.757, 0.587], [0.0, 0.757, 0.587]],
  "charge": 0,
  "spin": 0,
  "basis": "sto-3g"
}
```

Campos:
- `symbols`: Lista de símbolos químicos (requerido)
- `coordinates`: Lista de coordenadas 3D en Angstrom (requerido)
- `charge`: Carga total (opcional, default: 0)
- `spin`: Spin total (opcional, default: 0, singlete/RHF)
- `basis`: Base a usar (opcional, default: "sto-3g")

## Ejemplos predefinidos

En la carpeta `examples/` encontrarás archivos de ejemplo:
- `water.json`: Molécula de agua
- `hydrogen.json`: Molécula de hidrógeno
- `methane.json`: Molécula de metano

## Salida

El sistema devuelve un diccionario/JSON con:
- `energia`: Energía total en Hartree
- `convergio`: Boolean indicando si convergió
- `iteraciones`: Número de iteraciones SCF
- `metodo`: "RHF" o "UHF"




