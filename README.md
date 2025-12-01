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

> Nota para Windows: PySCF se instala más fácilmente en Linux.  
> Una opción muy cómoda es usar WSL (Ubuntu) y trabajar sobre el repo en `/mnt/c/...`.

## Uso

### Desde línea de comandos

#### Usando archivo de input JSON o XYZ:
```bash
python -m scf_lite.cli -f examples/water.json
```
o bien, si tienes un `water.xyz`:
```bash
python -m scf_lite.cli -f examples/water.xyz --basis sto-3g
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

### Formato JSON estándar

El archivo JSON debe tener el siguiente formato mínimo:

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

También se acepta el campo opcional:
- `name`: Nombre descriptivo de la molécula (solo informativo)

### Formato XYZ

También puedes usar archivos `.xyz` sencillos:
```text
3
water
O  0.0000  0.0000  0.0000
H  0.0000 -0.7570  0.5870
H  0.0000  0.7570  0.5870
```

- La primera línea es el número de átomos.
- La segunda es un comentario (se ignora).
- El resto son líneas `Símbolo x y z`.
- Por defecto, el sistema asume `charge = 0`, `spin = 0`, `basis = "sto-3g"`
  (puedes sobreescribirlos con `--charge`, `--spin`, `--basis` en la CLI).

## Ejemplos predefinidos

En la carpeta `examples/` encontrarás archivos de ejemplo:
- `water.json` / `water.xyz`: Molécula de agua neutra
- `hydrogen.json`: Molécula de H₂
- `methane.json`: Molécula de CH₄
- `co2.json` / `co2.xyz`: Molécula lineal de CO₂
- `nh3.json` / `nh3.xyz`: Molécula de NH₃ (pirámide trigonal)
- `oh_radical.json`: Radical OH (spin 1)
- `h2plus.json`: Ion H₂⁺ (carga +1, spin 1)
- `benzene.json`: Molécula de benceno C₆H₆ (estructura plana)

## Salida

El sistema devuelve un diccionario/JSON con:
- `energia`: Energía total en Hartree
- `convergio`: Boolean indicando si convergió
- `iteraciones`: Número de iteraciones SCF
- `metodo`: "RHF" o "UHF"
 - `tiempo_segundos`: Tiempo total de SCF en segundos

## Escaneos y gráficas incorporadas

La CLI incluye algunos modos de escaneo con gráficos interactivos (requiere `matplotlib`):

- Escaneo del enlace H–H en H₂:
  ```bash
  python -m scf_lite.cli --scan-h2
  ```
- Escaneo del enlace O–H en H₂O (manteniendo fijo el resto de la geometría):
  ```bash
  python -m scf_lite.cli --scan-oh
  ```

Ambos modos muestran:
- Una tabla con `R` (distancia), `E` (Hartree) y `ΔE` (kcal/mol), con el mínimo marcado.
- Una curva de energía y un ajuste tipo Morse anotado en la figura.




