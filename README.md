# SCF Lite

Sistema minimalista para cálculos SCF (Self-Consistent Field) usando PySCF.

---

## 1. Setup: entorno, librerías y uso básico

### 1.1 Entorno recomendado

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
> Una opción cómoda es usar WSL (Ubuntu) y trabajar sobre el repo en `/mnt/c/...`.

### 1.2 Librerías principales

- **PySCF**: motor de química cuántica. Calcula integrales electrónicas y ejecuta SCF (Hartree–Fock, UHF).  
- **matplotlib**: genera las gráficas de energía vs distancia.  
- **scipy (opcional)**: se usa, si está disponible, para ajustar curvas tipo Morse a los escaneos.

### 1.3 Uso desde línea de comandos

#### Cálculo rápido con archivo JSON o XYZ

```bash
python -m scf_lite.cli -f examples/water.json --format json
```

o bien, usando un archivo `water.xyz`:

```bash
python -m scf_lite.cli -f examples/water.xyz --basis sto-3g --format json
```

#### Cálculo rápido con argumentos directos

```bash
python -m scf_lite.cli \
  -s O H H \
  -c 0.0 0.0 0.0  0.0 -0.757 0.587  0.0 0.757 0.587 \
  --charge 0 --spin 0 --basis sto-3g --format json
```

#### Guardar resultado en archivo

```bash
python -m scf_lite.cli -f examples/water.json -o resultado.json
```

### 1.4 Escaneos y gráficas incorporadas

La CLI incluye modos de escaneo con gráficos interactivos (requiere `matplotlib`):

- **H₂ (enlace H–H):**
  ```bash
  python -m scf_lite.cli --scan-h2
  ```

- **H₂O (un enlace O–H):**
  ```bash
  python -m scf_lite.cli --scan-oh
  ```

- **Escaneo genérico de un enlace en tu propia molécula**  
  (índices de átomos 1‑based):
  ```bash
  python -m scf_lite.cli -f examples/co2.json --scan-bond 1 2
  ```
  Esto usa la geometría del archivo, escanea la distancia entre los átomos 1 y 2
  y muestra tabla + curva + ajuste tipo Morse.

En todos los escaneos se muestra:
- Una tabla con `R` (distancia), `E` (Hartree) y `ΔE` (kcal/mol), marcando el mínimo.  
- Una gráfica con los puntos SCF, el mínimo y una curva ajustada tipo Morse.

### 1.5 Uso desde Python

```python
from scf_lite import calculate_scf, validate_input

# Definir molécula (agua)
symbols = ["O", "H", "H"]
coordinates = [
    [0.0, 0.0, 0.0],
    [0.0, -0.757, 0.587],
    [0.0,  0.757, 0.587],
]

es_valido, mensaje = validate_input(symbols, coordinates)
if es_valido:
    resultados = calculate_scf(symbols, coordinates, charge=0, spin=0, basis="sto-3g")
    print(resultados)
```

### 1.6 Formato de archivo de entrada

#### Formato JSON estándar

```json
{
  "symbols": ["O", "H", "H"],
  "coordinates": [[0.0, 0.0, 0.0], [0.0, -0.757, 0.587], [0.0, 0.757, 0.587]],
  "charge": 0,
  "spin": 0,
  "basis": "sto-3g"
}
```

- `symbols`: lista de símbolos químicos (obligatorio).  
- `coordinates`: lista de coordenadas 3D en Å (obligatorio).  
- `charge`: carga total (opcional, por defecto 0).  
- `spin`: \(N_\alpha - N_\beta\) (opcional, por defecto 0 → singlete/RHF).  
- `basis`: base de cálculo (opcional, por defecto `"sto-3g"`).  
- `name`: nombre descriptivo de la molécula (opcional, solo informativo).

#### Formato XYZ

```text
3
water
O  0.0000  0.0000  0.0000
H  0.0000 -0.7570  0.5870
H  0.0000  0.7570  0.5870
```

- Primera línea: número de átomos.  
- Segunda línea: comentario (se ignora).  
- Resto: líneas `Símbolo x y z`.  
- Por defecto se asume `charge = 0`, `spin = 0`, `basis = "sto-3g"`, que puedes
  sobreescribir con `--charge`, `--spin`, `--basis`.

### 1.7 Ejemplos incluidos

En `examples/`:
- `water.json` / `water.xyz`: H₂O neutra.  
- `hydrogen.json`: H₂.  
- `methane.json`: CH₄.  
- `co2.json` / `co2.xyz`: CO₂ lineal.  
- `nh3.json` / `nh3.xyz`: NH₃ (pirámide trigonal).  
- `oh_radical.json`: radical OH (spin 1).  
- `h2plus.json`: ion H₂⁺ (carga +1, spin 1).  
- `benzene.json`: benceno C₆H₆ (estructura plana).

---

## 2. Física y química que implementa SCF Lite

Esta sección explica **qué problema físico resuelve** el código, sin entrar aún
en detalles de implementación.

### 2.1 Modelo de la molécula

- Usamos la **aproximación de Born–Oppenheimer**:
  - Los **núcleos** se consideran fijos en las coordenadas que proporcionas.  
  - Solo resolvemos el problema de los **electrones** en el campo de esos núcleos.
- Cada cálculo toma:
  - Símbolos atómicos (identidad de cada núcleo).  
  - Coordenadas en Å.  
  - Carga total.  
  - Spin total \(S_z = N_\alpha - N_\beta\).

### 2.2 Bases atómicas (STO‑3G y otras)

- Los orbitales electrónicos no se representan como funciones continuas arbitrarias,
  sino como combinaciones lineales de **funciones base centradas en los átomos**.  
- La base por defecto es **STO‑3G**:
  - Es una base **mínima**: un conjunto pequeño de funciones por átomo.  
  - Es rápida y suficiente para ver tendencias y hacer demos.  
- También se aceptan otras bases estándar (`"6-31g"`, `"cc-pvdz"`, `"def2-svp"`),
  que dan energías más realistas a costa de más tiempo de cómputo.

### 2.3 Hartree–Fock y campo autoconsistente (SCF)

- El método principal es **Hartree–Fock (HF)**:
  - Aproxima la función de onda electrónica como un **determinante de Slater**
    construido con **orbitales moleculares**.  
  - Cada electrón se mueve en el **campo promedio** de los demás
    (de ahí “campo autoconsistente”).

- Tipos de SCF que usamos:
  - **RHF (Restricted HF)** si `spin = 0`:
    - Los pares de electrones con espines opuestos comparten el mismo orbital espacial.  
    - Apropiado para sistemas cerrados (singletes).  
  - **UHF (Unrestricted HF)** si `spin ≠ 0`:
    - Orbitales distintos para electrones α y β.  
    - Útil para radicales o especies con spin abierto.

### 2.4 Qué hace una iteración SCF (conceptualmente)

En cada ciclo SCF, internamente PySCF:

1. Parte de una densidad electrónica inicial (conjetura).  
2. Construye la **matriz de Fock**, que es el Hamiltoniano efectivo de un
   electrón en el campo promedio de todos los demás y los núcleos.  
3. Resuelve una ecuación de valores propios:
   - Obtiene nuevas **energías orbitales** y **orbitales moleculares**.  
4. Con esos orbitales construye una **nueva densidad electrónica**.  
5. Compara densidades nueva y anterior:
   - Si la diferencia es pequeña, dice que el SCF ha **convergido**.  
   - Si no, repite el ciclo.

Nuestro código se apoya en PySCF para estos pasos, y solo mide:
- Energía total final.  
- Número de iteraciones empleadas.  
- Tiempo de cómputo.

### 2.5 Energía que reporta el programa

La `energia` que devuelve SCF Lite es la energía total HF (incluye interacción
electrón–electrón y repulsión núcleo–núcleo tal como la maneja PySCF), en
**unidades de Hartree**:

- \(1\ \text{Hartree} \approx 27.2\ \text{eV} \approx 627.5\ \text{kcal/mol}\).

Lo físicamente relevante suelen ser las **diferencias de energía**:
- Al estirar un enlace, miramos \(\Delta E = E(R) - E_\text{min}\).  
- En los escaneos se convierten estas diferencias a **kcal/mol** para una
  interpretación más química.

### 2.6 Escaneos de energía vs distancia y potencial de Morse

En los modos `--scan-h2`, `--scan-oh` y `--scan-bond`:

1. Se genera una familia de geometrías variando una **distancia de enlace**:
   - H₂: distancia H–H.  
   - H₂O: un enlace O–H, manteniendo fijo el resto.  
   - `--scan-bond`: enlace entre átomos `I` y `J` de tu molécula.
2. Para cada distancia se ejecuta un cálculo HF independiente.  
3. Se construye la curva **E(R)** con los puntos SCF.  
4. Se busca el mínimo numérico → distancia de equilibrio aproximada.  
5. Se ajusta (cuando es posible) un **potencial tipo Morse**:

\[
E(R) \approx E_0 + D_e \left(1 - e^{-a (R - R_e)}\right)^2
\]

Donde:
- \(R_e\): distancia de equilibrio.  
- \(D_e\): profundidad del pozo (energía de enlace aproximada).  
- \(a\): parámetro de rigidez del enlace.  
- \(E_0\): energía de referencia a grandes distancias.

La gráfica muestra:
- Puntos discretos SCF.  
- Mínimo marcado con un punto rojo y línea vertical.  
- Curva suave de Morse ajustada y un recuadro con la ecuación y los parámetros.

---

## 3. Cómo se conecta el código con la física

Esta sección describe, módulo por módulo, qué hace cada parte del código y
cómo se relaciona con la física explicada arriba.

### 3.1 `scf_lite.calculator.calculate_scf`

Responsable de **armar y ejecutar el cálculo SCF**:

- Recibe:
  - `symbols`, `coordinates`, `charge`, `spin`, `basis`.  
- Construye la descripción atómica para PySCF:
  - Crea cadenas del tipo `\"O 0.0 0.0 0.0\"`, `\"H 0.0 -0.757 0.587\"`, etc.  
- Crea el objeto molecular:
  - `mol = gto.M(atom=atom_string, basis=basis, charge=charge, spin=spin)`  
  - Esto fija el modelo de la molécula (Sección 2.1 y 2.2).
- Selecciona el método SCF:
  - `RHF` si `spin == 0` (Sección 2.3, singlete).  
  - `UHF` si `spin != 0` (radicales o spins abiertos).
- Cuenta iteraciones mediante un `callback`:
  - Cada vez que PySCF completa un ciclo SCF llama al callback, que suma 1.  
  - Así obtenemos un conteo realista de iteraciones (Sección 2.4).
- Mide el tiempo de cómputo con `time.perf_counter()`.  
- Llama a `mf.kernel()`:
  - PySCF ejecuta todo el procedimiento SCF y devuelve la energía total HF.
- Devuelve un diccionario con:
  - `energia`, `convergio`, `iteraciones`, `metodo`, `spin`, `charge`,
    `basis`, `natom`, `nelec`, `tiempo_segundos`.

Relación con la física:
- Es el punto donde la **aplicación se conecta con PySCF** para hacer
  Hartree–Fock siguiendo el esquema de la Sección 2.3–2.5.

### 3.2 `scf_lite.input_validator`

Encargado de que el input tenga **sentido químico básico**:

- `validate_input(...)`:
  - Verifica que `len(symbols) == len(coordinates)`.  
  - Comprueba que los símbolos estén en una lista de elementos soportados
    (hasta Ca, para mantenerlo simple).  
  - Revisa que cada coordenada tenga 3 componentes numéricos.  
  - Revisa que la base elegida esté en un conjunto permitido
    (`"sto-3g"`, `"6-31g"`, `"cc-pvdz"`, `"def2-svp"`).
- Carga de archivos:
  - `_load_json_input`:
    - Lee JSON con `symbols`, `coordinates`, `charge`, `spin`, `basis`.  
  - `_load_xyz_input`:
    - Lee `.xyz`, interpreta el formato estándar y construye un diccionario
      compatible con el JSON.  
  - `load_input_file(filepath)`:
    - Decide según la extensión (`.json` o `.xyz`) qué cargador usar.

Relación con la física:
- Garantiza que la geometría y los elementos son razonables antes de intentar
  construir el objeto `mol` de PySCF.  
- Ayuda a evitar errores numéricos obvios (como vectores mal formados o
  elementos inventados).

### 3.3 `scf_lite.output_formatter`

Separa la parte de **formateo de salida**:

- `format_results(resultados, formato="dict", output_file=None)`:
  - Toma el diccionario completo devuelto por `calculate_scf`.  
  - Construye una salida simple con:
    - `energia`, `convergio`, `iteraciones`, `metodo`,
      y opcionalmente `tiempo_segundos`.  
  - Puede devolver un `dict` de Python o un string JSON.  
  - Si se da `output_file`, guarda el JSON en disco.

Relación con la física:
- No cambia el contenido físico; solo lo empaqueta para que sea fácil de leer,
  guardar o pasar a otros programas (por ejemplo, scripts de graficado).

### 3.4 `scf_lite.cli`

Proporciona la **interfaz de línea de comandos** y los modos de escaneo:

- Parseo de argumentos:
  - `-f/--file`: archivo de entrada (`.json` o `.xyz`).  
  - `-s/--symbols` y `-c/--coordinates`: entrada manual desde la CLI.  
  - `--charge`, `--spin`, `--basis`: parámetros físicos del sistema.  
  - `--format`: `dict` o `json` para la salida.  
  - `-o/--output`: archivo donde guardar el resultado.  
  - Modos especiales:
    - `--scan-h2`: escaneo preconfigurado para H₂.  
    - `--scan-oh`: escaneo preconfigurado para un O–H en H₂O.  
    - `--scan-bond I J`: escaneo genérico de un enlace entre átomos `I` y `J`
      usando la geometría del usuario (índices 1‑based).

- Flujo:
  1. Carga el input desde archivo o argumentos.  
  2. Si se pidió `--scan-h2` o `--scan-oh`, ejecuta los escaneos predefinidos.  
  3. Si se pidió `--scan-bond I J`:
     - Valida el input.  
     - Llama a `run_bond_scan(...)` con la geometría del usuario.  
  4. Si no hay modo de escaneo, hace un cálculo SCF único:
     - Llama a `calculate_scf(...)`.  
     - Formatea y muestra el resultado con `format_results(...)`.

Relación con la física:
- Es la **capa de orquestación**:
  - Decide si vamos a hacer un **cálculo único** (energía de una geometría)
    o un **escaneo de energía vs distancia** (Sección 2.6).  
  - En los escaneos, controla cómo se generan las geometrías al variar un
    **enlace concreto**.

### 3.5 Funciones de escaneo: `run_h2_scan`, `run_oh_scan`, `run_bond_scan`

Estas funciones implementan la lógica de los escaneos de la Sección 2.6:

- **`run_h2_scan`**:
  - Construye una familia de geometrías H₂ cambiando la distancia H–H.  
  - Llama a `calculate_scf` para cada distancia.  
  - Encuentra el mínimo, imprime la tabla y dibuja la curva + ajuste de Morse.

- **`run_oh_scan`**:
  - Usa una geometría fija de agua.  
  - Escanea solo uno de los enlaces O–H, manteniendo fijo el resto.  
  - Igual que `run_h2_scan`, produce tabla y gráfica con ajuste tipo Morse.

- **`run_bond_scan`**:
  - Trabaja con **cualquier molécula** que des (JSON/XYZ o entrada manual).  
  - Identifica las posiciones de los átomos `i` y `j`.  
  - Calcula la distancia inicial \(R_0\) y la dirección del enlace.  
  - Genera nuevas geometrías variando `j` entre \(0.7 R_0\) y \(1.3 R_0\),
    dejando `i` fijo.  
  - Para cada geometría:
    - Ejecuta HF (RHF o UHF, según el `spin`).  
  - Obtiene la curva E(R), busca el mínimo y ajusta un potencial de Morse
    cuando es posible.  
  - Imprime tabla y dibuja la gráfica con anotaciones.

En conjunto, estas funciones muestran de forma visual la **relación entre
geometría y energía** que se deriva del modelo de Hartree–Fock y conectan
directamente con la interpretación química de la fuerza y estabilidad del
enlace.


