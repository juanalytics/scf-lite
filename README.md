# SCF-Lite

A minimal wrapper around PySCF for electronic energy calculations. This application takes molecular input (atoms, coordinates) and returns the electronic energy using quantum chemistry libraries.

## Features

- Simple API for SCF (Self-Consistent Field) calculations
- Support for Hartree-Fock (HF) and DFT methods
- Input validation and error handling
- Clean, structured output format

## Installation

### 1. Create and activate virtual environment

**Windows (PowerShell):**
```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
```

**Windows (Command Prompt):**
```cmd
python -m venv venv
venv\Scripts\activate.bat
```

**Linux/Mac:**
```bash
python -m venv venv
source venv/bin/activate
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

## Usage

### Minimal Example

```python
from src.scf_lite.calculator import calculate_h2o_energy

energy = calculate_h2o_energy()
print(f"Energy: {energy:.6f} Hartree")
```

### Running the test

```bash
python -m src.scf_lite.calculator
```

## Project Structure

```
scf-lite/
├── src/
│   └── scf_lite/
│       ├── __init__.py
│       ├── calculator.py      # Main SCF calculation logic
│       ├── input_validator.py # Input validation
│       └── output_formatter.py # Result formatting
├── tests/
│   └── test_calculator.py
├── requirements.txt
├── ROADMAP.md
└── README.md
```

## Development Status

🚧 **In Development** - Following the roadmap in `ROADMAP.md`

Current phase: Step 1 - Minimal Working Example

## Requirements

- Python 3.7+
- PySCF >= 2.0.0
- NumPy >= 1.20.0

## License

[To be determined]

