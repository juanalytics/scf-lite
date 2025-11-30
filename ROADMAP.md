# CORRECTED ROADMAP: SCF-Lite Application

## Overview
A minimal wrapper around PySCF/Psi4 that takes molecular input and returns electronic energy. The library handles all quantum chemistry; we handle input validation, execution, and output formatting.

---

## PHASE 1: Project Setup & Dependencies

### 1.1 Initialize Project Structure
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
├── setup.py (optional)
└── README.md
```

### 1.2 Define Dependencies
- **PySCF** (primary): `pyscf>=2.0.0`
- **Psi4** (optional alternative): `psi4`
- **NumPy**: For array operations
- **Pydantic** (optional): For input validation

---

## PHASE 2: Core Input Specification

### 2.1 Required Input Fields
```python
{
  "atoms": ["O", "H", "H"],           # List of atomic symbols
  "coords": [[0.0, 0.0, 0.0],         # Coordinates in Angstroms
             [0.957, 0.0, 0.0],
             [-0.24, 0.93, 0.0]],
  "charge": 0,                        # Molecular charge (default: 0)
  "multiplicity": 1,                  # Spin multiplicity (default: 1)
  "basis": "6-31g",                   # Basis set (default: "6-31g")
  "method": "HF"                      # Method: "HF" or "DFT" (default: "HF")
}
```

### 2.2 Input Validation Rules
- `len(atoms) == len(coords)` (must match)
- All coordinates must be 3D arrays `[x, y, z]`
- Atomic symbols must be valid (H, He, Li, ...)
- Charge must be integer
- Multiplicity must be positive integer (1 = singlet, 2 = doublet, ...)
- Basis set must be valid PySCF basis name
- Method must be "HF" or "DFT"

### 2.3 Charge & Multiplicity Calculation (if not provided)
- **Charge**: `charge = sum(atomic_numbers) - electrons`
  - If `electrons` provided: use it
  - If not: assume neutral (charge = 0)
- **Multiplicity**: `multiplicity = 2*S + 1` where S = total spin
  - Default: 1 (singlet) for even electrons, 2 (doublet) for odd
  - **Note**: User should provide this explicitly for accuracy

---

## PHASE 3: Core Calculation Engine

### 3.1 Molecule Construction (PySCF)
```python
from pyscf import gto, scf

# Build molecule object
mol = gto.Mole()
mol.atom = build_atom_string(atoms, coords)  # Format: "O 0 0 0; H 0.957 0 0; ..."
mol.basis = basis_set
mol.charge = charge
mol.spin = multiplicity - 1  # PySCF uses spin (0,1,2...) not multiplicity
mol.build()
```

### 3.2 SCF Calculation
```python
# For Hartree-Fock
if method == "HF":
    mf = scf.RHF(mol)  # Restricted HF (closed-shell)
    # For open-shell: scf.ROHF(mol) or scf.UHF(mol)
    
# For DFT
elif method == "DFT":
    mf = scf.RKS(mol)  # Restricted Kohn-Sham
    mf.xc = 'b3lyp'    # Functional
    
# Run calculation
energy = mf.kernel()  # Returns energy in Hartree
converged = mf.converged
```

### 3.3 Error Handling
- **Convergence failures**: Increase `mf.max_cycle` (default 50 → 200)
- **DIIS**: Enabled by default in PySCF
- **Damping**: Add `mf.damp = 0.5` if needed for difficult cases
- **Invalid inputs**: Raise clear error messages

---

## PHASE 4: API Structure

### 4.1 Main Function Signature
```python
def calculate_energy(
    atoms: List[str],
    coords: List[List[float]],
    charge: int = 0,
    multiplicity: int = 1,
    basis: str = "6-31g",
    method: str = "HF",
    functional: str = "b3lyp"  # Only for DFT
) -> Dict[str, Any]:
    """
    Calculate electronic energy using SCF.
    
    Returns:
        {
            "energy_au": float,      # Energy in Hartree
            "energy_ev": float,       # Energy in eV (optional)
            "converged": bool,
            "method": str,
            "basis": str,
            "n_iterations": int,     # Number of SCF cycles
            "charge": int,
            "multiplicity": int
        }
    """
```

### 4.2 Input/Output Format
- **Input**: JSON or Python dict
- **Output**: JSON-serializable dict
- **Units**: 
  - Coordinates: Angstroms (convert to Bohr if needed internally)
  - Energy: Hartree (primary), eV (converted: 1 Ha = 27.2114 eV)

---

## PHASE 5: Implementation Order

### Step 1: Minimal Working Example
- Hardcode a simple molecule (H₂O)
- Get PySCF working end-to-end
- Return energy value

### Step 2: Input Parsing
- Accept atoms + coords as input
- Validate input format
- Convert to PySCF format

### Step 3: Method Selection
- Implement HF path
- Implement DFT path
- Add basis set selection

### Step 4: Output Formatting
- Structure return dictionary
- Add unit conversions
- Include metadata (method, basis, convergence)

### Step 5: Error Handling
- Validate all inputs
- Handle convergence failures gracefully
- Provide informative error messages

### Step 6: Testing
- Test with known molecules (H₂O, H₂, CH₄)
- Verify energies match literature values
- Test edge cases (invalid input, convergence failures)

---

## PHASE 6: Optional Enhancements (Post-MVP)

### 6.1 Input Format Support
- **XYZ file parser**: Read standard XYZ format
- **MOL file parser**: Read MOL format
- **Inline string**: Accept XYZ as string

### 6.2 Geometry Preprocessing
- **Center of mass**: Optionally center molecule (doesn't affect energy)
- **Unit conversion**: Accept Bohr or Angstroms

### 6.3 Output Enhancements
- **Energy breakdown**: Electronic + nuclear repulsion
- **Orbital energies**: HOMO/LUMO (if requested)
- **Dipole moment**: (if available from method)

### 6.4 Visualization (Optional)
- Simple 3D structure plot (matplotlib/plotly)
- Molecular orbital visualization

---

## PHASE 7: What NOT to Implement

**Explicitly excluded** (library handles these):
- ❌ Manual integral calculation
- ❌ Building S, H, F matrices manually
- ❌ SCF convergence algorithm implementation
- ❌ Basis set construction
- ❌ Post-HF methods (MP2, CCSD, etc.)
- ❌ Custom functionals
- ❌ Symmetry handling
- ❌ Geometry optimization
- ❌ Frequency calculations

---

## TECHNICAL NOTES

### Charge Calculation Correction
- **Correct**: `charge = total_protons - electrons`
- **Incorrect**: `charge = electrons - total_protons` (from original roadmap)

### Multiplicity vs Spin
- PySCF uses `spin` (0, 1, 2, ...) = unpaired electrons
- Standard chemistry uses `multiplicity` (1, 2, 3, ...) = 2S + 1
- **Conversion**: `spin = multiplicity - 1`

### Default Method Selection
- **HF**: Fast, reliable, good for testing
- **DFT/B3LYP**: More accurate, standard for production
- **Basis sets**: 
  - `STO-3G`: Minimal, fast
  - `6-31G`: Standard, balanced
  - `6-31G(d)`: With polarization, better accuracy

### Convergence Strategy
1. Default settings (usually works)
2. Increase `max_cycle` to 200
3. Add damping: `mf.damp = 0.5`
4. Use level shifting if needed: `mf.level_shift = 0.3`
5. If still fails, return error with diagnostic info

---

## DEVELOPMENT WORKFLOW

1. **Start with PySCF tutorial**: Get familiar with basic usage
2. **Build minimal example**: H₂O molecule, hardcoded
3. **Generalize input**: Accept atoms + coords
4. **Add validation**: Ensure robust error handling
5. **Test with multiple molecules**: H₂, H₂O, CH₄, NH₃
6. **Add method selection**: HF vs DFT
7. **Polish output**: Format, units, metadata
8. **Document**: README with examples

---

## EXPECTED OUTPUT FORMAT

```json
{
  "energy_au": -76.027,
  "energy_ev": -2070.8,
  "converged": true,
  "method": "HF",
  "basis": "6-31G",
  "functional": null,
  "n_iterations": 12,
  "charge": 0,
  "multiplicity": 1,
  "n_atoms": 3,
  "n_electrons": 10
}
```

---

## ERROR RESPONSE FORMAT

```json
{
  "error": true,
  "error_type": "convergence_failure",
  "message": "SCF calculation did not converge after 200 cycles",
  "energy_au": null,
  "converged": false
}
```

