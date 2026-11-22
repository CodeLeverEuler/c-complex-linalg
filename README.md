# c-complex-linalg

Dependency-free complex linear algebra library in C99 featuring advanced matrix calculations and applications in AC circuit analysis and quantum entropy.

Built from scratch using native `<complex.h>` (`double complex`) and a flat 1D row-major contiguous memory layout.

---

## Features

* **Memory Management:** Flat 1D row-major matrix allocation (`create_matrix`, `free_matrix`) with explicit caller-owned memory contracts.
* **Core Arithmetic:** Matrix addition, subtraction, multiplication, real scalar multiplication, and integer matrix exponentiation (`matrix_pow`).
* **Transforms & Products:** Transpose ($A^T$), element-wise complex conjugate (\(A^*\)), Hermitian adjoint (\(A^H\)), Kronecker identity generator, Hermitian inner product (\(\langle A, B \rangle\)), and Kronecker tensor product (\(A \otimes B\)).
* **Matrix Invariants:** Determinant (via partial pivoting), trace, and Frobenius norm.
* **Decompositions & Solvers:**
  * Complex Gauss-Jordan matrix inversion (\(A^{-1}\)) with partial pivoting
  * Complex Householder QR decomposition (`matrix_qdecomposition`, `matrix_rdecomposition`)
  * Iterative QR eigenvalue solver (`matrix_eigenvalue`)
  * Complex linear system solver (\(Ax = b\)) via Gaussian elimination with partial pivoting (`matrix_solver`)

---

## Repository Structure

```text
├── linear_algebra.h                 # Public API & Doxygen documentation
├── linear_algebra.c                 # Core numerical implementations
├── Examples/                        # 15 standalone unit examples covering every API function
│   ├── matrix_add.c
│   ├── matrix_eigenvalues.c
│   ├── matrix_qrdecomp.c
│   ├── matrix_solver.c
│   └── ...
└── Applications/                    # Real-world physics & ECE domain solvers
    ├── von_neumann_entropy.c        # Quantum density matrix & Von Neumann entropy calculator
    └── ac_mesh_solver.c             # Multi-loop AC circuit KVL & complex power analyzer
```

---

## Domain Applications

### 1. Quantum Mixed-State Density Matrix & Von Neumann Entropy Calculator
Constructs a density matrix $\rho = \sum_i p_i |\psi_i\rangle\langle\psi_i|$ from an ensemble of pure quantum states, verifies Hermiticity ($\rho = \rho^H$), extracts eigenvalues via iterative QR decomposition, and computes Von Neumann entropy in bits:
$$S(\rho) = -\sum_i \lambda_i \log_2(\lambda_i)$$

### 2. AC Circuit Mesh & Complex Power Analyzer
Models multi-loop alternating-current circuits using Kirchhoff's Voltage Law in matrix form ($Z \cdot I = V$) to solve for complex loop currents and calculates total complex power ($S = V \cdot I^*$), separating active power (Watts) and reactive power (VAR).

---

## Build & Run

Compile any example or application with `gcc` (linking the standard math library `-lm`):

```bash
# Compile an API example
gcc linear_algebra.c Examples/matrix_solver.c -o solver_demo -lm
./solver_demo

# Compile the AC Mesh & Power Analyzer
gcc linear_algebra.c Applications/ac_mesh_solver.c -o ac_solver -lm
./ac_solver

# Compile the Von Neumann Entropy Calculator
gcc linear_algebra.c Applications/von_neumann_entropy.c -o entropy_calc -lm
./entropy_calc
```
