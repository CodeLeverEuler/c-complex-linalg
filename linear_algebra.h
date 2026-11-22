#pragma once
#include <complex.h>
#include <stddef.h>

/* =========================================================================
 * Memory Management
 * ========================================================================= */

/**
 * Allocates a flat row-major complex matrix of size (rows x columns).
 * Returns NULL on allocation failure. Caller must free with free_matrix().
 */
double complex* create_matrix(size_t rows, size_t columns);

/**
 * Frees a matrix allocated by this library. Safely ignores NULL pointers.
 */
void free_matrix(double complex* matrix);


/* =========================================================================
 * Basic Arithmetic & Matrix Operations
 * ========================================================================= */

/**
 * Computes element-wise addition (A + B) for two (rows x columns) matrices.
 * Caller owns and must free the returned matrix.
 */
double complex* matrix_addition(double complex* a, double complex* b, size_t rows, size_t columns);

/**
 * Computes element-wise subtraction (A - B) for two (rows x columns) matrices.
 * Caller owns and must free the returned matrix.
 */
double complex* matrix_subtraction(double complex* a, double complex* b, size_t rows, size_t columns);

/**
 * Multiplies matrix A (rows_1 x columns_1) by matrix B (columns_1 x columns_2).
 * Returns a newly allocated (rows_1 x columns_2) matrix, or NULL on invalid input.
 */
double complex* matrix_multiplication(double complex* a, double complex* b, size_t rows_1, size_t columns_1, size_t columns_2);

/**
 * Scales every element of matrix A (rows x columns) by a real scalar x.
 * Caller owns and must free the returned matrix.
 */
double complex* matrix_scalarprod(double complex* a, double x, size_t rows, size_t columns);

/**
 * Raises a square matrix A (rows == columns) to an integer power >= 0.
 * Returns the identity matrix when pow == 0. Caller must free the result.
 */
double complex* matrix_pow(double complex* a, size_t pow, size_t rows, size_t columns);


/* =========================================================================
 * Structural Transforms & Products
 * ========================================================================= */

/**
 * Returns the transpose (A^T) of size (columns x rows). Caller must free.
 */
double complex* matrix_transpose(double complex* a, size_t rows, size_t columns);

/**
 * Returns the element-wise complex conjugate (A*) of size (rows x columns).
 * Caller must free the returned matrix.
 */
double complex* complex_conjugate(double complex* a, size_t rows, size_t columns);

/**
 * Returns the conjugate transpose / Hermitian adjoint (A^H) of size (columns x rows).
 * Caller must free the returned matrix.
 */
double complex* hermitian_conjugate(double complex* a, size_t rows, size_t columns);

/**
 * Generates a square identity matrix of size (rows x columns).
 * Returns NULL if rows != columns or rows == 0. Caller must free.
 */
double complex* kronecker_delta(size_t rows, size_t columns);

/**
 * Computes the Hermitian inner product <A, B> = sum(conj(A_ij) * B_ij).
 * Returns 0 on NULL or zero-dimension inputs.
 */
double complex inner_product(double complex* a, double complex* b, size_t rows, size_t columns);

/**
 * Computes the Kronecker (tensor) product (A ⊗ B) of dimensions
 * ((rows_a * rows_b) x (columns_a * columns_b)). Caller must free.
 */
double complex* tensor_product(double complex* a, double complex* b, size_t rows_a, size_t columns_a, size_t rows_b, size_t columns_b);


/* =========================================================================
 * Matrix Metrics & Invariants
 * ========================================================================= */

/**
 * Computes the determinant of a square matrix using Gaussian elimination
 * with partial pivoting. Returns 0 if singular or non-square.
 */
double complex matrix_determinant(double complex* a, size_t rows, size_t columns);

/**
 * Computes the trace (sum of main diagonal elements) of a square matrix.
 * Returns 0 if non-square or NULL.
 */
double complex matrix_trace(double complex* a, size_t rows, size_t columns);

/**
 * Computes the Frobenius norm sqrt(sum(|A_ij|^2)) of matrix A.
 */
double complex matrix_frobenius_norm(double complex* a, size_t rows, size_t columns);


/* =========================================================================
 * Decompositions, Solvers & Eigenvalues
 * ========================================================================= */

/**
 * Computes the inverse (A^-1) of a square matrix via Gauss-Jordan elimination
 * with partial pivoting. Returns NULL if singular or non-square. Caller must free.
 */
double complex* matrix_inverse(double complex* a, size_t rows, size_t columns);

/**
 * Computes the upper-triangular R factor of QR decomposition using
 * complex Householder reflections. Caller must free.
 */
double complex* matrix_rdecomposition(double complex* a, size_t rows, size_t columns);

/**
 * Computes the unitary Q factor of QR decomposition using
 * complex Householder reflections. Caller must free.
 */
double complex* matrix_qdecomposition(double complex* a, size_t rows, size_t columns);

/**
 * Approximates all eigenvalues of a square matrix using iterative QR decomposition.
 * Returns a newly allocated 1D array of length 'rows'. Caller must free.
 */
double complex* matrix_eigenvalue(double complex* a, size_t rows, size_t columns, size_t max_iterations, double tolerance);

/**
 * Solves the complex linear system A * x = b using Gaussian elimination with
 * partial pivoting. Returns solution vector x of size (rows x 1), or NULL if singular.
 */
double complex* matrix_solver(double complex* a, double complex* b, size_t rows, size_t columns);
