#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"linear_algebra.h"
/* =========================================================================
 * Memory Management
 * ========================================================================= */
double complex* create_matrix(size_t rows, size_t columns)
{
    double complex* matrix = (double complex*)malloc(rows*columns*sizeof(double complex));
    if(matrix==NULL||rows==0||columns==0) return NULL;
    return matrix;
}
void free_matrix(double complex* matrix)
{
    if(matrix==NULL) return ;
    free(matrix);
}
/* =========================================================================
 * Basic Arithmetic & Matrix Operations
 * ========================================================================= */
double complex* matrix_addition(double complex* a, double complex* b, size_t rows, size_t columns)
{
    double complex* result = create_matrix(rows,columns);
    if(result==NULL||rows==0||columns==0) return NULL;
    for(size_t i=0; i<rows; ++i)
    {
        for(size_t j=0; j<columns; ++j)
        {
            size_t idx = i*columns+j;
            result[idx]= a[idx]+b[idx];
        }
    }
return result;
}
double complex* matrix_subtraction(double complex* a, double complex* b, size_t rows, size_t columns)
{
    double complex* result = create_matrix(rows,columns);
    if(result==NULL||rows==0||columns==0) return NULL;
    for(size_t i=0; i<rows; ++i)
    {
        for(size_t j=0; j<columns; ++j)
        {
            size_t idx = i*columns+j;
            result[idx]= a[idx]-b[idx];
        }
    }
return result;
}
double complex* matrix_multiplication(double complex*a, double complex*b,size_t rows_1, size_t columns_1, size_t columns_2)
{
    if(a==NULL||b==NULL|| rows_1 == 0||columns_1==0||columns_2==0) return NULL;
    double complex* result = create_matrix(rows_1,columns_2);
    if(result==NULL) return NULL;
     for(size_t i=0; i<rows_1; ++i)
    {
        for(size_t j=0; j<columns_2; ++j)
        {
            size_t idx_result = i*columns_2+j;
            result[idx_result]= 0;
            for(size_t k=0; k<columns_1;++k)
            {
                size_t idx_a = i*columns_1+k;   //A_ik
                size_t idx_b = k*columns_2+j;   //B_kj
                result[idx_result] += a[idx_a]*b[idx_b];
            }
        }
    }
    return result;
}
double complex* matrix_scalarprod(double complex*a, double x, size_t rows,size_t columns)
{
    if(a==NULL||rows==0||columns==0) return NULL;
    double complex* result = create_matrix(rows, columns);
    if(result==NULL) return NULL;
    for(size_t i=0; i<rows; ++i)
    {
        for(size_t j=0; j<columns; ++j)
        {
            result[i*columns+j] = x * a[i*columns+j];
        }
    }
    return result;
}
double complex* matrix_pow(double complex*a,size_t pow, size_t rows, size_t columns)
{
    if(a==NULL||rows==0||columns!=rows) return NULL;
    if(pow ==0)
    return kronecker_delta(rows,columns);
    double complex* result = create_matrix(rows, columns);
    if(result==NULL) return NULL;
    //Initialization
    for(size_t i=0; i<rows; ++i)
            {
                for(size_t j=0; j<columns;++j)
                result[i*columns+j] = a[i*columns+j];
            }
    for(size_t count = 1; count<pow;++count)
    {
        double complex* prod = matrix_multiplication(result,a,rows,columns,columns);
        if(prod==NULL)
        {
            free_matrix(result);
            return NULL;
        }
        free_matrix(result);
        result = prod;
    }
    return result;
}
/* =========================================================================
 * Structural Transforms & Products
 * ========================================================================= */
double complex* matrix_transpose(double complex* a,size_t rows, size_t columns)
{
    double complex* result = create_matrix(columns,rows);
    if(result==NULL||rows==0||columns==0) return NULL;
    for(size_t i =0; i<rows;++i)
    {
        for(size_t j=0; j<columns;++j)
        {
            size_t idx_original = i*columns+j;
            size_t idx_transpose = j*rows+i;
            result[idx_transpose] = a[idx_original];
        }
    }
    return result;
}
double complex* complex_conjugate(double complex* a,size_t rows,size_t columns)
{
    if(a==NULL||rows==0 ||columns ==0) return NULL;
    double complex* result = create_matrix(rows,columns);
    if(result == NULL) return NULL;
    for(size_t i =0; i<rows; ++i)
    {
        for(size_t j =0; j<columns; ++j)
        {
            result[i*columns+j]= conj(a[i*columns+j]);
        }
    }
    return result;
}
double complex* hermitian_conjugate(double complex* a,size_t rows,size_t columns)
{
    if(a==NULL||rows==0 ||columns ==0) return NULL;
    double complex* result = create_matrix(rows,columns);
    if(result == NULL) return NULL;
    result = complex_conjugate(a,rows,columns);
    result = matrix_transpose(result,rows,columns);
    return result;
}
double complex* kronecker_delta(size_t rows,size_t columns)
{
    if(rows!=columns||rows==0) return NULL;
    double complex* result = create_matrix(rows,columns);
    if(result==NULL) return NULL;
    for(size_t i=0;i<rows;++i)
    {
        for(size_t j=0; j< columns; ++j)
        {
            if(i==j)
            {
                result[i*columns+j] = 1;
            }
            else
            result[i*columns+j] = 0; 
        }
    }
    return result;
}
double complex inner_product(double complex* a, double complex* b, size_t rows, size_t columns) 
{
    if(a==NULL||b==NULL||rows ==0 || columns ==0) return 0;
    double complex result = 0;
    for(size_t i =0; i<rows;++i)
    {
        for(size_t j=0; j<columns;++j)
        {
            result += conj(a[i*columns+j])*b[i*columns+j];
        }
    }
    return result;
}
double complex* tensor_product(double complex*a, double complex* b, size_t rows_a, size_t columns_a, size_t rows_b, size_t columns_b)
{
    if (a == NULL || b == NULL || rows_a == 0 || columns_a == 0 || rows_b == 0 || columns_b == 0) return NULL;
    size_t rows = rows_a*rows_b; 
    size_t columns = columns_a*columns_b;
    double complex* result = create_matrix(rows, columns);
    if(result==NULL) return NULL;
    for(size_t i=0; i<rows_a; ++i )
    {
        for(size_t j=0; j<columns_a; ++j)
        {
            for(size_t k=0; k<rows_b; ++k)
            {
                for(size_t l=0;l<columns_b; ++l)
                {
                    size_t idx_result = (i*(rows_b)+k)*columns +(j*(columns_b)+l);
                    result[idx_result] = a[i*columns_a+j]*b[k*columns_b+l];
                }
            }
        }
    }
    return result;
}
/* =========================================================================
 * Matrix Metric and Invariants
 * ========================================================================= */
double complex matrix_determinant(double complex*a, size_t rows, size_t columns)
{
    if(a==NULL||rows!=columns||rows==0) return 0;
    double complex* matrix = create_matrix(rows,columns);
    if(matrix==NULL) return 0;
    for(size_t i =0; i<rows*columns; ++i)
        matrix[i] = a[i];
    double complex det = 1+0*I; int sign = 1;
    for(size_t i = 0; i<rows;++i)
    {
        size_t pivot_row = i;
        for(size_t j = i+1;j<rows;++j)
        {
            size_t idx_mat = j*columns +i;
            size_t idx_pivot = pivot_row*columns+i;
            if(cabs(matrix[idx_mat])>cabs(matrix[idx_pivot]))
            pivot_row=j;
        } 
        if(cabs(matrix[pivot_row*columns+i]) < 1e-9)
        {
            free_matrix(matrix);
            return 0;
        }
        if(pivot_row!=i)
        {
            for(size_t j = 0; j<columns;++j)
            {
                size_t idx_mat = i*columns +j;
                size_t idx_pivot = pivot_row*columns+j;
                double complex temp = matrix[idx_mat];
                matrix[idx_mat] = matrix[idx_pivot];
                matrix[idx_pivot]=temp;
            }
            sign = -sign;
        }
        for(size_t j = i+1;j<rows;++j)
        {
            double complex factor = matrix[j*columns+i]/matrix[i*columns+i];
            for(size_t k=i; k<columns;++k)
            {
                matrix[j*columns+k] -= factor*matrix[i*columns+k];
            }
        }
        det *= matrix[i*columns+i];
    }
    double complex result = det*sign;
    free_matrix(matrix); 
    return result;
}
double complex matrix_trace(double complex*a, size_t rows, size_t columns)
{
     if(a == NULL ||rows!=columns||rows==0) return 0;
     double complex result = 0;
     for(size_t i=0;i<rows; ++i)
     {
        result += a[i*columns+i]; 
     }
     return result;
}
double complex matrix_frobenius_norm(double complex*a, size_t rows, size_t columns)
{
   double result,norm=0;
   if(a == NULL ||columns==0||rows==0) return 0;
   for(size_t i =0; i<rows; ++i )
   {
        for(size_t j=0; j<columns; ++j)
        {
           norm += cabs(a[i*columns+j])*cabs(a[i*columns+j]);  
        }
   }
   result = sqrt(norm);
   return result;
}
/* =========================================================================
 * Decompositions and Solvers
 * ========================================================================= */
double complex* matrix_inverse(double complex*a,size_t rows, size_t columns)
{
    if(a==NULL||rows!=columns||rows==0) return NULL;
    double complex* result = create_matrix(rows,columns); 
    if(result==NULL) return NULL;
    size_t aug_columns = 2*columns;
    double complex* matrix_augmented = (double complex*)malloc(rows*aug_columns*sizeof(double complex));
    if (matrix_augmented == NULL)
    {
        free_matrix(result);
        return NULL;
    }
// Build augmented matrix [a | I]
    for(size_t i=0;i<rows;++i)
    {
        for(size_t j = 0; j<columns;++j)
        {
            matrix_augmented[i*aug_columns +j] = a[i*columns+j];
            if(i==j)
            matrix_augmented[i*aug_columns +(j+columns)] =1;
            else
            matrix_augmented[i*aug_columns +(j+columns)] =0;
        }   
    }
 // Gauss-Jordan elimination with partial pivoting
    for(size_t i=0;i<rows;++i)
    {
        
        size_t pivot_row = i;
        for(size_t j = i+1;j<rows;++j)
        {
            size_t idx_mat = j*aug_columns +i;
            size_t idx_pivot = pivot_row*aug_columns+i;
            if(cabs(matrix_augmented[idx_mat])>cabs(matrix_augmented[idx_pivot]))
            pivot_row=j;
        } 
        if(cabs(matrix_augmented[pivot_row*aug_columns+i]) < 1e-9)
        {
            free_matrix(matrix_augmented);
            free_matrix(result);
            return NULL;
        }
        if(pivot_row!=i)
        {
            for(size_t j = 0; j<aug_columns;++j)
            {
                size_t idx_mat = i*aug_columns +j;
                size_t idx_pivot = pivot_row*aug_columns+j;
                double complex temp = matrix_augmented[idx_mat];
                matrix_augmented[idx_mat] = matrix_augmented[idx_pivot];
                matrix_augmented[idx_pivot]=temp;
            }
        }
// Normalize pivot row and eliminate column entries above and below
        double complex pivot = matrix_augmented[i*aug_columns+i];
        for(size_t j =0; j<aug_columns;++j)
        matrix_augmented[i*aug_columns+j] /= pivot;
        for(size_t j = 0;j<rows;++j)
        {
            if(j!=i)
            {
                double complex factor = matrix_augmented[j*aug_columns+i]/matrix_augmented[i*aug_columns+i];
                for(size_t k=i; k<aug_columns;++k)
                matrix_augmented[j*aug_columns+k] -= factor*matrix_augmented[i*aug_columns+k];
            }
        }
    }
 // Extract inverted matrix from the right half [I | A^-1]
    for(size_t i=0; i<rows; ++i)
        {
            for(size_t j=0; j<columns;++j)
            {
                result[i*columns+j] = matrix_augmented[i*aug_columns+(j+columns)];
            }
        }
    free(matrix_augmented);
    return result;
}
double complex* matrix_rdecomposition(double complex*a, size_t rows, size_t columns)
{
    if(a == NULL ||columns==0||rows==0) return NULL;
    double complex* R = create_matrix(rows, columns);
    double complex* V = (double complex*)calloc(rows,sizeof(double complex));
    if(R==NULL || V==NULL)
    {
        free(R);
        free(V);
        return NULL;
    }
    for(size_t i=0; i<rows;++i)
    {
        for(size_t j=0; j<columns; ++j)
        R[i*columns+j] = a[i*columns+j];
    }   
    // Apply Householder reflections column-by-column to zero out sub-diagonal entries
    for(size_t i =0; i< rows-1 && i<columns;++i)
    {
        double norm =0;
        for(size_t j=i; j<rows; ++j)
        {
            norm += cabs(R[j*columns+i])*cabs(R[j*columns+i]);
        }
        norm=sqrt(norm);
        if(norm<1e-9) continue;

        for(size_t j=0; j<rows;++j)
        {
            V[j]= 0;
        }
        for(size_t j=i; j<rows;++j)
        {
            V[j]= R[j*columns+i];
        }
     // Phase-aligned reflection vector avoids catastrophic cancellation
        if(cabs(R[i*columns+i])>1e-9)
        {
            double complex phase = R[i*columns+i]/cabs(R[i*columns+i]);
            V[i]+=phase*norm;
        }
        else
        V[i] +=norm;
        double complex factor =0;
        for(size_t j=0; j<rows; ++j)
        { 
            factor += conj(V[j])*V[j];
        }
        if(cabs(factor)<1e-9)
        continue;
        for(size_t j=i;j<columns;++j)
        {
            double complex dot =0;
            for(size_t k=i; k<rows; ++k)
            dot += conj(V[k])*R[k*columns+j];
            for(size_t k=i; k<rows;++k)
            R[k*columns+j] -= 2 * V[k]*dot/factor; 
        }
    }
    free(V);
    return R;
}
double complex* matrix_qdecomposition(double complex*a, size_t rows, size_t columns)
{
    if(a == NULL ||columns==0||rows==0) return NULL;
    double complex* R = create_matrix(rows, columns);
    double complex* Q = kronecker_delta(rows, columns);
    double complex* V = (double complex*)calloc(rows,sizeof(double complex));
    if(R==NULL || V==NULL)
    {
        free(R);
        free(V);
        return NULL;
    }
    for(size_t i=0; i<rows;++i)
    {
        for(size_t j=0; j<columns; ++j)
        R[i*columns+j] = a[i*columns+j];
    }
// Accumulate Householder reflectors into Q while reducing R   
    for(size_t i =0; i< rows-1 && i<columns;++i)
    {
        double norm =0;
        for(size_t j=i; j<rows; ++j)
        {
            norm += cabs(R[j*columns+i])*cabs(R[j*columns+i]);
        }
        norm=sqrt(norm);
        if(norm<1e-9) continue;

        for(size_t j=0; j<rows;++j)
        {
            V[j]= 0;
        }
        for(size_t j=i; j<rows;++j)
        {
            V[j]= R[j*columns+i];
        }
        if(cabs(R[i*columns+i])>1e-9)
        {
            double complex phase = R[i*columns+i]/cabs(R[i*columns+i]);
            V[i]+=phase*norm;
        }
        else
        V[i] +=norm;
        double complex factor =0;
        for(size_t j=0; j<rows; ++j)
        { 
            factor += conj(V[j])*V[j];
        }
        if(cabs(factor)<1e-9)
        continue;
        for(size_t j=i;j<columns;++j)
        {
            double dot =0;
            for(size_t k=i; k<rows; ++k)
            dot += conj(V[k])*R[k*columns+j];
            for(size_t k=i; k<rows;++k)
            R[k*columns+j] -= 2 * V[k]*dot/factor; 
        }
        for(size_t j=0; j<rows; ++j)
        {
            double complex dot =0;
            for(size_t k=i; k<rows; ++k)
            dot += Q[j*rows+k]*V[k];
            for(size_t k=i; k<rows;++k)
            Q[j*rows+k] -= 2 * dot * conj(V[k])/factor;
        }
    }
    free_matrix(R);
    free(V);
    return Q;
}
double complex* matrix_eigenvalue(double complex* a, size_t rows, size_t columns, size_t max_loop, double tolerance)
{
    if (a == NULL || rows == 0 || rows != columns) return NULL;
    double complex* A_curr = create_matrix(rows, columns);
    if (A_curr == NULL) return NULL;
    for (size_t i = 0; i < rows * columns; ++i)
    {
        A_curr[i] = a[i];
    }
// Iterative QR algorithm: A_(k+1) = R_k * Q_k converges toward upper triangular Schur form
    for (size_t k = 0; k < max_loop; ++k)
    {
// 1. Compute Q and R 
        double complex* Q = matrix_qdecomposition(A_curr, rows, columns);
        double complex* R = matrix_rdecomposition(A_curr, rows, columns);

        if (Q == NULL || R == NULL)
        {
            free_matrix(Q);
            free_matrix(R);
            free_matrix(A_curr);
            return NULL;
        }
        double complex* A_next = matrix_multiplication(R, Q, rows, columns, columns);
        free_matrix(Q);
        free_matrix(R);
        if (A_next == NULL)
        {
            free_matrix(A_curr);
            return NULL;
        }
// Check Frobenius norm of sub-diagonal elements for convergence
        double off_diag_norm = 0.0;
        for (size_t i = 1; i < rows; ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                double val = cabs(A_next[i * columns + j]);
                off_diag_norm += val * val;
            }
        }

        free_matrix(A_curr);
        A_curr = A_next;

        if (sqrt(off_diag_norm) < tolerance)
        {
            break;
        }
    }

// 4. Extract the eigenvalues from the main diagonal of the converged Schur form
    double complex* result = (double complex*)malloc(rows * sizeof(double complex));
    if (result == NULL)
    {
        free_matrix(A_curr);
        return NULL;
    }

    for (size_t i = 0; i < rows; ++i)
    {
        result[i] = A_curr[i * columns + i];
    }

    free_matrix(A_curr);
    return result;
}
double complex* matrix_solver(double complex*a,double complex* b, size_t rows, size_t columns)
{
    if(a==NULL||b==NULL||rows!=columns||rows==0) return NULL;
    double complex* result = create_matrix(rows,1); 
    if(result==NULL) return NULL;
    size_t aug_columns = columns+1;
// Build augmented matrix [a | b]
    double complex* matrix_augmented = (double complex*)malloc(rows*aug_columns*sizeof(double complex));
    if (matrix_augmented == NULL)
    {
        free_matrix(result);
        return NULL;
    }

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < columns; ++j)
        {
            matrix_augmented[i * aug_columns + j] = a[i * columns + j];
        }
        matrix_augmented[i * aug_columns + columns] = b[i];
    } 
 // Forward elimination with partial pivoting
    for(size_t i=0;i<rows;++i)
    {
        
        size_t pivot_row = i;
        for(size_t j = i+1;j<rows;++j)
        {
            size_t idx_mat = j*aug_columns +i;
            size_t idx_pivot = pivot_row*aug_columns+i;
            if(cabs(matrix_augmented[idx_mat])>cabs(matrix_augmented[idx_pivot]))
            pivot_row=j;
        } 
        if(cabs(matrix_augmented[pivot_row*aug_columns+i]) < 1e-9)
        {
            free_matrix(matrix_augmented);
            free_matrix(result);
            return NULL;
        }
        if(pivot_row!=i)
        {
            for(size_t j = 0; j<aug_columns;++j)
            {
                size_t idx_mat = i*aug_columns +j;
                size_t idx_pivot = pivot_row*aug_columns+j;
                double complex temp = matrix_augmented[idx_mat];
                matrix_augmented[idx_mat] = matrix_augmented[idx_pivot];
                matrix_augmented[idx_pivot]=temp;
            }
        }
        for(size_t j = i+1;j<rows;++j)
        {
            if(j!=i)
            {
                double complex factor = matrix_augmented[j*aug_columns+i]/matrix_augmented[i*aug_columns+i];
                for(size_t k=i; k<aug_columns;++k)
                matrix_augmented[j*aug_columns+k] -= factor*matrix_augmented[i*aug_columns+k];
            }
        }
    }
// Back-substitution to solve for x
    for (size_t i = rows;i-->0;)
    {
        result[i]=matrix_augmented[i*aug_columns+columns];
        for(size_t j =i +1; j<rows;++j)
        {
            result[i]-=matrix_augmented[i*aug_columns+j]*result[j];
        }
         result[i]/=matrix_augmented[i*aug_columns+i];
    }
    free_matrix(matrix_augmented);
    return result;
}
