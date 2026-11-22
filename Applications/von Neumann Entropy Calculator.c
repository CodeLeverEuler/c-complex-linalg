/**
 * Quantum Mixed-State Density Matrix & Von Neumann Entropy Calculator
 *
 * Constructs a density matrix ρ = Σ p_i |ψ_i><ψ_i| from an ensemble of n pure quantum states,
 * verifies Hermiticity (ρ = ρ^H), extracts the eigenvalue spectrum via iterative QR decomposition,
 * and computes the Von Neumann entropy in bits: S(ρ) = -Σ λ_i log2(λ_i).
 *
 * Input format: Enter complex wavefunction amplitudes as "<real> <imag>" pairs,
 * followed by the classical ensemble probabilities p_i (summing to 1).
 */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"linear_algebra.h"
int main()
{
    size_t row,n,d,flag=0;
    double sum =0;
    printf("Enter the number of states:\n");
    scanf("%zu",&n);
    if(n<=0)
    {
        printf("\nEnter valid number of states!\n");
        return 1;
    }
    printf("Enter the number of rows:\n");
    scanf("%zu",&row);
    if(row<=0)
    {
        printf("Enter valid row dimensions!");
        return 1;
    }
// Total stacked vector length d = (Hilbert space dimension) * (number of ensemble states)
    d = row*n;
    double* p = (double*)malloc(n*sizeof(double));
    double complex* psi = create_matrix(d,1);
    double complex* arr = create_matrix(d,1);
    if(p==NULL)
    {
        printf("Memory Allocation Failed!");
        return 1;
    }
    if(psi==NULL)
    {
        printf("Memory Allocation Failed!");
        return 1;
    }
    if(arr == NULL)
    {
        printf("Memory Allocation Failed!");
        return 1;
    }
// Read stacked pure state vectors |ψ_1>, |ψ_2>, ..., |ψ_n>
    printf("\nEnter the elements Wavefunction Amplitude Vector (%zu):\n",d*2);
    for(size_t i = 0;i<d;++i)
    {
        double psi_r,psi_i;
        if (scanf("%lf %lf",&psi_r,&psi_i)==2)
        {
            psi[i]= psi_r + psi_i*I; 
        }
    }
//Verify normalization across the n stacked unit state vectors (||psi||_F = sqrt(n))
    double norm = matrix_frobenius_norm(psi,d,1);
    if(cabs(norm-sqrt(n))>1e-6)
    printf("\nEnter valid values!\n"); 
// Read classical mixing probabilities p_i and verify Σ p_i = 1
    printf("\nEnter the elements of the Classical Probability Vector:\n");
    for(size_t i = 0;i<n;++i)
    {
        double p_r;
        if (scanf("%lf",&p_r)==1)
        {
            p[i]= p_r;
            sum +=p[i];
        }
    }
    if(cabs(sum-1)>1e-6)
    {
        printf("Enter valid values!\n");
        return 1;
    }
// Compute <ψ_i| (complex conjugate) and probability-weighted states p_i |ψ_i>
    double complex* phi = complex_conjugate(psi,d,1);
    static size_t j =0;
    for(size_t i=0;i<n;++i)
    {
        for(j; j<d;++j)
        {
            if(j<row*(i+1))
                arr[j]=p[i]*psi[j];
            else 
            break;
        }
    }
 // Construct the (row x row) mixed-state density matrix ρ = Σ p_i |ψ_i><ψ_i|
    double complex* rho = create_matrix(row,row);
    if(rho==NULL)
    {
        printf("Memory Allocation Failed for Density Matrix!\n");
        return 1;
    }
    
    for(size_t i = 0; i < row * row; ++i)
    {
        rho[i] = 0;
    }
    for(size_t i = 0; i < n; ++i)
    {
        for(size_t r = 0; r < row; ++r)
        {
            for(size_t c = 0; c < row; ++c)
            {
                rho[r * row + c] += arr[i * row + r] * phi[i * row + c];
            }
        }
    }
// Verify that the resulting density matrix is Hermitian (ρ == ρ^H)
    double complex* rho_bar = hermitian_conjugate(rho,row,row);
    for(size_t i=0;i<row*row;++i)
    {
        if(cabs(rho_bar[i]-rho[i])>1e-6 )
        {
            printf("Invalid Density Matrix!\n");
            return 1;
        }
    }
    printf("\nThe Density Matrix:\n");
    for(size_t i=0;i<row;++i)
    {
        for (size_t k = 0; k< row;++k)
        {
            if(cimag(rho[i*row+k])==0)
                printf("%6.2lf          ", creal(rho[i*row+k]));
            else if(cimag(rho[i*row+k])>0)
                printf("%6.2lf+%-6.2lfi  ", creal(rho[i*row+k]),cimag(rho[i*row+k]));
            else
                printf("%6.2lf-%-6.2lfi  ", creal(rho[i*row+k]),-cimag(rho[i*row+k]));
        }
        printf("\n");       
    }
// Extract eigenvalues of ρ and compute Von Neumann entropy S(ρ) = -Σ λ_i log2(λ_i)
    double complex* eig_vals = matrix_eigenvalue(rho,row,row,1000,1e-9);
    double entropy=0;
    for(size_t i = 0; i<row; ++i)
    {
        double temp = creal(eig_vals[i]);
        if(temp>1e-12)
        entropy += temp*log2(temp);
    }
    entropy=-entropy;
//Check if it obeys S_max = log(dimensions)
    if(entropy>log2(row))
    {
        printf("ERROR!");
        return 1;
    }
    printf("\nEntropy = %lf\n", entropy);
// Free all allocated arrays and matrices
    free(p);
    free_matrix(psi);
    free_matrix(arr);
    free_matrix(phi);
    free_matrix(rho);
    free_matrix(rho_bar);
    free_matrix(eig_vals);
    return 0;
}