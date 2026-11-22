/**
 * AC Circuit Mesh & Complex Power Analyzer
 *
 * Solves multi-loop AC networks using Kirchhoff's Voltage Law in matrix form (Z * I = V)
 * and computes total complex power (S = V · I*), separating Active Power (Watts)
 * and Reactive Power (VAR).
 *
 * Input format: Enter complex values as "<real> <imag>" pairs.
 */
#include<stdio.h>
#include"linear_algebra.h"
int main()
{
    size_t size_v;
    printf("Enter the number of voltage sources:\n");
    scanf("%zu",&size_v);
    // Allocate and populate the (N x 1) complex voltage vector V
    double complex* V =  create_matrix(size_v,1);
    if(V==NULL) return 1;
    printf("Enter the values of the voltage sources:\n");
    for(size_t i = 0; i<size_v; ++i)
    {
        double complex V_r,V_i;
        if(scanf("%lf %lf", &V_r, &V_i)!=2)
        {
            printf("Wrong input format!");
            free_matrix(V);
            return 1;
        }
        V[i] = V_r + V_i*I;
    }
    // Allocate and populate the (N x N) complex mesh impedance matrix Z
    double complex* Z =  create_matrix(size_v,size_v);
    if(Z==NULL) return 1;
    printf("Enter the values of the impedances:\n");
    for(size_t i = 0; i<size_v; ++i)
    {
       for(size_t j =0; j<size_v;++j)
       {
         double complex Z_r,Z_i;
        if(scanf("%lf %lf", &Z_r, &Z_i)!=2)
        {
            printf("Wrong input format!");
            free_matrix(Z);
            return 1;
        }
        Z[i*size_v+j] = Z_r + Z_i*I;
       }
    }
// Solve Z * C = V for mesh currents C, then compute total complex power S = <C, V> = sum(conj(C_i) * V_i)
    double complex* C = matrix_solver(Z,V,size_v, size_v);
    double complex S = inner_product(C,V,size_v,1);
// Display per-loop mesh currents in rectangular form
    printf("\nMesh Currents:\n");
    for(size_t i=0;i<size_v;++i)
        {
            if(cimag(C[i]) ==0)
                printf("%lf\n",creal(C[i]));
            else if(cimag(C[i])<0)
                printf("%lf %lf\n", creal(C[i]),cimag(C[i]));
            else
            printf("%lf+%lf\n", creal(C[i]),cimag(C[i]));
        }
// Display per-loop mesh currents in rectangular form
    printf("Power:\n");
    if(cimag(S)==0)
        printf("%.2lf\t", creal(S));
    else if(cimag(S)<0)
        printf("%.2lf%.2lfi\t", creal(S),cimag(S));
    else
        printf("%.2lf+ %.2lfi\t", creal(S),cimag(S));
    printf("\nActive Power:\n");
    printf("%lf", creal(S));
    printf("\nReactive Power:\n");
    printf("%lf", cimag(S));
// Clean up allocated matrices
    free_matrix(V);
    free_matrix(C);
    return 0;
}