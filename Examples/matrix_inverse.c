#include <stdio.h>
#include "linear_algebra.h"
#define MAX 2
int main()
{
    //Matrix Creation
    printf("Enter elements of a:\n");
    double complex* a = create_matrix(MAX,MAX);
    if(a==NULL)
        {
            printf("Memory allocation failed for matrix a!\n");
            return 1;
        }
    for(size_t i=0;i<MAX;++i)
    {
        for (size_t j = 0; j< MAX;++j)
        {
            double a_r, a_i;
            if(scanf("%lf %lf", &a_r, &a_i)!=2)
            {
                printf("Wrong input!");
                free_matrix(a);
                return 1;
            }
            a[i*MAX+j] = a_r + a_i*I;
        }   
    }
    //Inversion
    double complex* inv = matrix_inverse(a,2,2);
    if(inv==NULL)
    {
        printf("Matrix is singular!");
        return 1;
    }
    printf("Inverse:\n");
    for(size_t i=0;i<MAX;++i)
    {
        for (size_t j = 0; j< MAX;++j)
        {
            if(fabs(cimag(inv[i*MAX+j])) < 1e-9)
                printf("%.2lf",creal(inv[i*MAX+j]) + 0);
            else if(fabs(creal(inv[i*MAX+j])) < 1e-9)
                printf("%.2lfi",cimag(inv[i*MAX+j]) + 0);
            else if(cimag(inv[i*MAX+j]) < 0)
                printf("%.2lf%.2lfi",creal(inv[i*MAX+j]),cimag(inv[i*MAX+j]));
            else 
            printf("%.2lf+%.2lfi",creal(inv[i*MAX+j]),cimag(inv[i*MAX+j]));
        }
        
    }
    free_matrix(a);
    free_matrix(inv);
    return 0;
    
}