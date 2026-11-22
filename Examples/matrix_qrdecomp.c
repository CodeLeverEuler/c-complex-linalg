#include<stdio.h>
#include"linear_algebra.h"
#define MAX 2
int main()
{
    //Matrix Creation
    double complex* a = create_matrix(MAX,MAX);
    if(a==NULL)
    {
        printf("Memory allocation failed for matrix a\n");
        return 1;
    }
    printf("Enter elements for a:\n");
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
    //R Decomposition
    double complex* R = matrix_rdecomposition(a,MAX,MAX);
    printf("R:\n");
    for(size_t i=0;i<MAX;++i)
    {
        for (size_t j = 0; j< MAX;++j)
        {
            if(fabs(cimag(R[i*MAX+j])) < 1e-9)
                printf("%.2lf",creal(R[i*MAX+j]) + 0);
            else if(fabs(creal(R[i*MAX+j])) < 1e-9)
                printf("%.2lfi",cimag(R[i*MAX+j]) + 0);
            else if(cimag(R[i*MAX+j]) < 0)
                printf("%.2lf%.2lfi",creal(R[i*MAX+j]),cimag(R[i*MAX+j]));
            else 
            printf("%.2lf+%.2lfi",creal(R[i*MAX+j]),cimag(R[i*MAX+j]));
        }
    }
    //Q Decomposition
    double complex* Q = matrix_qdecomposition(a,MAX,MAX);
    printf("Q:\n");
    for(size_t i=0;i<MAX;++i)
    {
        for (size_t j = 0; j< MAX;++j)
        {
            if(fabs(cimag(Q[i*MAX+j])) < 1e-9)
                printf("%.2lf",creal(Q[i*MAX+j]) + 0);
            else if(fabs(creal(Q[i*MAX+j])) < 1e-9)
                printf("%.2lfi",cimag(Q[i*MAX+j]) + 0);
            else if(cimag(Q[i*MAX+j]) < 0)
                printf("%.2lf%.2lfi",creal(Q[i*MAX+j]),cimag(Q[i*MAX+j]));
            else 
            printf("%.2lf+%.2lfi",creal(Q[i*MAX+j]),cimag(Q[i*MAX+j]));
        }
        
    }

    free_matrix(a);
    free_matrix(R);
    free_matrix(Q);
    return 0;
}