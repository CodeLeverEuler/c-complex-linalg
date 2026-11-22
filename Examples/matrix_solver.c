#include <stdio.h>
#include "linear_algebra.h"
#define MAX 2
int main()
{
    //Matrix Creation
    printf("Enter elements of a:\n");
    double complex* a = create_matrix(MAX,MAX);
    double complex* b = create_matrix(MAX,1);
    if(a==NULL)
        {
            printf("Memory allocation failed for matrix a!\n");
            return 1;
        }
    if(b==NULL)
        {
            printf("Memory allocation failed for matrix b!\n");
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
                free_matrix(b);
                return 1;
            }
            a[i*MAX+j] = a_r + a_i*I;
        }   
    }
    printf("Enter elements for b:\n");
    for(size_t i=0;i<MAX;++i)
    {
        double b_r, b_i;
        if(scanf("%lf %lf", &b_r, &b_i)!=2)
        {
            printf("Wrong input!");
            free_matrix(a);
            free_matrix(b);
            return 1;
        }
        b[i] = b_r + b_i*I;
    }
    //solve
    double complex* solve = matrix_solver(a,b,MAX,MAX);
    if (solve==NULL)
    {
        printf("System has no unique solution!\n");
        free_matrix(a);
        free_matrix(b);
        return 1;
    }
    printf("Solution:\n");
    for(size_t i=0;i<MAX;++i)
    {
        if(fabs(cimag(solve[i])) < 1e-9)
            printf("%.2lf",creal(solve[i]) + 0);
        else if(fabs(creal(solve[i])) < 1e-9)
            printf("%.2lfi",cimag(solve[i]) + 0);
        else if(cimag(solve[i]) < 0)
            printf("%.2lf%.2lfi",creal(solve[i]),cimag(solve[i]));
        else 
            printf("%.2lf+%.2lfi",creal(solve[i]),cimag(solve[i]));
    }
    free_matrix(a);
    free_matrix(b);
    free_matrix(solve);
    return 0;
}