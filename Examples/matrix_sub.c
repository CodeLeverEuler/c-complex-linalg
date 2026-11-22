#include <stdio.h>
#include "linear_algebra.h"
#define MAX 2
int main()
{
    //Matrix Creation
    printf("Enter elements of a:\n");
    double complex* a = create_matrix(MAX,MAX);
    double complex* b = create_matrix(MAX,MAX);
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
                return 1;
            }
            a[i*MAX+j] = a_r + a_i*I;
        }   
    }
    printf("Enter elements for b:\n");
    for(size_t i=0;i<MAX;++i)
    {
        for (size_t j = 0; j< MAX;++j)
        {
            double b_r, b_i;
            if(scanf("%lf %lf", &b_r, &b_i)!=2)
            {
                printf("Wrong input!");
                free_matrix(b);
                return 1;
            }
            b[i*MAX+j] = b_r + b_i*I;
        }   
    }
    //Subtraction
    double complex* sub = matrix_subtraction(a,b,2,2);
    printf("Difference:\n");
    for(size_t i=0;i<MAX;++i)
    {
        for (size_t j = 0; j< MAX;++j)
        {
            if(fabs(cimag(sub[i*MAX+j])) < 1e-9)
                printf("%.2lf",creal(sub[i*MAX+j]) + 0);
            else if(fabs(creal(sub[i*MAX+j])) < 1e-9)
                printf("%.2lfi",cimag(sub[i*MAX+j]) + 0);
            else if(cimag(sub[i*MAX+j]) < 0)
                printf("%.2lf%.2lfi",creal(sub[i*MAX+j]),cimag(sub[i*MAX+j]));
            else 
            printf("%.2lf+%.2lfi",creal(sub[i*MAX+j]),cimag(sub[i*MAX+j]));
        }
        
    }
    free_matrix(a);
    free_matrix(b);
    free_matrix(sub);
    return 0;
}