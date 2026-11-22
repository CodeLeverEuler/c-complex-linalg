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
            double temp_r, temp_i;
            if(scanf("%lf %lf", &temp_r, &temp_i)!=2)
            {
                printf("Wrong input!");
                free_matrix(a);
                return 1;
            }
            a[i*MAX+j] = temp_r + temp_i*I;
        }   
    }
    //Exponentiation
    double complex* exp = matrix_pow(a,2,MAX,MAX);
    printf("Result:\n");
    for(size_t i=0;i<MAX;++i)
    {
        for (size_t j = 0; j< MAX;++j)
        {
            if(fabs(cimag(exp[i*MAX+j])) < 1e-9)
                printf("%.2lf",creal(exp[i*MAX+j]) + 0);
            else if(fabs(creal(exp[i*MAX+j])) < 1e-9)
                printf("%.2lfi",cimag(exp[i*MAX+j]) + 0);
            else if(cimag(exp[i*MAX+j]) < 0)
                printf("%.2lf%.2lfi",creal(exp[i*MAX+j]),cimag(exp[i*MAX+j]));
            else 
            printf("%.2lf+%.2lfi",creal(exp[i*MAX+j]),cimag(exp[i*MAX+j]));
        }
        
    }
    free_matrix(a);
    free_matrix(exp);
    return 0;
}