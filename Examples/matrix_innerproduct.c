#include <stdio.h>
#include "linear_algebra.h"
#define MAX 2
int main()
{
    //Matrix Memory Allocation
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
    //Multiplication
    double complex ip = inner_product(a,b,MAX,MAX);
    printf("Inner Product:\n");
        if(fabs(cimag(ip)) < 1e-9)
            printf("%.2lf",creal(ip) + 0);
        else if(fabs(creal(ip)) < 1e-9)
            printf("%.2lfi",cimag(ip) + 0);
        else if(cimag(ip) < 0)
            printf("%.2lf%.2lfi",creal(ip),cimag(ip));
        else 
            printf("%.2lf+%.2lfi",creal(ip),cimag(ip));     
    free_matrix(a);
    free_matrix(b);
    return 0; 
}