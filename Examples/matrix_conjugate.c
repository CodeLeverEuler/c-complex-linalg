#include<stdio.h>
#include"linear_algebra.h"
#define MAX 2
int main()
{
    //Matrix Creation
    printf("Enter matrix\n");
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
    //Transposition
    double complex* a_bar = matrix_conjugate(a,MAX,MAX);
    printf("Complex Conjugate:\n");
    for(size_t i=0;i<MAX;++i)
    {
        for (size_t j = 0; j< MAX;++j)
        {
            if(fabs(cimag(a_bar[i*MAX+j])) < 1e-9)
                printf("%.2lf",creal(a_bar[i*MAX+j]) + 0);
            else if(fabs(creal(a_bar[i*MAX+j])) < 1e-9)
                printf("%.2lfi",cimag(a_bar[i*MAX+j]) + 0);
            else if(cimag(a_bar[i*MAX+j]) < 0)
                printf("%.2lf%.2lfi",creal(a_bar[i*MAX+j]),cimag(a_bar[i*MAX+j]));
            else 
            printf("%.2lf+%.2lfi",creal(a_bar[i*MAX+j]),cimag(a_bar[i*MAX+j]));
        }
        
    }
    free_matrix(a);
    free_matrix(a_bar);
    return 0;
}