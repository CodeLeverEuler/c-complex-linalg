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
    //Transposition
    double complex* trp = complex_conjugate(a,MAX,MAX);
    printf("Transpose:\n");
    for(size_t i=0;i<MAX;++i)
    {
        for (size_t j = 0; j< MAX;++j)
        {
            if(fabs(cimag(trp[i*MAX+j])) < 1e-9)
                printf("%.2lf",creal(trp[i*MAX+j]) + 0);
            else if(fabs(creal(trp[i*MAX+j])) < 1e-9)
                printf("%.2lfi",cimag(trp[i*MAX+j]) + 0);
            else if(cimag(trp[i*MAX+j]) < 0)
                printf("%.2lf%.2lfi",creal(trp[i*MAX+j]),cimag(trp[i*MAX+j]));
            else 
            printf("%.2lf+%.2lfi",creal(trp[i*MAX+j]),cimag(trp[i*MAX+j]));
        }
        
    }
    free_matrix(a);
    free_matrix(trp);
    return 0;
}