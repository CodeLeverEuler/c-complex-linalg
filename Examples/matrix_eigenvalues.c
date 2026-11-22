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
    //Eigenvalue Solving
    double complex* eigval = matrix_eigenvalue(a,MAX,MAX,1000,1e-9);
    printf("Eigenvalues:\n");
    for(size_t i=0;i<MAX;++i)
    {
        if(fabs(cimag(eigval[i])) < 1e-9)
            printf("%.2lf",creal(eigval[i]) + 0);
        else if(fabs(creal(eigval[i])) < 1e-9)
            printf("%.2lfi",cimag(eigval[i]) + 0);
        else if(cimag(eigval[i]) < 0)
            printf("%.2lf%.2lfi",creal(eigval[i]),cimag(eigval[i]));
        else 
            printf("%.2lf+%.2lfi",creal(eigval[i]),cimag(eigval[i]));  
    }
    free_matrix(a);
    free_matrix(eigval);
    return 0;
}