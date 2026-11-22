#include<stdio.h>
#include<math.h>
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
    //Determinant
    printf("Determinant:\n");
    double complex det = matrix_determinant(a,MAX,MAX);
    if(fabs(cimag(det)) < 1e-9)
        printf("%.2lf",creal(det) + 0);
    else if(fabs(creal(det)) < 1e-9)
        printf("%.2lfi",cimag(det) + 0);
    else if(cimag(det) < 0)
        printf("%.2lf%.2lfi",creal(det),cimag(det));
    else 
        printf("%.2lf+%.2lfi",creal(det),cimag(det));
    return 0;
}