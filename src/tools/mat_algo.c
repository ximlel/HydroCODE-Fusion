/**
 * @file  mat_algo.c
 * @brief There are some matrix algorithms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/**
 * @brief A function of matrix addition.
 * @details OUTPUT=> C[]: C_{m×n} = A_{m×n} + B_{m×n}.
 */
void mat_add(const double A[], const double B[], double C[], const int m, const int n) 
{
	int i,j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			C[i*n+j] = A[i*n+j] + B[i*n+j];
}


/**
 * @brief A function of matrix subtraction.
 * @details OUTPUT=> C[]: C_{m×n} = A_{m×n} - B_{m×n}.
 */
void mat_sub(const double A[], const double B[], double C[], const int m, const int n)
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			C[i*n+j] = A[i*n+j] - B[i*n+j];
}


/**
 * @brief A function of matrix multiplication.
 * @details OUTPUT=> C[]: C_{m×p} = A_{m×n} * B_{n×p}.
 */
void mat_mul(const double A[], const double B[], double C[], const int m, const int p, const int n)
{
	int i,j,k;
	for (i = 0; i < m; i++)
		for (j = 0; j < p; j++)
			{
				C[i*p+j] = 0.0;
				for (k = 0; k < n; k++)
					C[i*p+j] += A[i*n+k] * B[k*p+j];
			}
}


/**
 * @brief A function to caculate the inverse of the input square matrix.
 * @param[in,out] a: The pointer of the input/output square matrix.
 * @param[in]     n: The order of the input/output square matrix.
 * @return Matrix is invertible or not.
 *    @retval 0: No inverse matrix
 *    @retval 1: Invertible matrix
 * @attention This function is only valid for the matrix with consecutive array.
 */
int rinv(double a[], const int n)
{
    int *is,*js,i,j,k,l,u,v;
    double d,p;
    is=(int*)malloc(n*sizeof(int));
    if(is == NULL)
	{
	    printf("NOT enough memory! INT variable 'is' in 'rinv()' function 'rinv()'.\n");
	    return 0;
	}
    js=(int*)malloc(n*sizeof(int));
    if(js == NULL)
	{
	    free(is);
	    is = NULL;
	    printf("NOT enough memory! INT variable 'js' in 'rinv()' function 'rinv()'.\n");
	    return 0;
	}
    for (k=0; k<=n-1; k++)
		{
			d=0.0;
			for (i=k; i<=n-1; i++)
				for (j=k; j<=n-1; j++)
					{
						l=i*n+j;
						p=fabs(a[l]);
						if (p>d)
							{
								d=p;
								is[k]=i;
								js[k]=j;
							}
					}
			if (d+1.0==1.0)
				{
					free(is);  free(js);
					is = NULL; js = NULL;
					fprintf(stderr, "Error: no inverse matrix!\n");
					return 0;
				}
			if (is[k]!=k)
				for (j=0; j<=n-1; j++)
					{
						u=k*n+j;
						v=is[k]*n+j;
						p=a[u];
						a[u]=a[v];
						a[v]=p;
					}
			if (js[k]!=k)
				for (i=0; i<=n-1; i++)
					{
						u=i*n+k;
						v=i*n+js[k];
						p=a[u];
						a[u]=a[v];
						a[v]=p;
					}
			l=k*n+k;
			a[l]=1.0/a[l];
			for (j=0; j<=n-1; j++)
				if (j!=k)
					{
						u=k*n+j;
						a[u]=a[u]*a[l];
					}
			for (i=0; i<=n-1; i++)
				if (i!=k)
					for (j=0; j<=n-1; j++)
						if (j!=k)
							{
								u=i*n+j;
								a[u]=a[u]-a[i*n+k]*a[k*n+j];
							}
			for (i=0; i<=n-1; i++)
				if (i!=k)
					{
						u=i*n+k;
						a[u]=-a[u]*a[l];
					}
		}
    for (k=n-1; k>=0; k--)
		{
			if (js[k]!=k)
				for (j=0; j<=n-1; j++)
					{
						u=k*n+j;
						v=js[k]*n+j;
						p=a[u];
						a[u]=a[v];
						a[v]=p;
					}
			if (is[k]!=k)
				for (i=0; i<=n-1; i++)
					{
						u=i*n+k;
						v=i*n+is[k];
						p=a[u];
						a[u]=a[v];
						a[v]=p;
					}
		}
    free(is);  free(js);
    is = NULL; js = NULL;
    return 1;
}
