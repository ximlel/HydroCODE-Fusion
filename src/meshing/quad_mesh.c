#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../include/var_struc.h"
#include "../include/meshing.h"

#ifndef M_PI
#define M_PI acos(-1.0)
#endif



static int quad_mesh(struct mesh_var * mv, int n_x_add, int n_y_add)
{
	if(isinf(config[13]) || isinf(config[14]))
		{
			fprintf(stderr, "The initial data is not mentioned in a structural mesh!\n");
			exit(2);
		}
	if(isinf(config[10]) || isinf(config[11]) || config[10] < 0.0 || config[11] < 0.0)
		{
			fprintf(stderr, "Without a proper spatial grid size!\n");
			exit(2);
		}
	const int n_x = (int)config[13] + n_x_add, n_y = (int)config[14] + n_y_add;

	const int num_cell = n_x * n_y;
	const int num_border = 2*n_x + 2*n_y;
	if (num_cell > (int)config[3])
		printf("There are %d ghost cell!\n", mv->num_ghost = num_cell - (int)config[3]);

	mv->num_pt = (n_x+1)*(n_y+1);
	mv->X = (double*)malloc(mv->num_pt * sizeof(double));
	mv->Y = (double*)malloc(mv->num_pt * sizeof(double));
	if(mv->X == NULL || mv->Y == NULL)
		{
			printf("Not enough memory in quadrilateral mesh constructed!\n");
			goto return_NULL;
		}
	int k;
	for(k = 0; k < mv->num_pt; k++)
		{
			mv->X[k] = (k%(n_x+1)-n_x/2)*config[10];
			mv->Y[k] = (k/(n_x+1)-n_y/2)*config[11];
		}	
	
	mv->cell_pt = (int**)malloc(num_cell * sizeof(int *));
	if(mv->cell_pt == NULL)
		{
			fprintf(stderr, "Not enough memory in quadrilateral mesh constructed!\n");
			goto return_NULL;
		}
	for(k = 0; k < num_cell; k++)
		{
			mv->cell_pt[k] = (int*)malloc(5 * sizeof(int));
			if(mv->cell_pt[k] == NULL)
				{
					fprintf(stderr, "Not enough memory in CELL_POINT[%d]!\n", k);
					goto return_NULL;
				}			
			mv->cell_pt[k][0] = 4;
			mv->cell_pt[k][1] = k + k/n_x;
			mv->cell_pt[k][2] = mv->cell_pt[k][1] + 1;
			mv->cell_pt[k][3] = k + k/n_x + n_x + 2;
			mv->cell_pt[k][4] = mv->cell_pt[k][3] - 1;
		}

	mv->num_border[0] = 1;	
	mv->num_border[1] = num_border;	
	mv->border_pt = (int*)malloc((num_border+1) * sizeof(int));	
	if(mv->border_pt == NULL)
		{
			printf("Not enough memory in quadrilateral mesh constructed!\n");
			goto return_NULL;
		}	
	for(k = 0; k < n_x+1; k++)	
		mv->border_pt[k] = k; 
	for(k = n_x+1; k < n_x+n_y+1; k++)	
		mv->border_pt[k] = mv->border_pt[k-1] + n_x + 1;
	for(k = n_x+n_y+1; k < n_x*2 + n_y + 1; k++)	
		mv->border_pt[k] = mv->border_pt[k-1] - 1;
	for(k = n_x*2 + n_y + 1; k < num_border; k++)	
		mv->border_pt[k] = mv->border_pt[k-1] - n_x - 1;
	mv->border_pt[num_border] = 0;

	return 1;

 return_NULL:
	free(mv->X);
	mv->X = NULL;
	free(mv->Y);
	mv->Y = NULL;	
	free(mv->border_pt);
	mv->border_pt = NULL;
	for(k = 0; k < num_cell; k++)
		{
			if (mv->cell_pt[k] == NULL)
				break;
			free(mv->cell_pt[k]);
			mv->cell_pt[k] = NULL;
		}
	free(mv->cell_pt);
	mv->cell_pt = NULL;	
	exit(5);	
}

static int quad_border_cond
(struct mesh_var * mv, int n_x_add, int n_y_add,
 int down, int right, int up, int left)
{
	if (down >= 0 || right >= 0 || up >= 0|| left >= 0)
		{
			fprintf(stderr, "Input wrong boundary condition in quadrilateral mesh!\n");
			exit(2);
		}
	else if ((up == -7) - (down == -7) != 0 || (left == -7) - (right == -7) != 0)
		{
			fprintf(stderr, "Periodic boundary condition error!\n");
			exit(2);
		}
	else if ((up == -70) - (down == -70) != 0 || (up == -70 && down == -70 && fabs(config[70]) > config[13]))
		{
			fprintf(stderr, "Periodic boundary condition error!\n");
			exit(2);
		}
	else if ((up == -71) - (down == -71) != 0 || (up == -71 && down == -71 && fabs(config[70]) > config[13]))
		{
			fprintf(stderr, "Periodic boundary condition error!\n");
			exit(2);
		}

	const int n_x = (int)config[13] + n_x_add, n_y = (int)config[14] + n_y_add;
	const int num_cell = n_x * n_y;
	const int num_border = mv->num_border[1];
	int k;

	mv->border_cond = (int*)malloc(num_border * sizeof(int));
	if(mv->border_cond == NULL)
		{
			printf("Not enough memory in quadrilateral boundary constructed!\n");
			goto return_NULL;
		}	
	if(mv->num_ghost > 0)
		{					
			mv->peri_cell = (int*)malloc(num_cell * sizeof(int));
			if(mv->peri_cell == NULL)
				{
					printf("Not enough memory in quadrilateral periodic boundary constructed!\n");
					goto return_NULL;
				}
			for(k = 0; k < num_cell; k++)
				mv->peri_cell[k] = -1;
		}

	for(k = 0; k < n_x; k++)
		{
			if (down == -7)
				mv->peri_cell[k] = k + n_x * (n_y-2);
			else if (down == -70)
				{
					if (k-(int)config[70] < 0)
						mv->peri_cell[k] = n_x*(n_y-2);
					else if (k-(int)config[70] >= n_x)
						mv->peri_cell[k] = n_x*(n_y-1) - 1;
					else
						mv->peri_cell[k] = n_x*(n_y-2)+k - (int)config[70];
				}
			else if (down == -71)
				{
					if (k-(int)config[70] < 0)
						mv->peri_cell[k] = n_x*(n_y-1)+k - (int)config[70];
					else if (k-(int)config[70] >= n_x)
						mv->peri_cell[k] = n_x*(n_y-3)+k - (int)config[70];
					else
						mv->peri_cell[k] = n_x*(n_y-2)+k - (int)config[70];
				}
			else if (down == -77)
				{
					if (k < n_x/2)
						mv->peri_cell[k] = k + n_x + 1;
					else
						mv->peri_cell[k] = k + n_x - 1;
				}
			mv->border_cond[k] = down;
		}
	for(k = n_x; k < n_x+n_y; k++)
		{
			if (right == -7)				
				mv->peri_cell[n_x * (k-n_x+1) - 1] = n_x * (k-n_x) + 1;			
			else if (right == -77)
				{
					if (k < n_x + n_y/2)
						mv->peri_cell[n_x * (k-n_x+1) - 1] = n_x * (k-n_x+2) - 2;
					else
						mv->peri_cell[n_x * (k-n_x+1) - 1] = n_x * (k-n_x) - 2;
				}	
			mv->border_cond[k] = right;
		}
	for(k = n_x + n_y; k < n_x*2 + n_y; k++)
		{
			if (up == -7)
				mv->peri_cell[n_x*(n_y+1) + n_y - k - 1] = n_x*3 + n_y - k - 1;
			else if (up == -70)
				{					
					if (n_x*2+n_y-k-1 + (int)config[70] < 0)
						mv->peri_cell[n_x*(n_y+1)+n_y-k-1] = n_x;
					else if (n_x*2+n_y-k-1 + (int)config[70] >= n_x)
						mv->peri_cell[n_x*(n_y+1)+n_y-k-1] = n_x*2 - 1;
					else
						mv->peri_cell[n_x*(n_y+1)+n_y-k-1] = n_x*3+n_y-k-1 + (int)config[70];
				}
			else if (up == -71)
				{					
					if (n_x*2+n_y-k-1 + (int)config[70] < 0)
						mv->peri_cell[n_x*(n_y+1)+n_y-k-1] = n_x*4+n_y-k-1 + (int)config[70];
					else if (n_x*2+n_y-k-1 + (int)config[70] >= n_x)
						mv->peri_cell[n_x*(n_y+1)+n_y-k-1] = n_x*2+n_y-k-1 + (int)config[70];
					else
						mv->peri_cell[n_x*(n_y+1)+n_y-k-1] = n_x*3+n_y-k-1 + (int)config[70];
				}
			else if (up == -77)
				{
					if (k < n_x + n_y + n_x/2)
						mv->peri_cell[n_x*(n_y+1) + n_y - k - 1] = n_x * n_y + n_y-k-2;
					else
						mv->peri_cell[n_x*(n_y+1) + n_y - k - 1] = n_x * n_y + n_y-k;
				}
			mv->border_cond[k] = up;
		}
	for(k = n_x*2 + n_y; k < num_border; k++)
		{
			if (left == -7)
				mv->peri_cell[(num_border-k-1) * n_x] = (num_border-k) * n_x - 2;				
			else if (left == -77)
				{
					if (k < n_x*2 + n_y + n_y/2)
						mv->peri_cell[(num_border-k-1) * n_x] = (num_border-k-2) * n_x + 1;
					else
						mv->peri_cell[(num_border-k-1) * n_x] = (num_border-k) * n_x + 1;
				}
			mv->border_cond[k] = left;
		}

	if(mv->num_ghost > 0)
		period_cell_modi(mv);

	return 1;

 return_NULL:
	free(mv->border_cond);
	mv->border_cond = NULL;	
	free(mv->peri_cell);
	mv->peri_cell = NULL;
	exit(5);	
}

void Sod_mesh(struct mesh_var * mv)
{
	const int n_x_a = 0, n_y_a = 0;
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -2, -4, -2, -4);
}

void Shock_Bubble_mesh(struct mesh_var * mv)
{
	const int n_x_a = 0, n_y_a = 0;
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -4, -4, -4, -4);
}

void Shear_mesh(struct mesh_var * mv)
{	
	const int n_x_a = 0, n_y_a = 0;
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -1, -4, -4, -4);
}

void free_mesh(struct mesh_var * mv)
{		
	const int n_x_a = 0, n_y_a = 0;	
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -3, -3, -3, -3);
}

void RMI_mesh(struct mesh_var * mv)
{	
	const int n_x_a = 2, n_y_a = 0;		
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -4, -7, -4, -7);
}

void RMI_S_mesh(struct mesh_var * mv)
{	
	const int n_x_a = 0, n_y_a = 0;		
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -4, -4, -4, -4);
}

void R2D_mesh(struct mesh_var * mv)
{	
	const int n_x_a = 2, n_y_a = 2;		
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -77, -77, -77, -77);
}

void Vortex_mesh(struct mesh_var * mv)
{	
	const int n_x_a = 2, n_y_a = 2;		
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -7, -7, -7, -7);
}

void Shell_mesh(struct mesh_var * mv)
{	
	const int n_x_a = 0, n_y_a = 0;		
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -2, -4, -4, -2);
}

void cylinder_mesh(struct mesh_var * mv)
{
	const int n_x_a = 2, n_y_a = 0;
	quad_mesh(mv, n_x_a, n_y_a);
	
	const int n_x = (int)config[13] + n_x_a, n_y = (int)config[14] + n_y_a;
	const double R = config[11]*n_y/M_PI*9.0/8.0;	
	for(int k = 0; k < (n_x+1)*(n_y+1); k++)
		{
			mv->X[k] = (R+(n_x-k%(n_x+1))*config[10])*cos(13.0/9.0*M_PI-(k/(n_x+1))*8.0/9.0*M_PI/n_y);
			mv->Y[k] = (R+(n_x-k%(n_x+1))*config[10])*sin(13.0/9.0*M_PI-(k/(n_x+1))*8.0/9.0*M_PI/n_y);
		}	
	quad_border_cond(mv, n_x_a, n_y_a, -3, -2, -3, -1);
}

void odd_even_mesh(struct mesh_var * mv)
{
	const int n_x_a = 0, n_y_a = 0;
	quad_mesh(mv, n_x_a, n_y_a);
	
	const int n_x = (int)config[13] + n_x_a, n_y = (int)config[14] + n_y_a;
	for(int k = 0; k < (n_y+1); k++)
		mv->Y[(n_x/2)*(n_y+1)+k] += ((k%2)-0.5)*0.002*config[11];

	quad_border_cond(mv, n_x_a, n_y_a, -2, -3, -2, -1);
}

void odd_even_periodic_mesh(struct mesh_var * mv)
{
	const int n_x_a = 2, n_y_a = 0;
	quad_mesh(mv, n_x_a, n_y_a);
	
	const int n_x = (int)config[13] + n_x_a, n_y = (int)config[14] + n_y_a;
	for(int k = 0; k < (n_y+1); k++)
		mv->Y[(n_x/2)*(n_y+1)+k] += ((k%2)-0.5)*0.002*config[11];

	quad_border_cond(mv, n_x_a, n_y_a, -2, -7, -2, -7);
}

void odd_even_inflow_mesh(struct mesh_var * mv)
{
	const int n_x_a = 0, n_y_a = 2;
	quad_mesh(mv, n_x_a, n_y_a);
	
	const int n_x = (int)config[13] + n_x_a, n_y = (int)config[14] + n_y_a;
	for(int k = 0; k < (n_y+1); k++)
		mv->Y[(n_x/2)*(n_y+1)+k] += ((k%2)-0.5)*0.002*config[11];

	quad_border_cond(mv, n_x_a, n_y_a, -7, -3, -7, -1);
}

void rand_disturb_inflow_mesh(struct mesh_var * mv)
{
	const int n_x_a = 0, n_y_a = 2;

	quad_mesh(mv, n_x_a, n_y_a);

	const int n_x = (int)config[13] + n_x_a, n_y = (int)config[14] + n_y_a;
	srand((unsigned) time(NULL)); //seed--time.
	for(int k = n_x+1; k < n_y*(n_x+1); k++)
		{
			if(k%(n_x+1)&&(k%(n_x+1))!=n_x+1)
				{							
					mv->X[k] += (0.5-(rand()%10001)/10000.0)*0.000001*config[11];
					mv->Y[k] += (0.5-(rand()%10001)/10000.0)*0.000001*config[11];
				}
		}

	quad_border_cond(mv, n_x_a, n_y_a, -7, -3, -7, -1);
}

void oblique_periodic_mesh(struct mesh_var * mv)
{
	const int n_x_a = 0, n_y_a = 2;
	quad_mesh(mv, n_x_a, n_y_a);
	quad_border_cond(mv, n_x_a, n_y_a, -70, -3, -70, -3);
}

static int quad_border_normal_velocity
(struct mesh_var * mv, int n_x_add, int n_y_add,
 double down, double right, double up, double left)
{
	const int n_x = (int)config[13] + n_x_add, n_y = (int)config[14] + n_y_add;
	const int num_border = mv->num_border[1];
	int k;
	
	mv->normal_v = (double*)malloc(num_border * sizeof(double));
	if(mv->normal_v == NULL)
		{
			printf("Not enough memory in quadrilateral boundary constructed!\n");
			goto return_NULL;
		}

	for(k = 0; k < n_x; k++)
		{
			mv->normal_v[k] = down;
		}
	for(k = n_x; k < n_x+n_y; k++)
		{
			mv->normal_v[k] = right;
		}
	for(k = n_x + n_y; k < n_x*2 + n_y; k++)
		{
			mv->normal_v[k] = up;
		}
	for(k = n_x*2 + n_y; k < num_border; k++)
		{			
			mv->normal_v[k] = left;
		}

	return 1;

 return_NULL:
	free(mv->normal_v);
	mv->normal_v = NULL;
	exit(5);	
}

void Saltzman_mesh_Lag(struct mesh_var * mv)
{
	const int n_x_a = 0, n_y_a = 0;
	quad_mesh(mv, n_x_a, n_y_a);

	for(int k = 0; k < mv->num_pt; k++)		
		mv->X[k] += (0.1 - mv->Y[k])*sin(M_PI * mv->X[k]);
	
	quad_border_cond(mv, n_x_a, n_y_a, -2, -2, -2, -2);
	quad_border_normal_velocity(mv, n_x_a, n_y_a, 0.0, 0.0, 0.0, -1.0);
}
