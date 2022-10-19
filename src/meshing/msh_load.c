#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "../include_cii/mem.h"
#include "../include/var_struc.h"

//! Define default boundary condition.
#define DEFAULT_BC -3


static int msh_border_cell_dis(int type)
{
    if (type == 1)
	return 1;
    else if (type == 2 || type == 3)
	return 0;
    else
	return -1;
}


static int msh_border_build(struct mesh_var * mv, int num_bc_all)
{
	if (mv->num_border[0] != 1)
		{
			fprintf(stderr,"Not only one boundary!\n");
			return 0;
		}

	const int num_border = mv->num_border[1];
	if (num_border < 2)
		{
			fprintf(stderr,"Number of border point is less than 2!\n");
			return 0;
		}

	mv->border_pt   = (int*)ALLOC((num_border+1)* sizeof(int));
	mv->border_cond = (int*)ALLOC(num_border* sizeof(int));

	int  * bp = mv->border_pt;
	int ** cp = mv->cell_pt;
	int  * bc = mv->border_cond;
	int  * ct = mv->cell_type;

	bp[0] = cp[num_bc_all-1][1];
	bc[0] = ct[num_bc_all-1];

	bp[1] = cp[num_bc_all-1][2];

	// not num_border-1, the lastest one is used to test loop.
	for (int j = 1; j < num_border; j++)
	    for (int i = num_bc_all-2; i >= num_bc_all-num_border; i--)
		{
		    if (cp[i][0] != 2)
			{
			    fprintf(stderr,"There are not 2 border point on a boundary cell in 2-D case!\n");
			    return 0;
			}
		    if(bp[j] == cp[i][1] && bp[j-1] != cp[i][2])
			{
			    bp[j+1] = cp[i][2];
			    bc[j] = ct[i];
			}
		    else if(bp[j] == cp[i][2] && bp[j-1] != cp[i][1])
			{
			    bp[j+1] = cp[i][1];
			    bc[j] = ct[i];
			}
		}
	if (bp[num_border] != bp[0])
	    {
		fprintf(stderr,"The boundary isn't a loop!\n");
		return 0;
	    }

	return 1;	
}


// only one boundary loop is supported.
int msh_read(FILE * fp, struct mesh_var * mv)
{
	int s_now = -1, s_max = 0; //section index
	int num_of = 0, num_tag; // number of section or tags.
	int *idx_N = NULL; // order of NODE;
	int num_cell = (int)config[3], num_border = 0;
	int b_or_c, n_bc = 0, num_bc_all = 0; // boundary_or_cell
	int temp[30]; // store cell_point data
	int n_c, type, phy_entity;

	const char *section[] = {
		"MeshFormat", "PhysicalNames", "Nodes",
		"Elements" , "Periodic", "NodeData", "ElementData",
		"ElementNodeData", "InterpolationScheme" };
	char one_line[1000]; // input one line
	char *headptr, *endptr; // PTR of read data
	char Endsection[25];

	while (fgets(one_line, sizeof(one_line), fp) != NULL)
		{
			for (headptr = one_line; isspace(*headptr); )
				headptr++;

			if (*headptr == '$')
				{
					for (endptr = (++headptr); !isspace(*endptr); )
						endptr++;
					*endptr = '\0';

					if (s_now >= 0)
						{
							strcpy(Endsection, "End");
							strcat(Endsection, section[s_now]);
							if (strcmp(headptr, Endsection) != 0)
								{
									fprintf(stderr, "End of the section in .msh file doesn't match!\n");
									goto return_0;
								}
							else if (num_of > 0)
								{
									fprintf(stderr, "The count is not complete in the section in .msh file!\n");
									goto return_0;
								}
							else
								s_now = -1;

							continue;
						}

					for (int s = s_max; s < 9; s++)
						if (strcmp(headptr, section[s]) == 0)
							{
								s_now = s;
								s_max = s;
								break;
							}
					if (s_now == 0 && num_of <= 0)
						num_of = 1;
				}
			else if (strlen(headptr) <= 0)
				continue;
			else if (s_now == 0 && num_of == 1)
				{
					if (strtod(headptr, &headptr) - 2.2 > EPS)
						{
							fprintf(stderr, "Version-number isn't not equal to 2.2 in .msh file!\n");
							goto return_0;
						}
					if (strtol(headptr, &headptr, 10) != 0)
						{
							fprintf(stderr, "The .msh file isn't ASCII file format!\n");
							goto return_0;
						}
					if (strtol(headptr, NULL, 10) != 8)
						{
							fprintf(stderr, "Currently only data-size = sizeof(double) is supported in .msh file!\n");
							goto return_0;
						}
					num_of--;
				}
			else if (s_now == 0)
				{
					fprintf(stderr, "Not only one line in Meshformat!\n");
					goto return_0;
				}
			else if (s_now == 2 && num_of <= 0)
				{
					num_of = strtol(headptr, NULL, 10);
					if (num_of > 0)
						{
							mv->num_pt = num_of;
							idx_N = (int*)   ALLOC(mv->num_pt * sizeof(int));
							mv->X = (double*)ALLOC(mv->num_pt * sizeof(double));
							mv->Y = (double*)ALLOC(mv->num_pt * sizeof(double));
							if(mv->X == NULL || mv->Y == NULL)
								{
									printf("Not enough memory in msh_read!\n");
									goto return_0;
								}
						}
				}
			else if (s_now == 2)
				{
					idx_N[--num_of] = strtol(headptr, &headptr, 10);
					mv->X[num_of]   = strtod(headptr, &headptr);
					mv->Y[num_of]   = strtod(headptr, &headptr);
				}
			else if (s_now == 3 && num_of <= 0 && idx_N != NULL)
				{
					num_of = strtol(headptr, NULL, 10);
					if (num_of <= 0)
						continue;

					num_cell = 0;
					num_border = 0;
					num_bc_all = num_of;
					mv->cell_type = (int*) ALLOC(num_bc_all * sizeof(int));
					mv->cell_pt   = (int**)ALLOC(num_bc_all * sizeof(int *));
					if(mv->cell_type == NULL || mv->cell_pt == NULL)
						{
							fprintf(stderr, "Not enough memory in msh_read!\n");
							goto return_0;
						}
				}
			else if (s_now == 3 && idx_N != NULL)
				{
					--num_of;
					strtol(headptr, &headptr, 10);

					type = strtol(headptr, &headptr, 10);

					b_or_c = msh_border_cell_dis(type);
					if (b_or_c == -1)
						continue;
					if (b_or_c == 0)
						{
							n_bc = num_cell++;
							mv->cell_type[n_bc] = type;
						}

					num_tag = strtol(headptr, &headptr, 10);
					if (num_tag < 2)
						{
							fprintf(stderr, "Using the MSH2 format require at least the first two tags!\n");
							continue;
						}
					phy_entity = strtol(headptr, &headptr, 10);
					if (b_or_c == 1)
						{
							n_bc = num_bc_all - (++num_border);
							mv->cell_type[n_bc] = phy_entity ? -phy_entity : DEFAULT_BC;
						}

					while(--num_tag > 0)
						strtol(headptr, &headptr, 10);

					for(n_c = 0; (temp[n_c] = strtol(headptr, &headptr, 10)) > 0; )
						n_c++;


					mv->cell_pt[n_bc] = (int*)ALLOC((n_c+1) * sizeof(int));
					if(mv->cell_pt[n_bc] == NULL)
						{
							fprintf(stderr, "Not enough memory in msh_read!\n");
							goto return_0;
						}

					mv->cell_pt[n_bc][0] = n_c;
					for(int j = 0; j < n_c; j++)
						for(int i = 0; i < mv->num_pt; i++)
							if(temp[j] == idx_N[i])
								mv->cell_pt[n_bc][j+1] = i;
				}

		}

	if (ferror(fp))
		{
			fprintf(stderr, "Read error occurrs in .msh file!\n");
			goto return_0;
		}

	if (num_cell > (int)config[3])
		printf("There are %d ghost cell!\n", mv->num_ghost = num_cell - (int)config[3]);
	else if (num_cell < (int)config[3])
		{
			fprintf(stderr,"There are not enough cell in .msh file!\n");
			goto return_0;
		}

	mv->num_border[0] = 1;
	mv->num_border[1] = num_border;

	if(msh_border_build(mv, num_bc_all) == 0)
		{
			FREE(mv->border_pt);
			FREE(mv->border_cond);
			goto return_0;
		}

	FREE(idx_N);
	return 1;

 return_0:

	FREE(idx_N);
	FREE(mv->X);
	FREE(mv->Y);
	FREE(mv->cell_type);
	if (mv->cell_pt != NULL)
		{
			for(int i = 0; i < num_bc_all; ++i)
				FREE(mv->cell_pt[i]);
			FREE(mv->cell_pt);
		}
	return 0;
}
