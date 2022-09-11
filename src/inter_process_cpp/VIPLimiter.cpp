///////////////////////////////////////////////////
/// @file VIPLimiter.cpp
/// @brief  The subroutin implements the VIP limiter for simulations of 2D flows on structured grids.
/// @note   Note, this is only a limiter to limit the velocity vector V=(u,v) for 2D flows. 
/// @sa     The specific implementation is mainly based on section 2.1-2.3 of the following paper: \n
///         G.Luttwak, J.Falcovitz, Slope limiting for vectors, A novel vector limiting algorithm,
///         Int. J. Numer. Meth. Fluids., 2011 (65)1365-1375.
/// @date   Apr 11, 2018
/// @author Jian Cheng @ IAPCM
///////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

#include "../include/var_struc.h"


///////////////////////////////////////////////////
//some subroutines called by useViPLimiter...
///////////////////////////////////////////////////
static double getTriArea(double x0, double y0, double x1, double y1, double xp, double yp);
static void getPerpendFoot(double x0, double y0, double x1, double y1, double xc, double yc, double* pf);
static bool obtuseAngle(double x0, double y0, double xa, double ya, double xb, double yb);
static bool insideSegment(double x0, double y0, double x1, double y1, double xp, double yp);
static double insectionPoint(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double* Vp);
static bool insideTriCH(std::vector<std::vector<double> >& CH, bool flag, double* Vp);
static bool insideQuadCH(std::vector<std::vector<double> >& CH, bool flag, double* Vp);
static bool insideTriCH(std::vector<std::vector<double> >& CH, bool flag, double* V0, double* Vp, double& lambda);
static bool insideQuadCH(std::vector<std::vector<double> >& CH, bool flag, double* V0, double* Vp, double& lambda);


///////////////////////////////////////////////////
/// @brief Subroutine of using VIP limiter for 2D velocity vector
/// @param[in] neigh_cell_num: number of neighbor cells
///                            (Note: only 2D face neighbor cells are used here, thus, 0<=neigh_cell_num<=4)
/// @param[in] Vave: matrix for the average velocity vectors of neighbor cells
///                 i.e. [u^0_ave,v^0_ave]
///                      [u^1_ave,v^1_ave]
///                       [u^2_ave,v^2_ave]
///                       ...
/// @param[in] V0: vector for the cell average velocity of the target cell
/// @param[in,out] Vp: vector for the velocity of a given location in the target cell where VIP limiter needs to be applied
/// @return the limiting coefficient lambda ( in [0, 1] ) for gradient vector
///////////////////////////////////////////////////
double useVIPLimiter(int neigh_cell_num, double Vave[][2], double* V0, double* Vp)
{
	double const Alpha = config[41];
	Vp[0] = V0[0]+(Vp[0]-V0[0])*2.0/Alpha;
	Vp[1] = V0[1]+(Vp[0]-V0[0])*2.0/Alpha;

	bool colinear(true);
	int count(0);
	double area(0);
	std::vector<std::vector<double> > CH(2, std::vector<double>(2));

	//Temporally, we choose to do noting in these cases...
	if(neigh_cell_num != 3 && neigh_cell_num != 4)
		return 1.0;

	//Set initial guess for the CH...
	if(Vave[0][0] > Vave[1][0])
	{
		CH[0][0] = Vave[0][0];
		CH[0][1] = Vave[0][1];
		CH[1][0] = Vave[1][0];
		CH[1][1] = Vave[1][1];
	}
	else
	{
		CH[0][0] = Vave[1][0];
		CH[0][1] = Vave[1][1];
		CH[1][0] = Vave[0][0];
		CH[1][1] = Vave[0][1];
	}

	//Now, we try to find the initial triangle for the CH
	for(int e = 2; e < neigh_cell_num; ++e)
	{
		area = getTriArea(CH[0][0], CH[0][1], CH[1][0], CH[1][1], Vave[e][0], Vave[e][1]);

		if( fabs(area) < EPS ) //node[e] is colinear with node[0] and node[1]
		{
			if( Vave[e][0] > CH[0][0] )
			{
				CH[0][0] = Vave[e][0];
				CH[0][1] = Vave[e][1];
			}
			else if( Vave[e][0] < CH[1][0] )
			{
				CH[1][0] = Vave[e][0];
				CH[1][1] = Vave[e][1];
			}
			else
			{
				//do nothing here...
			}
		}
		else //build a triangle
		{
			count = e+1;

			std::vector<double> new_node(2);
			new_node[0] = Vave[e][0];
			new_node[1] = Vave[e][1];
			//CH 0->1->2 counterclockwise
			if( area > 0 )
				CH.push_back(new_node);
			else
				CH.insert(CH.begin()+1,new_node);

			colinear = false;

			break;
		}
	}

	//Case-1: the given cell velocities are all colinear...
	if(colinear)
	{
		if ( insideSegment(CH[0][0], CH[0][1], CH[1][0], CH[1][1], V0[0], V0[1]) )
		{
			area = getTriArea(CH[0][0], CH[0][1], CH[1][0], CH[1][1], Vp[0], Vp[1]);

			if( fabs(area) > EPS )
			{
				Vp[0] = V0[0];
				Vp[1] = V0[1];
				return 0.0;
			}
			else
			{
				if ( obtuseAngle(CH[0][0],CH[0][1], CH[1][0],CH[1][1], Vp[0],Vp[1]) )
				{
					double len0=sqrt( (CH[0][0]-V0[0])*(CH[0][0]-V0[0]) + (CH[0][1]-V0[1])*(CH[0][1]-V0[1]) );
					double lenp=sqrt( (Vp[0]-V0[0])*(Vp[0]-V0[0]) + (Vp[1]-V0[1])*(Vp[1]-V0[1]) );

					Vp[0] = CH[0][0];
					Vp[1] = CH[0][1];
					return len0/lenp;
				}
				else if ( obtuseAngle(CH[1][0],CH[1][1], CH[0][0],CH[0][1], Vp[0],Vp[1]) )
				{
					double len1=sqrt( (CH[1][0]-V0[0])*(CH[1][0]-V0[0]) + (CH[1][1]-V0[1])*(CH[1][1]-V0[1]) );
					double lenp=sqrt( (Vp[0]-V0[0])*(Vp[0]-V0[0]) + (Vp[1]-V0[1])*(Vp[1]-V0[1]) );

					Vp[0] = CH[1][0];
					Vp[1] = CH[1][1];
					return len1/lenp;
				}
				else
				{
					return 1.0;
				}
			}
		}
		else
		{
			Vp[0] = V0[0];
			Vp[1] = V0[1];
			return 0.0;
		}
	}

	//Case-2: check if the CH can be further extended using the last node...
	if( !colinear && count < neigh_cell_num )
	{
		bool face0(false), face1(false), face2(false);
		//Outward normal of CH
		//face0: 1->2
		double n0x =   CH[2][1]-CH[1][1];
		double n0y = -(CH[2][0]-CH[1][0]);

		//face1: 2->0
		double n1x =   CH[0][1]-CH[2][1];
		double n1y = -(CH[0][0]-CH[2][0]);

		//check the node v.s. face position
		if ((Vave[count][0] - CH[1][0])*n0x + (Vave[count][1] - CH[1][1])*n0y > EPS)
			face0 = true;

		if ((Vave[count][0] - CH[2][0])*n1x + (Vave[count][1] - CH[2][1])*n1y > EPS)
			face1 = true;

		if ((Vave[count][0] - CH[0][0])*n1x + (Vave[count][1] - CH[0][1])*n1y > EPS)
			face2 = true;

		//there are seven possible cases
		if (face0 && face1)
		{
			CH[2][0] = Vave[count][0];
			CH[2][1] = Vave[count][1];
		}
		else if (face0 && face2)
		{
			CH[1][0] = Vave[count][0];
			CH[1][1] = Vave[count][1];
		}
		else if (face1 && face2)
		{
			CH[0][0] = Vave[count][0];
			CH[0][1] = Vave[count][1];
		}
		else if (face0)
		{
			std::vector<double> new_node(2);
			new_node[0] = Vave[count][0];
			new_node[1] = Vave[count][1];
			CH.insert(CH.begin() + 2, new_node);
		}
		else if (face1)
		{
			std::vector<double> new_node(2);
			new_node[0] = Vave[count][0];
			new_node[1] = Vave[count][1];
			CH.insert(CH.begin() + 3, new_node);
		}
		else if (face2)
		{
			std::vector<double> new_node(2);
			new_node[0] = Vave[count][0];
			new_node[1] = Vave[count][1];
			CH.insert(CH.begin() + 1, new_node);
		}
		else
		{
			//do nothing here, the CH is not extended by the extra node
		}
	}

	//Now, we do the VIP limiter using the CH we find above...
	//(1) If V0 lies outside the CH, we set Vp=V0...
	//(2) If V0 lies inside the CH, then we check whether Vp lies inside the CH or not and limit Vp if necessary...

	double lambda(0.0);

	if ( CH.size() == 3 )
	{
		if ( insideTriCH(CH, false, V0) )
		{
			insideTriCH(CH, true, V0, Vp, lambda);
		}
		else
		{
			Vp[0] = V0[0];
			Vp[1] = V0[1];
		}
	}
	else if ( CH.size() == 4 )
	{
		if( insideQuadCH(CH, false, V0) )
		{
			insideQuadCH(CH, true, V0, Vp, lambda);
		}
		else
		{
			Vp[0] = V0[0];
			Vp[1] = V0[1];
		}
	}

	return lambda;
}

///////////////////////////////////////////////////
//some subroutines called by useViPLimiter...
///////////////////////////////////////////////////
static double getTriArea(double x0, double y0, double x1, double y1, double xp, double yp)
{
	return (xp - x0)*(yp - y1) - (xp - x1)*(yp - y0);
}

static void getPerpendFoot(double x0, double y0, double x1, double y1, double xc, double yc, double* pf)
{
	double k(0);

	k = ((xc-x0)*(x1-x0) + (yc-y0)*(y1-y0)) / ((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));

	pf[0] = x0 + k*(x1-x0);
	pf[1] = y0 + k*(y1-y0);
}

static bool obtuseAngle(double x0, double y0, double xa, double ya, double xb, double yb)
{
	if( (xa-x0)*(xb-x0) + (ya-y0)*(yb-y0) > EPS )
		return false;
	else
		return true;
}

static bool insideSegment(double x0, double y0, double x1, double y1, double xp, double yp)
{
	double area(0);

	area = getTriArea(x0, y0, x1, y1, xp, yp);

	if ( fabs(area) > EPS )
	{
		return false;
	}
	else if ( xp > std::max(x0,x1) || xp < std::min(x0,x1) )
	{
		return false;
	}

	return true;
}

static double insectionPoint(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double* Vp)
{
	double r[2], s[2], pq[2];
	r[0] = x1-x0;
	r[1] = y1-y0;
	s[0] = x3-x2;
	s[1] = y3-y2;
	pq[0] = x2-x0;
	pq[1] = y2-y0;

	double pqxs = pq[0]*s[1] - pq[1]*s[0];
	double rxs  = r[0]*s[1] - r[1]*s[0];
	double t = fabs(pqxs/rxs);

	Vp[0] = x0 + t*r[0];
	Vp[1] = y0 + t*r[1];

	return t;
}

static bool insideTriCH(std::vector<std::vector<double> >& CH, bool flag, double* Vp)
{
	bool face0(false), face1(false), face2(false);

	//face0: 1->2
	double n0x =   CH[2][1] - CH[1][1];
	double n0y = -(CH[2][0] - CH[1][0]);

	//face1: 2->0
	double n1x =   CH[0][1] - CH[2][1];
	double n1y = -(CH[0][0] - CH[2][0]);

	//face2: 0->1
	double n2x =   CH[1][1] - CH[0][1];
	double n2y = -(CH[1][0] - CH[0][0]);

	//check the node v.s. face position
	if ((Vp[0] - CH[1][0])*n0x + (Vp[1] - CH[1][1])*n0y > EPS)
		face0 = true;

	if ((Vp[0] - CH[2][0])*n1x + (Vp[1] - CH[2][1])*n1y > EPS)
		face1 = true;

	if ((Vp[0] - CH[0][0])*n2x + (Vp[1] - CH[0][1])*n2y > EPS)
		face2 = true;

	if (!flag)
	{
		if(face0 || face1 || face2)
			return false;
		else
			return true;
	}

	//if flag, we may need to limit Vp
	if (face0)
	{
		if(obtuseAngle(CH[2][0],CH[2][1], CH[1][0],CH[1][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[2][0];
			Vp[1] = CH[2][1];
		}
		else if(obtuseAngle(CH[1][0],CH[1][1], CH[2][0],CH[2][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[1][0];
			Vp[1] = CH[1][1];
		}
		else
		{
			getPerpendFoot(CH[1][0],CH[1][1],CH[2][0],CH[2][1],Vp[0],Vp[1],Vp);
		}

		return false;
	}
	else if (face1)
	{
		if(obtuseAngle(CH[2][0],CH[2][1], CH[0][0],CH[0][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[2][0];
			Vp[1] = CH[2][1];
		}
		else if(obtuseAngle(CH[0][0],CH[0][1], CH[2][0],CH[2][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[0][0];
			Vp[1] = CH[0][1];
		}
		else
		{
			getPerpendFoot(CH[2][0],CH[2][1],CH[0][0],CH[0][1],Vp[0],Vp[1],Vp);
		}

		return false;
	}
	else if (face2)
	{
		if(obtuseAngle(CH[0][0],CH[0][1], CH[1][0],CH[1][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[0][0];
			Vp[1] = CH[0][1];
		}
		else if(obtuseAngle(CH[1][0],CH[1][1], CH[0][0],CH[0][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[1][0];
			Vp[1] = CH[1][1];
		}
		else
		{
			getPerpendFoot(CH[0][0],CH[0][1],CH[1][0],CH[1][1],Vp[0],Vp[1],Vp);
		}

		return false;
	}
	else
	{
		//do nothing here...
	}

	return true;
}

static bool insideQuadCH(std::vector<std::vector<double> >& CH, bool flag, double* Vp)
{
	bool face0(false), face1(false), face2(false), face3(false);

	//face0: 0->1
	double n0x =   CH[1][1] - CH[0][1];
	double n0y = -(CH[1][0] - CH[0][0]);

	//face1: 1->2
	double n1x =   CH[2][1] - CH[1][1];
	double n1y = -(CH[2][0] - CH[1][0]);

	//face2: 2->3
	double n2x =   CH[3][1] - CH[2][1];
	double n2y = -(CH[3][0] - CH[2][0]);

	//face3: 3->0
	double n3x =   CH[0][1] - CH[3][1];
	double n3y = -(CH[0][0] - CH[3][0]);

	//check the node v.s. face position
	if ((Vp[0] - CH[0][0])*n0x + (Vp[1] - CH[0][1])*n0y > EPS)
		face0 = true;

	if ((Vp[0] - CH[1][0])*n1x + (Vp[1] - CH[1][1])*n1y > EPS)
		face1 = true;

	if ((Vp[0] - CH[2][0])*n2x + (Vp[1] - CH[2][1])*n2y > EPS)
		face2 = true;

	if ((Vp[0] - CH[3][0])*n3x + (Vp[1] - CH[3][1])*n3y > EPS)
		face3 = true;

	if (!flag)
	{
		if(face0 || face1 || face2 || face3)
			return false;
		else
			return true;
	}

	//if flag, we may need to limit Vp
	if (face0)
	{
		if(obtuseAngle(CH[0][0],CH[0][1], CH[1][0],CH[1][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[0][0];
			Vp[1] = CH[0][1];
		}
		else if(obtuseAngle(CH[1][0],CH[1][1], CH[0][0],CH[0][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[1][0];
			Vp[1] = CH[1][1];
		}
		else
		{
			getPerpendFoot(CH[0][0],CH[0][1],CH[1][0],CH[1][1],Vp[0],Vp[1],Vp);
		}
		return false;
	}
	else if (face1)
	{
		if(obtuseAngle(CH[2][0],CH[2][1], CH[1][0],CH[1][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[2][0];
			Vp[1] = CH[2][1];
		}
		else if(obtuseAngle(CH[1][0],CH[1][1], CH[2][0],CH[2][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[1][0];
			Vp[1] = CH[1][1];
		}
		else
		{
			getPerpendFoot(CH[1][0],CH[1][1],CH[2][0],CH[2][1],Vp[0],Vp[1],Vp);
		}
		return false;
	}
	else if (face2)
	{
		if(obtuseAngle(CH[2][0],CH[2][1], CH[3][0],CH[3][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[2][0];
			Vp[1] = CH[2][1];
		}
		else if(obtuseAngle(CH[3][0],CH[3][1], CH[2][0],CH[2][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[3][0];
			Vp[1] = CH[3][1];
		}
		else
		{
			getPerpendFoot(CH[3][0],CH[3][1],CH[2][0],CH[2][1],Vp[0],Vp[1],Vp);
		}
		return false;
	}
	else if (face3)
	{
		if(obtuseAngle(CH[0][0],CH[0][1], CH[3][0],CH[3][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[0][0];
			Vp[1] = CH[0][1];
		}
		else if(obtuseAngle(CH[3][0],CH[3][1], CH[0][0],CH[0][1], Vp[0],Vp[1]))
		{
			Vp[0] = CH[3][0];
			Vp[1] = CH[3][1];
		}
		else
		{
			getPerpendFoot(CH[3][0],CH[3][1],CH[0][0],CH[0][1],Vp[0],Vp[1],Vp);
		}
		return false;
	}
	else
	{
		//do nothing here...
	}

	return true;
}

static bool insideTriCH(std::vector<std::vector<double> >& CH, bool flag, double* V0, double* Vp, double& lambda)
{
	bool face0(false), face1(false), face2(false);

	//face0: 1->2
	double n0x =   CH[2][1] - CH[1][1];
	double n0y = -(CH[2][0] - CH[1][0]);

	//face1: 2->0
	double n1x =   CH[0][1] - CH[2][1];
	double n1y = -(CH[0][0] - CH[2][0]);

	//face2: 0->1
	double n2x =   CH[1][1] - CH[0][1];
	double n2y = -(CH[1][0] - CH[0][0]);

	//check the node v.s. face position
	if ( (Vp[0] - CH[1][0])*n0x + (Vp[1] - CH[1][1])*n0y > EPS )
		face0 = true;

	if ( (Vp[0] - CH[2][0])*n1x + (Vp[1] - CH[2][1])*n1y > EPS )
		face1 = true;

	if ( (Vp[0] - CH[0][0])*n2x + (Vp[1] - CH[0][1])*n2y > EPS )
		face2 = true;

	//if unflag, we only need to check whether the point is inside or not
	if (!flag)
	{
		if(face0 || face1 || face2)
			return false;
		else
			return true;
	}

	//if flag, we may need to limit Vp and compute lambda
	if( !face0 && !face1 && !face2 ) //Vp inside the given CH
	{
		lambda = 1.0;
		return true;
	}
	else //otherwise, we need to check which face the line V0-Vp will insect with, and comput the insection point and lambda
	{
		double area0 = getTriArea(V0[0],V0[1],CH[0][0],CH[0][1],Vp[0],Vp[1]);
		double area1 = getTriArea(V0[0],V0[1],CH[1][0],CH[1][1],Vp[0],Vp[1]);
		double area2 = getTriArea(V0[0],V0[1],CH[2][0],CH[2][1],Vp[0],Vp[1]);

		double len = sqrt( (Vp[0]-V0[0])*(Vp[0]-V0[0]) + (Vp[1]-V0[1])*(Vp[1]-V0[1]) );

		if( area0 > EPS && area1 < EPS )     //insect with face2
		{
			lambda = insectionPoint(V0[0],V0[1], Vp[0],Vp[1], CH[0][0],CH[0][1], CH[1][0],CH[1][1], Vp);
		}
		else if( area1 > EPS && area2 < EPS ) //insect with face0
		{
			lambda = insectionPoint(V0[0],V0[1], Vp[0],Vp[1], CH[1][0],CH[1][1], CH[2][0],CH[2][1], Vp);
		}
		else if( area2 > EPS && area0 < EPS ) //insect with face1
		{
			lambda = insectionPoint(V0[0],V0[1], Vp[0],Vp[1], CH[2][0],CH[2][1], CH[0][0],CH[0][1], Vp);
		}
		else if( fabs(area0) < EPS )//along V0-node0
		{
			Vp[0] = CH[0][0];
			Vp[1] = CH[0][1];

			lambda = sqrt( (CH[0][0]-V0[0])*(CH[0][0]-V0[0]) + (CH[0][1]-V0[1])*(CH[0][1]-V0[1]) )/len;
		}
		else if( fabs(area1) < EPS )//along V0-node1
		{
			Vp[0] = CH[1][0];
			Vp[1] = CH[1][1];

			lambda = sqrt( (CH[1][0]-V0[0])*(CH[1][0]-V0[0]) + (CH[1][1]-V0[1])*(CH[1][1]-V0[1]) )/len;
		}
		else if( fabs(area2) < EPS )//along V0-node2
		{
			Vp[0] = CH[2][0];
			Vp[1] = CH[2][1];

			lambda = sqrt( (CH[2][0]-V0[0])*(CH[2][0]-V0[0]) + (CH[2][1]-V0[1])*(CH[2][1]-V0[1]) )/len;
		}
		else
		{
			std::cout<<"Error: it should not be here for quadrilateral CH..."<<std::endl;
			lambda = 1.0;
			//exit(1);
		}

		return false;
	}
}

static bool insideQuadCH(std::vector<std::vector<double> >& CH, bool flag, double* V0, double* Vp, double& lambda)
{
	bool face0(false), face1(false), face2(false), face3(false);

	//face0: 0->1
	double n0x =   CH[1][1] - CH[0][1];
	double n0y = -(CH[1][0] - CH[0][0]);

	//face1: 1->2
	double n1x =   CH[2][1] - CH[1][1];
	double n1y = -(CH[2][0] - CH[1][0]);

	//face2: 2->3
	double n2x =   CH[3][1] - CH[2][1];
	double n2y = -(CH[3][0] - CH[2][0]);

	//face3: 3->0
	double n3x =   CH[0][1] - CH[3][1];
	double n3y = -(CH[0][0] - CH[3][0]);

	//check the node v.s. face position
	if ( (Vp[0] - CH[0][0])*n0x + (Vp[1] - CH[0][1])*n0y > EPS )
		face0 = true;

	if ( (Vp[0] - CH[1][0])*n1x + (Vp[1] - CH[1][1])*n1y > EPS )
		face1 = true;

	if ( (Vp[0] - CH[2][0])*n2x + (Vp[1] - CH[2][1])*n2y > EPS )
		face2 = true;

	if ( (Vp[0] - CH[3][0])*n3x + (Vp[1] - CH[3][1])*n3y > EPS )
		face3 = true;

	//if unflag, we only need to check whether the point is inside or not
	if ( !flag )
	{
		if( face0 || face1 || face2 || face3 )
			return false;
		else
			return true;
	}

	//if flag, we may need to limit Vp and compute lambda
	if( !face0 && !face1 && !face2 && !face3 ) //Vp inside the given CH
	{
		lambda = 1.0;
		return true;
	}
	else
	{
		double area0 = getTriArea(V0[0],V0[1],CH[0][0],CH[0][1],Vp[0],Vp[1]);
		double area1 = getTriArea(V0[0],V0[1],CH[1][0],CH[1][1],Vp[0],Vp[1]);
		double area2 = getTriArea(V0[0],V0[1],CH[2][0],CH[2][1],Vp[0],Vp[1]);
		double area3 = getTriArea(V0[0],V0[1],CH[3][0],CH[3][1],Vp[0],Vp[1]);

		double len = sqrt( (Vp[0]-V0[0])*(Vp[0]-V0[0]) + (Vp[1]-V0[1])*(Vp[1]-V0[1]) );

		if( area0 > EPS && area1 < EPS )     //insect with face0
		{
			lambda = insectionPoint(V0[0],V0[1], Vp[0],Vp[1], CH[0][0],CH[0][1], CH[1][0],CH[1][1], Vp);
		}
		else if( area1 > EPS && area2 < EPS ) //insect with face1
		{
			lambda = insectionPoint(V0[0],V0[1], Vp[0],Vp[1], CH[1][0],CH[1][1], CH[2][0],CH[2][1], Vp);
		}
		else if( area2 > EPS && area0 < EPS ) //insect with face2
		{
			lambda = insectionPoint(V0[0],V0[1], Vp[0],Vp[1], CH[2][0],CH[2][1], CH[3][0],CH[3][1], Vp);
		}
		else if( area3 > EPS && area0 < EPS)  //insect with face3
		{
			lambda = insectionPoint(V0[0],V0[1], Vp[0],Vp[1], CH[3][0],CH[3][1], CH[0][0],CH[0][1], Vp);
		}
		else if( fabs(area0) < EPS )//along V0-node0
		{
			Vp[0] = CH[0][0];
			Vp[1] = CH[0][1];

			lambda = sqrt( (CH[0][0]-V0[0])*(CH[0][0]-V0[0]) + (CH[0][1]-V0[1])*(CH[0][1]-V0[1]) )/len;
		}
		else if( fabs(area1) < EPS )//along V0-node1
		{
			Vp[0] = CH[1][0];
			Vp[1] = CH[1][1];

			lambda = sqrt( (CH[1][0]-V0[0])*(CH[1][0]-V0[0]) + (CH[1][1]-V0[1])*(CH[1][1]-V0[1]) )/len;
		}
		else if( fabs(area2) < EPS )//along V0-node2
		{
			Vp[0] = CH[2][0];
			Vp[1] = CH[2][1];

			lambda = sqrt( (CH[2][0]-V0[0])*(CH[2][0]-V0[0]) + (CH[2][1]-V0[1])*(CH[2][1]-V0[1]) )/len;
		}
		else if( fabs(area3) < EPS )//along V0-node3
		{
			Vp[0] = CH[3][0];
			Vp[1] = CH[3][1];

			lambda = sqrt( (CH[3][0]-V0[0])*(CH[3][0]-V0[0]) + (CH[3][1]-V0[1])*(CH[3][1]-V0[1]) )/len;
		}
		else
		{
			std::cout<<"Error: it should not be here for triangualr CH..."<<std::endl;
			lambda = 1.0;
			//exit(1);
		}

		return false;
	}
}
////////////////////////////////////////////////////
