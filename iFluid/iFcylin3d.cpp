/*******************************************************************
 * 			iFcartsn3d.cpp	
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"
#include <vector>
//--------------------------------------------------------------------------
// 		   Incompress_Solver_Smooth_3D_Cylindrical
//--------------------------------------------------------------------------
//
//
//
//
//-----------------------------------------------------------------------------------------------------
//  utility functions related to the computation of advection source term
//-------------------------------------------------------------------------------------------------------


void Incompress_Solver_Smooth_3D_Cylindrical::compAdvectionTerm_decoupled(void)
{
        int index;
    	int i,j,k,icoords[MAXD];
    	int I;

	setIndexMap();

    	for (k = kmin; k <= kmax; k++)
    	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	     I  = ijk_to_I[i][j][k];
	     if (I == -1) 
	     {
		 continue;
	     }

	     index = d_index3d(i,j,k,top_gmax);

	     icoords[0] = i;
	     icoords[1] = j;
	     icoords[2] = k;

	     double convectionTerm[3];
	     getAdvectionTerm_decoupled(icoords, convectionTerm);

	     cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	     cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	     cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}

}

void Incompress_Solver_Smooth_3D_Cylindrical::compAdvectionTerm_coupled(void)
{
        int index;
    	int i,j,k,icoords[MAXD];
    	int I;

	setIndexMap();

    	for (k = kmin; k <= kmax; k++)
    	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	     I  = ijk_to_I[i][j][k];
	     if (I == -1) 
	     {
		 continue;
	     }

	     index = d_index3d(i,j,k,top_gmax);

	     icoords[0] = i;
	     icoords[1] = j;
	     icoords[2] = k;

	     double convectionTerm[3];
	     getAdvectionTerm_coupled(icoords, convectionTerm);

	     cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	     cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	     cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}

}


void Incompress_Solver_Smooth_3D_Cylindrical::getAdvectionTerm_decoupled(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west, state_east, state_south, state_north, state_lower, state_upper;

    double dtheta = top_h[0];
    double dz 	  = top_h[1];
    double dr     = top_h[2];
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double r = cell_center[index].m_coords[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,EAST,sl);
	getFaceVelocity_middleStep(icoords,WEST,sr);
	getRiemannSolution(COORD_X,sl,sr,state_west);
    }
    else
	state_west = sl;

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,WEST,sr);
    	getFaceVelocity_middleStep(icoords,EAST,sl);
	getRiemannSolution(COORD_X,sl,sr,state_east);
    }
    else
	state_east = sr;

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,NORTH,sl);
	getFaceVelocity_middleStep(icoords,SOUTH,sr);
	getRiemannSolution(COORD_Y,sl,sr,state_south);
    }
    else
	state_south = sl;

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,SOUTH,sr);
	getFaceVelocity_middleStep(icoords,NORTH,sl);
	getRiemannSolution(COORD_Y,sl,sr,state_north);
    }
    else
	state_north = sr;

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep(ICoords,UPPER,sl);
	getFaceVelocity_middleStep(icoords,LOWER,sr);
	getRiemannSolution(COORD_Z,sl,sr,state_lower);
    }
    else
	state_lower = sl;

    // UPPER
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep(ICoords,LOWER,sr);
	getFaceVelocity_middleStep(icoords,UPPER,sl);
	getRiemannSolution(COORD_Z,sl,sr,state_upper);
    }
    else
	state_upper = sr;

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[0]-state_west.m_U[0])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[0]-state_south.m_U[0])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[0]-state_lower.m_U[0])/dr +
	    1/r * 1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * 1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]);

    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[1]-state_west.m_U[1])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[1]-state_south.m_U[1])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[1]-state_lower.m_U[1])/dr;

    convectionTerm[2] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[2]-state_west.m_U[2])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[2]-state_south.m_U[2])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[2]-state_lower.m_U[2])/dr -
	    1/r * sqr(1.0/2*(state_east.m_U[0]+state_west.m_U[0]));
}


void Incompress_Solver_Smooth_3D_Cylindrical::getAdvectionTerm_coupled(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west, state_east, state_south, state_north, state_lower, state_upper;

    double dtheta = top_h[0];
    double dz 	  = top_h[1];
    double dr     = top_h[2];
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double r = cell_center[index].m_coords[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,EAST,sl);
	getFaceVelocity_middleStep_coupled(icoords,WEST,sr);
	getRiemannSolution(COORD_X,sl,sr,state_west);
    }
    else
	state_west = sl;

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,WEST,sr);
    	getFaceVelocity_middleStep_coupled(icoords,EAST,sl);
	getRiemannSolution(COORD_X,sl,sr,state_east);
    }
    else
	state_east = sr;

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,NORTH,sl);
	getFaceVelocity_middleStep_coupled(icoords,SOUTH,sr);
	getRiemannSolution(COORD_Y,sl,sr,state_south);
    }
    else
	state_south = sl;

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,SOUTH,sr);
	getFaceVelocity_middleStep_coupled(icoords,NORTH,sl);
	getRiemannSolution(COORD_Y,sl,sr,state_north);
    }
    else
	state_north = sr;

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_coupled(ICoords,UPPER,sl);
	getFaceVelocity_middleStep_coupled(icoords,LOWER,sr);
	getRiemannSolution(COORD_Z,sl,sr,state_lower);
    }
    else
	state_lower = sl;

    // UPPER
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_coupled(ICoords,LOWER,sr);
	getFaceVelocity_middleStep_coupled(icoords,UPPER,sl);
	getRiemannSolution(COORD_Z,sl,sr,state_upper);
    }
    else
	state_upper = sr;

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[0]-state_west.m_U[0])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[0]-state_south.m_U[0])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[0]-state_lower.m_U[0])/dr +
	    1/r * 1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * 1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]);

    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[1]-state_west.m_U[1])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[1]-state_south.m_U[1])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[1]-state_lower.m_U[1])/dr;

    convectionTerm[2] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[2]-state_west.m_U[2])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[2]-state_south.m_U[2])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[2]-state_lower.m_U[2])/dr -
	    1/r * sqr(1.0/2*(state_east.m_U[0]+state_west.m_U[0]));
}


void Incompress_Solver_Smooth_3D_Cylindrical::getFaceVelocity_middleStep(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    L_STATE state;
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state = cell_center[index].m_state;
    double r = cell_center[index].m_coords[2];
    double rho = cell_center[index].m_state.m_rho;

    double dx = 0.0, dy = 0.0 ,dz = 0.0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};

    switch(dir)
    {
    case WEST:
	dx = -top_h[0];
	break;
    case EAST:
	dx =  top_h[0];
	break;
    case SOUTH:
	dy = -top_h[1];
	break;
    case NORTH:
	dy =  top_h[1];
	break;
    case LOWER:
	dz = -top_h[2];
	break;
    case UPPER:
	dz =  top_h[2];
	break;
    default:
	assert(false);
    }

    state_face.m_U[0] = state.m_U[0];
    state_face.m_U[1] = state.m_U[1];
    state_face.m_U[2] = state.m_U[2];

//    return;

    getLimitedSlope(icoords,COORD_X,slope_x_limited);
    getLimitedSlope(icoords,COORD_Y,slope_y_limited);
    getLimitedSlope(icoords,COORD_Z,slope_z_limited);

    // dx/2, dy/2, dz/2
    state_face.m_U[0] += dx/2 * slope_x_limited[0] + dy/2 * slope_y_limited[0] + dz/2 * slope_z_limited[0];
    state_face.m_U[1] += dx/2 * slope_x_limited[1] + dy/2 * slope_y_limited[1] + dz/2 * slope_z_limited[1];
    state_face.m_U[2] += dx/2 * slope_x_limited[2] + dy/2 * slope_y_limited[2] + dz/2 * slope_z_limited[2];

    //    return;
    // dt/2
    double diffusion[3];
    getDifffusion(icoords,diffusion);
    state_face.m_U[0] += m_dt/2 *
	    ( diffusion[0]/rho - (1/r*state.m_U[0]*slope_x_limited[0] + state.m_U[1]*slope_y_limited[0] + state.m_U[2]*slope_z_limited[0] + 1/r*(state.m_U[0]*state.m_U[2]) ) - cell_center[index].m_state.grad_q[0]/rho );
    state_face.m_U[1] += m_dt/2 *
	    ( diffusion[1]/rho - (1/r*state.m_U[0]*slope_x_limited[1] + state.m_U[1]*slope_y_limited[1] + state.m_U[2]*slope_z_limited[1]                                   ) - cell_center[index].m_state.grad_q[1]/rho );
    state_face.m_U[2] += m_dt/2 *
	    ( diffusion[2]/rho - (1/r*state.m_U[0]*slope_x_limited[2] + state.m_U[1]*slope_y_limited[2] + state.m_U[2]*slope_z_limited[2] - 1/r*(state.m_U[0]*state.m_U[0]) ) - cell_center[index].m_state.grad_q[2]/rho );

    //    return;
    // rhs
    double coords[3];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    // the source term does not contain density.
    computeSourceTerm_Adv(coords, source_term);

    state_face.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_face.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_face.m_U[2] += m_dt/2 * source_term.m_U[2];
}

void Incompress_Solver_Smooth_3D_Cylindrical::getFaceVelocity_middleStep_coupled(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    L_STATE state;
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state = cell_center[index].m_state;
    double r = cell_center[index].m_coords[2];
    double rho = cell_center[index].m_state.m_rho;

    double dx = 0.0, dy = 0.0 ,dz = 0.0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};

    switch(dir)
    {
    case WEST:
	dx = -top_h[0];
	break;
    case EAST:
	dx =  top_h[0];
	break;
    case SOUTH:
	dy = -top_h[1];
	break;
    case NORTH:
	dy =  top_h[1];
	break;
    case LOWER:
	dz = -top_h[2];
	break;
    case UPPER:
	dz =  top_h[2];
	break;
    default:
	assert(false);
    }

    state_face.m_U[0] = state.m_U[0];
    state_face.m_U[1] = state.m_U[1];
    state_face.m_U[2] = state.m_U[2];

//    return;

    getLimitedSlope(icoords,COORD_X,slope_x_limited);
    getLimitedSlope(icoords,COORD_Y,slope_y_limited);
    getLimitedSlope(icoords,COORD_Z,slope_z_limited);

    // dx/2, dy/2, dz/2
    state_face.m_U[0] += dx/2 * slope_x_limited[0] + dy/2 * slope_y_limited[0] + dz/2 * slope_z_limited[0];
    state_face.m_U[1] += dx/2 * slope_x_limited[1] + dy/2 * slope_y_limited[1] + dz/2 * slope_z_limited[1];
    state_face.m_U[2] += dx/2 * slope_x_limited[2] + dy/2 * slope_y_limited[2] + dz/2 * slope_z_limited[2];

    //    return;
    // dt/2
    double diffusion[3];
    getDiffusion_coupled(icoords,diffusion);
    state_face.m_U[0] += m_dt/2 *
	    ( diffusion[0]/rho - (1/r*state.m_U[0]*slope_x_limited[0] + state.m_U[1]*slope_y_limited[0] + state.m_U[2]*slope_z_limited[0] + 1/r*(state.m_U[0]*state.m_U[2]) ) - cell_center[index].m_state.grad_q[0]/rho );
    state_face.m_U[1] += m_dt/2 *
	    ( diffusion[1]/rho - (1/r*state.m_U[0]*slope_x_limited[1] + state.m_U[1]*slope_y_limited[1] + state.m_U[2]*slope_z_limited[1]                                   ) - cell_center[index].m_state.grad_q[1]/rho );
    state_face.m_U[2] += m_dt/2 *
	    ( diffusion[2]/rho - (1/r*state.m_U[0]*slope_x_limited[2] + state.m_U[1]*slope_y_limited[2] + state.m_U[2]*slope_z_limited[2] - 1/r*(state.m_U[0]*state.m_U[0]) ) - cell_center[index].m_state.grad_q[2]/rho );

    //    return;
    // rhs
    double coords[3];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    // the source term does not contain density.
    computeSourceTerm_Adv(coords, source_term);

    state_face.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_face.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_face.m_U[2] += m_dt/2 * source_term.m_U[2];
}


/**
* compute mu * (Uxx+Uyy+Uzz)
* @param icoords
* @param diffusion
* @param gradP
*/

void Incompress_Solver_Smooth_3D_Cylindrical::getDifffusion(
	int *icoords,
	double diffusion[3])
{
    double dU2_theta[3], dU2_z[3], dU2_r[3];
    double dU_theta[3], dU_z[3], dU_r[3];
    getDU2(icoords,COORD_X, dU2_theta, dU_theta);
    getDU2(icoords,COORD_Y, dU2_z,     dU_z);
    getDU2(icoords,COORD_Z, dU2_r,     dU_r);

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double mu = cell_center[index].m_state.m_mu;

    double r = cell_center[index].m_coords[2];
    L_STATE U = cell_center[index].m_state;

    diffusion[0] = mu * (1/(r*r)*dU2_theta[0] + 1/r*dU_r[0] + dU2_z[0] + dU2_r[0]);
    diffusion[0]+= mu * ( 2/(r*r)*dU_theta[2] - 1/(r*r)*U.m_U[0]);

    diffusion[1] = mu * (1/(r*r)*dU2_theta[1] + 1/r*dU_r[1] + dU2_z[1] + dU2_r[1]);

    diffusion[2] = mu * (1/(r*r)*dU2_theta[2] + 1/r*dU_r[2] + dU2_z[2] + dU2_r[2]);
    diffusion[2]+= mu * (-2/(r*r)*dU_theta[0] - 1/(r*r)*U.m_U[2]);
}

void Incompress_Solver_Smooth_3D_Cylindrical::getDiffusion_coupled(
	int *icoords,
	double diffusion[3])
{
    int index, index_nb[18];
    double mu[6], mu_edge[6], mu0;
    double r0, rr, r[6], r_edge[6];
    L_STATE Unb;
    double U_nb[3][18], U_center[3], U_face[3][6];
    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    bool bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];
    double dtheta,dz,dr;

    dtheta = top_h[0];
    dz = top_h[1];
    dr = top_h[2];

    diffusion[0] = 0.0;
    diffusion[1] = 0.0;
    diffusion[2] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];

    index = d_index3d(i,j,k,top_gmax);

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);
    
    //theta-z cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	
    //z-r cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //theta-r cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    mu0 = cell_center[index].m_state.m_mu;
    r0 = cell_center[index].m_coords[2];
    rr = r0*r0;

    for (int l = 0; l < 3; l++)
	U_center[l] = cell_center[index].m_state.m_U[l];

    for (nb = 0; nb < 6; nb++)
    {
	bNoBoundary[nb] = getNeighborOrBoundaryState(icoords, dir[nb], Unb, m_t_old);
	U_nb[0][nb] = Unb.m_U[0];
	U_nb[1][nb] = Unb.m_U[1];
	U_nb[2][nb] = Unb.m_U[2];

	if(!bNoBoundary[nb])
	{
	    mu[nb] = mu0;
	    mu_edge[nb] = mu0;
	    U_face[0][nb] = U_nb[0][nb];
	    U_face[1][nb] = U_nb[1][nb];
	    U_face[2][nb] = U_nb[2][nb];
	}
	else
	{
	    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
	    mu_edge[nb] = 0.5*(mu[nb] + mu0);
	    U_face[0][nb] = 0.5*(U_nb[0][nb] + U_center[0]);
	    U_face[1][nb] = 0.5*(U_nb[1][nb] + U_center[1]);
	    U_face[2][nb] = 0.5*(U_nb[2][nb] + U_center[2]);
	}
	r[nb] = cell_center[index_nb[nb]].m_coords[2];
	r_edge[nb] = 0.5 * (r[nb] + r0);
    }

    //traverse the corners on 3 cut planes

    //corner (i-1/2,j-1/2,k)

    if (!bNoBoundary[0] && bNoBoundary[2])
    {
	U_nb[0][6] = U_nb[0][0];
	U_nb[1][6] = U_nb[1][0];
	U_nb[2][6] = U_nb[2][0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[2])
    {
	U_nb[0][6] = U_nb[0][2];
	U_nb[1][6] = U_nb[1][2];
	U_nb[2][6] = U_nb[2][2];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[2])
    {
	U_nb[0][6] = U_nb[0][0];
	U_nb[1][6] = U_nb[1][0];
	U_nb[2][6] = U_nb[2][0];
    }
    else
    {
	U_nb[0][6] = (U_nb[0][0]+U_nb[0][2]+cell_center[index_nb[6]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][6] = (U_nb[1][0]+U_nb[1][2]+cell_center[index_nb[6]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][6] = (U_nb[2][0]+U_nb[2][2]+cell_center[index_nb[6]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i+1/2,j-1/2,k)

    if (!bNoBoundary[1] && bNoBoundary[2])
    {
	U_nb[0][7] = U_nb[0][1];
	U_nb[1][7] = U_nb[1][1];
	U_nb[2][7] = U_nb[2][1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[2])
    {
	U_nb[0][7] = U_nb[0][2];
	U_nb[1][7] = U_nb[1][2];
	U_nb[2][7] = U_nb[2][2];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[2])
    {
	U_nb[0][7] = U_nb[0][1];
	U_nb[1][7] = U_nb[1][1];
	U_nb[2][7] = U_nb[2][1];
    }
    else
    {
	U_nb[0][7] = (U_nb[0][1]+U_nb[0][2]+cell_center[index_nb[7]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][7] = (U_nb[1][1]+U_nb[1][2]+cell_center[index_nb[7]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][7] = (U_nb[2][1]+U_nb[2][2]+cell_center[index_nb[7]].m_state.m_U[2]+U_center[2])/4.0;
    }


    //corner (i+1/2,j+1/2,k)

    if (!bNoBoundary[1] && bNoBoundary[3])
    {
	U_nb[0][8] = U_nb[0][1];
	U_nb[1][8] = U_nb[1][1];
	U_nb[2][8] = U_nb[2][1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[3])
    {
	U_nb[0][8] = U_nb[0][3];
	U_nb[1][8] = U_nb[1][3];
	U_nb[2][8] = U_nb[2][3];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[3])
    {
	U_nb[0][8] = U_nb[0][1];
	U_nb[1][8] = U_nb[1][1];
	U_nb[2][8] = U_nb[2][1];
    }
    else
    {
	U_nb[0][8] = (U_nb[0][1]+U_nb[0][3]+cell_center[index_nb[8]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][8] = (U_nb[1][1]+U_nb[1][3]+cell_center[index_nb[8]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][8] = (U_nb[2][1]+U_nb[2][3]+cell_center[index_nb[8]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i-1/2,j+1/2,k)

    if (!bNoBoundary[0] && bNoBoundary[3])
    {
	U_nb[0][9] = U_nb[0][0];
	U_nb[1][9] = U_nb[1][0];
	U_nb[2][9] = U_nb[2][0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[3])
    {
	U_nb[0][9] = U_nb[0][3];
	U_nb[1][9] = U_nb[1][3];
	U_nb[2][9] = U_nb[2][3];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[3])
    {
	U_nb[0][9] = U_nb[0][0];
	U_nb[1][9] = U_nb[1][0];
	U_nb[2][9] = U_nb[2][0];
    }
    else
    {
	U_nb[0][9] = (U_nb[0][0]+U_nb[0][3]+cell_center[index_nb[9]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][9] = (U_nb[1][0]+U_nb[1][3]+cell_center[index_nb[9]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][9] = (U_nb[2][0]+U_nb[2][3]+cell_center[index_nb[9]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i,j-1/2,k-1/2)

    if (!bNoBoundary[2] && bNoBoundary[4])
    {
	U_nb[0][10] = U_nb[0][2];
	U_nb[1][10] = U_nb[1][2];
	U_nb[2][10] = U_nb[2][2];
    }
    else if(bNoBoundary[2] && !bNoBoundary[4])
    {
	U_nb[0][10] = U_nb[0][4];
	U_nb[1][10] = U_nb[1][4];
	U_nb[2][10] = U_nb[2][4];
    }
    else if(!bNoBoundary[2] && !bNoBoundary[4])
    {
	U_nb[0][10] = U_nb[0][2];
	U_nb[1][10] = U_nb[1][2];
	U_nb[2][10] = U_nb[2][2];
    }
    else
    {
	U_nb[0][10] = (U_nb[0][2]+U_nb[0][4]+cell_center[index_nb[10]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][10] = (U_nb[1][2]+U_nb[1][4]+cell_center[index_nb[10]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][10] = (U_nb[2][2]+U_nb[2][4]+cell_center[index_nb[10]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i,j+1/2,k-1/2)

    if (!bNoBoundary[3] && bNoBoundary[4])
    {
	U_nb[0][11] = U_nb[0][3];
	U_nb[1][11] = U_nb[1][3];
	U_nb[2][11] = U_nb[2][3];
    }
    else if(bNoBoundary[3] && !bNoBoundary[4])
    {
	U_nb[0][11] = U_nb[0][4];
	U_nb[1][11] = U_nb[1][4];
	U_nb[2][11] = U_nb[2][4];
    }
    else if(!bNoBoundary[3] && !bNoBoundary[4])
    {
	U_nb[0][11] = U_nb[0][3];
	U_nb[1][11] = U_nb[1][3];
	U_nb[2][11] = U_nb[2][3];
    }
    else
    {
	U_nb[0][11] = (U_nb[0][3]+U_nb[0][4]+cell_center[index_nb[11]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][11] = (U_nb[1][3]+U_nb[1][4]+cell_center[index_nb[11]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][11] = (U_nb[2][3]+U_nb[2][4]+cell_center[index_nb[11]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i,j+1/2,k+1/2)

    if (!bNoBoundary[3] && bNoBoundary[5])
    {
	U_nb[0][12] = U_nb[0][3];
	U_nb[1][12] = U_nb[1][3];
	U_nb[2][12] = U_nb[2][3];
    }
    else if(bNoBoundary[3] && !bNoBoundary[5])
    {
	U_nb[0][12] = U_nb[0][5];
	U_nb[1][12] = U_nb[1][5];
	U_nb[2][12] = U_nb[2][5];
    }
    else if(!bNoBoundary[3] && !bNoBoundary[5])
    {
	U_nb[0][12] = U_nb[0][3];
	U_nb[1][12] = U_nb[1][3];
	U_nb[2][12] = U_nb[2][3];
    }
    else
    {
	U_nb[0][12] = (U_nb[0][3]+U_nb[0][5]+cell_center[index_nb[12]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][12] = (U_nb[1][3]+U_nb[1][5]+cell_center[index_nb[12]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][12] = (U_nb[2][3]+U_nb[2][5]+cell_center[index_nb[12]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i,j-1/2,k+1/2)

    if (!bNoBoundary[2] && bNoBoundary[5])
    {
	U_nb[0][13] = U_nb[0][2];
	U_nb[1][13] = U_nb[1][2];
	U_nb[2][13] = U_nb[2][2];
    }
    else if(bNoBoundary[2] && !bNoBoundary[5])
    {
	U_nb[0][13] = U_nb[0][5];
	U_nb[1][13] = U_nb[1][5];
	U_nb[2][13] = U_nb[2][5];
    }
    else if(!bNoBoundary[2] && !bNoBoundary[5])
    {
	U_nb[0][13] = U_nb[0][2];
	U_nb[1][13] = U_nb[1][2];
	U_nb[2][13] = U_nb[2][2];
    }
    else
    {
	U_nb[0][13] = (U_nb[0][2]+U_nb[0][5]+cell_center[index_nb[13]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][13] = (U_nb[1][2]+U_nb[1][5]+cell_center[index_nb[13]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][13] = (U_nb[2][2]+U_nb[2][5]+cell_center[index_nb[13]].m_state.m_U[2]+U_center[2])/4.0;
    }


    //corner (i-1/2,j,k-1/2)

    if (!bNoBoundary[0] && bNoBoundary[4])
    {
	U_nb[0][14] = U_nb[0][0];
	U_nb[1][14] = U_nb[1][0];
	U_nb[2][14] = U_nb[2][0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[4])
    {
	U_nb[0][14] = U_nb[0][4];
	U_nb[1][14] = U_nb[1][4];
	U_nb[2][14] = U_nb[2][4];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[4])
    {
	U_nb[0][14] = U_nb[0][0];
	U_nb[1][14] = U_nb[1][0];
	U_nb[2][14] = U_nb[2][0];
    }
    else
    {
	U_nb[0][14] = (U_nb[0][0]+U_nb[0][4]+cell_center[index_nb[14]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][14] = (U_nb[1][0]+U_nb[1][4]+cell_center[index_nb[14]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][14] = (U_nb[2][0]+U_nb[2][4]+cell_center[index_nb[14]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i+1/2,j,k-1/2)

    if (!bNoBoundary[1] && bNoBoundary[4])
    {
	U_nb[0][15] = U_nb[0][1];
	U_nb[1][15] = U_nb[1][1];
	U_nb[2][15] = U_nb[2][1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[4])
    {
	U_nb[0][15] = U_nb[0][4];
	U_nb[1][15] = U_nb[1][4];
	U_nb[2][15] = U_nb[2][4];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[4])
    {
	U_nb[0][15] = U_nb[0][1];
	U_nb[1][15] = U_nb[1][1];
	U_nb[2][15] = U_nb[2][1];
    }
    else
    {
	U_nb[0][15] = (U_nb[0][1]+U_nb[0][4]+cell_center[index_nb[15]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][15] = (U_nb[1][1]+U_nb[1][4]+cell_center[index_nb[15]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][15] = (U_nb[2][1]+U_nb[2][4]+cell_center[index_nb[15]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i+1/2,j,k+1/2)

    if (!bNoBoundary[1] && bNoBoundary[5])
    {
	U_nb[0][16] = U_nb[0][1];
	U_nb[1][16] = U_nb[1][1];
	U_nb[2][16] = U_nb[2][1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[5])
    {
	U_nb[0][16] = U_nb[0][5];
	U_nb[1][16] = U_nb[1][5];
	U_nb[2][16] = U_nb[2][5];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[5])
    {
	U_nb[0][16] = U_nb[0][1];
	U_nb[1][16] = U_nb[1][1];
	U_nb[2][16] = U_nb[2][1];
    }
    else
    {
	U_nb[0][16] = (U_nb[0][1]+U_nb[0][5]+cell_center[index_nb[16]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][16] = (U_nb[1][1]+U_nb[1][5]+cell_center[index_nb[16]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][16] = (U_nb[2][1]+U_nb[2][5]+cell_center[index_nb[16]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i-1/2,j,k+1/2)

    if (!bNoBoundary[0] && bNoBoundary[5])
    {
	U_nb[0][17] = U_nb[0][0];
	U_nb[1][17] = U_nb[1][0];
	U_nb[2][17] = U_nb[2][0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[5])
    {
	U_nb[0][17] = U_nb[0][5];
	U_nb[1][17] = U_nb[1][5];
	U_nb[2][17] = U_nb[2][5];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[5])
    {
	U_nb[0][17] = U_nb[0][0];
	U_nb[1][17] = U_nb[1][0];
	U_nb[2][17] = U_nb[2][0];
    }
    else
    {
	U_nb[0][17] = (U_nb[0][0]+U_nb[0][5]+cell_center[index_nb[17]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][17] = (U_nb[1][0]+U_nb[1][5]+cell_center[index_nb[17]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][17] = (U_nb[2][0]+U_nb[2][5]+cell_center[index_nb[17]].m_state.m_U[2]+U_center[2])/4.0;
    }

    for (int l = 0; l < 3; l++)
    {
	dh[l] = top_h[l];

	if (bNoBoundary[l*2])
	    dh0[l] = top_h[l];
	else
	    dh0[l] = top_h[l]/2.0;

	if (bNoBoundary[l*2+1])
	    dh1[l] = top_h[l];
	else
	    dh1[l] = top_h[l]/2.0;
    }


    /***************** Diffusion term for Equation 1  ***************/

    ////////////////   Tensor term 1 //////////////
    //first term

    diffusion[0] += ( 
	             mu_edge[3]*(U_nb[0][3]-U_center[0])/dh1[1] 
	            -mu_edge[2]*(U_center[0]-U_nb[0][2])/dh0[1]
		    ) / dh[1];
    //second term

    diffusion[0] += (mu_edge[2]/r_edge[2]*U_nb[1][6] - mu_edge[2]/r_edge[2]*U_nb[1][7]
	            +mu_edge[3]/r_edge[3]*U_nb[1][8] - mu_edge[3]/r_edge[3]*U_nb[1][9])/(dz*dtheta);

    //////////// Tensor term 2 /////////////
    //first term

    diffusion[0] += ( 
	             mu_edge[5]*(U_nb[0][5]-U_center[0])/dh1[2] 
	            -mu_edge[4]*(U_center[0]-U_nb[0][4])/dh0[2]
		    ) / dh[2];

    //second term

    diffusion[0] += (mu_edge[4]/r_edge[4]*U_face[0][4] - mu_edge[5]/r_edge[5]*U_face[0][5]) / dr;

    //third term

    diffusion[0] += (mu_edge[4]/r_edge[4]*U_nb[2][14] - mu_edge[4]/r_edge[4]*U_nb[2][15]
	            +mu_edge[5]/r_edge[5]*U_nb[2][16] - mu_edge[5]/r_edge[5]*U_nb[2][17])/(dr*dtheta);

    ///////////////  Tensor term 3  /////////////////
    //first term

    diffusion[0] += ( 
	             2.0/r0*mu_edge[1]/r_edge[1]*(U_nb[0][1]-U_center[0])/dh1[0] 
	            -2.0/r0*mu_edge[0]/r_edge[0]*(U_center[0]-U_nb[0][0])/dh0[0]
		    ) / dh[0];

    //second term

    diffusion[0] += (-2.0/r0*mu_edge[0]/r_edge[0]*U_face[2][0] 
	             +2.0/r0*mu_edge[1]/r_edge[1]*U_face[2][1]) / dtheta;

    /////////////////  Tensor term 4     //////////////
    //first term

    diffusion[0] += (-2.0*mu0/r0*U_face[0][4] 
	             +2.0*mu0/r0*U_face[0][5]) / dr;

    //second term

    diffusion[0] += (-2.0*mu0/rr*U_center[0]);

    //third term
    
    diffusion[0] += (-2.0*mu0/rr*U_face[2][0] 
	             +2.0*mu0/rr*U_face[2][1]) / dtheta;

    /******************* Diffusion term for Equation 2 ************/

    //////////////  Tensor term 1  /////////////
    //first term

    diffusion[1] += ( 
	             2.0*mu_edge[3]*(U_nb[1][3]-U_center[1])/dh1[1] 
	            -2.0*mu_edge[2]*(U_center[1]-U_nb[1][2])/dh0[1]
		    ) / dh[1];

    //////////////  Tensor term 2 /////////////////
    //first term

    diffusion[1] += (mu_edge[4]*U_nb[2][10] - mu_edge[4]*U_nb[2][11]
	            +mu_edge[5]*U_nb[2][12] - mu_edge[5]*U_nb[2][13])/(dr*dz);

    //second term

    diffusion[1] += ( 
	             mu_edge[5]*(U_nb[1][5]-U_center[1])/dh1[2] 
	            -mu_edge[4]*(U_center[1]-U_nb[1][4])/dh0[2]
		    ) / dh[2];

    /////////////// Tensor term 3 /////////////////
    //first term

    diffusion[1] += (mu_edge[0]/r0*U_nb[0][6] - mu_edge[1]/r0*U_nb[0][7]
	            +mu_edge[1]/r0*U_nb[0][8] - mu_edge[0]/r0*U_nb[0][9])/(dtheta*dz);

    //second term

    diffusion[1] += ( 
	             mu_edge[1]/r_edge[1]/r0*(U_nb[1][1]-U_center[1])/dh1[0] 
	            -mu_edge[0]/r_edge[0]/r0*(U_center[1]-U_nb[1][0])/dh0[0]
		    ) / dh[0];

    //////////////// Tensor term 4  ///////////////////////
    //first term
 
    diffusion[1] += (-mu0/r0*U_face[2][2] 
	             +mu0/r0*U_face[2][3]) / dz;

    //second term
 
    diffusion[1] += (-mu0/r0*U_face[1][4] 
	             +mu0/r0*U_face[1][5]) / dr;


    /**************** Diffusion term for Equation 3  *****************/

    ///////////// Tensor term 1  ////////////////
    //first term

    diffusion[2] += ( 
	             mu_edge[3]*(U_nb[2][3]-U_center[2])/dh1[1] 
	            -mu_edge[2]*(U_center[2]-U_nb[2][2])/dh0[1]
		    ) / dh[1];

    //second term

    diffusion[2] += (mu_edge[2]*U_nb[1][10] - mu_edge[3]*U_nb[1][11]
	            +mu_edge[3]*U_nb[1][12] - mu_edge[2]*U_nb[1][13])/(dr*dz);

    ////////////// Tensor term 2  ////////////////
    //first term

    diffusion[2] += ( 
	             2.0*mu_edge[5]*(U_nb[2][5]-U_center[2])/dh1[2] 
	            -2.0*mu_edge[4]*(U_center[2]-U_nb[2][4])/dh0[2]
		    ) / dh[2];


    ////////////// Tensor term 3  //////////////////
    //first term

    diffusion[2] += (mu_edge[0]/r0*U_nb[0][14] - mu_edge[1]/r0*U_nb[0][15]
	            +mu_edge[1]/r0*U_nb[0][16] - mu_edge[0]/r0*U_nb[0][17])/(dr*dtheta);

    //second term

    diffusion[2] += (mu_edge[0]/r_edge[0]/r0*U_face[0][0] 
	            -mu_edge[1]/r_edge[1]/r0*U_face[0][1]) / dtheta;

    //third term

    diffusion[2] += ( 
	             mu_edge[1]/r_edge[1]/r0*(U_nb[2][1]-U_center[2])/dh1[0] 
	            -mu_edge[0]/r_edge[0]/r0*(U_center[2]-U_nb[2][0])/dh0[0]
		    ) / dh[0];

    ////////// Tensor term 4 ///////////
    //first term

    diffusion[2] += (-2.0*mu0/r0*U_face[2][4] 
	             +2.0*mu0/r0*U_face[2][5]) / dr;

    ///////////// Tensor term 5  /////////////////
    //first term

    diffusion[2] += (2.0*mu0/rr*U_face[0][0] 
	            -2.0*mu0/rr*U_face[0][1]) / dtheta;

    //second term

    diffusion[2] += (-2.0*mu0/rr*U_center[2]);

}

/**
* TODO: the calculation near the boundary might not be correct.
* @param dir
* @param icoords
* @param dU2
* @param dP
*/
void Incompress_Solver_Smooth_3D_Cylindrical::getDU2(
	int *icoords,
	EBM_COORD xyz,
	double dU2[3],
	double dU[3])
{
    double dh0, dh1, dh;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state;

    bool bNoBoundary[2];

    if(xyz==COORD_X)
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,WEST,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[0];
	else
	    dh0 = top_h[0]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,EAST,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[0];
	else
	    dh1 = top_h[0]/2;

	dh = top_h[0];
    }
    else if(xyz==COORD_Y)	//
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,SOUTH,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[1];
	else
	    dh0 = top_h[1]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,NORTH,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[1];
	else
	    dh1 = top_h[1]/2;

	dh = top_h[1];
    }
    else
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,LOWER,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[2];
	else
	    dh0 = top_h[2]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,UPPER,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[2];
	else
	    dh1 = top_h[2]/2;

	dh = top_h[2];
    }

    // second order derivative
    dU2[0] = ((U2.m_U[0] - U1.m_U[0])/dh1 - (U1.m_U[0] - U0.m_U[0])/dh0) / dh;
    dU2[1] = ((U2.m_U[1] - U1.m_U[1])/dh1 - (U1.m_U[1] - U0.m_U[1])/dh0) / dh;
    dU2[2] = ((U2.m_U[2] - U1.m_U[2])/dh1 - (U1.m_U[2] - U0.m_U[2])/dh0) / dh;

    // first order derivative
    dU[0] = ((U2.m_U[0] - U1.m_U[0])/dh1 + (U1.m_U[0] - U0.m_U[0])/dh0) / 2;
    dU[1] = ((U2.m_U[1] - U1.m_U[1])/dh1 + (U1.m_U[1] - U0.m_U[1])/dh0) / 2;
    dU[2] = ((U2.m_U[2] - U1.m_U[2])/dh1 + (U1.m_U[2] - U0.m_U[2])/dh0) / 2;

}


void Incompress_Solver_Smooth_3D_Cylindrical::getLimitedSlope(
	int *icoords,
	EBM_COORD xyz,
	double slope[3])
{
    double dh0, dh1;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state;

    bool bNoBoundary[2];

    if(xyz==COORD_X)
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,WEST,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[0];
	else
	    dh0 = top_h[0]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,EAST,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[0];
	else
	    dh1 = top_h[0]/2;
    }
    else if(xyz==COORD_Y)	//
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,SOUTH,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[1];
	else
	    dh0 = top_h[1]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,NORTH,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[1];
	else
	    dh1 = top_h[1]/2;
    }
    else
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,LOWER,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[1];
	else
	    dh0 = top_h[1]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,UPPER,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[1];
	else
	    dh1 = top_h[1]/2;
    }
    slope[0] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dh0, (U2.m_U[0]-U1.m_U[0])/dh1);
    slope[1] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dh0, (U2.m_U[1]-U1.m_U[1])/dh1);
    slope[2] = EBM_minmod((U1.m_U[2]-U0.m_U[2])/dh0, (U2.m_U[2]-U1.m_U[2])/dh1);
}


double Incompress_Solver_Smooth_3D_Cylindrical::EBM_minmod(
	double x,
	double y)
{

    double sign = x*y;

    if(sign<0)
	return 0;
    else if(sign>=0)
    {
	if(fabs(x)<fabs(y))
	    return x;
	else
	    return y;
    }
    return 0;
}

/**
* get the state from the neighbor cell or boundary.
* @param icoords
* @param dir
* @param comp
* @param state
* @return true,  valid state from neighbor cell
* 	   false, valid state from boundary
*/
bool Incompress_Solver_Smooth_3D_Cylindrical::getNeighborOrBoundaryState(
	int icoords[3],
	GRID_DIRECTION dir,
	L_STATE &state,
	double t)
{
    double crx_coords[MAXD];
    static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
    POINTER intfc_state;
    HYPER_SURF *hs;

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    int comp = cell_center[index].comp;

    int i,j,k;
    i = icoords[0];
    j = icoords[1];
    k = icoords[2];
    const int  nn = pp_numnodes();
    int        myid = pp_mynode();
    int   *ppgmax = front->pp_grid->gmax;
    int   ppx = myid % ppgmax[0];
    int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
    int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir,
	    comp,&intfc_state,&hs,crx_coords,t) &&
	    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
    {
        state.m_U[0] = getStateVel[0](intfc_state);
        state.m_U[1] = getStateVel[1](intfc_state);
        state.m_U[2] = getStateVel[2](intfc_state);
	return false;
    }
    else
    {
        if ( (ppz == 0 && k == kmin && dir == LOWER)  )
        {
            state.m_U[0] = iFparams->bvel[0][0];
            state.m_U[1] = iFparams->bvel[0][1];
            state.m_U[2] = iFparams->bvel[0][2];
	    return false;
        }
        else if ( (ppz == ppgmax[2]-1 && k == kmax && dir == UPPER) )
        {
            state.m_U[0] = iFparams->bvel[1][0];
            state.m_U[1] = iFparams->bvel[1][1];
            state.m_U[2] = iFparams->bvel[1][2];
            return false;
        }
        else
        {
	    int index_nb;
	    switch(dir)
	    {
	    case WEST:
	        index_nb = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
	        break;
	    case EAST:
	        index_nb = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
	        break;
	    case SOUTH:
	        index_nb = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
	        break;
	    case NORTH:
	        index_nb = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
	        break;
	    case LOWER:
	        index_nb = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
	        break;
	    case UPPER:
	        index_nb = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
	        break;
	    default:
	        assert(false);
	    }
	    state = cell_center[index_nb].m_state;
	    return true;
        }
    }
}


/**
*
* @param state_left
* @param state_right
* @param ans
*/
void Incompress_Solver_Smooth_3D_Cylindrical::getRiemannSolution(
	EBM_COORD xyz,
	L_STATE &state_left,
	L_STATE &state_right,
	L_STATE &ans)
{
    L_STATE sl, sr;


    // rotate state
    if(xyz==COORD_X)
    {
	sl.m_U[0] = state_left.m_U[0];
	sl.m_U[1] = state_left.m_U[1];
	sl.m_U[2] = state_left.m_U[2];

	sr.m_U[0] = state_right.m_U[0];
	sr.m_U[1] = state_right.m_U[1];
	sr.m_U[2] = state_right.m_U[2];
    }
    else if(xyz==COORD_Y)
    {
	sl.m_U[0] = state_left.m_U[1];
	sl.m_U[1] = state_left.m_U[0];
	sl.m_U[2] = state_left.m_U[2];

	sr.m_U[0] = state_right.m_U[1];
	sr.m_U[1] = state_right.m_U[0];
	sr.m_U[2] = state_right.m_U[2];
    }
    else
    {
	sl.m_U[0] = state_left.m_U[2];
	sl.m_U[1] = state_left.m_U[1];
	sl.m_U[2] = state_left.m_U[0];

	sr.m_U[0] = state_right.m_U[2];
	sr.m_U[1] = state_right.m_U[1];
	sr.m_U[2] = state_right.m_U[0];
    }

    // calculate the Riemann solution
    double uL = sl.m_U[0];
    double uR = sr.m_U[0];

    // BCG, JCP 85, 257-283 (1989)
    // ut + uux = 0
    if(uL>=0 && (uL+uR)>=0)
	ans.m_U[0] = uL;
    else if(uL<0 && uR>0)
	ans.m_U[0] = 0;
    else
	ans.m_U[0] = uR;

    // vt + uvx = 0
    if(ans.m_U[0]>0)
	ans.m_U[1] = sl.m_U[1];
    else if(ans.m_U[0]<0)
	ans.m_U[1] = sr.m_U[1];
    else
	ans.m_U[1] = 1.0/2*(sl.m_U[1]+sr.m_U[1]);

    if(ans.m_U[0]>0)
	ans.m_U[2] = sl.m_U[2];
    else if(ans.m_U[0]<0)
	ans.m_U[2] = sr.m_U[2];
    else
	ans.m_U[2] = 1.0/2*(sl.m_U[2]+sr.m_U[2]);

    // rotate state
    if(xyz==COORD_X)
	; // do nothing
    else if(xyz==COORD_Y)
	std::swap(ans.m_U[0],ans.m_U[1]);
    else
	std::swap(ans.m_U[0],ans.m_U[2]);
}

void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
        double coords[MAXD],crx_coords[MAXD];
	double coeff[18],mu0,rho,rhs,U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	int num_iter;
	double rel_residual;
	double r, rr, redge[2];

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 9, 9);
	// 7theta + 2r  for the first equation
	// 7z for the second equation
	// 7r + 2theta for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

        const int  nn = pp_numnodes();
        int        myid = pp_mynode();
        int   *ppgmax = front->pp_grid->gmax;
        int   ppx = myid % ppgmax[0];
        int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
        int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

	printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
	

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];
	    U2_center = cell_center[index].m_state.m_U[2];

	  
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                {
                    U0_nb[nb] = getStateVel[0](intfc_state);
                    U1_nb[nb] = getStateVel[1](intfc_state);
                    U2_nb[nb] = getStateVel[2](intfc_state);
                }
                else
                {
                if ( (ppz == 0 && k == kmin && dir[nb] == LOWER) )
	       	{
                    U0_nb[nb] = iFparams->bvel[0][0];
                    U1_nb[nb] = iFparams->bvel[0][1];
                    U2_nb[nb] = iFparams->bvel[0][2];
		}
                else if ( (ppz == ppgmax[2]-1 && k == kmax && dir[nb] == UPPER) )
                {
                    U0_nb[nb] = iFparams->bvel[1][0];
                    U1_nb[nb] = iFparams->bvel[1][1];
                    U2_nb[nb] = iFparams->bvel[1][2];
                }
                else
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
                }
	    }

            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


    //Setting the coeffecients for the first equation

	    rr = cell_center[index].m_coords[2] * cell_center[index].m_coords[2];
	    r = cell_center[index].m_coords[2];
	    redge[0] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[4]].m_coords[2]);
	    redge[1] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[5]].m_coords[2]);


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++) //For the second order boundary
	    {
		if (I_nb[nb] == -1)
		{
		    coeff[nb] = 2.0*coeff[nb];
		}
	    }

	    solver.Set_A(I*3,I*3,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+0.5*m_dt*mu0/(rho*rr));
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5]-0.5*m_dt*mu0/(rho*rr))*cell_center[index].m_state.m_U[0];


	    coeff[6] = -0.5*m_dt/rho * mu0/(rr*top_h[0]);
	    coeff[7] =  0.5*m_dt/rho * mu0/(rr*top_h[0]);


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }

	    if (I_nb[0] != -1)
	    {
		solver.Set_A(I*3,I_nb[0]*3+2,-coeff[6]);
		rhs += coeff[6]*U2_nb[0];
	    }
	    else
	    {
		solver.Set_A(I*3,I*3+2,coeff[6]);
		rhs += -coeff[6]*U2_center;
		coeff[6] = coeff[6] * 2.0;
		rhs += 2.0*coeff[6]*U2_nb[0];
	    }

	    if (I_nb[1] != -1)
	    {
		solver.Set_A(I*3,I_nb[1]*3+2,-coeff[7]);
		rhs += coeff[7]*U2_nb[1];
	    }
	    else
	    {
		solver.Set_A(I*3,I*3+2,coeff[7]);
		rhs += -coeff[7]*U2_center;
		coeff[7] = coeff[7] * 2.0;
		rhs += 2.0*coeff[7]*U2_nb[1];
	    }

	    //rhs -= m_dt*cell_center[index].m_state.m_U[0]*cell_center[index].m_state.m_U[2]/r; //Source term in advection step
	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[0];



	    solver.Set_b(I*3, rhs);

	    /************************************************************************/


    //Setting the coeffecients for the second equation


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+1,I*3+1,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[1];


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[1];


	    solver.Set_b(I*3+1, rhs);

	    /************************************************************************/

    //Setting the coeffecients for the third equation

	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+2,I*3+2,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+0.5*m_dt*mu0/(rho*rr));
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5]-0.5*m_dt*mu0/(rho*rr))*cell_center[index].m_state.m_U[2];


	    coeff[6] =  0.5*m_dt/rho * mu0/(rr*top_h[0]);
	    coeff[7] = -0.5*m_dt/rho * mu0/(rr*top_h[0]);


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }


	    if (I_nb[0] != -1)
	    {
		solver.Set_A(I*3+2, I_nb[0]*3, -coeff[6]);
		rhs += coeff[6]*U0_nb[0];
	    }
	    else
	    {
		solver.Set_A(I*3+2, I*3, coeff[6]);
		rhs += -coeff[6]*U0_center;
		coeff[6] = coeff[6] * 2.0;
		rhs += 2.0*coeff[6]*U0_nb[0];
	    }
	    
	    if (I_nb[1] != -1)
	    {
		solver.Set_A(I*3+2, I_nb[1]*3, -coeff[7]);
		rhs += coeff[7]*U0_nb[1];
	    }
	    else
	    {
		solver.Set_A(I*3+2, I*3, coeff[7]);
		rhs += -coeff[7]*U0_center;
		coeff[7] = coeff[7] * 2.0;
		rhs += 2.0*coeff[7]*U0_nb[1];
	    }

	    //rhs += m_dt* cell_center[index].m_state.m_U[0] * cell_center[index].m_state.m_U[0] / r; //Source term in advection step
	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[2];

	    solver.Set_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
/*
	if (rel_residual > 1)
	{
	    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	    solver.Reset_x();
	    solver.Solve_GMRES();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	}
*/
	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }

        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);
}       /* end 3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled */

void Incompress_Solver_Smooth_3D_Cylindrical::compDiff_CellFace(
	PETSc *pSolver,
	int I,
	int I_nb[18],
	double U_center[3],
	double U_nb[3][18],
	int flag[6],
	int equation_index,
	int vel_comp,
	int face_index,
	double coeff)
{
    double rhs;

    if (flag[face_index] == 1)
	rhs = coeff*U_nb[vel_comp][face_index];
    else 
    {
	rhs = coeff*(U_nb[vel_comp][face_index]+U_center[vel_comp])/4.0;

	pSolver->Add_A(I*3+equation_index,I*3+vel_comp,                -coeff/4.0);
	pSolver->Add_A(I*3+equation_index,I_nb[face_index]*3+vel_comp, -coeff/4.0);
    }
    pSolver->Add_b(I*3+equation_index,rhs);

}

void Incompress_Solver_Smooth_3D_Cylindrical::compDiff_CellCorner(
	PETSc *pSolver,
	int I, 
	int I_nb[18],
	double U_center[3],
	double U_nb[3][18], 
	int flag[6],
	int equation_index,
	int vel_comp, 
	int corner_index, 
	double coeff)
{
    double rhs;

    if (corner_index == 6)
    {
	if      (flag[0] == 1 && flag[2] == 0)
	    rhs = coeff*U_nb[vel_comp][0];
	else if (flag[0] == 0 && flag[2] == 1)
	    rhs = coeff*U_nb[vel_comp][2];
	else if (flag[0] == 1 && flag[2] == 1)
	    rhs = coeff*U_nb[vel_comp][0];
	else {
	    rhs = coeff*(U_nb[vel_comp][0]
		        +U_nb[vel_comp][2]
			+U_nb[vel_comp][6]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[0]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[2]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[6]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 7)
    {
	if      (flag[1] == 1 && flag[2] == 0)
	    rhs = coeff*U_nb[vel_comp][1];
	else if (flag[1] == 0 && flag[2] == 1)
	    rhs = coeff*U_nb[vel_comp][2];
	else if (flag[1] == 1 && flag[2] == 1)
	    rhs = coeff*U_nb[vel_comp][1];
	else {
	    rhs = coeff*(U_nb[vel_comp][1]
		        +U_nb[vel_comp][2]
			+U_nb[vel_comp][7]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[1]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[2]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[7]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 8)
    {
	if      (flag[1] == 1 && flag[3] == 0)
	    rhs = coeff*U_nb[vel_comp][1];
	else if (flag[1] == 0 && flag[3] == 1)
	    rhs = coeff*U_nb[vel_comp][3];
	else if (flag[1] == 1 && flag[3] == 1)
	    rhs = coeff*U_nb[vel_comp][1];
	else {
	    rhs = coeff*(U_nb[vel_comp][1]
		        +U_nb[vel_comp][3]
			+U_nb[vel_comp][8]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[1]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[3]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[8]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 9)
    {
	if      (flag[0] == 1 && flag[3] == 0)
	    rhs = coeff*U_nb[vel_comp][0];
	else if (flag[0] == 0 && flag[3] == 1)
	    rhs = coeff*U_nb[vel_comp][3];
	else if (flag[0] == 1 && flag[3] == 1)
	    rhs = coeff*U_nb[vel_comp][0];
	else {
	    rhs = coeff*(U_nb[vel_comp][0]
		        +U_nb[vel_comp][3]
			+U_nb[vel_comp][9]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[0]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[3]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[9]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 10)
    {
	if      (flag[2] == 1 && flag[4] == 0)
	    rhs = coeff*U_nb[vel_comp][2];
	else if (flag[2] == 0 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][4];
	else if (flag[2] == 1 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][2];
	else {
	    rhs = coeff*(U_nb[vel_comp][2]
		        +U_nb[vel_comp][4]
			+U_nb[vel_comp][10]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[2]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[4]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[10]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 11)
    {
	if      (flag[3] == 1 && flag[4] == 0)
	    rhs = coeff*U_nb[vel_comp][3];
	else if (flag[3] == 0 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][4];
	else if (flag[3] == 1 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][3];
	else {
	    rhs = coeff*(U_nb[vel_comp][3]
		        +U_nb[vel_comp][4]
			+U_nb[vel_comp][11]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[3]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[4]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[11]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 12)
    {
	if      (flag[3] == 1 && flag[5] == 0)
	    rhs = coeff*U_nb[vel_comp][3];
	else if (flag[3] == 0 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][5];
	else if (flag[3] == 1 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][3];
	else {
	    rhs = coeff*(U_nb[vel_comp][3]
		        +U_nb[vel_comp][5]
			+U_nb[vel_comp][12]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[3]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[5]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[12]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 13)
    {
	if      (flag[2] == 1 && flag[5] == 0)
	    rhs = coeff*U_nb[vel_comp][2];
	else if (flag[2] == 0 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][5];
	else if (flag[2] == 1 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][2];
	else {
	    rhs = coeff*(U_nb[vel_comp][2]
		        +U_nb[vel_comp][5]
			+U_nb[vel_comp][13]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[2]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[5]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[13]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 14)
    {
	if      (flag[0] == 1 && flag[4] == 0)
	    rhs = coeff*U_nb[vel_comp][0];
	else if (flag[0] == 0 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][4];
	else if (flag[0] == 1 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][0];
	else {
	    rhs = coeff*(U_nb[vel_comp][0]
		        +U_nb[vel_comp][4]
			+U_nb[vel_comp][14]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[0]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[4]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[14]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 15)
    {
	if      (flag[1] == 1 && flag[4] == 0)
	    rhs = coeff*U_nb[vel_comp][1];
	else if (flag[1] == 0 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][4];
	else if (flag[1] == 1 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][1];
	else {
	    rhs = coeff*(U_nb[vel_comp][1]
		        +U_nb[vel_comp][4]
			+U_nb[vel_comp][15]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[1]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[4]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[15]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 16)
    {
	if      (flag[1] == 1 && flag[5] == 0)
	    rhs = coeff*U_nb[vel_comp][1];
	else if (flag[1] == 0 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][5];
	else if (flag[1] == 1 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][1];
	else {
	    rhs = coeff*(U_nb[vel_comp][1]
		        +U_nb[vel_comp][5]
			+U_nb[vel_comp][16]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[1]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[5]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[16]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 17)
    {
	if      (flag[0] == 1 && flag[5] == 0)
	    rhs = coeff*U_nb[vel_comp][0];
	else if (flag[0] == 0 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][5];
	else if (flag[0] == 1 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][0];
	else {
	    rhs = coeff*(U_nb[vel_comp][0]
		        +U_nb[vel_comp][5]
			+U_nb[vel_comp][17]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[0]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[5]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[17]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
}


void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_coupled(void)
{
        COMPONENT comp;
        int index,index_nb[18],size;
        int I,I_nb[18];
        double coords[MAXD],crx_coords[MAXD];
	double coeff_temp, coeff0, coeff1;
	double mu0,rho,rhs,mu[6],mu_edge[6];
	int flag[6];
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;
	double r0, rr, r_edge[6],r[6];
	double dr = top_h[2];
	double dz = top_h[1];
	double dtheta = top_h[0];

	double U_center[3];
	double U_nb[3][18];

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 30, 30);
	// 7theta + 9z + 9r  for the first equation
	// 7z + 9theta + 9r  for the second equation
	// 7r + 9theta + 9z  for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

        const int  nn = pp_numnodes();
        int        myid = pp_mynode();
        int   *ppgmax = front->pp_grid->gmax;
        int   ppx = myid % ppgmax[0];
        int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
        int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

	printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

	start_clock("SetMatrixForDiffusion");
	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	//theta-z cut neighbours
	    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
	    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
	    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
	    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	
	//z-r cut neighbours
	    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
	    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
	    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
	    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

	//theta-r cut neighbours
	    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
	    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
	    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
	    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);



        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
		
	//theta-z cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
	    I_nb[7] = ijk_to_I[i+1][j-1][k];
	    I_nb[8] = ijk_to_I[i+1][j+1][k];
	    I_nb[9] = ijk_to_I[i-1][j+1][k];
	//z-r cut neighbours
	    I_nb[10] = ijk_to_I[i][j-1][k-1];
	    I_nb[11] = ijk_to_I[i][j+1][k-1];
	    I_nb[12] = ijk_to_I[i][j+1][k+1];
	    I_nb[13] = ijk_to_I[i][j-1][k+1];
	//theta-r cut neighbours
	    I_nb[14] = ijk_to_I[i-1][j][k-1];
	    I_nb[15] = ijk_to_I[i+1][j][k-1];
	    I_nb[16] = ijk_to_I[i+1][j][k+1];
	    I_nb[17] = ijk_to_I[i-1][j][k+1];



	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U_center[0] = cell_center[index].m_state.m_U[0];
	    U_center[1] = cell_center[index].m_state.m_U[1];
	    U_center[2] = cell_center[index].m_state.m_U[2];


	    r0 = cell_center[index].m_coords[2];
	    rr = r0*r0;

	  
            for (nb = 0; nb < 6; nb++) //dealing with 6 neighbours
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                {
                    flag[nb] = 1;
                    U_nb[0][nb] = getStateVel[0](intfc_state);
                    U_nb[1][nb] = getStateVel[1](intfc_state);
                    U_nb[2][nb] = getStateVel[2](intfc_state);
                    mu[nb] = mu0;
                    mu_edge[nb] = mu0;

                    //printf("\n FrontWholeStateAtGridCrossing_tmp_U i = %d, j = %d, k = %d, nb = %d U_theta = %.16g, U_z = %.16g, U_r = %.16g dir %d\n", i,j,k,nb,U0_nb[nb], U1_nb[nb], U2_nb[nb],dir[nb]);
                }
		else
		{
		    if ( (ppz == 0 && k == kmin && dir[nb] == LOWER) )
		    {
			flag[nb] = 1;
                        U_nb[0][nb] = iFparams->bvel[0][0];
                        U_nb[1][nb] = iFparams->bvel[0][1];
                        U_nb[2][nb] = iFparams->bvel[0][2];
		        mu[nb] = mu0;
		        mu_edge[nb] = mu0;
		    }
                    else if ( (ppz == ppgmax[2]-1 && k == kmax && dir[nb] == UPPER) )
                    {
                        flag[nb] = 1;
                        U_nb[0][nb] = iFparams->bvel[1][0];
                        U_nb[1][nb] = iFparams->bvel[1][1];
                        U_nb[2][nb] = iFparams->bvel[1][2];
                        mu[nb] = mu0;
                        mu_edge[nb] = mu0;
                    }
		    else
		    {
		        flag[nb] = 0;
		        U_nb[0][nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		        U_nb[1][nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		        U_nb[2][nb] = cell_center[index_nb[nb]].m_state.m_U[2];

		        mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		        mu_edge[nb] = 0.5*(mu0 + mu[nb]);
		    }
		}

		r[nb] = cell_center[index_nb[nb]].m_coords[2];
		r_edge[nb] = 0.5 * (r[nb] + r0);
	    }

	    for (nb = 6; nb < 18; nb++) //corner values for interior
	    {
		if (I_nb[nb] != -1)
		{
		    U_nb[0][nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U_nb[1][nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U_nb[2][nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }

	    //source term
            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


	    /****************** The first equation about theta *************/

	    solver.Add_A(I*3,I*3,1.0);
	    rhs = U_center[0];

	    /////////////////    term d(tao_z0)/dz      //////////////////
	    //first term in tao_z0: mu*du_0/dz

	    if (flag[2] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][2];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[2]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][2];
		rhs -= coeff_temp*U_center[0];
	    }

	    if (flag[3] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][3];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho *mu_edge[3]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[3]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][3];
		rhs -= coeff_temp*U_center[0];
	    }

	    //second term mu/r * du_z/d0
	    // 4 corners
	    //
	    coeff0 = m_dt/rho * mu_edge[2]/r_edge[2] / (dtheta*dz);
	    coeff1 = m_dt/rho * mu_edge[3]/r_edge[3] / (dtheta*dz);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,1,6,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,1,7, -coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,1,8,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,1,9, -coeff1);

	    ////////////      term d(tao_r0)/dr       /////////////////
	    //first term in tao_r0: mu*du_0/dr

	    if (flag[4] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][4];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[4]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][4];
		rhs -= coeff_temp*U_center[0];
	    }

	    if (flag[5] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[5]/(dr*dr);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][5];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho *mu_edge[5]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[5]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][5];
		rhs -= coeff_temp*U_center[0];
	    }

	    //second term -u_0/r
	    //
	    coeff0 = m_dt/rho * mu_edge[4]/r_edge[4]/dr;
	    coeff1 = m_dt/rho * mu_edge[5]/r_edge[5]/dr;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,0,4, coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,0,5,-coeff1);

	    //third term 
	    coeff0 = m_dt/rho * mu_edge[4]/r_edge[4]/(dtheta*dr);
	    coeff1 = m_dt/rho * mu_edge[5]/r_edge[5]/(dtheta*dr);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,14,   coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,15,  -coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,16,   coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,17,  -coeff1);

	    /////////////////   term 1/r * dtao_00/d0   //////////////
	    //first term 2*mu/r*du_0/d0
	    
	    if (flag[0] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][0];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[0]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][0];
		rhs -= coeff_temp*U_center[0];
	    }

	    if (flag[1] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][1];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[1]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][1];
		rhs -= coeff_temp*U_center[0];
	    }

	    //second term
	    coeff0 = 2.0*m_dt/rho/r0 * mu_edge[0]/r_edge[0] / dtheta;
	    coeff1 = 2.0*m_dt/rho/r0 * mu_edge[1]/r_edge[1] / dtheta;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,0,  -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,1,   coeff1);

	    /////////////////////   term 2*tao_r0/r  //////////////
	    //first term
	    coeff0 = 2.0*m_dt/rho/r0 * mu0 / dr;
	    coeff1 = 2.0*m_dt/rho/r0 * mu0 / dr;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,0,4, -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,0,5,  coeff1);

	    //second term
	    coeff_temp = -m_dt/rho * mu0/rr;
	    rhs += coeff_temp*U_center[0];
	    solver.Add_A(I*3,I*3, -coeff_temp);

	    //third term
	    coeff0 = 2.0*m_dt/rho/rr * mu0 / dtheta;
	    coeff1 = 2.0*m_dt/rho/rr * mu0 / dtheta;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,0, -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,1,  coeff1);

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[0];

	    solver.Add_b(I*3, rhs);

	    /***************** Second Equation ***************/

	    solver.Add_A(I*3+1,I*3+1,1.0);
	    rhs = U_center[1];

	    ////////////////     term d(tao_zz)/dz   ///////////
	    //first term

	    if (flag[2] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][2];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[2]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][2];
		rhs -= coeff_temp*U_center[1];
	    }

	    if (flag[3] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][3];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[3]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][3];
		rhs -= coeff_temp*U_center[1];
	    }

	    //////////////   term d(tar_rz)/dr     /////////
	    //first term

	    coeff0 = m_dt/rho * mu_edge[4] / (dz*dr);
	    coeff1 = m_dt/rho * mu_edge[5] / (dz*dr);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,10,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,11, -coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,12,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,13, -coeff1);

	    //second term

	    if (flag[4] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][4];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[4]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][4];
		rhs -= coeff_temp*U_center[1];
	    }
	    
	    if (flag[5] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[5]/(dr*dr);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][5];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[5]/(dr*dr);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[5]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][5];
		rhs -= coeff_temp*U_center[1];
	    }

	    ///////////////    term  1/r * d tao_0z/d0    //////////////
	    //first term

	    coeff0 = m_dt/rho/r0 * mu_edge[0] / (dtheta*dz);
	    coeff1 = m_dt/rho/r0 * mu_edge[1] / (dtheta*dz);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,0,6,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,0,7, -coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,0,8,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,0,9, -coeff0);

	    //second term

	    if (flag[0] == 1)
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][0];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[0]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][0];
		rhs -= coeff_temp*U_center[1];
	    }

	    if (flag[1] == 1)
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][1];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[1]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][1];
		rhs -= coeff_temp*U_center[1];
	    }

	    //////////    term  tao_rz/r    ////////////
	    //first term

	    coeff0 = m_dt/rho/r0 * mu0 / dz;
	    coeff1 = m_dt/rho/r0 * mu0 / dz;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,2, -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,3,  coeff1);

	    //second term

	    coeff0 = m_dt/rho/r0 * mu0 / dr;
	    coeff1 = m_dt/rho/r0 * mu0 / dr;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    1,1,4, -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    1,1,5,  coeff1);

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[1];

	    solver.Add_b(I*3+1, rhs);

	    /****************** Third Equation *****************/

	    solver.Add_A(I*3+2,I*3+2,1.0);
	    rhs = U_center[2];

	    /////////////    term d tao_zr / dz    //////////////////
	    //first term

	    if (flag[2] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][2];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[2]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][2];
		rhs -= coeff_temp*U_center[2];
	    }

	    if (flag[3] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][3];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[3]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][3];
		rhs -= coeff_temp*U_center[2];
	    }

	    //second term

	    coeff0 = m_dt/rho * mu_edge[2] / (dz*dr);
	    coeff1 = m_dt/rho * mu_edge[3] / (dz*dr);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,1,10,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,1,11, -coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,1,12,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,1,13, -coeff0);


	    ///////////////    term d tao_rr / dr   /////////////

	    if (flag[4] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][4];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[4]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][4];
		rhs -= coeff_temp*U_center[2];
	    }

	    if (flag[5] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho * mu_edge[5]/(dr*dr);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][5];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = m_dt/rho * mu_edge[5]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[5]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][5];
		rhs -= coeff_temp*U_center[2];
	    }

	    ///////////     term 1/r * d tao_0r / d0   ///////////
	    //first term

	    coeff0 = m_dt/rho/r0 * mu_edge[0] / (dtheta*dr);
	    coeff1 = m_dt/rho/r0 * mu_edge[1] / (dtheta*dr);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,14,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,15, -coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,16,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,17, -coeff0);

	    //second term

	    coeff0 = m_dt/rho/r0 * mu_edge[0] / r_edge[0] / dtheta;
	    coeff1 = m_dt/rho/r0 * mu_edge[1] / r_edge[1] / dtheta;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,0,  coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,1, -coeff1);

	    //third term
	    
	    if (flag[0] == 1)
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][0];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[0]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][0];
		rhs -= coeff_temp*U_center[2];
	    }

	    if (flag[1] == 1)
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][1];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[1]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][1];
		rhs -= coeff_temp*U_center[2];
	    }

	    /////////////////    term  tao_rr / r    /////////////////////

	    coeff0 = 2.0*m_dt/rho * mu0/r0 / dr;
	    coeff1 = 2.0*m_dt/rho * mu0/r0 / dr;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,2,4,  -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,2,5,   coeff1);


	    //////////////////   term  -tao_00 / r    ///////////////////
	    //first term

	    coeff0 = 2.0*m_dt/rho * mu0/rr / dtheta;
	    coeff1 = 2.0*m_dt/rho * mu0/rr / dtheta;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,0,   coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,1,  -coeff1);

	    //second term

	    coeff_temp = -m_dt/rho * mu0/rr;
	    rhs += coeff_temp*U_center[2];
	    solver.Add_A(I*3+2,I*3+2, -coeff_temp);


	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[2];

	    solver.Add_b(I*3+2, rhs);
        }

	stop_clock("SetMatrixForDiffusion");

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

	start_clock("Petsc_Solve_DIFF");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
/*
	if (rel_residual > 1)
	{
	    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	    solver.Reset_x();
	    solver.Solve_GMRES();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	}
*/
	stop_clock("Petsc_Solve_DIFF");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_coupled: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);


	start_clock("GetSolAndScat");
	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }

        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);

	stop_clock("GetSolAndScat");
}       /* end 3D_Cylindrical::compDiffWithSmoothProperty_2nd_coupled */


//-------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void Incompress_Solver_Smooth_3D_Cylindrical::computeAdvection_test(void)
{
	int i,j,k,l;
	int index,index00,index01,index10,index11,index20,index21,size;
	L_STATE state;
	COMPONENT comp;
	double speed;
	double *u, *v, *w;
	double u0,u00,u01,u10,u11,u20,u21;
	double v0,v00,v01,v10,v11,v20,v21;
	double w0,w00,w01,w10,w11,w20,w21;
	double crx_coords[MAXD];
	int icoords[MAXD];
	POINTER intfc_state;
	HYPER_SURF *hs;

	double rho, mu0;

	max_speed = 0.0;

	size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    u[index] = cell_center[index].m_state.m_U[0];
	    v[index] = cell_center[index].m_state.m_U[1];
	    w[index] = cell_center[index].m_state.m_U[2];
	}

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (comp == SOLID_COMP)
	    {
	    	cell_center[index].m_state.m_U[0] = 0.0; 
	    	cell_center[index].m_state.m_U[1] = 0.0; 
	    	cell_center[index].m_state.m_U[2] = 0.0; 
		continue;
	    }
	    u0 = u[index];
	    v0 = v[index];
	    w0 = w[index];

	    rho = cell_center[index].m_state.m_rho;
	    mu0 = cell_center[index].m_state.m_mu;


	    // To continue
	    index00 = d_index3d(i-1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u00 = getStateXvel(intfc_state);
		v00 = getStateYvel(intfc_state);
		w00 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u00 = u[index00];
		v00 = v[index00];
		w00 = w[index00];
	    }
	    index01 = d_index3d(i+1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,EAST,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u01 = getStateXvel(intfc_state);
		v01 = getStateYvel(intfc_state);
		w01 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u01 = u[index01];
		v01 = v[index01];
		w01 = w[index01];
	    }
	    index10 = d_index3d(i,j-1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,SOUTH,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u10 = getStateXvel(intfc_state);
		v10 = getStateYvel(intfc_state);
		w10 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u10 = u[index10];
		v10 = v[index10];
		w10 = w[index10];
	    }
	    index11 = d_index3d(i,j+1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,NORTH,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u11 = getStateXvel(intfc_state);
		v11 = getStateYvel(intfc_state);
		w11 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u11 = u[index11];
		v11 = v[index11];
		w11 = w[index11];
	    }
	    index20 = d_index3d(i,j,k-1,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u20 = getStateXvel(intfc_state);
		v20 = getStateYvel(intfc_state);
		w20 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u20 = u[index20];
		v20 = v[index20];
		w20 = w[index20];
	    }
	    index21 = d_index3d(i,j,k+1,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u21 = getStateXvel(intfc_state);
		v21 = getStateYvel(intfc_state);
		w21 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u21 = u[index21];
		v21 = v[index21];
		w21 = w[index21];
	    }


	    /*
	    cell_center[index].m_state.m_U[0] += -m_dt*(
				burger_flux(u00,u0,u01)/top_h[0] +
				linear_flux(v0,u10,u0,u11)/top_h[1] +
				linear_flux(w0,u20,u0,u21)/top_h[2]);

	    cell_center[index].m_state.m_U[1] += - m_dt*(
				linear_flux(u0,v00,v0,v01)/top_h[0] +
			 	burger_flux(v10,v0,v11)/top_h[1] +
				linear_flux(w0,v20,v0,v21)/top_h[2]);

	    cell_center[index].m_state.m_U[2] += - m_dt*(
				linear_flux(u0,w00,w0,w01)/top_h[0] +
				linear_flux(v0,w10,w0,w11)/top_h[1] +
			 	burger_flux(w20,w0,w21)/top_h[2]);
	    */

	    double r = cell_center[index].m_coords[2];

            cell_center[index].m_state.m_U[0] += - m_dt*(
                                burger_flux(u00,u0,u01)/(r*top_h[0]) +
                                linear_flux(v0,u10,u0,u11)/top_h[1] +
                                linear_flux(w0,u20,u0,u21)/top_h[2]);

            cell_center[index].m_state.m_U[1] += - m_dt*(
                                linear_flux(u0,v00,v0,v01)/(r*top_h[0]) +
                                burger_flux(v10,v0,v11)/top_h[1] +
                                linear_flux(w0,v20,v0,v21)/top_h[2]);

            cell_center[index].m_state.m_U[2] += - m_dt*(
                                linear_flux(u0,w00,w0,w01)/(r*top_h[0]) +
                                linear_flux(v0,w10,w0,w11)/top_h[1] +
                                burger_flux(w20,w0,w21)/top_h[2]);
          
            //SOURCE TERM
	    //
	    //How about calculating it in the diffusion solver?
	    
            cell_center[index].m_state.m_U[0] += -(m_dt*u0*w0/r);

            cell_center[index].m_state.m_U[2] += (m_dt*u0*u0/r);

	    // First derivative term in the diffusion solver

	    cell_center[index].m_state.m_U[0] += mu0/rho * 2.0/(r*r) * (w01 - w00)/ (2.0*top_h[0]);
	    cell_center[index].m_state.m_U[2] -= mu0/rho * 2.0/(r*r) * (u01 - u00)/ (2.0*top_h[0]);


	    
	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]) +
		    fabs(cell_center[index].m_state.m_U[2]);
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.m_U[l] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	FT_FreeThese(3,u,v,w);
/*
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    //printf("After advection, u_theta = %.16g, u_z = %.16g, u_r = %.16g\n",cell_center[index].m_state.m_U[0], cell_center[index].m_state.m_U[1], cell_center[index].m_state.m_U[2]);
	}
	*/
}	/* end computeAdvection3d */



void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_1st_decoupled_test(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
        double coords[MAXD],crx_coords[MAXD];
	double coeff[18],mu0,rho,rhs,U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	int num_iter;
	double rel_residual;
	double r, rr, redge[2],rnb[2];

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 9, 9);
	// 7theta + 2r  for the first equation
	// 7z for the second equation
	// 7r + 2theta for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
	

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];
	    U2_center = cell_center[index].m_state.m_U[2];

	  
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);
		}
                else
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }

            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


    //Setting the coeffecients for the first equation

	    rr = cell_center[index].m_coords[2] * cell_center[index].m_coords[2];
	    r = cell_center[index].m_coords[2];
	    redge[0] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[4]].m_coords[2]);
	    redge[1] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[5]].m_coords[2]);
	    rnb[0] = cell_center[index_nb[4]].m_coords[2];
	    rnb[1] = cell_center[index_nb[5]].m_coords[2];


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    //coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    //coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    coeff[4] = 0.5*m_dt/rho * mu0 / (redge[0]*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0 / (redge[1]*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++) //For the second order boundary
	    {
		if (I_nb[nb] == -1)
		{
		    //printf("\n i = %d, j = %d, k = %d, neighbour = %d\n",i,j,k,nb);
		    coeff[nb] = 2.0*coeff[nb];
		}
	    }

	    solver.Set_A(I*3,I*3,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]*r+coeff[5]*r);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]*r-coeff[5]*r)*cell_center[index].m_state.m_U[0];

	    if (I_nb[4] == -1)
		coeff[4] = coeff[4] * redge[0];
	    else
		coeff[4] = coeff[4] * rnb[0];

	    if (I_nb[5] == -1)
		coeff[5] = coeff[5] * redge[1];
	    else
		coeff[5] = coeff[5] * rnb[1];

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }


	    //rhs -= m_dt*cell_center[index].m_state.m_U[0]*cell_center[index].m_state.m_U[2]/r; //Source term in advection step
	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;

	    solver.Set_b(I*3, rhs);

	    /************************************************************************/


    //Setting the coeffecients for the second equation


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+1,I*3+1,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[1];


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;

	    solver.Set_b(I*3+1, rhs);

	    /************************************************************************/

    //Setting the coeffecients for the third equation

	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0/(redge[0]*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0/(redge[1]*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+2,I*3+2,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]*r+coeff[5]*r);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]*r-coeff[5]*r)*cell_center[index].m_state.m_U[2];

	    if (I_nb[4] == -1)
		coeff[4] = coeff[4] * redge[0];
	    else
		coeff[4] = coeff[4] * rnb[0];

	    if (I_nb[5] == -1)
		coeff[5] = coeff[5] * redge[1];
	    else
		coeff[5] = coeff[5] * rnb[1];

	    
	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }

	    //rhs += m_dt* cell_center[index].m_state.m_U[0] * cell_center[index].m_state.m_U[0] / r; //Source term in advection step
	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;

	    solver.Set_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
/*
	if (rel_residual > 1)
	{
	    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	    solver.Reset_x();
	    solver.Solve_GMRES();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	}
*/
	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_1st_decoupled_test: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }
        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);
	/*
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    //printf("After diffusion, u_theta = %.16g, u_z = %.16g, u_r = %.16g\n",cell_center[index].m_state.m_U[0], cell_center[index].m_state.m_U[1], cell_center[index].m_state.m_U[2]);
	}
*/
}       /* end compDiffWithSmoothProperty3d */


void Incompress_Solver_Smooth_3D_Cylindrical::computeAdvection(void)
{
	int i,j,k,l;
	int index,index00,index01,index10,index11,index20,index21,size;
	L_STATE state;
	COMPONENT comp;
	double speed;
	double *u, *v, *w;
	double u0,u00,u01,u10,u11,u20,u21;
	double v0,v00,v01,v10,v11,v20,v21;
	double w0,w00,w01,w10,w11,w20,w21;
	double crx_coords[MAXD];
	int icoords[MAXD];
	POINTER intfc_state;
	HYPER_SURF *hs;

	max_speed = 0.0;

	size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    u[index] = cell_center[index].m_state.m_U[0];
	    v[index] = cell_center[index].m_state.m_U[1];
	    w[index] = cell_center[index].m_state.m_U[2];
	}

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (comp == SOLID_COMP)
	    {
	    	cell_center[index].m_state.m_U[0] = 0.0; 
	    	cell_center[index].m_state.m_U[1] = 0.0; 
	    	cell_center[index].m_state.m_U[2] = 0.0; 
		continue;
	    }
	    u0 = u[index];
	    v0 = v[index];
	    w0 = w[index];
	    // To continue
	    index00 = d_index3d(i-1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u00 = getStateXvel(intfc_state);
		v00 = getStateYvel(intfc_state);
		w00 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u00 = u[index00];
		v00 = v[index00];
		w00 = w[index00];
	    }
	    index01 = d_index3d(i+1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,EAST,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u01 = getStateXvel(intfc_state);
		v01 = getStateYvel(intfc_state);
		w01 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u01 = u[index01];
		v01 = v[index01];
		w01 = w[index01];
	    }
	    index10 = d_index3d(i,j-1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,SOUTH,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u10 = getStateXvel(intfc_state);
		v10 = getStateYvel(intfc_state);
		w10 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u10 = u[index10];
		v10 = v[index10];
		w10 = w[index10];
	    }
	    index11 = d_index3d(i,j+1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,NORTH,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u11 = getStateXvel(intfc_state);
		v11 = getStateYvel(intfc_state);
		w11 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u11 = u[index11];
		v11 = v[index11];
		w11 = w[index11];
	    }
	    index20 = d_index3d(i,j,k-1,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u20 = getStateXvel(intfc_state);
		v20 = getStateYvel(intfc_state);
		w20 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u20 = u[index20];
		v20 = v[index20];
		w20 = w[index20];
	    }
	    index21 = d_index3d(i,j,k+1,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u21 = getStateXvel(intfc_state);
		v21 = getStateYvel(intfc_state);
		w21 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u21 = u[index21];
		v21 = v[index21];
		w21 = w[index21];
	    }


	    /*
	    cell_center[index].m_state.m_U[0] += -m_dt*(
				burger_flux(u00,u0,u01)/top_h[0] +
				linear_flux(v0,u10,u0,u11)/top_h[1] +
				linear_flux(w0,u20,u0,u21)/top_h[2]);

	    cell_center[index].m_state.m_U[1] += - m_dt*(
				linear_flux(u0,v00,v0,v01)/top_h[0] +
			 	burger_flux(v10,v0,v11)/top_h[1] +
				linear_flux(w0,v20,v0,v21)/top_h[2]);

	    cell_center[index].m_state.m_U[2] += - m_dt*(
				linear_flux(u0,w00,w0,w01)/top_h[0] +
				linear_flux(v0,w10,w0,w11)/top_h[1] +
			 	burger_flux(w20,w0,w21)/top_h[2]);
	    */

	    double r = cell_center[index].m_coords[2];

            cell_center[index].m_state.m_U[0] += - m_dt*(
                                burger_flux(u00,u0,u01)/(r*top_h[0]) +
                                linear_flux(v0,u10,u0,u11)/top_h[1] +
                                linear_flux(w0,u20,u0,u21)/top_h[2]);

            cell_center[index].m_state.m_U[1] += - m_dt*(
                                linear_flux(u0,v00,v0,v01)/(r*top_h[0]) +
                                burger_flux(v10,v0,v11)/top_h[1] +
                                linear_flux(w0,v20,v0,v21)/top_h[2]);

            cell_center[index].m_state.m_U[2] += - m_dt*(
                                linear_flux(u0,w00,w0,w01)/(r*top_h[0]) +
                                linear_flux(v0,w10,w0,w11)/top_h[1] +
                                burger_flux(w20,w0,w21)/top_h[2]);
          
            //SOURCE TERM
	    //
	    //How about calculating it in the diffusion solver?
	   /* 
            cell_center[index].m_state.m_U[0] += -(m_dt*u0*w0/r);

            cell_center[index].m_state.m_U[2] += (m_dt*u0*u0/r);
*/
	    
	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]) +
		    fabs(cell_center[index].m_state.m_U[2]);
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.m_U[l] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	FT_FreeThese(3,u,v,w);

}	/* end computeAdvection3d */



void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_1st_decoupled(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
        double coords[MAXD],crx_coords[MAXD];
	double coeff[18],mu0,rho,rhs,U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	int num_iter;
	double rel_residual;
	double r, rr, redge[2];

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 9, 9);
	// 7theta + 2r  for the first equation
	// 7z for the second equation
	// 7r + 2theta for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];
	    U2_center = cell_center[index].m_state.m_U[2];

	  
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);
		}
                else
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }

            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


    //Setting the coeffecients for the first equation

	    rr = cell_center[index].m_coords[2] * cell_center[index].m_coords[2];
	    r = cell_center[index].m_coords[2];
	    redge[0] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[4]].m_coords[2]);
	    redge[1] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[5]].m_coords[2]);


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++) //For the second order boundary
	    {
		if (I_nb[nb] == -1)
		{
		    coeff[nb] = 2.0*coeff[nb];
		}
	    }

	    solver.Set_A(I*3,I*3,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+0.5*m_dt*mu0/(rho*rr));
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5]-0.5*m_dt*mu0/(rho*rr))*cell_center[index].m_state.m_U[0];


	    coeff[6] = -0.5*m_dt/rho * mu0/(rr*top_h[0]);
	    coeff[7] =  0.5*m_dt/rho * mu0/(rr*top_h[0]);


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }

	    if (I_nb[0] != -1)
	    {
		solver.Set_A(I*3,I_nb[0]*3+2,-coeff[6]);
		rhs += coeff[6]*U2_nb[0];
	    }
	    else
	    {
		solver.Set_A(I*3,I*3+2,coeff[6]);
		rhs += -coeff[6]*U2_center;
		coeff[6] = coeff[6] * 2.0;
		rhs += 2.0*coeff[6]*U2_nb[0];
	    }

	    if (I_nb[1] != -1)
	    {
		solver.Set_A(I*3,I_nb[1]*3+2,-coeff[7]);
		rhs += coeff[7]*U2_nb[1];
	    }
	    else
	    {
		solver.Set_A(I*3,I*3+2,coeff[7]);
		rhs += -coeff[7]*U2_center;
		coeff[7] = coeff[7] * 2.0;
		rhs += 2.0*coeff[7]*U2_nb[1];
	    }

	    rhs -= m_dt*cell_center[index].m_state.m_U[0]*cell_center[index].m_state.m_U[2]/r; //Source term in advection step
	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;

	    solver.Set_b(I*3, rhs);

	    /************************************************************************/


    //Setting the coeffecients for the second equation


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+1,I*3+1,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[1];


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;

	    solver.Set_b(I*3+1, rhs);

	    /************************************************************************/

    //Setting the coeffecients for the third equation

	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+2,I*3+2,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+0.5*m_dt*mu0/(rho*rr));
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5]-0.5*m_dt*mu0/(rho*rr))*cell_center[index].m_state.m_U[2];


	    coeff[6] =  0.5*m_dt/rho * mu0/(rr*top_h[0]);
	    coeff[7] = -0.5*m_dt/rho * mu0/(rr*top_h[0]);


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }


	    if (I_nb[0] != -1)
	    {
		solver.Set_A(I*3+2, I_nb[0]*3, -coeff[6]);
		rhs += coeff[6]*U0_nb[0];
	    }
	    else
	    {
		solver.Set_A(I*3+2, I*3, coeff[6]);
		rhs += -coeff[6]*U0_center;
		coeff[6] = coeff[6] * 2.0;
		rhs += 2.0*coeff[6]*U0_nb[0];
	    }
	    
	    if (I_nb[1] != -1)
	    {
		solver.Set_A(I*3+2, I_nb[1]*3, -coeff[7]);
		rhs += coeff[7]*U0_nb[1];
	    }
	    else
	    {
		solver.Set_A(I*3+2, I*3, coeff[7]);
		rhs += -coeff[7]*U0_center;
		coeff[7] = coeff[7] * 2.0;
		rhs += 2.0*coeff[7]*U0_nb[1];
	    }

	    rhs += m_dt* cell_center[index].m_state.m_U[0] * cell_center[index].m_state.m_U[0] / r; //Source term in advection step
	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;

	    solver.Set_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
/*
	if (rel_residual > 1)
	{
	    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	    solver.Reset_x();
	    solver.Solve_GMRES();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	}
*/
	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_1st_decoupled: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }
        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);

}       /* end compDiffWithSmoothProperty3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computeProjection(void)
{
	int index, index_nb[6], size;
	double rhs, coeff[6], rho[6], rho0;
	int I,I_nb[6];
	int i,j,k,l,icoords[MAXD];
	INTERFACE *intfc = front->interf;
	double P_max,P_min;
	int icrds_Pmax[MAXD],icrds_Pmin[MAXD];
	COMPONENT comp;
	double aII,rho_nb[6];
	double coords[MAXD],crx_coords[MAXD];
	double **vel = iFparams->field->vel;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	
	max_value = 0.0;
	double value;
	double sum_div;
	sum_div = 0.0;
	int num_iter = 0;
	double rel_residual = 0.0;

	PETSc solver;
	solver.Create(ilower, iupper-1, 7, 7);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;

	start_clock("SetMatrixForProjection");
	setIndexMap();

	for (l = 0; l < dim; ++l)
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    vel[l][index] = cell_center[index].m_state.m_U[l];
	}

	/* Compute velocity divergence */
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index3d(i,j,k,top_gmax);
	    array[index] = computeFieldPointDiv(icoords,vel);
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.div_U = array[index];    
	}


	if(debugging("step_size"))
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        value = fabs(cell_center[index].m_state.div_U);
		sum_div = sum_div + cell_center[index].m_state.div_U * cell_center[index].m_coords[2] * top_h[0]*top_h[1]*top_h[2];
	        if(value > max_value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_div,1);
	    printf("\nThe summation of divergence of U is %.16g\n",sum_div);
	    pp_global_max(&max_value,1);
	    printf("\nThe max value of divergence of U is %.16g\n",max_value);
	    max_value = 0.0;
	}
	

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);
	    I_nb[0] = ijk_to_I[i-1][j][k];
	    I_nb[1] = ijk_to_I[i+1][j][k];
	    I_nb[2] = ijk_to_I[i][j-1][k];
	    I_nb[3] = ijk_to_I[i][j+1][k];
	    I_nb[4] = ijk_to_I[i][j][k-1];
	    I_nb[5] = ijk_to_I[i][j][k+1];
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	
	    rho0   = cell_center[index].m_state.m_rho;
	    for (l = 0; l < 6; ++l)
	    {
		if (I_nb[l] == -1)
		    index_nb[l] = index;
		rho[l] = 1.0/2*(rho0 + cell_center[index_nb[l]].m_state.m_rho);
		//coeff[l] = 1.0/rho[l]/sqr(top_h[l/2]);
	    }

	    double rr = sqr(cell_center[index].m_coords[2]);
	    double r = cell_center[index].m_coords[2];
	    double redge[2];
	    redge[0] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[4]].m_coords[2]);
	    redge[1] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[5]].m_coords[2]);


	    coeff[0] = 1.0/rho[0] * 1.0/(rr*top_h[0]*top_h[0]);
            coeff[1] = 1.0/rho[1] * 1.0/(rr*top_h[0]*top_h[0]);
            coeff[2] = 1.0/rho[2] * 1.0/(top_h[1]*top_h[1]);
            coeff[3] = 1.0/rho[3] * 1.0/(top_h[1]*top_h[1]);
            //coeff[4] = 1.0/rho[4] * ((1.0/(top_h[2]*top_h[2])) - (1.0/(2.0*r*top_h[2])));
            //coeff[5] = 1.0/rho[5] * ((1.0/(top_h[2]*top_h[2])) + (1.0/(2.0*r*top_h[2])));

	    coeff[4] = redge[0]/rho[4] / (r*top_h[2]*top_h[2]);
	    coeff[5] = redge[1]/rho[5] / (r*top_h[2]*top_h[2]);

            rhs = cell_center[index].m_state.div_U/accum_dt;

	    aII = 0.0;
	    for (l = 0; l < 6; ++l)
	    {
	    	if (I_nb[l] != -1)
		{
		    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
		}
	    }
            if (aII != 0.0)
	    {
                solver.Set_A(I,I,aII);
	    }
            else
            {
		printf("\nUsing the original pressure!\n");
                solver.Set_A(I,I,1.0);
                rhs = cell_center[index].m_state.m_P;
            }
            solver.Set_b(I,rhs);
	}
	stop_clock("SetMatrixForProjection");
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Petsc_Solver_Projection");
	solver.Solve_withPureNeumann();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
	if(rel_residual > 1e-10)
	{
	    printf("\n The solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	    solver.Reset_x();
	    solver.Solve_withPureNeumann_GMRES();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	}
	stop_clock("Petsc_Solver_Projection");

	start_clock("GetSolandScat");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("Incompress_Solver_Smooth_3D_Cylindrical::"
			"computeProjection: "
	       		"num_iter = %d, rel_residual = %le \n", 
			num_iter, rel_residual);
	
	P_max = -HUGE;		P_min = HUGE;
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    I = ijk_to_I[i][j][k];
	    array[index] = x[I-ilower];
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_phi = array[index];
	}

	if(debugging("step_size"))
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		value = fabs(cell_center[index].m_state.m_phi);
		if (value > max_value)
		    max_value = value;
	    }
	    pp_global_max(&max_value,1);
	    printf("\nThe max value of phi is %.16g\n",max_value);
	}

	FT_FreeThese(1,x);

	stop_clock("GetSolandScat");
}	/* end computeProjection3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computeNewVelocity(void)
{
	int i, j, k, l, index;
	double grad_phi[3], rho;
	COMPONENT comp;
	double speed;
	int icoords[MAXD];
	double r;

	max_speed = 0.0;

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_phi;
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
		for (l = 0; l < 3; ++l)
		    cell_center[index].m_state.m_U[l] = 0.0;
		continue;
	    }
	    rho = cell_center[index].m_state.m_rho;
	    r = cell_center[index].m_coords[2];
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGradPhi(icoords,array,grad_phi);
	    speed = 0.0;
	    for (l = 0; l < 3; ++l)
	    {
	    	cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
		speed += fabs(cell_center[index].m_state.m_U[l]);
	    }

	    if (speed > max_speed)
		max_speed = speed;
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.m_U[l] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
}	/* end computeNewVelocity3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computeSourceTerm_Adv(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
}

void Incompress_Solver_Smooth_3D_Cylindrical::computeSourceTerm(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
}
void Incompress_Solver_Smooth_3D_Cylindrical::computeSourceTerm(double *coords, double t, L_STATE &state) 
{
	computeSourceTerm(coords, state);
}

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void Incompress_Solver_Smooth_3D_Cylindrical::solve(double dt)
{
        printf("\nEntering solve, the dt for this solve is : %.16g\n",dt);

	m_t_old = front->time;
	m_t_int = front->time + dt/2.0;
	m_t_new = front->time + dt;

	printf("\nIn solve, old time = %.16g, new time = %.16g, step = %d\n", m_t_old, m_t_new, front->step);

	static boolean first = YES;
	if (first)
	{
	    accum_dt = 0.0;
	    first = NO;
	}
	m_dt = dt;
	max_speed = 0.0;

	setDomain();

	start_clock("SetCompAndIndex");

	setComponent();
	if (debugging("trace"))
	    printf("Passed setComponent()\n");
	setCornerComp();
	if (debugging("trace"))
	    printf("Passed setCornerComp()\n");

	setGlobalIndex();
	if (debugging("trace"))
	    printf("Passed setGlobalIndex()\n");

/*
	if (iFparams->movie_option->output_droplet_analysis == YES)
	{
	    start_clock("setTempFrontAndRedist");
	    setTempFrontAndRedist();	
	    if (debugging("trace"))
		printf("Passed setTempFrontAndRedist()\n");
	    stop_clock("setTempFrontAndRedist");

	    start_clock("SettingTriGlobalIndex");

	    setGlobalTriIndex();
	    if (debugging("trace"))
		printf("Passed setGlobalTriIndex()\n");

	    commGlobalTriIndex();
	    if (debugging("trace"))
		printf("Passed commGlobalTriIndex()\n");

	    commGlobalTriIndex2();
	    if (debugging("trace"))
		printf("Passed commGlobalTriIndex2()\n");


	    stop_clock("SettingTriGlobalIndex");
	}
	else
	    tempfront = NULL;
*/
	stop_clock("SetCompAndIndex");

	//start_clock("setSmoProOnePhase()\n");
	//setSmoProOnePhase();
	//stop_clock("setSmoProOnePhase()\n");

	//start_clock("setSmoProExact()\n");
	//setSmoProExact();
	//stop_clock("setSmoProExact()\n");

	setRhoMuOld();
	start_clock("setSmoothedProperties");
	setSmoothedProperties();
	stop_clock("setSmoothedProperties");

	if (front->step == 0 || front->regrid_restart == YES)
	    setRhoMuOld();


	if (debugging("trace"))
	    printf("Passed setSmoothedProperties()\n");
	
        if(front->turbulence_model == YES)
        {
            start_clock("compSGS");
            computeSubgridModel();
            stop_clock("compSGS");
        }

	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	//computeAdvection_test(); //Discretize the equation using the cylindrical Paper, result seems to be similar (First Order)
	//computeAdvection();  //First order discretization
	//compAdvectionTerm_decoupled(); //Second order convection term using slope limiter for first derivative (BELL 1989)
	compAdvectionTerm_coupled();
	stop_clock("computeAdvection");

	if (debugging("trace"))
	    printf("max_speed after computeAdvection(): %20.14f\n",
				max_speed);

	start_clock("compDiffWithSmoothProperty");
	//compDiffWithSmoothProperty_1st_decoupled_source();
	//compDiffWithSmoothProperty_1st_decoupled_test(); //Discretize the equatino using the cylindrical Paper, result seems to be similar
	//compDiffWithSmoothProperty_1st_decoupled();
	//compDiffWithSmoothProperty_2nd_decoupled(); //2nd order diffusion solver with the advection source terms
	compDiffWithSmoothProperty_2nd_coupled();
	//compDiffWithSmoothProperty_2nd_decoupled_Shuqiang(); //2nd order diffusion solver by Shuqiang
	stop_clock("compDiffWithSmoothProperty");

	if (debugging("trace"))
	    printf("max_speed after compDiffWithSmoothProperty(): %20.14f\n",
				max_speed);

	// 2) projection step
	accum_dt += m_dt;
	if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    computeProjection();
	    //computeProjection_Shuqiang(); //Projection step by Shuqiang
	    stop_clock("computeProjection");

	    start_clock("computePressure");
	    computePressure();
	    stop_clock("computePressure");

	    start_clock("computeNewVelocity");
	    computeNewVelocity();
	    stop_clock("computeNewVelocity");
	    accum_dt = 0.0;
	}

	if (debugging("sample_velocity"))
	{
	    sampleVelocity();
	}

	if (debugging("trace"))
	    printf("max_speed after computeNewVelocity(): %20.14f\n",
				max_speed);

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
}	/* end solve */


double Incompress_Solver_Smooth_3D_Cylindrical::getVorticityX(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dy,dz;
	double vorticity;

	dy = top_h[1];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j-1,k,top_gmax);
	index01 = d_index3d(i,j+1,k,top_gmax);
	index10 = d_index3d(i,j,k-1,top_gmax);
	index11 = d_index3d(i,j,k+1,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[2];
	v01 =  cell_center[index01].m_state.m_U[2];
	v10 =  cell_center[index10].m_state.m_U[1];
	v11 = -cell_center[index11].m_state.m_U[1];

	vorticity = (v00 + v01)/2.0/dz + (v10 + v11)/2.0/dy;
	return vorticity;
}	/* end getVorticityX */

double Incompress_Solver_Smooth_3D_Cylindrical::getVorticityY(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dz;
	double vorticity;

	dx = top_h[0];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j,k-1,top_gmax);
	index01 = d_index3d(i,j,k+1,top_gmax);
	index10 = d_index3d(i-1,j,k,top_gmax);
	index11 = d_index3d(i+1,j,k,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[0];
	v01 =  cell_center[index01].m_state.m_U[0];
	v10 =  cell_center[index10].m_state.m_U[2];
	v11 = -cell_center[index11].m_state.m_U[2];

	vorticity = (v00 + v01)/2.0/dx + (v10 + v11)/2.0/dz;
	return vorticity;
}	/* end getVorticityY */

double Incompress_Solver_Smooth_3D_Cylindrical::getVorticityZ(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i-1,j,k,top_gmax);
	index01 = d_index3d(i+1,j,k,top_gmax);
	index10 = d_index3d(i,j-1,k,top_gmax);
	index11 = d_index3d(i,j+1,k,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[1];
	v01 =  cell_center[index01].m_state.m_U[1];
	v10 =  cell_center[index10].m_state.m_U[0];
	v11 = -cell_center[index11].m_state.m_U[0];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticityZ */


void Incompress_Solver_Smooth_3D_Cylindrical::copyMeshStates()
{
	int i,j,k,d,index;
	double **vel = field->vel;
	double *pres = field->pres;
	double *vort = field->vort;
	double **vort3d = field->vort3d;

	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	for (k = kmin; k <= kmax; ++k)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    if (ifluid_comp(top_comp[index]))
	    {
		pres[index] = cell_center[index].m_state.m_P;
	    	vel[0][index] = cell_center[index].m_state.m_U[0];
	    	vel[1][index] = cell_center[index].m_state.m_U[1];
	    	vel[2][index] = cell_center[index].m_state.m_U[2];
		vort3d[0][index] = getVorticityX(i,j,k);
		vort3d[1][index] = getVorticityY(i,j,k);
		vort3d[2][index] = getVorticityZ(i,j,k);
	    }
	    else
	    {
	    	pres[index] = 0.0;
		for (d = 0; d < 3; ++d)
		{
		    vel[d][index] = 0.0;
		    vort3d[d][index] = 0.0;
		}
	    }
	}
	FT_ParallelExchGridArrayBuffer(pres,front);
	FT_ParallelExchGridArrayBuffer(vel[0],front);
	FT_ParallelExchGridArrayBuffer(vel[1],front);
	FT_ParallelExchGridArrayBuffer(vel[2],front);
	FT_ParallelExchGridArrayBuffer(vort3d[0],front);
	FT_ParallelExchGridArrayBuffer(vort3d[1],front);
	FT_ParallelExchGridArrayBuffer(vort3d[2],front);
}	/* end copyMeshStates */


void Incompress_Solver_Smooth_3D_Cylindrical::
	compDiffWithSmoothProperty_1st_decoupled_source(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
        L_STATE state;
        INTERFACE *intfc = front->interf;
	double coords[MAXD], crx_coords[MAXD];
	double coeff[6],mu[6],mu0,rho,corner[6],rhs,U_nb[6];
        double speed;
        double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	double (*getStateVel[3])(POINTER) = 
			{getStateXvel,getStateYvel,getStateZvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	int num_iter;
	double rel_residual;

        setIndexMap();

	max_speed = 0.0;

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
            PETSc solver;
            solver.Create(ilower, iupper-1, 7, 7);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I  = ijk_to_I[i][j][k];
            	if (I == -1) continue;

            	index  = d_index3d(i,j,k,top_gmax);
            	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		comp = top_comp[index];

            	I_nb[0] = ijk_to_I[i-1][j][k]; //west
            	I_nb[1] = ijk_to_I[i+1][j][k]; //east
            	I_nb[2] = ijk_to_I[i][j-1][k]; //south
            	I_nb[3] = ijk_to_I[i][j+1][k]; //north
            	I_nb[4] = ijk_to_I[i][j][k-1]; //lower
            	I_nb[5] = ijk_to_I[i][j][k+1]; //upper


            	mu0   = cell_center[index].m_state.m_mu;
            	rho   = cell_center[index].m_state.m_rho;

            	for (nb = 0; nb < 6; nb++)
            	{
                    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                                comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			U_nb[nb] = getStateVel[l](intfc_state);
			if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			    wave_type(hs) == NEUMANN_BOUNDARY)
			    mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
		    }
                    else
		    {
                    	U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
			mu[nb] = 1.0/2*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
		    }
            	}

                /*
            	coeff[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = 0.5*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = 0.5*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);
                */

                double rr = sqr(cell_center[index].m_coords[2]);
                double r = cell_center[index].m_coords[2];
                coeff[0] = 0.5*m_dt/rho * mu[0]/(rr*top_h[0]*top_h[0]);
                coeff[1] = 0.5*m_dt/rho * mu[1]/(rr*top_h[0]*top_h[0]);
                coeff[2] = 0.5*m_dt/rho * mu[2]/(top_h[1]*top_h[1]);
                coeff[3] = 0.5*m_dt/rho * mu[3]/(top_h[1]*top_h[1]);
                coeff[4] = 0.5*m_dt/rho * ((mu[4]/(top_h[2]*top_h[2]))-(mu[4]/(2.0*top_h[2]*r)));
                coeff[5] = 0.5*m_dt/rho * ((mu[5]/(top_h[2]*top_h[2]))+(mu[5]/(2.0*top_h[2]*r)));


            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, state);

                //SOURCE TERM
                double sour[3];
                double temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7,
                       temp8, temp9, temp10, temp11, temp12, temp13, temp14,
                       temp15, temp16, temp17, temp18, temp19, temp20;

                    temp0 = cell_center[index].m_state.m_U[0];
                    temp1 = cell_center[index].m_state.m_U[1];
                    temp2 = cell_center[index].m_state.m_U[2];
                    temp3 = cell_center[index_nb[0]].m_state.m_U[0];
                    temp4 = cell_center[index_nb[0]].m_state.m_U[1];
                    temp5 = cell_center[index_nb[0]].m_state.m_U[2];
                    temp6 = cell_center[index_nb[1]].m_state.m_U[0];
                    temp7 = cell_center[index_nb[1]].m_state.m_U[1];
                    temp8 = cell_center[index_nb[1]].m_state.m_U[2];
                    temp9 = cell_center[index_nb[2]].m_state.m_U[0];
                    temp10 = cell_center[index_nb[2]].m_state.m_U[1];
                    temp11 = cell_center[index_nb[2]].m_state.m_U[2];
                    temp12 = cell_center[index_nb[3]].m_state.m_U[0];
                    temp13 = cell_center[index_nb[3]].m_state.m_U[1];
                    temp14 = cell_center[index_nb[3]].m_state.m_U[2];
                    temp15 = cell_center[index_nb[4]].m_state.m_U[0];
                    temp16 = cell_center[index_nb[4]].m_state.m_U[1];
                    temp17 = cell_center[index_nb[4]].m_state.m_U[2];
                    temp18 = cell_center[index_nb[5]].m_state.m_U[0];
                    temp19 = cell_center[index_nb[5]].m_state.m_U[1];
                    temp20 = cell_center[index_nb[5]].m_state.m_U[2];

                    if(l == 0)
                    {
                        sour[0] = (m_dt*cell_center[index].m_state.m_mu/cell_center[index].m_state.m_rho
                                                  * (((2.0/sqr(cell_center[index].m_coords[2]))
                                                  * (temp8
                                                  - temp5) / top_h[0])
                                                  - (temp0 / sqr(cell_center[index].m_coords[2]))));
                    }
                    if(l == 1)
                    {
                        sour[1] = 0.0;
                    }
                    if(l == 2)
                    {
                        sour[2] = - (m_dt*cell_center[index].m_state.m_mu/cell_center[index].m_state.m_rho
                                                  * (((2.0/sqr(cell_center[index].m_coords[2]))
                                                  *(temp6
                                                  - temp3) / top_h[0])
                                                  + (temp2 / sqr(cell_center[index].m_coords[2]))));
                    }

        	//first equation
            	solver.Set_A(I,I,1+coeff[0]+coeff[1]+coeff[2]+coeff[3]+
				coeff[4]+coeff[5]);
		rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*
		      		cell_center[index].m_state.m_U[l];

		for(nb = 0; nb < 6; nb++)
		{
		    if(I_nb[nb] != -1)
		    {
			solver.Set_A(I,I_nb[nb],-coeff[nb]);
			rhs += coeff[nb]*U_nb[nb];
		    }
		    else
			rhs += 2.0*coeff[nb]*U_nb[nb];
		}
		rhs += m_dt*state.m_U[l];
                rhs += sour[l];
		rhs += m_dt*cell_center[index].m_state.f_surf[l];
		rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
 
		solver.Set_b(I, rhs);

            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-10);

	    start_clock("Befor Petsc solve");
            solver.Solve_GMRES();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);
/*
	    if(rel_residual > 1)
	    {
		printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
            	solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }
	    stop_clock("After Petsc solve");
*/
            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
			"compDiffWithSmoothProperty_1st_decoupled_source: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ijk_to_I[i][j][k];
                index = d_index3d(i,j,k,top_gmax);
                if (I >= 0)
                {
                    cell_center[index].m_state.m_U[l] = x[I-ilower];
                }
                else
                {
                    cell_center[index].m_state.m_U[l] = 0.0;
                }
            }

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]) +
                    fabs(cell_center[index].m_state.m_U[2]);
            if (speed > max_speed)
                    max_speed = speed;
	}
        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);
}       /* end compDiffWithSmoothProperty3d_decoupled */


void Incompress_Solver_Smooth_3D_Cylindrical::computePressurePmI(void)
{
        int i,j,k,index;

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P += cell_center[index].m_state.m_phi;
	    array[index] = cell_center[index].m_state.m_P;
	}
	scatMeshArray();

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}

}        /* end computePressurePmI3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computePressurePmII(void)
{
        int i,j,k,index;
        double mu0;

        for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P = (cell_center[index].m_state.m_q +
				cell_center[index].m_state.m_phi -
                        	mu0*cell_center[index].m_state.div_U);
	    array[index] = cell_center[index].m_state.m_P;
	}
	scatMeshArray();

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}
}        /* end computePressurePmII3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computePressurePmIII(void)
{
        int i,j,k,index;
        double mu0;

        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P = 
				cell_center[index].m_state.m_phi -
                        	accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = 0.0;
	}
}        /* end computePressurePmIII3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computePressure(void)
{
	switch (iFparams->num_scheme)
	{
	case BELL_COLELLA:
	    computePressurePmII();
	    break;
	case KIM_MOIN:
	    computePressurePmII();
	    break;
	case SIMPLE:
	case PEROT_BOTELLA:
	    computePressurePmIII();
	    break;
	}
	computeGradientQ();
}	/* end computePressure */

void Incompress_Solver_Smooth_3D_Cylindrical::computeGradientQ(void)
{
    //printf("\nEnter computeGradientQ in cylindrical coordinate\n");
	int i,j,k,l,index;
	double *grad_q;
	int icoords[MAXD];

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_q;

	}
	scatMeshArray();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    grad_q = cell_center[index].m_state.grad_q;
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGrad(icoords,array,grad_q); //grad_q is the Gradient Operator in Cylindrical Coordinate
	}
	for (l = 0; l < dim; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		array[index] = cell_center[index].m_state.grad_q[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.grad_q[l] = array[index];
	    }
	}
}	/* end computeGradientQ3d */

#define		MAX_TRI_FOR_INTEGRAL		100
void Incompress_Solver_Smooth_3D_Cylindrical::surfaceTension_Fedkiw(
	double delta,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
	double sigma,
	double rho,
	double *weight)
{
	TRI *tri;
	POINT *p;
	double nor[MAXD];
	double kappa = 0.0;
	for (int m = 0; m < 3; m++)
	    nor[m] = 0.0;

	tri = Tri_of_hse(hse);

	for (int m = 0; m < 3; m++)
	{
	    p = Point_of_tri(tri)[m];
	    kappa += p->curvature * weight[m];
	    for (int n = 0; n < 3; n++)
		nor[n] += p->_nor[n]* weight[m];
	}
	kappa *= 2.0;

	for (int m = 0; m < 3; m++)
	    force[m] = -delta*sigma*kappa*nor[m]/rho;

	/*
	TriAndFirstRing(hse,hs,&num_tris,tri_list);
	for (i = 0; i < num_tris; ++i)
	{
	    kappa = 0.0;
	    tri = tri_list[i];
	    for (j = 0; j < 3; ++j) median[j] = 0.0;
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		for (k = 0; k < 3; ++k) 
		    median[k] += Coords(p)[k];
	    	GetFrontCurvature(p,Hyper_surf_element(tri),hs,
				&kappa_tmp,front);
		kappa += kappa_tmp;
		nor[j] = Tri_normal(tri)[j];
	    }
	    kappa /= 3.0;
	    kappa *= 2.0;
	    mag_nor = mag_vector(nor,3);
	    area = 0.5*mag_nor; //use the triangle area in cartesian coordinate
	    for (j = 0; j < 3; ++j)  
	    {
		nor[j] /= mag_nor;
		median[j] /= 3.0;
	    }
	    delta = smoothedDeltaFunction(coords,median);
	    if (delta == 0.0) continue;
	    for (j = 0; j < dim; ++j) 
	    {
		force[j] += delta*sigma*area*kappa*nor[j];
	    }
	}

	for (i = 0; i < dim; i++)
	    force[i] /= cellVolume;
	*/
}	/* end surfaceTension3d */

void Incompress_Solver_Smooth_3D_Cylindrical::setInitialCondition()
{
	int i,j,k,index,l;
	COMPONENT comp;
	double coords[MAXD];
	double speed;

	FT_MakeGridIntfc(front);
	setDomain();
	setComponent();
	setCornerComp();

        m_rho[0] = iFparams->rho1;
        m_rho[1] = iFparams->rho2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	m_sigma = iFparams->surf_tension;
	mu_min = rho_min = HUGE;

	max_speed = 0.0;

	for (i = 0; i < 2; ++i)
	{
	    if (ifluid_comp(m_comp[i]))
	    {
        	mu_min = std::min(mu_min,m_mu[i]);
        	rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

       	double rho_1 = iFparams->rho1;
        double rho_2 = iFparams->rho2;
        double mu_1 = iFparams->mu1;
        double mu_2 = iFparams->mu2;
        //double Omega1 = 0.314;//0.314rad/ms
	//double Omega1 = 0.1;//0.1rad/s Re = 0.1*2.538*0.628/0.0102 = 15.626 laminar flow
        double R_1 = 2.538;
        double R_2 = 3.166;
        double Omega1 = iFparams->bvel[0][0]/R_1;//0.314rad/ms
        double a = R_1;
        double b = R_2;
        double c,A1,A2,B1,B2,term_3,term_4,term_5;
        double R_1s = sqr(R_1);
        double R_2s = sqr(R_2);

	// Initialize state at cell_center
        for (i = 0; i < cell_center.size(); i++)
        {
            getRectangleCenter(i, coords);
	    cell_center[i].m_state.setZero();
	    comp = top_comp[i];
	    if (getInitialState != NULL)
	    	(*getInitialState)(comp,coords,cell_center[i].m_state,dim,
						iFparams);

	    // Special initialization for contactor problem

	    /////////////// Two phase initialization with no surface tension//////////
	    /*
            c = 2.9; //Interface position at r = 2.9
            B1 = - Omega1*c/((1.0/c) - (c/R_1s) + (mu_1*c/(R_2s*mu_2)) - (mu_1/(c*mu_2)));
            A1 = Omega1 - (B1/R_1s);
            B2 = mu_1*B1/mu_2;
            A2 = - mu_1*B1/(R_2s*mu_2);

            term_3 = (sqr(A1)*sqr(cell_center[i].m_coords[2])/2.0)
                          + (2.0*A1*B1*log(cell_center[i].m_coords[2]))
                          - (sqr(B1)/(2.0*sqr(cell_center[i].m_coords[2])))
                          - (sqr(A1)*R_1s/2.0) - (2.0*A1*B1*log(R_1)) + (sqr(B1)/(2.0*R_1s));
            term_4 = (sqr(A1)*sqr(c)/2.0)
                          + (2.0*A1*B1*log(c))
                          - (sqr(B1)/(2.0*sqr(c)))
                          - (sqr(A1)*R_1s/2.0) - (2.0*A1*B1*log(R_1)) + (sqr(B1)/(2.0*R_1s));
            term_5 = (sqr(A2)*sqr(cell_center[i].m_coords[2])/2.0)
                          + (2.0*A2*B2*log(cell_center[i].m_coords[2]))
                          - (sqr(B2)/(2.0*sqr(cell_center[i].m_coords[2])))
                          - (sqr(A2)*sqr(c)/2.0) - (2.0*A2*B2*log(c)) + (sqr(B2)/(2.0*sqr(c)));

            if(cell_center[i].m_coords[2] <= c)
            {
                cell_center[i].m_state.m_U[0] = (A1*cell_center[i].m_coords[2])
                                                + (B1/cell_center[i].m_coords[2]);
                cell_center[i].m_state.m_U[1] = 0.0;
                cell_center[i].m_state.m_U[2] = 0.0;
                cell_center[i].m_state.m_P = (rho_1*term_3) + 1.0;
                cell_center[i].m_state.m_exactU[0] = cell_center[i].m_state.m_U[0];
                cell_center[i].m_state.m_exactU[1] = cell_center[i].m_state.m_U[1];
                cell_center[i].m_state.m_exactU[2] = cell_center[i].m_state.m_U[2];
                cell_center[i].m_state.m_exactP = cell_center[i].m_state.m_P;

            }
	    else if(cell_center[i].m_coords[2] > c)
            {
                cell_center[i].m_state.m_U[0] = (A2*cell_center[i].m_coords[2])
                                                + (B2/cell_center[i].m_coords[2]);
                cell_center[i].m_state.m_U[1] = 0.0;
                cell_center[i].m_state.m_U[2] = 0.0;
                cell_center[i].m_state.m_P = (rho_1*term_4) + (rho_2*term_5) + 1.0;
                cell_center[i].m_state.m_exactU[0] = cell_center[i].m_state.m_U[0];
                cell_center[i].m_state.m_exactU[1] = cell_center[i].m_state.m_U[1];
                cell_center[i].m_state.m_exactU[2] = cell_center[i].m_state.m_U[2];
                cell_center[i].m_state.m_exactP = cell_center[i].m_state.m_P;

            }
	    */
	    

	    /////////////// End of Two phase initialization with no surfac tension //////////
	
	    /////////////// Two phase initialization with surface tension//////////

            if(iFparams->use_couette_init_vel == YES)
            {
            c = 2.9; //Interface position at r = 2.9
            B1 = - Omega1*c/((1.0/c) - (c/R_1s) + (mu_1*c/(R_2s*mu_2)) - (mu_1/(c*mu_2)));
            A1 = Omega1 - (B1/R_1s);
            B2 = mu_1*B1/mu_2;
            A2 = - mu_1*B1/(R_2s*mu_2);

            term_3 = (sqr(A1)*sqr(cell_center[i].m_coords[2])/2.0)
                          + (2.0*A1*B1*log(cell_center[i].m_coords[2]))
                          - (sqr(B1)/(2.0*sqr(cell_center[i].m_coords[2])))
                          - (sqr(A1)*R_1s/2.0) - (2.0*A1*B1*log(R_1)) + (sqr(B1)/(2.0*R_1s));
            term_4 = (sqr(A1)*sqr(c)/2.0)
                          + (2.0*A1*B1*log(c))
                          - (sqr(B1)/(2.0*sqr(c)))
                          - (sqr(A1)*R_1s/2.0) - (2.0*A1*B1*log(R_1)) + (sqr(B1)/(2.0*R_1s));
            term_5 = (sqr(A2)*sqr(cell_center[i].m_coords[2])/2.0)
                          + (2.0*A2*B2*log(cell_center[i].m_coords[2]))
                          - (sqr(B2)/(2.0*sqr(cell_center[i].m_coords[2])))
                          - (sqr(A2)*sqr(c)/2.0) - (2.0*A2*B2*log(c)) + (sqr(B2)/(2.0*sqr(c)));

            if(cell_center[i].m_coords[2] <= c)
            {
                cell_center[i].m_state.m_U[0] = (A1*cell_center[i].m_coords[2])
                                                + (B1/cell_center[i].m_coords[2]);
                cell_center[i].m_state.m_U[1] = 0.0;
                cell_center[i].m_state.m_U[2] = 0.0;

                cell_center[i].m_state.m_P = (rho_1*term_3) + 1.0;
                cell_center[i].m_state.m_exactU[0] = cell_center[i].m_state.m_U[0];
                cell_center[i].m_state.m_exactU[1] = cell_center[i].m_state.m_U[1];
                cell_center[i].m_state.m_exactU[2] = cell_center[i].m_state.m_U[2];
                cell_center[i].m_state.m_exactP = cell_center[i].m_state.m_P;
            }
	    else if(cell_center[i].m_coords[2] > c)
            {
                cell_center[i].m_state.m_U[0] = (A2*cell_center[i].m_coords[2])
                                                + (B2/cell_center[i].m_coords[2]);
                cell_center[i].m_state.m_U[1] = 0.0;
                cell_center[i].m_state.m_U[2] = 0.0;
                cell_center[i].m_state.m_P = (rho_1*term_4) + (rho_2*term_5) + 1.0 - m_sigma/c;
                cell_center[i].m_state.m_exactU[0] = cell_center[i].m_state.m_U[0];
                cell_center[i].m_state.m_exactU[1] = cell_center[i].m_state.m_U[1];
                cell_center[i].m_state.m_exactU[2] = cell_center[i].m_state.m_U[2];
                cell_center[i].m_state.m_exactP = cell_center[i].m_state.m_P;

            }
	    
	    /////////////// End of Two phase initialization with surface tension //////////
	    
    
	    ///////////// One phase Initialization /////////////
/*
	    A1 = (Omega1*R_1s)/(R_1s-R_2s);
            B1 = -A1*R_2s;
            term_3 = (sqr(A1)*sqr(cell_center[i].m_coords[2])/2.0)
                          + (2.0*A1*B1*log(cell_center[i].m_coords[2]))
                          - (sqr(B1)/(2.0*sqr(cell_center[i].m_coords[2])))
                          - (sqr(A1)*R_1s/2.0) - (2.0*A1*B1*log(R_1)) + (sqr(B1)/(2.0*R_1s));

            cell_center[i].m_state.m_U[0] = (A1*cell_center[i].m_coords[2])
                                            + (B1/cell_center[i].m_coords[2]);
            cell_center[i].m_state.m_U[1] = 0.0;
            cell_center[i].m_state.m_U[2] = 0.0;
            cell_center[i].m_state.m_P = (rho_1*term_3) + 1.0;
            cell_center[i].m_state.m_exactU[0] = cell_center[i].m_state.m_U[0];
            cell_center[i].m_state.m_exactU[1] = cell_center[i].m_state.m_U[1];
            cell_center[i].m_state.m_exactU[2] = cell_center[i].m_state.m_U[2];
            cell_center[i].m_state.m_exactP = cell_center[i].m_state.m_P;
*/
	    ///////////////// End of One phase initialization ////////

            }
            else
            {
            /* Initialization with zero velocity*/

            cell_center[i].m_state.m_U[0] = iFparams->init_vel[0];
            cell_center[i].m_state.m_U[1] = iFparams->init_vel[1];
            cell_center[i].m_state.m_U[2] = iFparams->init_vel[2];
            cell_center[i].m_state.m_P = 0.0;
            cell_center[i].m_state.m_q = 0.0;
            }
        }


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]) +
		    fabs(cell_center[index].m_state.m_U[2]);

	    if (speed > max_speed)
		max_speed = speed;
	}

	pp_global_max(&max_speed,1);

	for (l = 0; l < dim; l++) // Scatter the velocity after the initialization
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.m_U[l] = array[index];
	    }
	}

	double totalVolume, cellVolume;
	double averagep, averagep_test;
	double r;
	totalVolume = 0.0;
	averagep = 0.0;
	averagep_test = 0.0;
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    r = cell_center[index].m_coords[2];
	    cellVolume = r*top_h[0]*top_h[1]*top_h[2];
	    totalVolume += cellVolume;
	    averagep += cell_center[index].m_state.m_P * cellVolume;
	}

	printf("\ntotalVolume for the processor itself is %.16g\n", totalVolume);

	pp_global_sum(&totalVolume, 1);
	pp_global_sum(&averagep, 1);
	averagep = averagep / totalVolume;
	printf("\nThe averagep is %.16g\n",averagep);
	printf("\nThe totalVolume is %.16g\n",totalVolume);



	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = cell_center[index].m_state.m_P - averagep;
	    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
	    cell_center[index].m_state.m_exactP = cell_center[index].m_state.m_P;
	    array[index] = cell_center[index].m_state.m_P; //For scatter

	    r = cell_center[index].m_coords[2];
	    cellVolume = r*top_h[0]*top_h[1]*top_h[2];
	    averagep_test += cell_center[index].m_state.m_P * cellVolume;
	}

	pp_global_sum(&averagep_test,1);
	averagep_test = averagep_test / totalVolume;
	printf("\nThe average of p after adjusting is %.16g\n", averagep_test);
	printf("\nTotal volume is %.16g\n", totalVolume);

	scatMeshArray();

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}

	computeGradientQ();
        copyMeshStates();
	setAdvectionDt();

}       /* end setInitialCondition */

double Incompress_Solver_Smooth_3D_Cylindrical::computeFieldPointDiv(
        int *icoords,
        double **field)
{
        int index;
        COMPONENT comp;
        int i, j,k,nb;
	int index_nb[6];
        double div,rur_nb[2],utheta_nb[2],uz_nb[2],ur;
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
      	double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];
	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	ur = field[2][index];
	double r, r_nb[2];
	r = cell_center[index].m_coords[2];
	r_nb[0] = cell_center[index_nb[4]].m_coords[2];
	r_nb[1] = cell_center[index_nb[5]].m_coords[2];
	
        const int  nn = pp_numnodes();
        int        myid = pp_mynode();
        int   *ppgmax = front->pp_grid->gmax;
        int   ppx = myid % ppgmax[0];
        int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
        int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    utheta_nb[0] = getStateXvel(intfc_state);
	else
	    utheta_nb[0] = (field[0][index] + field[0][index_nb[0]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,EAST,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    utheta_nb[1] = getStateXvel(intfc_state);
	else
	    utheta_nb[1] = (field[0][index] + field[0][index_nb[1]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,SOUTH,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    uz_nb[0] = getStateYvel(intfc_state);
	else
	    uz_nb[0] = (field[1][index] + field[1][index_nb[2]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,NORTH,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    uz_nb[1] = getStateYvel(intfc_state);
	else
	    uz_nb[1] = (field[1][index] + field[1][index_nb[3]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            rur_nb[0] = getStateZvel(intfc_state) * ((r + r_nb[0])/2.0);
        else
        {
        if ( (ppz == 0 && k == kmin) )
	    rur_nb[0] = 0.0;
            //rur_nb[0] = iFparams->bvel[0][2] * ((r + r_nb[0])/2.0);
	else
	    rur_nb[0] = ((field[2][index] + field[2][index_nb[4]])/2.0) * ((r + r_nb[0])/2.0);
        }

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            rur_nb[1] = getStateZvel(intfc_state) * ((r + r_nb[1])/2.0);
        else
        {
        if ( (ppz == ppgmax[2]-1 && k == kmax) )
	    rur_nb[1] = 0.0;
            //rur_nb[1] = iFparams->bvel[1][2] * ((r + r_nb[1])/2.0);
	else
	    rur_nb[1] = ((field[2][index] + field[2][index_nb[5]])/2.0) * ((r + r_nb[1])/2.0);
        }

        div = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) + (uz_nb[1] - uz_nb[0])/top_h[1] + (rur_nb[1] - rur_nb[0])/(r*top_h[2]);
        return div;
}       /* end computeFieldPointDiv */

void Incompress_Solver_Smooth_3D_Cylindrical::computeFieldPointGradPhi(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
	int index_nb[6];
        double p_nbedge[6],p0;  //the p values on the cell edges and cell center
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];
	
	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];
	p0 = field[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        const int  nn = pp_numnodes();
        int        myid = pp_mynode();
        int   *ppgmax = front->pp_grid->gmax;
        int   ppx = myid % ppgmax[0];
        int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
        int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

	for (nb = 0; nb < 6; nb++)
	{
	    if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
                
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state(hs) == NULL)
		{
		    p_nbedge[nb] = 0.0;
		    printf("\n Now using a flow through boundary for pressure!!\n");
		}
		else
		{
		    p_nbedge[nb] = p0;
		    /*
		    if (nb == 0)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[1]] - p0);
		    else if (nb == 1)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[0]]);
		    else if (nb == 2)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[3]] - p0);
		    else if (nb == 3)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[2]]);
		    else if (nb == 4)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[5]] - p0);
		    else if (nb == 5)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[4]]);
		    */
		}
	    }
	    else
	    {
		if ( (ppz == 0 && k == kmin && dir[nb] == LOWER) || (ppz == ppgmax[2]-1 && k == kmax && dir[nb] == UPPER) )
		    p_nbedge[nb] = p0;
		else
		    p_nbedge[nb] = (p0 + field[index_nb[nb]])/2.0;
	    }
	}

	double r;
	r = cell_center[index].m_coords[2];
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/(r*top_h[0]);
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
	grad_field[2] = (p_nbedge[5] - p_nbedge[4])/top_h[2];
}



void Incompress_Solver_Smooth_3D_Cylindrical::computeFieldPointGrad(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
	int index_nb[6];
        double p_nbedge[6],p0;  //the p values on the cell edges and cell center
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];
	
	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];
	p0 = field[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        const int  nn = pp_numnodes();
        int        myid = pp_mynode();
        int   *ppgmax = front->pp_grid->gmax;
        int   ppx = myid % ppgmax[0];
        int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
        int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

	for (nb = 0; nb < 6; nb++)
	{
	    if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
                
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state(hs) == NULL)
		{
		    p_nbedge[nb] = 0.0;
		    printf("\n Now using a flow through boundary for pressure!!\n");
		}
		else
		{
		    //p_nbedge[nb] = p0;

		    if (nb == 0)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[1]] - p0);
		    else if (nb == 1)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[0]]);
		    else if (nb == 2)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[3]] - p0);
		    else if (nb == 3)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[2]]);
		    else if (nb == 4)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[5]] - p0);
		    else if (nb == 5)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[4]]);

		}
	    }
	    else
	    {
		if ( (ppz == 0 && k == kmin && dir[nb] == LOWER) || (ppz == ppgmax[2]-1 && k == kmax && dir[nb] == UPPER) )
		{
		    if (nb == 0)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[1]] - p0);
		    else if (nb == 1)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[0]]);
		    else if (nb == 2)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[3]] - p0);
		    else if (nb == 3)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[2]]);
		    else if (nb == 4)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[5]] - p0);
		    else if (nb == 5)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[4]]);
		}
		else
		    p_nbedge[nb] = (p0 + field[index_nb[nb]])/2.0;
	    }
	}

	double r;
	r = cell_center[index].m_coords[2];
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/(r*top_h[0]);
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
	grad_field[2] = (p_nbedge[5] - p_nbedge[4])/top_h[2];
}



void Incompress_Solver_Smooth_3D_Cylindrical::computeError_part(FILE* outfile)
{
    int i,j,k,l,index;
    double error0, error1, error2;

    double L1error_0_inner, L2error_0_inner, Linferror_0_inner;
    double L1error_1_inner, L2error_1_inner, Linferror_1_inner;
    double L1error_2_inner, L2error_2_inner, Linferror_2_inner;

    double L1error_0_outer, L2error_0_outer, Linferror_0_outer;
    double L1error_1_outer, L2error_1_outer, Linferror_1_outer;
    double L1error_2_outer, L2error_2_outer, Linferror_2_outer;

    double errorp;
    double L1error_p_inner, L2error_p_inner, Linferror_p_inner;
    double L1error_p_outer, L2error_p_outer, Linferror_p_outer;

    L1error_0_inner = L1error_1_inner = L1error_2_inner = 0.0;
    L1error_0_outer = L1error_1_outer = L1error_2_outer = 0.0;

    L2error_0_inner = L2error_1_inner = L2error_2_inner = 0.0;
    L2error_0_outer = L2error_1_outer = L2error_2_outer = 0.0;

    Linferror_0_inner = Linferror_1_inner = Linferror_2_inner = 0.0;
    Linferror_0_outer = Linferror_1_outer = Linferror_2_outer = 0.0;

    L1error_p_inner = L2error_p_inner = Linferror_p_inner = 0.0;
    L1error_p_outer = L2error_p_outer = Linferror_p_outer = 0.0;

    double cellVolume;
    double r;

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);
	r = cell_center[index].m_coords[2];

	error0 = fabs(cell_center[index].m_state.m_U[0] - cell_center[index].m_state.m_exactU[0]);
	error1 = fabs(cell_center[index].m_state.m_U[1] - cell_center[index].m_state.m_exactU[1]);
	error2 = fabs(cell_center[index].m_state.m_U[2] - cell_center[index].m_state.m_exactU[2]);

	errorp = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);

	cellVolume = top_h[0]*top_h[1]*top_h[2]*r;

	if(r <= 2.9)
	{
	L1error_0_inner = L1error_0_inner + error0 * cellVolume;
	L1error_1_inner = L1error_1_inner + error1 * cellVolume;
	L1error_2_inner = L1error_2_inner + error2 * cellVolume;
	L1error_p_inner = L1error_p_inner + errorp * cellVolume;

	L2error_0_inner = L2error_0_inner + error0*error0*cellVolume;
	L2error_1_inner = L2error_1_inner + error1*error1*cellVolume;
	L2error_2_inner = L2error_2_inner + error2*error2*cellVolume;
	L2error_p_inner = L2error_p_inner + errorp*errorp*cellVolume;

	if (error0 >= Linferror_0_inner)
	    Linferror_0_inner = error0;
	if (error1 >= Linferror_1_inner)
	    Linferror_1_inner = error1;
	if (error2 >= Linferror_2_inner)
	    Linferror_2_inner = error2;
	if (errorp >= Linferror_p_inner)
	    Linferror_p_inner = errorp;
	}
	if(r >= 2.9)
	{
	L1error_0_outer = L1error_0_outer + error0 * cellVolume;
	L1error_1_outer = L1error_1_outer + error1 * cellVolume;
	L1error_2_outer = L1error_2_outer + error2 * cellVolume;
	L1error_p_outer = L1error_p_outer + errorp * cellVolume;

	L2error_0_outer = L2error_0_outer + error0*error0*cellVolume;
	L2error_1_outer = L2error_1_outer + error1*error1*cellVolume;
	L2error_2_outer = L2error_2_outer + error2*error2*cellVolume;
	L2error_p_outer = L2error_p_outer + errorp*errorp*cellVolume;

	if (error0 >= Linferror_0_outer)
	    Linferror_0_outer = error0;
	if (error1 >= Linferror_1_outer)
	    Linferror_1_outer = error1;
	if (error2 >= Linferror_2_outer)
	    Linferror_2_outer = error2;
	if (errorp >= Linferror_p_outer)
	    Linferror_p_outer = errorp;
	}

    }

    pp_global_max(&Linferror_0_inner,1);
    pp_global_max(&Linferror_0_outer,1);
    pp_global_max(&Linferror_1_inner,1);
    pp_global_max(&Linferror_1_outer,1);
    pp_global_max(&Linferror_2_inner,1);
    pp_global_max(&Linferror_2_outer,1);
    pp_global_max(&Linferror_p_inner,1);
    pp_global_max(&Linferror_p_outer,1);


    pp_global_sum(&L1error_0_inner,1);
    pp_global_sum(&L1error_0_outer,1);
    pp_global_sum(&L1error_1_inner,1);
    pp_global_sum(&L1error_1_outer,1);
    pp_global_sum(&L1error_2_inner,1);
    pp_global_sum(&L1error_2_outer,1);
    pp_global_sum(&L1error_p_inner,1);
    pp_global_sum(&L1error_p_outer,1);

    pp_global_sum(&L2error_0_inner,1);
    pp_global_sum(&L2error_0_outer,1);
    pp_global_sum(&L2error_1_inner,1);
    pp_global_sum(&L2error_1_outer,1);
    pp_global_sum(&L2error_2_inner,1);
    pp_global_sum(&L2error_2_outer,1);
    pp_global_sum(&L2error_p_inner,1);
    pp_global_sum(&L2error_p_outer,1);

    L2error_0_inner = sqrt(L2error_0_inner);
    L2error_1_inner = sqrt(L2error_1_inner);
    L2error_2_inner = sqrt(L2error_2_inner);
    L2error_p_inner = sqrt(L2error_p_inner);


    L2error_0_outer = sqrt(L2error_0_outer);
    L2error_1_outer = sqrt(L2error_1_outer);
    L2error_2_outer = sqrt(L2error_2_outer);
    L2error_p_outer = sqrt(L2error_p_outer);

    if(pp_mynode() == 0)
    {
	fprintf(outfile, "\nAt time t = %.16g, step = %d\n", front->time, front->step);
	fprintf(outfile, "\nFor Inner Part\n");

	fprintf(outfile, "L1   error = { %.16g   %.16g   %.16g  %.16g}\n", L1error_0_inner, L1error_1_inner, L1error_2_inner, L1error_p_inner);
	fprintf(outfile, "L2   error = { %.16g   %.16g   %.16g  %.16g}\n", L2error_0_inner, L2error_1_inner, L2error_2_inner, L2error_p_inner);
	fprintf(outfile, "Linf error = { %.16g   %.16g   %.16g  %.16g}\n", Linferror_0_inner, Linferror_1_inner, Linferror_2_inner, Linferror_p_inner);

	fprintf(outfile, "\nFor Outer Part\n");

	fprintf(outfile, "L1   error = { %.16g   %.16g   %.16g  %.16g}\n", L1error_0_outer, L1error_1_outer, L1error_2_outer, L1error_p_outer);
	fprintf(outfile, "L2   error = { %.16g   %.16g   %.16g  %.16g}\n", L2error_0_outer, L2error_1_outer, L2error_2_outer, L2error_p_outer);
	fprintf(outfile, "Linf error = { %.16g   %.16g   %.16g  %.16g}\n", Linferror_0_outer, Linferror_1_outer, Linferror_2_outer, Linferror_p_outer);
    }

}

void Incompress_Solver_Smooth_3D_Cylindrical::computeError(FILE* outfile)
{
    int i,j,k,l,index;
    double error0, error1, error2;
    double L1error_0, L2error_0, Linferror_0, L1error_1, L2error_1, Linferror_1, L1error_2, L2error_2, Linferror_2;
    double errorp;
    double L1error_p, L2error_p, Linferror_p;
    L1error_0 = 0.0; L1error_1 = 0.0; L1error_2 = 0.0;
    L2error_0 = 0.0; L2error_1 = 0.0; L2error_2 = 0.0;
    Linferror_0 = 0.0; Linferror_1 = 0.0; Linferror_2 = 0.0;
    L1error_p = 0.0; L2error_p = 0.0; Linferror_p = 0.0;
    double cellVolume;
    double r0, r1;

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);

	error0 = fabs(cell_center[index].m_state.m_U[0] - cell_center[index].m_state.m_exactU[0]);
	error1 = fabs(cell_center[index].m_state.m_U[1] - cell_center[index].m_state.m_exactU[1]);
	error2 = fabs(cell_center[index].m_state.m_U[2] - cell_center[index].m_state.m_exactU[2]);

	errorp = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);


	r0 = cell_center[index].m_coords[2] - top_h[2]/2.0;
	r1 = cell_center[index].m_coords[2] + top_h[2]/2.0;


	cellVolume = top_h[0]*top_h[1]*(r1*r1 - r0*r0)/2.0;
	L1error_0 = L1error_0 + error0 * cellVolume;
	L1error_1 = L1error_1 + error1 * cellVolume;
	L1error_2 = L1error_2 + error2 * cellVolume;
	L1error_p = L1error_p + errorp * cellVolume;

	L2error_0 = L2error_0 + error0*error0*cellVolume;
	L2error_1 = L2error_1 + error1*error1*cellVolume;
	L2error_2 = L2error_2 + error2*error2*cellVolume;
	L2error_p = L2error_p + errorp*errorp*cellVolume;

	if (error0 >= Linferror_0)
	    Linferror_0 = error0;
	if (error1 >= Linferror_1)
	    Linferror_1 = error1;
	if (error2 >= Linferror_2)
	    Linferror_2 = error2;
	if (errorp >= Linferror_p)
	    Linferror_p = errorp;
    }

    pp_global_max(&Linferror_0,1);
    pp_global_max(&Linferror_1,1);
    pp_global_max(&Linferror_2,1);
    pp_global_max(&Linferror_p,1);

    pp_global_sum(&L1error_0,1);
    pp_global_sum(&L1error_1,1);
    pp_global_sum(&L1error_2,1);
    pp_global_sum(&L1error_p,1);

    pp_global_sum(&L2error_0,1);
    pp_global_sum(&L2error_1,1);
    pp_global_sum(&L2error_2,1);
    pp_global_sum(&L2error_p,1);

    L2error_0 = sqrt(L2error_0);
    L2error_1 = sqrt(L2error_1);
    L2error_2 = sqrt(L2error_2);
    L2error_p = sqrt(L2error_p);


    INTERFACE *intfc = front->interf;
    HYPER_SURF_ELEMENT *hse;
    HYPER_SURF *hs;
    POINT *p;
    double Linferror_intfc, L1error_intfc, L2error_intfc, error_intfc;
    L1error_intfc = Linferror_intfc = L2error_intfc = 0.0;

    double L[3], U[3];
    for(int m = 0; m < 3; m++)
    {
	L[m] = comp_grid->L[m];
	U[m] = comp_grid->U[m];
    }


    double num_points = 0.0;
    (void) next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
	if(   Coords(p)[0] >= L[0] && Coords(p)[0] <= U[0] 
	   && Coords(p)[1] >= L[1] && Coords(p)[1] <= U[1]
	   && Coords(p)[2] >= L[2] && Coords(p)[2] <= U[2])
	{
	    num_points += 1.0;
	    error_intfc = fabs(Coords(p)[2] - 2.9);
	    L1error_intfc += error_intfc;
	    L2error_intfc += error_intfc*error_intfc;
	    if (error_intfc > Linferror_intfc)
		Linferror_intfc = error_intfc;
	}

    }
    pp_global_max(&Linferror_intfc, 1);
    pp_global_sum(&L1error_intfc, 1);
    pp_global_sum(&L2error_intfc, 1);
    pp_global_sum(&num_points, 1);

    L1error_intfc = L1error_intfc/num_points;
    L2error_intfc = L2error_intfc/num_points;

    L2error_intfc = sqrt(L2error_intfc);

    if (pp_mynode() == 0)
    {
	fprintf(outfile, "\nAt time t = %.16g, step = %d\n", front->time, front->step);

	fprintf(outfile, "\nInterface error L1 = %.16g,  L2 = %.16g,  Linf = %.16g\n", 
		L1error_intfc, L2error_intfc, Linferror_intfc);

	fprintf(outfile, "L1   error = { %.16g   %.16g   %.16g  %.16g}\n", L1error_0, L1error_1, L1error_2, L1error_p);
	fprintf(outfile, "L2   error = { %.16g   %.16g   %.16g  %.16g}\n", L2error_0, L2error_1, L2error_2, L2error_p);
	fprintf(outfile, "Linf error = { %.16g   %.16g   %.16g  %.16g}\n", Linferror_0, Linferror_1, Linferror_2, Linferror_p);
    }
}

void Incompress_Solver_Smooth_3D_Cylindrical::printExpandedMesh(char *out_name, bool binary)
{
    if(binary == YES)
    {
	if(hardware_is_little_endian())
	    printExpandedMesh_little_endian(out_name);
	else
	    printExpandedMesh_big_endian(out_name);
    }
    else
	printExpandedMesh_ascii(out_name);
}

void Incompress_Solver_Smooth_3D_Cylindrical::printExpandedMesh_big_endian(char *out_name)
{
    	double val[3];
	double ival[8];
        int   i,j,k,l,index,totalpoints;
        double coord_theta,coord_z,coord_r;
        int first;
        char dirname[256];
        char filename1[200],filename2[200],filename3[200],filename4[200],filename5[200],filename6[200];
        FILE *outfile1,*outfile2,*outfile3,*outfile4,*outfile5,*outfile6;
        INTERFACE *intfc = front->interf;
        RECT_GRID *gr = &topological_grid(intfc);
        int *ppgmax = front->pp_grid->gmax;
        int pointsr, pointsz, pointstheta;
        int gmax[3];
        gmax[0] = front->gmax[0];
        gmax[1] = front->gmax[1];
        gmax[2] = front->gmax[2];

        if(pp_mynode() == 0)
        {
            sprintf(dirname,"%s/domain-bdry",out_name);
	    if (!create_directory(dirname,YES))
            {
                screen("Cannot create directory %s\n",dirname);
                clean_up(ERROR);
            }

            sprintf(filename1,"%s/surface1.vtk",dirname);
            outfile1 = fopen(filename1,"wb");
            fprintf(outfile1,"# vtk DataFile Version 3.0\n");
            fprintf(outfile1,"Mesh of the topological grid\n");
            fprintf(outfile1,"BINARY\n");
            fprintf(outfile1,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = gmax[1] + 1;
            pointstheta = 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile1, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile1, "POINTS %d double\n", totalpoints);
            i = 0;
            for(j = 0; j <= pointsz-1; ++j)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile1);
            }
            fclose(outfile1);

            sprintf(filename2,"%s/surface2.vtk",dirname);
            outfile2 = fopen(filename2,"wb");
            fprintf(outfile2,"# vtk DataFile Version 3.0\n");
            fprintf(outfile2,"Mesh of the topological grid\n");
            fprintf(outfile2,"BINARY\n");
            fprintf(outfile2,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile2, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile2, "POINTS %d double\n", totalpoints);
            j = 0;
            for(i = 0; i <= pointstheta-1; ++i)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
  		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile2);
            }
            fclose(outfile2);

            sprintf(filename3,"%s/surface3.vtk",dirname);
            outfile3 = fopen(filename3,"wb");
            fprintf(outfile3,"# vtk DataFile Version 3.0\n");
            fprintf(outfile3,"Mesh of the topological grid\n");
            fprintf(outfile3,"BINARY\n");
            fprintf(outfile3,"DATASET STRUCTURED_GRID\n");
            pointsr = 1;
            pointsz = gmax[1] + 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile3, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile3, "POINTS %d double\n", totalpoints);
            k = 0;
            for(i = 0; i <= pointstheta-1; ++i)
            for(j = 0; j <= pointsz-1; ++j)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile3);
            }
            fclose(outfile3);

            sprintf(filename4,"%s/surface4.vtk",dirname);
            outfile4 = fopen(filename4,"wb");
            fprintf(outfile4,"# vtk DataFile Version 3.0\n");
            fprintf(outfile4,"Mesh of the topological grid\n");
            fprintf(outfile4,"BINARY\n");
            fprintf(outfile4,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = gmax[1] + 1;
            pointstheta = 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile4, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile4, "POINTS %d double\n", totalpoints);
            i = gmax[0];
            for(j = 0; j <= pointsz-1; ++j)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile4);
            }
            fclose(outfile4);

            sprintf(filename5,"%s/surface5.vtk",dirname);
            outfile5 = fopen(filename5,"wb");
            fprintf(outfile5,"# vtk DataFile Version 3.0\n");
            fprintf(outfile5,"Mesh of the topological grid\n");
            fprintf(outfile5,"BINARY\n");
            fprintf(outfile5,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile5, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile5, "POINTS %d double\n", totalpoints);
            j = gmax[1];
            for(i = 0; i <= pointstheta-1; ++i)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile5);
            }
            fclose(outfile5);

            sprintf(filename6,"%s/surface6.vtk",dirname);
            outfile6 = fopen(filename6,"wb");
            fprintf(outfile6,"# vtk DataFile Version 3.0\n");
            fprintf(outfile6,"Mesh of the topological grid\n");
            fprintf(outfile6,"BINARY\n");
            fprintf(outfile6,"DATASET STRUCTURED_GRID\n");
            pointsr = 1;
            pointsz = gmax[1] + 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile6, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile6, "POINTS %d double\n", totalpoints);
            k = gmax[2];
            for(i = 0; i <= pointstheta-1; ++i)
            for(j = 0; j <= pointsz-1; ++j)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile6);
            }
            fclose(outfile6);
        }
}

void Incompress_Solver_Smooth_3D_Cylindrical::printExpandedMesh_little_endian(char *out_name)
{
    	double val[3];
	double ival[8];
        int   i,j,k,l,index,totalpoints;
        double coord_theta,coord_z,coord_r;
        int first;
        char dirname[256];
        char filename1[200],filename2[200],filename3[200],filename4[200],filename5[200],filename6[200];
        FILE *outfile1,*outfile2,*outfile3,*outfile4,*outfile5,*outfile6;
        INTERFACE *intfc = front->interf;
        RECT_GRID *gr = &topological_grid(intfc);
        int *ppgmax = front->pp_grid->gmax;
        int pointsr, pointsz, pointstheta;
        int gmax[3];
        gmax[0] = front->gmax[0];
        gmax[1] = front->gmax[1];
        gmax[2] = front->gmax[2];

        if(pp_mynode() == 0)
        {
            sprintf(dirname,"%s/domain-bdry",out_name);
	    if (!create_directory(dirname,YES))
            {
                screen("Cannot create directory %s\n",dirname);
                clean_up(ERROR);
            }

            sprintf(filename1,"%s/surface1.vtk",dirname);
            outfile1 = fopen(filename1,"wb");
            fprintf(outfile1,"# vtk DataFile Version 3.0\n");
            fprintf(outfile1,"Mesh of the topological grid\n");
            fprintf(outfile1,"BINARY\n");
            fprintf(outfile1,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = gmax[1] + 1;
            pointstheta = 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile1, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile1, "POINTS %d double\n", totalpoints);
            i = 0;
            for(j = 0; j <= pointsz-1; ++j)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile1);
            }
            fclose(outfile1);

            sprintf(filename2,"%s/surface2.vtk",dirname);
            outfile2 = fopen(filename2,"wb");
            fprintf(outfile2,"# vtk DataFile Version 3.0\n");
            fprintf(outfile2,"Mesh of the topological grid\n");
            fprintf(outfile2,"BINARY\n");
            fprintf(outfile2,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile2, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile2, "POINTS %d double\n", totalpoints);
            j = 0;
            for(i = 0; i <= pointstheta-1; ++i)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
  		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile2);
            }
            fclose(outfile2);

            sprintf(filename3,"%s/surface3.vtk",dirname);
            outfile3 = fopen(filename3,"wb");
            fprintf(outfile3,"# vtk DataFile Version 3.0\n");
            fprintf(outfile3,"Mesh of the topological grid\n");
            fprintf(outfile3,"BINARY\n");
            fprintf(outfile3,"DATASET STRUCTURED_GRID\n");
            pointsr = 1;
            pointsz = gmax[1] + 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile3, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile3, "POINTS %d double\n", totalpoints);
            k = 0;
            for(i = 0; i <= pointstheta-1; ++i)
            for(j = 0; j <= pointsz-1; ++j)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile3);
            }
            fclose(outfile3);

            sprintf(filename4,"%s/surface4.vtk",dirname);
            outfile4 = fopen(filename4,"wb");
            fprintf(outfile4,"# vtk DataFile Version 3.0\n");
            fprintf(outfile4,"Mesh of the topological grid\n");
            fprintf(outfile4,"BINARY\n");
            fprintf(outfile4,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = gmax[1] + 1;
            pointstheta = 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile4, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile4, "POINTS %d double\n", totalpoints);
            i = gmax[0];
            for(j = 0; j <= pointsz-1; ++j)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile4);
            }
            fclose(outfile4);

            sprintf(filename5,"%s/surface5.vtk",dirname);
            outfile5 = fopen(filename5,"wb");
            fprintf(outfile5,"# vtk DataFile Version 3.0\n");
            fprintf(outfile5,"Mesh of the topological grid\n");
            fprintf(outfile5,"BINARY\n");
            fprintf(outfile5,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile5, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile5, "POINTS %d double\n", totalpoints);
            j = gmax[1];
            for(i = 0; i <= pointstheta-1; ++i)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile5);
            }
            fclose(outfile5);

            sprintf(filename6,"%s/surface6.vtk",dirname);
            outfile6 = fopen(filename6,"wb");
            fprintf(outfile6,"# vtk DataFile Version 3.0\n");
            fprintf(outfile6,"Mesh of the topological grid\n");
            fprintf(outfile6,"BINARY\n");
            fprintf(outfile6,"DATASET STRUCTURED_GRID\n");
            pointsr = 1;
            pointsz = gmax[1] + 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile6, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile6, "POINTS %d double\n", totalpoints);
            k = gmax[2];
            for(i = 0; i <= pointstheta-1; ++i)
            for(j = 0; j <= pointsz-1; ++j)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile6);
            }
            fclose(outfile6);
        }
}

void Incompress_Solver_Smooth_3D_Cylindrical::printExpandedMesh_ascii(char *out_name)
{
        int   i,j,k,l,index,totalpoints;
        double coord_theta,coord_z,coord_r;
        int first;
        //char filename[200];
        //FILE *outfile;
        char dirname[256];
        char filename1[200],filename2[200],filename3[200],filename4[200],filename5[200],filename6[200];
        FILE *outfile1,*outfile2,*outfile3,*outfile4,*outfile5,*outfile6;
        INTERFACE *intfc = front->interf;
        RECT_GRID *gr = &topological_grid(intfc);
        int *ppgmax = front->pp_grid->gmax;
        int pointsr, pointsz, pointstheta;
        int gmax[3];
        gmax[0] = front->gmax[0];
        gmax[1] = front->gmax[1];
        gmax[2] = front->gmax[2];

/*
        sprintf(dirname,"%s/P-%s",out_name,
                right_flush(pp_mynode(),4));
        sprintf(dirname,"%s/interface",dirname);
        if (!create_directory(dirname,YES))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s/mesh-p%s",dirname,right_flush(pp_mynode(),4));
#endif

        sprintf(filename,"%s.vtk",filename);
        i = 0;

	    outfile = fopen(filename,"w");
            fprintf(outfile,"# vtk DataFile Version 3.0\n");
            fprintf(outfile,"Mesh of the topological grid\n");
            fprintf(outfile,"ASCII\n");
            fprintf(outfile,"DATASET STRUCTURED_GRID\n");
            int pointsr, pointsz, pointstheta;
	    double thetamin, thetamax, rmin, rmax, zmin, zmax;
	    index = d_index3d(0,0,0,top_gmax);
	    thetamin = cell_center[index].m_coords[0] + (3.5*top_h[0]);
	    zmin = cell_center[index].m_coords[1] + (3.5*top_h[1]);
	    rmin = cell_center[index].m_coords[2] + (3.5*top_h[2]);
	    //index = d_index3d(top_gmax[0],top_gmax[1],top_gmax[2],top_gmax);
	    //thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    //zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    //rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;
        
	    //pointsr = top_gmax[2] + 2;
            //pointsz = top_gmax[1] + 2;
            //pointstheta = top_gmax[0] + 2;
            pointsr = top_gmax[2] - 6;
            pointsz = top_gmax[1] - 6;
            pointstheta = top_gmax[0] - 6;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile, "POINTS %d double\n", totalpoints);
            for(i = 0; i <= top_gmax[0]-7; ++i)
            for(j = 0; j <= top_gmax[1]-7; ++j)
            for(k = 0; k <= top_gmax[2]-7; ++k)
            {
                coord_theta = thetamin + i*top_h[0];
                coord_z = zmin + j*top_h[1];
                coord_r = rmin + k*top_h[2];
                fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
            }
            fclose(outfile);
*/

        if(pp_mynode() == 0)
        {
            sprintf(dirname,"%s/domain-bdry",out_name);
           if (!create_directory(dirname,YES))
            {
                screen("Cannot create directory %s\n",dirname);
                clean_up(ERROR);
            }

            sprintf(filename1,"%s/surface1.vtk",dirname);
            outfile1 = fopen(filename1,"w");
            fprintf(outfile1,"# vtk DataFile Version 3.0\n");
            fprintf(outfile1,"Mesh of the topological grid\n");
            fprintf(outfile1,"ASCII\n");
            fprintf(outfile1,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = gmax[1] + 1;
            pointstheta = 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile1, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile1, "POINTS %d double\n", totalpoints);
            i = 0;
            for(j = 0; j <= pointsz-1; ++j)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
                fprintf(outfile1, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
            }
            fclose(outfile1);

            sprintf(filename2,"%s/surface2.vtk",dirname);
            outfile2 = fopen(filename2,"w");
            fprintf(outfile2,"# vtk DataFile Version 3.0\n");
            fprintf(outfile2,"Mesh of the topological grid\n");
            fprintf(outfile2,"ASCII\n");
            fprintf(outfile2,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile2, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile2, "POINTS %d double\n", totalpoints);
            j = 0;
            for(i = 0; i <= pointstheta-1; ++i)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
                fprintf(outfile2, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
            }
            fclose(outfile2);

            sprintf(filename3,"%s/surface3.vtk",dirname);
            outfile3 = fopen(filename3,"w");
            fprintf(outfile3,"# vtk DataFile Version 3.0\n");
            fprintf(outfile3,"Mesh of the topological grid\n");
            fprintf(outfile3,"ASCII\n");
            fprintf(outfile3,"DATASET STRUCTURED_GRID\n");
            pointsr = 1;
            pointsz = gmax[1] + 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile3, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile3, "POINTS %d double\n", totalpoints);
            k = 0;
            for(i = 0; i <= pointstheta-1; ++i)
            for(j = 0; j <= pointsz-1; ++j)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
                fprintf(outfile3, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
            }
            fclose(outfile3);

            sprintf(filename4,"%s/surface4.vtk",dirname);
            outfile4 = fopen(filename4,"w");
            fprintf(outfile4,"# vtk DataFile Version 3.0\n");
            fprintf(outfile4,"Mesh of the topological grid\n");
            fprintf(outfile4,"ASCII\n");
            fprintf(outfile4,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = gmax[1] + 1;
            pointstheta = 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile4, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile4, "POINTS %d double\n", totalpoints);
            i = gmax[0];
            for(j = 0; j <= pointsz-1; ++j)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
                fprintf(outfile4, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
            }
            fclose(outfile4);

            sprintf(filename5,"%s/surface5.vtk",dirname);
            outfile5 = fopen(filename5,"w");
            fprintf(outfile5,"# vtk DataFile Version 3.0\n");
            fprintf(outfile5,"Mesh of the topological grid\n");
            fprintf(outfile5,"ASCII\n");
            fprintf(outfile5,"DATASET STRUCTURED_GRID\n");
            pointsr = gmax[2] + 1;
            pointsz = 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile5, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile5, "POINTS %d double\n", totalpoints);
            j = gmax[1];
            for(i = 0; i <= pointstheta-1; ++i)
            for(k = 0; k <= pointsr-1; ++k)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
                fprintf(outfile5, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
            }
            fclose(outfile5);

            sprintf(filename6,"%s/surface6.vtk",dirname);
            outfile6 = fopen(filename6,"w");
            fprintf(outfile6,"# vtk DataFile Version 3.0\n");
            fprintf(outfile6,"Mesh of the topological grid\n");
            fprintf(outfile6,"ASCII\n");
            fprintf(outfile6,"DATASET STRUCTURED_GRID\n");
            pointsr = 1;
            pointsz = gmax[1] + 1;
            pointstheta = gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile6, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile6, "POINTS %d double\n", totalpoints);
            k = gmax[2];
            for(i = 0; i <= pointstheta-1; ++i)
            for(j = 0; j <= pointsz-1; ++j)
            {
                coord_theta = gr->GL[0] + i*top_h[0];
                coord_z = gr->GL[1] + j*top_h[1];
                coord_r = gr->GL[2] + k*top_h[2];
                fprintf(outfile6, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
            }
            fclose(outfile6);
        }
}

/*
void Incompress_Solver_Smooth_3D_Cylindrical::printInteriorVelocity(char *out_name)
{
        int   i,j,k,l,index;
        double coord_theta,coord_z,coord_r;
        double vel_theta, vel_z, vel_r;
        int first;
        char filename[200];
        FILE *outfile;
        sprintf(filename,"%s/velo-visual-ts%s",out_name,
                right_flush(front->step,7));
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif
        sprintf(filename,"%s-liquid.vtk",filename);
        i = 0;
        if(iFparams->movie_option->plot_pres || iFparams->movie_option->plot_velo || iFparams->movie_option->plot_vort)
        {
	    outfile = fopen(filename,"w");
            fprintf(outfile,"# vtk DataFile Version 3.0\n");
            fprintf(outfile,"States of the whole computational domain\n");
            fprintf(outfile,"ASCII\n");
            fprintf(outfile,"DATASET STRUCTURED_GRID\n");
            int pointsr, pointsz, pointstheta;
        
	    pointsr = kmax - kmin + 2;
            pointsz = jmax - jmin + 2;
            pointstheta = imax - imin + 2;
            int totalpoints = pointsr*pointsz*pointstheta;
	    int totalcells = (pointsr-1)*(pointsz-1)*(pointstheta-1);

	    double thetamin, thetamax, zmin, zmax, rmin, rmax;

	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


            fprintf(outfile, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile, "POINTS %d double\n", totalpoints);
            for(i = 0; i <= imax - imin + 1; ++i)
            for(j = 0; j <= jmax - jmin + 1; ++j)
            for(k = 0; k <= kmax - kmin + 1; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
                fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
            }
            fprintf(outfile, "CELL_DATA %i\n", totalcells);
            if(iFparams->movie_option->plot_velo)
            {
                fprintf(outfile, "VECTORS velocity double\n");
    		for(i = imin; i <= imax; ++i)
	     	for(j = jmin; j <= jmax; ++j)
		for(k = kmin; k <= kmax; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    coord_theta = cell_center[index].m_coords[0];
                    coord_z = cell_center[index].m_coords[1];
                    coord_r = cell_center[index].m_coords[2];
                    //vel_theta = 0.0;
                    vel_theta = cell_center[index].m_state.m_U[0];
                    vel_z = cell_center[index].m_state.m_U[1];
                    vel_r = cell_center[index].m_state.m_U[2];
                    fprintf(outfile, "%.16g %.16g %16g\n",-vel_theta*sin(coord_theta)+vel_r*cos(coord_theta),vel_theta*cos(coord_theta)+vel_r*sin(coord_theta), vel_z);
                }
            }
            if(iFparams->movie_option->plot_pres)
            {
                fprintf(outfile, "SCALARS pressure double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
    		for(i = imin; i <= imax; ++i)
	     	for(j = jmin; j <= jmax; ++j)
		for(k = kmin; k <= kmax; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_P);
                } 
            }
            if(iFparams->movie_option->plot_vort)
            {
            }

	    if(true) // phase
	    {
		fprintf(outfile, "SCALARS phase double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
	
                for(i = imin; i <= imax; ++i)
                for(j = jmin; j <= jmax; ++j)
                for(k = kmin; k <= kmax; ++k)
		{

		    index = d_index3d(i,j,k,top_gmax);
		    fprintf(outfile,"%.16g\n", (double) cell_center[index].comp);
		}
	    }
	    if(true) // density
	    {
		fprintf(outfile, "SCALARS density double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
	
                for(i = imin; i <= imax; ++i)
                for(j = jmin; j <= jmax; ++j)
                for(k = kmin; k <= kmax; ++k)
		{

		    index = d_index3d(i,j,k,top_gmax);
		    fprintf(outfile,"%.16g\n", cell_center[index].m_state.m_rho);
		}
	    }
	    if(true) // pressure error
	    {
		fprintf(outfile, "SCALARS preserror double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
	
                for(i = imin; i <= imax; ++i)
                for(j = jmin; j <= jmax; ++j)
                for(k = kmin; k <= kmax; ++k)
		{

		    index = d_index3d(i,j,k,top_gmax);
		    double errorp = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);
		    fprintf(outfile,"%.16g\n", errorp);
		}
	    }

	    if(true)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
	     	int *ppgmax = front->pp_grid->gmax;
	    	int ppx = myid % ppgmax[0];
	    	int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
	    	int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "SCALARS divU double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
	
                for(i = imin; i <= imax; ++i)
                for(j = jmin; j <= jmax; ++j)
                for(k = kmin; k <= kmax; ++k)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
		    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	    	    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	    	    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);
		    
		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;
	    
		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;
	    
		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;
	    
		    if ( (ppz == 0 && k == kmin) )
			rur_nb[0] = 0.0;
		        //rur_nb[0] = iFparams->bvel[0][2] * r_nb[0];
		    else
		        rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
	    
		    if ( (ppz == (ppgmax[2]-1) && k == kmax) )
			rur_nb[1] = 0.0;
		        //rur_nb[1] = iFparams->bvel[1][2] * r_nb[1];
		    else
		        rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
	    
		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
		           (uz_nb[1] - uz_nb[0])/top_h[1] +
		           (rur_nb[1] - rur_nb[0])/(r*top_h[2]);
	    
                    fprintf(outfile,"%.16g\n",fabs(divu));

		}
	    }

	    if(true)//rho_t
	    {
		fprintf(outfile, "SCALARS rho_t double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
	
                for(i = imin; i <= imax; ++i)
                for(j = jmin; j <= jmax; ++j)
                for(k = kmin; k <= kmax; ++k)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    	    
		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

                    fprintf(outfile,"%.16g\n",fabs(rhot));

		}
	    }

	    if(true)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
	     	int *ppgmax = front->pp_grid->gmax;
	    	int ppx = myid % ppgmax[0];
	    	int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
	    	int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "SCALARS mass_conserv double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
	
                for(i = imin; i <= imax; ++i)
                for(j = jmin; j <= jmax; ++j)
                for(k = kmin; k <= kmax; ++k)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
		    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	    	    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	    	    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);
		    
		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;
	    
		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;
	    
		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;
	    
		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;
	    
		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;
	    
		    if ( (ppz == 0 && k == kmin) )
		    {
			rur_nb[0] = 0.0;
		        //rur_nb[0] = iFparams->bvel[0][2] * r_nb[0];
		        rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
		        rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
		        rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }
	    
		    if ( (ppz == (ppgmax[2]-1) && k == kmax) )
		    {
			rur_nb[1] = 0.0;
		        //rur_nb[1] = iFparams->bvel[1][2] * r_nb[1];
		        rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
		        rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
		        rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }
	    
		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
		           (uz_nb[1] - uz_nb[0])/top_h[1] +
		           (rur_nb[1] - rur_nb[0])/(r*top_h[2]);
	    
		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
		        	  (uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
		    	  (rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);
	    
		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;
	    
		    double mass_conserv = rhot + divrhou;
                    fprintf(outfile,"%.16g\n",fabs(mass_conserv));

		}
	    }
            fclose(outfile);
        }
}
*/


void Incompress_Solver_Smooth_3D_Cylindrical::printInteriorVelocity_vtu(char *out_name, bool binary)
{
    if (binary == YES)
	printInteriorVelocity_vtu_binary(out_name);
    else
	printInteriorVelocity_vtu_ascii(out_name);
}


void Incompress_Solver_Smooth_3D_Cylindrical::printInteriorVelocity_vtu_binary(char *out_name)
{
    	std::vector<int> ph_index;
        int   i,j,k,l,index;
	int   ii,jj,kk;
        double coord_theta,coord_z,coord_r;
        double vel_theta, vel_z, vel_r;
	int num_points;
	int num_cells;
	int num_cell_list;
	int pointsr, pointsz, pointstheta;
	int icoords[3];
	int p_gmax[3];
	    
	double thetamin, thetamax, zmin, zmax, rmin, rmax;    	

	char dirname1[256];
        char dirname2[256];
        char filename1[200];
	char filename2[200];
        char filename3[200];
        char filename4[200];

	char str[100];
	double val[3];
	int ival[8];
	unsigned int uival[2];

	FILE *outfile;

	if (pp_numnodes() > 1)
	    sprintf(dirname1,"%s/P-%s",out_name,
		    right_flush(pp_mynode(),4));
	else
	    sprintf(dirname1,"%s",out_name);
        sprintf(dirname1,"%s/phase_1",dirname1);
        if (!create_directory(dirname1,YES))
        {
            screen("Cannot create directory %s\n",dirname1);
            clean_up(ERROR);
        }
        sprintf(filename1,"%s/ph1-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename1,"%s-p%s",filename1,right_flush(pp_mynode(),4));
#endif
        sprintf(filename1,"%s.vtu",filename1);

        sprintf(filename3,"%s/ph1_debug-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename3,"%s-p%s",filename3,right_flush(pp_mynode(),4));
#endif
        sprintf(filename3,"%s.vtu",filename3);

        if (pp_numnodes() > 1)
	    sprintf(dirname2,"%s/P-%s",out_name,
		    right_flush(pp_mynode(),4));
	else
	    sprintf(dirname2,"%s",out_name);
        sprintf(dirname2,"%s/phase_2",dirname2);
        if (!create_directory(dirname2,YES))
        {
            screen("Cannot create directory %s\n",dirname2);
            clean_up(ERROR);
        }
        sprintf(filename2,"%s/ph2-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename2,"%s-p%s",filename2,right_flush(pp_mynode(),4));
#endif
        sprintf(filename2,"%s.vtu",filename2);

        sprintf(filename4,"%s/ph2_debug-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename4,"%s-p%s",filename4,right_flush(pp_mynode(),4));
#endif
        sprintf(filename4,"%s.vtu",filename4);

        if(iFparams->movie_option->plot_phase_1)
        {
	    /***** Phase One Visualization *****/

	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();

	    outfile = fopen(filename1,"wb");
	    
	    fprintf(outfile, "<\?xml version=\"1.0\"\?>\n");
	    if (hardware_is_little_endian())
		fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	    else
		fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	    fprintf(outfile, "  <UnstructuredGrid>\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;

	    fprintf(outfile, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

	    fprintf(outfile, "      <Points>\n");
	    fprintf(outfile, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">");
	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "      </Points>\n");

	    fprintf(outfile, "      <Cells>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">");
	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);

		ival[0] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[1] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[2] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[3] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[4] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[5] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[6] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		ival[7] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);
		fwrite(ival, sizeof(int), 8, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">");
	    for (i = 1; i <= (int)ph_index.size(); i++)
	    {
		ival[0] = 8*i;
		fwrite(ival, sizeof(int), 1, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">");
	    for (i = 1; i <= (int)ph_index.size(); i++)
	    {
		uival[0] = 11;
		fwrite(uival, sizeof(unsigned int), 1, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "      </Cells>\n");

	    fprintf(outfile, "      <CellData>\n");
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];
		    //vel_theta = 0.0;
		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];

		    val[0] = -vel_theta*sin(coord_theta)+vel_r*cos(coord_theta);
		    val[1] = vel_theta*cos(coord_theta)+vel_r*sin(coord_theta);
		    val[2] = vel_z;
		    fwrite(val, sizeof(double), 3, outfile);
		}
		fprintf(outfile,"\n");
		fprintf(outfile, "        </DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    val[0] = cell_center[index].m_state.m_P;
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"phase\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = 1.0;
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"density\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = cell_center[index].m_state.m_rho;
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }
	    fprintf(outfile, "      </CellData>\n");
	    fprintf(outfile, "    </Piece>\n");
	    fprintf(outfile, "  </UnstructuredGrid>\n");
	    fprintf(outfile, "</VTKFile>\n");
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename3,"wb");
	    fprintf(outfile, "<\?xml version=\"1.0\"\?>\n");
	    if (hardware_is_little_endian())
		fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	    else
		fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");

	    fprintf(outfile, "  <UnstructuredGrid>\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

	    fprintf(outfile, "      <Points>\n");
	    fprintf(outfile, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">");
	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "      </Points>\n");

	    fprintf(outfile, "      <Cells>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">");
	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);

		ival[0] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[1] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[2] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[3] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[4] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[5] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[6] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		ival[7] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);
		fwrite(ival, sizeof(int), 8, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">");
	    for (i = 1; i <= (int)ph_index.size(); i++)
	    {
		ival[0] = 8*i;
		fwrite(ival, sizeof(int), 1, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">");
	    for (i = 1; i <= (int)ph_index.size(); i++)
	    {
		uival[0] = 11;
		fwrite(uival, sizeof(unsigned int), 1, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "      </Cells>\n");


	    fprintf(outfile, "      <CellData>\n");
	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"preserror\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];


		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"divU\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    val[0] = fabs(divu);
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"rho_t\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    val[0] = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"mass_conserv\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;
		    val[0] = fabs(mass_conserv);
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }
	    fprintf(outfile, "      </CellData>\n");
	    fprintf(outfile, "    </Piece>\n");
	    fprintf(outfile, "  </UnstructuredGrid>\n");
	    fprintf(outfile, "</VTKFile>\n");

	    fclose(outfile);
	    }
	}

        if(iFparams->movie_option->plot_phase_2)
        {
	    /***** Phase One Visualization *****/

	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();

	    outfile = fopen(filename2,"wb");
	    
	    fprintf(outfile, "<\?xml version=\"1.0\"\?>\n");
	    if (hardware_is_little_endian())
		fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	    else
		fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	    fprintf(outfile, "  <UnstructuredGrid>\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;

	    fprintf(outfile, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

	    fprintf(outfile, "      <Points>\n");
	    fprintf(outfile, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">");
	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "      </Points>\n");

	    fprintf(outfile, "      <Cells>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">");
	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);

		ival[0] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[1] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[2] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[3] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[4] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[5] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[6] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		ival[7] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);
		fwrite(ival, sizeof(int), 8, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">");
	    for (i = 1; i <= (int)ph_index.size(); i++)
	    {
		ival[0] = 8*i;
		fwrite(ival, sizeof(int), 1, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">");
	    for (i = 1; i <= (int)ph_index.size(); i++)
	    {
		uival[0] = 11;
		fwrite(uival, sizeof(unsigned int), 1, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "      </Cells>\n");

	    fprintf(outfile, "      <CellData>\n");
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];
		    //vel_theta = 0.0;
		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];

		    val[0] = -vel_theta*sin(coord_theta)+vel_r*cos(coord_theta);
		    val[1] = vel_theta*cos(coord_theta)+vel_r*sin(coord_theta);
		    val[2] = vel_z;
		    fwrite(val, sizeof(double), 3, outfile);
		}
		fprintf(outfile,"\n");
		fprintf(outfile, "        </DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    val[0] = cell_center[index].m_state.m_P;
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"phase\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = 2.0;
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"density\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = cell_center[index].m_state.m_rho;
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }
	    fprintf(outfile, "      </CellData>\n");
	    fprintf(outfile, "    </Piece>\n");
	    fprintf(outfile, "  </UnstructuredGrid>\n");
	    fprintf(outfile, "</VTKFile>\n");
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename4,"wb");
	    fprintf(outfile, "<\?xml version=\"1.0\"\?>\n");
	    if (hardware_is_little_endian())
		fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	    else
		fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");

	    fprintf(outfile, "  <UnstructuredGrid>\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

	    fprintf(outfile, "      <Points>\n");
	    fprintf(outfile, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">");
	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "      </Points>\n");

	    fprintf(outfile, "      <Cells>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">");
	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);

		ival[0] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[1] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[2] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[3] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[4] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[5] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[6] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		ival[7] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);
		fwrite(ival, sizeof(int), 8, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">");
	    for (i = 1; i <= (int)ph_index.size(); i++)
	    {
		ival[0] = 8*i;
		fwrite(ival, sizeof(int), 1, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">");
	    for (i = 1; i <= (int)ph_index.size(); i++)
	    {
		uival[0] = 11;
		fwrite(uival, sizeof(unsigned int), 1, outfile);
	    }
	    //fprintf(outfile, "\n");
	    fprintf(outfile, "</DataArray>\n");
	    fprintf(outfile, "      </Cells>\n");


	    fprintf(outfile, "      <CellData>\n");
	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"preserror\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];


		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"divU\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    val[0] = fabs(divu);
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"rho_t\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    val[0] = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"mass_conserv\" NumberOfComponents=\"1\" format=\"binary\">");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;
		    val[0] = fabs(mass_conserv);
		    fwrite(val, sizeof(double), 1, outfile);
		}
		//fprintf(outfile, "\n");
		fprintf(outfile, "</DataArray>\n");
	    }
	    fprintf(outfile, "      </CellData>\n");
	    fprintf(outfile, "    </Piece>\n");
	    fprintf(outfile, "  </UnstructuredGrid>\n");
	    fprintf(outfile, "</VTKFile>\n");

	    fclose(outfile);
	    }
	}
}
	
void Incompress_Solver_Smooth_3D_Cylindrical::printInteriorVelocity_vtu_ascii(char *out_name)
{
    	std::vector<int> ph_index;
        int   i,j,k,l,index;
	int   ii,jj,kk;
        double coord_theta,coord_z,coord_r;
        double vel_theta, vel_z, vel_r;
	int num_points;
	int num_cells;
	int num_cell_list;
	int pointsr, pointsz, pointstheta;
	int icoords[3];
	int p_gmax[3];
	    
	double thetamin, thetamax, zmin, zmax, rmin, rmax;    	

	char dirname1[256];
        char dirname2[256];
        char filename1[200];
	char filename2[200];
        char filename3[200];
        char filename4[200];

	FILE *outfile;

	if (pp_numnodes() > 1)
	    sprintf(dirname1,"%s/P-%s",out_name,
		    right_flush(pp_mynode(),4));
	else
	    sprintf(dirname1,"%s",out_name);
        sprintf(dirname1,"%s/phase_1",dirname1);
        if (!create_directory(dirname1,YES))
        {
            screen("Cannot create directory %s\n",dirname1);
            clean_up(ERROR);
        }
        sprintf(filename1,"%s/ph1-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename1,"%s-p%s",filename1,right_flush(pp_mynode(),4));
#endif
        sprintf(filename1,"%s.vtu",filename1);

        sprintf(filename3,"%s/ph1_debug-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename3,"%s-p%s",filename3,right_flush(pp_mynode(),4));
#endif
        sprintf(filename3,"%s.vtu",filename3);

        if (pp_numnodes() > 1)
	    sprintf(dirname2,"%s/P-%s",out_name,
		    right_flush(pp_mynode(),4));
	else
	    sprintf(dirname2,"%s",out_name);
        sprintf(dirname2,"%s/phase_2",dirname2);
        if (!create_directory(dirname2,YES))
        {
            screen("Cannot create directory %s\n",dirname2);
            clean_up(ERROR);
        }
        sprintf(filename2,"%s/ph2-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename2,"%s-p%s",filename2,right_flush(pp_mynode(),4));
#endif
        sprintf(filename2,"%s.vtu",filename2);

        sprintf(filename4,"%s/ph2_debug-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename4,"%s-p%s",filename4,right_flush(pp_mynode(),4));
#endif
        sprintf(filename4,"%s.vtu",filename4);

        if(iFparams->movie_option->plot_phase_1)
        {
	    /***** Phase One Visualization *****/

	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();

	    outfile = fopen(filename1,"w");
	    
	    fprintf(outfile, "<\?xml version=\"1.0\"\?>\n");
	    fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	    fprintf(outfile, "  <UnstructuredGrid>\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;

	    fprintf(outfile, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

	    fprintf(outfile, "      <Points>\n");
	    fprintf(outfile, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
	    }
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "      </Points>\n");

	    fprintf(outfile, "      <Cells>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		int index0 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		int index1 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		int index2 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		int index3 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		int index4 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		int index5 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		int index6 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		int index7 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fprintf(outfile, "%i %i %i %i %i %i %i %i\n",
			index0,index1,index2,index3,index4,index5,index6,index7);
	    }
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
	    for (i = 1; i <= (int)ph_index.size(); i++)
		fprintf(outfile, "%d\n", 8*i);
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
	    for (i = 1; i <= (int)ph_index.size(); i++)
		fprintf(outfile, "11\n");
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "      </Cells>\n");

	    fprintf(outfile, "      <CellData>\n");
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];
		    //vel_theta = 0.0;
		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];
		    fprintf(outfile, "%.16g %.16g %16g\n",-vel_theta*sin(coord_theta)+vel_r*cos(coord_theta),vel_theta*cos(coord_theta)+vel_r*sin(coord_theta), vel_z);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_P);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"phase\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    if(cell_center[index].comp == LIQUID_COMP1)
			fprintf(outfile,"%.16g\n", 1.0);
		    if(cell_center[index].comp == LIQUID_COMP2)
			fprintf(outfile,"%.16g\n", 2.0);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    fprintf(outfile,"%.16g\n", cell_center[index].m_state.m_rho);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    fprintf(outfile, "      </CellData>\n");
	    fprintf(outfile, "    </Piece>\n");
	    fprintf(outfile, "  </UnstructuredGrid>\n");
	    fprintf(outfile, "</VTKFile>\n");
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename3,"w");
	    fprintf(outfile, "<\?xml version=\"1.0\"\?>\n");
	    fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	    fprintf(outfile, "  <UnstructuredGrid>\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

	    fprintf(outfile, "      <Points>\n");
	    fprintf(outfile, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
	    }
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "      </Points>\n");

	    fprintf(outfile, "      <Cells>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		int index0 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		int index1 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		int index2 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		int index3 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		int index4 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		int index5 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		int index6 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		int index7 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fprintf(outfile, "%i %i %i %i %i %i %i %i\n",
			index0,index1,index2,index3,index4,index5,index6,index7);
	    }
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
	    for (i = 1; i <= (int)ph_index.size(); i++)
		fprintf(outfile, "%d\n", 8*i);
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
	    for (i = 1; i <= (int)ph_index.size(); i++)
		fprintf(outfile, "11\n");
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "      </Cells>\n");

	    fprintf(outfile, "      <CellData>\n");
	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"preserror\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    double errorp = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);
		    fprintf(outfile,"%.16g\n", errorp);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];


		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"divU\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    fprintf(outfile,"%.16g\n",fabs(divu));

		}
		fprintf(outfile, "        </DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"rho_t\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    fprintf(outfile,"%.16g\n",fabs(rhot));

		}
		fprintf(outfile, "        </DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"mass_conserv\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;
		    fprintf(outfile,"%.16g\n",fabs(mass_conserv));

		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    fprintf(outfile, "      </CellData>\n");
	    fprintf(outfile, "    </Piece>\n");
	    fprintf(outfile, "  </UnstructuredGrid>\n");
	    fprintf(outfile, "</VTKFile>\n");

	    fclose(outfile);
	    }
	}

        if(iFparams->movie_option->plot_phase_2)
        {
	    /***** Phase Two Visualization *****/

	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();

	    outfile = fopen(filename2,"w");
	    
	    fprintf(outfile, "<\?xml version=\"1.0\"\?>\n");
	    fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	    fprintf(outfile, "  <UnstructuredGrid>\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;

	    fprintf(outfile, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

	    fprintf(outfile, "      <Points>\n");
	    fprintf(outfile, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
	    }
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "      </Points>\n");

	    fprintf(outfile, "      <Cells>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		int index0 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		int index1 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		int index2 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		int index3 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		int index4 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		int index5 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		int index6 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		int index7 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fprintf(outfile, "%i %i %i %i %i %i %i %i\n",
			index0,index1,index2,index3,index4,index5,index6,index7);
	    }
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
	    for (i = 1; i <= (int)ph_index.size(); i++)
		fprintf(outfile, "%d\n", 8*i);
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
	    for (i = 1; i <= (int)ph_index.size(); i++)
		fprintf(outfile, "11\n");
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "      </Cells>\n");

	    fprintf(outfile, "      <CellData>\n");
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];
		    //vel_theta = 0.0;
		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];
		    fprintf(outfile, "%.16g %.16g %16g\n",-vel_theta*sin(coord_theta)+vel_r*cos(coord_theta),vel_theta*cos(coord_theta)+vel_r*sin(coord_theta), vel_z);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_P);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"phase\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    if(cell_center[index].comp == LIQUID_COMP1)
			fprintf(outfile,"%.16g\n", 1.0);
		    if(cell_center[index].comp == LIQUID_COMP2)
			fprintf(outfile,"%.16g\n", 2.0);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    fprintf(outfile,"%.16g\n", cell_center[index].m_state.m_rho);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    fprintf(outfile, "      </CellData>\n");
	    fprintf(outfile, "    </Piece>\n");
	    fprintf(outfile, "  </UnstructuredGrid>\n");
	    fprintf(outfile, "</VTKFile>\n");
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename4,"w");
	    fprintf(outfile, "<\?xml version=\"1.0\"\?>\n");
	    fprintf(outfile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	    fprintf(outfile, "  <UnstructuredGrid>\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

	    fprintf(outfile, "      <Points>\n");
	    fprintf(outfile, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
	    }
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "      </Points>\n");

	    fprintf(outfile, "      <Cells>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		int index0 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		int index1 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		int index2 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		int index3 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		int index4 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		int index5 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		int index6 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		int index7 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fprintf(outfile, "%i %i %i %i %i %i %i %i\n",
			index0,index1,index2,index3,index4,index5,index6,index7);
	    }
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
	    for (i = 1; i <= (int)ph_index.size(); i++)
		fprintf(outfile, "%d\n", 8*i);
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
	    for (i = 1; i <= (int)ph_index.size(); i++)
		fprintf(outfile, "11\n");
	    fprintf(outfile, "        </DataArray>\n");
	    fprintf(outfile, "      </Cells>\n");

	    fprintf(outfile, "      <CellData>\n");
	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"preserror\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    double errorp = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);
		    fprintf(outfile,"%.16g\n", errorp);
		}
		fprintf(outfile, "        </DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];


		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"divU\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    fprintf(outfile,"%.16g\n",fabs(divu));

		}
		fprintf(outfile, "        </DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"rho_t\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    fprintf(outfile,"%.16g\n",fabs(rhot));

		}
		fprintf(outfile, "        </DataArray>\n");
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "        <DataArray type=\"Float64\" Name=\"mass_conserv\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;
		    fprintf(outfile,"%.16g\n",fabs(mass_conserv));

		}
		fprintf(outfile, "        </DataArray>\n");
	    }
	    fprintf(outfile, "      </CellData>\n");
	    fprintf(outfile, "    </Piece>\n");
	    fprintf(outfile, "  </UnstructuredGrid>\n");
	    fprintf(outfile, "</VTKFile>\n");
	    fclose(outfile);
	    
	    }
	}
}
	

void Incompress_Solver_Smooth_3D_Cylindrical::printInteriorVelocity(char *out_name, bool binary)
{
    if(binary == YES)
    {
	if(hardware_is_little_endian())
	    printInteriorVelocity_little_endian(out_name);
	else
	    printInteriorVelocity_big_endian(out_name);
    }
    else
	printInteriorVelocity_ascii(out_name);
}

void Incompress_Solver_Smooth_3D_Cylindrical::printInteriorVelocity_ascii(char *out_name)
{
    	std::vector<int> ph_index;
        int   i,j,k,l,index;
	int   ii,jj,kk;
        double coord_theta,coord_z,coord_r;
        double vel_theta, vel_z, vel_r;
	int num_points;
	int num_cells;
	int num_cell_list;
	int pointsr, pointsz, pointstheta;
	int icoords[3];
	int p_gmax[3];
	    
	double thetamin, thetamax, zmin, zmax, rmin, rmax;
        char dirname1[256];
        char dirname2[256];
        char filename1[200];
	char filename2[200];
        char filename3[200];
        char filename4[200];

	FILE *outfile;

	if (pp_numnodes() > 1)
	    sprintf(dirname1,"%s/P-%s",out_name,
		    right_flush(pp_mynode(),4));
	else
	    sprintf(dirname1,"%s",out_name);
        sprintf(dirname1,"%s/phase_1",dirname1);
        if (!create_directory(dirname1,YES))
        {
            screen("Cannot create directory %s\n",dirname1);
            clean_up(ERROR);
        }
        sprintf(filename1,"%s/ph1-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename1,"%s-p%s",filename1,right_flush(pp_mynode(),4));
#endif
        sprintf(filename1,"%s.vtk",filename1);

        sprintf(filename3,"%s/ph1_debug-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename3,"%s-p%s",filename3,right_flush(pp_mynode(),4));
#endif
        sprintf(filename3,"%s.vtk",filename3);

        if (pp_numnodes() > 1)
	    sprintf(dirname2,"%s/P-%s",out_name,
		    right_flush(pp_mynode(),4));
	else
	    sprintf(dirname2,"%s",out_name);
        sprintf(dirname2,"%s/phase_2",dirname2);
        if (!create_directory(dirname2,YES))
        {
            screen("Cannot create directory %s\n",dirname2);
            clean_up(ERROR);
        }
        sprintf(filename2,"%s/ph2-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename2,"%s-p%s",filename2,right_flush(pp_mynode(),4));
#endif
        sprintf(filename2,"%s.vtk",filename2);

        sprintf(filename4,"%s/ph2_debug-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename4,"%s-p%s",filename4,right_flush(pp_mynode(),4));
#endif
        sprintf(filename4,"%s.vtk",filename4);

        if(iFparams->movie_option->plot_phase_1)
        {
	    /***** Phase One Visualization *****/

	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();
	    outfile = fopen(filename1,"w");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase one\n");
	    fprintf(outfile,"ASCII\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
	    }

	    fprintf(outfile, "CELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		int index0 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		int index1 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		int index2 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		int index3 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		int index4 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		int index5 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		int index6 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		int index7 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fprintf(outfile, "8 %i %i %i %i %i %i %i %i\n",
			index0,index1,index2,index3,index4,index5,index6,index7);
	    }

	    fprintf(outfile, "CELL_TYPES %i\n", num_cells);
	    for (i = 0; i < (int)ph_index.size();i++)
		fprintf(outfile, "11\n");

	    fprintf(outfile, "CELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "VECTORS velocity double\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];
		    //vel_theta = 0.0;
		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];
		    fprintf(outfile, "%.16g %.16g %16g\n",-vel_theta*sin(coord_theta)+vel_r*cos(coord_theta),vel_theta*cos(coord_theta)+vel_r*sin(coord_theta), vel_z);
		}
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "SCALARS pressure double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_P);
		} 
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "SCALARS phase double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    if(cell_center[index].comp == LIQUID_COMP1)
			fprintf(outfile,"%.16g\n", 1.0);
		    if(cell_center[index].comp == LIQUID_COMP2)
			fprintf(outfile,"%.16g\n", 2.0);
		}
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "SCALARS density double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    fprintf(outfile,"%.16g\n", cell_center[index].m_state.m_rho);
		}
	    }
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename3,"w");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase one\n");
	    fprintf(outfile,"ASCII\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
	    }

	    fprintf(outfile, "CELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		int index0 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		int index1 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		int index2 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		int index3 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		int index4 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		int index5 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		int index6 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		int index7 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fprintf(outfile, "8 %i %i %i %i %i %i %i %i\n",
			index0,index1,index2,index3,index4,index5,index6,index7);
	    }

	    fprintf(outfile, "CELL_TYPES %i\n", num_cells);
	    for (i = 0; i < (int)ph_index.size();i++)
		fprintf(outfile, "11\n");

	    fprintf(outfile, "CELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "SCALARS preserror double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    double errorp = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);
		    fprintf(outfile,"%.16g\n", errorp);
		}
	    }
	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "SCALARS divU double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    fprintf(outfile,"%.16g\n",fabs(divu));

		}
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "SCALARS rho_t double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    fprintf(outfile,"%.16g\n",fabs(rhot));

		}
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "SCALARS MB_residual double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;
		    fprintf(outfile,"%.16g\n",fabs(mass_conserv));

		}
	    }
	    fclose(outfile);
	    }
	}

        if(iFparams->movie_option->plot_phase_2)
        {
	    /***** Phase Two Visualization *****/
	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();
	    outfile = fopen(filename2,"w");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase two\n");
	    fprintf(outfile,"ASCII\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
	    }

	    fprintf(outfile, "CELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		int index0 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		int index1 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		int index2 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		int index3 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		int index4 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		int index5 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		int index6 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		int index7 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fprintf(outfile, "8 %i %i %i %i %i %i %i %i\n",
			index0,index1,index2,index3,index4,index5,index6,index7);
	    }

	    fprintf(outfile, "CELL_TYPES %i\n", num_cells);
	    for (i = 0; i < (int)ph_index.size();i++)
		fprintf(outfile, "11\n");

	    fprintf(outfile, "CELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "VECTORS velocity double\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];
		    //vel_theta = 0.0;
		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];
		    fprintf(outfile, "%.16g %.16g %16g\n",-vel_theta*sin(coord_theta)+vel_r*cos(coord_theta),vel_theta*cos(coord_theta)+vel_r*sin(coord_theta), vel_z);
		}
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "SCALARS pressure double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_P);
		} 
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "SCALARS phase double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    if(cell_center[index].comp == LIQUID_COMP1)
			fprintf(outfile,"%.16g\n", 1.0);
		    if(cell_center[index].comp == LIQUID_COMP2)
			fprintf(outfile,"%.16g\n", 2.0);
		}
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "SCALARS density double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    fprintf(outfile,"%.16g\n", cell_center[index].m_state.m_rho);
		}
	    }
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename4,"w");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase two\n");
	    fprintf(outfile,"ASCII\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }
	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
	    }

	    fprintf(outfile, "CELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		int index0 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		int index1 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		int index2 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		int index3 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		int index4 = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		int index5 = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		int index6 = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		int index7 = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fprintf(outfile, "8 %i %i %i %i %i %i %i %i\n",
			index0,index1,index2,index3,index4,index5,index6,index7);
	    }

	    fprintf(outfile, "CELL_TYPES %i\n", num_cells);
	    for (i = 0; i < (int)ph_index.size();i++)
		fprintf(outfile, "11\n");

	    fprintf(outfile, "CELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "SCALARS preserror double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    double errorp = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);
		    fprintf(outfile,"%.16g\n", errorp);
		}
	    }

	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "SCALARS divU double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    fprintf(outfile,"%.16g\n",fabs(divu));

		}
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "SCALARS rho_t double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    fprintf(outfile,"%.16g\n",fabs(rhot));

		}
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "SCALARS MB_residual double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);


		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;
		    fprintf(outfile,"%.16g\n",fabs(mass_conserv));
		}
	    }
	    fclose(outfile);
	    }
	}
}

void Incompress_Solver_Smooth_3D_Cylindrical::printInteriorVelocity_big_endian(char* out_name)
{
    	std::vector<int> ph_index;
        int   i,j,k,l,index;
	int   ii,jj,kk;
        double coord_theta,coord_z,coord_r;
        double vel_theta, vel_z, vel_r;
	int num_points;
	int num_cells;
	int num_cell_list;
	int pointsr, pointsz, pointstheta;
	int icoords[3];
	int p_gmax[3];
	    
	double thetamin, thetamax, zmin, zmax, rmin, rmax;

        char dirname1[256];
        char dirname2[256];
        char filename1[200];
	char filename2[200];
        char filename3[200];
        char filename4[200];

	char str[100];
	double val[3];
	int ival[8];

	FILE *outfile;

	if (pp_numnodes() > 1)
        sprintf(dirname1,"%s/P-%s",out_name,
                right_flush(pp_mynode(),4));
	else
	    sprintf(dirname1, "%s", out_name);
        sprintf(dirname1,"%s/phase_1",dirname1);
        if (!create_directory(dirname1,YES))
        {
            screen("Cannot create directory %s\n",dirname1);
            clean_up(ERROR);
        }
        sprintf(filename1,"%s/ph1-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename1,"%s-p%s",filename1,right_flush(pp_mynode(),4));
#endif
        sprintf(filename1,"%s.vtk",filename1);

        sprintf(filename3,"%s/ph1_debug-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename3,"%s-p%s",filename3,right_flush(pp_mynode(),4));
#endif
        sprintf(filename3,"%s.vtk",filename3);

	if(pp_numnodes() > 1)
        sprintf(dirname2,"%s/P-%s",out_name,
                right_flush(pp_mynode(),4));
	else
	    sprintf(dirname2,"%s",out_name);
        sprintf(dirname2,"%s/phase_2",dirname2);
        if (!create_directory(dirname2,YES))
        {
            screen("Cannot create directory %s\n",dirname2);
            clean_up(ERROR);
        }
        sprintf(filename2,"%s/ph2-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename2,"%s-p%s",filename2,right_flush(pp_mynode(),4));
#endif
        sprintf(filename2,"%s.vtk",filename2);

        sprintf(filename4,"%s/ph2_debug-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename4,"%s-p%s",filename4,right_flush(pp_mynode(),4));
#endif
        sprintf(filename4,"%s.vtk",filename4);

        if(iFparams->movie_option->plot_phase_1)
        {
	    /***** Phase One Visualization *****/
	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();
	    outfile = fopen(filename1,"wb");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase one\n");
	    fprintf(outfile,"BINARY\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile);
	    }

	    fprintf(outfile, "\nCELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		ival[0] = 8;
		fwrite(ival, sizeof(int), 1, outfile);
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		ival[0] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[1] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[2] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[3] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[4] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[5] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[6] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		ival[7] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fwrite(ival, sizeof(int), 8, outfile);
	    }

	    fprintf(outfile, "\nCELL_TYPES %i\n", num_cells);
	    for (i = 0; i < (int)ph_index.size();i++)
	    {
		ival[0] = 11;
		fwrite(ival, sizeof(int), 1, outfile);
	    }

	    fprintf(outfile, "\nCELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "VECTORS velocity double\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];

		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];
		    
		    val[0] = -vel_theta*sin(coord_theta)+vel_r*cos(coord_theta);
		    val[1] = vel_theta*cos(coord_theta)+vel_r*sin(coord_theta);
		    val[2] = vel_z;

		    fwrite(val, sizeof(double), 3, outfile);
		}
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "\nSCALARS pressure double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    val[0] = cell_center[index].m_state.m_P;
		    fwrite(val, sizeof(double), 1, outfile);
		} 
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "\nSCALARS phase double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    val[0] = 1.0;
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "\nSCALARS density double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = cell_center[index].m_state.m_rho;
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename3,"wb");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase one\n");
	    fprintf(outfile,"BINARY\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile);	    
	    }

	    fprintf(outfile, "\nCELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		ival[0] = 8;
		fwrite(ival, sizeof(int), 1, outfile);
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		ival[0] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[1] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[2] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[3] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[4] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[5] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[6] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		ival[7] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fwrite(ival, sizeof(int), 8, outfile);
	    }

	    fprintf(outfile, "\nCELL_TYPES %i\n", num_cells);

	    for (i = 0; i < (int)ph_index.size();i++)
	    {
		ival[0] = 11;
		fwrite(ival, sizeof(int), 1, outfile);
	    }

	    fprintf(outfile, "\nCELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "SCALARS preserror double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }

	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "\nSCALARS divU double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    val[0] = divu;
		    fwrite(val, sizeof(double), 1, outfile);

		}
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "\nSCALARS rho_t double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;
		    val[0] = rhot;
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "\nSCALARS MB_residual double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;

		    val[0] = mass_conserv;
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    fclose(outfile);
	    }
	}

        if(iFparams->movie_option->plot_phase_2)
        {

	    /***** Phase Two Visualization *****/
	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();
	    outfile = fopen(filename2,"wb");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase two\n");
	    fprintf(outfile,"BINARY\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile);
	    }

	    fprintf(outfile, "\nCELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		ival[0] = 8;
		fwrite(ival, sizeof(int), 1, outfile);
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		ival[0] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[1] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[2] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[3] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[4] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[5] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[6] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		ival[7] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fwrite(ival, sizeof(int), 8, outfile);
	    }

	    fprintf(outfile, "\nCELL_TYPES %i\n", num_cells);
	    for (i = 0; i < (int)ph_index.size();i++)
	    {
		ival[0] = 11;
		fwrite(ival, sizeof(int), 1, outfile);
	    }

	    fprintf(outfile, "\nCELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "VECTORS velocity double\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];

		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];
		    val[0] = -vel_theta*sin(coord_theta)+vel_r*cos(coord_theta);
		    val[1] = vel_theta*cos(coord_theta)+vel_r*sin(coord_theta);
		    val[2] = vel_z;

		    fwrite(val, sizeof(double), 3, outfile);		}
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "\nSCALARS pressure double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    val[0] = cell_center[index].m_state.m_P;
		    fwrite(val, sizeof(double), 1, outfile);
		} 
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "\nSCALARS phase double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    val[0] = 2.0;
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "\nSCALARS density double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = cell_center[index].m_state.m_rho;
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename4,"wb");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase two\n");
	    fprintf(outfile,"BINARY\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = coord_r*cos(coord_theta);
		val[1] = coord_r*sin(coord_theta);
		val[2] = coord_z;
		fwrite(val, sizeof(double), 3, outfile);
	    }

	    fprintf(outfile, "\nCELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		ival[0] = 8;
		fwrite(ival, sizeof(int), 1, outfile);
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		ival[0] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[1] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax);
		ival[2] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[3] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax);
		ival[4] = d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[5] = d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax);
		ival[6] = d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax);
		ival[7] = d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax);

		fwrite(ival, sizeof(int), 8, outfile);
	    }

	    fprintf(outfile, "\nCELL_TYPES %i\n", num_cells);

	    for (i = 0; i < (int)ph_index.size();i++)
	    {
		ival[0] = 11;
		fwrite(ival, sizeof(int), 1, outfile);
	    }

	    fprintf(outfile, "\nCELL_DATA %i\n", num_cells);

	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "SCALARS preserror double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }

	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "\nSCALARS divU double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    val[0] = divu;
		    fwrite(val, sizeof(double), 1, outfile);

		}
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "\nSCALARS rho_t double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;
		    val[0] = rhot;
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "\nSCALARS MB_residual double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;

		    val[0] = mass_conserv;
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    fclose(outfile);
	    }
	}
}


void Incompress_Solver_Smooth_3D_Cylindrical::printInteriorVelocity_little_endian(char* out_name)
{
    	std::vector<int> ph_index;
        int   i,j,k,l,index;
	int   ii,jj,kk;
        double coord_theta,coord_z,coord_r;
        double vel_theta, vel_z, vel_r;
	int num_points;
	int num_cells;
	int num_cell_list;
	int pointsr, pointsz, pointstheta;
	int icoords[3];
	int p_gmax[3];
	    
	double thetamin, thetamax, zmin, zmax, rmin, rmax;

        char dirname1[256];
        char dirname2[256];
        char filename1[200];
	char filename2[200];
        char filename3[200];
        char filename4[200];

	char str[100];
	double val[3];
	int ival[8];

	FILE *outfile;

	if (pp_numnodes() > 1)
        sprintf(dirname1,"%s/P-%s",out_name,
                right_flush(pp_mynode(),4));
	else
	    sprintf(dirname1,"%s",out_name);
        sprintf(dirname1,"%s/phase_1",dirname1);
        if (!create_directory(dirname1,YES))
        {
            screen("Cannot create directory %s\n",dirname1);
            clean_up(ERROR);
        }
        sprintf(filename1,"%s/ph1-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename1,"%s-p%s",filename1,right_flush(pp_mynode(),4));
#endif
        sprintf(filename1,"%s.vtk",filename1);

        sprintf(filename3,"%s/ph1_debug-t%s",dirname1,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename3,"%s-p%s",filename3,right_flush(pp_mynode(),4));
#endif
        sprintf(filename3,"%s.vtk",filename3);

	if (pp_numnodes() > 1)
        sprintf(dirname2,"%s/P-%s",out_name,
                right_flush(pp_mynode(),4));
	else
	    sprintf(dirname2,"%s",out_name);
        sprintf(dirname2,"%s/phase_2",dirname2);
        if (!create_directory(dirname2,YES))
        {
            screen("Cannot create directory %s\n",dirname2);
            clean_up(ERROR);
        }
        sprintf(filename2,"%s/ph2-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename2,"%s-p%s",filename2,right_flush(pp_mynode(),4));
#endif
        sprintf(filename2,"%s.vtk",filename2);

        sprintf(filename4,"%s/ph2_debug-t%s",dirname2,right_flush(front->step,7));

#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename4,"%s-p%s",filename4,right_flush(pp_mynode(),4));
#endif
        sprintf(filename4,"%s.vtk",filename4);

        if(iFparams->movie_option->plot_phase_1)
        {
	    /***** Phase One Visualization *****/
	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();
	    outfile = fopen(filename1,"wb");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase one\n");
	    fprintf(outfile,"BINARY\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile);
	    }

	    fprintf(outfile, "\nCELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		ival[0] = endian_int_swap(8);
		fwrite(ival, sizeof(int), 1, outfile);
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		ival[0] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax));
		ival[1] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax));
		ival[2] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax));
		ival[3] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax));
		ival[4] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax));
		ival[5] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax));
		ival[6] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax));
		ival[7] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax));


		fwrite(ival, sizeof(int), 8, outfile);
	    }

	    fprintf(outfile, "\nCELL_TYPES %i\n", num_cells);
	    for (i = 0; i < (int)ph_index.size();i++)
	    {
		ival[0] = endian_int_swap(11);
		fwrite(ival, sizeof(int), 1, outfile);
	    }

	    fprintf(outfile, "\nCELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "VECTORS velocity double\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];

		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];
		    
		    val[0] = endian_double_swap(-vel_theta*sin(coord_theta)+vel_r*cos(coord_theta));
		    val[1] = endian_double_swap(vel_theta*cos(coord_theta)+vel_r*sin(coord_theta));
		    val[2] = endian_double_swap(vel_z);

		    fwrite(val, sizeof(double), 3, outfile);
		}
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "\nSCALARS pressure double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    val[0] = endian_double_swap(cell_center[index].m_state.m_P);
		    fwrite(val, sizeof(double), 1, outfile);
		} 
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "\nSCALARS phase double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    val[0] = endian_double_swap(1.0);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "\nSCALARS density double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = endian_double_swap(cell_center[index].m_state.m_rho);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename3,"wb");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase one\n");
	    fprintf(outfile,"BINARY\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP1)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile);	    
	    }

	    fprintf(outfile, "\nCELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		ival[0] = endian_int_swap(8);
		fwrite(ival,sizeof(int),1,outfile);
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		ival[0] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax));
		ival[1] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax));
		ival[2] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax));
		ival[3] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax));
		ival[4] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax));
		ival[5] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax));
		ival[6] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax));
		ival[7] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax));

		fwrite(ival, sizeof(int), 8, outfile);
	    }

	    fprintf(outfile, "\nCELL_TYPES %i\n", num_cells);

	    for (i = 0; i < (int)ph_index.size();i++)
	    {
		ival[0] = endian_int_swap(11);
		fwrite(ival, sizeof(int), 1, outfile);
	    }

	    fprintf(outfile, "\nCELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "SCALARS preserror double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = endian_double_swap(fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP));
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }

	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "\nSCALARS divU double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    val[0] = endian_double_swap(divu);
		    fwrite(val, sizeof(double), 1, outfile);

		}
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "\nSCALARS rho_t double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;
		    val[0] = endian_double_swap(rhot);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "\nSCALARS MB_residual double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;

		    val[0] = endian_double_swap(mass_conserv);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    fclose(outfile);
	    }
	}

        if(iFparams->movie_option->plot_phase_2)
        {

	    /***** Phase Two Visualization *****/
	    if (iFparams->movie_option->plot_velo ||
		iFparams->movie_option->plot_pres ||
		iFparams->movie_option->plot_comp ||
		iFparams->movie_option->plot_dens)
	    {
	    ph_index.clear();
	    outfile = fopen(filename2,"wb");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase two\n");
	    fprintf(outfile,"BINARY\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile);
	    }

	    fprintf(outfile, "\nCELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		ival[0] = endian_int_swap(8);
		fwrite(ival, sizeof(int), 1, outfile);
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		ival[0] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax));
		ival[1] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax));
		ival[2] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax));
		ival[3] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax));
		ival[4] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax));
		ival[5] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax));
		ival[6] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax));
		ival[7] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax));

		fwrite(ival, sizeof(int), 8, outfile);
	    }

	    fprintf(outfile, "\nCELL_TYPES %i\n", num_cells);
	    for (i = 0; i < (int)ph_index.size();i++)
	    {
		ival[0] = endian_int_swap(11);
		fwrite(ival, sizeof(int), 1, outfile);
	    }

	    fprintf(outfile, "\nCELL_DATA %i\n", num_cells);
	    if(iFparams->movie_option->plot_velo)
	    {
		fprintf(outfile, "VECTORS velocity double\n");
		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    coord_theta = cell_center[index].m_coords[0];
		    coord_z = cell_center[index].m_coords[1];
		    coord_r = cell_center[index].m_coords[2];

		    vel_theta = cell_center[index].m_state.m_U[0];
		    vel_z = cell_center[index].m_state.m_U[1];
		    vel_r = cell_center[index].m_state.m_U[2];
		    val[0] = endian_double_swap(-vel_theta*sin(coord_theta)+vel_r*cos(coord_theta));
		    val[1] = endian_double_swap(vel_theta*cos(coord_theta)+vel_r*sin(coord_theta));
		    val[2] = endian_double_swap(vel_z);

		    fwrite(val, sizeof(double), 3, outfile);		}
	    }
	    if(iFparams->movie_option->plot_pres)
	    {
		fprintf(outfile, "\nSCALARS pressure double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");
		for(i = 0; i < (int)ph_index.size(); i++) 
		{
		    index = ph_index[i];
		    val[0] = endian_double_swap(cell_center[index].m_state.m_P);
		    fwrite(val, sizeof(double), 1, outfile);
		} 
	    }
	    if(iFparams->movie_option->plot_vort)
	    {
	    }

	    if(iFparams->movie_option->plot_comp) // phase
	    {
		fprintf(outfile, "\nSCALARS phase double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    val[0] = endian_double_swap(2.0);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    if(iFparams->movie_option->plot_dens) // density
	    {
		fprintf(outfile, "\nSCALARS density double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = endian_double_swap(cell_center[index].m_state.m_rho);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    fclose(outfile);
	    }

	    //debug
	    if(iFparams->movie_option->plot_pres_error ||
	       iFparams->movie_option->plot_divU ||
	       iFparams->movie_option->plot_rhot ||
	       iFparams->movie_option->plot_MB_residual)
	    {
	    ph_index.clear();
	    outfile = fopen(filename4,"wb");
	    fprintf(outfile,"# vtk DataFile Version 3.0\n");
	    fprintf(outfile,"States of the whole computational domain, phase two\n");
	    fprintf(outfile,"BINARY\n");
	    fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (cell_center[index].comp == LIQUID_COMP2)
		    ph_index.push_back(index);
	    }

	    num_cells = (int) ph_index.size();
	    num_cell_list = 9*num_cells;

	    pointsr = kmax - kmin + 2;
	    pointsz = jmax - jmin + 2;
	    pointstheta = imax - imin + 2;
	    num_points = pointsr*pointsz*pointstheta;

	    p_gmax[0] = pointstheta - 1;
	    p_gmax[1] = pointsz - 1;
	    p_gmax[2] = pointsr - 1;


	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    thetamin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    zmin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    rmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(imax,jmax,kmax,top_gmax);
	    thetamax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    zmax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    rmax = cell_center[index].m_coords[2] + top_h[2]/2.0;


	    fprintf(outfile, "POINTS %d double\n", num_points);

	    for (k = 0; k < pointsr; k++)
	    for (j = 0; j < pointsz; j++)
	    for (i = 0; i < pointstheta; i++)
	    {
		coord_theta = thetamin + i*top_h[0];
		coord_z = zmin + j*top_h[1];
		coord_r = rmin + k*top_h[2];
		val[0] = endian_double_swap(coord_r*cos(coord_theta));
		val[1] = endian_double_swap(coord_r*sin(coord_theta));
		val[2] = endian_double_swap(coord_z);
		fwrite(val, sizeof(double), 3, outfile);
	    }

	    fprintf(outfile, "\nCELLS %i %i\n", num_cells, num_cell_list);

	    for (i = 0; i < (int)ph_index.size(); i++)
	    {
		ival[0] = endian_int_swap(8);
		fwrite(ival, sizeof(int), 1, outfile);
		index = ph_index[i];
		getIcoords(index, icoords, top_gmax);
		ival[0] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-4, p_gmax));
		ival[1] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-4, p_gmax));
		ival[2] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-4, p_gmax));
		ival[3] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-4, p_gmax));
		ival[4] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-4, icoords[2]-3, p_gmax));
		ival[5] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-4, icoords[2]-3, p_gmax));
		ival[6] = endian_int_swap(d_index3d(icoords[0]-4, icoords[1]-3, icoords[2]-3, p_gmax));
		ival[7] = endian_int_swap(d_index3d(icoords[0]-3, icoords[1]-3, icoords[2]-3, p_gmax));

		fwrite(ival, sizeof(int), 8, outfile);
	    }

	    fprintf(outfile, "\nCELL_TYPES %i\n", num_cells);

	    for (i = 0; i < (int)ph_index.size();i++)
	    {
		ival[0] = endian_int_swap(11);
		fwrite(ival, sizeof(int), 1, outfile);
	    }

	    fprintf(outfile, "\nCELL_DATA %i\n", num_cells);

	    if(iFparams->movie_option->plot_pres_error) // pressure error
	    {
		fprintf(outfile, "SCALARS preserror double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    val[0] = endian_double_swap(fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP));
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }

	    if(iFparams->movie_option->plot_divU)//DivU
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "\nSCALARS divU double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
			rur_nb[0] = 0.0;
		    else
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
			rur_nb[1] = 0.0;
		    else
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    val[0] = endian_double_swap(divu);
		    fwrite(val, sizeof(double), 1, outfile);

		}
	    }

	    if(iFparams->movie_option->plot_rhot)//rho_t
	    {
		fprintf(outfile, "\nSCALARS rho_t double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;
		    val[0] = endian_double_swap(rhot);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }

	    if(iFparams->movie_option->plot_MB_residual)//Mass Conservation
	    {
		const int nn = pp_numnodes();
		int myid = pp_mynode();
		int *ppgmax = front->pp_grid->gmax;
		int ppx = myid % ppgmax[0];
		int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
		int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

		double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];
		double rur_nb[2], utheta_nb[2], uz_nb[2];
		int index_nb[6];

		fprintf(outfile, "\nSCALARS MB_residual double\n");
		fprintf(outfile, "LOOKUP_TABLE default\n");

		for(i = 0; i < (int)ph_index.size(); i++)
		{
		    index = ph_index[i];
		    getIcoords(index, icoords, top_gmax);
		    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
		    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
		    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
		    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
		    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
		    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

		    double r, r_nb[2];
		    r = cell_center[index].m_coords[2];
		    r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
		    r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

		    utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
		    utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

		    rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
		    rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

		    uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
		    uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

		    rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
		    rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

		    if ( (ppz == 0) && (icoords[2] == kmin) )
		    {
			rur_nb[0] = 0.0;
			rho_r_nb[0] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
			rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
		    }

		    if ( (ppz == (ppgmax[2]-1)) && (icoords[2] == kmax) )
		    {
			rur_nb[1] = 0.0;
			rho_r_nb[1] = cell_center[index].m_state.m_rho;
		    }
		    else
		    {
			rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
			rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
		    }

		    double divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
			(uz_nb[1] - uz_nb[0])/top_h[1] +
			(rur_nb[1] - rur_nb[0])/(r*top_h[2]);

		    double divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
			(uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
			(rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

		    double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

		    double mass_conserv = rhot + divrhou;

		    val[0] = endian_double_swap(mass_conserv);
		    fwrite(val, sizeof(double), 1, outfile);
		}
	    }
	    fclose(outfile);
	    }
	}
}

double Incompress_Solver_Smooth_3D_Cylindrical::getCellVolume(int *icoords)
{

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double r = cell_center[index].m_coords[2];
    double cellVolume = r*top_h[0]*top_h[1]*top_h[2];

    return cellVolume;
}

double Incompress_Solver_Smooth_3D_Cylindrical::getFaceArea(int *icoords, GRID_DIRECTION dir)
{

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double r = cell_center[index].m_coords[2];

    double faceArea;

    switch(dir)
    {
    case WEST:
    case EAST:
	faceArea = top_h[1]*top_h[2];
	break;
    case SOUTH:
    case NORTH:
	faceArea = r*top_h[2]*top_h[0];
	break;
    case LOWER:
	faceArea = (r-top_h[2]/2.0)*top_h[0]*top_h[1];
	break;
    case UPPER:
	faceArea = (r+top_h[2]/2.0)*top_h[0]*top_h[1];
	break;
    default:
	assert(false);
    }
    return faceArea;

}

void Incompress_Solver_Smooth_3D_Cylindrical::getFaceCenter(int *icoords, GRID_DIRECTION dir, double faceCenter[3])
{
    double dx = 0, dy = 0,dz = 0;

    switch(dir)
    {
    case WEST:
	dx = -top_h[0];
	break;
    case EAST:
	dx =  top_h[0];
	break;
    case SOUTH:
	dy = -top_h[1];
	break;
    case NORTH:
	dy =  top_h[1];
	break;
    case LOWER:
	dz = -top_h[2];
	break;
    case UPPER:
	dz =  top_h[2];
	break;
    default:
	assert(false);
    }

    faceCenter[0] = dx/2 + cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_coords[0];
    faceCenter[1] = dy/2 + cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_coords[1];
    faceCenter[2] = dz/2 + cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_coords[2];

}



void Incompress_Solver_Smooth_3D_Cylindrical::computeProjection_Shuqiang(void)
{
//    Incompress_Solver_Smooth_3D_Cylindrical::computeProjection();
//    return;

    int index, index_nb[6], size;
    double rhs, coeff[6], rho[6], rho0;
    int I,I_nb[6];
    int i,j,k,l,icoords[MAXD];
    INTERFACE *intfc = front->interf;
    double P_max,P_min;
    int icrds_Pmax[MAXD],icrds_Pmin[MAXD];
    COMPONENT comp;
    double aII,rho_nb[6];
    double coords[MAXD],crx_coords[MAXD];
    double **vel = iFparams->field->vel;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

    max_value = 0.0;
    double value;
    double sum_div;
    sum_div = 0.0;
    int num_iter = 0;
    double rel_residual = 0.0;

    PETSc solver;
    solver.Create(ilower, iupper-1, 7, 7);
    solver.Reset_A();
    solver.Reset_b();
    solver.Reset_x();
    size = iupper - ilower;

    setIndexMap();

    for (l = 0; l < dim; ++l)
	for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
		for (i = 0; i <= top_gmax[0]; i++)
		{
		    index  = d_index3d(i,j,k,top_gmax);
		    vel[l][index] = cell_center[index].m_state.m_U[l];
		}

    /* Compute velocity divergence */
    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		index  = d_index3d(i,j,k,top_gmax);
		array[index] = computeFieldPointDiv(icoords,vel);
	    }
    scatMeshArray();
    for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		index  = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.div_U = array[index];
	    }


    if(debugging("step_size"))
    {
	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    value = fabs(cell_center[index].m_state.div_U);
		    sum_div = sum_div + cell_center[index].m_state.div_U * cell_center[index].m_coords[2] * top_h[0]*top_h[1]*top_h[2];
		    if(value > max_value)
			max_value = value;
		}
	pp_global_sum(&sum_div,1);
	printf("\nThe summation of divergence of U is %.16g\n",sum_div);
	pp_global_max(&max_value,1);
	printf("\nThe max value of divergence of U is %.16g\n",max_value);
	max_value = 0.0;
    }

//    double normal[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
//    double flux[6];

    double dh[3] = {top_h[0],top_h[1], top_h[2]};

    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index  = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		I = ijk_to_I[i][j][k];
		if (I == -1) continue;

		index_nb[0] = d_index3d(i-1,j,k,top_gmax);
		index_nb[1] = d_index3d(i+1,j,k,top_gmax);
		index_nb[2] = d_index3d(i,j-1,k,top_gmax);
		index_nb[3] = d_index3d(i,j+1,k,top_gmax);
		index_nb[4] = d_index3d(i,j,k-1,top_gmax);
		index_nb[5] = d_index3d(i,j,k+1,top_gmax);
		I_nb[0] = ijk_to_I[i-1][j][k];
		I_nb[1] = ijk_to_I[i+1][j][k];
		I_nb[2] = ijk_to_I[i][j-1][k];
		I_nb[3] = ijk_to_I[i][j+1][k];
		I_nb[4] = ijk_to_I[i][j][k-1];
		I_nb[5] = ijk_to_I[i][j][k+1];
		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;

		rho0   = cell_center[index].m_state.m_rho;
		for (l = 0; l < 6; ++l)
		{
		    if (I_nb[l] == -1)
			index_nb[l] = index;
		    rho[l] = 1.0/2*(rho0 + cell_center[index_nb[l]].m_state.m_rho);
		    //coeff[l] = 1.0/rho[l]/sqr(top_h[l/2]);
		}

//		double rr = sqr(cell_center[index].m_coords[2]);
		double r = cell_center[index].m_coords[2];
		double cellVolume = getCellVolume(icoords);
		double faceArea;

		rhs = 0;
		// WEST, theta
		faceArea = getFaceArea(icoords,WEST);
		if(I_nb[0]>=0)
		{
		    solver.Add_A(I, I,       -1/r*1/dh[0]*faceArea /rho[0]);
		    solver.Add_A(I, I_nb[0],  1/r*1/dh[0]*faceArea /rho[0]);
		}
		else
		    ;
		// EAST
		faceArea = getFaceArea(icoords,EAST);
		if(I_nb[1]>=0)
		{
		    solver.Add_A(I, I,       -1/r*1/dh[0]*faceArea /rho[1]);
		    solver.Add_A(I, I_nb[1],  1/r*1/dh[0]*faceArea /rho[1]);
		}
		else
		    ;
		// SOUTH, z
		faceArea = getFaceArea(icoords,SOUTH);
		if(I_nb[2]>=0)
		{
		    solver.Add_A(I, I,       -1  *1/dh[1]*faceArea /rho[2]);
		    solver.Add_A(I, I_nb[2],  1  *1/dh[1]*faceArea /rho[2]);
		}
		else
		    ;
		// NORTH
		faceArea = getFaceArea(icoords,NORTH);
		if(I_nb[3]>=0)
		{
		    solver.Add_A(I, I,       -1  *1/dh[1]*faceArea /rho[3]);
		    solver.Add_A(I, I_nb[3],  1  *1/dh[1]*faceArea /rho[3]);
		}
		else
		    ;
		// LOWER, r
		faceArea = getFaceArea(icoords,LOWER);
		if(I_nb[4]>=0)
		{
		    solver.Add_A(I, I,       -1  *1/dh[2]*faceArea /rho[4]);
		    solver.Add_A(I, I_nb[4],  1  *1/dh[2]*faceArea /rho[4]);
		}
		else
		    ;
		// UPPER
		faceArea = getFaceArea(icoords,UPPER);
		if(I_nb[5]>=0)
		{
		    solver.Add_A(I, I,       -1  *1/dh[2]*faceArea /rho[5]);
		    solver.Add_A(I, I_nb[5],  1  *1/dh[2]*faceArea /rho[5]);
		}
		else
		    ;
		rhs = cell_center[index].m_state.div_U/accum_dt * cellVolume;
		solver.Add_b(I, rhs);


	    }

    solver.SetMaxIter(40000);
    solver.SetTol(1e-10);

    start_clock("Before Petsc Solver in Projection step");
    solver.Solve_withPureNeumann();

//    solver.Print_A(NULL);
//    solver.Print_b(NULL);
//    exit(0);


    solver.GetNumIterations(&num_iter);
    solver.GetFinalRelativeResidualNorm(&rel_residual);
    if(rel_residual > 1)
    {
	printf("\n The solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	solver.Reset_x();
	solver.Solve_withPureNeumann_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
    }
    stop_clock("After Petsc Solver in Projection step");

    double *x;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
    solver.Get_x(x);

    if (debugging("PETSc"))
	(void) printf("Incompress_Solver_Smooth_3D_Cylindrical::"
		"computeProjection: "
		"num_iter = %d, rel_residual = %le \n",
		num_iter, rel_residual);

    P_max = -HUGE;		P_min = HUGE;
    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		I = ijk_to_I[i][j][k];
		array[index] = x[I-ilower];
	    }
    scatMeshArray();
    for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		index  = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.m_phi = array[index];
	    }

    if(debugging("step_size"))
    {
	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    value = fabs(cell_center[index].m_state.m_phi);
		    if (value > max_value)
			max_value = value;
		}
	pp_global_max(&max_value,1);
	printf("\nThe max value of phi is %.16g\n",max_value);
    }

    FT_FreeThese(1,x);
}

void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled_Shuqiang(void)
{
    COMPONENT comp;
    int index,index_nb[6],size;
    int I,I_nb[6];
    double coords[MAXD],crx_coords[MAXD];
    double coeff[18],mu,rho; //,U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
    L_STATE source_term, U_nb[6], U_nb_new[6], U_center, rhs;

    int i,j,k,l,nb,icoords[MAXD];
    INTERFACE *intfc = front->interf;
    double speed;
    double *x;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
    POINTER intfc_state;
    HYPER_SURF *hs;
    int num_iter;
    double rel_residual;
    double r, rr;

    max_speed = 0.0;
    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

    PETSc solver;
    solver.Create(3*ilower, 3*iupper-1, 10, 10);
    // 7theta + 2r  for the first equation
    // 7z for the second equation
    // 7r + 2theta for the third equation

    solver.Reset_A();
    solver.Reset_b();
    solver.Reset_x();

    printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

    double dtheta, dz, dr; //, faceCenter[3];
    dtheta = top_h[0];
    dz	   = top_h[1];
    dr     = top_h[2];

    double dh[]  = {dtheta,dz,dr};

    static int debug_index = 0;

    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		I  = ijk_to_I[i][j][k];
		if (I == -1) continue;

		index  = d_index3d(i,j,k,top_gmax);
		//6 neighbours of the center cell
		index_nb[0] = d_index3d(i-1,j,k,top_gmax);
		index_nb[1] = d_index3d(i+1,j,k,top_gmax);
		index_nb[2] = d_index3d(i,j-1,k,top_gmax);
		index_nb[3] = d_index3d(i,j+1,k,top_gmax);
		index_nb[4] = d_index3d(i,j,k-1,top_gmax);
		index_nb[5] = d_index3d(i,j,k+1,top_gmax);


		//6 neighbours of the center cell
		I_nb[0] = ijk_to_I[i-1][j][k];
		I_nb[1] = ijk_to_I[i+1][j][k];
		I_nb[2] = ijk_to_I[i][j-1][k];
		I_nb[3] = ijk_to_I[i][j+1][k];
		I_nb[4] = ijk_to_I[i][j][k-1];
		I_nb[5] = ijk_to_I[i][j][k+1];


		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		comp = top_comp[index];

		mu = cell_center[index].m_state.m_mu;
		rho = cell_center[index].m_state.m_rho;
		U_center = cell_center[index].m_state;

		for (nb = 0; nb < 6; nb++)
		{
		    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			    comp,&intfc_state,&hs,crx_coords,m_t_old) &&
			    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			// old boundary condition
			U_nb[nb].m_U[0] = getStateVel[0](intfc_state);
			U_nb[nb].m_U[1] = getStateVel[1](intfc_state);
			U_nb[nb].m_U[2] = getStateVel[2](intfc_state);

			FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
				comp,&intfc_state,&hs,crx_coords,m_t_new);
			U_nb_new[nb].m_U[0] = getStateVel[0](intfc_state);
			U_nb_new[nb].m_U[1] = getStateVel[1](intfc_state);
			U_nb_new[nb].m_U[2] = getStateVel[2](intfc_state);

			assert(I_nb[nb]<0);
		    }
		    else
		    {
			U_nb[nb] = cell_center[index_nb[nb]].m_state;
		    }
		}

		double r = cell_center[index].m_coords[2];
		double cellVolume = getCellVolume(icoords);
		double faceArea;

		for(nb = 0; nb<6; nb++)
		{
		    faceArea = getFaceArea(icoords,dir[nb]);

		    if(I_nb[nb]>=0)
		    {
			compDiffWithSmoothProperty_cellFace(
				&solver,icoords,I,I_nb[nb],
				dir[nb],
				dh[nb/2],
				faceArea,cellVolume,r,mu,rho,
				U_nb[nb],U_nb_new[nb],U_center);
		    }
		    else
			compDiffWithSmoothProperty_Dirichlet(
				&solver,icoords,I,I_nb[nb],
				dir[nb],
				dh[nb/2],
				faceArea,cellVolume,r,mu,rho,
				U_nb[nb],U_nb_new[nb],U_center);
		}

		compDiffWithSmoothProperty_cellInterior(
			&solver,icoords,I,cellVolume,r,mu,rho,U_center);

	    }

    solver.SetMaxIter(40000);
    solver.SetTol(1e-10);

    start_clock("Before Petsc Solve");
    solver.Solve_GMRES();
    solver.GetNumIterations(&num_iter);
    solver.GetFinalRelativeResidualNorm(&rel_residual);
    /*
	if (rel_residual > 1)
	{
	    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	    solver.Reset_x();
	    solver.Solve_GMRES();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	}
     */
    stop_clock("After Petsc Solve");

    // get back the solution
    solver.Get_x(x);

    if (debugging("PETSc"))
	(void) printf("Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled_Shuqiang: "
		"num_iter = %d, rel_residual = %le. \n",
		num_iter,rel_residual);


    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		I = ijk_to_I[i][j][k];
		index = d_index3d(i,j,k,top_gmax);
		if (I >= 0)
		{
		    cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
		    cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		    cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
		    speed = fabs(cell_center[index].m_state.m_U[0]) +
			    fabs(cell_center[index].m_state.m_U[1]) +
			    fabs(cell_center[index].m_state.m_U[2]);
		    if (speed > max_speed)
			max_speed = speed;
		}
		else
		{
		    cell_center[index].m_state.m_U[0] = 0.0;
		    cell_center[index].m_state.m_U[1] = 0.0;
		    cell_center[index].m_state.m_U[2] = 0.0;
		}
	    }
    for (l = 0; l < 3; ++l)
    {
	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)

		{
		    index  = d_index3d(i,j,k,top_gmax);
		    array[index] = cell_center[index].m_state.m_U[l];
		}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
		for (i = 0; i <= top_gmax[0]; i++)

		{
		    index  = d_index3d(i,j,k,top_gmax);
		    cell_center[index].m_state.m_U[l] = array[index];
		}
    }
    pp_global_max(&max_speed,1);

    FT_FreeThese(1,x);
}       /* end compDiffWithSmoothProperty3d */

void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_cellFace(
	PETSc *pSolver,
	int *icoords,
	int I, int I_nb,
	GRID_DIRECTION dir,
	double dh,
	double faceArea,
	double cellVolume,
	double r,
	double mu,
	double rho,
	L_STATE &U_nb,
	L_STATE &U_nb_new,
	L_STATE &U_center)
{
    double mymu = 1.0/2*m_dt*mu/rho;

    double coeff = 1;
    if(dir==WEST || dir==EAST)
	coeff *= 1.0/r;

    // used for first order derivative with theta
    int sign = 0;
    if(dir==WEST)
	sign = -1;
    else if(dir==EAST)
	sign = 1;

    // 1st equ
    // u^{*}
    pSolver->Add_A(I*3, I*3,     - coeff/dh*faceArea*(-mymu));
    pSolver->Add_A(I*3, I_nb*3,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3,          - coeff/dh*faceArea*( mymu) * U_center.m_U[0]);
    pSolver->Add_b(I*3,            coeff/dh*faceArea*( mymu) * U_nb.m_U[0]);

    if(sign!=0)
    {
	pSolver->Add_A(I*3, I*3+2,    -sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_A(I*3, I_nb*3+2,  sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_b(I*3,           -sign/(r*r)/dh*cellVolume*( mymu) * U_center.m_U[2]);
	pSolver->Add_b(I*3,            sign/(r*r)/dh*cellVolume*( mymu) * U_nb.m_U[2]);
    }

    // 2nd equ
    pSolver->Add_A(I*3+1, I*3+1,     - coeff/dh*faceArea*(-mymu));
    pSolver->Add_A(I*3+1, I_nb*3+1,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3+1,            - coeff/dh*faceArea*( mymu) * U_center.m_U[1]);
    pSolver->Add_b(I*3+1,              coeff/dh*faceArea*( mymu) * U_nb.m_U[1]);

    // 3rd equ
    pSolver->Add_A(I*3+2, I*3+2,     - coeff/dh*faceArea*(-mymu));
    pSolver->Add_A(I*3+2, I_nb*3+2,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3+2,            - coeff/dh*faceArea*( mymu) * U_center.m_U[2]);
    pSolver->Add_b(I*3+2,              coeff/dh*faceArea*( mymu) * U_nb.m_U[2]);

    if(sign!=0)
    {
	sign = -sign;	// this makes the following code similar with 1st equ.
	pSolver->Add_A(I*3+2, I*3,    -sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_A(I*3+2, I_nb*3,  sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_b(I*3+2,         -sign/(r*r)/dh*cellVolume*( mymu) * U_center.m_U[0]);
	pSolver->Add_b(I*3+2,          sign/(r*r)/dh*cellVolume*( mymu) * U_nb.m_U[0]);
    }
}
/**
 * TODO: modify the following code for boundary.
 */
void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_Dirichlet(
	PETSc *pSolver,
	int *icoords,
	int I, int I_nb,
	GRID_DIRECTION dir,
	double dh,
	double faceArea,
	double cellVolume,
	double r,
	double mu,
	double rho,
	L_STATE &U_nb,
	L_STATE &U_nb_new,
	L_STATE &U_center)
{
    double mymu = 1.0/2*m_dt*mu/rho;

    double coeff = 1;
    if(dir==WEST || dir==EAST)
	coeff *= 1.0/r;

    // used for first order derivative with theta
    int sign = 0;
    if(dir==WEST)
	sign = -1;
    else if(dir==EAST)
	sign = 1;

    //-----------------------------------------------------
    // 			1st equ
    //-----------------------------------------------------
    // u^{*}
    pSolver->Add_A(I*3, I*3,     - coeff/dh*faceArea*(-mymu) * 2);
//  pSolver->Add_A(I*3, I_nb*3,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3,            coeff/dh*faceArea*( mymu) * 2*U_nb_new.m_U[0]);
    pSolver->Add_b(I*3,          - coeff/dh*faceArea*( mymu) * 2*U_center.m_U[0]);
    pSolver->Add_b(I*3,            coeff/dh*faceArea*( mymu) * 2*U_nb.m_U[0]);

    if(sign!=0)
    {
	pSolver->Add_A(I*3, I*3+2,    -sign/(r*r)/dh*cellVolume*(-mymu) * 2);
//	pSolver->Add_A(I*3, I_nb*3+2,  sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_b(I*3,            sign/(r*r)/dh*cellVolume*( mymu) * 2*U_nb_new.m_U[2]);
	pSolver->Add_b(I*3,           -sign/(r*r)/dh*cellVolume*( mymu) * 2*U_center.m_U[2]);
	pSolver->Add_b(I*3,            sign/(r*r)/dh*cellVolume*( mymu) * 2*U_nb.m_U[2]);
    }

    //-----------------------------------------------------
    // 			2nd equ
    //-----------------------------------------------------
    pSolver->Add_A(I*3+1, I*3+1,     - coeff/dh*faceArea*(-mymu) * 2);
//  pSolver->Add_A(I*3+1, I_nb*3+1,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3+1, 	       coeff/dh*faceArea*( mymu) * 2*U_nb_new.m_U[1]);
    pSolver->Add_b(I*3+1,            - coeff/dh*faceArea*( mymu) * 2*U_center.m_U[1]);
    pSolver->Add_b(I*3+1,              coeff/dh*faceArea*( mymu) * 2*U_nb.m_U[1]);

    //-----------------------------------------------------
    // 			3rd equ
    //-----------------------------------------------------
    pSolver->Add_A(I*3+2, I*3+2,     - coeff/dh*faceArea*(-mymu) * 2);
//  pSolver->Add_A(I*3+2, I_nb*3+2,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3+2,              coeff/dh*faceArea*( mymu) * 2*U_nb_new.m_U[2]);
    pSolver->Add_b(I*3+2,            - coeff/dh*faceArea*( mymu) * 2*U_center.m_U[2]);
    pSolver->Add_b(I*3+2,              coeff/dh*faceArea*( mymu) * 2*U_nb.m_U[2]);

    if(sign!=0)
    {
	sign = -sign;	// this makes the following code similar with 1st equ.
	pSolver->Add_A(I*3+2, I*3,    -sign/(r*r)/dh*cellVolume*(-mymu) * 2);
//	pSolver->Add_A(I*3+2, I_nb*3,  sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_b(I*3+2,          sign/(r*r)/dh*cellVolume*( mymu) * 2*U_nb_new.m_U[0]);
	pSolver->Add_b(I*3+2,         -sign/(r*r)/dh*cellVolume*( mymu) * 2*U_center.m_U[0]);
	pSolver->Add_b(I*3+2,          sign/(r*r)/dh*cellVolume*( mymu) * 2*U_nb.m_U[0]);
    }
}
void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_cellInterior(
	PETSc *pSolver,
	int *icoords,
	int I,
	double cellVolume,
	double r,
	double mu,
	double rho,
	L_STATE &U_center)
{
    int index  = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double coords[3], rhs;
    double mymu = 1.0/2*m_dt*mu/rho;


    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm(coords, source_term);


    //-----------------------------------------------------
    // 			1st equ
    //-----------------------------------------------------
    pSolver->Add_A(3*I, 3*I, 1*cellVolume);
    pSolver->Add_b(3*I,      1*cellVolume * U_center.m_U[0]);
    pSolver->Add_A(3*I, 3*I, 1/(r*r)*cellVolume*( mymu));
    pSolver->Add_b(3*I,      1/(r*r)*cellVolume*(-mymu) * U_center.m_U[0]);

    rhs  = m_dt*source_term.m_U[0]*cellVolume;
    rhs += m_dt*cell_center[index].m_state.f_surf[0]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.grad_q[0]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.m_adv[0] *cellVolume;
    pSolver->Add_b(3*I, rhs);

    //-----------------------------------------------------
    // 			2nd equ
    //-----------------------------------------------------
    pSolver->Add_A(3*I+1, 3*I+1, 1*cellVolume);
    pSolver->Add_b(3*I+1,        1*cellVolume * U_center.m_U[1]);

    rhs  = m_dt*source_term.m_U[1]*cellVolume;
    rhs += m_dt*cell_center[index].m_state.f_surf[1]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.grad_q[1]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.m_adv[1] *cellVolume;
    pSolver->Add_b(3*I+1, rhs);

    //-----------------------------------------------------
    // 			3rd equ
    //-----------------------------------------------------
    pSolver->Add_A(3*I+2, 3*I+2, 1*cellVolume);
    pSolver->Add_b(3*I+2,        1*cellVolume * U_center.m_U[2]);
    pSolver->Add_A(3*I+2, 3*I+2, 1/(r*r)*cellVolume*( mymu));
    pSolver->Add_b(3*I+2,        1/(r*r)*cellVolume*(-mymu) * U_center.m_U[2]);
    rhs  = m_dt*source_term.m_U[2]*cellVolume;
    rhs += m_dt*cell_center[index].m_state.f_surf[2]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.grad_q[2]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.m_adv[2] *cellVolume;
    pSolver->Add_b(3*I+2, rhs);
}


void Incompress_Solver_Smooth_3D_Cylindrical::computeSubgridModel(void)
{
        int i,j,k,l;
        int index,index0,index1,index2,index3,index4,index5,index6;
        int size;
        int index000,index100,index010,index001,index101,index110,index011,index111;
        L_STATE state;
        double *u, *v, *w;
        double ulx,uly,ulz,vlx,vly,vlz,wlx,wly,wlz;
        double urx,ury,urz,vrx,vry,vrz,wrx,wry,wrz;
        double ux,uy,uz,vx,vy,vz,wx,wy,wz;
        double *s, *s11, *s12, *s13, *s22, *s23, *s33;
        double *ss11, *ss12, *ss13, *ss22, *ss23, *ss33;
        double *vel_u, *vel_v, *vel_w, *vel_uu, *vel_vv, *vel_ww, *vel_uv, *vel_vw, *vel_uw;
        double sum_vel_u,sum_vel_v,sum_vel_w,sum_vel_uu,sum_vel_vv,sum_vel_ww,sum_vel_uv,sum_vel_vw,sum_vel_uw;
        double sum_s11,sum_s12,sum_s13,sum_s22,sum_s23,sum_s33,sum_ss11,sum_ss12,sum_ss13,sum_ss22,sum_ss23,sum_ss33,sum_s;
        double ma11, ma12, ma13, ma22, ma23, ma33, la11, la12, la13, la22, la23, la33;
        double cs, r;
        int    ii,jj,kk,iii,jjj,kkk,iiii,jjjj,kkkk;
        int    NB = 2, numr = 1;
        int    NBC = pow(NB,3);
        //int    numk = (2.0-1.0)/(2*top_h[2]);  //1.0: GL[2] 2.0: GU[2] 
        int    numk = front->gmax[2]/(numr*NB), numkk;   //32: mesh size of r
        //double *deno, *nume; 
        double gdeno, gnume;
        double ***ldeno, ***lnume;
        //int indexb[4][4][4];
        double delta2, tdelta2;
        delta2 = sqr(pow((top_h[0]*top_h[1]*top_h[2]),(1.0/3.0)));
        tdelta2 = sqr(pow((NB*top_h[0]*NB*top_h[1]*NB*top_h[2]),(1.0/3.0)));
 
        size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
        FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s13,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s23,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s33,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss13,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss23,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss33,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_w,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uu,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_vv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_ww,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uw,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_vw,size,sizeof(double));
        //FT_VectorMemoryAlloc((POINTER*)&deno,numk,sizeof(double));
        //FT_VectorMemoryAlloc((POINTER*)&nume,numk,sizeof(double));
        FT_TriArrayMemoryAlloc((POINTER*)&ldeno,(imax-imin+1)/NB,(jmax-jmin+1)/NB,(kmax-kmin+1)/NB,sizeof(double));
        FT_TriArrayMemoryAlloc((POINTER*)&lnume,(imax-imin+1)/NB,(jmax-jmin+1)/NB,(kmax-kmin+1)/NB,sizeof(double));

        const int  nn = pp_numnodes();
        int        myid = pp_mynode();
        int   *ppgmax = front->pp_grid->gmax;
        int   ppx = myid % ppgmax[0];
        int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
        int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);


        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            u[index] = cell_center[index].m_state.m_U[0];
            v[index] = cell_center[index].m_state.m_U[1];
            w[index] = cell_center[index].m_state.m_U[2];

            if ( (ppz == 0 && k == kmin-1) )
            {
                u[index] = iFparams->bvel[0][0];        
                v[index] = iFparams->bvel[0][1];
                w[index] = iFparams->bvel[0][2];
            }
            if ( (ppz == (ppgmax[2]-1) && k == kmax+1) )
            {
                u[index] = iFparams->bvel[1][0];
                v[index] = iFparams->bvel[1][1];
                w[index] = iFparams->bvel[1][2];
            }
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index0  = d_index3d(i,j,k,top_gmax);
            r = cell_center[index0].m_coords[2];
            index1  = d_index3d(i-1,j,k,top_gmax);
            ulx = u[index1];
            vlx = v[index1];
            wlx = w[index1];
            index2  = d_index3d(i+1,j,k,top_gmax);
            urx = u[index2];
            vrx = v[index2];
            wrx = w[index2];
            index3  = d_index3d(i,j-1,k,top_gmax);
            uly = u[index3];
            vly = v[index3];
            wly = w[index3];
            index4  = d_index3d(i,j+1,k,top_gmax);
            ury = u[index4];
            vry = v[index4];
            wry = w[index4];
            index5  = d_index3d(i,j,k-1,top_gmax);
            ulz = u[index5];
            vlz = v[index5];
            wlz = w[index5];
            index6  = d_index3d(i,j,k+1,top_gmax);
            urz = u[index6];
            vrz = v[index6];
            wrz = w[index6];

            ux = (urx - ulx) / (2.0*top_h[0]);
            uy = (ury - uly) / (2.0*top_h[1]);
            uz = (urz - ulz) / (2.0*top_h[2]);
            vx = (vrx - vlx) / (2.0*top_h[0]);
            vy = (vry - vly) / (2.0*top_h[1]);
            vz = (vrz - vlz) / (2.0*top_h[2]);
            wx = (wrx - wlx) / (2.0*top_h[0]);
            wy = (wry - wly) / (2.0*top_h[1]);
            wz = (wrz - wlz) / (2.0*top_h[2]);
            s11[index0] = (1.0/r)*(ux+w[index0]);
            s12[index0] = (uy + ((1.0/r)*vx))/2.0;
            s13[index0] = (((1.0/r)*wx) + uz - (u[index0]/r))/2.0;
            s22[index0] = vy;
            s23[index0] = (wy + vz)/2.0;
            s33[index0] = wz;
            s[index0] = sqrt(2*( sqr(s11[index0]) + sqr(s12[index0])
                         + sqr(s13[index0]) + sqr(s22[index0])
                         + sqr(s23[index0]) + sqr(s33[index0]) 
                         + sqr(s12[index0]) + sqr(s13[index0])
                         + sqr(s23[index0]) ));

            ss11[index0] = s[index0]*s11[index0];
            ss12[index0] = s[index0]*s12[index0];
            ss13[index0] = s[index0]*s13[index0];
            ss22[index0] = s[index0]*s22[index0];
            ss23[index0] = s[index0]*s23[index0];
            ss33[index0] = s[index0]*s33[index0];

            vel_u[index0] = u[index0];
            vel_v[index0] = v[index0];
            vel_w[index0] = w[index0];
            vel_uu[index0] = u[index0]*u[index0];
            vel_vv[index0] = v[index0]*v[index0];
            vel_ww[index0] = w[index0]*w[index0];
            vel_uv[index0] = u[index0]*v[index0];
            vel_vw[index0] = v[index0]*w[index0];
            vel_uw[index0] = u[index0]*w[index0];
        }

        /*
        gdeno = 0.0;
        gnume = 0.0;
        for (k = 0; k < numk; k++)
        {
            deno[k] = 0.0;
            nume[k] = 0.0;
        }
        */

        for (k = 0; k <= ((kmax-kmin+1)/NB)-1; k++)
        for (j = 0; j <= ((jmax-jmin+1)/NB)-1; j++)
        for (i = 0; i <= ((imax-imin+1)/NB)-1; i++)
        {
            kk = (NB*k)+kmin;
            jj = (NB*j)+jmin;
            ii = (NB*i)+imin;
            index000 = d_index3d(ii,jj,kk,top_gmax);

                sum_vel_u = sum_vel_v = sum_vel_w = 0.0;
                sum_vel_uu = sum_vel_vv = sum_vel_ww = 0.0;
                sum_vel_uv = sum_vel_vw = sum_vel_uw = 0.0;
                sum_s11 = sum_s12 = sum_s13 = sum_s22 = sum_s23 = sum_s33 = 0.0;
                sum_ss11 = sum_ss12 = sum_ss13 = sum_ss22 = sum_ss23 = sum_ss33 = 0.0;
                sum_s = 0.0;
                for(kkk = kk; kkk < kk+NB; kkk++)
                for(jjj = jj; jjj < jj+NB; jjj++)
                for(iii = ii; iii < ii+NB; iii++)
                {
                    index0  = d_index3d(iii,jjj,kkk,top_gmax);
                    sum_vel_u += vel_u[index0];
                    sum_vel_v += vel_v[index0];
                    sum_vel_w += vel_w[index0];
                    sum_vel_uu += vel_uu[index0];
                    sum_vel_vv += vel_vv[index0];
                    sum_vel_ww += vel_ww[index0];
                    sum_vel_uv += vel_uv[index0];
                    sum_vel_uw += vel_uw[index0];
                    sum_vel_vw += vel_vw[index0];
                    sum_s11 += s11[index0];
                    sum_s12 += s12[index0];
                    sum_s13 += s13[index0];
                    sum_s22 += s22[index0];
                    sum_s23 += s23[index0];
                    sum_s33 += s33[index0];
                    sum_ss11 += ss11[index0];
                    sum_ss12 += ss12[index0];
                    sum_ss13 += ss13[index0];
                    sum_ss22 += ss22[index0];
                    sum_ss23 += ss23[index0];
                    sum_ss33 += ss33[index0];
                    sum_s += s[index0];
                }
                ma11 = (2.0*delta2*(sum_ss11/(NBC)))
                        - (2.0*tdelta2*(sum_s/(NBC))*(sum_s11/(NBC)));
                ma12 = (2.0*delta2*(sum_ss12/(NBC)))
                        - (2.0*tdelta2*(sum_s/(NBC))*(sum_s12/(NBC)));
                ma13 = (2.0*delta2*(sum_ss13/(NBC)))
                        - (2.0*tdelta2*(sum_s/(NBC))*(sum_s13/(NBC)));
                ma22 = (2.0*delta2*(sum_ss22/(NBC)))
                        - (2.0*tdelta2*(sum_s/(NBC))*(sum_s22/(NBC)));
                ma23 = (2.0*delta2*(sum_ss23/(NBC)))
                        - (2.0*tdelta2*(sum_s/(NBC))*(sum_s23/(NBC)));
                ma33 = (2.0*delta2*(sum_ss33/(NBC)))
                        - (2.0*tdelta2*(sum_s/(NBC))*(sum_s33/(NBC)));
                la11 = (sum_vel_uu/(NBC))-((sum_vel_u/(NBC))*(sum_vel_u/(NBC)));
                la12 = (sum_vel_uv/(NBC))-((sum_vel_u/(NBC))*(sum_vel_v/(NBC)));
                la13 = (sum_vel_uw/(NBC))-((sum_vel_u/(NBC))*(sum_vel_w/(NBC)));
                la22 = (sum_vel_vv/(NBC))-((sum_vel_v/(NBC))*(sum_vel_v/(NBC)));
                la23 = (sum_vel_vw/(NBC))-((sum_vel_v/(NBC))*(sum_vel_w/(NBC)));
                la33 = (sum_vel_ww/(NBC))-((sum_vel_w/(NBC))*(sum_vel_w/(NBC)));

                ldeno[i][j][k] = ((ma11*ma11) + (ma12*ma12)
                       + (ma13*ma13) + (ma22*ma22)
                       + (ma23*ma23));
                lnume[i][j][k] = ((la11*ma11) + (la12*ma12)
                       + (la13*ma13) + (la22*ma22)
                       + (la23*ma23));

                /*
                numkk = (int)(((cell_center[index000].m_coords[2]-2.538)+(top_h[2]/2.0))/(numr*NB*top_h[2])); 

                deno[numkk] += ((ma11*ma11) + (ma12*ma12)
                       + (ma13*ma13) + (ma22*ma22)
                       + (ma23*ma23));
                nume[numkk] += ((la11*ma11) + (la12*ma12)
                       + (la13*ma13) + (la22*ma22)
                       + (la23*ma23));
                gdeno += ((ma11*ma11) + (ma12*ma12)
                       + (ma13*ma13) + (ma22*ma22)
                       + (ma23*ma23));
                gnume += ((la11*ma11) + (la12*ma12)
                       + (la13*ma13) + (la22*ma22)
                       + (la23*ma23));
                */
        }

        /*
        if (nn > 1)
        {
           for (k = 0; k < numk; k++)
           {
              pp_global_sum(&deno[k],1L);
              pp_global_sum(&nume[k],1L);
           }
           pp_global_sum(&gdeno,1L);
           pp_global_sum(&gnume,1L);
        }

        cs = gnume/gdeno;
        if(cs < 0.0 || deno < 1e-16)
            cs = 0.0;
        */

        double max_cs = 0.0;
        double max_turbulent_kinematic_viscosity = 0.0;

        for (k = 0; k <= ((kmax-kmin+1)/NB)-1; k++)
        for (j = 0; j <= ((jmax-jmin+1)/NB)-1; j++)
        for (i = 0; i <= ((imax-imin+1)/NB)-1; i++)
        {
            kk = (NB*k)+kmin;
            jj = (NB*j)+jmin;
            ii = (NB*i)+imin;
            cs = lnume[i][j][k]/ldeno[i][j][k];
            if(cs < 0.0 || ldeno[i][j][k] < 1e-16)
                cs = 0.0;
           
            if(cs > max_cs)
                max_cs = cs;

            index000 = d_index3d(ii,jj,kk,top_gmax);

            for(kkkk = kk; kkkk < kk+NB; kkkk++)
            for(jjjj = jj; jjjj < jj+NB; jjjj++)
            for(iiii = ii; iiii < ii+NB; iiii++)
            {
                index = d_index3d(iiii,jjjj,kkkk,top_gmax);
                //cell_center[index].m_state.m_mu += cell_center[index].m_state.m_rho*cs*delta2*s[index];
                cell_center[index].m_state.m_mu += cell_center[index].m_state.m_rho*cs*delta2*s[index];

                if(cs*delta2*s[index] > max_turbulent_kinematic_viscosity)
                    max_turbulent_kinematic_viscosity = cs*delta2*s[index];
            }
        }

        pp_global_max(&max_cs,1);
        pp_global_max(&max_turbulent_kinematic_viscosity,1);

        //FT_FreeThese(2,deno,nume);
        FT_FreeThese(2,ldeno,lnume);
        FT_FreeThese(3,u,v,w);
        FT_FreeThese(9,vel_u,vel_v,vel_w,vel_uu,vel_vv,vel_ww,vel_uv,vel_vw,vel_uw);
        FT_FreeThese(13,s,s11,s12,s13,s22,s23,s33,ss11,ss12,ss13,ss22,ss23,ss33);
}       /* end 3D_Cylindrical::computeSubgridModel */

void Incompress_Solver_Smooth_3D_Cylindrical::multiDiagnosisStep(char *out_name)
{
    char filename[500];
    FILE *outfile;


    if(pp_mynode() == 0)
    {
        sprintf(filename, "%s/global-quant.dat",out_name);
        if(front->step == 1)
            outfile = fopen(filename,"w");
        else
            outfile = fopen(filename,"a");
    }

    double area, tridiam_min, tridiam_max;
    double mass_conserv_max, mass_conserv_ave, divu_max, divu_ave;
    double vol_phone, vol_phtwo;

    getInterfaceArea(area);
    getMaxMinDiameterTri(tridiam_min, tridiam_max);
    getMassConservation(mass_conserv_max, mass_conserv_ave, divu_max, divu_ave);
    getPhaseVolume(vol_phone, vol_phtwo);


    if(pp_mynode() == 0)
    {
        if(front->step == 1)
        {
            (void) fprintf(outfile,"#time           intfc_area      phase_vol1      phase_vol2      max_MB_residual ave_MB_residual max_divU        ave_divU        min_tridia      max_tridia\n");
            (void) fprintf(outfile,"#(ms)           (cm^2)          (cm^3)          (cm^3)          (g/(cm^3 ms))   (g/(cm^3 ms))   (1/ms)          (1/ms)          (um)            (um)\n");
        }
        (void) fprintf(outfile,"%.5E     ",front->time);
        (void) fprintf(outfile,"%.5E     ",area);
        (void) fprintf(outfile,"%.5E     ",vol_phone);
        (void) fprintf(outfile,"%.5E     ",vol_phtwo);
        (void) fprintf(outfile,"%.5E     ",mass_conserv_max);
        (void) fprintf(outfile,"%.5E     ",mass_conserv_ave);
        (void) fprintf(outfile,"%.5E     ",divu_max);
        (void) fprintf(outfile,"%.5E     ",divu_ave);
        (void) fprintf(outfile,"%.5E     ",tridiam_min*10000);
        (void) fprintf(outfile,"%.5E     ",tridiam_max*10000);
        (void) fprintf(outfile,"\n");
        fclose(outfile);
    }
}

void Incompress_Solver_Smooth_3D_Cylindrical::recordMeshLength(void)
{
    int i;
    SURFACE **s;

    TRI *tri;
    POINT *p;


    INTERFACE *intfc = front->interf;

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    for (int k = 0; k < 3; ++k)
	    {

	    }
	}
    }
}

void Incompress_Solver_Smooth_3D_Cylindrical::printNormalAndCurvature(void)
{
    INTERFACE *intfc = front->interf;

    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
    POINT *p;
    double curv;
    int i;

    (void) next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
	GetFrontCurvature(p,hse,hs,&curv,front);

	double distance = sqrt( p->_nor[0]*p->_nor[0] + p->_nor[1]*p->_nor[1] + (p->_nor[2] - 1.0)* (p->_nor[2] -1.0) );

	if ( (fabs(p->curvature - 0.172413) > 0.05) || (distance > 0.05) )
	{
	printf("\nThe Point is (%10.8g, %10.8g, %10.8g)\n",p->_coords[0], p->_coords[1], p->_coords[2]);
	printf("Normal is (%10.8g, %10.8g, %10.8g) \n",p->_nor[0], p->_nor[1], p->_nor[2]);
	printf("Curvature is %10.8g\n",p->curvature);
	}
    }

}

void Incompress_Solver_Smooth_3D_Cylindrical::setAdvectionDt(void)
{
    //printf("\nEnter setAdvectionDt in 3D_Cylindrical \n");
    pp_global_max(&max_speed,1);
    if (max_speed != 0.0)
	max_dt = hmin/max_speed;
    else
	max_dt = HUGE;
    min_dt = 0.0000001*sqr(hmin)/mu_min;

    double tension_dt, average_rho;

    average_rho = (m_rho[0] + m_rho[1])/2.0;

    if(m_sigma != 0.0)
	tension_dt = sqrt(average_rho * hmin*hmin*hmin/(2.0*PI*m_sigma));
    else
	tension_dt = HUGE;

    max_dt = std::min(max_dt, tension_dt);

    if (debugging("trace"))
    {
	if (max_dt == HUGE)
	    printf("In setAdvectionDt: \n"
		    "max_dt = HUGE min_dt = %24.18g\n", min_dt);
	else
	    printf("In setAdvectionDt: \n"
		    "max_dt = %24.18g min_dt = %24.18g\n", max_dt, min_dt);
    }
}   /* end setAdvectionDt including surfacetesion effect */


void Incompress_Solver_Smooth_3D_Cylindrical::surfaceTension_Peskin(void)
{
    printf("\nEnter 3D_Cylindrical::surfaceTension_Peskin \n");

    INTERFACE* intfc = front->interf;

    //Go through all the triangles
    SURFACE **s;
    TRI *tri;

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    compSurfaceTension_Tri(tri);
	}
    }
}

void Incompress_Solver_Smooth_3D_Cylindrical::compSurfaceTension_Tri(TRI* triangle)
{
    POINT *point_tri;
    int i,j,k;
    //int range = m_smoothing_radius;
    int range = 2;
    int icoords[3];
    double tri_center[3], tri_nor[3];
    double grid_center[3];
    double force_on_tri[3], *force_on_cell;

    double tri_area = 0.0;
    double tri_curvature = 0.0;
    double point_curvature = 0.0;
    double mag_nor = 0.0;

    int cell_imin, cell_jmin, cell_kmin;
    int cell_imax, cell_jmax, cell_kmax;

    double prec[3][3];
    double edge[3][3];

    int index;

    for (int m = 0; m < 3; m++)
    {
	tri_center[m] = tri_nor[m] = grid_center[m]
	    = force_on_tri[m] = 0.0;
    }

    for (int p = 0; p < 3; p++)
    {
	point_tri = Point_of_tri(triangle)[p];

	prec[p][0] = Coords(point_tri)[2]*cos(Coords(point_tri)[0]);
	prec[p][1] = Coords(point_tri)[2]*sin(Coords(point_tri)[0]);
	prec[p][2] = Coords(point_tri)[1];
	

	for (int l = 0; l < 3; l++)
	{
	    tri_center[l] += Coords(point_tri)[l];
	}

	point_curvature = point_tri->curvature;
	tri_curvature += point_curvature;
    }
    tri_curvature /= 3.0;

    for (int l = 0; l < 3; l++)
    {
	tri_center[l] /= 3.0;
	tri_nor[l] = Tri_normal(triangle)[l];
    }

    mag_nor = mag_vector(tri_nor,3);

    //tri_area might need to be revised

    for (int m = 0; m < 3; m++)
    {
	edge[0][m] = prec[1][m] - prec[0][m];
	edge[1][m] = prec[2][m] - prec[1][m];
	edge[2][m] = prec[0][m] - prec[2][m];
    }
    double edge_len1 = mag_vector(edge[0], 3);
    double edge_len2 = mag_vector(edge[1], 3);
    double edge_len3 = mag_vector(edge[2], 3);
    double hallen = (edge_len1 + edge_len2 + edge_len3) / 2.0;

    tri_area = sqrt( hallen*(hallen - edge_len1)*(hallen - edge_len2)*(hallen - edge_len3) );

    /////////////////////////////////
    
    for (int l = 0; l < 3; l++)
    {
	tri_nor[l] /= mag_nor;
	force_on_tri[l] = 2.0 * m_sigma * tri_curvature * tri_nor[l] * tri_area;
        
        if(tri_area < 1e-16)
            force_on_tri[l] = 0.0;
    }



    for (int l = 0; l < 3; l++)
	icoords[l] = cell_index(tri_center[l],l,top_grid);

    cell_imin = std::max(imin, icoords[0] - range);
    cell_jmin = std::max(jmin, icoords[1] - range);
    cell_kmin = std::max(kmin, icoords[2] - range);

    cell_imax = std::min(imax, icoords[0] + range);
    cell_jmax = std::min(jmax, icoords[1] + range);
    cell_kmax = std::min(kmax, icoords[2] + range);

    for (k = cell_kmin; k <= cell_kmax; k++)
    for (j = cell_jmin; j <= cell_jmax; j++)
    for (i = cell_imin; i <= cell_imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);
	for (int m = 0; m < 3; m++)
	    grid_center[m] = cell_center[index].m_coords[m];

	double D = getDiscreteDelta(tri_center, grid_center);

	double rho = cell_center[index].m_state.m_rho;
	//force_on_cell = cell_center[index].m_state.f_surf;

	for (int m = 0; m < 3; m++)
	    cell_center[index].m_state.f_surf[m] += (-force_on_tri[m]*D/rho);
    }
}

double Incompress_Solver_Smooth_3D_Cylindrical::getDiscreteDelta(
	const double *tri_center,
	const double *grid_center)
{
    double x[3];

    for (int k = 0; k < 3; k++)
	x[k] = tri_center[k] - grid_center[k];

    return 1/(grid_center[2]*top_h[0]*top_h[1]*top_h[2]) *
	     getDiscreteDelta_Phi(x[0]/top_h[0]) *
	     getDiscreteDelta_Phi(x[1]/top_h[1]) *
	     getDiscreteDelta_Phi(x[2]/top_h[2]);

}

double Incompress_Solver_Smooth_3D_Cylindrical::getDiscreteDelta_Phi(
	double r)
{
    double R = fabs(r);
    double ret;

    if(R<1)
	ret = 1.0/8 * (3 - 2*R + sqrt(1+4*R-4*r*r));
    else if (R>=1 && R<2)
	ret = 1.0/8 * (5 - 2*R - sqrt(-7+12*R-4*r*r));
    else
	ret = 0;

    return ret;

}

void Incompress_Solver_Smooth_3D_Cylindrical::setRhoMuOld(void)
{
    int i,j,k,index;

    for (k = 0; k <= top_gmax[2]; k++)
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
	index = d_index3d(i,j,k,top_gmax);
	cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
	cell_center[index].m_state.m_mu_old = cell_center[index].m_state.m_mu;
    }
}


void Incompress_Solver_Smooth_3D_Cylindrical::getMassConservation(double &mass_conserv_max, double &mass_conserv_ave, double &divu_max, double &divu_ave)
{
    int i,j,k,index;
    int icoords[3];

    int index_nb[6];

    double rho_theta_nb[2], rho_r_nb[2], rho_z_nb[2];

    double rur_nb[2], utheta_nb[2], uz_nb[2];

    double num_cells = 0;
    double divu;
    double divrhou;
    double divu_sum = 0.0;
    divu_max = 0.0;
    double rhot_sum = 0.0;
    double rhot_max = 0.0;

    double mass_conserv;
    double mass_conserv_sum = 0.0;
    mass_conserv_max = 0.0;

    const int nn = pp_numnodes();
    int myid = pp_mynode();
    int *ppgmax = front->pp_grid->gmax;
    int ppx = myid % ppgmax[0];
    int ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
    int ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);
	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	double r, r_nb[2];
	r = cell_center[index].m_coords[2];
	r_nb[0] = (cell_center[index_nb[4]].m_coords[2] + r)/2.0;
	r_nb[1] = (cell_center[index_nb[5]].m_coords[2] + r)/2.0;

	utheta_nb[0] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0])/2.0;
	utheta_nb[1] = (cell_center[index].m_state.m_U[0] + cell_center[index_nb[1]].m_state.m_U[0])/2.0;

	rho_theta_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[0]].m_state.m_rho)/2.0;
	rho_theta_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[1]].m_state.m_rho)/2.0;

	uz_nb[0] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1])/2.0;
	uz_nb[1] = (cell_center[index].m_state.m_U[1] + cell_center[index_nb[3]].m_state.m_U[1])/2.0;

	rho_z_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[2]].m_state.m_rho)/2.0;
	rho_z_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[3]].m_state.m_rho)/2.0;

	if ( (ppz == 0 && k == kmin) )
	{
	    rur_nb[0] = 0.0;
	    rho_r_nb[0] = cell_center[index].m_state.m_rho;
	}
	else
	{
	    rur_nb[0] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2])/2.0) * r_nb[0];
	    rho_r_nb[0] = (cell_center[index].m_state.m_rho + cell_center[index_nb[4]].m_state.m_rho)/2.0;
	}

	if ( (ppz == (ppgmax[2]-1) && k == kmax) )
	{
	    rur_nb[1] = 0.0;
	    rho_r_nb[1] = cell_center[index].m_state.m_rho;
	}
	else
	{
	    rur_nb[1] = ((cell_center[index].m_state.m_U[2] + cell_center[index_nb[5]].m_state.m_U[2])/2.0) * r_nb[1];
	    rho_r_nb[1] = (cell_center[index].m_state.m_rho + cell_center[index_nb[5]].m_state.m_rho)/2.0;
	}

	divu = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) +
	       (uz_nb[1] - uz_nb[0])/top_h[1] +
	       (rur_nb[1] - rur_nb[0])/(r*top_h[2]);

	divrhou = (utheta_nb[1]*rho_theta_nb[1] - utheta_nb[0]*rho_theta_nb[0]) / (r*top_h[0]) +
	    	  (uz_nb[1]*rho_z_nb[1] - uz_nb[0]*rho_z_nb[0]) / top_h[1] +
		  (rur_nb[1]*rho_r_nb[1] - rur_nb[0]*rho_r_nb[0]) / (r*top_h[2]);

	double rhot = (cell_center[index].m_state.m_rho - cell_center[index].m_state.m_rho_old)/m_dt;

	mass_conserv = rhot + divrhou;

	if (fabs(divu) > divu_max)
	    divu_max = fabs(divu);
	if (fabs(mass_conserv) > mass_conserv_max)
	    mass_conserv_max = fabs(mass_conserv);
	if (fabs(rhot) > rhot_max)
	    rhot_max = fabs(rhot);

	num_cells += 1.0;

	divu_sum += fabs(divu);
	mass_conserv_sum += fabs(mass_conserv);
	rhot_sum += fabs(rhot);


    }

    pp_global_max(&divu_max , 1);
    pp_global_max(&mass_conserv_max, 1);
    pp_global_max(&rhot_max, 1);

    pp_global_sum(&num_cells , 1);
    pp_global_sum(&divu_sum , 1);
    pp_global_sum(&mass_conserv_sum, 1);
    pp_global_sum(&rhot_sum, 1);


    divu_ave = divu_sum/num_cells;
    mass_conserv_ave = mass_conserv_sum/num_cells;
}


void Incompress_Solver_Smooth_3D_Cylindrical::printTridia(char *out_name)
{
    INTERFACE* intfc = front->interf;
    SURFACE **s;
    TRI *tri;
    int p_in_domain = 0;
    POINT *p;
    double diameter_tri;
    double prec[3][3];
    double edge[3][3];
    double tot_num_tris = 0;

    char filename[500];
    FILE *outfile;

    
    sprintf(filename,"%s/tridia",out_name);
    if (pp_min_status(create_directory(filename,NO)) == NO)
    {
        screen("Cannot create directory %s\n",filename);
        clean_up(ERROR);
    }

    {
        sprintf(filename, "%s/tridia-ts%s",filename,right_flush(front->step,7));
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-p%s",filename,right_flush(pp_mynode(),4));
#endif
        sprintf(filename,"%s.dat",filename);

        if(front->step == 1)
            outfile = fopen(filename,"w");
        else
            outfile = fopen(filename,"a");
    }

    for (s = intfc->surfaces; s && *s; ++s)
    {
        for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
        {
            for (int i = 0; i < 3; i++)
            {
                p = Point_of_tri(tri)[i];
                if ( (Coords(p)[0] >= comp_grid->L[0]) &&
                     (Coords(p)[0] <= comp_grid->U[0]) &&
                     (Coords(p)[1] >= comp_grid->L[1]) &&
                     (Coords(p)[1] <= comp_grid->U[1]) &&
                     (Coords(p)[2] >= comp_grid->L[2]) &&
                     (Coords(p)[2] <= comp_grid->U[2]) )
                     p_in_domain++;
            }

            if (p_in_domain > 1)
            {
                tot_num_tris++;
                for (int i = 0; i < 3; i++)
                {
                    p = Point_of_tri(tri)[i];

                    prec[i][0] = Coords(p)[2]*cos(Coords(p)[0]);
                    prec[i][1] = Coords(p)[2]*sin(Coords(p)[0]);
                    prec[i][2] = Coords(p)[1];
                }

                for (int m = 0; m < 3; m++)
                {
                    edge[0][m] = prec[1][m] - prec[0][m];
                    edge[1][m] = prec[2][m] - prec[1][m];
                    edge[2][m] = prec[0][m] - prec[2][m];
                }
                double edge_len1 = mag_vector(edge[0], 3);
                double edge_len2 = mag_vector(edge[1], 3);
                double edge_len3 = mag_vector(edge[2], 3);
                double hallen = (edge_len1 + edge_len2 + edge_len3) / 2.0;

                diameter_tri = 2.0*sqrt( (hallen - edge_len1)*(hallen - edge_len2)*(hallen - edge_len3)/hallen );
                fprintf(outfile, "%.5E\n", diameter_tri);
            }
        }
    }
    fclose(outfile);
    pp_global_sum(&tot_num_tris,1);
    fprintf(stdout, "#total number of triangles %.5E\n", tot_num_tris);
}


void Incompress_Solver_Smooth_3D_Cylindrical::outputParallelVisitFile(char *out_name, bool binary)
{
    int nn = pp_numnodes();
    char filename_pvtu_intfc[200];
    char filename_pvtu_ph1[200];
    char filename_pvtu_ph2[200];
    char dir_name[256];

    FILE *outfile_intfc;
    FILE *outfile_ph1;
    FILE *outfile_ph2;

    sprintf(dir_name, "%s/parallel-visual", out_name);
    if (!create_directory(dir_name,YES))
    {
	screen("Cannot create directory %s\n",dir_name);
	clean_up(ERROR);
    }
    sprintf(filename_pvtu_intfc, "%s/parallel-intfc-t%s.pvtp", dir_name, right_flush(front->step,7));
    sprintf(filename_pvtu_ph1, "%s/parallel-ph1-t%s.pvtu", dir_name, right_flush(front->step,7));
    sprintf(filename_pvtu_ph2, "%s/parallel-ph2-t%s.pvtu", dir_name, right_flush(front->step,7));

    outfile_intfc = fopen(filename_pvtu_intfc, "w");
    outfile_ph1 = fopen(filename_pvtu_ph1, "w");
    outfile_ph2 = fopen(filename_pvtu_ph2, "w");

    fprintf(outfile_ph1,"<\?xml version=\"1.0\"\?>\n");

    if (binary == YES)
    {
	if(hardware_is_little_endian())
	    fprintf(outfile_ph1,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	else
	    fprintf(outfile_ph1,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
    }
    else
	fprintf(outfile_ph1,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");

    fprintf(outfile_ph1,"  <PUnstructuredGrid GhostLevel=\"0\">\n");

    fprintf(outfile_ph1,"    <PPoints>\n");
    fprintf(outfile_ph1,"      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n");
    fprintf(outfile_ph1,"    </PPoints>\n");

    fprintf(outfile_ph1,"    <PPointData>\n");
    fprintf(outfile_ph1,"    </PPointData>\n");

    fprintf(outfile_ph1,"    <PCellData>\n");
    if (iFparams->movie_option->plot_phase_1)
    {
	if (iFparams->movie_option->plot_velo)
	    fprintf(outfile_ph1, "      <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n");
	if (iFparams->movie_option->plot_pres)
	    fprintf(outfile_ph1, "      <PDataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\"/>\n");
	if (iFparams->movie_option->plot_comp)
	    fprintf(outfile_ph1, "      <PDataArray type=\"Int32\" Name=\"phase\" NumberOfComponents=\"1\"/>\n");
	if (iFparams->movie_option->plot_dens)
	    fprintf(outfile_ph1, "      <PDataArray type=\"Float64\" Name=\"density\" NumberOfComponents=\"1\"/>\n");
    }
    fprintf(outfile_ph1,"    </PCellData>\n");
    for (int i = 0; i < nn; i++)
    {
	char piece_node[10];
	sprintf(piece_node, "%s", right_flush(i,4));
	char piece_source[256];
	sprintf(piece_source, "    <Piece Source=\"../P-%s/phase_1/ph1-t%s-p%s.vtu\"/>\n",piece_node, right_flush(front->step, 7), piece_node);
	fprintf(outfile_ph1, piece_source);
    }

    fprintf(outfile_ph1,"  </PUnstructuredGrid>\n");
    fprintf(outfile_ph1,"</VTKFile>\n");


    fprintf(outfile_ph2,"<\?xml version=\"1.0\"\?>\n");

    if (binary == YES)
    {
	if(hardware_is_little_endian())
	    fprintf(outfile_ph2,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	else
	    fprintf(outfile_ph2,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
    }
    else
	fprintf(outfile_ph2,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");

    fprintf(outfile_ph2,"  <PUnstructuredGrid GhostLevel=\"0\">\n");

    fprintf(outfile_ph2,"    <PPointData>\n");
    fprintf(outfile_ph2,"    </PPointData>\n");

    fprintf(outfile_ph2,"    <PCellData>\n");
    if (iFparams->movie_option->plot_phase_2)
    {
	if (iFparams->movie_option->plot_velo)
	    fprintf(outfile_ph2,"      <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n");
	if (iFparams->movie_option->plot_pres)
	    fprintf(outfile_ph2,"      <PDataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\"/>\n");
	if (iFparams->movie_option->plot_comp)
	    fprintf(outfile_ph2,"      <PDataArray type=\"Int32\" Name=\"phase\" NumberOfComponents=\"1\"/>\n");
	if (iFparams->movie_option->plot_dens)
	    fprintf(outfile_ph2,"      <PDataArray type=\"Float64\" Name=\"density\" NumberOfComponents=\"1\"/>\n");
    }
    fprintf(outfile_ph2,"    </PCellData>\n");

    fprintf(outfile_ph2,"    <PPoints>\n");
    fprintf(outfile_ph2,"      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n");
    fprintf(outfile_ph2,"    </PPoints>\n");

    for (int i = 0; i < nn; i++)
    {
	char piece_node[10];
	sprintf(piece_node, "%s", right_flush(i,4));
	char piece_source[256];
	sprintf(piece_source, "    <Piece Source=\"../P-%s/phase_2/ph2-t%s-p%s.vtu\"/>\n",piece_node, right_flush(front->step, 7), piece_node);
	fprintf(outfile_ph2, piece_source);    
    }

    fprintf(outfile_ph2,"  </PUnstructuredGrid>\n");
    fprintf(outfile_ph2,"</VTKFile>\n");



    fprintf(outfile_intfc,"<\?xml version=\"1.0\"\?>\n");

    if (binary == YES)
    {
	if(hardware_is_little_endian())
	    fprintf(outfile_intfc,"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	else
	    fprintf(outfile_intfc,"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"BigEndian\">\n");
    }
    else
	fprintf(outfile_intfc,"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");

    fprintf(outfile_intfc,"  <PPolyData GhostLevel=\"0\">\n");

    fprintf(outfile_intfc,"    <PPoints>\n");
    fprintf(outfile_intfc,"      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n");
    fprintf(outfile_intfc,"    </PPoints>\n");

    fprintf(outfile_intfc,"    <PPointData>\n");
    fprintf(outfile_intfc,"      <PDataArray type=\"Float64\" Name=\"normal\" NumberOfComponents=\"3\"/>\n");
    fprintf(outfile_intfc,"      <PDataArray type=\"Float64\" Name=\"curvature\" NumberOfComponents=\"1\"/>\n");
    fprintf(outfile_intfc,"    </PPointData>\n");

    fprintf(outfile_intfc,"    <PCellData>\n");
    fprintf(outfile_intfc,"    </PCellData>\n");
    for (int i = 0; i < nn; i++)
    {
	char piece_node[10];
	sprintf(piece_node, "%s", right_flush(i,4));
	char piece_source[256];
	sprintf(piece_source, "    <Piece Source=\"../P-%s/interface/intfc.t%s-p%s.vtp\"/>\n",piece_node, right_flush(front->step, 7), piece_node);
	fprintf(outfile_intfc, piece_source);      
    }

    fprintf(outfile_intfc,"  </PPolyData>\n");
    fprintf(outfile_intfc,"</VTKFile>\n");

    fclose(outfile_intfc);
    fclose(outfile_ph1);
    fclose(outfile_ph2);

}

void Incompress_Solver_Smooth_3D_Cylindrical::multiDiagnosisInterval(char *out_name)
{
    if(iFparams->movie_option->output_tridia == YES)
	printTridia(out_name);

    if(iFparams->movie_option->output_error == YES)
    {
	static bool first_try = true;
	char filename1[500];//error output file
	char filename2[500];//for other interval outputs
	FILE *outfile1;
	FILE *outfile2;

	//if(pp_mynode() == 0)
	{
	    sprintf(filename1, "%s/error.dat",out_name);
	    if(first_try == true)
		outfile1 = fopen(filename1,"w");
	    else
		outfile1 = fopen(filename1,"a");

	    computeError(outfile1);
	    computeError_part(outfile1);
	    //printNormalAndCurvature();
    
	    first_try = false;
	    fclose(outfile1);
	}
    }
}

void Incompress_Solver_Smooth_3D_Cylindrical::adjustReflectionBuffer(void)
{
    const int  nn = pp_numnodes();
    int        myid = pp_mynode();
    int   *ppgmax = front->pp_grid->gmax;
    int   ppx = myid % ppgmax[0];
    int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
    int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

    int i,j,k,index;

    double **vel = field->vel;

    if (ppz == 0)
    {
	for (i = 0; i <= top_gmax[0]; i++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (k = 0; k < 4; k++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    vel[2][index] = -vel[2][index];
	}

    }

    if (ppz == (ppgmax[2] - 1) )
    {
	for (i = 0; i <= top_gmax[0]; i++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (k = 0; k < 4; k++)
	{
	    index = d_index3d(i,j,kmax+k+1,top_gmax);
	    vel[2][index] = -vel[2][index];
	}

    }
}



void Incompress_Solver_Smooth_3D_Cylindrical::getInterfaceArea(double &area)
{
    INTERFACE* intfc = front->interf;

    SURFACE **s;
    TRI *tri;

    area = 0.0;

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	   area += getTriArea(tri); 
	}
    }

    pp_global_sum(&area, 1);
}

void Incompress_Solver_Smooth_3D_Cylindrical::getMaxMinDiameterTri(double &min, double &max)
{
    INTERFACE* intfc = front->interf;

    SURFACE **s;
    TRI *tri;
    min = HUGE;
    max = 0.0;
    int p_in_domain = 0;
    POINT *p;
    double min_diameter_tri, max_diameter_tri;

    double prec[3][3];
    double edge[3][3];

    for (s = intfc->surfaces; s && *s; ++s)
    {
        for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
        {
            for (int i = 0; i < 3; i++)
            {
                p = Point_of_tri(tri)[i];
                if ( (Coords(p)[0] >= comp_grid->L[0]) &&
                     (Coords(p)[0] <= comp_grid->U[0]) &&
                     (Coords(p)[1] >= comp_grid->L[1]) &&
                     (Coords(p)[1] <= comp_grid->U[1]) &&
                     (Coords(p)[2] >= comp_grid->L[2]) &&
                     (Coords(p)[2] <= comp_grid->U[2]) )
                     p_in_domain++;
            }

            if (p_in_domain > 1)
            {
                for (int i = 0; i < 3; i++)
                {
                    p = Point_of_tri(tri)[i];

                    prec[i][0] = Coords(p)[2]*cos(Coords(p)[0]);
                    prec[i][1] = Coords(p)[2]*sin(Coords(p)[0]);
                    prec[i][2] = Coords(p)[1];
                }

                for (int m = 0; m < 3; m++)
                {
                    edge[0][m] = prec[1][m] - prec[0][m];
                    edge[1][m] = prec[2][m] - prec[1][m];
                    edge[2][m] = prec[0][m] - prec[2][m];
                }
                double edge_len1 = mag_vector(edge[0], 3);
                double edge_len2 = mag_vector(edge[1], 3);
                double edge_len3 = mag_vector(edge[2], 3);
                double hallen = (edge_len1 + edge_len2 + edge_len3) / 2.0;

                min_diameter_tri = 2.0*sqrt( (hallen - edge_len1)*(hallen - edge_len2)*(hallen - edge_len3)/hallen );
                max_diameter_tri = 2.0*sqrt( (hallen - edge_len1)*(hallen - edge_len2)*(hallen - edge_len3)/hallen );
            }
            else
            {
                min_diameter_tri = HUGE;
                max_diameter_tri = 0.0;
            }

           if(min > min_diameter_tri)
               min = max_diameter_tri;
           if(max < max_diameter_tri)
               max = max_diameter_tri;
        }
    }

    pp_global_min(&min, 1);
    pp_global_max(&max, 1);
}

double Incompress_Solver_Smooth_3D_Cylindrical::getTriArea(TRI* triangle)
{
    int p_in_domain = 0;
    POINT *p;
    double tri_area;

    double prec[3][3];
    double edge[3][3];

    for (int i = 0; i < 3; i++)
    {
	p = Point_of_tri(triangle)[i];

	if ( (Coords(p)[0] >= comp_grid->L[0]) &&
	     (Coords(p)[0] <= comp_grid->U[0]) &&
	     (Coords(p)[1] >= comp_grid->L[1]) &&
	     (Coords(p)[1] <= comp_grid->U[1]) &&
	     (Coords(p)[2] >= comp_grid->L[2]) &&
	     (Coords(p)[2] <= comp_grid->U[2]) )
	    p_in_domain++;
    }

    if (p_in_domain > 1)
    {
	for (int i = 0; i < 3; i++)
	{
	    p = Point_of_tri(triangle)[i];
	    
	    prec[i][0] = Coords(p)[2]*cos(Coords(p)[0]);
	    prec[i][1] = Coords(p)[2]*sin(Coords(p)[0]);
	    prec[i][2] = Coords(p)[1];
	}

	for (int m = 0; m < 3; m++)
	{
	    edge[0][m] = prec[1][m] - prec[0][m];
	    edge[1][m] = prec[2][m] - prec[1][m];
	    edge[2][m] = prec[0][m] - prec[2][m];
	}

	double edge_len1 = mag_vector(edge[0], 3);
	double edge_len2 = mag_vector(edge[1], 3);
	double edge_len3 = mag_vector(edge[2], 3);
	double hallen = (edge_len1 + edge_len2 + edge_len3) / 2.0;
	
	tri_area = sqrt( hallen*(hallen - edge_len1)*(hallen - edge_len2)*(hallen - edge_len3) );

	return tri_area;
    }
    else 
	return 0.0;
    
}

void Incompress_Solver_Smooth_3D_Cylindrical::getIcoords(int index, int *icoords, int *gmax)
{
    icoords[0] = index % (gmax[0]+1);
    icoords[1] = ((index-icoords[0])/(gmax[0]+1)) % (gmax[1]+1);
    icoords[2] = (index-icoords[0]-((gmax[0]+1)*icoords[1]))/((gmax[0]+1)*(gmax[1]+1));
}


void Incompress_Solver_Smooth_3D_Cylindrical::getPhaseVolume(double &vol_phone, double &vol_phtwo)
{
    int i,j,k,index;
    vol_phone = 0.0;
    vol_phtwo = 0.0;
    double r;
    double cellVolume;

    for (i = imin; i <= imax; i++)
    for (j = jmin; j <= jmax; j++)
    for (k = kmin; k <= kmax; k++)
    {
	index = d_index3d(i,j,k,top_gmax);
	r = cell_center[index].m_coords[2];
	cellVolume = r*top_h[0]*top_h[1]*top_h[2];
	if (cell_center[index].comp == LIQUID_COMP1)
	    vol_phone += cellVolume;
	else if(cell_center[index].comp == LIQUID_COMP2)
	    vol_phtwo += cellVolume;
	else
	{
	    screen("Unphysical comp in the computational domain!");
	    clean_up(ERROR);
	}
	
    }

    pp_global_sum(&vol_phone, 1);
    pp_global_sum(&vol_phtwo, 1);
}


