/*******************************************************************
 * 		         iFbasic.cpp
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

//-------------------------------------------------------------------
//		L_STATE
//-------------------------------------------------------------------
L_STATE::L_STATE(): m_P(0)
{
	int i;
	for (i = 0; i < MAXD; ++i)
	{
	    m_U[i] = 0;
	    grad_q[i] = 0;
	    f_surf[i] = 0;
	}
	m_P = 0;
	m_phi = 0;
	m_q = 0;
	m_mu = 0;
	m_rho = 0;
	div_U = 0;
}
void L_STATE::setZero()
{
	int i;
	for (i = 0; i < MAXD; ++i)
	{
	    m_U[i] = 0;
	    grad_q[i] = 0;
	    f_surf[i] = 0;
	}
	m_P = 0;
	m_q = 0;
	m_phi = 0;
	m_mu = 0;
	m_rho = 0;
	div_U = 0;
}


//----------------------------------------------------------------
//		L_RECTANGLE
//----------------------------------------------------------------
void L_POINT::setCoords(double *coords, int dim)
{
    int i;
    for (i = 0; i < dim; i++)
	p_coords[i] = coords[i];
}

L_RECTANGLE::L_RECTANGLE(): m_index(-1), comp(-1)
{
}

void L_RECTANGLE::setCoords(
	double *coords,
	int dim)
{
	int i;
	for (i = 0; i < dim; ++i)
	    m_coords[i] = coords[i];
}
//--------------------------------------------------------------------------
//               Incompress_Solver_Basis
//               Pure virtual class	
//--------------------------------------------------------------------------



//-------------------------------------------------------------------------------
//               Incompress_Solver_Smooth_Basis
//------------------------------------------------------------------------------
Incompress_Solver_Smooth_Basis::Incompress_Solver_Smooth_Basis(Front &front):front(&front)
{
}

boolean Incompress_Solver_Smooth_Basis::FT_StateStructAtGridCrossing_tmp(Front *front, int *icoords, GRID_DIRECTION dir, COMPONENT comp, Locstate *state, HYPER_SURF **hs, double *crx_coords, double t)
{
    return FT_StateStructAtGridCrossing(front,icoords,dir,comp,state,hs,crx_coords);
}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup cell_center
//---------------------------------------------------------------
void Incompress_Solver_Smooth_Basis::setCornerMesh(void)
{
    int index;
    double coords[MAXD];
    L_POINT m_point;
    int i,j,k,l,lindex;
    int num_corners;
    int icoords[MAXD];
    int loc_index[8];

    int num_nodes = pp_numnodes();
    int myid = pp_mynode();


    int ixlower, iylower, izlower;
    
    int Nxpps = pp_grid->gmax[0];
    int Nypps = pp_grid->gmax[1];
    int Nzpps = pp_grid->gmax[2];

    int *n_xdist;
    int *n_ydist;
    int *n_zdist;

    FT_VectorMemoryAlloc((POINTER*)&n_xdist, Nxpps, sizeof(int));
    FT_VectorMemoryAlloc((POINTER*)&n_ydist, Nypps, sizeof(int));
    FT_VectorMemoryAlloc((POINTER*)&n_zdist, Nzpps, sizeof(int));


    int my_xid = myid % Nxpps;
    int my_yid = ((myid - my_xid)/Nxpps) % Nypps;
    int my_zid = (myid-my_xid-(Nxpps*my_yid))/(Nxpps*Nypps);

    //printf("\nMy x id = %d, y id = %d, z id = %d\n",my_xid, my_yid, my_zid);


    for (i = 0; i < Nxpps; i++) n_xdist[i] = 0;
    for (i = 0; i < Nypps; i++) n_ydist[i] = 0;
    for (i = 0; i < Nzpps; i++) n_zdist[i] = 0;

    n_xdist[my_xid] = comp_grid->gmax[0];
    n_ydist[my_yid] = comp_grid->gmax[1];
    n_zdist[my_zid] = comp_grid->gmax[2];

    pp_global_imax(n_xdist, Nxpps);
    pp_global_imax(n_ydist, Nypps);
    pp_global_imax(n_zdist, Nzpps);

    ixlower = 0; iylower = 0; izlower = 0;

    for (i = 1; i <= my_xid; i++) ixlower += n_xdist[i-1];
    for (i = 1; i <= my_yid; i++) iylower += n_ydist[i-1];
    for (i = 1; i <= my_zid; i++) izlower += n_zdist[i-1];

    //printf("\nixlower = %d, iylower = %d, izlower= %d\n",ixlower, iylower, izlower);
    //printf("\ndim = %d\n", dim);

    num_corners = 1;
    for (i = 0; i < dim; i++)
	num_corners *= (comp_grid->gmax[i] + 1);
    cell_corner.insert(cell_corner.end(),num_corners,m_point);

    //printf("\nthe size of cell_corner:%d\n", (int) cell_corner.size());

    switch (dim)
    {
    case 2:
	for (j = 0; j <= comp_grid->gmax[1]; j++)
	for (i = 0; i <= comp_grid->gmax[0]; i++)
	{
	    coords[0] = comp_grid->L[0] + comp_grid->h[0]*i;
	    coords[1] = comp_grid->L[1] + comp_grid->h[1]*j;
	    lindex = d_index2d(i,j,comp_grid->gmax);
	    cell_corner[lindex].setCoords(coords,dim);
	    cell_corner[lindex].pl_index
		= d_index2d(i,j,comp_grid->gmax);
	    cell_corner[lindex].pg_index 
		=d_index2d(i + ixlower, j + iylower, pp_grid->Global_grid.gmax);
	}
	break;
    case 3:
	for (k = 0; k <= comp_grid->gmax[2]; k++)
	for (j = 0; j <= comp_grid->gmax[1]; j++)
	for (i = 0; i <= comp_grid->gmax[0]; i++)
	{
	    coords[0] = comp_grid->L[0] + comp_grid->h[0]*i;
	    coords[1] = comp_grid->L[1] + comp_grid->h[1]*j;
	    coords[2] = comp_grid->L[2] + comp_grid->h[2]*k;
	    lindex = d_index3d(i,j,k,comp_grid->gmax);
	    cell_corner[lindex].setCoords(coords,dim);
	    cell_corner[lindex].pl_index = lindex;

	    cell_corner[lindex].pg_index 
		= d_index3d(i+ixlower, j+iylower, k+izlower, pp_grid->Global_grid.gmax);

	}

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    for (l = 0; l < 8; l++)
		cell_center[index].corner[l] = 0;
	}

	// Buffer cells does not have pointers to corners

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    loc_index[0] = d_index3d(i-imin,  j-jmin,  k-kmin,  comp_grid->gmax);
	    loc_index[1] = d_index3d(i-imin+1,j-jmin,  k-kmin,  comp_grid->gmax);
	    loc_index[2] = d_index3d(i-imin+1,j-jmin+1,k-kmin,  comp_grid->gmax);
	    loc_index[3] = d_index3d(i-imin,  j-jmin+1,k-kmin,  comp_grid->gmax);
	    loc_index[4] = d_index3d(i-imin,  j-jmin,  k-kmin+1,comp_grid->gmax);
	    loc_index[5] = d_index3d(i-imin+1,j-jmin,  k-kmin+1,comp_grid->gmax);
	    loc_index[6] = d_index3d(i-imin+1,j-jmin+1,k-kmin+1,comp_grid->gmax);
	    loc_index[7] = d_index3d(i-imin,  j-jmin+1,k-kmin+1,comp_grid->gmax);

	    for (l = 0; l < 8; l++)
		cell_center[index].corner[l] = &cell_corner[loc_index[l]];
	}
	break;
    }

    FT_FreeThese(3,n_xdist,n_ydist,n_zdist);

}

void Incompress_Solver_Smooth_Basis::setCornerComp(void)
{

    int i,j,k,lindex;
    int comp;

    for (k = 0; k <= comp_grid->gmax[2]; k++)
    for (j = 0; j <= comp_grid->gmax[1]; j++)
    for (i = 0; i <= comp_grid->gmax[0]; i++)
    {
	lindex = d_index3d(i,j,k,comp_grid->gmax);
	comp = getComponent(cell_corner[lindex].p_coords);
	cell_corner[lindex].comp = comp;
    }

    /*
    
    printf("\n comp_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    comp_grid->L[0],comp_grid->U[0],comp_grid->L[1],comp_grid->U[1],
	    comp_grid->L[2],comp_grid->U[2]);
    printf("\n space is %10.5g, %10.5g, %10.5g \n", comp_grid->h[0], comp_grid->h[1],
	    comp_grid->h[2]);
    printf("\n gmax is %d, %d, %d\n", comp_grid->gmax[0],comp_grid->gmax[1],
	    comp_grid->gmax[2]);
    printf("\n Global mesh for comp_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    comp_grid->GL[0],comp_grid->GU[0],comp_grid->GL[1],comp_grid->GU[1],
	    comp_grid->GL[2],comp_grid->GU[2]);
    printf("\n Virtual mesh for comp_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    comp_grid->VL[0],comp_grid->VU[0],comp_grid->VL[1],comp_grid->VU[1],
	    comp_grid->VL[2],comp_grid->VU[2]);
    printf("\n lbuf for comp_grid \n");
    printf("\n %d, %d, %d \n", comp_grid->lbuf[0],comp_grid->lbuf[1],
	    comp_grid->lbuf[2]);
    printf("\n ubuf for comp_grid \n");
    printf("\n %d, %d, %d \n", comp_grid->ubuf[0],comp_grid->ubuf[1],
	    comp_grid->ubuf[2]);

    printf("\n top_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    top_grid->L[0],top_grid->U[0],top_grid->L[1],top_grid->U[1],
	    top_grid->L[2],top_grid->U[2]);
    printf("\n space is %10.5g, %10.5g, %10.5g \n", top_grid->h[0], top_grid->h[1],
	    top_grid->h[2]);
    printf("\n gmax is %d, %d, %d\n", top_grid->gmax[0],top_grid->gmax[1],
	    top_grid->gmax[2]);
    printf("\n Global mesh for top_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    top_grid->GL[0],top_grid->GU[0],top_grid->GL[1],top_grid->GU[1],
	    top_grid->GL[2],top_grid->GU[2]);
    printf("\n Virtual mesh for top_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    top_grid->VL[0],top_grid->VU[0],top_grid->VL[1],top_grid->VU[1],
	    top_grid->VL[2],top_grid->VU[2]);
    printf("\n lbuf for top_grid \n");
    printf("\n %d, %d, %d \n", top_grid->lbuf[0],top_grid->lbuf[1],
	    top_grid->lbuf[2]);
    printf("\n ubuf for top_grid \n");
    printf("\n %d, %d, %d \n", top_grid->ubuf[0],top_grid->ubuf[1],
	    top_grid->ubuf[2]);

    printf("\n global_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    pp_grid->Global_grid.L[0],pp_grid->Global_grid.U[0],pp_grid->Global_grid.L[1],pp_grid->Global_grid.U[1],
	    pp_grid->Global_grid.L[2],pp_grid->Global_grid.U[2]);
    printf("\n space is %10.5g, %10.5g, %10.5g \n", pp_grid->Global_grid.h[0], pp_grid->Global_grid.h[1],
	    pp_grid->Global_grid.h[2]);
    printf("\n gmax is %d, %d, %d\n", pp_grid->Global_grid.gmax[0],pp_grid->Global_grid.gmax[1],
	    pp_grid->Global_grid.gmax[2]);
    printf("\n Global mesh for top_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    pp_grid->Global_grid.GL[0],pp_grid->Global_grid.GU[0],pp_grid->Global_grid.GL[1],pp_grid->Global_grid.GU[1],
	    pp_grid->Global_grid.GL[2],pp_grid->Global_grid.GU[2]);
    printf("\n Virtual mesh for top_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    pp_grid->Global_grid.VL[0],pp_grid->Global_grid.VU[0],pp_grid->Global_grid.VL[1],pp_grid->Global_grid.VU[1],
	    pp_grid->Global_grid.VL[2],pp_grid->Global_grid.VU[2]);
    printf("\n lbuf for top_grid \n");
    printf("\n %d, %d, %d \n", pp_grid->Global_grid.lbuf[0],pp_grid->Global_grid.lbuf[1],
	    pp_grid->Global_grid.lbuf[2]);
    printf("\n ubuf for top_grid \n");
    printf("\n %d, %d, %d \n", pp_grid->Global_grid.ubuf[0],pp_grid->Global_grid.ubuf[1],
	    pp_grid->Global_grid.ubuf[2]);

    printf("\n Zoom_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    pp_grid->Zoom_grid.L[0],pp_grid->Zoom_grid.U[0],pp_grid->Zoom_grid.L[1],pp_grid->Zoom_grid.U[1],
	    pp_grid->Zoom_grid.L[2],pp_grid->Zoom_grid.U[2]);
    printf("\n space is %10.5g, %10.5g, %10.5g \n", pp_grid->Zoom_grid.h[0], pp_grid->Zoom_grid.h[1],
	    pp_grid->Zoom_grid.h[2]);
    printf("\n gmax is %d, %d, %d\n", pp_grid->Zoom_grid.gmax[0],pp_grid->Zoom_grid.gmax[1],
	    pp_grid->Zoom_grid.gmax[2]);
    printf("\n Zoom mesh for top_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    pp_grid->Zoom_grid.GL[0],pp_grid->Zoom_grid.GU[0],pp_grid->Zoom_grid.GL[1],pp_grid->Zoom_grid.GU[1],
	    pp_grid->Zoom_grid.GL[2],pp_grid->Zoom_grid.GU[2]);
    printf("\n Virtual mesh for top_grid:\n");
    printf("\n [%10.5g, %10.5g] * [%10.5g, %10.5g] * [%10.5g, %10.5g] \n",
	    pp_grid->Zoom_grid.VL[0],pp_grid->Zoom_grid.VU[0],pp_grid->Zoom_grid.VL[1],pp_grid->Zoom_grid.VU[1],
	    pp_grid->Zoom_grid.VL[2],pp_grid->Zoom_grid.VU[2]);
    printf("\n lbuf for top_grid \n");
    printf("\n %d, %d, %d \n", pp_grid->Zoom_grid.lbuf[0],pp_grid->Zoom_grid.lbuf[1],
	    pp_grid->Zoom_grid.lbuf[2]);
    printf("\n ubuf for top_grid \n");
    printf("\n %d, %d, %d \n", pp_grid->Zoom_grid.ubuf[0],pp_grid->Zoom_grid.ubuf[1],
	    pp_grid->Zoom_grid.ubuf[2]);


    printf("\n pp_grid:\n");
    printf("\n total number of nodes: %d\n", pp_grid->nn);
    printf("\n pp_gmax is %d, %d, %d\n", pp_grid->gmax[0], pp_grid->gmax[1], pp_grid->gmax[2]);
    printf("\n buffer extension is %d, %d, %d\n", pp_grid->buf[0], pp_grid->buf[1], pp_grid->buf[2]);


    int num_nodes = pp_numnodes();
    int myid = pp_mynode();

    printf("\nNumber of nodes : %d\n", num_nodes);
    printf("\nMy id number: %d\n", myid);
    */

}

void Incompress_Solver_Smooth_Basis::PrintMesh(void)
{
    int i,j,k,index;
    double coords[MAXD];
    switch (dim)
    {
    case 2:
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    printf("i = %d, j = %d, index = %d, I = %d, rho = %f, mu = %f\n"
		    ,i,j,index,ij_to_I[i][j],
		    cell_center[index].m_state.m_rho
		    ,cell_center[index].m_state.m_mu);
	}
	break;
    case 3:
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    printf("i = %d, j = %d, k = %d, index = %d, I = %d, rho = %f, mu = %f\n"
		    ,i,j,k,index,ijk_to_I[i][j][k],
		    cell_center[index].m_state.m_rho,
		    cell_center[index].m_state.m_mu);
	}
	break;
    }
    
}

void Incompress_Solver_Smooth_Basis::initMesh(void)
{
	int i,j,k, index;
	double coords[MAXD];
	int num_cells;
	int icoords[MAXD];
	int cell_index;

	// init cell_center
	L_RECTANGLE       rectangle;

	FT_MakeGridIntfc(front); 
	//Calculate the crossings ,etc. assuming that the interface is topological correct
	setDomain();

	num_cells = 1;
	for (i = 0; i < dim; ++i)
	{
	    num_cells *= (top_gmax[i] + 1);
	}
	cell_center.insert(cell_center.end(),num_cells,rectangle);
	
	// setup vertices
	// left to right, down to up
	switch (dim)
	{
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
		index = d_index2d(i,j,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;

	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
	    	coords[2] = top_L[2] + top_h[2]*k;
		index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    	cell_center[index].icoords[2] = k;
	    }
	}

	setCornerMesh();
	//Set cell_corner vector & cell_center[i].corner index and coordinates
	
	setCornerComp(); //Cell corner components
	setComponent();  //Cell center components
	FT_FreeGridIntfc(front);
}

void Incompress_Solver_Smooth_Basis::setComponent(void)
{
	int i;
	double coords[MAXD];
	
	// cell center components
	for (i = 0; i < cell_center.size(); i++)
	{
	    cell_center[i].comp = 
	    		getComponent(cell_center[i].icoords);
	}
}

void Incompress_Solver_Smooth_Basis::setIndexMap(void)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];

	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 2:
	    	FT_MatrixMemoryAlloc((POINTER*)&ij_to_I,top_gmax[0]+1,
				top_gmax[1]+1,INT);
	    	break;
	    case 3:
	    	FT_TriArrayMemoryAlloc((POINTER*)&ijk_to_I,top_gmax[0]+1,
				top_gmax[1]+1,top_gmax[2]+1,INT);
	    	break;
	    }
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[i] : 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[i] : 1;
	}
	switch (dim)
	{
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
		    ij_to_I[i][j] = -1;
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
                if (cell_center[ic].comp != SOLID_COMP)
                {
                    ij_to_I[i][j] = index + ilower;
                    index++;
                }
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
		    ijk_to_I[i][j][k] = -1;
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp != SOLID_COMP)
		{
                    ijk_to_I[i][j][k] = index + ilower;
                    index++;
                }
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
}	/* end setIndexMap */

void Incompress_Solver_Smooth_Basis::computeExactSolution(double *coords, L_STATE &state) 
{
	state.setZero();
}

void Incompress_Solver_Smooth_Basis::getVelocity(double *p, double *U)
{
	int i;
	double **vel = field->vel;

	if (!FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[0],getStateXvel,&U[0]))
	{
	    for (i = 0; i < dim; ++i) U[i] = 0.0;
	    return;
	}
	FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[1],getStateYvel,&U[1]);
	if (dim == 3)
	    FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[2],getStateZvel,&U[2]);
}

void Incompress_Solver_Smooth_Basis::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void Incompress_Solver_Smooth_Basis::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int Incompress_Solver_Smooth_Basis::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void Incompress_Solver_Smooth_Basis::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].m_coords[i];
}

double Incompress_Solver_Smooth_Basis::getDistance(double *c0, double *c1)
{
	int i;
	double dist;
	dist = 0.0;
	for (i = 0; i < dim; ++i)
	    dist += sqr(c0[i] - c1[i]);
	dist = sqrt(dist);
	return dist;
}


// input : p[]
// output: q[]

void Incompress_Solver_Smooth_Basis::getNearestInterfacePoint(
	COMPONENT comp,
	double *p, 
	double *q,
	double *nor,		// Normal vector
	double *kappa)		// curvature
{
	INTERFACE *intfc = front->interf;
	int i,j,dim = front->rect_grid->dim;
	double t[MAXD],mag_nor;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	POINT *pts[MAXD];
	BOND *bond;
	TRI *tri;
	static double **pts_nor,*pts_kappa;

	if (pts_kappa == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&pts_kappa,MAXD,FLOAT);
	    FT_MatrixMemoryAlloc((POINTER*)&pts_nor,MAXD,MAXD,FLOAT);
	}
	nearest_interface_point(p,comp,intfc,NO_BOUNDARIES,NULL,q,t,&hse,&hs);
	if (hse != NULL)
	{
	    switch (dim)
	    {
	    case 2:
		bond = Bond_of_hse(hse);
		pts[0] = bond->start;
		pts[1] = bond->end;
		for (i = 0; i < dim; ++i)
		{
		    GetFrontCurvature(pts[i],hse,hs,&pts_kappa[i],front);
		    GetFrontNormal(pts[i],hse,hs,pts_nor[i],front);
		}
		*kappa = (1.0 - t[0])*pts_kappa[0] + t[0]*pts_kappa[1];
		for (i = 0; i < dim; ++i)
		    nor[i] = (1.0 - t[0])*pts_nor[0][i] + t[0]*pts_nor[1][i];
		mag_nor = mag_vector(nor,dim);
		for (i = 0; i < dim; ++i)
		    nor[i] /= mag_nor;
		break;
	    case 3:
		tri = Tri_of_hse(hse);	
		for (i = 0; i < dim; ++i)
		{
		    pts[i] = Point_of_tri(tri)[i];
		    GetFrontCurvature(pts[i],hse,hs,&pts_kappa[i],front);
		    GetFrontNormal(pts[i],hse,hs,pts_nor[i],front);
		}
		for (i = 0; i < dim; ++i)
		{
		    *kappa = t[i]*pts_kappa[i];
		    for (j = 0; j < dim; ++j)
			nor[j] = t[i]*pts_nor[i][j];
		}
		mag_nor = mag_vector(nor,dim);
		for (i = 0; i < dim; ++i)
		    nor[i] /= mag_nor;
	    }
	}
	else
	{
	    *kappa = 0.0;
	    nor[0] = 1.0;
	    for (i = 1; i < dim; ++i) nor[i] = 0.0;
	}
}	/* end getNearestInterfacePoint */

int Incompress_Solver_Smooth_Basis::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int Incompress_Solver_Smooth_Basis::getComponent(
	int *icoords)
{
	int index;
	switch (dim)
	{
	case 2:
	    index = d_index2d(icoords[0],icoords[1],top_gmax);
	    return top_comp[index];
	case 3:
	    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
	    return top_comp[index];
	}
}

void Incompress_Solver_Smooth_Basis::save(char *filename)
{
	
	RECT_GRID *rect_grid = front->rect_grid;
	INTERFACE *intfc    = front->interf;
		
	int i, j;
	int xmax = rect_grid->gmax[0];
	int ymax = rect_grid->gmax[1];
	double x, y;
	
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in "
		       "SaveAsTecplot_rect_grid_and_interface().", filename);
		exit(0);
	}
	
	// secondly print out the interface
		
	if(exists_interface(intfc))
	{
	    CURVE		**cur;
	    CURVE		*curve;
	    BOND		*bond;
			
	    for(cur=intfc->curves; cur && *cur; cur++)
	    {
		curve = *cur;
		fprintf(hfile, "ZONE I=%d J=%d F=POINT \n", 
				curve->num_points, 1);
		bond=curve->first;
		fprintf(hfile, "%.4f %.4f \n",bond->start->_coords[0], 
				bond->start->_coords[1]);
		for(bond=curve->first; bond!=NULL; bond=bond->next)
		    fprintf(hfile, "%.4f %.4f \n",bond->end->_coords[0], 
		    		bond->end->_coords[1]);
		}					
	}		
	fclose(hfile);
}



void Incompress_Solver_Smooth_Basis::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i,size;

	comp_grid = front->rect_grid;
	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	comp_grid = computational_grid(grid_intfc);
	pp_grid = front->pp_grid;

	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;
	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	iFparams = (IF_PARAMS*)front->extra1;
	
	hmin = top_h[0];
	size = top_gmax[0]+1;
        for (i = 1; i < dim; ++i)
	{
            if (hmin > top_h[i]) hmin = top_h[i];
	    size *= (top_gmax[i]+1);
	}

	switch (dim)
	{
	case 2:
	    if (first)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(IF_FIELD));
		iFparams->field = field;
	    	FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&source,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&diff_coeff,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->pres,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->vort,size,
					sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&field->vel,2,size,
					sizeof(double));
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    break;
	case 3:
	    if (first)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(IF_FIELD));
		iFparams->field = field;
	    	FT_VectorMemoryAlloc((POINTER*)&array,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->pres,size,
					sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&field->vel,3,size,
					sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&field->vort3d,3,size,
					sizeof(double));
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    break;
	}
}

void Incompress_Solver_Smooth_Basis::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front);
}


bool Incompress_Solver_Smooth_Basis::interiorTri(TRI *tri)
{
    double center[3];

    center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0])/3.0;
    center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1])/3.0;
    center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2])/3.0;

    if ( 
	 ( (front->interf->rect_bdry_type[0][0] != PERIODIC_BOUNDARY) || (center[0] >= x_pp_bdry[x_pp_index])  ) && 
	 ( (front->interf->rect_bdry_type[0][1] != PERIODIC_BOUNDARY) || (center[0] < x_pp_bdry[x_pp_index+1]) ) &&
	 ( (front->interf->rect_bdry_type[1][0] != PERIODIC_BOUNDARY) || (center[1] >= y_pp_bdry[y_pp_index])  ) &&
	 ( (front->interf->rect_bdry_type[1][1] != PERIODIC_BOUNDARY) || (center[1] < y_pp_bdry[y_pp_index+1]) ) &&
	 ( (front->interf->rect_bdry_type[2][0] != PERIODIC_BOUNDARY) || (center[2] >= z_pp_bdry[z_pp_index])  ) &&
	 ( (front->interf->rect_bdry_type[2][1] != PERIODIC_BOUNDARY) || (center[2] < z_pp_bdry[z_pp_index+1]) ) 
       )
	return true;
    else
	return false;

}


void Incompress_Solver_Smooth_Basis::recvTrisFrom2(int dim, int dir, int rec_id)
{
    int *G = pp_grid->gmax;
    int i;
    int tag = 1;
    double *tris_pos0_rec;
    double *tris_pos1_rec;
    double *tris_pos2_rec;

    MPI_Status stat;

    int *tris_global_index_rec;
    int *num_tris_rec = new int;

    INTERFACE* intfc = tempfront->interf;

    SURFACE **s;
    TRI *tri;

    MPI_Recv(num_tris_rec, 1, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

    tris_pos0_rec = new double [(*num_tris_rec)];
    tris_pos1_rec = new double [(*num_tris_rec)];
    tris_pos2_rec = new double [(*num_tris_rec)];
    tris_global_index_rec = new int [(*num_tris_rec)];

    MPI_Recv(tris_pos0_rec, (*num_tris_rec), MPI_DOUBLE, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(tris_pos1_rec, (*num_tris_rec), MPI_DOUBLE, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(tris_pos2_rec, (*num_tris_rec), MPI_DOUBLE, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(tris_global_index_rec, (*num_tris_rec), MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);


	for (s = intfc->surfaces; s && *s; ++s)
	{
	    for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	    {
		if (tri->global_index == -1)
    		{
		    double center[3];
		    double p00, p10, p20;
		    double p01, p11, p21;
		    double p02, p12, p22;
		    double T;
		    if( (x_pp_index == 0) && (dir == 0) && (dim == 0) )
		    {
			T = tempfront->rect_grid->GU[0] - tempfront->rect_grid->GL[0];
			p00 = Coords(Point_of_tri(tri)[0])[0] + T;
			p10 = Coords(Point_of_tri(tri)[1])[0] + T;
			p20 = Coords(Point_of_tri(tri)[2])[0] + T;
			center[0] = (p00 + p10 + p20);
		    }
		    else
			center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);

		    if( (y_pp_index == 0) && (dir == 0) && (dim == 1) )
		    {
			T = tempfront->rect_grid->GU[1] - tempfront->rect_grid->GL[1];
			p01 = Coords(Point_of_tri(tri)[0])[1] + T;
			p11 = Coords(Point_of_tri(tri)[1])[1] + T;
			p21 = Coords(Point_of_tri(tri)[2])[1] + T;
			center[1] = (p01 + p11 + p21);
		    }
		    else
			center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);

		    if( (z_pp_index == 0) && (dir == 0) && (dim == 2) )
		    {
			T = tempfront->rect_grid->GU[2] - tempfront->rect_grid->GL[2];
			p02 = Coords(Point_of_tri(tri)[0])[2] + T;
			p12 = Coords(Point_of_tri(tri)[1])[2] + T;
			p22 = Coords(Point_of_tri(tri)[2])[2] + T;
			center[2] = (p02 + p12 + p22);
		    }
		    else
			center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);
		    /*
			center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);
			center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);
			center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);
			*/
		    for (int i = 0; i < (*num_tris_rec); i++)
		    {
			if ( (center[0] == tris_pos0_rec[i]) && (center[1] == tris_pos1_rec[i]) && (center[2] == tris_pos2_rec[i]) )
			//if ( (fabs(center[0] - tris_pos0_rec[i]) < 1e-15) && (fabs(center[1] - tris_pos1_rec[i]) < 1e-15) && (fabs(center[2] - tris_pos2_rec[i]) < 1e-15) )
			    tri->global_index = tris_global_index_rec[i];
		    }
		}
	    }
	}

    delete[] tris_pos0_rec;
    delete[] tris_pos1_rec;
    delete[] tris_pos2_rec;
    delete[] tris_global_index_rec;
    delete num_tris_rec;

    tris_pos0_rec = 0;
    tris_pos1_rec = 0;
    tris_pos2_rec = 0;
    tris_global_index_rec = 0;
    num_tris_rec = 0;
}


void Incompress_Solver_Smooth_Basis::sendTrisTo(int dim, int dir, int dst_id, int *&num_tris, double *&tris0, double *&tris1, double *&tris2, int *&tris_index)
{
    int *G = pp_grid->gmax;
    int i;
    int tag = 1;

    MPI_Request request_num;
    MPI_Request request_pos0;
    MPI_Request request_pos1;
    MPI_Request request_pos2;
    MPI_Request request_index;


    MPI_Isend(num_tris, 1, MPI_INT, dst_id, tag, MPI_COMM_WORLD, &request_num);
    
    MPI_Isend(tris0, (*num_tris), MPI_DOUBLE, dst_id, tag, MPI_COMM_WORLD, &request_pos0);
    MPI_Isend(tris1, (*num_tris), MPI_DOUBLE, dst_id, tag, MPI_COMM_WORLD, &request_pos1);
    MPI_Isend(tris2, (*num_tris), MPI_DOUBLE, dst_id, tag, MPI_COMM_WORLD, &request_pos2);
    MPI_Isend(tris_index, (*num_tris), MPI_INT, dst_id, tag, MPI_COMM_WORLD, &request_index);
}


void Incompress_Solver_Smooth_Basis::recvTrisFrom(int dim, int dir, int rec_id)
{
    int *G = pp_grid->gmax;
    int i;
    int tag = 1;
    double *tris_pos0_rec;
    double *tris_pos1_rec;
    double *tris_pos2_rec;

    MPI_Status stat;

    int *tris_global_index_rec;
    int *num_tris_rec = new int;

    INTERFACE* intfc = tempfront->interf;

    SURFACE **s;
    TRI *tri;

    MPI_Recv(num_tris_rec, 1, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);


    tris_pos0_rec = new double [(*num_tris_rec)];
    tris_pos1_rec = new double [(*num_tris_rec)];
    tris_pos2_rec = new double [(*num_tris_rec)];
    tris_global_index_rec = new int [(*num_tris_rec)];


    MPI_Recv(tris_pos0_rec, (*num_tris_rec), MPI_DOUBLE, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(tris_pos1_rec, (*num_tris_rec), MPI_DOUBLE, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(tris_pos2_rec, (*num_tris_rec), MPI_DOUBLE, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(tris_global_index_rec, (*num_tris_rec), MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	    {
		if (tri->global_index == -1)
    		{
		    double center[3];
		    double p00, p10, p20;
		    double p01, p11, p21;
		    double p02, p12, p22;
		    double T;
		    if( (x_pp_index == G[0] - 1) && (dir == 1) && (dim == 0) )
		    {
			T = tempfront->rect_grid->GL[0] - tempfront->rect_grid->GU[0];
			p00 = Coords(Point_of_tri(tri)[0])[0] + T;
			p10 = Coords(Point_of_tri(tri)[1])[0] + T;
			p20 = Coords(Point_of_tri(tri)[2])[0] + T;
			center[0] = (p00 + p10 + p20);
		    }
		    else
			center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);

		    if( (y_pp_index == G[1] - 1) && (dir == 1) && (dim == 1) )
		    {
			T = tempfront->rect_grid->GL[1] - tempfront->rect_grid->GU[1];
			p01 = Coords(Point_of_tri(tri)[0])[1] + T;
			p11 = Coords(Point_of_tri(tri)[1])[1] + T;
			p21 = Coords(Point_of_tri(tri)[2])[1] + T;
			center[1] = (p01 + p11 + p21);
		    }
		    else
			center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);

		    if( (z_pp_index == G[2] - 1) && (dir == 1) && (dim == 2) )
		    {
			T = tempfront->rect_grid->GL[2] - tempfront->rect_grid->GU[2];
			p02 = Coords(Point_of_tri(tri)[0])[2] + T;
			p12 = Coords(Point_of_tri(tri)[1])[2] + T;
			p22 = Coords(Point_of_tri(tri)[2])[2] + T;
			center[2] = (p02 + p12 + p22);
		    }
		    else
			center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);

		    for (int i = 0; i < (*num_tris_rec); i++)
		    {
			if ( (center[0] == tris_pos0_rec[i]) && (center[1] == tris_pos1_rec[i]) && (center[2] == tris_pos2_rec[i]) )
			    tri->global_index = tris_global_index_rec[i];
		    }
		}
	    }
	}


    delete[] tris_pos0_rec;
    delete[] tris_pos1_rec;
    delete[] tris_pos2_rec;
    delete[] tris_global_index_rec;
    delete num_tris_rec;

    tris_pos0_rec = 0;
    tris_pos1_rec = 0;
    tris_pos2_rec = 0;
    tris_global_index_rec = 0;
    num_tris_rec = 0;

}


void Incompress_Solver_Smooth_Basis::checkUnsetGlobalIndex(void)
{
    printf("x_pp_index = %d, y_pp_index = %d, z_pp_index = %d\n", x_pp_index, y_pp_index, z_pp_index);
    printf("lbound = %10.8g ubound = %10.8g\n", x_pp_bdry[x_pp_index], x_pp_bdry[x_pp_index+1]);

    INTERFACE *intfc = tempfront->interf;
    SURFACE **s;
    TRI *tri;
    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{

		if(tri->global_index == -1)
		{
		    double centerp[3];
		    centerp[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);
		    centerp[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);
		    centerp[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);
		    printf("\nGlobal index is still -1! local index = %d", tri->local_index);
		    printf("\ncenter = (%10.8g, %10.8g, %10.8g)\n", centerp[0], centerp[1], centerp[2]);
		    //clean_up(ERROR);
		}

	}
    }

}

void Incompress_Solver_Smooth_Basis::commGlobalTriIndex2()
{
    int i,j,k;
    int bdry_type[MAXD][2];
    MPI_Request request;
    MPI_Status status;

    int *G = pp_grid->gmax;

    int me[MAXD], him[MAXD];

    int my_id, dst_id, rcv_id;

    INTERFACE *intfc = tempfront->interf;

    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < 2; j++)
	{
	    bdry_type[i][j] = rect_boundary_type(intfc,i,j);

	}
    }

    my_id = pp_mynode();
    find_Cartesian_coordinates(my_id, pp_grid, me);

    for (i = 0; i < dim; i++)
    {
	for(j = 0; j < 2; j++)
	{

	    int *tris_num_send = 0;

	    double *tris_send0 = 0;
	    double *tris_send1 = 0;
	    double *tris_send2 = 0;
	    int *tris_send_index = 0;


	    for (k = 0; k < dim; k++)
		him[k] = me[k];

	    pp_gsync();
	    if (bdry_type[i][j] == SUBDOMAIN_BOUNDARY)
	    {
		makeTrisForSend2(i, j, tris_num_send, tris_send0, tris_send1, tris_send2, tris_send_index);

		him[i] = me[i] + 2*j - 1;
		him[i] = (him[i] + G[i]) % G[i];
		dst_id = domain_id(him,G,dim);
		sendTrisTo(i, j, dst_id, tris_num_send, tris_send0, tris_send1, tris_send2, tris_send_index);
	    }

	    if (bdry_type[i][((j+1)%2)] == SUBDOMAIN_BOUNDARY)
	    {
		him[i] = me[i] - 2*j + 1;
		him[i] = (him[i] + G[i]) % G[i];
		rcv_id = domain_id(him,G,dim);
		recvTrisFrom2(i, ((j+1)%2), rcv_id);
	    }

	    pp_gsync();

	    delete tris_num_send;

	    delete[] tris_send0;
	    delete[] tris_send1;
	    delete[] tris_send2;
	    delete[] tris_send_index;
	}
    }

    printf("\nSecond round finished\n");
/*
    SURFACE **s;
    TRI *tri;
    int nn = 1;
    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{

		if(tri->global_index == -1)
		{
		    double centerp[3];
		    centerp[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);
		    centerp[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);
		    centerp[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);
		    printf("\nGlobal index is still -1! local index = %d", tri->local_index);
		    printf("\ncenter = (%10.8g, %10.8g, %10.8g)\n", centerp[0], centerp[1], centerp[2]);
		    clean_up(ERROR);
		}

	}
    }
*/
}


void Incompress_Solver_Smooth_Basis::commGlobalTriIndex()
{
    int i,j,k;
    int bdry_type[MAXD][2];
    MPI_Request request;
    MPI_Status status;

    int *G = pp_grid->gmax;

    int me[MAXD], him[MAXD];

    int my_id, dst_id, rcv_id;

    INTERFACE *intfc = tempfront->interf;

    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < 2; j++)
	{
	    bdry_type[i][j] = rect_boundary_type(intfc,i,j);

	}
    }

    my_id = pp_mynode();
    find_Cartesian_coordinates(my_id, pp_grid, me);


    for (i = 0; i < dim; i++)
    {
	for(j = 0; j < 2; j++)
	{


	    int *tris_num_send = 0;

	    double *tris_send0 = 0;
	    double *tris_send1 = 0;
	    double *tris_send2 = 0;
	    int *tris_send_index = 0;

	    for (k = 0; k < dim; k++)
		him[k] = me[k];

	    pp_gsync();

	    if (bdry_type[i][j] == SUBDOMAIN_BOUNDARY)
	    {
		makeTrisForSend(i, j, tris_num_send, tris_send0, tris_send1, tris_send2, tris_send_index);
		him[i] = me[i] + 2*j - 1;
		him[i] = (him[i] + G[i]) % G[i];
		dst_id = domain_id(him,G,dim);
		sendTrisTo(i, j, dst_id, tris_num_send, tris_send0, tris_send1, tris_send2, tris_send_index);
	    }

	    if (bdry_type[i][((j+1)%2)] == SUBDOMAIN_BOUNDARY)
	    {
		him[i] = me[i] - 2*j + 1;
		him[i] = (him[i] + G[i]) % G[i];
		rcv_id = domain_id(him,G,dim);
		recvTrisFrom(i, ((j+1)%2), rcv_id);

	    }

	    pp_gsync();
	    delete tris_num_send;

	    delete[] tris_send0;
	    delete[] tris_send1;
	    delete[] tris_send2;
	    delete[] tris_send_index;

	}
    }

    printf("\nFirst round finished\n");
/*
    SURFACE **s;
    TRI *tri;
    int nn = 1;
    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    if (!(tri->interior))
	    {
		printf("\nThe %d-th exterior tri", nn++);

		if(tri->global_index == -1)
		{
		    printf("\nGlobal index is still -1!\n");
		}

	    }
	}
    }
*/
}

void Incompress_Solver_Smooth_Basis::makeTrisForSend2(int dim, int dir, int *&num_tris, double *&tris0, double *&tris1, double *&tris2, int *&tris_index)
{
    int *G = pp_grid->gmax;
    int i;
    int num_tris_send;

    double lBound, uBound;

    std::vector <double> pos0_temp;
    std::vector <double> pos1_temp;
    std::vector <double> pos2_temp;
    std::vector <int> gindex_temp;
/*
    pos0_temp.reserve(Ntris_all);
    pos1_temp.reserve(Ntris_all);
    pos2_temp.reserve(Ntris_all);
    gindex_temp.reserve(Ntris_all);
*/
    INTERFACE* intfc = tempfront->interf;
    SURFACE **s;
    TRI *tri;

    switch(dim)
    {
	case 0:
	    if (dir == 0)
	    {
		lBound = 3.0*(x_pp_bdry[x_pp_index] - 6.0*top_h[0]);
		uBound = 3.0*(x_pp_bdry[x_pp_index] + 10.0*top_h[0]);

	    }
	    else if (dir == 1)
	    {
		lBound = 3.0*(x_pp_bdry[x_pp_index+1] - 10.0*top_h[0]);
		uBound = 3.0*(x_pp_bdry[x_pp_index+1] + 6.0*top_h[0]);
	    }
	    else
	    {
		printf("\nWrong Direction!\n");
		clean_up(ERROR);
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
		{
		    double center[3];
		    center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);
		    center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);
		    center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);
		    //if( (center[0] >= lBound) && (center[0] <= uBound) )
		    {
			double T;
			double p00, p10, p20;

			if( (x_pp_index == 0) && (dir == 0) )
			{
			    T = tempfront->rect_grid->GU[0] - tempfront->rect_grid->GL[0];
			    p00 = Coords(Point_of_tri(tri)[0])[0] + T;
			    p10 = Coords(Point_of_tri(tri)[1])[0] + T;
			    p20 = Coords(Point_of_tri(tri)[2])[0] + T;
			    center[0] = (p00 + p10 + p20);
			}
/*
			if( (x_pp_index == G[0] - 1) && (dir == 1) )
			{
			    T = front->rect_grid->GL[0] - front->rect_grid->GU[0];
			    p00 = Coords(Point_of_tri(tri)[0])[0] + T;
			    p10 = Coords(Point_of_tri(tri)[1])[0] + T;
			    p20 = Coords(Point_of_tri(tri)[2])[0] + T;
			    center[0] = (p00 + p10 + p20);
			}
*/

			pos0_temp.push_back(center[0]);
			pos1_temp.push_back(center[1]);
			pos2_temp.push_back(center[2]);
			gindex_temp.push_back(tri->global_index);
		    }
		}
	    }
	    break;
	case 1:
	    if (dir == 0)
	    {
		lBound = 3.0*(y_pp_bdry[y_pp_index] - 6.0*top_h[1]);
		uBound = 3.0*(y_pp_bdry[y_pp_index] + 10.0*top_h[1]);

	    }
	    else if (dir == 1)
	    {
		lBound = 3.0*(y_pp_bdry[y_pp_index+1] - 10.0*top_h[1]);
		uBound = 3.0*(y_pp_bdry[y_pp_index+1] + 6.0*top_h[1]);
	    }
	    else
	    {
		printf("\nWrong Direction!\n");
		clean_up(ERROR);
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
		{
		    double center[3];
		    center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);
		    center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);
		    center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);

		    //if( (center[1] >= lBound) && (center[1] <= uBound) )
		    {
			double T;
			double p01, p11, p21;
			
			if( (y_pp_index == 0) && (dir == 0) )
			{
			    T = tempfront->rect_grid->GU[1] - tempfront->rect_grid->GL[1];
			    p01 = Coords(Point_of_tri(tri)[0])[1] + T;
			    p11 = Coords(Point_of_tri(tri)[1])[1] + T;
			    p21 = Coords(Point_of_tri(tri)[2])[1] + T;
			    center[1] = (p01 + p11 + p21);
			}
			/*
			if( (y_pp_index == G[1] - 1) && (dir == 1) )
			{
			    T = front->rect_grid->GL[1] - front->rect_grid->GU[1];
			    p01 = Coords(Point_of_tri(tri)[0])[1] + T;
			    p11 = Coords(Point_of_tri(tri)[1])[1] + T;
			    p21 = Coords(Point_of_tri(tri)[2])[1] + T;
			    center[1] = (p01 + p11 + p21);
			}
			*/

			pos0_temp.push_back(center[0]);
			pos1_temp.push_back(center[1]);
			pos2_temp.push_back(center[2]);
			gindex_temp.push_back(tri->global_index);
		    }
		}
	    }
	    break;
	case 2:
	    if (dir == 0)
	    {
		lBound = 3.0*(z_pp_bdry[z_pp_index] - 6.0*top_h[2]);
		uBound = 3.0*(z_pp_bdry[z_pp_index] + 10.0*top_h[2]);

	    }
	    else if (dir == 1)
	    {
		lBound = 3.0*(z_pp_bdry[z_pp_index+1] - 10.0*top_h[2]);
		uBound = 3.0*(z_pp_bdry[z_pp_index+1] + 6.0*top_h[2]);
	    }
	    else
	    {
		printf("\nWrong Direction!\n");
		clean_up(ERROR);
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
		{
		    double center[3];
		    center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);
		    center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);
		    center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);

		    //if( (center[2] >= lBound) && (center[2] <= uBound) )
		    {
			double T;
			double p02, p12, p22;
			if( (z_pp_index == 0) && (dir == 0) )
			{
			    T = tempfront->rect_grid->GU[2] - tempfront->rect_grid->GL[2];
			    p02 = Coords(Point_of_tri(tri)[0])[2] + T;
			    p12 = Coords(Point_of_tri(tri)[1])[2] + T;
			    p22 = Coords(Point_of_tri(tri)[2])[2] + T;
			    center[2] = (p02 + p12 + p22);
			}
			/*
			if( (z_pp_index == G[2] - 1) && (dir == 1) )
			{
			    T = front->rect_grid->GL[2] - front->rect_grid->GU[2];
			    p02 = Coords(Point_of_tri(tri)[0])[2] + T;
			    p12 = Coords(Point_of_tri(tri)[1])[2] + T;
			    p22 = Coords(Point_of_tri(tri)[2])[2] + T;
			    center[2] = (p02 + p12 + p22);
			}
			*/
		
			pos0_temp.push_back(center[0]);
			pos1_temp.push_back(center[1]);
			pos2_temp.push_back(center[2]);
			gindex_temp.push_back(tri->global_index);
		    }
		}
	    }
	    break;
	default:
	    printf("\nWrong Dimension!\n");
	    clean_up(ERROR);
	    break;
    }

    
    num_tris_send = (int) gindex_temp.size();

    num_tris = new int (num_tris_send);
    tris0 = new double [num_tris_send];
    tris1 = new double [num_tris_send];
    tris2 = new double [num_tris_send];

    tris_index = new int [num_tris_send];


    for (i = 0; i < num_tris_send; i++)
    {
	tris0[i] = pos0_temp[i];
	tris1[i] = pos1_temp[i];
	tris2[i] = pos2_temp[i];

	tris_index[i] = gindex_temp[i];
    }

    pos0_temp.clear();
    pos1_temp.clear();
    pos2_temp.clear();
    gindex_temp.clear();

}

void Incompress_Solver_Smooth_Basis::makeTrisForSend(int dim, int dir, int *&num_tris, double *&tris0, double *&tris1, double *&tris2, int *&tris_index)
{
    int *G = pp_grid->gmax;
    int i;
    int num_tris_send;

    double lBound, uBound;

    std::vector <double> pos0_temp;
    std::vector <double> pos1_temp;
    std::vector <double> pos2_temp;
    std::vector <int> gindex_temp;
/*
    pos0_temp.reserve(Ntris_all);
    pos1_temp.reserve(Ntris_all);
    pos2_temp.reserve(Ntris_all);
    gindex_temp.reserve(Ntris_all);
*/
    INTERFACE* intfc = tempfront->interf;
    SURFACE **s;
    TRI *tri;

    switch(dim)
    {
	case 0:
	    if (dir == 0)
	    {
		lBound = 3.0*(x_pp_bdry[x_pp_index] - 6.0*top_h[0]);
		uBound = 3.0*(x_pp_bdry[x_pp_index] + 10.0*top_h[0]);

	    }
	    else if (dir == 1)
	    {
		lBound = 3.0*(x_pp_bdry[x_pp_index+1] - 10.0*top_h[0]);
		uBound = 3.0*(x_pp_bdry[x_pp_index+1] + 6.0*top_h[0]);
	    }
	    else
	    {
		printf("\nWrong Direction!\n");
		clean_up(ERROR);
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
		{
		    double center[3];
		    center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);
		    center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);
		    center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);
		    //if( (center[0] >= lBound) && (center[0] <= uBound) )
		    {
			double T;
			double p00, p10, p20;
/*
			if( (x_pp_index == 0) && (dir == 0) )
			{
			    T = front->rect_grid->GU[0] - front->rect_grid->GL[0];
			    p00 = Coords(Point_of_tri(tri)[0])[0] + T;
			    p10 = Coords(Point_of_tri(tri)[1])[0] + T;
			    p20 = Coords(Point_of_tri(tri)[2])[0] + T;
			    center[0] = (p00 + p10 + p20);
			}
*/
			if( (x_pp_index == G[0] - 1) && (dir == 1) )
			{
			    T = tempfront->rect_grid->GL[0] - tempfront->rect_grid->GU[0];
			    p00 = Coords(Point_of_tri(tri)[0])[0] + T;
			    p10 = Coords(Point_of_tri(tri)[1])[0] + T;
			    p20 = Coords(Point_of_tri(tri)[2])[0] + T;
			    center[0] = (p00 + p10 + p20);
			}


			pos0_temp.push_back(center[0]);
			pos1_temp.push_back(center[1]);
			pos2_temp.push_back(center[2]);
			gindex_temp.push_back(tri->global_index);
		    }
		}
	    }
	    break;
	case 1:
	    if (dir == 0)
	    {
		lBound = 3.0*(y_pp_bdry[y_pp_index] - 6.0*top_h[1]);
		uBound = 3.0*(y_pp_bdry[y_pp_index] + 10.0*top_h[1]);

	    }
	    else if (dir == 1)
	    {
		lBound = 3.0*(y_pp_bdry[y_pp_index+1] - 10.0*top_h[1]);
		uBound = 3.0*(y_pp_bdry[y_pp_index+1] + 6.0*top_h[1]);
	    }
	    else
	    {
		printf("\nWrong Direction!\n");
		clean_up(ERROR);
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
		{
		    double center[3];
		    center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);
		    center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);
		    center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);

		    //if( (center[1] >= lBound) && (center[1] <= uBound) )
		    {
			double T;
			double p01, p11, p21;
		/*	
			if( (y_pp_index == 0) && (dir == 0) )
			{
			    T = front->rect_grid->GU[1] - front->rect_grid->GL[1];
			    p01 = Coords(Point_of_tri(tri)[0])[1] + T;
			    p11 = Coords(Point_of_tri(tri)[1])[1] + T;
			    p21 = Coords(Point_of_tri(tri)[2])[1] + T;
			    center[1] = (p01 + p11 + p21);
			}
*/
			if( (y_pp_index == G[1] - 1) && (dir == 1) )
			{
			    T = tempfront->rect_grid->GL[1] - tempfront->rect_grid->GU[1];
			    p01 = Coords(Point_of_tri(tri)[0])[1] + T;
			    p11 = Coords(Point_of_tri(tri)[1])[1] + T;
			    p21 = Coords(Point_of_tri(tri)[2])[1] + T;
			    center[1] = (p01 + p11 + p21);
			}

			pos0_temp.push_back(center[0]);
			pos1_temp.push_back(center[1]);
			pos2_temp.push_back(center[2]);
			gindex_temp.push_back(tri->global_index);
		    }
		}
	    }
	    break;
	case 2:
	    if (dir == 0)
	    {
		lBound = 3.0*(z_pp_bdry[z_pp_index] - 6.0*top_h[2]);
		uBound = 3.0*(z_pp_bdry[z_pp_index] + 10.0*top_h[2]);

	    }
	    else if (dir == 1)
	    {
		lBound = 3.0*(z_pp_bdry[z_pp_index+1] - 10.0*top_h[2]);
		uBound = 3.0*(z_pp_bdry[z_pp_index+1] + 6.0*top_h[2]);
	    }
	    else
	    {
		printf("\nWrong Direction!\n");
		clean_up(ERROR);
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
		{
		    double center[3];
		    center[0] = (Coords(Point_of_tri(tri)[0])[0] + Coords(Point_of_tri(tri)[1])[0] + Coords(Point_of_tri(tri)[2])[0]);
		    center[1] = (Coords(Point_of_tri(tri)[0])[1] + Coords(Point_of_tri(tri)[1])[1] + Coords(Point_of_tri(tri)[2])[1]);
		    center[2] = (Coords(Point_of_tri(tri)[0])[2] + Coords(Point_of_tri(tri)[1])[2] + Coords(Point_of_tri(tri)[2])[2]);

		    //if( (center[2] >= lBound) && (center[2] <= uBound) )
		    {
			double T;
			double p02, p12, p22;
/*
			if( (z_pp_index == 0) && (dir == 0) )
			{
			    T = front->rect_grid->GU[2] - front->rect_grid->GL[2];
			    p02 = Coords(Point_of_tri(tri)[0])[2] + T;
			    p12 = Coords(Point_of_tri(tri)[1])[2] + T;
			    p22 = Coords(Point_of_tri(tri)[2])[2] + T;
			    center[2] = (p02 + p12 + p22);
			}
*/	
			if( (z_pp_index == G[2] - 1) && (dir == 1) )
			{
			    T = tempfront->rect_grid->GL[2] - tempfront->rect_grid->GU[2];
			    p02 = Coords(Point_of_tri(tri)[0])[2] + T;
			    p12 = Coords(Point_of_tri(tri)[1])[2] + T;
			    p22 = Coords(Point_of_tri(tri)[2])[2] + T;
			    center[2] = (p02 + p12 + p22);
			}
		
			pos0_temp.push_back(center[0]);
			pos1_temp.push_back(center[1]);
			pos2_temp.push_back(center[2]);
			gindex_temp.push_back(tri->global_index);
		    }
		}
	    }
	    break;
	default:
	    printf("\nWrong Dimension!\n");
	    clean_up(ERROR);
	    break;
    }

    
    num_tris_send = (int) gindex_temp.size();


    num_tris = new int (num_tris_send);
    tris0 = new double [num_tris_send];
    tris1 = new double [num_tris_send];
    tris2 = new double [num_tris_send];

    tris_index = new int [num_tris_send];


    for (i = 0; i < num_tris_send; i++)
    {
	tris0[i] = pos0_temp[i];
	tris1[i] = pos1_temp[i];
	tris2[i] = pos2_temp[i];

	tris_index[i] = gindex_temp[i];
    }

    pos0_temp.clear();
    pos1_temp.clear();
    pos2_temp.clear();
    gindex_temp.clear();

}


void Incompress_Solver_Smooth_Basis::checkTempFront(void)
{
    has_small_edge = 0;
    INTERFACE *intfc = tempfront->interf;

    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
    POINT *p0, *p1, *p2;

    SURFACE **s;
    TRI *tri;

    double p00, p01, p02;
    double p10, p11, p12;
    double p20, p21, p22;

    double L01, L02, L12;

    int i = 0;

    double tol = 0.02;

    printf ("\nhmin = %10.8g, tol for edge = %10.8g\n", hmin, tol*hmin);

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    p0 = Point_of_tri(tri)[0];
	    p1 = Point_of_tri(tri)[1];
	    p2 = Point_of_tri(tri)[2];

	    p00 = p0->_coords[0];
	    p01 = p0->_coords[1];
	    p02 = p0->_coords[2];

	    p10 = p1->_coords[0];
	    p11 = p1->_coords[1];
	    p12 = p1->_coords[2];

	    p20 = p2->_coords[0];
	    p21 = p2->_coords[1];
	    p22 = p2->_coords[2];

	    L01 = sqrt((p00 - p10)*(p00 - p10) + (p01 - p11)*(p01 - p11) + (p02 - p12)*(p02 - p12));
	    L02 = sqrt((p00 - p20)*(p00 - p20) + (p01 - p21)*(p01 - p21) + (p02 - p22)*(p02 - p22));
	    L12 = sqrt((p20 - p10)*(p20 - p10) + (p21 - p11)*(p21 - p11) + (p22 - p12)*(p22 - p12));

    	    if ( (L01 < tol*hmin) || (L02 < tol*hmin) || (L12 < tol*hmin) )
	    {
		has_small_edge = 1;
		printf("\n Small edge! Index of tri: %d, edge length L01 = %10.8g L02 = %10.8g L12 = %10.8g \n", i, L01, L02, L12);
		printf(" Point 1: (%10.8g, %10.8g, %10.8g)\n", p00, p01, p02);
		printf(" Point 2: (%10.8g, %10.8g, %10.8g)\n", p10, p11, p12);
		printf(" Point 3: (%10.8g, %10.8g, %10.8g)\n", p20, p21, p22);
	    }

	    
	    i++;
	}

	printf("\nNumber of tris in checkTempFront : %d\n", i);

    }


}


void Incompress_Solver_Smooth_Basis::setTempFrontAndRedist(void)
{

            tempfront = copy_front(front);/* f_copy_front called */
            tempfront->redis_flag = front->redis_flag;
            Interface_redistributed(tempfront) = NO;
            set_size_of_intfc_state(size_of_state(front->interf));
            set_copy_intfc_states(NO);
            tempfront->interf = pp_copy_interface(front->interf);
            if (tempfront->interf == NULL)
            {
                (void) printf("WARNING in use_redistribution, "
                              "unable to copy interface\n");
                clean_up(ERROR);
            }
	    interpolate_intfc_states(tempfront->interf) = YES;


	    checkTempFront();

	    pp_global_imax(&has_small_edge,1);

	    printf("\n has_small_edge = %d \n", has_small_edge);

	    if (has_small_edge == 1)
		FT_RedistMesh_temp(tempfront);
	
	    printf("\nAfter redistribution \n");
	    checkTempFront();

	    scatter_front(tempfront);

}

void Incompress_Solver_Smooth_Basis::setGlobalTriIndex()
{
    /*
    printf("\n Top grid: \n");
    printf("Dimension 0: [%10.8g, %10.8g]\n", top_grid->L[0], top_grid->U[0]);
    printf("Dimension 1: [%10.8g, %10.8g]\n", top_grid->L[1], top_grid->U[1]);
    printf("Dimension 2: [%10.8g, %10.8g]\n", top_grid->L[2], top_grid->U[2]);

        printf("\n Comp grid: \n");
    printf("Dimension 0: [%10.8g, %10.8g]\n", comp_grid->L[0], comp_grid->U[0]);
    printf("Dimension 1: [%10.8g, %10.8g]\n", comp_grid->L[1], comp_grid->U[1]);
    printf("Dimension 2: [%10.8g, %10.8g]\n", comp_grid->L[2], comp_grid->U[2]);

    printf("\n PP grid: \n");
    printf("Dimension 0: [%10.8g, %10.8g]\n", pp_grid->Global_grid.L[0], pp_grid->Global_grid.U[0]);
    printf("Dimension 1: [%10.8g, %10.8g]\n", pp_grid->Global_grid.L[1], pp_grid->Global_grid.U[1]);
    printf("Dimension 2: [%10.8g, %10.8g]\n", pp_grid->Global_grid.L[2], pp_grid->Global_grid.U[2]);
*/

    INTERFACE* intfc = tempfront->interf;

    SURFACE **s;
    TRI *tri;

    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
    POINT *p;

    int i;
    int num_nodes = pp_numnodes();
    int myid = pp_mynode();
    int   *ppgmax = tempfront->pp_grid->gmax;
    x_pp_index = myid % ppgmax[0];
    y_pp_index = ((myid-x_pp_index)/ppgmax[0]) % ppgmax[1];
    z_pp_index = (myid-x_pp_index-(ppgmax[0]*y_pp_index))/(ppgmax[0]*ppgmax[1]);


    static boolean first = YES;

    if (first)
    {
	first = NO;
	FT_VectorMemoryAlloc((POINTER*)&n_dist_tri,num_nodes,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&x_pp_bdry,ppgmax[0]+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&y_pp_bdry,ppgmax[1]+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&z_pp_bdry,ppgmax[2]+1,sizeof(double));

	for (i = 0; i <= ppgmax[0]; i++)
	    x_pp_bdry[i] = pp_grid->Global_grid.L[0] + i*(pp_grid->Global_grid.U[0] - pp_grid->Global_grid.L[0])/ppgmax[0];

	for (i = 0; i <= ppgmax[1]; i++)
	    y_pp_bdry[i] = pp_grid->Global_grid.L[1] + i*(pp_grid->Global_grid.U[1] - pp_grid->Global_grid.L[1])/ppgmax[1];

	for (i = 0; i <= ppgmax[2]; i++)
	    z_pp_bdry[i] = pp_grid->Global_grid.L[2] + i*(pp_grid->Global_grid.U[2] - pp_grid->Global_grid.L[2])/ppgmax[2];

    }
    Ntris = 0; 

    Ntris_all = 0;

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    tri->local_index = Ntris_all;
	    if (interiorTri(tri))
	    {
		tri->global_index = Ntris;
		tri->interior = TRUE;
		Ntris++;
	    }
	    else
	    {
		tri->global_index = -1;
		tri->interior = FALSE;
	    }
	    Ntris_all++;
	}
    }

    for (i = 0; i < num_nodes; ++i) 
	n_dist_tri[i] = 0;
    n_dist_tri[myid] = Ntris;

    pp_global_imax(n_dist_tri,num_nodes);

    tri_ilower = 0;
    tri_iupper = n_dist_tri[0];

    for (i = 1; i <= myid; i++)
    {
	tri_ilower += n_dist_tri[i-1];
	tri_iupper += n_dist_tri[i];
    }

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    if (tri->interior)
	    {
		tri->global_index = tri->global_index + tri_ilower;
	    }
	}
    }


    Npoints = 0;
    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    p = Point_of_tri(tri)[0];
	    p->local_index = -1;
	    p = Point_of_tri(tri)[1];
	    p->local_index = -1;
	    p = Point_of_tri(tri)[2];
	    p->local_index = -1;
	}
    }

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    p = Point_of_tri(tri)[0];
	    if (p->local_index == -1)
	    {
		p->local_index = Npoints;
		Npoints++;
	    }
	    p = Point_of_tri(tri)[1];
	    if (p->local_index == -1)
	    {
		p->local_index = Npoints;
		Npoints++;
	    }
	    p = Point_of_tri(tri)[2];
	    if (p->local_index == -1)
	    {
		p->local_index = Npoints;
		Npoints++;
	    }
	}
    }

    /*
    next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
	p->local_index = Npoints;
	if (p->local_index == 0 || p->local_index == 154)
	    printf("address of p is %p\n", p);
	Npoints++;
    }
    */
    printf("\n Num of triangles inside this processor: %d \n", Ntris);
    printf("\n Num of triangles for this processor: %d \n", Ntris_all);
    printf("Tri ilower ---> Tri iupper is %d ---> %d\n", tri_ilower, tri_iupper);

    int Ntris_global = Ntris;

    pp_global_isum(&Ntris_global, 1);

    printf("\n Num of triangles for the entire domain: %d \n", Ntris_global);
    printf("\n Num of points for this processor: %d \n", Npoints);

}

void Incompress_Solver_Smooth_Basis::setGlobalIndex()
{
	int i,j,k,ic;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(int));
	}
	NLblocks = 0;
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (cell_center[ic].comp == SOLID_COMP) continue;
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (cell_center[ic].comp == SOLID_COMP) continue;
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp == SOLID_COMP) continue;
		NLblocks++;
	    }
	    break;
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = NLblocks;
	pp_global_imax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];

        for (i = 1; i <= myid; i++)
        {
            ilower += n_dist[i-1];
            iupper += n_dist[i];
        }	
}

void Incompress_Solver_Smooth_Basis::printFrontInteriorStates(char *out_name, bool binary)
{
	int i,j,k,l,index;
	char filename[100];
	FILE *outfile;
        char dirname[256];

	double val[3];
	char str[100];
	int ival[8];
	char blank[2];
	sprintf(blank, "\n");

        sprintf(dirname,"%s/restart",out_name);
        if (pp_min_status(create_directory(dirname,NO)) == NO)
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }

	sprintf(filename,"%s/state.t%s",dirname,
			right_flush(front->step,7));
#if defined(__MPI__)
	if (pp_numnodes() > 1)
            sprintf(filename,"%s-p%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

	if (binary == YES)
	{
	    outfile = fopen(filename,"wb");

	    fluid_print_front_states(outfile,front, YES);

	    sprintf(str,"\nInterior ifluid states:\n");
	    fwrite(str, sizeof(char), 25, outfile);
	    switch (dim)
	    {
		case 2:
		    if(hardware_is_little_endian())
		    {
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			{
			    index = d_index2d(i,j,top_gmax);
			    val[0] = endian_double_swap(cell_center[index].m_state.m_rho);
			    val[1] = endian_double_swap(cell_center[index].m_state.m_P);
			    val[2] = endian_double_swap(cell_center[index].m_state.m_mu);
			    fwrite(val, sizeof(double), 3, outfile);
			    val[0] = endian_double_swap(cell_center[index].m_state.m_U[0]);
			    val[1] = endian_double_swap(cell_center[index].m_state.m_U[1]);
			    fwrite(val, sizeof(double), 2, outfile);
			    ival[0] = endian_int_swap(top_comp[index]);
			    fwrite(ival, sizeof(int), 1, outfile);
			}
		    }
		    else //big endian system
		    {
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			{
			    index = d_index2d(i,j,top_gmax);
			    val[0] = cell_center[index].m_state.m_rho;
			    val[1] = cell_center[index].m_state.m_P;
			    val[2] = cell_center[index].m_state.m_mu;
			    fwrite(val, sizeof(double), 3, outfile);
			    val[0] = cell_center[index].m_state.m_U[0];
			    val[1] = cell_center[index].m_state.m_U[1];
			    fwrite(val, sizeof(double), 2, outfile);
			    ival[0] = top_comp[index];
			    fwrite(ival, sizeof(int), 1, outfile);
			}
		    }
		    break;
		case 3:
		    if(hardware_is_little_endian())
		    {
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			for (k = 0; k <= top_gmax[2]; ++k)
			{
			    index = d_index3d(i,j,k,top_gmax);
			    val[0] = endian_double_swap(cell_center[index].m_state.m_rho);
			    val[1] = endian_double_swap(cell_center[index].m_state.m_P);
			    val[2] = endian_double_swap(cell_center[index].m_state.m_mu);
			    fwrite(val, sizeof(double), 3, outfile);
			    val[0] = endian_double_swap(cell_center[index].m_state.m_U[0]);
			    val[1] = endian_double_swap(cell_center[index].m_state.m_U[1]);
			    val[2] = endian_double_swap(cell_center[index].m_state.m_U[2]);
			    fwrite(val, sizeof(double), 3, outfile);
			    ival[0] = endian_int_swap(top_comp[index]);
			    fwrite(ival, sizeof(int), 1, outfile);

			}
		    }
		    else //big endian system
		    {
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			for (k = 0; k <= top_gmax[2]; ++k)
			{
			    index = d_index3d(i,j,k,top_gmax);
			    val[0] = cell_center[index].m_state.m_rho;
			    val[1] = cell_center[index].m_state.m_P;
			    val[2] = cell_center[index].m_state.m_mu;
			    fwrite(val, sizeof(double), 3, outfile);
			    val[0] = cell_center[index].m_state.m_U[0];
			    val[1] = cell_center[index].m_state.m_U[1];
			    val[2] = cell_center[index].m_state.m_U[2];
			    fwrite(val, sizeof(double), 3, outfile);
			    ival[0] = top_comp[index];
			    fwrite(ival, sizeof(int), 1, outfile);
			}
		    }
		    break;
	    }
	    fclose(outfile);
	}
	else//ASCII output
	{
	    outfile = fopen(filename,"w");

	    fluid_print_front_states(outfile,front,NO);

	    fprintf(outfile,"\nInterior ifluid states:\n");
	    switch (dim)
	    {
		case 2:
		    for (i = 0; i <= top_gmax[0]; ++i)
		    for (j = 0; j <= top_gmax[1]; ++j)
		    {
			index = d_index2d(i,j,top_gmax);
			fprintf(outfile,"%24.18g\n",
				cell_center[index].m_state.m_rho);
			fprintf(outfile,"%24.18g\n",
				cell_center[index].m_state.m_P);
			fprintf(outfile,"%24.18g\n",
				cell_center[index].m_state.m_mu);
			for (l = 0; l < dim; ++l)
			    fprintf(outfile,"%24.18g\n",
				    cell_center[index].m_state.m_U[l]);
			fprintf(outfile, "%d\n", top_comp[index]);
		    }
		    break;
		case 3:
		    for (i = 0; i <= top_gmax[0]; ++i)
		    for (j = 0; j <= top_gmax[1]; ++j)
		    for (k = 0; k <= top_gmax[2]; ++k)
		    {
			index = d_index3d(i,j,k,top_gmax);
			fprintf(outfile,"%24.18g\n",
				cell_center[index].m_state.m_rho);
			fprintf(outfile,"%24.18g\n",
				cell_center[index].m_state.m_P);
			fprintf(outfile,"%24.18g\n",
				cell_center[index].m_state.m_mu);
			for (l = 0; l < dim; ++l)
			    fprintf(outfile,"%24.18g\n",
				    cell_center[index].m_state.m_U[l]);
			fprintf(outfile, "%d\n", top_comp[index]);
		    }
	    }
	    fclose(outfile);
	}
}

void Incompress_Solver_Smooth_Basis::printFrontInteriorStatesRegridRep(char *out_name, bool binary, bool RegridRun)
{
	int i,j,k,l,index,temp;
        int index_i,index_j,index_k;
	char filename[100];
	FILE *outfile;
        char dirname[256];

        iFparams = (IF_PARAMS*)front->extra1;
        double **vel = iFparams->field->vel;
        double *pres = iFparams->field->pres;
        int rep_top_gmin[MAXD];
        int rep_top_gmax[MAXD];
        int top_gmax_regrid[MAXD];
        double coords[MAXD];

        double vel_out[MAXD];
        double pres_out;
        double mu_out;
        double rho_out;

	double val[3];
	char str[100];
	int ival[8];
	char blank[2];
	sprintf(blank, "\n");

        PP_GRID *pp_grid = front->pp_grid;
        PP_GRID *new_pp_grid = front->new_pp_grid;
        Front           *tempfront;
        RECT_GRID       dual_grid;
        RECT_GRID       *zoom_grid = &(new_pp_grid->Zoom_grid);
        char            fname[1024];
        double           L[MAXD], U[MAXD], *GL, *GU;
        int             ii,jj,kk;
        int             ratio[MAXD],ratio_regrid[MAXD];
        int             dim = front->rect_grid->dim;
        int             icrds_old[MAXD],icrds_new[MAXD];
        int             gmax[MAXD],lbuf[MAXD], ubuf[MAXD];

        GL = pp_grid->Global_grid.L;
        GU = pp_grid->Global_grid.U;

        find_Cartesian_coordinates(pp_mynode(),pp_grid,icrds_old);
        for (ii = 0; ii < dim; ++ii)
            ratio[ii] = (new_pp_grid->gmax[ii])/(pp_grid->gmax[ii]);

        for (ii = 0; ii < new_pp_grid->nn; ++ii)
        {
            find_Cartesian_coordinates(ii,new_pp_grid,icrds_new);

            kk = 1;
            for (jj = 0; jj < dim; ++jj)
            {
                kk = kk && ((ratio[jj]*icrds_old[jj] <= icrds_new[jj])
                            &&(icrds_new[jj] <= ratio[jj]*icrds_old[jj] +
                               ratio[jj]-1));
            }

            if (kk)
            {
               //set the zoom_grid of new_pp_grid here
                for (jj = 0; jj < dim; ++jj)
                {
                    L[jj] = new_pp_grid->dom[jj][icrds_new[jj]];
                    U[jj] = new_pp_grid->dom[jj][icrds_new[jj] + 1];
                    gmax[jj] = irint((U[jj] - L[jj])/(front->rect_grid->h[jj]));

                    switch (dim) // TODO Unify 2 and 3 D 
                    {
                    case 1:
                    case 2:
                        lbuf[jj] = (icrds_new[jj] > 0) ?
                                     new_pp_grid->buf[jj] : 0;
                        ubuf[jj] = (icrds_new[jj] <(new_pp_grid->gmax[jj]-1))
                                   ? new_pp_grid->buf[jj]:0;
                        break;
                    case 3:
                        lbuf[jj] = (icrds_new[jj] > 0)
                                ? new_pp_grid->buf[jj] :
                                front->rect_grid->lbuf[jj];

                        ubuf[jj] = (icrds_new[jj] <(new_pp_grid->gmax[jj]-1))
                                ? new_pp_grid->buf[jj] :
                                front->rect_grid->ubuf[jj];
                        break;
                    }
                }
                set_rect_grid(L,U,GL,GU,lbuf,ubuf,gmax,dim,
                    &front->rect_grid->Remap,&new_pp_grid->Zoom_grid);

                //clip the front with the new zoom_grid
                tempfront = copy_front(front);
                set_copy_intfc_states(YES);
                set_size_of_intfc_state(size_of_state(front->interf));
                tempfront->interf = copy_interface(front->interf);
                for (jj = 0; jj < dim; ++jj)
                {
                    if (icrds_new[jj] > 0)
                        rect_boundary_type(tempfront->interf,jj,0)
                                = SUBDOMAIN_BOUNDARY;
                    if (icrds_new[jj] <(new_pp_grid->gmax[jj]-1))
                        rect_boundary_type(tempfront->interf,jj,1)
                                = SUBDOMAIN_BOUNDARY;
                }

                //set topological grid here
                set_dual_grid(&dual_grid,zoom_grid);
                for (jj = 0; jj < dim; ++jj)
                    gmax[jj] = dual_grid.gmax[jj]+dual_grid.lbuf[jj]+
                                       dual_grid.ubuf[jj];
                set_rect_grid(dual_grid.VL,dual_grid.VU,dual_grid.GL,
                      dual_grid.GU,NOBUF,NOBUF,gmax,dim,&zoom_grid->Remap,
                                     &topological_grid(tempfront->interf));
                set_computational_grid(tempfront->interf,zoom_grid);

                //clip the interface to subdomain
                clip_front_for_output(tempfront,zoom_grid);

                //set relative coordinates of new nodes 
                for(jj = 0; jj < dim; ++jj)
                    icrds_new[jj] = icrds_new[jj] - icrds_old[jj]*ratio[jj];

                for (i = 0; i < dim; ++i)
                {
                    top_gmax_regrid[i] = top_gmax[i] + ((front->gmax_regrid[i] - front->gmax[i])/pp_grid->gmax[i]);
                    ratio_regrid[i] = front->gmax_regrid[i]/front->gmax[i];
                }

                rep_top_gmin[0]=(((top_gmax_regrid[0]-7)/ratio[0]))*(icrds_new[0]);
                rep_top_gmin[1]=(((top_gmax_regrid[1]-7)/ratio[1]))*(icrds_new[1]);
                rep_top_gmin[2]=(((top_gmax_regrid[2]-7)/ratio[2]))*(icrds_new[2]);
                rep_top_gmax[0]=(((top_gmax_regrid[0]-7)/ratio[0]))*(icrds_new[0]+1)+7;
                rep_top_gmax[1]=(((top_gmax_regrid[1]-7)/ratio[1]))*(icrds_new[1]+1)+7;
                rep_top_gmax[2]=(((top_gmax_regrid[2]-7)/ratio[2]))*(icrds_new[2]+1)+7;


                sprintf(dirname,"%s/restart",out_name);
                if (pp_min_status(create_directory(dirname,NO)) == NO)
                {
                    screen("Cannot create directory %s\n",dirname);
                    clean_up(ERROR);
                }

                sprintf(filename,"%s/state.t%s",dirname,
                                right_flush(front->step,7));
#if defined(__MPI__)
                if (pp_numnodes() > 1)
                    sprintf(filename,"%s-p%s",filename,right_flush(ii,4));
#endif /* defined(__MPI__) */

	    if (binary == YES)
	    {
	    outfile = fopen(filename,"wb");

            //fluid_print_front_states(outfile,front, YES);

            sprintf(str,"\nInterior ifluid states:\n");
            fwrite(str, sizeof(char), 25, outfile);
            switch (dim)
            {
		case 2:
		    if(hardware_is_little_endian())
		    {
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			{
			    index = d_index2d(i,j,top_gmax);
			    val[0] = endian_double_swap(cell_center[index].m_state.m_rho);
			    val[1] = endian_double_swap(cell_center[index].m_state.m_P);
			    val[2] = endian_double_swap(cell_center[index].m_state.m_mu);
			    fwrite(val, sizeof(double), 3, outfile);
			    val[0] = endian_double_swap(cell_center[index].m_state.m_U[0]);
			    val[1] = endian_double_swap(cell_center[index].m_state.m_U[1]);
			    fwrite(val, sizeof(double), 2, outfile);
			    ival[0] = endian_int_swap(top_comp[index]);
			    fwrite(ival, sizeof(int), 1, outfile);
			}
		    }
		    else //big endian system
		    {
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			{
			    index = d_index2d(i,j,top_gmax);
			    val[0] = cell_center[index].m_state.m_rho;
			    val[1] = cell_center[index].m_state.m_P;
			    val[2] = cell_center[index].m_state.m_mu;
			    fwrite(val, sizeof(double), 3, outfile);
			    val[0] = cell_center[index].m_state.m_U[0];
			    val[1] = cell_center[index].m_state.m_U[1];
			    fwrite(val, sizeof(double), 2, outfile);
			    ival[0] = top_comp[index];
			    fwrite(ival, sizeof(int), 1, outfile);
			}
		    }
		    break;
		    case 3:
		    if(hardware_is_little_endian())
		    {
                        for (i = rep_top_gmin[0]; i <= rep_top_gmax[0]; ++i)
                        for (j = rep_top_gmin[1]; j <= rep_top_gmax[1]; ++j)
                        for (k = rep_top_gmin[2]; k <= rep_top_gmax[2]; ++k)
                        {
                            if (RegridRun == YES)
                            { 
                            index_i = rep_top_gmin[0] + lbuf[0];
                            index_j = rep_top_gmin[1] + lbuf[1];
                            index_k = rep_top_gmin[2] + lbuf[2];
                            index = d_index3d(index_i,index_j,index_k,top_gmax);
                            coords[0] = cell_center[index].m_coords[0]-(0.5*top_h[0])-((lbuf[0]-0.5)*top_h[0]/ratio_regrid[0])+(i*top_h[0]/ratio_regrid[0]);
                            coords[1] = cell_center[index].m_coords[1]-(0.5*top_h[1])-((lbuf[1]-0.5)*top_h[1]/ratio_regrid[1])+(j*top_h[1]/ratio_regrid[1]);
                            coords[2] = cell_center[index].m_coords[2]-(0.5*top_h[2])-((lbuf[2]-0.5)*top_h[2]/ratio_regrid[2])+(k*top_h[2]/ratio_regrid[2]);

                            FT_IntrpStateVarAtCoords(front,NO_COMP,coords,vel[0],getStateXvel,&vel_out[0]);
                            FT_IntrpStateVarAtCoords(front,NO_COMP,coords,vel[1],getStateYvel,&vel_out[1]);
                            FT_IntrpStateVarAtCoords(front,NO_COMP,coords,vel[2],getStateZvel,&vel_out[2]);
                            FT_IntrpStateVarAtCoords(front,NO_COMP,coords,pres,getStatePres,&pres_out);
                            rho_out = 1.0;
                            mu_out = 1.0;
                            temp = 1;
                            val[0] = endian_double_swap(rho_out);
                            val[1] = endian_double_swap(pres_out);
                            val[2] = endian_double_swap(mu_out);
                            fwrite(val, sizeof(double), 3, outfile);
                            val[0] = endian_double_swap(vel_out[0]);
                            val[1] = endian_double_swap(vel_out[1]);
                            val[2] = endian_double_swap(vel_out[2]);
                            fwrite(val, sizeof(double), 3, outfile);
                            ival[0] = endian_int_swap(temp);
                            fwrite(ival, sizeof(int), 1, outfile);
                            }
                            else
                            {
                            index = d_index3d(i,j,k,top_gmax);
                            val[0] = endian_double_swap(cell_center[index].m_state.m_rho);
                            val[1] = endian_double_swap(cell_center[index].m_state.m_P);
                            val[2] = endian_double_swap(cell_center[index].m_state.m_mu);
                            fwrite(val, sizeof(double), 3, outfile);
                            val[0] = endian_double_swap(cell_center[index].m_state.m_U[0]);
                            val[1] = endian_double_swap(cell_center[index].m_state.m_U[1]);
                            val[2] = endian_double_swap(cell_center[index].m_state.m_U[2]);
                            fwrite(val, sizeof(double), 3, outfile);
                            ival[0] = endian_int_swap(top_comp[index]);
                            fwrite(ival, sizeof(int), 1, outfile);
                            }
                        }
		    }
		    else //big endian system
		    {
                        for (i = rep_top_gmin[0]; i <= rep_top_gmax[0]; ++i)
                        for (j = rep_top_gmin[1]; j <= rep_top_gmax[1]; ++j)
                        for (k = rep_top_gmin[2]; k <= rep_top_gmax[2]; ++k)
                        {
                            if (RegridRun == YES)
                            {
                            index_i = rep_top_gmin[0] + lbuf[0];
                            index_j = rep_top_gmin[1] + lbuf[1];
                            index_k = rep_top_gmin[2] + lbuf[2];
                            index = d_index3d(index_i,index_j,index_k,top_gmax);
                            coords[0] = cell_center[index].m_coords[0]-(0.5*top_h[0])-((lbuf[0]-0.5)*top_h[0]/ratio_regrid[0])+(i*top_h[0]/ratio_regrid[0]);
                            coords[1] = cell_center[index].m_coords[1]-(0.5*top_h[1])-((lbuf[1]-0.5)*top_h[1]/ratio_regrid[1])+(j*top_h[1]/ratio_regrid[1]);
                            coords[2] = cell_center[index].m_coords[2]-(0.5*top_h[2])-((lbuf[2]-0.5)*top_h[2]/ratio_regrid[2])+(k*top_h[2]/ratio_regrid[2]);

                            FT_IntrpStateVarAtCoords(front,NO_COMP,coords,vel[0],getStateXvel,&vel_out[0]);
                            FT_IntrpStateVarAtCoords(front,NO_COMP,coords,vel[1],getStateYvel,&vel_out[1]);
                            FT_IntrpStateVarAtCoords(front,NO_COMP,coords,vel[2],getStateZvel,&vel_out[2]);
                            FT_IntrpStateVarAtCoords(front,NO_COMP,coords,pres,getStatePres,&pres_out);
                            rho_out = 1.0;
                            mu_out = 1.0;
                            temp = 1;

                            val[0] = rho_out;
                            val[1] = pres_out;
                            val[2] = mu_out;
                            fwrite(val, sizeof(double), 3, outfile);
                            val[0] = vel_out[0];
                            val[1] = vel_out[1];
                            val[2] = vel_out[2];
                            fwrite(val, sizeof(double), 3, outfile);
                            ival[0] = temp;
                            fwrite(ival, sizeof(int), 1, outfile);
                            }
                            else
                            {
                            index = d_index3d(i,j,k,top_gmax);
                            val[0] = cell_center[index].m_state.m_rho;
                            val[1] = cell_center[index].m_state.m_P;
                            val[2] = cell_center[index].m_state.m_mu;
                            fwrite(val, sizeof(double), 3, outfile);
                            val[0] = cell_center[index].m_state.m_U[0];
                            val[1] = cell_center[index].m_state.m_U[1];
                            val[2] = cell_center[index].m_state.m_U[2];
                            fwrite(val, sizeof(double), 3, outfile);
                            ival[0] = top_comp[index];
                            fwrite(ival, sizeof(int), 1, outfile);
                            }
                        }
		    }
		    break;
	    }
	    fclose(outfile);
	    }
	    else//ASCII output
	    {
	    outfile = fopen(filename,"w");

	    //fluid_print_front_states(outfile,front,NO);

	    fprintf(outfile,"\nInterior ifluid states:\n");
	    switch (dim)
	    {
		case 2:
		    for (i = 0; i <= top_gmax[0]; ++i)
		    for (j = 0; j <= top_gmax[1]; ++j)
		    {
			index = d_index2d(i,j,top_gmax);
			fprintf(outfile,"%24.18g\n",
				cell_center[index].m_state.m_rho);
			fprintf(outfile,"%24.18g\n",
				cell_center[index].m_state.m_P);
			fprintf(outfile,"%24.18g\n",
				cell_center[index].m_state.m_mu);
			for (l = 0; l < dim; ++l)
			    fprintf(outfile,"%24.18g\n",
				    cell_center[index].m_state.m_U[l]);
			fprintf(outfile, "%d\n", top_comp[index]);
		    }
		    break;
		case 3:
                    for (i = rep_top_gmin[0]; i <= rep_top_gmax[0]; ++i)
                    for (j = rep_top_gmin[1]; j <= rep_top_gmax[1]; ++j)
                    for (k = rep_top_gmin[2]; k <= rep_top_gmax[2]; ++k)

                    {
                        if (RegridRun == YES)
                        {
                        index_i = rep_top_gmin[0] + lbuf[0];
                        index_j = rep_top_gmin[1] + lbuf[1];
                        index_k = rep_top_gmin[2] + lbuf[2];
                        index = d_index3d(index_i,index_j,index_k,top_gmax);
                        coords[0] = cell_center[index].m_coords[0]-(0.5*top_h[0])-((lbuf[0]-0.5)*top_h[0]/ratio_regrid[0])+(i*top_h[0]/ratio_regrid[0]);
                        coords[1] = cell_center[index].m_coords[1]-(0.5*top_h[1])-((lbuf[1]-0.5)*top_h[1]/ratio_regrid[1])+(j*top_h[1]/ratio_regrid[1]);
                        coords[2] = cell_center[index].m_coords[2]-(0.5*top_h[2])-((lbuf[2]-0.5)*top_h[2]/ratio_regrid[2])+(k*top_h[2]/ratio_regrid[2]);

                        FT_IntrpStateVarAtCoords(front,NO_COMP,coords,vel[0],getStateXvel,&vel_out[0]);
                        FT_IntrpStateVarAtCoords(front,NO_COMP,coords,vel[1],getStateYvel,&vel_out[1]);
                        FT_IntrpStateVarAtCoords(front,NO_COMP,coords,vel[2],getStateZvel,&vel_out[2]);
                        FT_IntrpStateVarAtCoords(front,NO_COMP,coords,pres,getStatePres,&pres_out);
                        rho_out = 1.0;
                        mu_out = 1.0;
                        temp = 1;

                        fprintf(outfile,"%24.18g\n",
                                rho_out);
                        fprintf(outfile,"%24.18g\n",
                                pres_out);
                        fprintf(outfile,"%24.18g\n",
                                mu_out);
                        for (l = 0; l < dim; ++l)
                            fprintf(outfile,"%24.18g\n",
                                    vel_out[l]);
                        fprintf(outfile, "%d\n", temp);
                        }
                        else
                        {
                        index = d_index3d(i,j,k,top_gmax);
                        fprintf(outfile,"%24.18g\n",
                                cell_center[index].m_state.m_rho);
                        fprintf(outfile,"%24.18g\n",
                                cell_center[index].m_state.m_P);
                        fprintf(outfile,"%24.18g\n",
                                cell_center[index].m_state.m_mu);
                        for (l = 0; l < dim; ++l)
                            fprintf(outfile,"%24.18g\n",
                                    cell_center[index].m_state.m_U[l]);
                        fprintf(outfile, "%d\n", top_comp[index]);
                        }
                    }
            }
            fclose(outfile);
            free_front(tempfront);
	    }  //else 
        }  //kk
        }  //ii
}


void Incompress_Solver_Smooth_Basis::readFrontInteriorStates(char *restart_name, bool binary, bool RegridRestart)
{
    	printf("\nEnter readFrontInteriorStates \n");
	FILE *infile;
	int i,j,k,l,index,temp;
	char fname[100];

	double val[3];
	int ival[3];

	m_rho[0] = iFparams->rho1;		
	m_rho[1] = iFparams->rho2;		
	m_mu[0] = iFparams->mu1;		
	m_mu[1] = iFparams->mu2;		
	m_sigma = iFparams->surf_tension;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	mu_min = HUGE;
	for(i = 0; i < 2; i++)
	{
	    if(ifluid_comp(m_comp[i]))
	    {
		mu_min = std::min(mu_min,m_mu[i]);
		rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

	if (binary == YES)
	    infile = fopen(restart_name,"rb");
	else
	    infile = fopen(restart_name,"r");
	
	/* Initialize states at the interface */
        if (RegridRestart != YES)
	    fluid_read_front_states(infile,front,binary);

	FT_MakeGridIntfc(front);
	setDomain();

	next_output_line_containing_string(infile,"Interior ifluid states:");

	if (binary == YES)
	{
	    if (hardware_is_little_endian())
	    {
		switch (dim)
		{
		    case 2:
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			{
			    index = d_index2d(i,j,top_gmax);
			    fread(val, sizeof(double), 3, infile);
			    cell_center[index].m_state.m_rho = endian_double_swap(val[0]);
			    cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
			    cell_center[index].m_state.m_P = endian_double_swap(val[1]);
			    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
			    cell_center[index].m_state.m_mu = endian_double_swap(val[2]);
			    cell_center[index].m_state.m_mu_old = cell_center[index].m_state.m_mu;

			    fread(val, sizeof(double), 2, infile);
			    cell_center[index].m_state.m_U[0] = endian_double_swap(val[0]);
			    cell_center[index].m_state.m_U[1] = endian_double_swap(val[1]);

			    fread(ival, sizeof(int), 1, infile);
                            if (RegridRestart != YES)
			        top_comp[index] = endian_int_swap(ival[0]);
                            else
                                temp = endian_int_swap(ival[0]);

			}
			break;
		    case 3:
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			for (k = 0; k <= top_gmax[2]; ++k)
			{
			    index = d_index3d(i,j,k,top_gmax);
			    fread(val, sizeof(double), 3, infile);
			    cell_center[index].m_state.m_rho = endian_double_swap(val[0]);
			    cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
			    cell_center[index].m_state.m_P = endian_double_swap(val[1]);
			    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
			    cell_center[index].m_state.m_mu = endian_double_swap(val[2]);
			    cell_center[index].m_state.m_mu_old = cell_center[index].m_state.m_mu;

			    fread(val, sizeof(double), 3, infile);
			    cell_center[index].m_state.m_U[0] = endian_double_swap(val[0]);
			    cell_center[index].m_state.m_U[1] = endian_double_swap(val[1]);
			    cell_center[index].m_state.m_U[2] = endian_double_swap(val[2]);

			    fread(ival, sizeof(int), 1, infile);
                            if (RegridRestart != YES)
                                top_comp[index] = endian_int_swap(ival[0]);
                            else
                                temp = endian_int_swap(ival[0]);
			}
			break;
		}

	    }
	    else //big endian system
	    {
		switch (dim)
		{
		    case 2:
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			{
			    index = d_index2d(i,j,top_gmax);
			    fread(val, sizeof(double), 3, infile);
			    cell_center[index].m_state.m_rho = val[0];
			    cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
			    cell_center[index].m_state.m_P = val[1];
			    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
			    cell_center[index].m_state.m_mu = val[2];
			    cell_center[index].m_state.m_mu_old = cell_center[index].m_state.m_mu;

			    fread(val, sizeof(double), 2, infile);
			    cell_center[index].m_state.m_U[0] = val[0];
			    cell_center[index].m_state.m_U[1] = val[1];

			    fread(ival, sizeof(int), 1, infile);
                            if (RegridRestart != YES)
			        top_comp[index] = ival[0];
                            else
                                temp = ival[0];
			}
			break;
		    case 3:
			for (i = 0; i <= top_gmax[0]; ++i)
			for (j = 0; j <= top_gmax[1]; ++j)
			for (k = 0; k <= top_gmax[2]; ++k)
			{
			    index = d_index3d(i,j,k,top_gmax);
			    fread(val, sizeof(double), 3, infile);
			    cell_center[index].m_state.m_rho = val[0];
			    cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
			    cell_center[index].m_state.m_P = val[1];
			    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
			    cell_center[index].m_state.m_mu = val[2];
			    cell_center[index].m_state.m_mu_old = cell_center[index].m_state.m_mu;

			    fread(val, sizeof(double), 3, infile);
			    cell_center[index].m_state.m_U[0] = val[0];
			    cell_center[index].m_state.m_U[1] = val[1];
			    cell_center[index].m_state.m_U[2] = val[2];

			    fread(ival, sizeof(int), 1, infile);
                            if (RegridRestart != YES)
			        top_comp[index] = ival[0];
                            else
                                temp = ival[0];
			}
			break;
		}
	    }
	}
	else //reading ASCII file
	{
	    switch (dim)
	    {
		case 2:
		    for (i = 0; i <= top_gmax[0]; ++i)
		    for (j = 0; j <= top_gmax[1]; ++j)
		    {
			index = d_index2d(i,j,top_gmax);
			fscanf(infile,"%lf",&cell_center[index].m_state.m_rho);
			cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
			fscanf(infile,"%lf",&cell_center[index].m_state.m_P);
			cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
			fscanf(infile,"%lf",&cell_center[index].m_state.m_mu);
			cell_center[index].m_state.m_mu_old = cell_center[index].m_state.m_mu;
			for (l = 0; l < dim; ++l)
			    fscanf(infile,"%lf",&cell_center[index].m_state.m_U[l]);
                        if (RegridRestart != YES)
			    fscanf(infile, "%d", &top_comp[index]);
                        else
                            fscanf(infile, "%d", &temp);
		    }
		    break;
		case 3:
		    for (i = 0; i <= top_gmax[0]; ++i)
		    for (j = 0; j <= top_gmax[1]; ++j)
		    for (k = 0; k <= top_gmax[2]; ++k)
		    {
			index = d_index3d(i,j,k,top_gmax);
			fscanf(infile,"%lf",&cell_center[index].m_state.m_rho);
			cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
			fscanf(infile,"%lf",&cell_center[index].m_state.m_P);
			cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
			fscanf(infile,"%lf",&cell_center[index].m_state.m_mu);
			cell_center[index].m_state.m_mu_old = cell_center[index].m_state.m_mu;
			for (l = 0; l < dim; ++l)
			    fscanf(infile,"%lf",&cell_center[index].m_state.m_U[l]);
                        if (RegridRestart != YES)
			    fscanf(infile, "%d", &top_comp[index]);
                        else
                            fscanf(infile, "%d", &temp);
		    }
		    break;
	    }
	}
	fclose(infile);

	computeGradientQ();
	copyMeshStates();
	setAdvectionDt();
}


void Incompress_Solver_Smooth_Basis::setAdvectionDt()
{
	pp_global_max(&max_speed,1);
	if (max_speed != 0.0)
	    max_dt = hmin/max_speed;
	else
	    max_dt = HUGE;
	min_dt = 0.0000001*sqr(hmin)/mu_min;
	if (debugging("trace"))
	{
	    if (max_dt == HUGE)
	    	(void) printf("In setAdvectionDt: \n"
			"max_dt = HUGE min_dt = %24.18g\n",min_dt);
	    else
	    	(void) printf("In setAdvectionDt:\n"
			"max_dt = %24.18g min_dt = %24.18g\n",max_dt,min_dt);
	}
}	/* end setAdvectionDt */


void Incompress_Solver_Smooth_Basis::augmentMovieVariables()
{
	int i,j,k,index;
	static HDF_MOVIE_VAR *hdf_movie_var;
	static VTK_MOVIE_VAR *vtk_movie_var;
	static double *pres,*vort,*xvel,*yvel,*zvel;
	int n,offset,num_var;
	IF_MOVIE_OPTION *movie_option = iFparams->movie_option;

	offset = front->hdf_movie_var->num_var;
	num_var = (dim == 2) ? offset + 4 : offset + 12;

	if (hdf_movie_var == NULL)
	{
	    /* Added for vtk movie of vector field */
	    if (movie_option->plot_vel_vector)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&vtk_movie_var,
				sizeof(VTK_MOVIE_VAR));
		vtk_movie_var->num_var = 1;
		FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->top_var,1,
					sizeof(double**));
		FT_MatrixMemoryAlloc((POINTER*)&vtk_movie_var->var_name,1,100,
					sizeof(char));
	    	sprintf(vtk_movie_var->var_name[0],"velo");
		vtk_movie_var->top_var[0] = field->vel;
		front->vtk_movie_var = vtk_movie_var;
	    }
	    else
		vtk_movie_var = NULL;

	    /* Begin hdf movie */
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,num_var,
					100,sizeof(char));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,num_var,
					sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->idir,num_var,
					sizeof(int));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,
					num_var,sizeof(COMPONENT));
	    for (i = 0; i < front->hdf_movie_var->num_var; ++i)
	    {
		strcpy(hdf_movie_var->var_name[i],
			front->hdf_movie_var->var_name[i]);
		hdf_movie_var->get_state_var[i] =
			front->hdf_movie_var->get_state_var[i];
		hdf_movie_var->top_var[i] = 
			front->hdf_movie_var->top_var[i];
		hdf_movie_var->obstacle_comp[i] = 
			front->hdf_movie_var->obstacle_comp[i];
	    }
	    
	    switch (dim)
	    {
	    case 2:
	        hdf_movie_var->plot_comp = YES;
		hdf_movie_var->num_var = n = offset;
		if (movie_option->plot_pres)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"pres");
	    	    hdf_movie_var->get_state_var[n] = getStatePres;
	    	    hdf_movie_var->top_var[n] = field->pres;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_vort)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"vort");
	    	    hdf_movie_var->get_state_var[n] = getStateVort;
	    	    hdf_movie_var->top_var[n] = field->vort;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_velo)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"xvel");
	    	    hdf_movie_var->get_state_var[n] = getStateXvel;
	    	    hdf_movie_var->top_var[n] = field->vel[0];
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
	    	    sprintf(hdf_movie_var->var_name[n],"yvel");
	    	    hdf_movie_var->get_state_var[n] = getStateYvel;
	    	    hdf_movie_var->top_var[n] = field->vel[1];
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		break;
	    case 3:
		if (!movie_option->plot_comp)
	            hdf_movie_var->plot_comp = NO;
		hdf_movie_var->num_var = n = offset;
		if (movie_option->plot_cross_section[0])
		{
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = field->pres;
	    		hdf_movie_var->idir[n] = 0;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-yz-y");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvel;
	    		hdf_movie_var->top_var[n] = field->vel[1];
	    		hdf_movie_var->idir[n] = 0;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-yz-z");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvel;
	    		hdf_movie_var->top_var[n] = field->vel[2];
	    		hdf_movie_var->idir[n] = 0;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_vort)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"vort-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvort;
	    		hdf_movie_var->top_var[n] = field->vort3d[0];
	    		hdf_movie_var->idir[n] = 0;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
		if (movie_option->plot_cross_section[1])
		{
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-xz");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = field->pres;
	    		hdf_movie_var->idir[n] = 1;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xz-x");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvel;
	    		hdf_movie_var->top_var[n] = field->vel[0];
	    		hdf_movie_var->idir[n] = 1;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xz-z");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvel;
	    		hdf_movie_var->top_var[n] = field->vel[2];
	    		hdf_movie_var->idir[n] = 1;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_vort)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"vort-xz");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvort;
	    		hdf_movie_var->top_var[n] = field->vort3d[1];
	    		hdf_movie_var->idir[n] = 1;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
		if (movie_option->plot_cross_section[2])
		{
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = field->pres;
	    		hdf_movie_var->idir[n] = 2;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xy-x");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvel;
	    		hdf_movie_var->top_var[n] = field->vel[0];
	    		hdf_movie_var->idir[n] = 2;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xy-y");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvel;
	    		hdf_movie_var->top_var[n] = field->vel[1];
	    		hdf_movie_var->idir[n] = 2;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_vort)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"vort-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvort;
	    		hdf_movie_var->top_var[n] = field->vort3d[2];
	    		hdf_movie_var->idir[n] = 2;
	    	        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
	    }
	    FT_FreeThese(3,front->hdf_movie_var->var_name,
			front->hdf_movie_var->top_var,
			front->hdf_movie_var->obstacle_comp);
	    FT_FreeThese(1,front->hdf_movie_var);
	    front->hdf_movie_var = hdf_movie_var;
	    front->vtk_movie_var = vtk_movie_var;
	}

}	/* end augmentMovieVariables */

void Incompress_Solver_Smooth_Basis::initMovieVariables()
{
	int i,j,k,index;
	static HDF_MOVIE_VAR *hdf_movie_var;
	static VTK_MOVIE_VAR *vtk_movie_var;
	static double *pres,*vort,*xvel,*yvel,*zvel;
	int n;
	IF_MOVIE_OPTION *movie_option = iFparams->movie_option;

	if (hdf_movie_var == NULL)
	{
	    /* Added for vtk movie of vector field */
	    if (movie_option->plot_vel_vector)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&vtk_movie_var,
				sizeof(VTK_MOVIE_VAR));
		vtk_movie_var->num_var = 1;
		FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->top_var,1,
					sizeof(double**));
		FT_MatrixMemoryAlloc((POINTER*)&vtk_movie_var->var_name,1,100,
					sizeof(char));
	    	sprintf(vtk_movie_var->var_name[0],"velo");
		vtk_movie_var->top_var[0] = field->vel;
		front->vtk_movie_var = vtk_movie_var;
	    }
	    else
		front->vtk_movie_var = NULL;

	    /* Begin hdf movies */
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    switch (dim)
	    {
	    case 2:
	    	hdf_movie_var->plot_comp = YES;
		hdf_movie_var->num_var = n = 0;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,4,100,
					sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,4,
					sizeof(double*));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,4,
					sizeof(COMPONENT));
		if (movie_option->plot_pres)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"pres");
	    	    hdf_movie_var->get_state_var[n] = getStatePres;
	    	    hdf_movie_var->top_var[n] = field->pres;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_vort)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"vort");
	    	    hdf_movie_var->get_state_var[n] = getStateVort;
	    	    hdf_movie_var->top_var[n] = field->vort;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_velo)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"xvel");
	    	    hdf_movie_var->get_state_var[n] = getStateXvel;
	    	    hdf_movie_var->top_var[n] = field->vel[0];
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
	    	    sprintf(hdf_movie_var->var_name[n],"yvel");
	    	    hdf_movie_var->get_state_var[n] = getStateYvel;
	    	    hdf_movie_var->top_var[n] = field->vel[1];
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		break;
	    case 3:
		hdf_movie_var->num_var = n = 0;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,
					12,100,sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,
					12,sizeof(double*));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->idir,
					12,sizeof(int));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,
					12,sizeof(COMPONENT));
		if (!movie_option->plot_comp)
		    hdf_movie_var->plot_comp = NO;
		if (movie_option->plot_cross_section[0])
		{
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = field->pres;
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-yz-y");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvel;
	    		hdf_movie_var->top_var[n] = field->vel[1];
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-yz-z");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvel;
	    		hdf_movie_var->top_var[n] = field->vel[2];
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_vort)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"vort-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvort;
	    		hdf_movie_var->top_var[n] = field->vort3d[0];
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
		if (movie_option->plot_cross_section[1])
		{
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-xz");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = field->pres;
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xz-x");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvel;
	    		hdf_movie_var->top_var[n] = field->vel[0];
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xz-z");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvel;
	    		hdf_movie_var->top_var[n] = field->vel[2];
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_vort)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"vort-xz");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvort;
	    		hdf_movie_var->top_var[n] = field->vort3d[1];
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
		if (movie_option->plot_cross_section[2])
		{
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = field->pres;
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xy-x");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvel;
	    		hdf_movie_var->top_var[n] = field->vel[0];
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xy-y");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvel;
	    		hdf_movie_var->top_var[n] = field->vel[1];
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_vort)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"vort-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvort;
	    		hdf_movie_var->top_var[n] = field->vort3d[2];
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
	    }
	    front->hdf_movie_var = hdf_movie_var;
	}
}	/* end initMovieVariables */



//-------------------------------------------------------------------------------
//               Incompress_Solver_Smooth_2D_Basis
//------------------------------------------------------------------------------

double Incompress_Solver_Smooth_2D_Basis::getSmoothingFunction(double phi)
{
	// Heaviside function [1]
	if (phi < -m_smoothing_radius)	
	    return 0;
	else if (phi > m_smoothing_radius)
	    return 1;
	else
	    return 1.0/2 + phi/(2*m_smoothing_radius) + 
		   1/(2*PI)*sin(PI*phi/m_smoothing_radius);
}

double Incompress_Solver_Smooth_2D_Basis::getSmoothingFunctionD(double *center, double *point)
{
        if (fabs(center[0]-point[0]) < 2*top_h[0] && 
		fabs(center[1]-point[1]) < 2*top_h[1])
	    return ((1.0 + cos((PI*(center[0]-point[0]))/(2.0*top_h[0]))) 
		    *(1.0 + cos((PI*(center[1]-point[1]))/(2.0*top_h[1])))) 
		    /(16.0*top_h[0]*top_h[1]);
	else
            return 0.0;
}	/* end getSmoothingFunctionD in 2D */

double Incompress_Solver_Smooth_2D_Basis::smoothedDeltaFunction(double *p, double *center)
{
	int i;
	double len,d[MAXD],r,delta;

	for (i = 0; i < dim; ++i) d[i] = p[i] - center[i];
	len = mag_vector(d,dim);
	r = len/hmin/m_smoothing_radius;	
	delta = (fabs(r) > 1.0) ? 0.0 : 0.5*(1.0 + cos(PI*r));
	return delta;
}	/* end smoothedDeltaFunction */

double Incompress_Solver_Smooth_2D_Basis::smoothedStepFunction(double *p, double *center, int sign)
{
	int i;
	double dist,dp[MAXD],x,H;

	for (i = 0; i < dim; ++i) dp[i] = p[i] - center[i];
	dist = mag_vector(dp,dim);
	x = sign*dist/hmin/m_smoothing_radius;	
	if (x < -1.0) H = 0.0;
	else if (x > 1.0) H = 1.0;
	else H = 0.5*(1.0 + x) + 0.5/PI*sin(PI*x);
	return H;
}	/* end smoothedStepFunction */

void Incompress_Solver_Smooth_2D_Basis::sampleVelocity()
{
        int i,j,index;
	SAMPLE *sample = front->sample;
	char *sample_type = sample->sample_type;
	double *line = sample->sample_coords;
	char *out_name = front->out_name;
        double coords[MAXD];
        double velo1,velo2,velo;
        FILE *sfile;
        char sname[100];
        static int count = 0;
        static int step = 0;
        static int l = -1;
        static double lambda;

	if (front->step < sample->start_step || front->step > sample->end_step)
	    return;
	if ((front->step - sample->start_step)%sample->step_interval)
	    return;
        if (step != front->step)
        {
            step = front->step;
            count = 0;
        }
        switch (sample_type[0])
        {
        case 'x':
            sprintf(sname, "%s/vertical-x-%d-%d.xg",out_name,step,count);
            sfile = fopen(sname,"w");
            if (l == -1)
            {
                double x1,x2;
                do
                {
                    ++l;
                    index = d_index2d(l,0,top_gmax);
                    getRectangleCenter(index, coords);
                } while(line[0] >= coords[0]);
                --l;
                index = d_index2d(l,0,top_gmax);
                getRectangleCenter(index,coords);
                x1 = coords[0];
                index = d_index2d(l+1,0,top_gmax);
                getRectangleCenter(index,coords);
                x2 = coords[0];
                lambda = (line[0] - x1) / (x2 - line[0]);
            }
            i = l;
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                velo1 = cell_center[index].m_state.m_U[0];
                index = d_index2d(i+1,j,top_gmax);
                velo2 = cell_center[index].m_state.m_U[0];
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
            }
            fclose(sfile);
            sprintf(sname,"%s/vertical-y-%d-%d.xg",out_name,step,count++);
            sfile = fopen(sname,"w");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                velo1 = cell_center[index].m_state.m_U[1];
                index = d_index2d(i+1,j,top_gmax);
                velo2 = cell_center[index].m_state.m_U[1];
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
            }
            fclose(sfile);
            break;
        case 'y':
            sprintf(sname, "%s/horizontal-x-%d-%d.xg",out_name,step,count);
            sfile = fopen(sname,"w");
            if (l == -1)
            {
                double y1,y2;
                do
                {
                    ++l;
                    index = d_index2d(0,l,top_gmax);
                    getRectangleCenter(index, coords);
                } while (line[0] >= coords[1]);
                --l;
                index = d_index2d(0,l,top_gmax);
                getRectangleCenter(index,coords);
                y1 = coords[1];
                index = d_index2d(0,l+1,top_gmax);
                getRectangleCenter(index,coords);
                y2 = coords[1];
               lambda = (line[0] - y1) / (y2 - line[0]);
            }
            j = l;
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                velo1 = cell_center[index].m_state.m_U[0];
                index = d_index2d(i,j+1,top_gmax);
                velo2 = cell_center[index].m_state.m_U[0];
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);
            sprintf(sname,"%s/horizontal-y-%d-%d.xg",out_name,step,count++);
            sfile = fopen(sname,"w");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                velo1 = cell_center[index].m_state.m_U[1];
                index = d_index2d(i,j+1,top_gmax);
                velo2 = cell_center[index].m_state.m_U[1];
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);
            break;
        }
}	/* end sampleVelocity2d */

void Incompress_Solver_Smooth_2D_Basis::setSmoProOnePhase(void)
{
	boolean status;
	int i,j,l,index,sign;
	COMPONENT comp;
	double t[MAXD],*force;

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);			
	    comp  = cell_center[index].comp;
	    if (!ifluid_comp(comp)) continue;

	    force = cell_center[index].m_state.f_surf;
	    for (l = 0; l < dim; ++l) 
		force[l] = 0.0; //set the surface tension to 0

	
	    cell_center[index].m_state.m_mu = m_mu[0]; //set viscosity to any of the two (should be the same in input file)
	    cell_center[index].m_state.m_rho = m_rho[0]; //set density to any of the two (should be the same in input file)
	}

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = cell_center[index].m_state.m_mu;
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_mu = array[index];
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = cell_center[index].m_state.m_rho;
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_rho = array[index];
	}
	for (l = 0; l < dim; ++l)
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = cell_center[index].m_state.f_surf[l];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	cell_center[index].m_state.f_surf[l] = array[index];
	    }
	}
}

void Incompress_Solver_Smooth_2D_Basis::setSmoothedProperties(void)
{
	boolean status;
	int i,j,l,index,sign;
	COMPONENT comp;
	double t[MAXD],*force;
	double center[MAXD],point[MAXD],nor[MAXD],phi,H,D;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	int range = (int)(m_smoothing_radius+1);
	double rho;

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);			
	    comp  = cell_center[index].comp;
	    rho = cell_center[index].m_state.m_rho;
	    if (!ifluid_comp(comp)) continue;

	    getRectangleCenter(index, center);
	    status = FT_FindNearestIntfcPointInRange(front,comp,center,point,
				t,&hse,&hs,range);

	    force = cell_center[index].m_state.f_surf;
	    for (l = 0; l < dim; ++l) force[l] = 0.0;

	    if (status == YES && 
		ifluid_comp(positive_component(hs)) &&
		ifluid_comp(negative_component(hs)))
	    {
		sign = (comp == m_comp[0]) ? -1 : 1;
		D = smoothedDeltaFunction(center,point);
		H = smoothedStepFunction(center,point,sign);
		cell_center[index].m_state.m_mu = m_mu[0]  + 
					(m_mu[1]-m_mu[0])*H;
		cell_center[index].m_state.m_rho = m_rho[0] + 
					(m_rho[1]-m_rho[0])*H; 
		
		if (m_sigma != 0.0 && D != 0.0)
		{
		    surfaceTension_Fedkiw(D,hse,hs,force,m_sigma,rho,t);
		}
	    }
	    else
	    {
		switch (comp)
		{
		case LIQUID_COMP1:
		    cell_center[index].m_state.m_mu = m_mu[0];
		    cell_center[index].m_state.m_rho = m_rho[0];
		    break;
		case LIQUID_COMP2:
		    cell_center[index].m_state.m_mu = m_mu[1];
		    cell_center[index].m_state.m_rho = m_rho[1];
		    break;
		}
	    }
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = cell_center[index].m_state.m_mu;
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_mu = array[index];
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = cell_center[index].m_state.m_rho;
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_rho = array[index];
	}
	for (l = 0; l < dim; ++l)
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = cell_center[index].m_state.f_surf[l];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	cell_center[index].m_state.f_surf[l] = array[index];
	    }
	}
}	/* end setSmoothedProperties2d */

//-------------------------------------------------------------------------------
//               Incompress_Solver_Smooth_3D_Basis
//------------------------------------------------------------------------------

double Incompress_Solver_Smooth_3D_Basis::getSmoothingFunction(double phi)
{
	// Heaviside function [1]
	if (phi < -m_smoothing_radius)	
	    return 0;
	else if (phi > m_smoothing_radius)
	    return 1;
	else
	    return 1.0/2 + phi/(2*m_smoothing_radius) + 
		   1/(2*PI)*sin(PI*phi/m_smoothing_radius);
}

double Incompress_Solver_Smooth_3D_Basis::getSmoothingFunctionD(double *center, double *point)
{
        if (fabs(center[0]-point[0]) < 2*top_h[0] && 
		fabs(center[1]-point[1]) < 2*top_h[1] && 
		fabs(center[2]-point[2]) < 2*top_h[2])
	    return ((1.0 + cos((PI*(center[0]-point[0]))/(2.0*top_h[0]))) 
		    *(1.0 + cos((PI*(center[1]-point[1]))/(2.0*top_h[1]))) 
		    *(1.0 + cos((PI*(center[2]-point[2]))/(2.0*top_h[2])))) 
		/(64.0*top_h[0]*top_h[1]*top_h[2]);
	else
            return 0.0;
}	/* end getSmoothingFunctionD in 3D */

double Incompress_Solver_Smooth_3D_Basis::smoothedDeltaFunction(double *p, double *center)
{
	int i;
	double len,d[MAXD],r,delta;

	for (i = 0; i < dim; ++i) d[i] = p[i] - center[i];
	len = mag_vector(d,dim);
	r = len/hmin/m_smoothing_radius;	
	delta = (fabs(r) > 1.0) ? 0.0 : (0.5/hmin/m_smoothing_radius)*(1.0 + cos(PI*r));
	return delta;
}	/* end smoothedDeltaFunction */

double Incompress_Solver_Smooth_3D_Basis::smoothedStepFunction(double *p, double *center, int sign)
{
	int i;
	double dist,dp[MAXD],x,H;

	for (i = 0; i < dim; ++i) dp[i] = p[i] - center[i];
	dist = mag_vector(dp,dim);
	x = sign*dist/hmin/m_smoothing_radius;	
	if (x < -1.0) H = 0.0;
	else if (x > 1.0) H = 1.0;
	else H = 0.5*(1.0 + x) + 0.5/PI*sin(PI*x);
	return H;
}	/* end smoothedStepFunction */

void Incompress_Solver_Smooth_3D_Basis::sampleVelocity()
{
        int i,j,k,index;
        double coords[MAXD];
        double velo1,velo2,velo_tmp1,velo_tmp2,velo;
        FILE *sfile;
        char sname[100];
        static int count = 0;
        static int step = 0;
        static int l=-1,m=-1;
        static double lambda1,lambda2;
	SAMPLE *sample = front->sample;
	char *sample_type = sample->sample_type;
	double *sample_line = sample->sample_coords;
	char *out_name = front-> out_name;

	if (front->step < sample->start_step || front->step > sample->end_step)
	    return;
	if ((front->step - sample->start_step)%sample->step_interval)
	    return;
        if (step != front->step)
        {
            step = front->step;
            count = 0;
        }
        switch (sample_type[0])
        {
        case 'x':
            if (l == -1)
            {
                double x1,x2;
                do
                {
                    ++l;
                    index = d_index3d(l,0,0,top_gmax);
                    getRectangleCenter(index, coords);
                }while(sample_line[0]>=coords[0]);
                --l;
                index = d_index3d(l,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x1 = coords[0];
                index = d_index3d(l+1,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x2 = coords[0];
                lambda1 = (sample_line[0] - x1) / (x2 - sample_line[0]);
            }

            switch (sample_type[1])
            {
                case 'y':
                    if (m == -1)
                    {
                        double y1,y2;
                        do
                        {
                            ++m;
                            index = d_index3d(0,m,0,top_gmax);
                            getRectangleCenter(index,coords);
                        }while(sample_line[1]>=coords[1]);
                        --m;
                        index = d_index3d(0,m,0,top_gmax);
                        getRectangleCenter(index,coords);
                        y1 = coords[1];
                        index = d_index3d(0,m+1,0,top_gmax);
                        getRectangleCenter(index,coords);
                        y2 = coords[1];
                        lambda2 = (sample_line[1] - y1)/(y2 - sample_line[1]);
                    }
                    i = l;
                    j = m;
                    sprintf(sname, "%s/x-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (k = kmin; k <= kmax; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[0];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[0];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[0];
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[0];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (k = kmin; k <= kmax; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[1];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[1];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[1];
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[1];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (k = kmin; k <= kmax; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[2];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[2];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[2];
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[2];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    printf("sample line: x = %20.14f, y = %20.14f\n",coords[0],
                        coords[1]);

                    break;

                case 'z':
                    if (m == -1)
                    {
                        double z1,z2;
                        do
                        {
                            ++m;
                            index = d_index3d(0,0,m,top_gmax);
                            getRectangleCenter(index,coords);
                        }while(sample_line[1]>=coords[2]);
                        --m;
                        index = d_index3d(0,0,m,top_gmax);
                        getRectangleCenter(index,coords);
                        z1 = coords[2];
                        index = d_index3d(0,0,m+1,top_gmax);
                        getRectangleCenter(index,coords);
                        z2 = coords[2];
                        lambda2 = (sample_line[1] - z1)/(z2 - sample_line[1]);
                    }
                    i = l;
                    k = m;
                    sprintf(sname, "%s/x-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (j = jmin; j <= jmax; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[0];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[0];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[0];
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[0];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (j = jmin; j <= jmax; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[1];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[1];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[1];
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[1];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (j = jmin; j <= jmax; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[2];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[2];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[2];
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[2];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    printf("sample line: x = %20.14f, z = %20.14f\n",coords[0],
                        coords[2]);

                    break;

                    default:
                        printf("Incorrect input for sample velocity!\n");
                        break;

            }
            break;

        case 'y':
            if (l == -1)
            {
                double y1,y2;
                do
                {
                    ++l;
                    index = d_index3d(0,l,0,top_gmax);
                    getRectangleCenter(index, coords);
                }while(sample_line[0]>=coords[1]);
                --l;
                index = d_index3d(0,l,0,top_gmax);
                getRectangleCenter(index,coords);
                y1 = coords[1];
                index = d_index3d(0,l+1,0,top_gmax);
                getRectangleCenter(index,coords);
                y2 = coords[1];
                lambda1 = (sample_line[0] - y1)/(y2 - sample_line[0]);
            }

            switch (sample_type[1])
            {
                case 'z':
                    if (m == -1)
                    {
                        double z1,z2;
                        do
                        {
                            ++m;
                            index = d_index3d(0,0,m,top_gmax);
                            getRectangleCenter(index,coords);
                        }while(sample_line[1]>=coords[2]);
                        --m;
                        index = d_index3d(0,0,m,top_gmax);
                        getRectangleCenter(index,coords);
                        z1 = coords[2];
                        index = d_index3d(0,0,m+1,top_gmax);
                        getRectangleCenter(index,coords);
                        z2 = coords[2];
                        lambda2 = (sample_line[1] - z1)/(z2 - sample_line[1]);
                    }
                    j = l;
                    k = m;
                    sprintf(sname, "%s/x-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[0];
                        index = d_index3d(i,j+1,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[0];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[0];
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[0];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[1];
                        index = d_index3d(i,j+1,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[1];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[1];
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[1];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[2];
                        index = d_index3d(i,j+1,k,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[2];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = cell_center[index].m_state.m_U[2];
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        velo2 = cell_center[index].m_state.m_U[2];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    printf("sample line: y = %20.14f, z = %20.14f\n",coords[1],
                        coords[2]);

                    break;

                default:
                    printf("Incorrect input for sample velocity!\n");
                    break;
            }
        default:
            printf("Incorrect input for sample velocity!\n");
            break;
        }
}	/* end sampleVelocity in 3D */

void Incompress_Solver_Smooth_Basis::initSampleVelocity(char *in_name)
{
        FILE *infile;
	static SAMPLE *sample;
	char *sample_type;
	double *sample_line;

	infile = fopen(in_name,"r");
	FT_ScalarMemoryAlloc((POINTER*)&sample,sizeof(SAMPLE));
	sample_type = sample->sample_type;
	sample_line = sample->sample_coords;

	if (dim == 2)
	{
            CursorAfterString(infile,"Enter the sample line type:");
            fscanf(infile,"%s",sample_type);
            (void) printf(" %s\n",sample_type);
            CursorAfterString(infile,"Enter the sample line coordinate:");
            fscanf(infile,"%lf",sample_line);
            (void) printf(" %f\n",sample_line[0]);
	}
	else if (dim == 3)
        {
            CursorAfterString(infile,"Enter the sample line type:");
            fscanf(infile,"%s",sample_type);
            (void) printf(" %s\n",sample_type);
            CursorAfterString(infile,"Enter the sample line coordinate:");
            fscanf(infile,"%lf %lf",sample_line,sample_line+1);
            (void) printf(" %lf %lf\n",sample_line[0],sample_line[1]);
        }
        CursorAfterString(infile,"Enter the start step for sample: ");
        fscanf(infile,"%d",&sample->start_step);
        (void) printf("%d\n",sample->start_step);
        CursorAfterString(infile,"Enter the end step for sample: ");
        fscanf(infile,"%d",&sample->end_step);
        (void) printf("%d\n",sample->end_step);
        CursorAfterString(infile,"Enter the step interval for sample: ");
        fscanf(infile,"%d",&sample->step_interval);
        (void) printf("%d\n",sample->step_interval);
	front->sample = sample;
        fclose(infile);
}	/* end initSampleVelocity */

void Incompress_Solver_Smooth_3D_Basis::setSmoProOnePhase(void)
{
	int i,j,k,l,index,sign; 
	COMPONENT comp;
        double *force;

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);			
	    comp  = cell_center[index].comp;
	    if (!ifluid_comp(comp)) continue;

	    force = cell_center[index].m_state.f_surf;
            for (l = 0; l < dim; ++l) 
		force[l] = 0.0;

                cell_center[index].m_state.m_mu = m_mu[0];
                cell_center[index].m_state.m_rho = m_rho[0];
	}

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    index  = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_mu;
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{	
	    index  = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_mu = array[index];
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_rho;
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_rho = array[index];
	}
	for (l = 0; l < dim; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.f_surf[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.f_surf[l] = array[index];
	    }
	}
}

void Incompress_Solver_Smooth_3D_Basis::setSmoProExact(void)
//setting exact discrete heaviside function
{

}


void Incompress_Solver_Smooth_3D_Basis::setSmoothedProperties(void)
{
    // A boundary condition capturing method for multiphase incomrpessible flow,
    // M.Kang, R.P. Fedkiw, X-D, Liu, JSC, 2000
    // formula (63) & (91)
	boolean status;
	int i,j,k,l,index,sign; 
	COMPONENT comp;
        double t[MAXD],*force;
	double center[MAXD],point[MAXD],nor[MAXD],phi,H,D;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF *hs;
	double rho;

	//Make sure the curvature and point normal are calculated
	INTERFACE* intfc = front->interf;
	double curv;
	HYPER_SURF *hs_curv;
	HYPER_SURF_ELEMENT *hse_curv;
	POINT *p;
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse_curv,&hs_curv))
	    GetFrontCurvature(p,hse_curv,hs_curv,&curv,front);
	//////////////////////////////////////////////////
	

	start_clock("SetRhoAndMu");
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);			
	    comp  = cell_center[index].comp;
	    D = 0.0;
	    H = 0.0;

	    if (!ifluid_comp(comp)) continue;

	    getRectangleCenter(index, center);
            status = FT_FindNearestIntfcPointInRange(front,comp,center,point,
				t,&hse,&hs,(int)m_smoothing_radius);

	    force = cell_center[index].m_state.f_surf;
            for (l = 0; l < dim; ++l) force[l] = 0.0;

	    if (status  == YES && 
		ifluid_comp(positive_component(hs)) &&
                ifluid_comp(negative_component(hs)))
	    {
		sign = (comp == m_comp[0]) ? -1 : 1;
                D = smoothedDeltaFunction(center,point);
                H = smoothedStepFunction(center,point,sign);
                cell_center[index].m_state.m_mu = m_mu[0]  +
                                        (m_mu[1]-m_mu[0])*H;
                cell_center[index].m_state.m_rho = m_rho[0] +
                                        (m_rho[1]-m_rho[0])*H;
	    }
	    else
	    {
		switch (comp)
		{
		case LIQUID_COMP1:
		    cell_center[index].m_state.m_mu = m_mu[0];
		    cell_center[index].m_state.m_rho = m_rho[0];
		    break;
		case LIQUID_COMP2:
		    cell_center[index].m_state.m_mu = m_mu[1];
		    cell_center[index].m_state.m_rho = m_rho[1];
		    break;
		}
	    }
	    rho = cell_center[index].m_state.m_rho;

	/*    
	    if (m_sigma != 0.0 && D != 0.0)
	    {
                surfaceTension_Fedkiw(D,hse,hs,force,m_sigma,rho,t);
            }
	  */  
	}
	stop_clock("SetRhoAndMu");

	start_clock("SurfaceTension");
	surfaceTension_Peskin();
	stop_clock("SurfaceTension");

	start_clock("ScatStates");

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    index  = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_mu;
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{	
	    index  = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_mu = array[index];
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_rho;
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_rho = array[index];
	}
	for (l = 0; l < dim; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.f_surf[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.f_surf[l] = array[index];
	    }
	}
	stop_clock("ScatStates");
}
/* end setSmoothedProperties in 3D */



// Flux of Riemann solution of Burgers equation u_t + uu_x = 0

double burger_flux(	
	double ul,
	double um,
	double ur)
{
	double u_Rl,u_Rr;
	if (ul < um)
	{
	    if (ul > 0.0) u_Rl = ul;
	    else if (um < 0.0) u_Rl = um;
	    else u_Rl = 0.0;
	}
	else
	{
	    if (ul + um > 0.0) u_Rl = ul;
	    else u_Rl = um;
	}

	if (um < ur)
	{
	    if (um > 0.0) u_Rr = um;
	    else if (ur < 0.0) u_Rr = ur;
	    else u_Rr = 0.0;
	}
	else
	{
	    if (um + ur > 0.0) u_Rr = um;
	    else u_Rr = ur;
	}
	return 0.5*(u_Rr*u_Rr - u_Rl*u_Rl);
}	/* end flux */

// Flux of Riemann solution of linear equation u_t + au_x = 0

double linear_flux(	
	double a,
	double ul,
	double um,
	double ur)
{
	double u;
	
	if (a > 0.0)
	    return a*(um - ul);
	else
	    return a*(ur - um);
}	/* end net_uwind_flux */

/*
void Incompress_Solver_Smooth_3D_Basis::getLengthScale(void)
{
// Finding the distance of a random point from the interface in 3D 

	float dphi,drad,dz,radii,phi,z,start_rad,start_phi,start_z,dr,next_rad,next_z;
	float alpha;
	float distance_nor_pos,distance_nor_neg,distance_nor;
	float distance_nor_pos_z,distance_nor_neg_z,dist_z;
	float distance_norpos_z,distance_norneg_z,distance_tanpos_z,distance_tanneg_z;
	float distance_tan_pos,distance_tan_neg,distance_tan;
	float distance_1,distance_2,distance_3,distance_4,distance,distance_min;
	float distance1_z,distance2_z,distance3_z,distance4_z;
	COMPONENT comp,next_comp;
	float coords[MAXD], next_coords[MAXD];
	int Nphi,Nrad,Nz;
	char file_name[512];
	FILE *CLS;

	Nphi = 500;
	dphi = PI/Nphi;
	start_phi = -0.5*PI;
	Nrad = 20;
	drad = 0.628/Nrad;
	start_rad = 2.538;
	Nz = 100;
	dz = 10/Nz;
	start_z = 0;
	dr = rgrid->h[0]/2;

	for (i = 0; i <= Nrad; i++)
	{
	    radii = start_rad+i*drad;
	    sprintf(file_name,"%s_CLS_t%g_r%g_%d",output_filename(grid->init),
		front->time,radii,myid);
	    CLS = fopen(file_name,"w");
		
	    for (j = 0; j < Nphi; j++)
	    {
	    	phi = start_phi+j*dphi;
	    	
	    	for (k = 0; k < Nz; k++)
	    	{
	    	   z = start_z+k*dz;
	    	   
	    	coords[0] = radii*cos(phi);
	    	coords[1] = radii*sin(phi);
	    	coords[2] = z;
	    	comp = component(coords,front->interf);
	    	
		   if( (coords[0] >= L[0]) && (coords[0] <= U[0])
		    && (coords[1] >= L[1]) && (coords[1] <= U[1])
		    && (coords[2] >= L[2]) && (coords[2] <= U[2]))
		    
		    {
		        // normal positive side 
		          next_rad = radii+dr;
		          next_coords[0] = next_rad*cos(phi);
		          next_coords[1] = next_rad*sin(phi);
		          next_coords[2] = z;
		          next_comp = component(next_coords,front->interf);

		          while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		          && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		          && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		          && (next_coords[2] <= GU[2]))
		          {
		             next_z = z+dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z+dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_pos_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		                  
		             next_z = z-dr;                     
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z-dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_neg_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		             dist_z = min(distance_pos_z,distance_neg_z);   
		             distance_norpos_z = min(distance_norpos_z,dist_z);
		             next_rad = next_rad+dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = z;
		             next_comp = component(next_coords,front->interf);
		          }
		          distance_nor_pos = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));

				// normal negative side 
				next_rad = radii-dr;
				next_coords[0] = next_rad*cos(phi);
				next_coords[1] = next_rad*sin(phi);
				next_coords[2] = z;
				next_comp = component(next_coords,front->interf);

				while((next_comp == comp) && (next_coords[0]>=GL[0]) 
				&& (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
				&& (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
				&& (next_coords[2] <= GU[2]))
				{
					
		             next_z = z+dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z+dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_pos_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		                                  
		             next_z = z-dr;                     
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z-dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_neg_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		             dist_z = min(distance_pos_z,distance_neg_z);   
		             distance_norneg_z = min(distance_norneg_z,dist_z);
		             next_rad = next_rad-dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = z;
		             next_comp = component(next_coords,front->interf);
		      }
		        distance_nor_neg = sqrt(sqr(next_coords[0]-coords[0])+
		                                sqr(next_coords[1]-coords[1])+
		                                sqr(next_coords[2]-coords[2]));
	
		distance_nor = min(min(distance_nor_pos,distance_nor_neg),min(distance_norpos_z,distance_norneg_z));
	
		// tangential positive side 
		next_rad = sqrt(sqr(radii)+sqr(dr));
		next_coords[0] = next_rad*cos(phi-atan2(dr,radii));
		next_coords[1] = next_rad*sin(phi-atan2(dr,radii));
		next_coords[2] = z;
		next_comp = component(next_coords,front->interf);
		k = 1;

		while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		   && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		   && (next_coords[1] <= GU[1]) && (next_coords[2] >=GL[2])
		   && (next_coords[2] <= GU[2]))
		{
		             next_z = z+dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z+dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_pos_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		                  
		             next_z = z-dr;                     
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z-dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_neg_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		             dist_z = min(distance_pos_z,distance_neg_z);   
		             distance_tanpos_z = min(distance_tanpos_z,dist_z);
		     k++;
		    next_rad = sqrt(sqr(radii)+sqr(k*dr));
		    next_coords[0] = next_rad*cos(phi-atan2(k*dr,radii));
		    next_coords[1] = next_rad*sin(phi-atan2(k*dr,radii));
		    next_coords[2] = z;
		    next_comp = component(next_coords,front->interf);
		}
		distance_tan_pos = sqrt(sqr(next_coords[0]-coords[0])+
		                        sqr(next_coords[1]-coords[1])+
		                        sqr(next_coords[2]-coords[2]));

		// tangential negative side 
		next_rad = sqrt(sqr(radii)+sqr(dr));
		next_coords[0] = next_rad*cos(phi+atan2(dr,radii));
		next_coords[1] = next_rad*sin(phi+atan2(dr,radii));
		next_coords[2] = z;
		next_comp = component(next_coords,front->interf);
		k = 1;

		while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		   && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		   && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		   && (next_coords[2] <= GU[2]))
		{
		    
		             next_z = z+dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z+dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_pos_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		                                  
		             next_z = z-dr;                     
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z-dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_neg_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		             dist_z = min(distance_pos_z,distance_neg_z);   
		             distance_tanneg_z = min(distance_norneg_z,dist_z);
		    k++;
		    next_rad = sqrt(sqr(radii)+sqr(k*dr));
		    next_coords[0] = next_rad*cos(phi+atan2(k*dr,radii));
		    next_coords[1] = next_rad*sin(phi+atan2(k*dr,radii));
		    next_coords[2] = z;
		    next_comp = component(next_coords,front->interf);
		}
		distance_tan_neg = sqrt(sqr(next_coords[0]-coords[0])+
		                        sqr(next_coords[1]-coords[1])+
		                        sqr(next_coords[2]-coords[2]));
	
		distance_tan = min(min(distance_tan_pos,distance_tan_neg),min(distance_tanpos_z,distance_tanneg_z));
		
		distance_min = min(distance_tan,distance_nor);

		// random direction 
		for (kk = 1; kk < 5; kk++)
		{
		    alpha = kk*(PI/10);
		
		// side1 
		next_rad = sqrt(sqr(radii+dr*cos(alpha))+sqr(dr*sin(alpha)));
		next_coords[0] = next_rad*cos(phi+atan2(dr*sin(alpha),radii+dr*cos(alpha)));
		next_coords[1] = next_rad*sin(phi+atan2(dr*sin(alpha),radii+dr*cos(alpha)));
		next_coords[2] = z;
		next_comp = component(next_coords,front->interf);
		k = 1;

		while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		   && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		   && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		   && (next_coords[2]))
		{
		             next_z = z+dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z+dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_pos_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		                  
		             next_z = z-dr;                     
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z-dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_neg_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		             dist_z = min(distance_pos_z,distance_neg_z);   
		             distance1_z = min(distance1_z,dist_z);
		    k++;
		    next_rad = sqrt(sqr(radii+k*dr*cos(alpha))+sqr(k*dr*sin(alpha)));
		    next_coords[0] = next_rad*cos(phi+atan2(k*dr*sin(alpha),radii+k*dr*cos(alpha)));
		    next_coords[1] = next_rad*sin(phi+atan2(k*dr*sin(alpha),radii+k*dr*cos(alpha)));
		    next_coords[2] = z;
		    next_comp = component(next_coords,front->interf);
		}
		distance_1 = sqrt(sqr(next_coords[0]-coords[0])+
		                  sqr(next_coords[1]-coords[1])+
		                  sqr(next_coords[2]-coords[2]));
		
		// side2 
		next_rad = sqrt(sqr(radii+dr*cos(alpha))+sqr(dr*sin(alpha)));
		next_coords[0] = next_rad*cos(phi-atan2(dr*sin(alpha),radii+dr*cos(alpha)));
		next_coords[1] = next_rad*sin(phi-atan2(dr*sin(alpha),radii+dr*cos(alpha)));
		next_coords[2] = z;
		next_comp = component(next_coords,front->interf);
		k = 1;

		while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		   && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		   && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		   && (next_coords[2] <= GU[2]))
		{
		             next_z = z+dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z+dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_pos_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		                  
		             next_z = z-dr;                     
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z-dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_neg_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		             dist_z = min(distance_pos_z,distance_neg_z);   
		             distance2_z = min(distance2_z,dist_z);
		      k++;
		    next_rad = sqrt(sqr(radii+k*dr*cos(alpha))+sqr(k*dr*sin(alpha)));
		    next_coords[0] = next_rad*cos(phi-atan2(k*dr*sin(alpha),radii+k*dr*cos(alpha)));
		    next_coords[1] = next_rad*sin(phi-atan2(k*dr*sin(alpha),radii+k*dr*cos(alpha)));
		    next_coords[2] = z;
		    next_comp = component(next_coords,front->interf);
		}
		distance_2 = sqrt(sqr(next_coords[0]-coords[0])+
		                  sqr(next_coords[1]-coords[1])+
		                  sqr(next_coords[2]-coords[2]));
		
		// side3 
		next_rad = sqrt(sqr(radii-dr*cos(alpha))+sqr(dr*sin(alpha)));
		next_coords[0] = next_rad*cos(phi+atan2(dr*sin(alpha),radii-dr*cos(alpha)));
		next_coords[1] = next_rad*sin(phi+atan2(dr*sin(alpha),radii-dr*cos(alpha)));
		next_coords[2] = z;
		next_comp = component(next_coords,front->interf);
		k = 1;

		while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		   && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		   && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		   && (next_coords[2] <= GU[2]))
		{
		             next_z = z+dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z+dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_pos_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		                  
		             next_z = z-dr;                     
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z-dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_neg_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		             dist_z = min(distance_pos_z,distance_neg_z);   
		             distance3_z = min(distance3_z,dist_z);
		     k++;
		    next_rad = sqrt(sqr(radii-k*dr*cos(alpha))+sqr(k*dr*sin(alpha)));
		    next_coords[0] = next_rad*cos(phi+atan2(k*dr*sin(alpha),radii-k*dr*cos(alpha)));
		    next_coords[1] = next_rad*sin(phi+atan2(k*dr*sin(alpha),radii-k*dr*cos(alpha)));
		    next_coords[2] = z;
		    next_comp = component(next_coords,front->interf);
		}
		distance_3 = sqrt(sqr(next_coords[0]-coords[0])+
		                  sqr(next_coords[1]-coords[1])+
		                  sqr(next_coords[2]-coords[2]));
		
		// side4 
		next_rad = sqrt(sqr(radii-dr*cos(alpha))+sqr(dr*sin(alpha)));
		next_coords[0] = next_rad*cos(phi-atan2(dr*sin(alpha),radii-dr*cos(alpha)));
		next_coords[1] = next_rad*sin(phi-atan2(dr*sin(alpha),radii-dr*cos(alpha)));
	    next_coords[2] = z;
		next_comp = component(next_coords,front->interf);
		k = 1;

		while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		   && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		   && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		   && (next_coords[2] <= GU[2]))
		{
	        	     next_z = z+dr;
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z+dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_pos_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		                  
		             next_z = z-dr;                     
		             next_coords[0] = next_rad*cos(phi);
		             next_coords[1] = next_rad*sin(phi);
		             next_coords[2] = next_z;
		             next_comp = component(next_coords,front->interf);
		             while((next_comp == comp) && (next_coords[0]>=GL[0]) 
		                  && (next_coords[0] <= GU[0]) && (next_coords[1] >= GL[1])
		                  && (next_coords[1] <= GU[1]) && (next_coords[2] >= GL[2])
		                  && (next_coords[2] <= GU[2]))
		            {
		              next_z = next_z-dr;
		              next_coords[0] = next_rad*cos(phi);
		              next_coords[1] = next_rad*sin(phi);
		              next_coords[2] = next_z;
		              next_comp = component(next_coords,front->interf);
		            }
		             distance_neg_z = sqrt(sqr(next_coords[0]-coords[0])+
		                                  sqr(next_coords[1]-coords[1])+
		                                  sqr(next_coords[2]-coords[2]));
		             dist_z = min(distance_pos_z,distance_neg_z);   
		             distance4_z = min(distance4_z,dist_z);
		    k++;
		    next_rad = sqrt(sqr(radii-k*dr*cos(alpha))+sqr(k*dr*sin(alpha)));
		    next_coords[0] = next_rad*cos(phi-atan2(k*dr*sin(alpha),radii-k*dr*cos(alpha)));
		    next_coords[1] = next_rad*sin(phi-atan2(k*dr*sin(alpha),radii-k*dr*cos(alpha)));
		    next_coords[2] = z;
		    next_comp = component(next_coords,front->interf);
		}
		distance_4 = sqrt(sqr(next_coords[0]-coords[0])+
		                  sqr(next_coords[1]-coords[1])+
		                  sqr(next_coords[2]-coords[2]));
		distance = min(min(distance1_z,distance2_z),min(distance3_z,distance4_z));
		distance = min(min(distance_1,distance_2),distance);
		distance = min(min(distance_3,distance_4),distance);

		distance_min = min(distance,distance_min);
		}

		fprintf(CLS,"%- 22.16f",coords[0]);
		fprintf(CLS,"%- 22.16f",coords[1]);
		fprintf(CLS,"%- 22.16f",coords[2]);
		fprintf(CLS,"%- 22.16f",distance_min);
		fprintf(CLS,"\n");
		}
		}
	    }
	    fclose(CLS);
	}
}	
*/

