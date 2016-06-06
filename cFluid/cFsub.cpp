/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


#include "cFluid.h"

	/*  Function Declarations */
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static boolean g_intrp_state(STATE*,STATE*,POINT*,HYPER_SURF*,Front*);
LOCAL boolean find_block(double*,int*,double*,RECT_GRID*);
LOCAL void sample_state_1d(Locstate,double*,int,double*,int*,Front*);
LOCAL void sample_state_2d(Locstate,double*,int,double*,int*,Front*);
LOCAL void sample_state_3d(Locstate,double*,int,double*,int*,Front*);
LOCAL void g_intrp_state_1d(int*,double*,double*,POINT*,HYPER_SURF*,Front*);
LOCAL void g_intrp_state_2d(int*,double*,double*,POINT*,HYPER_SURF*,Front*);
LOCAL void g_intrp_state_3d(int*,double*,double*,POINT*,HYPER_SURF*,Front*);
static void promptForDirichletBdryState(FILE*,Front*,HYPER_SURF**,int,int);

static void cfluid_compute_force_and_torque2d(Front*,HYPER_SURF*,double,
                        double*,double*);
static void cfluid_compute_force_and_torque3d(Front*,HYPER_SURF*,double,
                        double*,double*);
static boolean force_on_hse(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,double*,
                                        double*,double*,boolean);
static boolean force_on_hse2d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static boolean force_on_hse3d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static double intrp_between(double,double,double,double,double);

static double (*getStateVel[MAXD])(Locstate) =
               {getStateXvel,getStateYvel,getStateZvel};
static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};
static void set_state_max_speed(Front*,STATE*,double*);
static void get_variable_bdry_params(int,FILE*,POINTER*);
static void cF_variableBoundaryState2d(double*,HYPER_SURF*,Front*,
					POINTER,POINTER);
static void cF_variableBoundaryState3d(double*,HYPER_SURF*,Front*,
					POINTER,POINTER);

double getStateDens(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->dens;
}	/* end getStateDens */

double getStatePdens0(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->pdens[0];
}       /* end getStatePdens */

double getStatePdens1(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->pdens[1];
}       /* end getStatePdens */

double getStateEngy(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->engy;
}	/* end getStateEngy */

double getStatePres(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->pres;
}	/* end getStatePres */

double getStateVort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort;
}	/* end getStateVort */

double getStateXmom(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->momn[0];
}	/* end getStateXmom */

double getStateYmom(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->momn[1];
}	/* end getStateYmom */

double getStateZmom(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->momn[2];
}	/* end getStateZmom */

double getStateXvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[0];
}	/* end getStateXvel */

double getStateYvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[1];
}	/* end getStateYvel */

double getStateZvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[2];
}	/* end getStateZvel */

double getStateXvort(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->vort3d[0];
}       /* end getStateXvort */

double getStateYvort(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->vort3d[1];
}       /* end getStateYvort */

double getStateZvort(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->vort3d[2];
}       /* end getStateZvort */

void read_dirichlet_bdry_data(
	char *inname,
	Front *front)
{
	char msg[100];
	int i,j,k,nhs,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	HYPER_SURF *hs,**hss;
	INTERFACE *intfc = front->interf;
	int i_hs = 0;

	for (i = 0; i < dim; ++i)
	for (j = 0; j < 2; ++j)
	{
	    if (rect_boundary_type(intfc,i,j) == DIRICHLET_BOUNDARY)
	    {
		hs = FT_RectBoundaryHypSurf(intfc,DIRICHLET_BOUNDARY,i,j);
		if (hs == NULL)
		{
		    printf("ERROR: cannot find Dirichlet boundary"
			   " in dimension %d direction %d\n",i,j);
		    clean_up(ERROR);
		}
		if (j == 0)
		    sprintf(msg,"For lower boundary in %d-th dimension",i);
		else
		    sprintf(msg,"For upper boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		promptForDirichletBdryState(infile,front,&hs,1,i_hs);
		i_hs++;
	    }
	    else if (rect_boundary_type(intfc,i,j) == MIXED_TYPE_BOUNDARY)
	    {
		hss = FT_MixedBoundaryHypSurfs(intfc,i,j,DIRICHLET_BOUNDARY,
					&nhs);
		printf("Number of Dirichlet boundaries on dir %d side %d: %d\n",
					i,j,nhs);
		if (dim == 2)
		{
		    for (k = 0; k < nhs; ++k)
		    {
			CURVE *c = Curve_of_hs(hss[k]);
		    	(void) printf("Curve %d start and end at: ",k+1);
		    	(void) printf("(%f %f)->(%f %f)\n",
				  Coords(c->start->posn)[0],
				  Coords(c->start->posn)[1],
				  Coords(c->end->posn)[0],
				  Coords(c->end->posn)[1]);
			promptForDirichletBdryState(infile,front,hss+k,1,i_hs);
			i_hs++;
		    }
		}
	    }
	}
	hss = FT_InteriorHypSurfs(intfc,DIRICHLET_BOUNDARY,&nhs);
	if (nhs == 0 || hss == NULL) return;

	sprintf(msg,"For interior Dirichlet boundary:");
	CursorAfterString(infile,msg);
	(void) printf("\n");
	promptForDirichletBdryState(infile,front,hss,nhs,i_hs);
	i_hs++;
}	/* end read_dirichlet_bdry_data */

void cF_variableBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    cF_variableBoundaryState2d(p0,hs,front,params,state);
	    return;
	case 3:
	    cF_variableBoundaryState3d(p0,hs,front,params,state);
	    return;
	}
}	/* end cF_variableBoundaryState */

void cF_variableBoundaryState2d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	VAR_BDRY_PARAMS *bdry_params;
	STATE *newst = (STATE*) state;
	int i=0, dim, nbr_pist;
	double *angles,half_angular_width,*center,jet_duration_time;
	double radius,theta,vec[MAXD];
	boolean within_piston = NO;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;

	bdry_params = (VAR_BDRY_PARAMS*)boundary_state_function_params(hs);
	jet_duration_time = bdry_params->jet_duration_time;
	if (front->time > jet_duration_time)
	{
	    cF_flowThroughBoundaryState(p0,hs,front,params,state);
	    return;
	}

	dim = bdry_params->dim;
	center = bdry_params->center;
	nbr_pist = bdry_params->number_pistons;
	half_angular_width = bdry_params->half_angular_width;
	angles = bdry_params->angles_pistons;

	radius = 0.0;
	for (i = 0; i < dim; ++i) 
	{
	    vec[i] = p0[i] - center[i];
	    radius += sqr(vec[i]);
	}
	radius = sqrt(radius);
	for (i = 0; i < dim; ++i) 
	    vec[i] /= -radius;
	theta = asin(fabs(p0[1] - center[1])/radius);
	if (p0[0]-center[0] < 0 && p0[1]-center[1] > 0)
            theta = PI - theta;
	else if (p0[0]-center[0] < 0 && p0[1]-center[1] < 0)
            theta = PI + theta;
	else if (p0[0]-center[0] > 0 && p0[1]-center[1] < 0)
            theta = 2*PI - theta;
	for (i = 0; i < nbr_pist; ++i)
	{
	    if (theta > angles[i] - half_angular_width &&
		theta < angles[i] + half_angular_width)
	    {
		within_piston = YES;
	    }
	}
	if (within_piston)
	{
	    POINT *oldp = ft_params->oldp;
	    HYPER_SURF *oldhs = oldp->hs;
	    HYPER_SURF_ELEMENT *oldhse = oldp->hse;
	    STATE *sl,*sr;
	    COMPONENT comp;
	    slsr(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	    if (gas_comp(negative_component(oldhs)))
	    {
	        newst = (STATE*)state;
	        comp = negative_component(oldhs);
	        newst->eos = sl->eos;
	    }
	    else if (gas_comp(positive_component(oldhs)))
	    {
	        newst = (STATE*)state;
	        comp = positive_component(oldhs);
	        newst->eos = sr->eos;
	    }

	    newst->dens = bdry_params->bdry_dens;
	    newst->pres = bdry_params->bdry_pres;
	    for (i = 0; i < dim; ++i)
	    {
		newst->vel[i] = bdry_params->bdry_vel*vec[i];
		newst->momn[i] = (newst->dens)*(newst->vel[i]);
	    }
	    newst->engy = EosEnergy(newst);
	    set_state_max_speed(front,newst,p0);
	}
	else
	{
	    cF_flowThroughBoundaryState(p0,hs,front,params,state);
	}
}	/* end cF_variableBoundaryState2d */

void cF_variableBoundaryState3d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
}	/* end cF_variableBoundaryState3d */

void cF_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	Tan_stencil **tsten;
	Nor_stencil *nsten;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	EQN_PARAMS *eqn_params = ft_params->eqn_params;
	static SWEEP *st_stencil;
	static FSWEEP *st_flux;
	double dir[MAXD];
	double u[3];		/* velocity in the sweeping direction */
	double v[3][MAXD];	/* velocity in the orthogonal direction */
	double vort[3];		/* vorticity stencil */
	double pres[3];		/* pressure stencil */
	double dens[3];		/* pressure stencil */
	double f_u;		/* u flux in the sweeping direction */
	double f_v[MAXD];	/* v flux in the orthogonal direction */
	double f_vort;		/* vort flux */
	double f_pres;		/* pressure flux */
	double f_dens;		/* density flux */
	double dn,dt = front->dt;
	STATE *newst = (STATE*)state;
	STATE  *s0,*sl,*sr,**sts;
	static STATE *s1;
	int i,j,dim = front->rect_grid->dim;
	int nrad = 3;
	int size = 2*nrad + 1;
	
	if (debugging("flow_through"))
	    printf("Entering cF_flowThroughBoundaryState()\n");
	if (s1 == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&s1,sizeof(STATE));
	    FT_ScalarMemoryAlloc((POINTER*)&st_stencil,sizeof(SWEEP));
	    FT_ScalarMemoryAlloc((POINTER*)&st_flux,sizeof(FSWEEP));
	    FT_VectorMemoryAlloc((POINTER*)&st_stencil->dens,size,
					sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&st_stencil->pdens,2,size,
                                        sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&st_stencil->engy,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&st_stencil->pres,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&st_stencil->momn,MAXD,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&st_flux->dens_flux,size,
					sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&st_flux->pdens_flux,2,size,
                                        sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&st_flux->engy_flux,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&st_flux->momn_flux,MAXD,size,
					sizeof(double));
	}

	tsten = FrontGetTanStencils(front,oldp,nrad);
	if (tsten != NULL)
	{
	    for (i = 0; i < dim; ++i)
	    	dir[i] = tsten[0]->dir[i];
	    dn = FT_GridSizeInDir(dir,front);

	    if (comp == negative_component(hs))  
	    {
	        sts = (STATE**)tsten[0]->leftst;
		s0 = sts[0];
	    }
	    else 
	    {
	        sts = (STATE**)tsten[0]->rightst;
		s0 = sts[0];
	    }

	    if (debugging("flow_through"))
	    {
	    	(void) printf("Ambient component: %d\n",comp);
	    	(void) printf("hs = %p  oldp->hs = %p\n",
					(POINTER)hs,(POINTER)oldp->hs);
	    	(void) printf("Time step = %f  Tangential grid size = %f\n",
					dt,dn);
	    	(void) printf("Tangential direction: ");
	    	for (j = 0; j < dim; ++j)
		    (void) printf("%f ",tsten[0]->dir[j]);
	    	(void) printf("\n");
	    	(void) printf("Tan_stencil at point p(%f %f)\n",Coords(oldp)[0],
				Coords(oldp)[1]);
	    	(void) printf("Left points:\n");
	    	for (i = 0; i < nrad; ++i)
	    	{
		    for (j = 0; j < dim; ++j)
	    	    	(void) printf("%f ",Coords(tsten[0]->p[-i])[j]);
		    (void) printf("\n");
	    	}
	    	(void) printf("Right points:\n");
	    	for (i = 0; i < nrad; ++i)
	    	{
		    for (j = 0; j < dim; ++j)
	    	    	(void) printf("%f ",Coords(tsten[0]->p[i])[j]);
		    (void) printf("\n");
	    	}
	    }

	    for (j = 0; j < 3; ++j)
	    	u[j] = 0.0;
	    for (j = 0; j < 3; ++j)
	    {
//	    	vort[j] = sts[j-1]->vort;
	    	pres[j] = sts[j-1]->pres;
	    	dens[j] = sts[j-1]->dens;
	    	for (i = 0; i < dim; ++i)
	    	{
		    u[j] += sts[j-1]->vel[i]*dir[i];
		    v[j][i] = sts[j-1]->vel[i]*(1.0 - dir[i]);
	    	}
	    }

	    f_u = burger_flux(u[0],u[1],u[2]);
	    for (i = 0; i < dim; ++i)
	    	f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
//	    f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	    f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	    f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

	    for (i = 0; i < dim; ++i)
	    	newst->vel[i] = sts[0]->vel[i] - dt/dn*(f_u*dir[i] + f_v[i]) ;
//	    newst->vort = sts[0]->vort - dt/dn*f_vort;
	    newst->pres = sts[0]->pres - dt/dn*f_pres;
	    newst->dens = sts[0]->dens - dt/dn*f_dens;
	}
	else
	{
	    slsr(oldp,oldp->hse,oldp->hs,(POINTER*)&sl,(POINTER*)&sr);
	    if (comp == negative_component(hs))  
		s0 = sl;
	    else
		s0 = sr;
	}
	
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(dir,front);

	if (debugging("flow_through"))
	{
	    printf("Time step = %f  Normal grid size = %f\n",dt,dn);
	    printf("Normal direction: ");
	    for (j = 0; j < dim; ++j)
		printf("%f ",nsten->nor[j]);
	    printf("\n");
	    printf("Nor_stencil at point p(%f %f)\n",Coords(oldp)[0],
				Coords(oldp)[1]);
	    printf("Nor_stencil:\n");
	    for (i = 0; i < nrad; ++i)
	    {
		for (j = 0; j < dim; ++j)
	    	    printf("%f ",nsten->pts[i][j]);
		printf("\n");
	    }
	}

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;
	for (j = 0; j < 2; ++j)
	{
//	    vort[j] = s0->vort;
	    pres[j] = s0->pres;
	    dens[j] = s0->dens;
	    for (i = 0; i < dim; ++i)
	    {
		u[j] += s0->vel[i]*dir[i];
		v[j][i] = s0->vel[i]*(1.0 - dir[i]);
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			eqn_params->vel[i],getStateVel[i],&vtmp,&s0->vel[i]);
	    s1->vel[i] = vtmp;
	}
//	if (dim != 1)
//	{
//	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
//			eqn_params->vort,getStateVort,&vort[2],&s0->vort);
//	    s1->vort = vort[2];
//	}
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->pres,
                            getStatePres,&pres[2],&s0->pres);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->dens,
                            getStateDens,&dens[2],&s0->dens);
	s1->pres = pres[2];
	s1->dens = dens[2];
	for (i = 0; i < dim; ++i)
	{
	    u[2] += s1->vel[i]*dir[i];
	    v[2][i] = s1->vel[i] - s1->vel[i]*dir[i];
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
//	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
//	newst->vort += - dt/dn*f_vort;
	newst->pres += - dt/dn*f_pres;
	newst->dens += - dt/dn*f_dens;
//	set_state_max_speed(front,newst,p0);
	if (debugging("flow_through"))
	{
	    printf("flow through boundary state:\n");
	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    printf("Pressure: %f\n",newst->pres);
//	    printf("Vorticity: %f\n",newst->vort);
	}
}       /* end cF_flowThroughBoundaryState */

void cFluid_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	switch(wave_type(oldhs))
	{
	case NEUMANN_BOUNDARY:
	case MOVABLE_BODY_BOUNDARY:
	    return neumann_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case DIRICHLET_BOUNDARY:
	    return dirichlet_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case SUBDOMAIN_BOUNDARY:
	case PASSIVE_BOUNDARY:
            return;
	default:
	    return contact_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	}
}       /* cFluid_point_propagate */

static  void neumann_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double *m_pres = eqn_params->pres;
	double *m_dens = eqn_params->dens;
	double *m_engy = eqn_params->engy;
        double **m_pdens = eqn_params->pdens;
	double nor[MAXD],tan[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (gas_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
            oldst->eos = &eqn_params->eos[comp];
	    newst = (STATE*)left_state(newp);
	}
	else if (gas_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
            oldst->eos = &eqn_params->eos[comp];
	    newst = (STATE*)right_state(newp);
	}
	FT_NormalAtPoint(oldp,front,nor,comp);

	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;
	tan[0] = -nor[1]; 	tan[1] = nor[0];

	if (wave_type(oldhs) == MOVABLE_BODY_BOUNDARY)
	{
            double omega_dt,crds_com[MAXD];
            omega_dt = angular_velo(oldhs)*dt;
            for (i = 0; i < dim; ++i)
            {
                vel[i] = center_of_mass_velo(oldhs)[i];
                crds_com[i] = Coords(oldp)[i] +dt*vel[i] - 
			center_of_mass(oldhs)[i];
            }
            vel[0] += -angular_velo(oldhs)*crds_com[1]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[0]*sin(omega_dt);
            vel[1] +=  angular_velo(oldhs)*crds_com[0]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[1]*sin(omega_dt);
	}
	else
	{
            for (i = 0; i < dim; ++i)
	    	vel[i] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	{
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
	    newst->vel[i] = vel[i];
//            newst->momn[i] = newst->dens*vel[i];  //HK FIX
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
	}
	FT_IntrpStateVarAtCoords(front,comp,p1,m_dens,
			getStateDens,&newst->dens,&oldst->dens);
	FT_IntrpStateVarAtCoords(front,comp,p1,m_engy,
			getStateEngy,&newst->engy,&oldst->engy);
        if(eqn_params->multi_comp_non_reactive == YES)
        {
            //int ii;
            //for(ii = 0; ii < eqn_params->n_comps; ii++)
            {
                FT_IntrpStateVarAtCoords(front,comp,p1,m_pdens[0],
                        getStatePdens0,&newst->pdens[0],&oldst->pdens[0]);
                FT_IntrpStateVarAtCoords(front,comp,p1,m_pdens[1],
                        getStatePdens1,&newst->pdens[1],&oldst->pdens[1]);
            }
        }
        newst->eos = oldst->eos;
	newst->pres = EosPressure(newst);
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	    newst->momn[i] = newst->dens*vel[i];
	}
	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
	set_state_max_speed(front,newst,Coords(newp));
}	/* end neumann_point_propagate */

static  void dirichlet_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        int i, dim = front->rect_grid->dim;
	STATE *sl,*sr,*newst = NULL;
	STATE *bstate;
	FLOW_THROUGH_PARAMS ft_params;
	COMPONENT comp;

	if (debugging("dirichlet_bdry"))
	{
	    printf("Entering dirichlet_point_propagate()\n");
	    print_general_vector("oldp:  ",Coords(oldp),dim,"\n");
	}

	slsr(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	if (gas_comp(negative_component(oldhs)))
	{
	    newst = (STATE*)left_state(newp);
	    comp = negative_component(oldhs);
	    newst->eos = sl->eos;
	}
	else if (gas_comp(positive_component(oldhs)))
	{
	    newst = (STATE*)right_state(newp);
	    comp = positive_component(oldhs);
	    newst->eos = sr->eos;
	}
	if (newst == NULL) return;	// node point

	if (boundary_state(oldhs) != NULL)
	{
	    bstate = (STATE*)boundary_state(oldhs);
	    newst->dens = bstate->dens;
	    newst->engy = bstate->engy;
            for (i = 0; i < dim; ++i)
	    	newst->vel[i] = bstate->vel[i];
            newst->pres = bstate->pres;
            newst->eos = bstate->eos;
            newst->vort = 0.0;
	    set_state_max_speed(front,newst,Coords(oldp));

	    if (debugging("dirichlet_bdry"))
	    {
		printf("Preset boundary state:\n");
		print_general_vector("Velocity: ",newst->vel,dim,"\n");
		printf("Density: %f\n",newst->dens);
		printf("Energy: %f\n",newst->engy);
		printf("Pressure: %f\n",newst->pres);
		printf("Vorticity: %f\n",newst->vort);
	    }
	}
	else if (boundary_state_function(oldhs))
	{
	    oldp->hse = oldhse;
	    oldp->hs = oldhs;
	    ft_params.oldp = oldp;
	    ft_params.eqn_params = eqn_params;
	    ft_params.comp = comp;
	    (*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
			(POINTER)&ft_params,(POINTER)newst);	
	}
	if (debugging("dirichlet_bdry"))
	    printf("Leaving dirichlet_point_propagate()\n");
        return;
}	/* end dirichlet_point_propagate */

static  void contact_point_propagate(
        Front			*front,
        POINTER			wave,
        POINT			*oldp,
        POINT			*newp,
        HYPER_SURF_ELEMENT	*oldhse,
        HYPER_SURF		*oldhs,
        double			dt,
        double			*V)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        int i, dim = front->rect_grid->dim;
	double **m_mom = eqn_params->mom;
	double *m_dens = eqn_params->dens;
	double *m_engy = eqn_params->engy;
	double *m_pres = eqn_params->pres;
	double *p0;
	STATE *oldst,*newst;
	STATE *newstl, *newstr;
	POINTER sl,sr;
	EOS_PARAMS *eos = eqn_params->eos;
	double default_var;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);

        newstl = (STATE*)left_state(newp);
        newstl->eos = &eos[negative_component(oldhs)]; 

        newstr = (STATE*)right_state(newp);
        newstr->eos = &eos[positive_component(oldhs)];

	STATE *stl, *str;
	STATE *ssl, *ssr;
	FT_ScalarMemoryAlloc((POINTER*)&stl,sizeof(STATE));
	FT_ScalarMemoryAlloc((POINTER*)&str,sizeof(STATE));
	FT_ScalarMemoryAlloc((POINTER*)&ssl,sizeof(STATE));
	FT_ScalarMemoryAlloc((POINTER*)&ssr,sizeof(STATE));
        *stl = *(STATE*)sl;
        *str = *(STATE*)sr;
        *ssl = *(STATE*)sl;
        *ssr = *(STATE*)sr;
        stl->eos = newstl->eos;
        str->eos = newstr->eos;

	if (!g_intrp_state(stl,str,oldp,oldhs,front))
	{
	    //if can't find block for the point, keep the old point
	    *newp = *oldp;
	    return;
	}

        double lvnor, rvnor;
        double lvtang[MAXD], rvtang[MAXD];
        double nor[MAXD];
        double velp[MAXD],s;
        const double epsilon = 1.e-10;
        const double delta = 1.e-10;

        normal(oldp,oldhse,oldhs,nor,front);

        lvnor = 0;
        rvnor = 0;
        for (i=0; i<dim; i++)
        {
                lvnor += stl->vel[i]*nor[i];
                rvnor += str->vel[i]*nor[i];
        }
        for (i=0; i<dim; i++)
        {
                lvtang[i] = stl->vel[i]-lvnor*nor[i];
                rvtang[i] = str->vel[i]-rvnor*nor[i];
        }

	double pml, pmr, uml, umr, ml, mr;
	RIEMANN_SOLVER_WAVE_TYPE l_wave,r_wave;
	STATE *stl2, *str2;
	STATE *ansl, *ansr;

	FT_ScalarMemoryAlloc((POINTER*)&stl2,sizeof(STATE));
	FT_ScalarMemoryAlloc((POINTER*)&str2,sizeof(STATE));
	*stl2 = *stl;
	*str2 = *str;
        stl2->eos = newstl->eos;
        str2->eos = newstr->eos;

	stl2->vel[0] = lvnor;
	str2->vel[0] = rvnor;
	for (i = 1; i < dim; i++)
	    stl2->vel[i] = str2->vel[i] = 0.0;

	FT_ScalarMemoryAlloc((POINTER*)&ansl,sizeof(STATE));
	FT_ScalarMemoryAlloc((POINTER*)&ansr,sizeof(STATE));

	if (!find_mid_state(stl2,str2,0.0/*pjump*/,&pml,&pmr,&uml,&umr,&ml,&mr,
		    &l_wave,&r_wave))
	{
	    printf("In contact_point_propagate(), find_mid_state() failed at (%lf, %lf, %lf).\n", 
		    Coords(oldp)[0], Coords(oldp)[1], Coords(oldp)[2]);
	    printf("nor: (%lf, %lf, %lf).\n", nor[0], nor[1], nor[2]);
	    printf("ssl: dens = %e, pres = %e, vel = %e.\n", 
		    ssl->dens, ssl->pres, ssl->vel[0]);
	    printf("ssr: dens = %e, pres = %e, vel = %e.\n", 
		    ssr->dens, ssr->pres, ssr->vel[0]);
	    printf("stl: dens = %e, pres = %e, vel = %e.\n", 
		    stl2->dens, stl2->pres, stl2->vel[0]);
	    printf("str: dens = %e, pres = %e, vel = %e.\n", 
		    str2->dens, str2->pres, str2->vel[0]);
	    clean_up(ERROR);
	}
	midstate(stl2,ansl,ml,uml,pml,l_wave,1);
	midstate(str2,ansr,mr,umr,pmr,r_wave,-1);	//1 for left, -1 for right

        for (i=0; i<dim; i++)
            velp[i] = ansl->vel[0]*nor[i] + lvtang[i];

        newstl->dens = ansl->dens;
        newstl->pres = ansl->pres;
        newstr->dens = ansr->dens;
        newstr->pres = ansr->pres;
        if(eqn_params->multi_comp_non_reactive == YES)
        {
            int ii;
            for(ii = 0; ii < eqn_params->n_comps; ii++)
            {
                newstl->pdens[ii] = ansl->pdens[ii];
                newstr->pdens[ii] = ansr->pdens[ii];
            }
        }
        for (i=0; i<dim; i++)
        {
            newstl->vel[i] = velp[i];
            newstl->momn[i] = newstl->vel[i]*newstl->dens;
            newstr->vel[i] = velp[i];
            newstr->momn[i] = newstr->vel[i]*newstr->dens;
        }
        newstl->engy = EosEnergy(newstl);
        newstr->engy = EosEnergy(newstr);

	FT_FreeThese(2,stl2,str2);
	FT_FreeThese(2,ansl,ansr);

        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(oldp)[i] + dt*velp[i];
            set_max_front_speed(i,fabs(velp[i]),NULL,Coords(newp),front);
        }
        s = mag_vector(velp,dim);

	double cl, cr, c;
	cl = EosSoundSpeed(newstl);
	cr = EosSoundSpeed(newstr);
	c = std::max(cl, cr);
	set_max_front_speed(dim,s+c,NULL,Coords(newp),front);
	FT_FreeThese(4,stl,str,ssl,ssr);
}	/* end contact_point_propagate */

static boolean g_intrp_state(
	STATE		*stl, 
	STATE		*str, 
	POINT		*p, 
	HYPER_SURF	*hs,
	Front		*front)
{
	RECT_GRID *gr = computational_grid(front->interf);
	int i, j, k, dim = gr->dim;
	double *L = gr->L;
	double *U = gr->U;
	double *h = gr->h;
	double *coords = Coords(p);
	int icrds[MAXD];
	double f[MAXD];
	double V[MAXD];

	if (!find_block(f, icrds, coords, gr))
	{
	    return NO;
	}

	for (i = 0; i < dim; i++)
		icrds[i] += (gr->lbuf[i] == 0) ? 1 : gr->lbuf[i];

	switch(dim)
	{
		case 1:
			g_intrp_state_1d(icrds, V, f, p, hs, front);
			break;
		case 2:
			g_intrp_state_2d(icrds, V, f, p, hs, front);
			break;
		case 3:
			g_intrp_state_3d(icrds, V, f, p, hs, front);
			break;
	}

	*stl = *(STATE*)left_state(p);
	*str = *(STATE*)right_state(p);

	return YES;
}

LOCAL	void g_intrp_state_1d(
	int		*icrds,
	double		*V,
	double		*f,
	POINT		*p,
	HYPER_SURF	*hs,
	Front		*front)
{
	if(negative_component(hs) == 2)
	    sample_state_1d(left_state(p),V,0,f,icrds,front);
	else
	    sample_state_1d(left_state(p),V,1,f,icrds,front);
	
	if(positive_component(hs) == 3)
	    sample_state_1d(right_state(p),V,1,f,icrds,front);
	else
	    sample_state_1d(right_state(p),V,0,f,icrds,front);
}

LOCAL	void g_intrp_state_2d(
	int		*icrds,
	double		*V,
	double		*f,
	POINT		*p,
	HYPER_SURF	*hs,
	Front		*front)
{
	if(negative_component(hs) == 2)
	    sample_state_2d(left_state(p),V,0,f,icrds,front);
	else
	    sample_state_2d(left_state(p),V,1,f,icrds,front);
	
	if(positive_component(hs) == 3)
	    sample_state_2d(right_state(p),V,1,f,icrds,front);
	else
	    sample_state_2d(right_state(p),V,0,f,icrds,front);
}

LOCAL	void g_intrp_state_3d(
	int		*icrds,
	double		*V,
	double		*f,
	POINT		*p,
	HYPER_SURF	*hs,
	Front		*front)
{
	if(negative_component(hs) == 2)
	    sample_state_3d(left_state(p),V,0,f,icrds,front);
	else
	    sample_state_3d(left_state(p),V,1,f,icrds,front);
	
	if(positive_component(hs) == 3)
	    sample_state_3d(right_state(p),V,1,f,icrds,front);
	else
	    sample_state_3d(right_state(p),V,0,f,icrds,front);
}

LOCAL	void sample_state_1d(
	Locstate	ans,
	double		*v,
	int		ind,
	double		*f,
	int		*icrds,
	Front		*front)
{
	EQN_PARAMS	*eqn_params = (EQN_PARAMS*)front->extra1;
	double		***Gv = eqn_params->Gvel;
	double		**Gd = eqn_params->Gdens;
	double		**Gp = eqn_params->Gpres;
	int		index1, index2;

	index1 = icrds[0];
	index2 = icrds[0]+1;
	((STATE*)ans)->dens = Gd[ind][index1]*(1.0-f[0]) + Gd[ind][index2]*f[0];
	((STATE*)ans)->pres = Gp[ind][index1]*(1.0-f[0]) + Gp[ind][index2]*f[0];
	v[0] = Gv[ind][0][index1]*(1.0-f[0]) + Gv[ind][0][index2]*f[0];
	((STATE*)ans)->vel[0] = v[0];
}

LOCAL	void sample_state_2d(
	Locstate	ans,
	double		*v,
	int		ind,
	double		*f,
	int		*icrds,
	Front		*front)
{
	RECT_GRID	*gr = &topological_grid(front->grid_intfc);
	int		dim = gr->dim;
	int		*gmax = gr->gmax;
	EQN_PARAMS	*eqn_params = (EQN_PARAMS*)front->extra1;
	double		***Gv = eqn_params->Gvel;
	double		**Gd = eqn_params->Gdens;
	double		**Gp = eqn_params->Gpres;
	int		index1, index2;
	int		i, ic[2], k;
	double		densz[2], presz[2], velz[2][2];

	for (i=0; i<2; i++)
	{
		ic[0] = icrds[0]+i;
		ic[1] = icrds[1];
		index1 = d_index(ic,gmax,dim);
		ic[1]++;
		index2 = d_index(ic,gmax,dim);

		densz[i] = Gd[ind][index1]*(1.0-f[1]) + Gd[ind][index2]*f[1];
		presz[i] = Gp[ind][index1]*(1.0-f[1]) + Gp[ind][index2]*f[1];
                for(k=0; k<2; k++)
                    velz[i][k] = Gv[ind][k][index1]*(1.0-f[1]) + Gv[ind][k][index2]*f[1];

	}
	((STATE*)ans)->dens = densz[0]*(1.0-f[0]) + densz[1]*f[0];
	((STATE*)ans)->pres = presz[0]*(1.0-f[0]) + presz[1]*f[0];
        for(k=0; k<3; k++)
	{
		v[k] = velz[0][k]*(1.0-f[0]) + velz[1][k]*f[0];
		((STATE*)ans)->vel[k] = v[k];
	}
}

LOCAL	void sample_state_3d(
	Locstate	ans,
	double		*v,
	int		ind,
	double		*f,
	int		*icrds,
	Front		*front)
{
	int		i, j, k, ic[3],iii;
	Locstate	s1, s2;
	double		velz[2][2][3], presz[2][2], densz[2][2], pdensz[2][2][10], sum;
	int		st_type;
	RECT_GRID	*gr = &topological_grid(front->grid_intfc);
	int		dim = gr->dim;
	int		*gmax = gr->gmax;
	EQN_PARAMS	*eqn_params = (EQN_PARAMS*)front->extra1;
	double		***Gv = eqn_params->Gvel;
	double		**Gd = eqn_params->Gdens;
        double          ***Gpd = eqn_params->Gpdens;
	double		**Gp = eqn_params->Gpres;
	int		index1, index2;

	for(i=0; i<2; i++)
	    for(j=0; j<2; j++)
	    {
		ic[0] = icrds[0] + i;
		ic[1] = icrds[1] + j;
		ic[2] = icrds[2];
		index1 = d_index(ic,gmax,dim);
		
		ic[2]++;
		index2 = d_index(ic,gmax,dim);
	
                for(k=0; k<3; k++)
                    velz[i][j][k] = Gv[ind][k][index1]*(1.0-f[2]) + Gv[ind][k][index2]*f[2];

		presz[i][j] = Gp[ind][index1]*(1.0-f[2]) + Gp[ind][index2]*f[2];
		densz[i][j] = Gd[ind][index1]*(1.0-f[2]) + Gd[ind][index2]*f[2];	       
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        pdensz[i][j][ii] = Gpd[ind][ii][index1]*(1.0-f[2]) + Gpd[ind][ii][index2]*f[2];
                    }
                }
	    }

        for(k=0; k<3; k++)
        {
            velz[0][0][k] = velz[0][0][k]*(1.0-f[1]) + velz[0][1][k]*f[1];
            velz[1][0][k] = velz[1][0][k]*(1.0-f[1]) + velz[1][1][k]*f[1];
            v[k] = velz[0][0][k]*(1.0-f[0]) + velz[1][0][k]*f[0];
        }

	presz[0][0] = presz[0][0]*(1.0-f[1]) + presz[0][1]*f[1];
	presz[1][0] = presz[1][0]*(1.0-f[1]) + presz[1][1]*f[1];
	((STATE*)ans)->pres = presz[0][0]*(1.0-f[0]) + presz[1][0]*f[0];
	
	densz[0][0] = densz[0][0]*(1.0-f[1]) + densz[0][1]*f[1];
	densz[1][0] = densz[1][0]*(1.0-f[1]) + densz[1][1]*f[1];
	((STATE*)ans)->dens = densz[0][0]*(1.0-f[0]) + densz[1][0]*f[0];

        if(eqn_params->multi_comp_non_reactive == YES)
        {
            int ii;
            for(ii = 0; ii < eqn_params->n_comps; ii++)
            {
                pdensz[0][0][ii] = pdensz[0][0][ii]*(1.0-f[1]) + pdensz[0][1][ii]*f[1];
                pdensz[1][0][ii] = pdensz[1][0][ii]*(1.0-f[1]) + pdensz[1][1][ii]*f[1];
            }
            for(ii = 0; ii < eqn_params->n_comps; ii++)
            {
                ((STATE*)ans)->pdens[ii] = pdensz[0][0][ii]*(1.0-f[0]) + pdensz[1][0][ii]*f[0];
            }
            //normalize
            double sum = 0.0;
            for(ii = 0; ii < eqn_params->n_comps; ii++)
                sum += ((STATE*)ans)->pdens[ii];
            for(ii = 0; ii < eqn_params->n_comps; ii++)
                ((STATE*)ans)->pdens[ii] = ((STATE*)ans)->dens/sum*((STATE*)ans)->pdens[ii];
        }

	((STATE*)ans)->vel[0] = v[0];
	((STATE*)ans)->vel[1] = v[1];
	((STATE*)ans)->vel[2] = v[2];

}

#define TOL 1e-10
#define	corner_index(p,i,gr)	irint(floor(((p)-(gr)->L[i])/(gr)->h[i]-0.5))
LOCAL	boolean find_block(
	double		*f,
	int		*icrds,
	double		*p,
	RECT_GRID	*gr)
{
	int	i;
	int	dim=gr->dim;

	for(i=0; i<dim; i++)
	{
	    icrds[i] = corner_index(p[i],i,gr);
	    if(icrds[i] < -gr->lbuf[i] || icrds[i] >= gr->gmax[i]+gr->ubuf[i]-1)
		return  NO;
	    f[i] = p[i] - (gr->L[i]+(0.5+icrds[i])*gr->h[i]);
	    f[i] /= gr->h[i];
	}
	return  YES;
}

static void set_state_max_speed(
	Front *front,
	STATE *state,
	double *coords)
{
	int i,dim = front->rect_grid->dim;
	double c,s;
	s = 0.0;
	for (i = 0; i < dim; ++i)
	    s += sqr(state->momn[i]/state->dens);
	s = sqrt(s);
	c = EosSoundSpeed(state);
	set_max_front_speed(dim,s+c,NULL,coords,front);
}	/* end set_state_max_speed */

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
	if (a > 0.0)
	    return a*(um - ul);
	else
	    return a*(ur - um);
}	/* end net_uwind_flux */

void readFrontStates(
	Front		*front,
	char		*restart_name)
{
	FILE 		*infile;
	EQN_PARAMS 	*eqn_params = (EQN_PARAMS*)front->extra1;
	INTERFACE 	*intfc = front->interf;
        STATE 		*sl,*sr;
        POINT 		*p;
        HYPER_SURF 	*hs;
        HYPER_SURF_ELEMENT *hse;
	STATE 		*lstate,*rstate;
	char 		fname[100];
	int 		i,dim = front->rect_grid->dim;
	int		comp;
	EOS_PARAMS	*eos = eqn_params->eos;

	sprintf(fname,"%s-gas",restart_name);
	infile = fopen(fname,"r");
	
	/* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface gas states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    lstate = (STATE*)sl;	rstate = (STATE*)sr;
            fscanf(infile,"%lf %lf",&lstate->dens,&rstate->dens);
            fscanf(infile,"%lf %lf",&lstate->engy,&rstate->engy);
	    for (i = 0; i < dim; ++i)
            	fscanf(infile,"%lf %lf",&lstate->momn[i],&rstate->momn[i]);	
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                int ii;
                for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    fscanf(infile,"%lf %lf",&lstate->pdens[ii],&rstate->pdens[ii]);
                }
            }
    
	    comp = negative_component(hs);
	    lstate->eos = &eos[comp];
	    lstate->dim = dim;
	    if(gas_comp(comp))
	    	lstate->pres = EosPressure(lstate);
		
	    comp = positive_component(hs);

	    rstate->eos = &eos[comp];
	    rstate->dim = dim;
	    if(gas_comp(comp))
	    	rstate->pres = EosPressure(rstate);
	    lstate->dim = rstate->dim = dim;
//            fprintf(stdout, "In the readFrontStates %e %e\n",lstate->pres,rstate->pres);
        }
	FT_MakeGridIntfc(front);
	fclose(infile);
}

extern void reflectVectorThroughPlane(
	double *vec,
	double *nor,
	double *vec_ref,
	int dim)
{
	int i;
	double vec_nor[MAXD];
	for (i = 0; i < dim; ++i)
	{
	    vec_nor[i] = vec[i]*fabs(nor[i]);
	    vec_ref[i] = vec[i] - 2.0*vec_nor[i];
	}	
}	/* end reflectVectorThroughPlane */

extern boolean reflectNeumannState(
	Front *front,
	HYPER_SURF *hs,
	double *coords,
	COMPONENT comp,
	SWEEP *m_vst,
	STATE *state)
{
	int i,dim = front->rect_grid->dim;
	double coordsref[MAXD],nor[MAXD];
	double momn[MAXD];

	if (!FrontReflectPointViaNeumannBdry(coords,coordsref,nor,comp,
				hs,front))
	{
	    printf("ERROR: in appendGhostBuffer(), cannot reflect point!\n");
	    return NO;
	}
	FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->dens,getStateDens,
					&state->dens,NULL);
	FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->engy,getStateEngy,
					&state->engy,NULL);
	FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->pres,getStatePres,
					&state->pres,NULL);
	for (i = 0; i < dim; ++i)
	{
	    FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->momn[i],
				getStateMom[i],&momn[i],NULL);
	}
        reflectVectorThroughPlane(momn,nor,state->momn,dim);
	return YES;
}	/* end reflectNeumannState */	

extern void findGhostState(
	STATE intfc_st,
	STATE inter_st,
	STATE *ghost_st)
{
	double vel[MAXD];
	vel[0] = inter_st.momn[0]/inter_st.dens;
	vel[1] = inter_st.momn[1]/inter_st.dens;
	vel[2] = inter_st.momn[2]/inter_st.dens;


	ghost_st->dens = intfc_st.dens;
	ghost_st->pres = intfc_st.pres;
	ghost_st->momn[0] = intfc_st.dens*vel[0];
	ghost_st->momn[1] = intfc_st.dens*vel[1];
	ghost_st->momn[2] = intfc_st.dens*vel[2];
	ghost_st->engy = EosEnergy(ghost_st);
}	/* end findGhostState */

static void promptForDirichletBdryState(
	FILE *infile,
	Front *front,
	HYPER_SURF **hs,
	int nhs,
	int i_hs)
{
	static STATE *state;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	char s[100];
	COMPONENT comp;
	int i,k,dim = front->rect_grid->dim;
	POINTER func_params;

	FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
	state->dim = dim;

	CursorAfterString(infile,"Enter type of Dirichlet boundary:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	switch (s[0])
	{
	case 'c':			// Constant state
	case 'C':
	    comp = gas_comp(positive_component(hs[0])) ? 
				positive_component(hs[0]) :
				negative_component(hs[0]);
	    state->eos = &(eqn_params->eos[comp]);
	    CursorAfterString(infile,"Enter velocity:");
	    for (k = 0; k < dim; ++k)
	    {
		fscanf(infile,"%lf",&state->vel[k]);
		(void) printf("%f ",state->vel[k]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter pressure:");
	    fscanf(infile,"%lf",&state->pres);
	    (void) printf("%f\n",state->pres);
	    CursorAfterString(infile,"Enter density:");
	    fscanf(infile,"%lf",&state->dens);
	    (void) printf("%f\n",state->dens);
	    for (k = 0; k < dim; ++k)
                       state->momn[k] = state->dens*state->vel[k];
	    state->engy = EosEnergy(state);
	    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
			(POINTER)state,hs[0],i_hs);
	    for (i = 1; i < nhs; ++i)
		bstate_index(hs[i]) = bstate_index(hs[0]);
	    break;
	case 'f':			// Flow through state
	case 'F':
	    FT_InsertDirichletBoundary(front,cF_flowThroughBoundaryState,
			"cF_flowThroughBoundaryState",NULL,NULL,hs[0],i_hs);
	    for (i = 1; i < nhs; ++i)
		bstate_index(hs[i]) = bstate_index(hs[0]);
	    break;
	case 'v':			// Flow through state
	case 'V':
	    get_variable_bdry_params(dim,infile,&func_params);
	    FT_InsertDirichletBoundary(front,cF_variableBoundaryState,
			"cF_variableBoundaryState",func_params,NULL,hs[0],i_hs);
	    for (i = 1; i < nhs; ++i)
		bstate_index(hs[i]) = bstate_index(hs[0]);
	    break;
	}
} 	/* end  promptForDirichletBdryState */

/* 	Function computes velocity of center of mass and
 *  	angular velocity of a regid body, must be a closed curve. 
*/

extern	void cfluid_compute_force_and_torque(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	switch (fr->rect_grid->dim)
	{
	case 2:
	    return cfluid_compute_force_and_torque2d(fr,hs,dt,force,torque);
	case 3:
	    return cfluid_compute_force_and_torque3d(fr,hs,dt,force,torque);
	}
}	/* end cfluid_compute_force_and_torque */

static	void cfluid_compute_force_and_torque2d(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	RECT_GRID *gr = computational_grid(fr->interf);
	double f[MAXD],rr[MAXD];
	double t,pres;
	double area[MAXD],posn[MAXD];
	BOND *b;
	boolean pos_side;
	int i,dim = gr->dim;
	EQN_PARAMS *cFparams = (EQN_PARAMS*)fr->extra1;
	double *gravity = cFparams->gravity;
	CURVE *curve = Curve_of_hs(hs);

	if (debugging("rigid_body"))
	    (void) printf("Entering cfluid_compute_force_and_torque2d()\n");

	if (gas_comp(negative_component(curve)))
	    pos_side = NO;
	else 
	    pos_side = YES;

	for (i = 0; i < dim; ++i)
	{
	    force[i] = 0.0;
	}
	*torque = 0.0;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    if (force_on_hse(Hyper_surf_element(b),Hyper_surf(curve),gr,
			&pres,area,posn,pos_side))
	    {
	    	for (i = 0; i < dim; ++i)
	    	{
		    f[i] = pres*area[i];
	    	    rr[i] = 0.5*(Coords(b->start)[i] + Coords(b->end)[i])
				- rotation_center(curve)[i];
	    	    force[i] += f[i];
	    	}
	    	Cross2d(rr,f,t);
	    	*torque += t;
	    }
	}
	 /* Add gravity to the total force */
	if (motion_type(curve) != ROTATION)
	{
	    for (i = 0; i < dim; ++i)
	    	force[i] += gravity[i]*total_mass(curve);
	}
	if (debugging("rigid_body"))
	{
	    (void) printf("Leaving cfluid_compute_force_and_torque2d()\n");
	    (void) printf("total_force = %f %f\n",force[0],force[1]);
	    (void) printf("torque = %f\n",*torque);
	}
}	/* end cfluid_compute_force_and_torque2d */

#define         MAX_TRI_FOR_INTEGRAL            100
static	void cfluid_compute_force_and_torque3d(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	RECT_GRID *gr = computational_grid(fr->interf);
	double f[MAXD],rr[MAXD];
	double t[MAXD],tdir,pres;
	double area[MAXD],posn[MAXD];
	TRI *tri;
	boolean pos_side;
	int i,dim = gr->dim;
	EQN_PARAMS *cFparams = (EQN_PARAMS*)fr->extra1;
	double *gravity = cFparams->gravity;
	SURFACE *surface = Surface_of_hs(hs);

	if (gas_comp(negative_component(surface)))
	    pos_side = NO;
	else 
	    pos_side = YES;

	for (i = 0; i < dim; ++i)
	{
	    force[i] = 0.0;
	    torque[i] = 0.0;
	}
	for (tri = first_tri(surface); !at_end_of_tri_list(tri,surface); 
			tri = tri->next)
	{
	    if (force_on_hse(Hyper_surf_element(tri),Hyper_surf(surface),gr,
			&pres,area,posn,pos_side))
	    {
	    	for (i = 0; i < dim; ++i)
	    	{
		    f[i] = pres*area[i];
	    	    force[i] += f[i];
		    rr[i] = posn[i] - rotation_center(surface)[i];
		}
		Cross3d(rr,f,t);
		tdir = Dot3d(t,(rotation_direction(hs)));
	    	for (i = 0; i < dim; ++i)
		{
		    t[i] = tdir*rotation_direction(hs)[i];
		    torque[i] += t[i];
		}
	    }
	}
	 /* Add gravity to the total force */
	if (motion_type(surface) != ROTATION)
	{
	    for (i = 0; i < dim; ++i)
	    	force[i] += gravity[i]*total_mass(surface);
	}
	if (debugging("rigid_body"))
	{
	    printf("In cfluid_compute_force_and_torque3d()\n");
	    printf("total_force = %f %f %f\n",force[0],force[1],force[2]);
	    printf("torque = %f %f %f\n",torque[0],torque[1],torque[2]);
	}
}	/* end cfluid_compute_force_and_torque3d */


static boolean force_on_hse(
	HYPER_SURF_ELEMENT *hse,	/* Bond (2D) or tri (3D) */
	HYPER_SURF *hs,			/* Curve (2D) or surface (3D) */
	RECT_GRID *gr,			/* Rectangular grid */
	double *pres,		/* Average pressure */
	double *area,		/* Area as a vector, pointing onto body */
	double *posn,		/* Position of the pressure */
	boolean pos_side)	/* Is the body on the positive side of hs? */
{
	int dim = gr->dim;
	switch (dim)
	{
	case 2: 
	    return force_on_hse2d(hse,hs,gr,pres,area,posn,pos_side);
	case 3: 
	    return force_on_hse3d(hse,hs,gr,pres,area,posn,pos_side);
	default: 
	    return NO; 
	}
	
}	/* end force_on_hse */

static boolean force_on_hse2d(
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	RECT_GRID *gr,
	double *pres,
	double *area,
	double *posn,
	boolean pos_side)
{
	double crds1[MAXD],crds2[MAXD];
	double p1,p2;
	Locstate s1,s2;
	BOND *b = Bond_of_hse(hse);
	CURVE *c = Curve_of_hs(hs);
	double *L = gr->L;
	double *U = gr->U;
	int i;
	
	/* Get pressure at two end points of the bond */
	if (b->start == c->start->posn)
	    s1 = pos_side ? right_start_state(c) : left_start_state(c);
	else
	    s1 = pos_side ? right_state(b->start) : left_state(b->start);
	if (b->end == c->end->posn)
	    s2 = pos_side ? right_end_state(c) : left_end_state(c);
	else
	    s2 = pos_side ? right_state(b->end) : left_state(b->end);

	p1 = getStatePres(s1);	p2 = getStatePres(s2);
	for (i = 0; i < 2; ++i)
	{
	    crds1[i] = Coords(b->start)[i];
	    crds2[i] = Coords(b->end)[i];
	}

	/* Cut and interpolate if one end is outside the domain */
	for (i = 0; i < 2; ++i)
	{
	    if (crds1[i] <= L[i])
	    {
		if (crds2[i] <= L[i]) return NO; // both ends out
		else
		{
		    crds1[(i+1)%2] = intrp_between(crds1[i],crds2[i],L[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p1 = intrp_between(crds1[i],crds2[i],L[i],p1,p2);
		    crds1[i] = L[i];
		}
	    }
	    if (crds1[i] >= U[i])
	    {
		if (crds2[i] >= U[i]) return NO; // both ends out
		else
		{
		    crds1[(i+1)%2] = intrp_between(crds1[i],crds2[i],U[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p1 = intrp_between(crds1[i],crds2[i],U[i],p1,p2);
		    crds1[i] = U[i];
		}
	    }
	}
	for (i = 0; i < 2; ++i)
	{
	    if (crds2[i] <= L[i])
	    {
		if (crds1[i] <= L[i]) return NO; // both ends out
		else
		{
		    crds2[(i+1)%2] = intrp_between(crds1[i],crds2[i],L[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p2 = intrp_between(crds1[i],crds2[i],L[i],p1,p2);
		    crds2[i] = L[i];
		}
	    }
	    if (crds2[i] >= U[i])
	    {
		if (crds1[i] >= U[i]) return NO; // both ends out
		else
		{
		    crds2[(i+1)%2] = intrp_between(crds1[i],crds2[i],U[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p2 = intrp_between(crds1[i],crds2[i],U[i],p1,p2);
		    crds2[i] = U[i];
		}
	    }
	}
	area[0] = pos_side ? crds1[1] - crds2[1] : crds2[1] - crds1[1];
	area[1] = pos_side ? crds2[0] - crds1[0] : crds1[0] - crds2[0];
	*pres = 0.5*(p1 + p2);
	posn[0] = 0.5*(crds1[0] + crds2[0]);
	posn[1] = 0.5*(crds1[1] + crds2[1]);
	return YES;
}	/* end force_on_hse2d */

static boolean force_on_hse3d(
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	RECT_GRID *gr,
	double *pres,
	double *area,
	double *posn,
	boolean pos_side)
{
        TRI *t = Tri_of_hse(hse);
	POINT *point;
	Locstate sl,sr;
        int i,j,dim = gr->dim;

	*pres = 0.0;
	for (i = 0; i < 3; ++i)
	    posn[i] = 0.0;
	for (i = 0; i < 3; ++i)
	{
	    point = Point_of_tri(t)[i];
	    for (j = 0; j < dim; ++j)
		posn[j] += Coords(point)[j];
	    FT_GetStatesAtPoint(point,hse,hs,&sl,&sr);
	    if (pos_side)
		*pres += getStatePres(sr);
	    else
		*pres += getStatePres(sl);
	}
	*pres /= 3.0;
	for (i = 0; i < dim; ++i)
	{
	    area[i] = pos_side ? -Tri_normal(t)[i] : Tri_normal(t)[i];
	    posn[i] /= 3.0;
	}
	/* Need to treat subdomain boundary */
	return YES;
}	/* end force_on_hse3d */

static double intrp_between(
	double x1,
	double x2,
	double x,
	double y1,
	double y2)
{
	double y;
	if (x1 == x2) return y1;
	y = y1 + (y2 - y1)/(x2 - x1)*(x - x1);
	return y;
}

static void get_variable_bdry_params(
	int dim,
	FILE *infile,
	POINTER *func_params)
{
	static VAR_BDRY_PARAMS params;
	int nbr_pistons;
	double half_angular_width;
	int i;

	params.dim=dim;
	CursorAfterString(infile,"Enter the center of the circle:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&params.center[i]);
            (void) printf("%f ",params.center[i]);
	}
	(void) printf("\n");

	CursorAfterString(infile,"Enter number of pistons:");
	fscanf(infile,"%d",&nbr_pistons);
	(void) printf("%d\n",nbr_pistons);

	params.number_pistons = nbr_pistons;

	FT_VectorMemoryAlloc((POINTER*)&params.angles_pistons,nbr_pistons+1,
				sizeof(double));
	for (i = 0; i < nbr_pistons+1; ++i)
		params.angles_pistons[i] = 2*PI*i/nbr_pistons;
	params.number_pistons = nbr_pistons + 1;

	CursorAfterString(infile,"Enter angular width of pistons:");
	fscanf(infile,"%lf",&half_angular_width);
	(void) printf("%f\n",half_angular_width);
	params.half_angular_width = half_angular_width;

	CursorAfterString(infile,
		"Enter radial velocity, density and pressure at piston:");
	fscanf(infile,"%lf %lf %lf",&params.bdry_vel,&params.bdry_dens,
					&params.bdry_pres);
	(void) printf("%f %f %f\n",params.bdry_vel,params.bdry_dens,
					params.bdry_pres);
	CursorAfterString(infile,
		"Enter time duration of the piston:");
	fscanf(infile,"%lf",&params.jet_duration_time);
	(void) printf("%f\n",params.jet_duration_time);

	*func_params = (POINTER)&params;
}	/* end get_variable_bdry_params */

extern void restart_set_dirichlet_bdry_function(Front *front)
{
        INTERFACE *intfc = front->interf;
        int i;
        BOUNDARY_STATE  *bstate;
        const char *s;
        for (i = 0; i < num_bstates(intfc); ++i)
        {
            bstate = bstate_list(intfc)[i];
            if (bstate == NULL) continue;
            s = bstate->_boundary_state_function_name;
            if (s == NULL) continue;
            if (strcmp(s,"cF_flowThroughBoundaryState") == 0)
                bstate->_boundary_state_function = cF_flowThroughBoundaryState;
	    else if (strcmp(s,"cF_variableBoundaryState") == 0)
                bstate->_boundary_state_function = cF_variableBoundaryState;
        }
}       /* end restart_set_dirichlet_bdry_function */
