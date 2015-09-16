/*******************************************************************
 * 		G_CARTESIAN.c
 *******************************************************************/
#include "cFluid.h"

class ToFill{
public:
int icoords[3];
};

EXPORT  void    tecplot_interface_states(const char*, INTERFACE	*);

static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};

//----------------------------------------------------------------
//		L_RECTANGLE
//----------------------------------------------------------------

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
// 		G_CARTESIAN
//--------------------------------------------------------------------------

G_CARTESIAN::~G_CARTESIAN()
{
}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup cell_center
//---------------------------------------------------------------
void G_CARTESIAN::initMesh(void)
{
	int i,j,k, index;
	double coords[2];
	int num_cells;

	// init cell_center
	L_RECTANGLE       rectangle;

	/*TMP*/
	min_dens = 0.0001;
	min_pres = 0.0001;
	FT_MakeGridIntfc(front);
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
	case 1:
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
		index = d_index1d(i,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    }
	    break;
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
	
	setComponent();
	FT_FreeGridIntfc(front);
}

void G_CARTESIAN::setComponent_old(void)
{
	int		i,j, ind;
	double 		*coords;
	int 		*icoords;
	COMPONENT 	old_comp,new_comp;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
        double          ***Gpdens = eqn_params->Gpdens;
	double		**Gpres = eqn_params->Gpres;
	static STATE 	*state = NULL;
	double		*dens = field.dens;
        double          **pdens = field.pdens;
	double		*engy = field.engy;
	double		*pres = field.pres;
	double		**momn = field.momn;
	int		size = (int)cell_center.size();
	
	// cell center components
	if(state == NULL)
	    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));

	for (i = 0; i < size; i++)
	{
	    icoords = cell_center[i].icoords;
	    coords = cell_center[i].m_coords;
	    old_comp = cell_center[i].comp;
	    new_comp = top_comp[i];
	    if (eqn_params->tracked && cell_center[i].comp != -1 &&
		cell_center[i].comp != top_comp[i] && gas_comp(new_comp))
	    {
		if (!FrontNearestIntfcState(front,coords,new_comp,
				(POINTER)state))
		{
		    (void) printf("In setComponent()\n");
		    (void) printf("FrontNearestIntfcState() failed\n");
		    (void) printf("old_comp = %d new_comp = %d\n",
					old_comp,new_comp);
		    clean_up(ERROR);
		}

		//GFM
		state->dim = dim;
		state->eos = &eqn_params->eos[new_comp];
		if (gas_comp(old_comp) && gas_comp(new_comp))
		{
		    if(new_comp == GAS_COMP1)
			ind = 0;
		    else
			ind = 1;

		    state->dens = Gdens[ind][i];
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            state->pdens[ii] = Gpdens[ind][ii][i];
                        }
                    }
		    state->pres = Gpres[ind][i];
		    for(j = 0; j < dim; ++j)
			state->momn[j] = Gvel[ind][j][i]*Gdens[ind][i];
		    state->engy = EosEnergy(state);
		}

		dens[i] = state->dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        pdens[ii][i] = state->pdens[ii];
                    }
                }
		pres[i] = state->pres;
		engy[i] = state->engy;
		for (j = 0; j < dim; ++j)
		    momn[j][i] = state->momn[j];
	    }
	    cell_center[i].comp = top_comp[i];
	}
}	/* end setComponent_old() */

void G_CARTESIAN::setComponent(void)
{
        int             i,j, ind;
        double          *coords;
        int             *icoords;
        COMPONENT       old_comp,new_comp;
        double          ***Gvel = eqn_params->Gvel;
        double          **Gdens = eqn_params->Gdens;
        double          **Gpres = eqn_params->Gpres;
        static STATE    *state = NULL;
        double          *dens = field.dens;
        double          **pdens = field.pdens;
        double          *engy = field.engy;
        double          *pres = field.pres;
        double          **momn = field.momn;
	double		*gamma = field.gamma;
        int             size = (int)cell_center.size();

        // cell center components
        if(state == NULL)
            FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));

        for (i = 0; i < size; i++)
        {
            icoords = cell_center[i].icoords;
            coords = cell_center[i].m_coords;
            old_comp = cell_center[i].comp;
            new_comp = top_comp[i];

            if (eqn_params->tracked && cell_center[i].comp != -1 &&
                cell_center[i].comp != top_comp[i] && gas_comp(new_comp))
            {
                if (!FrontNearestIntfcState(front,coords,new_comp,
                                (POINTER)state))
                {
                    (void) printf("In setComponent()\n");
                    (void) printf("FrontNearestIntfcState() failed\n");
                    (void) printf("old_comp = %d new_comp = %d\n",
                                        old_comp,new_comp);
                    clean_up(ERROR);
                }

                state->dim = dim;
                state->eos = &eqn_params->eos[new_comp];

                if (gas_comp(old_comp) && gas_comp(new_comp))
                {
                    if(new_comp == GAS_COMP1)
                        ind = 0;
                    else
                        ind = 1;
                    setRPGhost(i,state,ind);
                }

                dens[i] = state->dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        pdens[ii][i] = state->pdens[ii];
                    }
                }
                pres[i] = state->pres;
                engy[i] = state->engy;
                for (j = 0; j < dim; ++j)
                    momn[j][i] = state->momn[j];
//		gamma[i] = EosGamma(state);	//Dan
            }

            cell_center[i].comp = top_comp[i];
        }
}       /* end setComponent() */

void G_CARTESIAN::setRPGhost(
        int     index,
        STATE   *state,
        int     ind)
{
        int             i, dim = state->dim;
        double          **vel = eqn_params->vel;
        double          *dens = eqn_params->dens;
        double          **pdens = eqn_params->pdens;
        double          *pres = eqn_params->pres;
        EOS_PARAMS      *eos = eqn_params->eos;
        double          ***Gvel = eqn_params->Gvel;
        double          **Gdens = eqn_params->Gdens;
        double          ***Gpdens = eqn_params->Gpdens;
        double          **Gpres = eqn_params->Gpres;
        double          **Gnor = eqn_params->gnor;
        double          nor[MAXD];
        double          lvnor, rvnor;
        double          lvtang[MAXD], rvtang[MAXD];

        double pml, pmr, uml, umr, ml, mr;
        RIEMANN_SOLVER_WAVE_TYPE l_wave,r_wave;
        STATE *stl, *str;
        STATE *ansl;

        FT_ScalarMemoryAlloc((POINTER*)&stl,sizeof(STATE));
        FT_ScalarMemoryAlloc((POINTER*)&str,sizeof(STATE));
        FT_ScalarMemoryAlloc((POINTER*)&ansl,sizeof(STATE));

        if (ind == 0)
        {
            stl->eos = &(eos[GAS_COMP1]);
            str->eos = &(eos[GAS_COMP2]);
        }
        else
        {
            stl->eos = &(eos[GAS_COMP2]);
            str->eos = &(eos[GAS_COMP1]);
        }
        stl->dim = dim;
        str->dim = dim;
        stl->dens = Gdens[ind][index];
        stl->pdens[0] = Gpdens[ind][0][index];
        stl->pdens[1] = Gpdens[ind][1][index];
        str->dens = dens[index];
        str->pdens[0] = pdens[0][index];
        str->pdens[1] = pdens[1][index];
        stl->pres = Gpres[ind][index];
        str->pres = pres[index];
        lvnor = 0;
        rvnor = 0;
        for (i=0; i<dim; i++)
                rvnor += vel[i][index]*Gnor[i][index];
        if (rvnor < 0)
        {
                rvnor = -rvnor;
                for (i=0; i<dim; i++)
                        nor[i] = -Gnor[i][index];
        }
        else
                for (i=0; i<dim; i++)
                        nor[i] = Gnor[i][index];
        for (i=0; i<dim; i++)
                lvnor += Gvel[ind][i][index]*nor[i];
        for (i=0; i<dim; i++)
        {
                lvtang[i] = Gvel[ind][i][index] - lvnor*nor[i];
                rvtang[i] = vel[i][index] - rvnor*nor[i];
        }
        stl->vel[0] = lvnor;
        str->vel[0] = rvnor;
        for (i = 1; i < dim; i++)
            stl->vel[i] = str->vel[i] = 0.0;

        find_mid_state(stl,str,0.0/*pjump*/,&pml,&pmr,&uml,&umr,&ml,&mr,&l_wave,&r_wave);
        midstate(stl,ansl,ml,uml,pml,l_wave,1);

        state->pres = ansl->pres;
        state->dens = ansl->dens;
	if(eqn_params->multi_comp_non_reactive == YES)
	{
	    int ii;
	    for(ii = 0; ii < eqn_params->n_comps; ii++)
		state->pdens[ii] = ansl->pdens[ii];
	}
        for (i=0; i<dim; i++)
        {
                state->vel[i] = ansl->vel[0]*nor[i] + lvtang[i];
                state->momn[i] = state->dens*state->vel[i];
        }
        state->engy = EosEnergy(state);

        FT_FreeThese(2,stl,str);
        FT_FreeThese(1,ansl);
}

void G_CARTESIAN::setInitialIntfc(
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	dim = front->rect_grid->dim;
	eqn_params = (EQN_PARAMS*)front->extra1;
	switch (eqn_params->prob_type)
	{
	case TWO_FLUID_RT:
	case TWO_FLUID_RM:
	    initSinePertIntfc(level_func_pack,inname);
	    break;
	case TWO_FLUID_VST_RM:
//	    initSinePertIntfc(level_func_pack,inname);
	    initVSTRMIntfc(level_func_pack,inname);
	    break;
	case TWO_FLUID_BUBBLE:
	    initCirclePlaneIntfc(level_func_pack,inname);
	    break;
	case IMPLOSION:
	    initImplosionIntfc(level_func_pack,inname);
	    break;
	case MT_FUSION:
	    initMTFusionIntfc(level_func_pack,inname);
	    break;
	case PROJECTILE:
	    initProjectileIntfc(level_func_pack,inname);
	    break;
	case RIEMANN_PROB:
	case ONED_BLAST:
	case ONED_SSINE:
	    initRiemannProb(level_func_pack,inname);
	    break;
	default:
	    (void) printf("Problem type not implemented, code needed!\n");
	    clean_up(ERROR);
	}
}	/* end setInitialIntfc */

void G_CARTESIAN::setProbParams(char *inname)
{
	dim = front->rect_grid->dim;
	eqn_params = (EQN_PARAMS*)front->extra1;
	switch (eqn_params->prob_type)
	{
	case TWO_FLUID_RT:
	    setRayleiTaylorParams(inname);
	    break;
	case TWO_FLUID_RM:
	    setRichtmyerMeshkovParams(inname);
	    break;
	case TWO_FLUID_VST_RM:
//	    setRichtmyerMeshkovParams(inname);
	    setVSTRMParams(inname);
	    break;
	case TWO_FLUID_BUBBLE:
	    setBubbleParams(inname);
	    break;
	case IMPLOSION:
	    setImplosionParams(inname);
	    break;
	case MT_FUSION:
	    setMTFusionParams(inname);
	    break;
	case PROJECTILE:
	    setProjectileParams(inname);
	    break;
	case RIEMANN_PROB:
	    setRiemProbParams(inname);
	    break;
	case ONED_BLAST:
	case ONED_SSINE:
	    setOnedParams(inname);
	    break;
	default:
	    printf("In setProbParams(), unknown problem type!\n");
	    clean_up(ERROR);
	}
}	/* end setProbParams */

void G_CARTESIAN::setInitialStates()
{
	switch (eqn_params->prob_type)
	{
	case TWO_FLUID_RT:
	    initRayleiTaylorStates();
	    break;
	case TWO_FLUID_RM:
	    initRichtmyerMeshkovStates();
	    break;
	case TWO_FLUID_VST_RM:
//	    initRichtmyerMeshkovStates();
	    initVSTRMStates();
	    break;
	case TWO_FLUID_BUBBLE:
	    initBubbleStates();
	    break;
	case IMPLOSION:
	    initImplosionStates();
	    break;
	case MT_FUSION:
	    initMTFusionStates();
	    break;
	case PROJECTILE:
	    initProjectileStates();
	    break;
	case RIEMANN_PROB:
	    initRiemProbStates();
	    break;
	case ONED_BLAST:
	    initBlastWaveStates();
	    break;
	case ONED_SSINE:
	    initShockSineWaveStates();
	    break;
	default:
	    (void) printf("In setInitialStates(), case not implemented!\n");
	    clean_up(ERROR);
	}
	copyMeshStates();
}	/* end setInitialStates */

void G_CARTESIAN::computeAdvection(void)
{
	int order;
	switch (eqn_params->num_scheme)
	{
	case TVD_FIRST_ORDER:
	case WENO_FIRST_ORDER:
	    nrad = 3;
	    order = 1;
	    break;
	case TVD_SECOND_ORDER:
	case WENO_SECOND_ORDER:
	    nrad = 3;
	    order = 2;
	    break;
	case TVD_FOURTH_ORDER:
	case WENO_FOURTH_ORDER:
	    nrad = 3;
	    order = 4;
	    break;
	default:
	    order = -1;
	}
	solveRungeKutta(order);
}	/* end computeAdvection */


void G_CARTESIAN::solveRungeKutta(int order)
{
	static SWEEP *st_field,st_tmp;
	static FSWEEP *st_flux;
	static double **a,*b;
	double delta_t;
	int i,j;

	/* Allocate memory for Runge-Kutta of order */
	start_clock("solveRungeKutta");
	if (st_flux == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&b,order,sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&a,order,order,sizeof(double));

	    FT_VectorMemoryAlloc((POINTER*)&st_field,order,sizeof(SWEEP));
	    FT_VectorMemoryAlloc((POINTER*)&st_flux,order,sizeof(FSWEEP));
	    for (i = 0; i < order; ++i)
	    {
	    	allocMeshVst(&st_tmp);
	    	allocMeshVst(&st_field[i]);
	    	allocMeshFlux(&st_flux[i]);
	    }
	    /* Set coefficient a, b, c for different order of RK method */
	    switch (order)
	    {
	    case 1:
		b[0] = 1.0;
	    	break;
	    case 2:
	    	a[0][0] = 1.0;
	    	b[0] = 0.5;  b[1] = 0.5;
	    	break;
	    case 4:
	    	a[0][0] = 0.5;
	    	a[1][0] = 0.0;  a[1][1] = 0.5;
	    	a[2][0] = 0.0;  a[2][1] = 0.0;  a[2][2] = 1.0;
	    	b[0] = 1.0/6.0;  b[1] = 1.0/3.0;
	    	b[2] = 1.0/3.0;  b[3] = 1.0/6.0;
	    	break;
	    default:
	    	(void)printf("ERROR: %d-th order RK method not implemented\n",
					order);
	    	clean_up(ERROR);
	    }
	}
	delta_t = m_dt;

	/* Compute flux and advance field */

	copyToMeshVst(&st_field[0]);
	computeMeshFlux(st_field[0],&st_flux[0],delta_t);
	
	for (i = 0; i < order-1; ++i)
	{
	    copyMeshVst(st_field[0],&st_field[i+1]);
	    for (j = 0; j <= i; ++j)
	    {
		if (a[i][j] != 0.0)
		{
		    addMeshFluxToVst(&st_field[i+1],st_flux[j],a[i][j]);
		}
	    }
	    computeMeshFlux(st_field[i+1],&st_flux[i+1],delta_t);
	}
	for (i = 0; i < order; ++i)
	{
	    if (b[i] != 0.0)
	    {
		addMeshFluxToVst(&st_field[0],st_flux[i],b[i]);
	    }
	}
	copyFromMeshVst(st_field[0]);
	stop_clock("solveRungeKutta");
}	/* end solveRungeKutta */

void G_CARTESIAN::computeMeshFlux(
	SWEEP m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int dir;

	if(eqn_params->tracked)
	{
	    get_ghost_state(m_vst, 2, 0);
	    get_ghost_state(m_vst, 3, 1);
	    scatMeshGhost();
	    solve_exp_value();
	}

	resetFlux(m_flux);
	for (dir = 0; dir < dim; ++dir)
	{
	    addFluxInDirection(dir,&m_vst,m_flux,delta_t);
	}
	addSourceTerm(&m_vst,m_flux,delta_t);
}	/* end computeMeshFlux */

void G_CARTESIAN::resetFlux(FSWEEP *m_flux)
{
	int i,j;
	int size = (int)cell_center.size();
	for (i = 0; i < size; i++)
	{
	    m_flux->dens_flux[i] = 0.0;
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                int ii;
                for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    m_flux->pdens_flux[ii][i] = 0.0;
                }
            }
	    m_flux->engy_flux[i] = 0.0;
	    for (j = 0; j < MAXD; ++j)
	    	m_flux->momn_flux[j][i] = 0.0;
	}
}	/* resetFlux */

void G_CARTESIAN::addFluxInDirection(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	switch (dim)
	{
	case 1:
	    return addFluxInDirection1d(dir,m_vst,m_flux,delta_t);
	case 2:
	    return addFluxInDirection2d(dir,m_vst,m_flux,delta_t);
	case 3:
	    return addFluxInDirection3d(dir,m_vst,m_flux,delta_t);
	}
}	/* end addFluxInDirection */

void G_CARTESIAN::addFluxInDirection1d(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int i,n,index;
	SCHEME_PARAMS scheme_params;
	EOS_PARAMS	*eos;
	static SWEEP vst;
	static FSWEEP vflux;
	static boolean first = YES;
	COMPONENT comp;
	int seg_min,seg_max;
	static int icoords[MAXD];
	STATE st;
	
	start_clock("addFluxInDirection1d");
	if (first)
	{
	    first = NO;
	    allocDirVstFlux(&vst,&vflux);
	}

	scheme_params.lambda = delta_t/top_h[dir];
	scheme_params.beta = 0.0;

	seg_min = imin[0];
	while (seg_min <= imax[0])
	{
	    for (; seg_min <= imax[0]; ++seg_min)
	    {
		i = seg_min;
	    	index = d_index1d(i,top_gmax);
	    	comp = top_comp[index];
	    	if (gas_comp(comp)) break;
	    }
	    if (seg_min > imax[0]) break;
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
	    	vst.dens[i] = 0.0; 
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        vst.pdens[ii][i] = 0.0;
                    }
                }
	    	vst.pres[i] = 0.0; 
	    	vst.engy[i] = 0.0;
		vst.gamma[i] = 0.0;
	    	vst.momn[0][i] = vst.momn[1][i] = vst.momn[2][i] = 0.0;
	    }
	    i = seg_min;
	    index = d_index1d(i,top_gmax);
	    comp = top_comp[index];
	    n = 0;
	    vst.dens[n+nrad] = m_vst->dens[index];
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                int ii;
                for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                }
            }
            vst.engy[n+nrad] = m_vst->engy[index];
            vst.pres[n+nrad] = m_vst->pres[index];
            vst.momn[0][n+nrad] = m_vst->momn[0][index];
            vst.momn[1][n+nrad] = 0.0;
            vst.momn[2][n+nrad] = 0.0;
	    seg_max = i;
	    n++;
	    for (i = seg_min+1; i <= imax[0]; i++)
	    {
		index = d_index1d(i,top_gmax);
		if (needBufferFromIntfc(comp,top_comp[index]))
		    break;
		else
		{
	    	    vst.dens[n+nrad] = m_vst->dens[index];
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                        }
                    }
	    	    vst.engy[n+nrad] = m_vst->engy[index];
	    	    vst.pres[n+nrad] = m_vst->pres[index];
	    	    vst.momn[0][n+nrad] = m_vst->momn[0][index];
	    	    vst.momn[1][n+nrad] = 0.0;
	    	    vst.momn[2][n+nrad] = 0.0;
		    n++;
		}
		seg_max = i;
	    }
	    icoords[0] = seg_min;
	    appendGhostBuffer(&vst,m_vst,n,icoords,0,0);
	    icoords[0] = seg_max;
	    appendGhostBuffer(&vst,m_vst,n,icoords,0,1);

	    //Dan	FIXME
	    for (i = 0; i <= n+2*nrad-1; i++)
	    {
		st.eos = &(eqn_params->eos[comp]);
		st.dens = vst.dens[i];
		if(eqn_params->multi_comp_non_reactive == YES)
		{
		    int ii;
		    for(ii = 0; ii < eqn_params->n_comps; ii++)
		    {
			st.pdens[ii] = vst.pdens[ii][n+nrad];
		    }
		}
		vst.gamma[i] = EosGamma(&st);
	    }
	    //Dan	FIXME
	    
	    eos = &(eqn_params->eos[comp]);
	    EosSetTVDParams(&scheme_params, eos);
	    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
	    
	    n = 0;
	    for (i = seg_min; i <= seg_max; ++i)
	    {
	    	index = d_index1d(i,top_gmax);
	    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_flux->pdens_flux[ii][index] += vflux.pdens_flux[ii][n+nrad];
                    }
                }
	    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
	    	m_flux->momn_flux[0][index] += 
				vflux.momn_flux[0][n+nrad];
	    	m_flux->momn_flux[1][index] = 0.0;
	    	m_flux->momn_flux[2][index] = 0.0;
		n++;
	    }
	    seg_min = seg_max + 1;
	}
	stop_clock("addFluxInDirection1d");
}	/* end addFluxInDirection1d */

void G_CARTESIAN::addFluxInDirection2d(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int i,j,n,index;
	SCHEME_PARAMS scheme_params;
	EOS_PARAMS	*eos;
	static SWEEP vst;
	static FSWEEP vflux;
	static boolean first = YES;
	COMPONENT comp;
	int seg_min,seg_max;
	static int icoords[MAXD];
	STATE st;
	
	start_clock("addFluxInDirection2d");
	if (first)
	{
	    first = NO;
	    allocDirVstFlux(&vst,&vflux);
	}

	scheme_params.lambda = delta_t/top_h[dir];
	scheme_params.beta = 0.0;
	switch (dir)
	{
	case 0:
	    for (j = imin[1]; j <= imax[1]; j++)
	    {
		seg_min = imin[0];
		while (seg_min <= imax[0])
		{
		    for (; seg_min <= imax[0]; ++seg_min)
		    {
			i = seg_min;
		    	index = d_index2d(i,j,top_gmax);
		    	comp = top_comp[index];
		    	if (gas_comp(comp)) break;
		    }
		    if (seg_min > imax[0]) break;
		    for (i = 0; i <= top_gmax[0]; ++i)
		    {
		    	vst.dens[i] = 0.0; 
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;
                            for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {
                                vst.pdens[ii][i] = 0.0;
                            }
                        }
		    	vst.pres[i] = 0.0; 
		    	vst.engy[i] = 0.0; 
			vst.gamma[i] = 0.0;
		    	vst.momn[0][i] = vst.momn[1][i] = vst.momn[2][i] = 0.0;
		    }
		    i = seg_min;
		    index = d_index2d(i,j,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                        }
                    }
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[0][index];
                    vst.momn[1][n+nrad] = m_vst->momn[1][index];
                    vst.momn[2][n+nrad] = 0.0;
		    seg_max = i;
		    n++;
		    for (i = seg_min+1; i <= imax[0]; i++)
		    {
			index = d_index2d(i,j,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
			    break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
                            if(eqn_params->multi_comp_non_reactive == YES)
                            {
                                int ii;
                                for(ii = 0; ii < eqn_params->n_comps; ii++)
                                {
                                    vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                                }
                            }
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[0][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[1][index];
		    	    vst.momn[2][n+nrad] = 0.0;
			    n++;
			}
			seg_max = i;
		    }
		    icoords[1] = j;
		    icoords[0] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,0);
		    icoords[0] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,1);

		    //Dan	FIXME
		    for (i = 0; i <= n+2*nrad-1; i++)
		    {
			st.eos = &(eqn_params->eos[comp]);
			st.dens = vst.dens[i];
			if(eqn_params->multi_comp_non_reactive == YES)
			{
			    int ii;
			    for(ii = 0; ii < eqn_params->n_comps; ii++)
			    {
				st.pdens[ii] = vst.pdens[ii][n+nrad];
			    }
			}
			vst.gamma[i] = EosGamma(&st);
		    }
		    //Dan	FIXME

		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (i = seg_min; i <= seg_max; ++i)
		    {
		    	index = d_index2d(i,j,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;                            
                            for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {                                
                                m_flux->pdens_flux[ii][index] += vflux.pdens_flux[ii][n+nrad];
                            }
                        }
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[1][n+nrad];
		    	m_flux->momn_flux[2][index] = 0.0;
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	case 1:
	    for (i = imin[0]; i <= imax[0]; i++)
	    {
		seg_min = imin[1];
		while (seg_min <= imax[1])
		{
		    for (; seg_min <= imax[1]; ++seg_min)
		    {
			j = seg_min;
		    	index = d_index2d(i,j,top_gmax);
		    	comp = top_comp[index];
		    	if (gas_comp(comp)) break;
		    }
		    if (seg_min > imax[1]) break;
		    for (j = 0; j <= top_gmax[1]; ++j)
		    {
		    	vst.dens[j] = 0.0; 
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;
			    for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {
                                vst.pdens[ii][i] = 0.0;
                            }
                        }
		    	vst.pres[j] = 0.0; 
		    	vst.engy[j] = 0.0; 
			vst.gamma[i] = 0.0;
		    	vst.momn[0][j] = vst.momn[1][j] = vst.momn[2][j] = 0.0;
		    }
		    j = seg_min;
		    index = d_index2d(i,j,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    if(eqn_params->multi_comp_non_reactive == YES)
                    { 
                       int ii;           
                       for(ii = 0; ii < eqn_params->n_comps; ii++)
                       {
                           vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                       }
                    }
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[1][index];
                    vst.momn[1][n+nrad] = m_vst->momn[0][index];
                    vst.momn[2][n+nrad] = 0.0;
		    seg_max = j;
		    n++;
		    for (j = seg_min+1; j <= imax[1]; j++)
		    {
			index = d_index2d(i,j,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
			    break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
                            if(eqn_params->multi_comp_non_reactive == YES)
                            {
                                int ii;
                                for(ii = 0; ii < eqn_params->n_comps; ii++)
                                {
                                    vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                                }
                            }
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[1][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[0][index];
		    	    vst.momn[2][n+nrad] = 0.0;
			    n++;
			}
			seg_max = j;
		    }
		    icoords[0] = i;
		    icoords[1] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,0);
		    icoords[1] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,1);
		    
		    //Dan	FIXME
		    for (j = 0; j <= n+2*nrad-1; j++)
		    {
			st.eos = &(eqn_params->eos[comp]);
			st.dens = vst.dens[j];
			if(eqn_params->multi_comp_non_reactive == YES)
			{
			    int ii;
			    for(ii = 0; ii < eqn_params->n_comps; ii++)
			    {
				st.pdens[ii] = vst.pdens[ii][n+nrad];
			    }
			}
			vst.gamma[j] = EosGamma(&st);
		    }
		    //Dan	FIXME

		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (j = seg_min; j <= seg_max; ++j)
		    {
		    	index = d_index2d(i,j,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;
                            for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {
                                m_flux->pdens_flux[ii][index] += vflux.pdens_flux[ii][n+nrad];
                            }
                        }
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[1][n+nrad];
		    	m_flux->momn_flux[2][index] = 0.0;
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	}
	stop_clock("addFluxInDirection2d");
}	/* end addFluxInDirection2d */

void G_CARTESIAN::addFluxInDirection3d(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int		i,j,k,n,index;
	SCHEME_PARAMS	scheme_params;
	EOS_PARAMS	*eos;
	static SWEEP 	vst;
	static FSWEEP 	vflux;
	static boolean 	first = YES;
	COMPONENT 	comp;
	int 		seg_min,seg_max;
	int 		icoords[3];
	STATE st;
	
	start_clock("addFluxInDirection3d");
	if (first)
	{
	    first = NO;
	    allocDirVstFlux(&vst,&vflux);
	}
	
	scheme_params.lambda = delta_t/top_h[dir];
	scheme_params.beta = 0.0;
	
	switch (dir)
	{
	case 0:
	    for (k = imin[2]; k <= imax[2]; k++)
	    for (j = imin[1]; j <= imax[1]; j++)
	    {
		seg_min = imin[0];
		while (seg_min <= imax[0])
		{
		    for (; seg_min <= imax[0]; ++seg_min)
                    {
                        i = seg_min;
                        index = d_index3d(i,j,k,top_gmax);
                        comp = top_comp[index];
                        if (gas_comp(comp)) break;
                    }
                    if (seg_min > imax[0]) break;
		    for (i = 0; i <= top_gmax[1]; ++i)
		    {
		    	vst.dens[i] = 0.0; 
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;
                            for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {
                                vst.pdens[ii][i] = 0.0;
                            }
                        }
		    	vst.pres[i] = 0.0; 
		    	vst.engy[i] = 0.0; 
			vst.gamma[i] = 0.0;
		    	vst.momn[0][i] = vst.momn[1][i] = vst.momn[2][i] = 0.0;
		    }
		    i = seg_min;
		    index = d_index3d(i,j,k,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                        }
                    }
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[0][index];
                    vst.momn[1][n+nrad] = m_vst->momn[1][index];
                    vst.momn[2][n+nrad] = m_vst->momn[2][index];
		    seg_max = i;
		    n++;
		    for (i = seg_min+1; i <= imax[0]; i++)
		    {
			index = d_index3d(i,j,k,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
                            break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
                            if(eqn_params->multi_comp_non_reactive == YES)
                            {
                                int ii;
                                for(ii = 0; ii < eqn_params->n_comps; ii++)
                                {
                                    vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                                }
                            }
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[0][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[1][index];
		    	    vst.momn[2][n+nrad] = m_vst->momn[2][index];
			    n++;
			}
			seg_max = i;
		    }
		    
		    icoords[1] = j;
		    icoords[2] = k;
		    icoords[0] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,0);
		    icoords[0] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,1);
		    
		    //Dan	FIXME
		    for (i = 0; i <= n+2*nrad-1; i++)
		    {
			st.eos = &(eqn_params->eos[comp]);
			st.dens = vst.dens[i];
			if(eqn_params->multi_comp_non_reactive == YES)
			{
			    int ii;
			    for(ii = 0; ii < eqn_params->n_comps; ii++)
			    {
				st.pdens[ii] = vst.pdens[ii][n+nrad];
			    }
			}
			vst.gamma[i] = EosGamma(&st);
		    }
		    //Dan	FIXME

		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (i = seg_min; i <= seg_max; ++i)
		    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;
                            for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {
                                m_flux->pdens_flux[ii][index] += vflux.pdens_flux[ii][n+nrad];
                            }
                        }
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[1][n+nrad];
		    	m_flux->momn_flux[2][index] +=
					vflux.momn_flux[2][n+nrad];
			n++;
		    }

		    seg_min = seg_max + 1;
		}
	    }
	    break;
	case 1:
	    for (k = imin[2]; k <= imax[2]; k++)
	    for (i = imin[0]; i <= imax[0]; i++)
	    {
		seg_min = imin[1];
		while (seg_min <= imax[1])
		{
		    for (; seg_min <= imax[1]; ++seg_min)
                    {
                        j = seg_min;
                        index = d_index3d(i,j,k,top_gmax);
                        comp = top_comp[index];
                        if (gas_comp(comp)) break;
                    }
                    if (seg_min > imax[1]) break;
		    for (j = 0; j <= top_gmax[1]; ++j)
		    {
		    	vst.dens[j] = 0.0; 
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;
                            for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {
                                vst.pdens[ii][j] = 0.0;
                            }
                        }
		    	vst.pres[j] = 0.0; 
		    	vst.engy[j] = 0.0; 
			vst.gamma[i] = 0.0;
		    	vst.momn[0][j] = vst.momn[1][j] = vst.momn[2][j] = 0.0;
		    }
		    j = seg_min;
		    index = d_index3d(i,j,k,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                        }
                    }
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[1][index];
                    vst.momn[1][n+nrad] = m_vst->momn[2][index];
                    vst.momn[2][n+nrad] = m_vst->momn[0][index];
		    seg_max = j;
		    n++;
		    
		    for (j = seg_min+1; j <= imax[1]; j++)
		    {
			index = d_index3d(i,j,k,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
			    break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
                            if(eqn_params->multi_comp_non_reactive == YES)
                            {
                                int ii;
                                for(ii = 0; ii < eqn_params->n_comps; ii++)
                                {
                                    vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                                }
                            }
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[1][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[2][index];
		    	    vst.momn[2][n+nrad] = m_vst->momn[0][index];
			    n++;
			}
			seg_max = j;
		    }
		    icoords[0] = i;
		    icoords[2] = k;
		    icoords[1] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,0);
		    icoords[1] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,1);
		    
		    //Dan	FIXME
		    for (j = 0; j <= n+2*nrad-1; j++)
		    {
			st.eos = &(eqn_params->eos[comp]);
			st.dens = vst.dens[j];
			if(eqn_params->multi_comp_non_reactive == YES)
			{
			    int ii;
			    for(ii = 0; ii < eqn_params->n_comps; ii++)
			    {
				st.pdens[ii] = vst.pdens[ii][n+nrad];
			    }
			}
			vst.gamma[j] = EosGamma(&st);
		    }
		    //Dan	FIXME

		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (j = seg_min; j <= seg_max; ++j)
		    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;
                            for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {
                                m_flux->pdens_flux[ii][index] += vflux.pdens_flux[ii][n+nrad];
                            }
                        }
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[2][n+nrad];
		    	m_flux->momn_flux[2][index] += 
					vflux.momn_flux[1][n+nrad];
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; j++)
	    for (i = imin[0]; i <= imax[0]; i++)
	    {
		seg_min = imin[2];
		while (seg_min <= imax[2])
		{
		    for (; seg_min <= imax[2]; ++seg_min)
                    {
                        k = seg_min;
                        index = d_index3d(i,j,k,top_gmax);
                        comp = top_comp[index];
                        if (gas_comp(comp)) break;
                    }
                    if (seg_min > imax[2]) break;
		    for (k = 0; k <= top_gmax[2]; ++k)
		    {
		    	vst.dens[k] = 0.0; 
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;
                            for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {
                                vst.pdens[ii][k] = 0.0;
                            }
                        }
		    	vst.pres[k] = 0.0; 
		    	vst.engy[k] = 0.0; 
			vst.gamma[i] = 0.0;
		    	vst.momn[0][k] = vst.momn[1][k] = vst.momn[2][k] = 0.0;
		    }
		    k = seg_min;
		    index = d_index3d(i,j,k,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                        }
                    }
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[2][index];
                    vst.momn[1][n+nrad] = m_vst->momn[0][index];
                    vst.momn[2][n+nrad] = m_vst->momn[1][index];
		    seg_max = k;
		    n++;
		    
		    for (k = seg_min+1; k <= imax[2]; k++)
		    {
			index = d_index3d(i,j,k,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
			    break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
                            if(eqn_params->multi_comp_non_reactive == YES)
                            {
                                int ii;
                                for(ii = 0; ii < eqn_params->n_comps; ii++)
                                {
                                    vst.pdens[ii][n+nrad] = m_vst->pdens[ii][index];
                                }
                            }
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[2][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[0][index];
		    	    vst.momn[2][n+nrad] = m_vst->momn[1][index];
			    n++;
			}
			seg_max = k;
		    }
		    icoords[0] = i;
		    icoords[1] = j;
		    icoords[2] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,2,0);
		    icoords[2] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,2,1);
		    
		    //Dan	FIXME
		    for (k = 0; k <= n+2*nrad-1; k++)
		    {
			st.eos = &(eqn_params->eos[comp]);
			st.dens = vst.dens[k];
			if(eqn_params->multi_comp_non_reactive == YES)
			{
			    int ii;
			    for(ii = 0; ii < eqn_params->n_comps; ii++)
			    {
				st.pdens[ii] = vst.pdens[ii][n+nrad];
			    }
			}
			vst.gamma[k] = EosGamma(&st);
		    }
		    //Dan	FIXME

		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (k = seg_min; k <= seg_max; ++k)
		    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
                        if(eqn_params->multi_comp_non_reactive == YES)
                        {
                            int ii;
                            for(ii = 0; ii < eqn_params->n_comps; ii++)
                            {
                                m_flux->pdens_flux[ii][index] += vflux.pdens_flux[ii][n+nrad];
                            }
                        }
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[2][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[1][n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[2][n+nrad];
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	}
	stop_clock("addFluxInDirection3d");
}

void G_CARTESIAN::scatMeshFlux(FSWEEP *m_flux)
{
	int i,j,k,l,index;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_flux->dens_flux[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_flux->dens_flux[index] = array[index];
	    }
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_flux->engy_flux[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_flux->engy_flux[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index1d(i,top_gmax);
		    array[index] = m_flux->momn_flux[l][index];
	    	}
	    	scatMeshArray();
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index1d(i,op_gmax);
		    m_flux->momn_flux[l][index] = array[index];
	    	}
	    }
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_flux->dens_flux[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_flux->dens_flux[index] = array[index];
	    }
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_flux->engy_flux[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_flux->engy_flux[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index2d(i,j,top_gmax);
		    array[index] = m_flux->momn_flux[l][index];
	    	}
	    	scatMeshArray();
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index2d(i,j,top_gmax);
		    m_flux->momn_flux[l][index] = array[index];
	    	}
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		array[index] = m_flux->dens_flux[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
		index = d_index3d(i,j,k,top_gmax);
		m_flux->dens_flux[index] = array[index];
	    }
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		array[index] = m_flux->engy_flux[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_flux->engy_flux[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (k = imin[2]; k <= imax[2]; ++k)
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index3d(i,j,k,top_gmax);
		    array[index] = m_flux->momn_flux[l][index];
	    	}
	    	scatMeshArray();
            	for (k = 0; k <= top_gmax[2]; k++)
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
		    index = d_index3d(i,j,k,top_gmax);
		    m_flux->momn_flux[l][index] = array[index];
	    	}
	    }
	}
}	/* end scatMeshFlux */

void G_CARTESIAN::addSourceTerm(
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int i,j,k,l,index;
	double *gravity = eqn_params->gravity;

	switch (dim)
	{
	case 1:
            for (i = imin[0]; i <= imax[0]; i++)
            {
		index = d_index1d(i,top_gmax);
		if (!gas_comp(top_comp[index]))
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] = 0.0; 
		    	m_flux->engy_flux[index] = 0.0; 
		    }
		}
		else
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] += 
				delta_t*gravity[l]*m_vst->dens[index];
		    	m_flux->engy_flux[index] += 
				delta_t*gravity[l]*m_vst->momn[l][index];
		    }
		}
	    }
	    break;
	case 2:
            for (j = imin[1]; j <= imax[1]; j++)
            for (i = imin[0]; i <= imax[0]; i++)
            {
		index = d_index2d(i,j,top_gmax);
		if (!gas_comp(top_comp[index]))
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] = 0.0; 
		    	m_flux->engy_flux[index] = 0.0; 
		    }
		}
		else
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] += 
				delta_t*gravity[l]*m_vst->dens[index];
		    	m_flux->engy_flux[index] += 
				delta_t*gravity[l]*m_vst->momn[l][index];
		    }
		}
	    }
	    break;
	case 3:
            for (k = imin[2]; k <= imax[2]; k++)
            for (j = imin[1]; j <= imax[1]; j++)
            for (i = imin[0]; i <= imax[0]; i++)
            {
		index = d_index3d(i,j,k,top_gmax);
		if (!gas_comp(top_comp[index]))
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] = 0.0; 
		    	m_flux->engy_flux[index] = 0.0; 
		    }
		}
		else
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] += 
				delta_t*gravity[l]*m_vst->dens[index];
		    	m_flux->engy_flux[index] += 
				delta_t*gravity[l]*m_vst->momn[l][index];
		    }
		}
	    }
	}
	
}	/* end addSourceTerm */

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void G_CARTESIAN::solve(double dt)
{
	m_dt = dt;
	max_speed = 0.0;

	start_clock("solve");
	setDomain();

	setComponent();
	
	if (debugging("trace"))
	    printf("Passed setComponent()\n");

	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	computeAdvection();
	if (debugging("trace"))
	    printf("max_speed after computeAdvection(): %20.14f\n",max_speed);
	stop_clock("computeAdvection");

	//turbulence boundary layer	FIXME
	//for now, it's only set for Kathy's simulation
	if (eqn_params->TBL)
	{
	    computeTBL();	//Dan
	}
	else
	    printf("TBL is not called.\n");

        /* parabolic step added by PRAO */
        if (eqn_params->parabolic_step == true)
        {
            computeParab();
        }

	if (debugging("sample_velocity"))
	{
	    sampleVelocity();
	}

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
	stop_clock("solve");

}	/* end solve */


// check http://en.wikipedia.org/wiki/Bilinear_interpolation
void G_CARTESIAN::getVelocity(double *p, double *U)
{
        double **vel = eqn_params->vel;

        FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[0],getStateXvel,&U[0],
					NULL);
        if (dim > 1)
            FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[1],getStateYvel,&U[1],
					NULL);
        if (dim > 2)
            FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[2],getStateZvel,&U[2],
					NULL);
}

void G_CARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void G_CARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int G_CARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void G_CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].m_coords[i];
}

void G_CARTESIAN::getRectangleCenter(
	int index0, 
	int index1, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = 0.5*(cell_center[index0].m_coords[i] +
	    		     cell_center[index1].m_coords[i]);
	}
}


double G_CARTESIAN::getDistance(double *c0, double *c1)
{
	return sqrt( (c0[0]-c1[0])*(c0[0]-c1[0])
		    +(c0[1]-c1[1])*(c0[1]-c1[1]) );
}


// input : p[]
// output: q[]

void G_CARTESIAN::getNearestInterfacePoint(
	double *p, 
	double *q)
{
	INTERFACE *intfc = front->interf;
	double t;
	HYPER_SURF_ELEMENT *phse;
	HYPER_SURF *phs;
	nearest_interface_point(p,getComponent(p),intfc,NO_BOUNDARIES,
				NULL,q,&t,&phse,&phs);
}

int G_CARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int G_CARTESIAN::getComponent(
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
	default:
	    return NO_COMP;
	}
}

int G_CARTESIAN::getComponent(
	int index)
{
	switch (dim)
	{
	case 2:
	    return top_comp[index];
	case 3:
	    return top_comp[index];
	default:
	    return NO_COMP;
	}
}

void G_CARTESIAN::save(char *filename)
{
	
	INTERFACE *intfc    = front->interf;
		
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in "
		       "SaveAsTecplot_rect_grid_and_interface().", filename);
		clean_up(ERROR);
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

G_CARTESIAN::G_CARTESIAN(Front &front):front(&front)
{
}

void G_CARTESIAN::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i,size;

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	eqn_params = (EQN_PARAMS*)front->extra1;
	
	if (first)
	{
	    first = NO;
	    dim = grid_intfc->dim;

	    hmin = HUGE;
	    size = 1;
	    
            for (i = 0; i < 3; ++i)
	    	top_gmax[i] = 0;

            for (i = 0; i < dim; ++i)
	    {
	    	lbuf[i] = front->rect_grid->lbuf[i];
	    	ubuf[i] = front->rect_grid->ubuf[i];
	    	top_gmax[i] = top_grid->gmax[i];
	    	top_L[i] = top_grid->L[i];
	    	top_U[i] = top_grid->U[i];
	    	top_h[i] = top_grid->h[i];

                if (hmin > top_h[i]) hmin = top_h[i];
	        size *= (top_gmax[i]+1);
	    	imin[i] = (lbuf[i] == 0) ? 1 : lbuf[i];
	    	imax[i] = (ubuf[i] == 0) ? top_gmax[i] - 1 : 
				top_gmax[i] - ubuf[i];
	    }

	    FT_VectorMemoryAlloc((POINTER*)&eqn_params->dens,size,
					sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&eqn_params->pdens,2,size,
                                        sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&eqn_params->pres,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&eqn_params->engy,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->vel,dim,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->mom,dim,size,
					sizeof(double));
	    //GFM
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->gnor,dim,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->Gdens,2,size,
					sizeof(double));
            FT_TriArrayMemoryAlloc((POINTER*)&eqn_params->Gpdens,2,2,size,
                                        sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->Gpres,2,size,
					sizeof(double));
	    FT_TriArrayMemoryAlloc((POINTER*)&eqn_params->Gvel,2,dim,size,
					sizeof(double));

            /* fields added by PRAO for SGS */
            FT_TriArrayMemoryAlloc((POINTER*)&field.tau,dim,dim,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field.qt,size,sizeof(double));
            
            /* end  of fields added by PRAO */ 

	    FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
	    if (dim == 2)
	    	FT_VectorMemoryAlloc((POINTER*)&eqn_params->vort,size,
					sizeof(double));
	    else if (dim == 3)
	    	FT_MatrixMemoryAlloc((POINTER*)&eqn_params->vort3d,dim,size,
					sizeof(double));
	    field.dens = eqn_params->dens;
            field.pdens = eqn_params->pdens;
	    field.engy = eqn_params->engy;
	    field.pres = eqn_params->pres;
	    field.momn = eqn_params->mom;
	    field.vel = eqn_params->vel;
	}
}

void G_CARTESIAN::allocMeshVst(
	SWEEP *vst)
{
	int i,size;

	size = 1;
        for (i = 0; i < dim; ++i)
	    size *= (top_gmax[i]+1);

	FT_VectorMemoryAlloc((POINTER*)&vst->dens,size,sizeof(double));
        FT_MatrixMemoryAlloc((POINTER*)&vst->pdens,2,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->engy,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->pres,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&vst->momn,MAXD,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->comp,size,sizeof(COMPONENT));
	FT_VectorMemoryAlloc((POINTER*)&vst->gamma,size,sizeof(double));
}	/* end allocMeshVstFlux */

void G_CARTESIAN::allocMeshFlux(
	FSWEEP *flux)
{
	int i,size;

	size = 1;
        for (i = 0; i < dim; ++i)
	    size *= (top_gmax[i]+1);

	FT_VectorMemoryAlloc((POINTER*)&flux->dens_flux,size,sizeof(double));
        FT_MatrixMemoryAlloc((POINTER*)&flux->pdens_flux,2,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux->engy_flux,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&flux->momn_flux,MAXD,size,sizeof(double));
}	/* end allocMeshVstFlux */

void G_CARTESIAN::allocDirVstFlux(
        SWEEP *vst,
        FSWEEP *flux)
{
	int i,size;

	size = 1;
        for (i = 0; i < dim; ++i)
	    if (size < top_gmax[i]+7) 
		size = top_gmax[i]+7;
	FT_VectorMemoryAlloc((POINTER*)&vst->dens,size,sizeof(double));
        FT_MatrixMemoryAlloc((POINTER*)&vst->pdens,2,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->engy,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->pres,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&vst->momn,MAXD,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->gamma,size,sizeof(double));

	FT_VectorMemoryAlloc((POINTER*)&flux->dens_flux,size,sizeof(double));
        FT_MatrixMemoryAlloc((POINTER*)&flux->pdens_flux,2,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux->engy_flux,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&flux->momn_flux,MAXD,size,sizeof(double));
}	/* end allocDirMeshVstFlux */

void G_CARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front);
}

void G_CARTESIAN::checkVst(SWEEP *vst)
{
	int i,j,index;
	for (j = imin[1]; j < imax[1]; j++)
	for (i = imin[0]; i < imax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    if (isnan(vst->dens[index]))
		printf("At %d %d: dens is nan\n",i,j);
	    if (vst->dens[index] < 0.0)
		printf("At %d %d: dens is negative\n",i,j);
	}
}

void G_CARTESIAN::checkFlux(FSWEEP *flux)
{
	int i,j,index;
	//for (j = imin[1]; j < imax[1]; j++)
	j = 140;
	for (i = imin[0]; i <= imax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    printf("%d %f  %f\n",i,flux->momn_flux[1][index],
				flux->engy_flux[index]);
	}
}

void G_CARTESIAN::printFrontInteriorStates(char *out_name)
{
	int i,j,k,l,index;
	char filename[100];
	FILE *outfile;
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double *dens = field.dens;
	double *engy = field.engy;
	double **momn = field.momn;
        double **pdens = field.pdens;

	sprintf(filename,"%s/state.ts%s",out_name,
			right_flush(front->step,7));
	if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
	sprintf(filename,"%s-gas",filename);
	outfile = fopen(filename,"w");

        /* Initialize states at the interface */
        fprintf(outfile,"Interface gas states:\n");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",getStateDens(sl),
				getStateDens(sr));
            fprintf(outfile,"%24.18g %24.18g\n",getStateEngy(sl),
				getStateEngy(sr));
	    for (i = 0; i < dim; ++i)
            	fprintf(outfile,"%24.18g %24.18g\n",getStateMom[i](sl),
				getStateMom[i](sr));
            if(eqn_params->multi_comp_non_reactive == YES)
            { 
                //int ii;
                //for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    fprintf(outfile,"%24.18g %24.18g\n",getStatePdens0(sl),
                                getStatePdens0(sr));
                    fprintf(outfile,"%24.18g %24.18g\n",getStatePdens1(sl),
                                getStatePdens1(sr));
                }
            }
        }
	
	fprintf(outfile,"\nInterior gas states:\n");
	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	        fprintf(outfile,"%24.18g\n",dens[index]);
	        fprintf(outfile,"%24.18g\n",engy[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",momn[l][index]);
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        fprintf(outfile,"%24.18g\n",pdens[ii][index]);
                    }
                }
	    }
	    break;
	case 2:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%24.18g\n",dens[index]);
	        fprintf(outfile,"%24.18g\n",engy[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",momn[l][index]);
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        fprintf(outfile,"%24.18g\n",pdens[ii][index]);
                    }
                }
	    }
	    break;
	case 3:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (k = imin[2]; k <= imax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%24.18g\n",dens[index]);
	        fprintf(outfile,"%24.18g\n",engy[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",momn[l][index]);
                if(eqn_params->multi_comp_non_reactive == YES)
                { 
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        fprintf(outfile,"%24.18g\n",pdens[ii][index]);
                    }
                }
	    }
	}
	fclose(outfile);
}

void G_CARTESIAN::printOneDimStates(char *out_name)
{
	int i, index;
	char filename1[100], filename2[100], filename3[100];
	FILE *outfile;
	double *dens = field.dens;
	double *pres = field.pres;
	double **vel = field.vel;

	sprintf(filename1,"%s/dens%s",out_name,right_flush(front->step,7));
	sprintf(filename1,"%s-dat",filename1);
	sprintf(filename2,"%s/pres%s",out_name,right_flush(front->step,7));
	sprintf(filename2,"%s-dat",filename2);
	sprintf(filename3,"%s/vel%s",out_name,right_flush(front->step,7));
	sprintf(filename3,"%s-dat",filename3);

	outfile = fopen(filename1,"w");

	for (i=0; i<=top_gmax[0]; ++i)
	{
	    index = d_index1d(i,top_gmax);
	    fprintf(outfile,"%24.18g\n", dens[index]);
	}
	fclose(outfile);

	outfile = fopen(filename2,"w");

	for (i=0; i<=top_gmax[0]; ++i)
	{
	    index = d_index1d(i,top_gmax);
	    fprintf(outfile,"%24.18g\n", pres[index]);
	}
	fclose(outfile);

	outfile = fopen(filename3,"w");

	for (i=0; i<=top_gmax[0]; ++i)
	{
	    index = d_index1d(i,top_gmax);
	    fprintf(outfile,"%24.18g\n", vel[0][index]);
	}
	fclose(outfile);
}

void G_CARTESIAN::printInteriorVtk(char *out_name)
{
    	int i,j,k,index,l;
	char filename[100];
	FILE *outfile;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

	sprintf(filename,"%s/state.ts%07d", out_name, front->step);
	if (pp_numnodes() > 1)
		sprintf(filename,"%s-nd%07d", filename, pp_mynode());
	sprintf(filename,"%s.vtk",filename);
	outfile = fopen(filename,"w");

	int data_size = (top_gmax[0]+1)* (top_gmax[1]+1)* (dim==3?(top_gmax[2]+1):1);
	fprintf(outfile, 
		"# vtk DataFile Version 2.0\n"
		"comment line\n"
		"ASCII\n"
		"DATASET STRUCTURED_POINTS\n"
		"DIMENSIONS %d %d %d\n"
		"SPACING %f %f %f\n"
		"ORIGIN %f %f %f\n"
		"CELL_DATA %d\n",
		top_gmax[0]+2, top_gmax[1]+2, dim==3?top_gmax[2]+2:1,
		top_h[0], top_h[1], dim==3?top_h[2]:0,
		top_L[0]-top_h[0]/2, top_L[1]-top_h[1]/2, dim==3?top_L[2]:0,
		data_size);

	fprintf(outfile,"SCALARS DENSITY double 1\nLOOKUP_TABLE default\n");
	for(int i=0; i < data_size; ++i)
	{
	    fprintf(outfile,"%24.18g\n",dens[i]);
	}

	fprintf(outfile,"SCALARS XVEL double 1\nLOOKUP_TABLE default\n");
	for(int i=0; i < data_size; ++i)
	{
	    fprintf(outfile,"%24.18g\n",momn[0][i]/dens[i]);
	}

	fprintf(outfile,"SCALARS YVEL double 1\nLOOKUP_TABLE default\n");
        for(int i=0; i < data_size; ++i)
        {
            fprintf(outfile,"%24.18g\n",momn[1][i]/dens[i]);
        }

	fclose(outfile);
}

void G_CARTESIAN::readInteriorStates(char *restart_name)
{
	FILE *infile;
	int i,j,k,l,index;
	STATE st_tmp;
	char fname[100];
	int		comp;
	EOS_PARAMS	*eos = eqn_params->eos;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
        double **pdens = field.pdens;

	setDomain();
	m_dens[0] = eqn_params->rho1;		
	m_dens[1] = eqn_params->rho2;		
	m_mu[0] = eqn_params->mu1;		
	m_mu[1] = eqn_params->mu2;		
	if (eqn_params->prob_type == FLUID_SOLID_CIRCLE ||
	    eqn_params->prob_type == FLUID_RIGID_BODY ||
	    eqn_params->prob_type == FLUID_CRYSTAL)
	    m_comp[0] = SOLID_COMP;
	else
	    m_comp[0] = GAS_COMP1;
	m_comp[1] = GAS_COMP2;
	m_smoothing_radius = top_h[0] < top_h[1] ? top_h[1] : top_h[0];
	m_smoothing_radius *= 2.0;
	
	st_tmp.dim = eqn_params->dim;

	sprintf(fname,"%s-gas",restart_name);
	infile = fopen(fname,"r");
	
	next_output_line_containing_string(infile,"Interior gas states:");

	switch (dim)
	{
	case 2:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		st_tmp.eos = &(eos[comp]);
	    	
		fscanf(infile,"%lf",&dens[index]);
	    	fscanf(infile,"%lf",&engy[index]);
		st_tmp.dens = dens[index];
		st_tmp.engy = engy[index];
		for (l = 0; l < dim; ++l)
		{
	    	    fscanf(infile,"%lf",&momn[l][index]);
		    st_tmp.momn[l] = momn[l][index];
		}
		pres[index] = EosPressure(&st_tmp);
	    }
	    break;
	case 3:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (k = imin[2]; k <= imax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		st_tmp.eos = &(eos[comp]);

	    	fscanf(infile,"%lf",&dens[index]);
	    	fscanf(infile,"%lf",&engy[index]);
		st_tmp.dens = dens[index];
		st_tmp.engy = engy[index];
		for (l = 0; l < dim; ++l)
		{
	    	    fscanf(infile,"%lf",&momn[l][index]);
		    st_tmp.momn[l] = momn[l][index];
		}
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        fscanf(infile,"%lf",&pdens[ii][index]);
                        st_tmp.pdens[ii] = pdens[ii][index];
                    }
                }
		pres[index] = EosPressure(&st_tmp);
	    }
	}
	fclose(infile);
	scatMeshStates();
	copyMeshStates();
}


void G_CARTESIAN::setAdvectionDt()
{
	double d = (double)dim;
	pp_global_max(&max_speed,1);
	if (max_speed != 0.0)
//	    max_dt = hmin/max_speed/d;
	    max_dt = hmin/max_speed;
	else
	    max_dt = 0.0;
	if (debugging("trace"))
	    printf("In setAdvectionDt: max_dt = %24.18g\n",max_dt);
}	/* end setAdvectionDt */


void G_CARTESIAN::augmentMovieVariables()
{
	int i;
	static HDF_MOVIE_VAR *hdf_movie_var;
	int offset,num_var;

	hdf_movie_var = front->hdf_movie_var;
	offset = front->hdf_movie_var->num_var;
	if (hdf_movie_var == NULL)
	    return initMovieVariables();
	else
	{
	    num_var = offset + dim + 3;
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    hdf_movie_var->num_var = num_var;
	    FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,
				num_var,100,sizeof(char));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,
				num_var,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,
				num_var,sizeof(COMPONENT));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,
				num_var,sizeof(boolean));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,
				num_var,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,
				num_var,sizeof(double));
	    for (i = 0; i < front->hdf_movie_var->num_var; ++i)
	    {
		strcpy(hdf_movie_var->var_name[i],
			front->hdf_movie_var->var_name[i]);
		hdf_movie_var->get_state_var[i] =
			front->hdf_movie_var->get_state_var[i];
		hdf_movie_var->top_var[i] = 
			front->hdf_movie_var->top_var[i];
	    }
	    sprintf(hdf_movie_var->var_name[offset+0],"dens");
	    sprintf(hdf_movie_var->var_name[offset+1],"pres");
	    sprintf(hdf_movie_var->var_name[offset+2],"vort");
	    sprintf(hdf_movie_var->var_name[offset+3],"xvel");
	    sprintf(hdf_movie_var->var_name[offset+4],"yvel");
	    hdf_movie_var->get_state_var[offset+0] = getStateDens;
	    hdf_movie_var->get_state_var[offset+1] = getStatePres;
	    hdf_movie_var->get_state_var[offset+2] = getStateVort;
	    hdf_movie_var->get_state_var[offset+3] = getStateXvel;
	    hdf_movie_var->get_state_var[offset+4] = getStateYvel;
	    if (dim == 3)
	    {
	    	sprintf(hdf_movie_var->var_name[offset+5],"zvel");
	    	hdf_movie_var->get_state_var[offset+5] = getStateZvel;
	    }
	}
	hdf_movie_var->top_var[offset+0] = eqn_params->dens;
	hdf_movie_var->top_var[offset+1] = eqn_params->pres;
	hdf_movie_var->top_var[offset+2] = eqn_params->vort;
	hdf_movie_var->top_var[offset+3] = eqn_params->vel[0];
	hdf_movie_var->top_var[offset+4] = eqn_params->vel[1];
	if (dim == 3)
	    hdf_movie_var->top_var[offset+5] = eqn_params->vel[2];
	FT_FreeThese(2,front->hdf_movie_var->var_name,
			front->hdf_movie_var->top_var);
	FT_FreeThese(1,front->hdf_movie_var);
	front->hdf_movie_var = hdf_movie_var;
	front->hdf_movie_var->num_var = num_var;

}	/* end augmentMovieVariables */

void G_CARTESIAN::initMovieVariables()
{
	static HDF_MOVIE_VAR *hdf_movie_var;
	int n;
	MOVIE_OPTION *movie_option = eqn_params->movie_option;

	if (hdf_movie_var == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    if (eqn_params->tracked == NO)
		hdf_movie_var->untracked = YES;
	    switch (dim)
	    {
	    case 1:
		hdf_movie_var->num_var = n = 0;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,3,100,
				sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,3,
				sizeof(double*));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,3,
				sizeof(COMPONENT));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,3,
				sizeof(boolean));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,3,
				sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,3,
				sizeof(double));
		if (movie_option->plot_dens)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"density");
	    	    hdf_movie_var->get_state_var[n] = getStateDens;
	    	    hdf_movie_var->top_var[n] = eqn_params->dens;
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_dens;
			hdf_movie_var->var_max[n] = movie_option->max_dens;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_pres)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"pressure");
	    	    hdf_movie_var->get_state_var[n] = getStatePres;
	    	    hdf_movie_var->top_var[n] = eqn_params->pres;
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_pres;
			hdf_movie_var->var_max[n] = movie_option->max_pres;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_velo)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"velocity");
	    	    hdf_movie_var->get_state_var[n] = getStateXvel;
	    	    hdf_movie_var->top_var[n] = eqn_params->vel[0];
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_velo;
			hdf_movie_var->var_max[n] = movie_option->max_velo;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
		}
		break;
	    case 2:
		hdf_movie_var->num_var = n = 0;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,5,100,
					sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,5,
					sizeof(double*));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,5,
					sizeof(COMPONENT));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,5,
				sizeof(boolean));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,5,
				sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,5,
				sizeof(double));
		if (movie_option->plot_dens)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"dens");
	    	    hdf_movie_var->get_state_var[n] = getStateDens;
	    	    hdf_movie_var->top_var[n] = eqn_params->dens;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_pres)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"pres");
	    	    hdf_movie_var->get_state_var[n] = getStatePres;
	    	    hdf_movie_var->top_var[n] = eqn_params->pres;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_vort)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"vort");
	    	    hdf_movie_var->get_state_var[n] = getStateVort;
	    	    hdf_movie_var->top_var[n] = eqn_params->vort;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_velo)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"xvel");
	    	    hdf_movie_var->get_state_var[n] = getStateXvel;
	    	    hdf_movie_var->top_var[n] = eqn_params->vel[0];
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
	    	    sprintf(hdf_movie_var->var_name[n],"yvel");
	    	    hdf_movie_var->get_state_var[n] = getStateYvel;
	    	    hdf_movie_var->top_var[n] = eqn_params->vel[1];
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		break;
	    case 3:
		hdf_movie_var->num_var = n = 0;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,15,100,
					sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,15,
					sizeof(double*));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->idir,15,
					sizeof(int));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,15,
					sizeof(COMPONENT));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,15,
				sizeof(boolean));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,15,
				sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,15,
				sizeof(double));
		if (movie_option->plot_cross_section[0])
		{
		    if (movie_option->plot_dens)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"dens-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStateDens;
	    		hdf_movie_var->top_var[n] = eqn_params->dens;
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = eqn_params->pres;
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-yz-y");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[1];
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-yz-z");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[2];
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_vort)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"vort-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvort;
	    		hdf_movie_var->top_var[n] = eqn_params->vort3d[0];
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
		if (movie_option->plot_cross_section[1])
		{
		    if (movie_option->plot_dens)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"dens-xz");
	    	    	hdf_movie_var->get_state_var[n] = getStateDens;
	    		hdf_movie_var->top_var[n] = eqn_params->dens;
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-xz");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = eqn_params->pres;
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xz-x");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[0];
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xz-z");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[2];
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_vort)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"vort-xz");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvort;
	    		hdf_movie_var->top_var[n] = eqn_params->vort3d[1];
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
		if (movie_option->plot_cross_section[2])
		{
		    if (movie_option->plot_dens)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"dens-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStateDens;
	    		hdf_movie_var->top_var[n] = eqn_params->dens;
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = eqn_params->pres;
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xy-x");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[0];
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xy-y");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[1];
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_vort)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"vort-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvort;
	    		hdf_movie_var->top_var[n] = eqn_params->vort3d[2];
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
	    }
	}
	front->hdf_movie_var = hdf_movie_var;
}	/* end initMovieVariables */

double G_CARTESIAN::getVorticityX(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dy,dz;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dy = top_h[1];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j-1,k,top_gmax);
	index01 = d_index3d(i,j+1,k,top_gmax);
	index10 = d_index3d(i,j,k-1,top_gmax);
	index11 = d_index3d(i,j,k+1,top_gmax);
	v00 = -momn[2][index00]/dens[index00];
	v01 =  momn[2][index01]/dens[index01];
	v10 =  momn[1][index10]/dens[index10];
	v11 = -momn[1][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dz + (v10 + v11)/2.0/dy;
	return vorticity;
}	/* end getVorticityX */

double G_CARTESIAN::getVorticityY(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dz;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dx = top_h[0];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j,k-1,top_gmax);
	index01 = d_index3d(i,j,k+1,top_gmax);
	index10 = d_index3d(i-1,j,k,top_gmax);
	index11 = d_index3d(i+1,j,k,top_gmax);
	v00 = -momn[0][index00]/dens[index00];
	v01 =  momn[0][index01]/dens[index01];
	v10 =  momn[2][index10]/dens[index10];
	v11 = -momn[2][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dx + (v10 + v11)/2.0/dz;
	return vorticity;
}	/* end getVorticityY */

double G_CARTESIAN::getVorticityZ(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i-1,j,k,top_gmax);
	index01 = d_index3d(i+1,j,k,top_gmax);
	index10 = d_index3d(i,j-1,k,top_gmax);
	index11 = d_index3d(i,j+1,k,top_gmax);
	v00 = -momn[1][index00]/dens[index00];
	v01 =  momn[1][index01]/dens[index01];
	v10 =  momn[0][index10]/dens[index10];
	v11 = -momn[0][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticityZ */

double G_CARTESIAN::getVorticity(int i, int j)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index2d(i,j,top_gmax);
	index00 = d_index2d(i-1,j,top_gmax);
	index01 = d_index2d(i+1,j,top_gmax);
	index10 = d_index2d(i,j-1,top_gmax);
	index11 = d_index2d(i,j+1,top_gmax);
	v00 = -momn[1][index00]/dens[index00];
	v01 =  momn[1][index01]/dens[index01];
	v10 =  momn[0][index10]/dens[index10];
	v11 = -momn[0][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticity */

void G_CARTESIAN::copyMeshStates()
{
	int i,j,k,l,index;
	double **vel = eqn_params->vel;
	double **mom = eqn_params->mom;
	double *dens = eqn_params->dens;
        double **pdens = eqn_params->pdens;
	double *pres = eqn_params->pres;
	double *engy = eqn_params->engy;
	double *vort = eqn_params->vort;
	double **vort3d = eqn_params->vort3d;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		for (l = 0; l < dim; ++l)
		    vel[l][index] = mom[l][index]/dens[index];	
	    }
	    FT_ParallelExchGridArrayBuffer(dens,front);
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                int ii;
                for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    FT_ParallelExchGridArrayBuffer(pdens[ii],front);
                }
            }
	    FT_ParallelExchGridArrayBuffer(pres,front);
	    FT_ParallelExchGridArrayBuffer(engy,front);
	    FT_ParallelExchGridArrayBuffer(mom[0],front);
	    FT_ParallelExchGridArrayBuffer(vel[0],front);
	    break;
	case 2:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		for (l = 0; l < dim; ++l)
		    vel[l][index] = mom[l][index]/dens[index];	
	    }
	    FT_ParallelExchGridArrayBuffer(dens,front);
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                int ii;
                for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    FT_ParallelExchGridArrayBuffer(pdens[ii],front);
                }
            }
	    FT_ParallelExchGridArrayBuffer(pres,front);
	    FT_ParallelExchGridArrayBuffer(engy,front);
	    FT_ParallelExchGridArrayBuffer(vort,front);
	    for (l = 0; l < dim; ++l)
	    {
	    	FT_ParallelExchGridArrayBuffer(mom[l],front);
	    	FT_ParallelExchGridArrayBuffer(vel[l],front);
	    }
	    break;
	case 3:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (k = imin[2]; k <= imax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
		for (l = 0; l < dim; ++l)
		    vel[l][index] = mom[l][index]/dens[index];	
	    }
	    FT_ParallelExchGridArrayBuffer(dens,front);
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                int ii;
                for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    FT_ParallelExchGridArrayBuffer(pdens[ii],front);
                }
            }
	    FT_ParallelExchGridArrayBuffer(pres,front);
	    FT_ParallelExchGridArrayBuffer(engy,front);
            FT_ParallelExchGridArrayBuffer(vort3d[0],front);
            FT_ParallelExchGridArrayBuffer(vort3d[1],front);
            FT_ParallelExchGridArrayBuffer(vort3d[2],front);
	    for (l = 0; l < dim; ++l)
	    {
	    	FT_ParallelExchGridArrayBuffer(mom[l],front);
		state_reflect(l,mom[l]);	//Dan	FIXME
	    	FT_ParallelExchGridArrayBuffer(vel[l],front);
		state_reflect(l,vel[l]);	//Dan	FIXME
	    }
	    break;
	}
}	/* end copyMeshStates */

void G_CARTESIAN::computeTBL()	//Dan
{
    	INTERFACE	*intfc = front->grid_intfc;
	int		dir, side;
	static SWEEP	st_field;
	static boolean	first = YES;
	double delta_t;

	delta_t = m_dt;

	if (first)
	{
	    allocMeshVst(&st_field);
	    first = NO;
	}

	copyToMeshVst(&st_field);

	oned_turbulence_boundary_layer(&st_field,dir,side);

	scatMeshVst(&st_field);
	copyFromMeshVst(st_field);
}

void G_CARTESIAN::oned_turbulence_boundary_layer(
	SWEEP	*m_vst,
	int	dir,
	int	side)
{
	if (dim == 2)
	    TBL2D(m_vst,dir,side);
	else if (dim == 3)
	    TBL3D(m_vst,dir,side);
	else
	{
	    printf("Turbulent boundary layer is not implemented for %dD.\n", dim);
	    exit(0);
	}
}

void G_CARTESIAN::TBL2D(
	SWEEP	*m_vst,
	int	dir,
	int	side)
{
    	printf("Turbulent boundary layer is not implemented for 2D yet.\n");
	clean_up(ERROR);
}

void G_CARTESIAN::TBL3D(
	SWEEP	*m_vst,
	int	dir,
	int	side)
{
    	INTERFACE	*intfc = front->grid_intfc;
    	int		i, j, k;
	int		ic[MAXD];

	if (dir == 0)
	{
	    for (k = imin[2]; k <= imax[2]; k++)
	    for (j = imin[1]; j <= imax[1]; j++)
	    {
		ic[1] = j;
		ic[2] = k;
		if (side == 0)
		    ic[0] = imin[0];
		else if (side == 1)
		    ic[0] = imax[0];
		else
		{
		    printf("In TBL3D, side should be 0 or 1.\n");
		    clean_up(ERROR);
		}
		if (rect_boundary_type(intfc,dir,side) == REFLECTION_BOUNDARY)
		    TBLsolver(m_vst,ic,dir,side);
	    }
	}

	if (dir == 1)
	{
	    for (k = imin[2]; k <= imax[2]; k++)
	    for (i = imin[0]; i <= imax[0]; i++)
	    {
		ic[0] = i;
		ic[2] = k;
		if (side == 0)
		    ic[1] = imin[1];
		else if (side == 1)
		    ic[1] = imax[1];
		else
		{
		    printf("In TBL3D, side should be 0 or 1.\n");
		    clean_up(ERROR);
		}
		if (rect_boundary_type(intfc,dir,side) == REFLECTION_BOUNDARY)
		    TBLsolver(m_vst,ic,dir,side);
	    }
	}

	if (dir == 0)
	{
	    for (j = imin[1]; j <= imax[1]; j++)
	    for (i = imin[0]; i <= imax[0]; i++)
	    {
		ic[0] = i;
		ic[1] = j;
		if (side == 0)
		    ic[2] = imin[2];
		else if (side == 1)
		    ic[2] = imax[2];
		else
		{
		    printf("In TBL3D, side should be 0 or 1.\n");
		    clean_up(ERROR);
		}
		if (rect_boundary_type(intfc,dir,side) == REFLECTION_BOUNDARY)
		    TBLsolver(m_vst,ic,dir,side);
	    }
	}
}

void G_CARTESIAN::TBLsolver(
	SWEEP	*m_vst,
	int	*icoords,
	int	idir,
	int	nb)
{
    	int		i;
	int		ic[MAXD];
	int		index_in, index_out;
	COMPONENT	comp_in, comp_out;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	GRID_DIRECTION	dir;
	STATE		*state;
	HYPER_SURF	*hs;
	double		coords[MAXD], crx_coords[MAXD];
	STATE		st_tmp;
	double		dt = m_dt;

	if (nb == 0)
	    dir = ldir[idir];
	else if (nb == 1)
	    dir = rdir[idir];
	else
	{
	    printf("Incorrect nb in TBLsolver.\n");
	    clean_up(ERROR);
	}

	for (i = 0; i < dim; i++)
	    ic[i] = icoords[i];
	index_in = d_index(ic,top_gmax,dim);
	comp_in = cell_center[index_in].comp;


	for (i = 0; i < dim; i++)
	    coords[i] = top_L[i] + icoords[i] * top_h[i];
	double y0, yw, tblflux[4];
	yw = nb == 0 ? coords[idir] - top_h[idir] : coords[idir] + top_h[idir];
	y0 = coords[idir];

	st_tmp.dens = m_vst->dens[index_in];
        st_tmp.pdens[0] = m_vst->pdens[0][index_in];
        st_tmp.pdens[1] = m_vst->pdens[1][index_in];
	st_tmp.engy = m_vst->engy[index_in];
	st_tmp.pres = m_vst->pres[index_in];
	for (i = 0; i < dim; i++)
	    st_tmp.momn[i] = m_vst->momn[i][index_in];
	st_tmp.eos = &(eqn_params->eos[comp_in]);

        if(fabs(st_tmp.momn[0]) < 1e-12 && fabs(st_tmp.momn[1]) < 1e-12 && fabs(st_tmp.momn[2]) < 1e-12)
            return;        

	computeTBLflux(st_tmp,tblflux,y0,yw,idir);

	m_vst->momn[0][index_in] += tblflux[2]*dt/top_h[idir];
	m_vst->momn[1][index_in] += tblflux[3]*dt/top_h[idir];
	m_vst->momn[2][index_in] += tblflux[1]*dt/top_h[idir];
	m_vst->engy[index_in] += tblflux[0]*dt/top_h[idir];
}

#define	WALLMODEL_STRETCHING	1.1	//FIXME
#define	WALLMODEL_APLUS		17.0	//FIXME
#define	WALLMODEL_KAPPA		0.41	//FIXME
#define	WALLMODEL_PRT		0.9	//FIXME
#define	Wall_Temp		0.1	//FIXME
void G_CARTESIAN::computeTBLflux(
	STATE	st,
	double	*flux,
	double	y0,
	double	yw,
	int	idir)	//Dan
{
    	int NY, i;
	double delta = fabs(y0 - yw);
	double ymatch = delta;
	double ycrans[200];	//FIXME
	double dyrans[200], yfrans[200];
	double wallunit = 0.0001;	//FIXME
	double dyrans0;
	double rnm1 = 1.0;	//??	FIXME
	double Umatch, Tw, Tmatch, rho0, C_P;
	double Urans[200], Trans[200], MUrans[200], Krans[200], Rrans[200], MUTrans[200];	//FIXME
	double ap[200], qp[200], aw[200], ae[200];	//FIXME
	double tauw, tauwold, qw, qwold, muw, rhow, utau, lv, kw;
	double damping;
	double intfac, muint, kint, mutint, factor, uint;
	double u, v, w;

	if (dim == 2)
	{
	    u = st.momn[0]/st.dens;
	    v = st.momn[1]/st.dens;
	    w = 0.0;
	    if (idir == 0)
		Umatch = fabs(v);
	    else
		Umatch = fabs(u);
	}
	else if (dim == 3)
	{
	    u = st.momn[0]/st.dens;
	    v = st.momn[1]/st.dens;
	    w = st.momn[2]/st.dens;
	    if (idir == 0)
		Umatch = sqrt(v*v + w*w);
	    else if (idir == 1)
		Umatch = sqrt(u*u + w*w);
	    else if (idir == 2)
		Umatch = sqrt(u*u + v*v);
	}
	else if (dim == 1)
	{
	    printf("TBL is not working for 1D simulation.\n");
	    clean_up(ERROR);
	}

	Tmatch = EosTemperature(&st);
	if(Tmatch >= Wall_Temp)
	    Tw = Wall_Temp;
	else Tw = Tmatch;
	rho0 = st.dens;
	C_P = EosCP(&st);	//FIXME

	NY = 1;
	ycrans[0] = 0.5;
	dyrans0 = wallunit;
	while (ycrans[NY-1] < ymatch/dyrans0)
	{
	    ycrans[NY] = ycrans[NY-1] + 0.5*rnm1*(1.0 + WALLMODEL_STRETCHING);	//??	FIXME
	    rnm1 *= WALLMODEL_STRETCHING;
	    NY++;
	}
	dyrans0 = ymatch/ycrans[NY-1];

	ycrans[0] *= dyrans0;
	dyrans[0] = dyrans0;
	yfrans[0] = 0.0;
	for (i = 1; i < NY; i++)
	{
	    ycrans[i] *= dyrans0;
	    dyrans[i] = dyrans[i-1]*WALLMODEL_STRETCHING;
	    yfrans[i] = yfrans[i-1] + dyrans[i-1];
	}

    	//initial guesses -- linear profiles
	for (i = 0; i < NY; i++)
	{
	    Urans[i] = Umatch/ymatch*ycrans[i];
	    Trans[i] = Tw + (Tmatch - Tw)/ymatch*ycrans[i];
	}
//	dynamic_viscosity_thermalconduct(statep,Tw,&muw,&kw);	//get muw and kw	FIXME
	EosViscTherm(&st,&muw,&kw);	//FIXME
	tauw = muw*Urans[0]/ycrans[0];
	rhow = rho0*Tmatch/Tw;
	utau = sqrt(tauw/rhow);
	lv = muw/(rhow*utau);
	qw = - kw*(Trans[0] - Tw)/ycrans[0];

	int iter = 0;
	double converged = HUGE_VAL;
	while (converged > 1e-8)	//FIXME
	{
	    tauwold = tauw;
	    qwold = qw;

	    for (i = 0; i < NY; i++)
	    {
//		dynamic_viscosity_thermalconduct(statep,Trans[j],&MUrans[j],&Krans[j]);
//		MUrans[i] = ??;	//FIXME
//		Krans[i] = ??;	//FIXME
		EosViscTherm(&st,&MUrans[i],&Krans[i]);	//FIXME
		damping = 1.0 - exp( - ycrans[i] / ( lv * WALLMODEL_APLUS)) ;
		damping = damping*damping;
		Rrans[i] =  rho0*Tmatch/Trans[i];
		MUTrans[i] = WALLMODEL_KAPPA * ycrans[i] * sqrt( Rrans[i] * tauw ) * damping ;
		ap[i] = 0.0 ;
		qp[i] = 0.0 ;   //  if LHS not zero, then should set qp to LHS * dyrans[j]...
	    }

	    //  Set up momentum equation using finite volume

	    ap[0] = -muw/ycrans[0];
	    for (i = 1 ; i < NY ; i++) 
	    {
		intfac = ( yfrans[i] - ycrans[i-1] ) / ( ycrans[i] - ycrans[i-1] );
		muint  = ( 1.0 - intfac ) * MUrans[i-1]  + intfac * MUrans[i];
		kint = (1.0 - intfac ) * Krans[i-1] + intfac*Krans[i];
		mutint = ( 1.0 - intfac ) * MUTrans[i-1] + intfac * MUTrans[i];
		factor = ( muint + mutint ) / ( ycrans[i] - ycrans[i-1] );
		aw[i] = factor;
		ap[i] -= factor;
		ae[i-1] = factor;
		ap[i-1] -= factor;
	    }
	    Urans[NY-1] = Umatch;

	    //  Solve momentum equation by single TDMA sweep
	    
	    for (i = 1 ; i < NY-1 ; i++) 
	    {
		factor = aw[i] / ap[i-1];
		ap[i] -= factor * ae[i-1];
		qp[i] -= factor * qp[i-1];
	    }
	    for (i = NY-2 ; i >= 0 ; i--)
		Urans[i] = ( qp[i] - ae[i] * Urans[i+1] ) / ap[i];
	    
	    //  Auxiliary stuff
	    tauw = muw * Urans[0] / ycrans[0];
	    utau = sqrt( tauw / rhow );
	    lv = muw / ( rhow * utau );
	    for (i = 0; i < NY; i++) 
	    {
		damping = 1.0 - exp(-ycrans[i] / ( lv * WALLMODEL_APLUS ));
		damping = damping*damping;
		MUTrans[i] = WALLMODEL_KAPPA * ycrans[i] * sqrt( Rrans[i] * tauw ) * damping;
		ap[i] = 0.0;
		qp[i] = 0.0;   //  if LHS not zero, then should set qp to LHS * dyrans[j]...
	    }
	
	    //  Set up energy equation using finite volume

	    ap[0] = - kw/ycrans[0];
	    qp[0] = - kw/ycrans[0]*Tw;

	    for (i = 1 ; i < NY ; i++) 
	    {
		intfac = ( yfrans[i] - ycrans[i-1] ) / ( ycrans[i] - ycrans[i-1] );
		muint  = ( 1.0 - intfac ) * MUrans[i-1]  + intfac * MUrans[i];
		mutint = ( 1.0 - intfac ) * MUTrans[i-1] + intfac * MUTrans[i];
		kint = (1.0 - intfac ) * Krans[i-1] + intfac*Krans[i];
		factor = ( kint + C_P/WALLMODEL_PRT * mutint ) / ( ycrans[i] - ycrans[i-1] );
		aw[i] = factor;
		ap[i] -= factor;
		ae[i-1] = factor;
		ap[i-1] -= factor;
		uint = ( 1.0 - intfac ) * Urans[i-1] + intfac * Urans[i];
		factor = ( muint + mutint ) * uint * ( Urans[i] - Urans[i-1] ) / ( ycrans[i] - ycrans[i-1] );
		qp[i] += factor;
		qp[i-1] -= factor;
	    }
	    Trans[NY-1] = Tmatch;

	    //  Solve energy equation by single TDMA sweep
	    for (i = 1 ; i < NY-1 ; i++) 
	    {
		factor = aw[i] / ap[i-1];
		ap[i] -= factor * ae[i-1];
		qp[i] -= factor * qp[i-1];

	    }
	    for (i = NY-2 ; i >= 0 ; i--)
            {
		Trans[i] = (qp[i] - ae[i] * Trans[i+1]) / ap[i];
            }

	    if(Tmatch < Wall_Temp)
		Tw = Trans[0];
	    tauw = muw * Urans[0] / ycrans[0];
	    rhow = rho0 * Tmatch / Tw;
	    utau = sqrt(tauw / rhow);
	    lv = muw / (rhow * utau);
	    qw = kw * (Trans[0] - Tw) / ycrans[0];   // wall heat-flux


	    converged = fabs((tauw - tauwold) / (tauwold + MACH_EPS));

	    if( fabs((qw - qwold) / (qwold + MACH_EPS)) < converged)
		converged = fabs((qw - qwold) / (qwold + MACH_EPS));
	    iter++; 
	}
	
	flux[0] = - qw;
	flux[1] = - tauw / Umatch * u;
	flux[2] = - tauw / Umatch * v;
	flux[3] = - muw * 4.0/3.0 * w/delta;

}

void G_CARTESIAN::compSGS2D(SWEEP *m_vst)
{

        int i, j;
        int index0, index1, index2, index3, index4;
        int index00, index10, index01, index11;
        double *u, *v;
        double ux, uy, vx, vy;
        double *tx, *ty, *cx, *cy, *cx0, *cy0;
        double *s, *s11, *s12, *s22;
        double *cp, *temp, *conc, *conc0;
        double sum_rho_u, sum_rho_v, sum_rho_uu, sum_rho_vv, sum_rho_uv;
        double sum_s11, sum_s12, sum_rhoss11, sum_rhoss12, sum_rhoss, sum_rho, \
               sum_s, sum_ss, sum_cp, sum_tx, sum_ty, sum_rho_t, sum_rho_cp,\
               sum_rho_cp_s_tx, sum_rho_cp_s_ty, sum_rho_cp_tu, sum_rho_cp_tv, \
               sum_cx, sum_cy, sum_cx0, sum_cy0, sum_rho_s_cx, sum_rho_s_cy, \
               sum_rho_s_cx0, sum_rho_s_cy0, sum_rho_uc,sum_rho_vc, sum_rho_uc0, \
               sum_rho_vc0, sum_rho_c, sum_rho_c0, sum_u;
        double MA11, MA12, L11, L12, L22, LA11, LI, MI, MH1, MH2, LH1, LH2;
        double MC1, MC2, MC10, MC20, LC1, LC2, LC10, LC20;
	double CS, CI, Prt, Sct, Sct0;
        double *cs, *ci, *prt, *sct, *sct0;
        int    ii, jj, iiii, jjjj;
        int    NB = 2.0;
        int    NBC = pow(NB, 2);
        double CS_deno, CS_nume, CI_deno, CI_nume, Prt_nume, Prt_deno, Sct_nume, \
               Sct_deno, Sct0_nume, Sct0_deno;
        double delta2, tdelta2;

        delta2 = sqr(pow((top_h[0]*top_h[1]),(1.0/2.0)));  // filter width 
        tdelta2 = sqr(pow((NB*top_h[0]*NB*top_h[1]),(1.0/2.0)));  // test Filter width

        int size = (top_gmax[0]+1)*(top_gmax[1]+1);
        FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tx,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ty,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&cx0,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&cx,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&cy0,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&cy,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cp,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&temp,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&conc0,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&conc,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cs,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ci,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&prt,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&sct0,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&sct,size,sizeof(double));

        STATE st;
        st.dim = eqn_params->dim;
        EOS_PARAMS *eos = eqn_params->eos;
	double *dens = m_vst->dens;
	double *pres = m_vst->pres;
	double **momn = m_vst->momn;
	double *engy = m_vst->engy;
	double **pdens = m_vst->pdens;

        int comp;
	bool subgrid_con = eqn_params->subgrid_con;
	bool subgrid_md = eqn_params->subgrid_md;
        bool subgrid_vis = eqn_params->subgrid_vis;
        
        for(j =imin[1]-1; j <= imax[1]+1; j++)
        for(i =imin[0]-1; i <= imax[0]+1; i++)

        {
            index0 = d_index2d(i,j,top_gmax);

            comp = getComponent(index0);
	    if (comp == 0 && m_vst->comp[index0] != 0)
		comp = m_vst->comp[index0];
	    else if (comp == 0 && m_vst->comp[index0] == 0)
		continue;

            st.eos = &(eos[comp]);
            st.dens = dens[index0];
	    if(eqn_params->multi_comp_non_reactive == YES)
	    {
		int ii;
		for(ii = 0; ii < eqn_params->n_comps; ii++)
		{
		    st.pdens[ii] = pdens[ii][index0];
		}
	    }
            st.engy = engy[index0];
            st.momn[0] = momn[0][index0];
            st.momn[1] = momn[1][index0];
            st.pres = EosPressure(&st);
	
	    u[index0] = momn[0][index0]/dens[index0];
	    v[index0] = momn[1][index0]/dens[index0];
	    conc0[index0] = pdens[0][index0]/dens[index0];
	    conc[index0]  = pdens[1][index0]/dens[index0];
            cp[index0] = EosCP(&st);
            temp[index0] = EosTemperature(&st);
        }

        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
        {
	   
            index0 = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);

            ux = (u[index2] - u[index1]) / (2.0*top_h[0]);
            uy = (u[index4] - u[index3]) / (2.0*top_h[1]);
            vx = (v[index2] - v[index1]) / (2.0*top_h[0]);
            vy = (v[index4] - v[index3]) / (2.0*top_h[1]);

            tx[index0] = (temp[index2] - temp[index1]) / (2.0*top_h[0]);
            ty[index0] = (temp[index4] - temp[index3]) / (2.0*top_h[1]);
		
	    cx0[index0] = (conc0[index2] - conc0[index1]) / (2.0*top_h[0]);
            cx[index0] = (conc[index2] - conc[index1]) / (2.0*top_h[0]);
	    cy0[index0] = (conc0[index4] - conc0[index3]) / (2.0*top_h[1]);
	    cy[index0] = (conc[index4] - conc[index3]) / (2.0*top_h[1]);

            s11[index0] = ux;
            s12[index0] = (uy + vx)/2.0;
            s22[index0] = vy;
            s[index0] = sqrt(2*( sqr(s11[index0]) + sqr(s12[index0])
                         + sqr(s22[index0]) + sqr(s12[index0]) ));

        }


        for (j = 0; j <= ((imax[1]-imin[1]+1)/NB)-1; j++)
        for (i = 0; i <= ((imax[0]-imin[0]+1)/NB)-1; i++)
        {

            jj = (NB*j) + imin[1];
            ii = (NB*i) + imin[0];

            sum_rho = 0.0;
            sum_u = 0.0;
            sum_cp = 0.0;
            sum_rho_u = sum_rho_v = 0.0;
            sum_rho_uu = sum_rho_vv = sum_rho_uv = 0.0;
            sum_s11 = sum_s12 = 0.0;
            sum_rhoss11 = sum_rhoss12 = sum_rhoss = 0.0;
            sum_s = sum_ss = 0.0;
            sum_tx = sum_ty = 0.0;
            sum_rho_t = 0.0;
            sum_rho_cp = sum_rho_cp_s_tx = sum_rho_cp_s_ty = 0.0;
            sum_cx = sum_cy = sum_cx0 = sum_cy0 = 0.0;
            sum_rho_s_cx = sum_rho_s_cy = sum_rho_s_cx0 = sum_rho_s_cy0 = 0.0;
            sum_rho_cp_tu = sum_rho_cp_tv = 0.0;
            sum_rho_uc =  sum_rho_vc = sum_rho_uc0 = sum_rho_vc0 = 0.0;
            sum_rho_c = sum_rho_c0 = 0.0;

	    double max_CS = 0.0;
            double max_CI = 0.0;
            double max_Prt = 0.0;
            double max_Sct = 0.0;
            double max_Sct0 = 0.0;


            for(jjjj = jj; jjjj < jj+NB; jjjj++)
            for(iiii = ii; iiii < ii+NB; iiii++)
            {

                index0 = d_index2d(iiii, jjjj, top_gmax);
                sum_u += u[index0];
                sum_rho_u += dens[index0]*u[index0];
                sum_rho_v += dens[index0]*v[index0];
                sum_rho_uu += dens[index0]*u[index0]*u[index0];
                sum_rho_vv += dens[index0]*v[index0]*v[index0];
                sum_rho_uv += dens[index0]*u[index0]*v[index0];
                sum_s11 += 0.5*(s11[index0]-s22[index0]);
                sum_s12 += s12[index0];
                sum_rhoss11 += dens[index0]*s[index0]*0.5*(s11[index0]-s22[index0]);
                sum_rhoss12 += dens[index0]*s[index0]*s12[index0];
                sum_rhoss += dens[index0]*s[index0]*s[index0];
                sum_s += s[index0];
                sum_ss += s[index0]*s[index0];
                sum_rho += dens[index0];

                sum_cp += cp[index0];
                sum_tx += tx[index0];
                sum_ty += ty[index0];
                sum_rho_t += dens[index0]*temp[index0];
                sum_rho_cp += dens[index0]*cp[index0];
                sum_rho_cp_s_tx += dens[index0]*cp[index0]*s[index0]*tx[index0];
                sum_rho_cp_s_ty += dens[index0]*cp[index0]*s[index0]*ty[index0];
                sum_rho_cp_tu += dens[index0]*cp[index0]*temp[index0]*u[index0];
                sum_rho_cp_tv += dens[index0]*cp[index0]*temp[index0]*v[index0];

		sum_cx += cx[index0];
                sum_cy += cy[index0];
                sum_cx0 += cx0[index0];
                sum_cy0 += cy0[index0];
                sum_rho_s_cx += dens[index0]*s[index0]*cx[index0];
                sum_rho_s_cy += dens[index0]*s[index0]*cy[index0];
                sum_rho_s_cx0 += dens[index0]*s[index0]*cx0[index0];
                sum_rho_s_cy0 += dens[index0]*s[index0]*cy0[index0];
                sum_rho_uc += dens[index0]*u[index0]*conc[index0];
                sum_rho_vc += dens[index0]*v[index0]*conc[index0];
                sum_rho_uc0 += dens[index0]*u[index0]*conc0[index0];
                sum_rho_vc0 += dens[index0]*v[index0]*conc0[index0];
                sum_rho_c += dens[index0]*conc[index0];
                sum_rho_c0 += dens[index0]*conc0[index0];

            }

            MA11 = 2.0*delta2*(sum_rhoss11/NBC)
                        - 2.0*tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_s11/NBC);
            MA12 = 2.0*delta2*sum_rhoss12/(NBC)
                        - 2.0*tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_s12/NBC);

            MI = -2.0*delta2*(sum_rhoss/NBC) + 2.0*tdelta2*(sum_rho/NBC)*(sum_ss/NBC);

            L11 = sum_rho_uu/NBC - (sum_rho_u/NBC)*(sum_rho_u/NBC)/(sum_rho/NBC);
            L12 = sum_rho_uv/NBC - (sum_rho_u/NBC)*(sum_rho_v/NBC)/(sum_rho/NBC);
            L22 = sum_rho_vv/NBC - (sum_rho_v/NBC)*(sum_rho_v/NBC)/(sum_rho/NBC);
            LA11 = 0.5*(L11 - L22);
            LI = 0.5*(L11 + L22);
            MH1 = delta2*(sum_rho_cp_s_tx/NBC) \
                    - tdelta2*(sum_rho/NBC)*(sum_cp/NBC)*(sum_s/NBC)*(sum_tx/NBC); 
            MH2 = delta2*(sum_rho_cp_s_ty/NBC) \
                    - tdelta2*(sum_rho/NBC)*(sum_cp/NBC)*(sum_s/NBC)*(sum_ty/NBC);
            LH1 = sum_rho_cp_tu/NBC  \
                    - (sum_rho_cp/NBC)*(sum_rho_t/NBC)*(sum_rho_u/NBC)/((sum_rho/NBC)*(sum_rho/NBC));
            LH2 = sum_rho_cp_tv/NBC  \
                    - (sum_rho_cp/NBC)*(sum_rho_t/NBC)*(sum_rho_v/NBC)/((sum_rho/NBC)*(sum_rho/NBC));

	    MC1 = delta2*(sum_rho_s_cx/NBC)  \
                    - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cx/NBC);
            MC2 = delta2*(sum_rho_s_cy/NBC)  \
                    - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cy/NBC);
            MC10 = delta2*(sum_rho_s_cx0/NBC)  \
                    - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cx0/NBC);
            MC20 = delta2*(sum_rho_s_cy0/NBC)  \
                    - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cy0/NBC);

	    LC1 = sum_rho_uc/NBC - (sum_rho_u/NBC)*(sum_rho_c/NBC)/(sum_rho/NBC);
            LC2 = sum_rho_vc/NBC - (sum_rho_v/NBC)*(sum_rho_c/NBC)/(sum_rho/NBC);
            LC10 = sum_rho_uc0/NBC - (sum_rho_u/NBC)*(sum_rho_c0/NBC)/(sum_rho/NBC);
            LC20 = sum_rho_vc0/NBC - (sum_rho_v/NBC)*(sum_rho_c0/NBC)/(sum_rho/NBC);

            CS_deno = MA11*MA11 + MA12*MA12;
            CS_nume = LA11*MA11 + L12*MA12;
            CI_deno = MI;
            CI_nume = LI;
            Prt_nume = LH1*MH1 + LH2*MH2;
            Prt_deno = MH1*MH1 + MH2*MH2;
	    Sct_nume  = LC1*MC1 + LC2*MC2;
            Sct_deno = MC1*MC1 + MC2*MC2;
            Sct0_nume = LC10*MC10 + LC20*MC20;
            Sct0_deno = MC10*MC10 + MC20*MC20;

            CS = CS_nume/CS_deno;
            CI = CI_nume/CI_deno;
            Prt = Prt_nume/Prt_deno;
	    Sct = Sct_nume/Sct_deno;
            Sct0 = Sct0_nume/Sct0_deno;

            if(CS < 0.0 || CS_deno < 1e-16)
                CS = 0.0;
            if(CS > max_CS)
                max_CS = CS;
            if(CI < 0.0 || CI_deno < 1e-16)
                CI = 0.0;
            if(CI > max_CI)
                max_CI = CI;
            if(Prt < 0.0 || Prt_deno < 1e-16)
                Prt = 0.0;
            if(Prt > max_Prt)
                max_Prt = Prt;
	    if(Sct < 0.0 || Sct_deno < 1e-16)
                Sct = 0.0;
            if(Sct > max_Sct)
                max_Sct = Sct;
            if(Sct0 < 0.0 || Sct0_deno < 1e-16)
                Sct0 = 0.0;
            if(Sct0 > max_Sct0)
                max_Sct0 = Sct0;

            
            index00 = d_index2d(ii, jj, top_gmax);
            index10 = d_index2d(ii+1, jj, top_gmax);
            index01 = d_index2d(ii, jj+1, top_gmax);
            index11 = d_index2d(ii+1, jj+1, top_gmax);


            cs[index00] = CS;
            cs[index10] = CS;
            cs[index01] = CS;
            cs[index11] = CS;
            ci[index0] = CI;
            ci[index10] = CI;
            ci[index01] = CI;
            ci[index11] = CI;
            prt[index00] = Prt;
            prt[index10] = Prt;
            prt[index01] = Prt;
            prt[index11] = Prt;
	    sct[index00] = Sct;
            sct[index10] = Sct;
            sct[index01] = Sct;
            sct[index11] = Sct;
            sct0[index00] = Sct0;
            sct0[index10] = Sct0;
            sct0[index01] = Sct0;
            sct0[index11] = Sct0;

        }

        double ***tau;
        double *qt1, *qt2, *qp1, *qp2, *qp3, *qp4;
        double dtemx, dtemy, concx, concy, conc0x, conc0y;

        FT_VectorMemoryAlloc((POINTER*)&qt1,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qt2,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&qp1,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qp2,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qp3,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qp4,size,sizeof(double));
	FT_TriArrayMemoryAlloc((POINTER*)&tau,dim,dim,size,sizeof(double));

	for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
	{

            index0  = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);
 
            tau[0][0][index0] = -2.0*cs[index0]*delta2*dens[index0]*s[index0]*0.5*(s11[index0]-s22[index0]) + (2.0/2.0)*ci[index0]*delta2*dens[index0]*s[index0]*s[index0];
            tau[0][1][index0] = -2.0*cs[index0]*delta2*dens[index0]*s[index0]*s12[index0];
            tau[1][1][index0] = -2.0*cs[index0]*delta2*dens[index0]*s[index0]*0.5*(s22[index0]-s11[index0]) + (2.0/2.0)*ci[index0]*delta2*dens[index0]*s[index0]*s[index0];


            dtemx = (temp[index2] - temp[index1])/(2.0*top_h[0]);
            dtemy = (temp[index4] - temp[index3])/(2.0*top_h[1]);
	    concx = (conc[index2] - conc[index1])/(2.0*top_h[0]);
            concy = (conc[index4] - conc[index3])/(2.0*top_h[1]);
            conc0x = (conc0[index2] - conc0[index1])/(2.0*top_h[0]);
            conc0y = (conc0[index4] - conc0[index3])/(2.0*top_h[1]);

            qt1[index0] = -(dens[index0]*cp[index0]*prt[index0]*delta2*s[index0])*dtemx;
            qt2[index0] = -(dens[index0]*cp[index0]*prt[index0]*delta2*s[index0])*dtemy;
	 
	    qp1[index0] = -(dens[index0]*sct[index0]*delta2*s[index0])*concx;
            qp2[index0] = -(dens[index0]*sct[index0]*delta2*s[index0])*concy;
            qp3[index0] = -(dens[index0]*sct0[index0]*delta2*s[index0])*conc0x;
            qp4[index0] = -(dens[index0]*sct0[index0]*delta2*s[index0])*conc0y;


        }

        double *qt, *qc0, *qc;
	FT_VectorMemoryAlloc((POINTER*)&qt,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qc0,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qc,size,sizeof(double));

        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
        {

            index0 = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);

            qt[index0] = (qt1[index2]-qt1[index1])/(2.0*top_h[0])  \
                           + (qt2[index4]-qt2[index3])/(2.0*top_h[1]);
	    qc[index0] = (qp1[index2]-qp1[index1])/(2.0*top_h[0])  \
                           + (qp2[index4]-qp2[index3])/(2.0*top_h[1]);
            qc0[index0] = (qp3[index2]-qp3[index1])/(2.0*top_h[0])  \
                           + (qp4[index4]-qp4[index3])/(2.0*top_h[1]);
		
	    if(subgrid_con == true)
            {
               engy[index0] += -m_dt*qt[index0];
	       
            }

	    if(subgrid_md == true)
            {
               pdens[0][index0] += -m_dt*qc0[index0];
               pdens[1][index0] += -m_dt*qc[index0];
            }

	}

	double taukkr, taukkl, t1, t2;
		
        for (j = imin[1]; j <= imax[1]; j++) 
        for (i = imin[0]; i <= imax[0]; i++) 
        {    
	 
            index0 = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);

            if(subgrid_vis == true)
            {
               taukkr = tau[0][0][index2] + tau[1][1][index4];
               taukkl = tau[0][0][index1] + tau[1][1][index3];
               t1 = (taukkr*u[index2] - taukkl*u[index1])/(2*top_h[0]);
               t2 = (taukkr*v[index4] - taukkl*v[index3])/(2*top_h[1]);

	       // tau is a symmetric tensor
               momn[0][index0] += -m_dt*( (tau[0][0][index2]-tau[0][0][index1])/(2.0*top_h[0])  \
                                          + (tau[0][1][index4]-tau[0][1][index3])/(2.0*top_h[1]) );

               momn[1][index0] += -m_dt*( (tau[0][1][index2]-tau[0][1][index1])/(2.0*top_h[0]) \
                                          + (tau[1][1][index4]-tau[1][1][index3])/(2.0*top_h[1]) );

               engy[index0] += m_dt*0.5*(t1 + t2);
             }

                
	}  

        FT_FreeThese(2, u, v);
        FT_FreeThese(4, s, s11, s12, s22);
        FT_FreeThese(6, tx, ty, cx, cy, cx0, cy0);
        FT_FreeThese(4, cp, temp, conc, conc0);
        FT_FreeThese(5, cs, ci, prt, sct, sct0);
        FT_FreeThese(6, qt1, qt2, qp1, qp2, qp3, qp4);
        FT_FreeThese(1, tau);
        FT_FreeThese(3, qt, qc0, qc);


         
	 
} //end of compSGS2D() 



void G_CARTESIAN::compSGS3D(SWEEP *m_vst)
{

        int i, j, k;
        int index0, index1, index2, index3, index4, index5, index6;
        int index000, index100, index010, index110, index001, index101, index011, \
            index111;
        double *u, *v, *w;
        double ux, uy, uz, vx, vy, vz, wx, wy, wz;
        double *tx, *ty, *tz;
        double *cx, *cy, *cz, *cx0, *cy0, *cz0;
        double *s, *s11, *s12, *s22, *s13, *s23, *s33;
        double *cp, *temp, *conc, *conc0;
        double sum_rho_u, sum_rho_v, sum_rho_w, sum_rho_uu, sum_rho_vv, sum_rho_ww, \
               sum_rho_uv, sum_rho_uw, sum_rho_vw;
        double sum_s11, sum_s12, sum_s22, sum_s13, sum_s23, sum_s33, sum_rhoss11, \
               sum_rhoss12, sum_rhoss22, sum_rhoss13, sum_rhoss23, sum_rhoss, \
               sum_rho, sum_s, sum_ss, sum_cp, sum_tx, sum_ty, sum_tz, sum_rho_t, \
               sum_rho_cp, sum_rho_cp_s_tx, sum_rho_cp_s_ty, sum_rho_cp_s_tz, \
               sum_rho_cp_tu, sum_rho_cp_tv, sum_rho_cp_tw, sum_u;
	double sum_cx, sum_cy, sum_cz, sum_cx0, sum_cy0, sum_cz0, sum_rho_s_cx, \
               sum_rho_s_cy, sum_rho_s_cz, sum_rho_s_cx0, sum_rho_s_cy0, \
               sum_rho_s_cz0, sum_rho_uc, sum_rho_vc, sum_rho_wc, sum_rho_uc0, \
               sum_rho_vc0, sum_rho_wc0, sum_rho_c, sum_rho_c0;
        double MA11, MA12, MA22, MA13, MA23, L11, L12, L22, L13, L23, L33, \
               LA11, LA22, LI, MI, MH1, MH2, MH3, LH1, LH2, LH3;
        double MC1, MC2, MC3, MC01, MC02, MC03, LC1, LC2, LC3, LC01, LC02, LC03;
	double CS, CI, Prt, Sct0, Sct;
        double *cs, *ci, *prt, *sct0, *sct;
        int    ii, jj, kk, iiii, jjjj, kkkk;
        int    NB = 2.0;
        int    NBC = pow(NB, 3);
        double CS_deno, CS_nume, CI_deno, CI_nume, Prt_nume, Prt_deno, Sct0_nume, \
               Sct0_deno, Sct_nume, Sct_deno;
        double delta2, tdelta2;

	delta2 = sqr( pow( top_h[0]*top_h[1]*top_h[2], 1.0/3.0 ) ); //filter width 
        tdelta2 = sqr( pow( NB*top_h[0]*NB*top_h[1]*NB*top_h[2], 1.0/3.0) );  // test filter width

        int size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
        FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s22,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&s13,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&s23,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&s33,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tx,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ty,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&tz,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&cx0,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cy0,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cz0,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&cx,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cy,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cz,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cp,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&temp,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&conc0,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&conc,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cs,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ci,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&prt,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&sct0,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&sct,size,sizeof(double));

        STATE st;
        st.dim = eqn_params->dim;
        EOS_PARAMS *eos = eqn_params->eos;

	double *dens = m_vst->dens;
        double *pres = m_vst->pres;
        double **momn = m_vst->momn;
        double *engy = m_vst->engy;
        double **pdens = m_vst->pdens;

        int comp;
	bool subgrid_con = eqn_params->subgrid_con;
	bool subgrid_md = eqn_params->subgrid_md;
	bool subgrid_vis = eqn_params->subgrid_vis;        

	for(k =imin[2]-1; k <= imax[2]+1; k++)
        for(j =imin[1]-1; j <= imax[1]+1; j++)
        for(i =imin[0]-1; i <= imax[0]+1; i++)
        {
            index0 = d_index3d(i,j,k,top_gmax);
            comp = getComponent(index0);
	    if (comp == 0 && m_vst->comp[index0] != 0)
		comp = m_vst->comp[index0];
	    else if (comp == 0 && m_vst->comp[index0] == 0)
		continue;

            st.eos = &(eos[comp]);
            st.dens = dens[index0];
	    if(eqn_params->multi_comp_non_reactive == YES)
	    {
		int ii;
		for(ii = 0; ii < eqn_params->n_comps; ii++)
		{
		    st.pdens[ii] = pdens[ii][index0];
		}
	    }
            st.engy = engy[index0];
            //st.pres = pres[index0];
            st.momn[0] = momn[0][index0];
            st.momn[1] = momn[1][index0];
            st.momn[2] = momn[2][index0];
            st.pres = EosPressure(&st);

            u[index0] = momn[0][index0]/dens[index0];
            v[index0] = momn[1][index0]/dens[index0];
            w[index0] = momn[2][index0]/dens[index0];
            conc0[index0] = pdens[0][index0]/dens[index0];
            conc[index0]  = pdens[1][index0]/dens[index0];
            cp[index0] = EosCP(&st);
	    temp[index0] = EosTemperature(&st);
        }

	for (k = imin[2]; k <= imax[2]; k++)
        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
        {
	   
            index0 = d_index3d(i,j,k,top_gmax);
            index1 = d_index3d(i-1,j,k,top_gmax);
            index2 = d_index3d(i+1,j,k,top_gmax);
            index3 = d_index3d(i,j-1,k,top_gmax);
            index4 = d_index3d(i,j+1,k,top_gmax);
	    index5 = d_index3d(i,j,k-1,top_gmax);
	    index6 = d_index3d(i,j,k+1,top_gmax);
            
            ux = (u[index2] - u[index1]) / (2.0*top_h[0]);
            uy = (u[index4] - u[index3]) / (2.0*top_h[1]);
	    uz = (u[index6] - u[index5]) / (2.0*top_h[2]);
            vx = (v[index2] - v[index1]) / (2.0*top_h[0]);
            vy = (v[index4] - v[index3]) / (2.0*top_h[1]);
	    vz = (v[index6] - v[index5]) / (2.0*top_h[2]);
	    wx = (w[index2] - w[index1]) / (2.0*top_h[0]);
            wy = (w[index4] - w[index3]) / (2.0*top_h[1]);
            wz = (w[index6] - w[index5]) / (2.0*top_h[2]);

            tx[index0] = (temp[index2] - temp[index1]) / (2.0*top_h[0]);
            ty[index0] = (temp[index4] - temp[index3]) / (2.0*top_h[1]);
	    tz[index0] = (temp[index6] - temp[index5]) / (2.0*top_h[2]);
		
	    cx0[index0] = (conc0[index2] - conc0[index1]) / (2.0*top_h[0]);
            cx[index0] = (conc[index2] - conc[index1]) / (2.0*top_h[0]);
            cy0[index0] = (conc0[index4] - conc0[index3]) / (2.0*top_h[1]);
            cy[index0] = (conc[index4] - conc[index3]) / (2.0*top_h[1]);
	    cz0[index0] = (conc0[index6] - conc0[index5]) / (2.0*top_h[2]);
	    cz[index0] = (conc[index6] - conc[index5]) / (2.0*top_h[2]);

            s11[index0] = ux;
            s12[index0] = 0.5*(uy + vx);
            s22[index0] = vy;
	    s13[index0] = 0.5*(uz + wx);
	    s23[index0] = 0.5*(vz + wy);
            s33[index0] = wz;	
            s[index0] = sqrt( 2*( sqr(s11[index0]) + sqr(s22[index0])  \
                          + sqr(s33[index0]) + 2.0*(sqr(s12[index0])  \
                          + sqr(s13[index0]) + sqr(s23[index0])) ) );

        }

	for (k = 0; k <= ((imax[2]-imin[2]+1)/NB)-1; k++)
        for (j = 0; j <= ((imax[1]-imin[1]+1)/NB)-1; j++)
        for (i = 0; i <= ((imax[0]-imin[0]+1)/NB)-1; i++)
        {
	    kk = (NB*k) + imin[2];
            jj = (NB*j) + imin[1];
            ii = (NB*i) + imin[0];

            sum_rho = 0.0;
            sum_u = 0.0;
            sum_cp = 0.0;
            sum_rho_u = sum_rho_v = sum_rho_w = 0.0;
            sum_rho_uu = sum_rho_vv = sum_rho_ww = sum_rho_uv = sum_rho_uw  \
                       = sum_rho_vw = 0.0;
            sum_s11 = sum_s12 =  sum_s22 = sum_s23 = sum_s13 = sum_s33 = 0.0;
            sum_rhoss11 = sum_rhoss12 = sum_rhoss22 = sum_rhoss13 = sum_rhoss23 \
                        = sum_rhoss = 0.0;
            sum_s = sum_ss = 0.0;
            sum_tx = sum_ty = sum_tz = 0.0;
            sum_rho_t = 0.0;
            sum_rho_cp = sum_rho_cp_s_tx = sum_rho_cp_s_ty = sum_rho_cp_s_tz = 0.0;
            sum_rho_cp_tu = sum_rho_cp_tv = sum_rho_cp_tw = 0.0;
	    
	    sum_cx = sum_cy = sum_cz = sum_cx0 = sum_cy0 = sum_cz0 = 0.0;
            sum_rho_s_cx = sum_rho_s_cy = sum_rho_s_cz = sum_rho_s_cx0 \
                       = sum_rho_s_cy0 = sum_rho_s_cz0 = 0.0;
            sum_rho_uc =  sum_rho_vc = sum_rho_wc = sum_rho_uc0 = sum_rho_vc0 \
                       = sum_rho_wc0 = 0.0;
            sum_rho_c = sum_rho_c0 = 0.0;

            double max_CS = 0.0;
            double max_CI = 0.0;
            double max_Prt = 0.0;
	    double max_Sct = 0.0;
            double max_Sct0 = 0.0;


	    for(kkkk = kk; kkkk < kk+NB; kkkk++)
            for(jjjj = jj; jjjj < jj+NB; jjjj++)
            for(iiii = ii; iiii < ii+NB; iiii++)
            {

                index0 = d_index3d(iiii, jjjj, kkkk, top_gmax);
                sum_u += u[index0];
                sum_rho_u += dens[index0]*u[index0];
                sum_rho_v += dens[index0]*v[index0];
		sum_rho_w += dens[index0]*w[index0];
                sum_rho_uu += dens[index0]*u[index0]*u[index0];
                sum_rho_vv += dens[index0]*v[index0]*v[index0];
                sum_rho_uv += dens[index0]*u[index0]*v[index0];
		sum_rho_uw += dens[index0]*u[index0]*w[index0];
		sum_rho_vw += dens[index0]*v[index0]*w[index0];
		sum_rho_ww += dens[index0]*w[index0]*w[index0];
                sum_s11 += (1.0/3.0)*(2.0*s11[index0]-s22[index0]-s33[index0]); 
                sum_s12 += s12[index0];
		sum_s22 += (1.0/3.0)*(2.0*s22[index0]-s11[index0]-s33[index0]);
		sum_s13 += s13[index0];
		sum_s23 += s23[index0];
		sum_s33 += (1.0/3.0)*(2.0*s33[index0]-s11[index0]-s22[index0]);
                sum_rhoss11 += dens[index0]*s[index0]  \
                                  *(1.0/3.0)*(2.0*s11[index0]-s22[index0]-s33[index0]);
                sum_rhoss12 += dens[index0]*s[index0]*s12[index0];
                sum_rhoss22 += dens[index0]*s[index0]  \
                                  *(1.0/3.0)*(2.0*s22[index0]-s11[index0]-s33[index0]);
                sum_rhoss13 += dens[index0]*s[index0]*s13[index0];
                sum_rhoss23 += dens[index0]*s[index0]*s23[index0];
                sum_rhoss += dens[index0]*s[index0]*s[index0];
                sum_s += s[index0];
                sum_ss += s[index0]*s[index0];
                sum_rho += dens[index0];

                sum_cp += cp[index0];
                sum_tx += tx[index0];
                sum_ty += ty[index0];
		sum_tz += tz[index0];
                sum_rho_t += dens[index0]*temp[index0];
                sum_rho_cp += dens[index0]*cp[index0];
                sum_rho_cp_s_tx += dens[index0]*cp[index0]*s[index0]*tx[index0];
                sum_rho_cp_s_ty += dens[index0]*cp[index0]*s[index0]*ty[index0];
		sum_rho_cp_s_tz += dens[index0]*cp[index0]*s[index0]*tz[index0];
                sum_rho_cp_tu += dens[index0]*cp[index0]*temp[index0]*u[index0];
                sum_rho_cp_tv += dens[index0]*cp[index0]*temp[index0]*v[index0];
		sum_rho_cp_tw += dens[index0]*cp[index0]*temp[index0]*w[index0];

		sum_cx += cx[index0];
                sum_cy += cy[index0];
		sum_cz += cz[index0];
                sum_cx0 += cx0[index0];
                sum_cy0 += cy0[index0];
		sum_cz0 += cz0[index0];
                sum_rho_s_cx += dens[index0]*s[index0]*cx[index0];
                sum_rho_s_cy += dens[index0]*s[index0]*cy[index0];
	        sum_rho_s_cz += dens[index0]*s[index0]*cz[index0];
                sum_rho_s_cx0 += dens[index0]*s[index0]*cx0[index0];
                sum_rho_s_cy0 += dens[index0]*s[index0]*cy0[index0];
		sum_rho_s_cz0 += dens[index0]*s[index0]*cz0[index0];
                sum_rho_uc += dens[index0]*u[index0]*conc[index0];
                sum_rho_vc += dens[index0]*v[index0]*conc[index0];
		sum_rho_wc += dens[index0]*w[index0]*conc[index0];
                sum_rho_uc0 += dens[index0]*u[index0]*conc0[index0];
                sum_rho_vc0 += dens[index0]*v[index0]*conc0[index0];
		sum_rho_wc0 += dens[index0]*w[index0]*conc0[index0];
                sum_rho_c += dens[index0]*conc[index0];
                sum_rho_c0 += dens[index0]*conc0[index0];


            }

	// anistropic and independent part of M
            MA11 = 2.0*delta2*(sum_rhoss11/NBC)  \
                      - 2.0*tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_s11/NBC); 
            MA12 = 2.0*delta2*(sum_rhoss12/NBC)  \
                      - 2.0*tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_s12/NBC);
	    MA22 = 2.0*delta2*(sum_rhoss22/NBC)  \
                      - 2.0*tdelta2*(sum_rho/NBC)*(sum_s/(NBC))*(sum_s22/NBC); 
            MA13 = 2.0*delta2*(sum_rhoss13/NBC)  \
                      - 2.0*tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_s13/NBC);
            MA23 = 2.0*delta2*(sum_rhoss23/NBC)  \
                      - 2.0*tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_s23/NBC);

	// isotropic part of M
            MI = -2.0*delta2*(sum_rhoss/NBC) + 2.0*tdelta2*(sum_rho/NBC)*(sum_ss/NBC);

	// Leonard stress tensor
            L11 = sum_rho_uu/NBC - (sum_rho_u/NBC)*(sum_rho_u/NBC)/(sum_rho/NBC);
            L12 = sum_rho_uv/NBC - (sum_rho_u/NBC)*(sum_rho_v/NBC)/(sum_rho/NBC);
            L22 = sum_rho_vv/NBC - (sum_rho_v/NBC)*(sum_rho_v/NBC)/(sum_rho/NBC);
	    L13 = sum_rho_uw/NBC - (sum_rho_u/NBC)*(sum_rho_w/NBC)/(sum_rho/NBC);
	    L23 = sum_rho_vw/NBC - (sum_rho_v/NBC)*(sum_rho_w/NBC)/(sum_rho/NBC);
	    L33 = sum_rho_ww/NBC - (sum_rho_w/NBC)*(sum_rho_w/NBC)/(sum_rho/NBC);
	
	 // Breaking L into anistropic components(that will be a part of the 
         // five independent relations for C_s) and isotropic part.
            LA11 = (2.0*L11 - L22 - L33)/3.0;
            LA22 = (2.0*L22 - L11 - L22)/3.0;
            // isotropic part 
            LI = (L11 + L22 + L33)/3.0;

            MH1 = delta2*(sum_rho_cp_s_tx/NBC)  \
                    - tdelta2*(sum_rho/NBC)*(sum_cp/NBC)*(sum_s/NBC)*(sum_tx/NBC);
            MH2 = delta2*(sum_rho_cp_s_ty/NBC)  \
                    - tdelta2*(sum_rho/NBC)*(sum_cp/NBC)*(sum_s/NBC)*(sum_ty/NBC);
	    MH3 = delta2*(sum_rho_cp_s_tz/NBC)  \
                    - tdelta2*(sum_rho/NBC)*(sum_cp/NBC)*(sum_s/NBC)*(sum_tz/NBC);

            LH1 = sum_rho_cp_tu/NBC  \
                   - (sum_rho_cp/NBC)*(sum_rho_t/NBC)*(sum_rho_u/NBC)/((sum_rho/NBC)*(sum_rho/NBC));
            LH2 = sum_rho_cp_tv/NBC  \
                   - (sum_rho_cp/NBC)*(sum_rho_t/NBC)*(sum_rho_v/NBC)/((sum_rho/NBC)*(sum_rho/NBC));
	    LH3 = sum_rho_cp_tw/NBC  \
                   - (sum_rho_cp/NBC)*(sum_rho_t/NBC)*(sum_rho_w/NBC)/((sum_rho/NBC)*(sum_rho/NBC));

	    MC1 = delta2*(sum_rho_s_cx/NBC) - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cx/NBC);
            MC2 = delta2*(sum_rho_s_cy/NBC) - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cy/NBC);
	    MC3 = delta2*(sum_rho_s_cz/NBC) - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cz/NBC);
            MC01 = delta2*(sum_rho_s_cx0/NBC) - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cx0/NBC);
            MC02 = delta2*(sum_rho_s_cy0/NBC) - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cy0/NBC);
	    MC03 = delta2*(sum_rho_s_cz0/NBC) - tdelta2*(sum_rho/NBC)*(sum_s/NBC)*(sum_cz0/NBC);

            LC1 = sum_rho_uc/NBC - (sum_rho_u/NBC)*(sum_rho_c/NBC)/(sum_rho/NBC);
            LC2 = sum_rho_vc/NBC - (sum_rho_v/NBC)*(sum_rho_c/NBC)/(sum_rho/NBC);
	    LC3 = sum_rho_wc/NBC - (sum_rho_w/NBC)*(sum_rho_c/NBC)/(sum_rho/NBC);
            LC01 = sum_rho_uc0/NBC - (sum_rho_u/NBC)*(sum_rho_c0/NBC)/(sum_rho/NBC);
            LC02 = sum_rho_vc0/NBC - (sum_rho_v/NBC)*(sum_rho_c0/NBC)/(sum_rho/NBC);
	    LC03 = sum_rho_wc0/NBC - (sum_rho_w/NBC)*(sum_rho_c0/NBC)/(sum_rho/NBC);

            CS_deno = MA11*MA11 + MA12*MA12 + MA22*MA22 + MA13*MA13 + MA23*MA23;
            CS_nume = LA11*MA11 + L12*MA12 + LA22*MA22 + L13*MA13 + L23*MA23;
            CI_deno = MI;
            CI_nume = LI;
            Prt_nume = LH1*MH1 + LH2*MH2 + LH3*MH3; /* Cs is omitted to make the update to the energy equation simpler */
            Prt_deno = MH1*MH1 + MH2*MH2 + MH3*MH3;
	    Sct0_nume  = LC01*MC01 + LC02*MC02 + LC03*MC03;
            Sct0_deno = MC01*MC01 + MC02*MC02 + MC03*MC03;
            Sct_nume = LC1*MC1 + LC2*MC2 + LC3*MC3;
            Sct_deno = MC1*MC1 + MC2*MC2 + MC3*MC3;

            CS = CS_nume/CS_deno;
            CI = CI_nume/CI_deno;
            Prt = Prt_nume/Prt_deno;
	    Sct0 = Sct0_nume/Sct0_deno;
	    Sct = Sct_nume/Sct_deno;

            if(CS < 0.0 || CS_deno < 1e-16)
                CS = 0.0;
            if(CS > max_CS)
                max_CS = CS;
            if(CI < 0.0 || CI_deno < 1e-16)
                CI = 0.0;
            if(CI > max_CI)
                max_CI = CI;
            if(Prt < 0.0 || Prt_deno < 1e-16)
                Prt = 0.0;
            if(Prt > max_Prt)
                max_Prt = Prt;
	    if(Sct0 < 0.0 || Sct0_deno < 1e-16)
                Sct0 = 0.0;
            if(Sct0 > max_Sct0)
                max_Sct0 = Sct0;
	    if(Sct < 0.0 || Sct_deno < 1e-16)
                Sct = 0.0;
            if(Sct > max_Sct)
                max_Sct = Sct;
            
            index000 = d_index3d(ii, jj, kk, top_gmax);
            index100 = d_index3d(ii+1, jj, kk, top_gmax);
            index010 = d_index3d(ii, jj+1, kk, top_gmax);
            index110 = d_index3d(ii+1, jj+1, kk, top_gmax);
	    index001 = d_index3d(ii, jj, kk+1, top_gmax);
            index101 = d_index3d(ii+1, jj, kk+1, top_gmax);
            index011 = d_index3d(ii, jj+1, kk+1, top_gmax);
            index111 = d_index3d(ii+1, jj+1, kk+1, top_gmax);

            cs[index000] = CS;
            cs[index100] = CS;
            cs[index010] = CS;
            cs[index110] = CS;
	    cs[index001] = CS;
            cs[index101] = CS;
            cs[index011] = CS;
            cs[index111] = CS;
	
	    ci[index000] = CI;
            ci[index100] = CI;
            ci[index010] = CI;
            ci[index110] = CI;
            ci[index001] = CI;
            ci[index101] = CI;
            ci[index011] = CI;
            ci[index111] = CI;

	    prt[index000] = Prt;
            prt[index100] = Prt;
            prt[index010] = Prt;
            prt[index110] = Prt;
            prt[index001] = Prt;
            prt[index101] = Prt;
            prt[index011] = Prt;
            prt[index111] = Prt;

	    sct0[index000] = Sct0;
            sct0[index100] = Sct0;
            sct0[index010] = Sct0;
            sct0[index110] = Sct0;
            sct0[index001] = Sct0;
            sct0[index101] = Sct0;
            sct0[index011] = Sct0;
            sct0[index111] = Sct0;

            sct[index000] = Sct;
            sct[index100] = Sct;
            sct[index010] = Sct;
            sct[index110] = Sct;
            sct[index001] = Sct;
            sct[index101] = Sct;
            sct[index011] = Sct;
            sct[index111] = Sct;

        }

        double ***tau;
        double *qt1, *qt2, *qt3;
	double  *qp1, *qp2, *qp3, *qp4, *qp5, *qp6;
        double dtemx, dtemy, dtemz;
	double concx, concy, concz, conc0x, conc0y, conc0z;

        FT_VectorMemoryAlloc((POINTER*)&qt1,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qt2,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&qt3,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&qp1,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qp2,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qp3,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&qp4,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&qp5,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&qp6,size,sizeof(double));
	FT_TriArrayMemoryAlloc((POINTER*)&tau,dim,dim,size,sizeof(double));


	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
	{

            index0  = d_index3d(i,j,k,top_gmax);
            index1 = d_index3d(i-1,j,k,top_gmax);
            index2 = d_index3d(i+1,j,k,top_gmax);
            index3 = d_index3d(i,j-1,k,top_gmax);
            index4 = d_index3d(i,j+1,k,top_gmax);
	    index5 = d_index3d(i,j,k-1,top_gmax);
	    index6 = d_index3d(i,j,k+1,top_gmax);
 
            tau[0][0][index0] = -2.0*cs[index0]*delta2*dens[index0]*s[index0]*(2.0*s11[index0]-s22[index0]- s22[index0])/3.0 + (2.0/3.0)*ci[index0]*delta2*dens[index0]*s[index0]*s[index0];
            tau[0][1][index0] = -2.0*cs[index0]*delta2*dens[index0]*s[index0]*(s12[index0]);
            tau[1][1][index0] = -2.0*cs[index0]*delta2*dens[index0]*s[index0]*(2.0*s22[index0]-s11[index0]- s33[index0])/3.0 + (2.0/3.0)*ci[index0]*delta2*dens[index0]*s[index0]*s[index0];
	    tau[0][2][index0] = -2.0*cs[index0]*delta2*dens[index0]*s[index0]*(s13[index0]);
	    tau[1][2][index0] = -2.0*cs[index0]*delta2*dens[index0]*s[index0]*(s23[index0]);
            tau[2][2][index0] = -2.0*cs[index0]*delta2*dens[index0]*s[index0]*(2.0*s33[index0]-s11[index0]- s22[index0])/3.0 + (2.0/3.0)*ci[index0]*delta2*dens[index0]*s[index0]*s[index0];

            dtemx = (temp[index2]-temp[index1])/(2.0*top_h[0]);
            dtemy = (temp[index4]-temp[index3])/(2.0*top_h[1]);
	    dtemy = (temp[index6]-temp[index5])/(2.0*top_h[2]);
	
	    concx = (conc[index2] - conc[index1])/(2.0*top_h[0]);
            concy = (conc[index4] - conc[index3])/(2.0*top_h[1]);
	    concz = (conc[index6] - conc[index5])/(2.0*top_h[2]);
            conc0x = (conc0[index2] - conc0[index1])/(2.0*top_h[0]);
            conc0y = (conc0[index4] - conc0[index3])/(2.0*top_h[1]); 
	    conc0z = (conc[index6] - conc[index5])/(2.0*top_h[2]);

            qt1[index0] = -(dens[index0]*cp[index0]*prt[index0]*delta2*s[index0])*dtemx;
            qt2[index0] = -(dens[index0]*cp[index0]*prt[index0]*delta2*s[index0])*dtemy;
	    qt3[index0] = -(dens[index0]*cp[index0]*prt[index0]*delta2*s[index0])*dtemz;

	    qp1[index0] = -(dens[index0]*sct0[index0]*delta2*s[index0])*conc0x;
            qp2[index0] = -(dens[index0]*sct0[index0]*delta2*s[index0])*conc0y;
	    qp3[index0] = -(dens[index0]*sct0[index0]*delta2*s[index0])*conc0z;
            qp4[index0] = -(dens[index0]*sct[index0]*delta2*s[index0])*concx;
            qp5[index0] = -(dens[index0]*sct[index0]*delta2*s[index0])*concy;
	    qp6[index0] = -(dens[index0]*sct[index0]*delta2*s[index0])*concz;

        }

        double *qt, *qc0, *qc;
	FT_VectorMemoryAlloc((POINTER*)&qt,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qc0,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&qc,size,sizeof(double));


	for (k = imin[2]; k <= imax[2]; k++)
        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
        {

            index0 = d_index3d(i,j,k,top_gmax);
            index1 = d_index3d(i-1,j,k,top_gmax);
            index2 = d_index3d(i+1,j,k,top_gmax);
            index3 = d_index3d(i,j-1,k,top_gmax);
            index4 = d_index3d(i,j+1,k,top_gmax);
	    index5 = d_index3d(i,j,k-1,top_gmax);
	    index6 = d_index3d(i,j,k+1,top_gmax);

            qt[index0] = (qt1[index2]-qt1[index1])/(2.0*top_h[0])  \
                           + (qt2[index4]-qt2[index3])/(2.0*top_h[1])  \
                           + (qt3[index6]-qt3[index5])/(2.0*top_h[2]);
	
	    qc0[index0] = (qp1[index2]-qp1[index1])/(2.0*top_h[0])  \
                            + (qp2[index4]-qp2[index3])/(2.0*top_h[1])  \
                            + (qp3[index6]-qp3[index5])/(2.0*top_h[2]);
            qc[index0] = (qp4[index2]-qp4[index1])/(2.0*top_h[0])  \
                            + (qp5[index4]-qp5[index3])/(2.0*top_h[1])  \
                            + (qp6[index6]-qp6[index5])/(2.0*top_h[2]);

		
	    if(subgrid_con == true)
            {
               engy[index0] += -m_dt*qt[index0];
            }

	   if(subgrid_md == true)
            {
               pdens[0][index0] += -m_dt*qc0[index0];
               pdens[1][index0] += -m_dt*qc[index0];
            }


	}

	double t1, t2, t3;
	double t1r, t1l, t2r, t2l, t3r, t3l;
		
	for (k = imin[2]; k <= imax[2]; k++)
        for (j = imin[1]; j <= imax[1]; j++) 
        for (i = imin[0]; i <= imax[0]; i++) 
        {    
	 
	    index0 = d_index3d(i,j,k,top_gmax);
            index1 = d_index3d(i-1,j,k,top_gmax);
            index2 = d_index3d(i+1,j,k,top_gmax);
            index3 = d_index3d(i,j-1,k,top_gmax);
            index4 = d_index3d(i,j+1,k,top_gmax);
            index5 = d_index3d(i,j,k-1,top_gmax);
            index6 = d_index3d(i,j,k+1,top_gmax);
		

            if(subgrid_vis == true)
            {
		t1r = (tau[0][0][index2]+tau[1][1][index2]+tau[2][2][index2])*u[index2];
		t1l = (tau[0][0][index1]+tau[1][1][index1]+tau[2][2][index1])*u[index1];
		t1 = (t1r-t1l)/(2.0*top_h[0]);

		t2r = (tau[0][0][index4]+tau[1][1][index4]+tau[2][2][index4])*v[index4];
                t2l = (tau[0][0][index3]+tau[1][1][index3]+tau[2][2][index3])*v[index3];
		t2 = (t2r-t2l)/(2.0*top_h[1]);

		t3r = (tau[0][0][index6]+tau[1][1][index6]+tau[2][2][index6])*w[index6];
                t3l = (tau[0][0][index5]+tau[1][1][index5]+tau[2][2][index5])*w[index5];
		t3 = (t3r-t3l)/(2.0*top_h[2]);

               // tau is a symmetric tensor
               momn[0][index0] += -m_dt*( (tau[0][0][index2]-tau[0][0][index1])/(2.0*top_h[0]) \
                                             + (tau[0][1][index4]-tau[0][1][index3])/(2.0*top_h[1]) \
                                             + (tau[2][0][index6]-tau[2][0][index5])/(2.0*top_h[2]) );

               momn[1][index0] += -m_dt*( (tau[0][1][index2]-tau[0][1][index1])/(2.0*top_h[0]) \
                                             + (tau[1][1][index4]-tau[1][1][index3])/(2.0*top_h[1]) \
                                             + (tau[1][2][index6]-tau[1][2][index5])/(2.0*top_h[2]) );

	       momn[2][index0] += -m_dt*( (tau[0][2][index2]-tau[0][2][index1])/(2.0*top_h[0]) \
                                             + (tau[1][2][index4]-tau[1][2][index3])/(2.0*top_h[1]) \
                                             + (tau[2][2][index6]-tau[2][2][index5])/(2.0*top_h[2]) );

               engy[index0] += m_dt*0.5*(t1 + t2 + t3);
             }

                
	}  

	 
        FT_FreeThese(3, u, v, w);
        FT_FreeThese(7, s, s11, s12, s22, s13, s23, s33);
        FT_FreeThese(9, tx, ty, tz, cx0, cy0, cz0, cx, cy, cz);
        FT_FreeThese(4, cp, temp, conc0, conc);
        FT_FreeThese(5, cs, ci, prt, sct0, sct);
        FT_FreeThese(9, qt1, qt2, qt3, qp1, qp2, qp3, qp4, qp5, qp6);
        FT_FreeThese(3, qt, qc0, qc);
        FT_FreeThese(1, tau); 

         
	 
} //end of computeSGS3D() 


void G_CARTESIAN::computeParab(void)
{
	static SWEEP st_field;
	static boolean first = YES;

	if (first)
	{
	    allocMeshVst(&st_field);
	    first = NO;
	}

	copyToMeshVst(&st_field);

	if (dim == 2)
	    parab2D(&st_field);
	else if (dim == 3)
	    parab3D(&st_field);
	else
	{
	    printf("Parabolic step is not implemented for %d dimension.\n", dim);
	    exit(0);
	}

	scatMeshVst(&st_field);
	copyFromMeshVst(st_field);
}

void G_CARTESIAN::parab2D(SWEEP *m_vst) 
{
	int size;
	int i, j;
	int index0, index1, index2, index3, index4;
	double *u, *v;
        double *temp, *conc0, *conc;
		
        size = (top_gmax[0]+1)*(top_gmax[1]+1);
        FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&temp,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&conc0,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&conc,size,sizeof(double));

        STATE st;
        st.dim = eqn_params->dim;
        EOS_PARAMS *eos = eqn_params->eos;
	double *dens = m_vst->dens;
        double *pres = m_vst->pres;
        double **momn = m_vst->momn;
        double *engy = m_vst->engy;
        double **pdens = m_vst->pdens;
	
        int comp;
       
	bool parabolic = eqn_params->parabolic_step;
        double mu = eqn_params->mu;
        double kappa = eqn_params->kappa;
	double D = eqn_params->D;
        bool viscosity = eqn_params->viscosity;
        bool thermal_conduction  = eqn_params->thermal_conduction;
	bool mass_diffusion  = eqn_params->mass_diffusion;
        bool subgrid = eqn_params->subgrid;		

	appendGhostBufferforParab2D(m_vst);

        for(j = imin[1]-1; j <= imax[1]+1; j++)
        for(i = imin[0]-1; i <= imax[0]+1; i++)
        {
            index0 = d_index2d(i,j,top_gmax);
            comp = getComponent(index0);
	    if (comp == 0 && m_vst->comp[index0] != 0)
		comp = m_vst->comp[index0];
	    else if (comp == 0 && m_vst->comp[index0] == 0)
		continue;

            st.eos = &(eos[comp]);
            st.dens = dens[index0];
	    if(eqn_params->multi_comp_non_reactive == YES)
	    {
		int ii;
		for(ii = 0; ii < eqn_params->n_comps; ii++)
		{
		    st.pdens[ii] = pdens[ii][index0];
		}
	    }
            st.engy = engy[index0];
            st.momn[0] = momn[0][index0];
            st.momn[1] = momn[1][index0];
            st.pres = EosPressure(&st);

            u[index0] = momn[0][index0]/dens[index0];
            v[index0] = momn[1][index0]/dens[index0];
            conc0[index0] = pdens[0][index0]/dens[index0];
            conc[index0]  = pdens[1][index0]/dens[index0];
            temp[index0] = EosTemperature(&st);
	}

	if (subgrid == true)
	{ 
           compSGS2D(m_vst);
 	}
	
	double rhox, rhoy;
	double concx, concy, conc0x, conc0y, concxx, concyy, conc0xx, conc0yy;
	double txx, tyy;
        double ux, uy, vx, vy, uxx, uyy, vxx, vyy;
 
        for (j = imin[1]; j <= imax[1]; j++) 
        for (i = imin[0]; i <= imax[0]; i++) 
        {    

            index0 = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);

	 
           if (viscosity == true)
            {
		ux = ( u[index2] - u[index1] ) / (2.0*top_h[0]);
		uy = ( u[index4] - u[index3] ) / (2.0*top_h[1]);
		vx = ( v[index2] - v[index1] ) / (2.0*top_h[0]);
		vy = ( v[index4] - v[index3] ) / (2.0*top_h[1]);
		uxx = ( u[index2] - 2.0*u[index0] + u[index1]) / (top_h[0]*top_h[0]);
		uyy = ( u[index4] - 2.0*u[index0] + u[index3]) / (top_h[1]*top_h[1]);
		vxx = ( v[index2] - 2.0*v[index0] + v[index1]) / (top_h[0]*top_h[0]);
		vyy = ( v[index4] - 2.0*v[index0] + v[index3]) / (top_h[1]*top_h[1]);

                momn[0][index0] += m_dt*mu*(uxx+uyy);
                momn[1][index0] += m_dt*mu*(vxx+vyy);
                
		engy[index0] += m_dt*mu*( u[index0]*(uxx+uyy)  \
                                    + v[index0]*(vxx+vyy) + sqr(ux-vy) + sqr(vx+uy) );
	
            }

           if (thermal_conduction == true)
            {
	        txx = (temp[index2] - 2*temp[index0] + temp[index1]) / (top_h[0]*top_h[0]);
                tyy = (temp[index4] - 2*temp[index0] + temp[index3]) / (top_h[1]*top_h[1]);
                engy[index0] += m_dt* (kappa*(txx + tyy)); 
            }

	  if (mass_diffusion == true)
            {
                rhox = (dens[index2] - dens[index1]) / (2.0*top_h[0]);
                rhoy = (dens[index4] - dens[index3]) / (2.0*top_h[1]);
                concx = (conc[index2] - conc[index1]) / (2.0*top_h[0]);
                concy = (conc[index4] - conc[index3]) / (2.0*top_h[1]);
                conc0x = (conc0[index2] - conc0[index1]) / (2.0*top_h[0]);
                conc0y = (conc0[index4] - conc0[index3]) / (2.0*top_h[1]);
                concxx = (conc[index2] - 2.0*conc[index0] + conc[index1]) / (top_h[0]*top_h[0]);
                concyy = (conc[index4] - 2.0*conc[index0] + conc[index3]) / (top_h[1]*top_h[1]);
                conc0xx = (conc0[index2] - 2.0*conc0[index0] + conc0[index1]) / (top_h[0]*top_h[0]);
                conc0yy = (conc0[index4] - 2.0*conc0[index0] + conc0[index3]) / (top_h[1]*top_h[1]);
                pdens[0][index0] += m_dt*D* (rhox*conc0x + dens[index0]*conc0xx + rhoy*conc0y + dens[index0]*conc0yy);
                pdens[1][index0] += m_dt*D* (rhox*concx + dens[index0]*concxx + rhoy*concy + dens[index0]*concyy);
             }

       	}

        FT_FreeThese(5, u, v, temp, conc0, conc);

}

void G_CARTESIAN::parab3D(SWEEP *m_vst)
{
        int size;
        int i, j, k;
        int index0, index1, index2, index3, index4, index5, index6;
        double *u, *v, *w;
        double *temp, *conc0, *conc, *HDiff;
	double conc0x, conc0y, conc0z, concx, concy, concz;

        size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
        FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&temp,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&conc0,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&conc,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&HDiff,size,sizeof(double));


        STATE st;
        st.dim = eqn_params->dim;
        EOS_PARAMS *eos = eqn_params->eos;

	double *dens = m_vst->dens;
        double *pres = m_vst->pres;
        double **momn = m_vst->momn;
        double *engy = m_vst->engy;
        double **pdens = m_vst->pdens;

        int comp;
        bool parabolic = eqn_params->parabolic_step;
	double mu = eqn_params->mu;
        double kappa = eqn_params->kappa;
	double D = eqn_params->D;
        bool viscosity = eqn_params->viscosity;
        bool thermal_conduction  = eqn_params->thermal_conduction;
	bool mass_diffusion = eqn_params->mass_diffusion;
        bool subgrid = eqn_params->subgrid;

	appendGhostBufferforParab3D(m_vst);	

	for(k = imin[2]-1; k <= imax[2]+1; k++)
        for(j = imin[1]-1; j <= imax[1]+1; j++)
        for(i = imin[0]-1; i <= imax[0]+1; i++)
        {
            index0 = d_index3d(i,j,k,top_gmax);
            comp = getComponent(index0);
	    if (comp == 0 && m_vst->comp[index0] != 0)
		comp = m_vst->comp[index0];
	    else if (comp == 0 && m_vst->comp[index0] == 0)
		continue;

            st.eos = &(eos[comp]);
            st.dens = dens[index0];
	    if(eqn_params->multi_comp_non_reactive == YES)
	    {
		int ii;
		for(ii = 0; ii < eqn_params->n_comps; ii++)
		{
		    st.pdens[ii] = pdens[ii][index0];
		}
	    }
            st.engy = engy[index0];
            //st.pres = pres[index0];
            st.momn[0] = momn[0][index0];
            st.momn[1] = momn[1][index0];
            st.momn[2] = momn[2][index0];
            st.pres = EosPressure(&st);

            u[index0] = momn[0][index0]/dens[index0];
            v[index0] = momn[1][index0]/dens[index0];
            w[index0] = momn[2][index0]/dens[index0];
            conc0[index0] = pdens[0][index0]/dens[index0];
            conc[index0]  = pdens[1][index0]/dens[index0];
            temp[index0] = EosTemperature(&st);
	    HDiff[index0] = EosEnthalpyDifference(&st);
	}

	if (subgrid == true)
        {
           compSGS3D(m_vst);
        }

        double rhox, rhoy, rhoz;
        double txx, tyy, tzz;
        double ux, uy, uz, vx, vy, vz, wx, wy, wz;
        double uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz;
	double vyx, wxz, uyx, wyz, vzy, uxz;
	double term1, term2, term3, term4, term5;
	int index_r12, index_r13, index_r23, index_l12, index_l13, index_l23;
	
	double *d1, *d2, *d3, *d4, *d5, *d6;
	double *e1, *e2, *e3;
	FT_VectorMemoryAlloc((POINTER*)&d1,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&d2,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&d3,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&d4,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&d5,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&d6,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&e1,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&e2,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&e3,size,sizeof(double));
	
	
	for (k = imin[2]; k <= imax[2]; k++)
        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
        {

            index0 = d_index3d(i,j,k,top_gmax);
            index1 = d_index3d(i-1,j,k,top_gmax);
            index2 = d_index3d(i+1,j,k,top_gmax);
            index3 = d_index3d(i,j-1,k,top_gmax);
            index4 = d_index3d(i,j+1,k,top_gmax);
            index5 = d_index3d(i,j,k-1,top_gmax);
            index6 = d_index3d(i,j,k+1,top_gmax);

	    conc0x = (conc0[index2] - conc0[index1])/(2.0*top_h[0]);
            conc0y = (conc0[index4] - conc0[index3])/(2.0*top_h[1]);
            conc0z = (conc0[index6] - conc0[index5])/(2.0*top_h[2]);
	    concx = (conc[index2] - conc[index1])/(2.0*top_h[0]);
            concy = (conc[index4] - conc[index3])/(2.0*top_h[1]);
            concz = (conc[index6] - conc[index5])/(2.0*top_h[2]);

	    d1[index0] = dens[index0]*D*conc0x;
	    d2[index0] = dens[index0]*D*conc0y;
	    d3[index0] = dens[index0]*D*conc0z;
	    d4[index0] = dens[index0]*D*concx;
            d5[index0] = dens[index0]*D*concy;
            d6[index0] = dens[index0]*D*concz;
	    
	    e1[index0] = HDiff[index0]*dens[index0]*D*conc0x;
	    e2[index0] = HDiff[index0]*dens[index0]*D*conc0y;
	    e3[index0] = HDiff[index0]*dens[index0]*D*conc0z;
	}

 	for (k = imin[2]; k <= imax[2]; k++)
        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
        {

            index0 = d_index3d(i,j,k,top_gmax);
            index1 = d_index3d(i-1,j,k,top_gmax);
            index2 = d_index3d(i+1,j,k,top_gmax);
            index3 = d_index3d(i,j-1,k,top_gmax);
            index4 = d_index3d(i,j+1,k,top_gmax);
	    index5 = d_index3d(i,j,k-1,top_gmax);
            index6 = d_index3d(i,j,k+1,top_gmax);

	   // Indexing for mixed derivatives

	    index_r12 = d_index3d(i+1,j+1,k,top_gmax);
	    index_r13 = d_index3d(i+1,j,k+1,top_gmax);
	    index_r23 = d_index3d(i,j+1,k+1,top_gmax);
	    index_l12 = d_index3d(i-1,j-1,k,top_gmax);
	    index_l13 = d_index3d(i-1,j,k-1,top_gmax);
	    index_l23 = d_index3d(i,j-1,k-1,top_gmax);

           if (viscosity == true)
            {
                ux = ( u[index2] - u[index1] ) / (2.0*top_h[0]);
                uy = ( u[index4] - u[index3] ) / (2.0*top_h[1]);
		uz = ( u[index6] - u[index5] ) / (2.0*top_h[2]);
                vx = ( v[index2] - v[index1] ) / (2.0*top_h[0]);
                vy = ( v[index4] - v[index3] ) / (2.0*top_h[1]);
		vz = ( v[index6] - v[index5] ) / (2.0*top_h[2]);
		wx = ( w[index2] - w[index1] ) / (2.0*top_h[0]);
                wy = ( w[index4] - w[index3] ) / (2.0*top_h[1]);
                wz = ( w[index6] - w[index5] ) / (2.0*top_h[2]);
                uxx = ( u[index2] - 2.0*u[index0] + u[index1]) / (top_h[0]*top_h[0]);
                uyy = ( u[index4] - 2.0*u[index0] + u[index3]) / (top_h[1]*top_h[1]);
		uzz = ( u[index6] - 2.0*u[index0] + u[index5]) / (top_h[2]*top_h[2]);
                vxx = ( v[index2] - 2.0*v[index0] + v[index1]) / (top_h[0]*top_h[0]);
                vyy = ( v[index4] - 2.0*v[index0] + v[index3]) / (top_h[1]*top_h[1]);
		vzz = ( v[index6] - 2.0*v[index0] + v[index5]) / (top_h[2]*top_h[2]);
                wxx = ( w[index2] - 2.0*w[index0] + w[index1]) / (top_h[0]*top_h[0]);
                wyy = ( w[index4] - 2.0*w[index0] + w[index3]) / (top_h[1]*top_h[1]);
                wzz = ( w[index6] - 2.0*w[index0] + w[index5]) / (top_h[2]*top_h[2]);

		// Mixed derivatives
		//vyx, wxz, uyx, wyz, vzy, uxz
	
		uyx = (u[index_r12]-u[index4]-u[index2]+u[index_l12])/(4.0*top_h[1]*top_h[0]);
		uxz = (u[index_r13]-u[index2]-u[index6]+u[index_l13])/(4.0*top_h[0]*top_h[2]);
		vyx = (v[index_r12]-v[index4]-v[index2]+v[index_l12])/(4.0*top_h[1]*top_h[0]);
		vzy = (v[index_r23]-v[index6]-v[index4]+v[index_l23])/(4.0*top_h[2]*top_h[1]);
		wxz = (w[index_r13]-w[index2]-w[index6]+w[index_l13])/(4.0*top_h[0]*top_h[2]);
		wyz = (w[index_r23]-w[index4]-w[index6]+w[index_l23])/(4.0*top_h[1]*top_h[2]);

                momn[0][index0] += m_dt*mu*(uxx+uyy+uzz);
                momn[1][index0] += m_dt*mu*(vxx+vyy+vzz);
		momn[2][index0] += m_dt*mu*(wxx+wyy+wzz);

                //engy[index0] += m_dt*mu*(u[index0]*(uxx+uyy) + v[index0]*(vxx+vyy) + sqr(ux-vy) + sqr(vx+uy) );
		term1 = u[index0]* (4.0/3.0 *uxx + 1.0/3.0 * vyx + 1.0/3.0 * wxz + uyy + uzz); 
		term2 = v[index0]* (4.0/3.0 *vyy + 1.0/3.0 * uyx + 1.0/3.0 * wyz + vxx + vzz); 
		term3 = w[index0]* (4.0/3.0 *wzz + 1.0/3.0 * vzy + 1.0/3.0 * uxz + wxx + wyy);
		term4 = 4.0/3.0 * (ux*ux+vy*vy+wz*wz) + uy*uy + uz*uz + vx*vx + vz*vz + wx*wx + wy*wy;
		term5 = 2.0 * (wx*uz + vz*wy + vx*uy) - 4.0/3.0 * (ux*vy + vy*wz + ux*wz);
		engy[index0] += m_dt*mu*(term1 + term2 + term3 + term4 + term5); 
            }

	 if (thermal_conduction == true)
            {
                txx = (temp[index2] - 2*temp[index0] + temp[index1]) / (top_h[0]*top_h[0]);
                tyy = (temp[index4] - 2*temp[index0] + temp[index3]) / (top_h[1]*top_h[1]);
		tzz = (temp[index6] - 2*temp[index0] + temp[index5]) / (top_h[2]*top_h[2]);
                engy[index0] += m_dt* (kappa*(txx + tyy + tzz));

            }

	  double dc0, dc, ed;

	  if (mass_diffusion == true)
            {
		dc0 = (d1[index2]-d1[index1])/(2*top_h[0])   \
                       + (d2[index4]-d2[index3])/(2*top_h[1])  \
                        + (d3[index6]-d3[index5])/(2.0*top_h[2]);
		dc = (d4[index2]-d4[index1])/(2*top_h[0])   \
                       + (d5[index4]-d5[index3])/(2*top_h[1]) \
                        + (d6[index6]-d6[index5])/(2.0*top_h[2]);    
         	pdens[0][index0] += m_dt*dc0;
		pdens[1][index0] += m_dt*dc; 
                //fprintf(stdout, "PRINT pdens0 pdens1 %e %e\n",pdens[0][index0],pdens[1][index0]);

		ed = (e1[index2]-e1[index1])/(2*top_h[0])   \
                      + (e2[index4]-e2[index3])/(2*top_h[1])   \
                      + (e3[index6]-e3[index5])/(2.0*top_h[2]);
		engy[index0] += m_dt*ed;
	    }
        }

	FT_FreeThese(7, u, v, w, temp, conc, conc0, HDiff); 
        FT_FreeThese(9, d1, d2, d3, d4, d5, d6, e1, e2, e3);
}

void G_CARTESIAN::appendGhostBufferforParab2D(
	SWEEP	*m_vst)
{
    	int		i, j;
	int		ic[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	COMPONENT	bcomp;
	STATE		*state;
	HYPER_SURF	*hs;
	double		crx_coords[MAXD];
	int		index, bindex;

	for (j = imin[1]; j <= imax[1]; j++)
	{
	    ic[1] = j;
	    ic[0] = imin[0]-1;
	    index = d_index(ic,top_gmax,dim);
	    ic[0] = imin[0];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;
	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,ldir[0],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab2D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
		    setNeumannStatesforParab(m_vst,hs,state,ic,0,0,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab2D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }

	    ic[0] = imax[0]+1;
	    index = d_index(ic,top_gmax,dim);
	    ic[0] = imax[0];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;
	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,rdir[0],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab2D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
		    setNeumannStatesforParab(m_vst,hs,state,ic,0,1,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab2D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }
	}

	for (i = imin[0]; i <= imax[0]; i++)
	{
	    ic[0] = i;
	    ic[1] = imin[1]-1;
	    index = d_index(ic,top_gmax,dim);
	    ic[1] = imin[1];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;
	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,ldir[1],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab2D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
//		    setNeumannStatesForParab(ic,m_vst,1,0);
		    setNeumannStatesforParab(m_vst,hs,state,ic,1,0,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab2D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }

	    ic[1] = imax[1]+1;
	    index = d_index(ic,top_gmax,dim);
	    ic[1] = imax[1];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;
	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,rdir[1],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab2D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
		    setNeumannStatesforParab(m_vst,hs,state,ic,1,1,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab2D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }
	}
}

void G_CARTESIAN::appendGhostBufferforParab3D(
	SWEEP	*m_vst)
{
    	int		i, j, k;
	int		ic[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	COMPONENT	bcomp;
	STATE		*state;
	HYPER_SURF	*hs;
	double		crx_coords[MAXD];
	int		index, bindex;

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	{
	    ic[1] = j;
	    ic[2] = k;
	    ic[0] = imin[0]-1;
	    index = d_index(ic,top_gmax,dim);
	    ic[0] = imin[0];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;
	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,ldir[0],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab3D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
		    setNeumannStatesforParab(m_vst,hs,state,ic,0,0,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab3D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }
	}

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	{
	    ic[1] = j;
	    ic[2] = k;
	    ic[0] = imax[0]+1;
	    index = d_index(ic,top_gmax,dim);
	    ic[0] = imax[0];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;
	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,rdir[0],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab3D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
		    setNeumannStatesforParab(m_vst,hs,state,ic,0,1,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab3D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }
	}

	for (k = imin[2]; k <= imax[2]; k++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    ic[0] = i;
	    ic[2] = k;
	    ic[1] = imin[1]-1;
	    index = d_index(ic,top_gmax,dim);
	    ic[1] = imin[1];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;
	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,ldir[1],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab3D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
		    setNeumannStatesforParab(m_vst,hs,state,ic,1,0,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab3D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }
	}

	for (k = imin[2]; k <= imax[2]; k++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    ic[0] = i;
	    ic[2] = k;
	    ic[1] = imax[1]+1;
	    index = d_index(ic,top_gmax,dim);
	    ic[1] = imax[1];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;
	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,rdir[1],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab3D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
		    setNeumannStatesforParab(m_vst,hs,state,ic,1,1,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab3D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }
	}

	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    ic[0] = i;
	    ic[1] = j;
	    ic[2] = imin[2]-1;
	    index = d_index(ic,top_gmax,dim);
	    ic[2] = imin[2];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;

	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,ldir[2],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab3D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
		    setNeumannStatesforParab(m_vst,hs,state,ic,2,0,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab3D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }
	}

	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    ic[0] = i;
	    ic[1] = j;
	    ic[2] = imax[2]+1;
	    index = d_index(ic,top_gmax,dim);
	    ic[2] = imax[2];
	    bindex = d_index(ic,top_gmax,dim);
	    bcomp = cell_center[bindex].comp;
	    if (!needBufferFromIntfc(bcomp,cell_center[index].comp))
		continue;
	    if (!FT_StateStructAtGridCrossing(front,ic,rdir[2],
		bcomp,(POINTER*)&state,&hs,crx_coords))
	    {
		(void) printf("In parab3D():\n");
		(void) printf("ERROR: No crossing found!\n");
		(void) print_int_vector("icoords=", ic,3,"\n");
		clean_up(ERROR);
	    }
	    state->eos = &(eqn_params->eos[bcomp]);
	    switch (wave_type(hs))
	    {
		case NEUMANN_BOUNDARY:
		case MOVABLE_BODY_BOUNDARY:
		    setNeumannStatesforParab(m_vst,hs,state,ic,2,1,bcomp);
		    break;
		case DIRICHLET_BOUNDARY:
		    printf("Dirichlet boundary is not implemented in parab3D().\n");
		    exit(0);
		    break;
		default:
		    break;
	    }
	}
}

void G_CARTESIAN::setNeumannStatesforParab(
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	COMPONENT	comp)
{
	int 		i,j,index;
	int 		ic[MAXD];
	double		*vel_ref = state->vel;
	double		coords[MAXD],coords_ref[MAXD],crx_coords[MAXD];
	double		nor[MAXD],vn,v[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	GRID_DIRECTION  dir;
	STATE		st_tmp;

	st_tmp.eos = state->eos;
	st_tmp.dim = dim;
	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic[i] = icoords[i];
	}
	dir = (nb == 0) ? ldir[idir] : rdir[idir];
	FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);

	/* Find ghost point */
	ic[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;

	for (j = 0; j < dim; ++j)
	    coords_ref[j] = top_L[j] + ic[j]*top_h[j];

	/* Reflect ghost point through intfc-mirror at crossing */
	coords_ref[idir] = 2.0*crx_coords[idir] - coords_ref[idir];
	vn = 0.0;
	for (j = 0; j < dim; ++j)
	{
	    v[j] = coords_ref[j] - crx_coords[j];
	    vn += v[j]*nor[j];
	}
	for (j = 0; j < dim; ++j)
	    v[j] = 2.0*vn*nor[j] - v[j];
	for (j = 0; j < dim; ++j)
	    coords_ref[j] = crx_coords[j] + v[j];
			
	/* Interpolate the state at the reflected point */
	FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->dens,getStateDens,&st_tmp.dens,&m_vst->dens[index]);

	if(eqn_params->multi_comp_non_reactive == YES)
	{
	    int ii;
	    for(ii = 0; ii < eqn_params->n_comps; ii++)
		FT_IntrpStateVarAtCoords(front,comp,coords_ref,m_vst->pdens[ii],
			getStatePdens0,&st_tmp.pdens[ii],&m_vst->pdens[ii][index]);
	}

	FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->engy,getStateEngy,&st_tmp.engy,&m_vst->engy[index]);
	FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->pres,getStatePres,&st_tmp.pres,&m_vst->pres[index]);
	FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->momn[0],getStateXmom,&st_tmp.momn[0],&m_vst->momn[0][index]);
	if (dim > 1)
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->momn[1],getStateYmom,&st_tmp.momn[1],&m_vst->momn[1][index]);
	if (dim > 2)
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->momn[2],getStateZmom,&st_tmp.momn[2],&m_vst->momn[2][index]);

	/* Galileo Transformation */
	vn = 0.0;
	for (j = 0; j < dim; j++)
	{
	    v[j] = st_tmp.momn[j]/st_tmp.dens - vel_ref[j];
	    vn += v[j]*nor[j];
	}
	for (j = 0; j < dim; j++)
	{
	    v[j] += vel_ref[j] - 2.0*vn*nor[j];
	    st_tmp.momn[j] = v[j]*st_tmp.dens;
	}
	st_tmp.pres = EosPressure(&st_tmp);
	if (st_tmp.pres < min_pres) st_tmp.pres = min_pres;
	st_tmp.engy = EosEnergy(&st_tmp);

	index = d_index(ic,top_gmax,dim);
	m_vst->dens[index] = st_tmp.dens;
	if(eqn_params->multi_comp_non_reactive == YES)
	{
	    int ii;
	    for(ii = 0; ii < eqn_params->n_comps; ii++)
	    {
		m_vst->pdens[ii][index] = st_tmp.pdens[ii];
	    }
	}
	m_vst->engy[index] = st_tmp.engy;
	m_vst->pres[index] = st_tmp.pres;
	for (j = 0; j < 3; j++)
	    m_vst->momn[j][index] = 0.0;
	if (dim == 1)
	    m_vst->momn[0][index] = st_tmp.momn[0];
	else if (dim == 2)
	    for (j = 0; j < 2; j++)
		m_vst->momn[j][index] = st_tmp.momn[j];
	else if (dim == 3)
	    for (j = 0; j < 3; j++)
		m_vst->momn[j][index] = st_tmp.momn[j];
	m_vst->comp[index] = comp;
}	/* end setNeumannStatesforParab */

void G_CARTESIAN::sampleVelocity()
{
	switch (dim)
	{
	case 2:
	    return sampleVelocity2d();
	case 3:
	    return sampleVelocity3d();
	}
}	/* end sampleVelocity */

void G_CARTESIAN::sampleVelocity3d()
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
	double dens;

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
                    for (k = imin[2]; k <= imax[2]; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (k = imin[2]; k <= imax[2]; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (k = imin[2]; k <= imax[2]; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
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
                    for (j = imin[1]; j <= imax[1]; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (j = imin[1]; j <= imax[1]; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (j = imin[1]; j <= imax[1]; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
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
                    for (i = imin[0]; i <= imax[0]; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (i = imin[0]; i <= imax[0]; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (i = imin[0]; i <= imax[0]; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
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
}	/* end sampleVelocity3d */

void G_CARTESIAN::sampleVelocity2d()
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
	double dens;

	if (front->step < sample->start_step || front->step > sample->end_step)
            return;
        if ((front->step - sample->start_step)%sample->step_interval)
            return;
        if (step != front->step)
            step = front->step;
	
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
            for (j = imin[1]; j <= imax[1]; ++j)
            {
                index = d_index2d(i,j,top_gmax);
		dens = field.dens[index];
                velo1 = field.momn[0][index]/dens;
                index = d_index2d(i+1,j,top_gmax);
		dens = field.dens[index];
                velo2 = field.momn[0][index]/dens;
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
            }
            fclose(sfile);
            sprintf(sname,"%s/vertical-y-%d-%d.xg",out_name,step,count++);
            sfile = fopen(sname,"w");
            for (j = imin[1]; j <= imax[1]; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                dens = field.dens[index];
                velo1 = field.momn[1][index]/dens;
                index = d_index2d(i+1,j,top_gmax);
                dens = field.dens[index];
                velo2 = field.momn[1][index]/dens;
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
            for (i = imin[0]; i <= imax[0]; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                dens = field.dens[index];
                velo1 = field.momn[0][index]/dens;
                index = d_index2d(i,j+1,top_gmax);
                dens = field.dens[index];
                velo2 = field.momn[0][index]/dens;
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);
            sprintf(sname,"%s/horizontal-y-%d-%d.xg",out_name,step,count++);
            sfile = fopen(sname,"w");
            for (i = imin[0]; i <= imax[0]; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                dens = field.dens[index];
                velo1 = field.momn[1][index]/dens;
                index = d_index2d(i,j+1,top_gmax);
                dens = field.dens[index];
                velo2 = field.momn[1][index]/dens;
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);
            break;
        }
}	/* end sampleVelocity2d */

void G_CARTESIAN::numericalFlux(
	POINTER scheme_params,
	SWEEP *sweep,
	FSWEEP *fsweep,
	int n)
{
	switch (eqn_params->num_scheme)
	{
	case TVD_FIRST_ORDER:
	case TVD_SECOND_ORDER:
	case TVD_FOURTH_ORDER:
	    TVD_flux(scheme_params,sweep,fsweep,n);
	    break;
	case WENO_FIRST_ORDER:
	case WENO_SECOND_ORDER:
	case WENO_FOURTH_ORDER:
	    WENO_flux(scheme_params,sweep,fsweep,n);
	    break;
	default:
	    (void) printf("Unknow numerical scheme\n");
	    clean_up(ERROR);
	}
}	/* numericalFlux */


void G_CARTESIAN::scatMeshVst(SWEEP *m_vst)
{
	int i,j,k,l,index;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_vst->dens[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_vst->dens[index] = array[index];
	    }
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_vst->engy[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_vst->engy[index] = array[index];
	    }
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_vst->pres[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_vst->pres[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index1d(i,top_gmax);
		    array[index] = m_vst->momn[l][index];
	    	}
	    	scatMeshArray();
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index1d(i,top_gmax);
		    m_vst->momn[l][index] = array[index];
	    	}
	    }

            if(eqn_params->multi_comp_non_reactive == YES)
            {
            int ii;
            for(ii = 0; ii < eqn_params->n_comps; ii++)
            {
                for (i = imin[0]; i <= imax[0]; ++i)
                {
                    index = d_index1d(i,top_gmax);
                    array[index] = m_vst->pdens[ii][index];
                }
                scatMeshArray();
                for (i = 0; i <= top_gmax[0]; i++)
                {
                    index = d_index1d(i,top_gmax);
                    m_vst->pdens[ii][index] = array[index];
                }
            } //for ii 
            } //if
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_vst->dens[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_vst->dens[index] = array[index];
	    }
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_vst->engy[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_vst->engy[index] = array[index];
	    }
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_vst->pres[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_vst->pres[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index2d(i,j,top_gmax);
		    array[index] = m_vst->momn[l][index];
	    	}
	    	scatMeshArray();
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index2d(i,j,top_gmax);
		    m_vst->momn[l][index] = array[index];
	    	}
	    }

            if(eqn_params->multi_comp_non_reactive == YES)
            {
            int ii;
            for(ii = 0; ii < eqn_params->n_comps; ii++)
            {
                for (j = imin[1]; j <= imax[1]; ++j)
                for (i = imin[0]; i <= imax[0]; ++i)
                {
                    index = d_index2d(i,j,top_gmax);
                    array[index] = m_vst->pdens[ii][index];
                }
                scatMeshArray();
                for (j = 0; j <= top_gmax[1]; j++)
                for (i = 0; i <= top_gmax[0]; i++)
                {
                    index = d_index2d(i,j,top_gmax);
                    m_vst->pdens[ii][index] = array[index];
                }
            } //for ii 
            } //if
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
                array[index] = m_vst->dens[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
		index = d_index3d(i,j,k,top_gmax);
                m_vst->dens[index] = array[index];
	    }
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
                array[index] = m_vst->engy[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                m_vst->engy[index] = array[index];
	    }
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
                array[index] = m_vst->pres[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                m_vst->pres[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (k = imin[2]; k <= imax[2]; ++k)
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index3d(i,j,k,top_gmax);
                    array[index] = m_vst->momn[l][index];
	    	}
	    	scatMeshArray();
		state_reflect(l,array);	//Dan	FIXME
            	for (k = 0; k <= top_gmax[2]; k++)
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
		    index = d_index3d(i,j,k,top_gmax);
                    m_vst->momn[l][index] = array[index];
	    	}
	    }

            if(eqn_params->multi_comp_non_reactive == YES)
            {
            int ii;
            for(ii = 0; ii < eqn_params->n_comps; ii++)
            {
                for (k = imin[2]; k <= imax[2]; ++k)
                for (j = imin[1]; j <= imax[1]; ++j)
                for (i = imin[0]; i <= imax[0]; ++i)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    array[index] = m_vst->pdens[ii][index];
                }
                scatMeshArray();
                for (k = 0; k <= top_gmax[2]; k++)
                for (j = 0; j <= top_gmax[1]; j++)
                for (i = 0; i <= top_gmax[0]; i++)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    m_vst->pdens[ii][index] = array[index];
                }
            } //for ii 
            } //if
	}
}	/* end scatMeshStates */

void G_CARTESIAN::copyMeshVst(
	SWEEP m_vst_orig,
	SWEEP *m_vst)
{
	int i,j,k,l,index;
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		m_vst->dens[index] = m_vst_orig.dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] = m_vst_orig.pdens[ii][index];
                    }
                }
		m_vst->engy[index] = m_vst_orig.engy[index];
		m_vst->pres[index] = m_vst_orig.pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = m_vst_orig.momn[l][index];
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		m_vst->dens[index] = m_vst_orig.dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] = m_vst_orig.pdens[ii][index];
                    }
                }
		m_vst->engy[index] = m_vst_orig.engy[index];
		m_vst->pres[index] = m_vst_orig.pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = m_vst_orig.momn[l][index];
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		m_vst->dens[index] = m_vst_orig.dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] = m_vst_orig.pdens[ii][index];
                    }
                }
		m_vst->engy[index] = m_vst_orig.engy[index];
		m_vst->pres[index] = m_vst_orig.pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = m_vst_orig.momn[l][index];
	    }
	}
}	/* end copyMeshVst */

void G_CARTESIAN::copyToMeshVst(
	SWEEP *m_vst)
{
	int i,j,k,l,index;
	double *dens = field.dens;
        double **pdens = field.pdens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
	double *gamma = field.gamma;

	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		m_vst->dens[index] = dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] = pdens[ii][index];
                    }
                }
		m_vst->engy[index] = engy[index];
		m_vst->pres[index] = pres[index];
//		m_vst->gamma[index] = gamma[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = momn[l][index];
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		m_vst->dens[index] = dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] = pdens[ii][index];
                    }
                }
		m_vst->engy[index] = engy[index];
		m_vst->pres[index] = pres[index];
//		m_vst->gamma[index] = gamma[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = momn[l][index];
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		m_vst->dens[index] = dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] = pdens[ii][index];
                    }
                }
		m_vst->engy[index] = engy[index];
		m_vst->pres[index] = pres[index];
//		m_vst->gamma[index] = gamma[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = momn[l][index];
	    }
	}
}	/* end copyToMeshVst */

void G_CARTESIAN::copyFromMeshVst(
	SWEEP m_vst)
{
	int i,j,k,l,index;
	STATE state;
	COMPONENT comp;
	double *dens = field.dens;
        double **pdens = field.pdens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
	
	//GFM
	if(eqn_params->tracked)
	{
	    get_ghost_state(m_vst, 2, 0);
	    get_ghost_state(m_vst, 3, 1);
	    scatMeshGhost();
	}

	state.dim = dim;
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		state.dens = m_vst.dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        state.pdens[ii] = m_vst.pdens[ii][index];
                    }
                }
		state.engy = m_vst.engy[index];
		state.pres = m_vst.pres[index];
		for (l = 0; l < dim; ++l)
		    state.momn[l] = m_vst.momn[l][index];
		if (gas_comp(top_comp[index]))
		{
		    state.eos = &(eqn_params->eos[comp]);
		    checkCorrectForTolerance(&state);
		}
		dens[index] = state.dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        pdens[ii][index] = state.pdens[ii];
                    }
                }
		engy[index] = state.engy;
		pres[index] = state.pres;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		state.dens = m_vst.dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        state.pdens[ii] = m_vst.pdens[ii][index];
                    }
                }
		state.engy = m_vst.engy[index];
		state.pres = m_vst.pres[index];
		for (l = 0; l < dim; ++l)
		    state.momn[l] = m_vst.momn[l][index];
		if (gas_comp(top_comp[index]))
		{
		    state.eos = &(eqn_params->eos[comp]);
		    checkCorrectForTolerance(&state);
		}
		dens[index] = state.dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        pdens[ii][index] = state.pdens[ii];
                    }
                }
		engy[index] = state.engy;
		pres[index] = state.pres;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		state.dens = m_vst.dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        state.pdens[ii] = m_vst.pdens[ii][index];
                    }
                }
		state.engy = m_vst.engy[index];
		state.pres = m_vst.pres[index];
		for (l = 0; l < dim; ++l)
		    state.momn[l] = m_vst.momn[l][index];
		if (gas_comp(top_comp[index]))
		{
		    state.eos = &(eqn_params->eos[comp]);
		    checkCorrectForTolerance(&state);
		}
		dens[index] = state.dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        pdens[ii][index] = state.pdens[ii];
                    }
                }
		engy[index] = state.engy;
		pres[index] = state.pres;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	}
}	/* end copyFromMeshVst */

void G_CARTESIAN::appendStencilBuffer2d(
	SWEEP *vst,
	SWEEP *m_vst,
	int i,
	int dir)
{
	int		i1,i2,k,offset,index0,index;
	INTERFACE 	*intfc = front->interf;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	HYPER_SURF 	*hs;
	double 		crx_coords[MAXD];
	STATE 		*state;
	int		comp, icoords[3];

	switch (dir)
	{
	case 0:
	    i2 = i;
	    if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		i1 = imin[0];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imin[0] + k - 1;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = -m_vst->momn[1][index];
		    vst->momn[1][3-k] = m_vst->momn[0][index];
		    vst->momn[2][3-k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imin[0] - k;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = m_vst->momn[0][index];
		    vst->momn[1][3-k] = m_vst->momn[1][index];
		    vst->momn[2][3-k] = m_vst->momn[2][index];
		}
	    }
	    else
	    {
		i1 = imin[0];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,icoords,ldir[dir],
				comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[dir]), 0);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,0,0,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d: "
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }

	    if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		i1 = imax[0];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imax[0] - k + 1;
		    offset = imax[0] - imin[0] + 3;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = -m_vst->momn[1][index];
		    vst->momn[1][offset+k] = m_vst->momn[0][index];
		    vst->momn[2][offset+k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imax[0] + k;
		    index = d_index2d(i1,i2,top_gmax);
		    offset = imax[0] - imin[0] + 3;
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = m_vst->momn[0][index];
		    vst->momn[1][offset+k] = m_vst->momn[1][index];
		    vst->momn[2][offset+k] = m_vst->momn[2][index];
		}
	    }
	    else
	    {
		i1 = imax[0];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,icoords,rdir[dir],
				comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[0]), 1);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    offset = imax[dir] - imin[dir] + nrad;
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,1,
					offset,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d"
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }
	    break;
	case 1:
	    i1 = i;
	    if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		i2 = imin[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imin[1] + k - 1;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = -m_vst->momn[1][index];
		    vst->momn[1][3-k] = m_vst->momn[0][index];
		    vst->momn[2][3-k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		i2 = imin[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imin[1] - k;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = m_vst->momn[1][index];
		    vst->momn[1][3-k] = m_vst->momn[0][index];
		    vst->momn[2][3-k] = 0.0;
		}
	    }
	    else
	    {
		i2 = imin[1];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,icoords,ldir[dir],
				comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[dir]), 0);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,0,0,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d"
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }

	    if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		i2 = imax[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imax[1] - k + 1;
		    offset = imax[1] - imin[1] + 3;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = -m_vst->momn[1][index];
		    vst->momn[1][offset+k] = m_vst->momn[0][index];
		    vst->momn[2][offset+k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		i2 = imax[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imax[1] + k;
		    offset = imax[1] - imin[1] + 3;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = m_vst->momn[1][index];
		    vst->momn[1][offset+k] = m_vst->momn[0][index];
		    vst->momn[2][offset+k] = 0.0;
		}
	    }
	    else
	    {
		i2 = imax[1];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,icoords,rdir[dir],
				comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[dir]), 0);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    offset = imax[dir] - imin[dir] + nrad;
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,1,
					offset,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d"
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }
	}

}	/* end appendStencilBuffer2d */

void G_CARTESIAN::appendStencilBuffer3d(
	SWEEP *vst,
	SWEEP *m_vst,
	int i1,
	int i2,
	int dir)
{
	int i,j,k,l,offset,index;
	INTERFACE *intfc = front->interf;

	switch (dir)
	{
	case 0:
	    j = i1;	k = i2;
	    if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    i = imin[0] - l;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = m_vst->momn[0][index];
		    vst->momn[1][3-l] = m_vst->momn[1][index];
		    vst->momn[2][3-l] = m_vst->momn[2][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    i = imax[0] + l;
		    index = d_index3d(i,j,k,top_gmax);
		    offset = imax[0] - imin[0] + 3;
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = m_vst->momn[0][index];
		    vst->momn[1][offset+l] = m_vst->momn[1][index];
		    vst->momn[2][offset+l] = m_vst->momn[2][index];
		}
	    }
	    break;
	case 1:
	    k = i1;	i = i2;
	    if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imin[1] - l;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = m_vst->momn[1][index];
		    vst->momn[1][3-l] = m_vst->momn[2][index];
		    vst->momn[2][3-l] = m_vst->momn[0][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imin[1] + l - 1;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = -m_vst->momn[1][index];
		    vst->momn[1][3-l] = m_vst->momn[2][index];
		    vst->momn[2][3-l] = m_vst->momn[0][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imax[1] + l;
		    index = d_index3d(i,j,k,top_gmax);
		    offset = imax[1] - imin[1] + 3;
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = m_vst->momn[1][index];
		    vst->momn[1][offset+l] = m_vst->momn[2][index];
		    vst->momn[2][offset+l] = m_vst->momn[0][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imax[1] - l + 1;
		    offset = imax[1] - imin[1] + 3;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = -m_vst->momn[1][index];
		    vst->momn[1][offset+l] = m_vst->momn[2][index];
		    vst->momn[2][offset+l] = m_vst->momn[0][index];
		}
	    }
	    break;
	case 2:
	    i = i1;	j = i2;
	    if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imin[2] + l - 1;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = -m_vst->momn[2][index];
		    vst->momn[1][3-l] = m_vst->momn[0][index];
		    vst->momn[2][3-l] = m_vst->momn[1][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imax[2] - l + 1;
		    offset = imax[2] - imin[2] + 3;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = -m_vst->momn[2][index];
		    vst->momn[1][offset+l] = m_vst->momn[0][index];
		    vst->momn[2][offset+l] = m_vst->momn[1][index];
		}
	    }
	    
	    if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imin[2] - l; 
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = m_vst->momn[2][index];
		    vst->momn[1][3-l] = m_vst->momn[0][index];
		    vst->momn[2][3-l] = m_vst->momn[1][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imax[2] + l;
		    offset = imax[2] - imin[2] + 3;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = m_vst->momn[2][index];
		    vst->momn[1][offset+l] = m_vst->momn[0][index];
		    vst->momn[2][offset+l] = m_vst->momn[1][index];
		}
	    }
	}
}	/* end appendStencilBuffer3d */

void G_CARTESIAN::scatMeshStates()
{
	SWEEP vst;
	allocMeshVst(&vst);
	copyToMeshVst(&vst);
	scatMeshVst(&vst);
	copyFromMeshVst(vst);
	freeVst(&vst);
}	/* end scatMeshStates */

void G_CARTESIAN::freeVst(
	SWEEP *vst)
{
        FT_FreeThese(5,vst->dens,vst->pdens,vst->engy,vst->pres,vst->momn);
}	/* end freeVstFlux */

void G_CARTESIAN::freeFlux(
	FSWEEP *flux)
{
	FT_FreeThese(4,flux->dens_flux,flux->pdens_flux,flux->engy_flux,flux->momn_flux);
}

void G_CARTESIAN::addMeshFluxToVst(
	SWEEP *m_vst,
	FSWEEP m_flux,
	double chi)
{
	int 		i,j,k,l,index;
	double		ke,c,u;
	EOS_PARAMS	*eos;
	STATE		st;
	int		comp;
	double		temp;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		if (!gas_comp(comp))
		{
		    m_vst->dens[index] = 0.0;
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            m_vst->pdens[ii][index] = 0.0;
                        }
                    }
		    m_vst->engy[index] = 0.0;
		    for (l = 0; l < dim; ++l)
		    	m_vst->momn[l][index] = 0.0; 
		    continue;
		}
		eos = &(eqn_params->eos[comp]);

		m_vst->dens[index] += chi*m_flux.dens_flux[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] += chi*m_flux.pdens_flux[ii][index];
                    }
                }
		m_vst->engy[index] += chi*m_flux.engy_flux[index];
		ke = u = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    m_vst->momn[l][index] += 
			chi*m_flux.momn_flux[l][index];
		    ke += sqr(m_vst->momn[l][index]);
		    u += sqr(m_vst->momn[l][index]);
		}
		
		CovertVstToState(&st, m_vst, eos, index, dim);
		checkCorrectForTolerance(&st);
		m_vst->dens[index] = st.dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] = st.pdens[ii];
                    }
                }
		m_vst->pres[index] = st.pres;
		m_vst->engy[index] = st.engy;
		u = sqrt(u)/m_vst->dens[index];
		c = EosSoundSpeed(&st);
		temp = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
                if (max_speed < temp)
                    max_speed = temp;
	    }
            // scaling partial density
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                for (i = imin[0]; i <= imax[0]; ++i)
                {
                    index = d_index1d(i,top_gmax);
                    int ii;
                    double sum = 0.0;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                        sum += m_vst->pdens[ii][index];
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                        m_vst->pdens[ii][index] *= m_vst->dens[index]/sum;
                }
            }
	    scatMeshVst(m_vst);
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		if (!gas_comp(comp))
		{
		    m_vst->dens[index] = 0.0;
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            m_vst->pdens[ii][index] = 0.0;
                        }
                    }
		    m_vst->engy[index] = 0.0;
		    for (l = 0; l < dim; ++l)
		    	m_vst->momn[l][index] = 0.0; 
		    continue;
		}
		eos = &(eqn_params->eos[comp]);

		m_vst->dens[index] += chi*m_flux.dens_flux[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] += chi*m_flux.pdens_flux[ii][index];
                    }
                }
		m_vst->engy[index] += chi*m_flux.engy_flux[index];
		ke = u = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    m_vst->momn[l][index] += 
			chi*m_flux.momn_flux[l][index];
		    ke += sqr(m_vst->momn[l][index]);
		    u += sqr(m_vst->momn[l][index]);
		}
		
		CovertVstToState(&st, m_vst, eos, index, dim);
		checkCorrectForTolerance(&st);
		m_vst->dens[index] = st.dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] = st.pdens[ii];
                    }
                }
		m_vst->pres[index] = st.pres;
		m_vst->engy[index] = st.engy;
		u = sqrt(u)/m_vst->dens[index];
		c = EosSoundSpeed(&st);
		temp = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
                if (max_speed < temp)
                    max_speed = temp;
	    }
            // scaling partial density
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                for (j = imin[1]; j <= imax[1]; ++j)
                for (i = imin[0]; i <= imax[0]; ++i)
                {
                    index = d_index2d(i,j,top_gmax);
                    int ii;
                    double sum = 0.0;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                        sum += m_vst->pdens[ii][index];
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                        m_vst->pdens[ii][index] *= m_vst->dens[index]/sum;
                }
            }
	    scatMeshVst(m_vst);
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
                if (!gas_comp(comp))
                {
                    m_vst->dens[index] = 0.0;
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            m_vst->pdens[ii][index] = 0.0;
                        }
                    }
                    m_vst->engy[index] = 0.0;
                    for (l = 0; l < dim; ++l)
                        m_vst->momn[l][index] = 0.0;
                    continue;
                }
                eos = &(eqn_params->eos[comp]);

		m_vst->dens[index] += chi*m_flux.dens_flux[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] += chi*m_flux.pdens_flux[ii][index];
                    }
                }
		m_vst->engy[index] += chi*m_flux.engy_flux[index];
		ke = u = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    m_vst->momn[l][index] += 
				chi*m_flux.momn_flux[l][index];
		    ke += sqr(m_vst->momn[l][index]);
		    u += sqr(m_vst->momn[l][index]);
		}
		
		CovertVstToState(&st, m_vst, eos, index, dim);
		checkCorrectForTolerance(&st);
		m_vst->dens[index] = st.dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        m_vst->pdens[ii][index] = st.pdens[ii];
                    }
                }
		m_vst->pres[index] = st.pres;
		m_vst->engy[index] = st.engy;
		u = sqrt(u)/m_vst->dens[index];
		c = EosSoundSpeed(&st);
		temp = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
                if (max_speed < temp)
                    max_speed = temp;
	    }
            // scaling partial density
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                for (k = imin[2]; k <= imax[2]; ++k)
                for (j = imin[1]; j <= imax[1]; ++j)
                for (i = imin[0]; i <= imax[0]; ++i)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    int ii;
                    double sum = 0.0;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                        sum += m_vst->pdens[ii][index];
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                        m_vst->pdens[ii][index] *= m_vst->dens[index]/sum;
                }
            }
	    scatMeshVst(m_vst);
	}
}	/* end addMeshFluxToVst */


void G_CARTESIAN::appendGhostBuffer(
	SWEEP *vst,
	SWEEP *m_vst,
	int n,
	int *icoords,
	int idir,
	int nb)
{
	int		i,j,k,index,ic[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	HYPER_SURF 	*hs;
	COMPONENT 	comp;
	double 		crx_coords[MAXD];
	STATE 		*state,ghost_st;
	int		ind2[2][2] = {{0,1},{1,0}};
	int		ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int 		ic_next[MAXD];

	if (debugging("append_buffer"))
		printf("Entering appendGhostBuffer()\n");
	for (i = 0; i < dim; ++i) ic[i] = icoords[i];
	
	index = d_index(ic,top_gmax,dim);
	comp = cell_center[index].comp;
	
	switch(nb)
	{
	case 0:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] - i;
		index = d_index(ic,top_gmax,dim);
		    
		if (!needBufferFromIntfc(comp,cell_center[index].comp))
		{
		    vst->dens[nrad-i] = m_vst->dens[index];
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            vst->pdens[ii][nrad-i] = m_vst->pdens[ii][index];
                        }
                    }
		    vst->engy[nrad-i] = m_vst->engy[index];
		    vst->pres[nrad-i] = m_vst->pres[index];

		    for (j = 0; j < 3; j++)
			vst->momn[j][nrad-i] = 0.0;
		    if (dim == 1)
			vst->momn[0][nrad-i] = m_vst->momn[0][index];
		    else if (dim == 2)
			for(j = 0; j < 2; j++)
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3)
			for (j = 0; j < 3; j++)
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind3[idir][j]][index];
		}
		else
		{
		    for (k = 0; k < dim; ++k)
			ic_next[k] = ic[k];
		    ic_next[idir]++;
		    if (!FT_StateStructAtGridCrossing(front,ic_next,
			ldir[idir],comp,(POINTER*)&state,&hs,crx_coords))
		    {
		    	(void) printf("In appendGhostBuffer()\n");
		    	(void) printf("ERROR: No crossing found!\n");
		    	(void) print_int_vector("icoords=", ic_next,3,"\n");
		    	(void) printf("direction: %s side %d\n",
		           		grid_direction_name(rdir[idir]), nb);
		    	clean_up(ERROR);
		    }
                    state->eos = &(eqn_params->eos[comp]);
		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,0,i,comp);
		    	break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,
					idir,nb,0,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	GFMGhostState(ic_next,ic,comp,&ghost_st,-1,m_vst,idir);
		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[nrad-k] = ghost_st.dens;
                            if(eqn_params->multi_comp_non_reactive == YES)
                            {
                                int ii;
                                for(ii = 0; ii < eqn_params->n_comps; ii++)
                                {
                                    vst->pdens[ii][nrad-k] = ghost_st.pdens[ii];
                                }
                            }
		    	    vst->engy[nrad-k] = ghost_st.engy;
		    	    vst->pres[nrad-k] = ghost_st.pres;
			
			    for (j=0; j < 3; j++)
			    	    vst->momn[j][nrad-k] = 0.0;
			    if (dim == 1)
				vst->momn[0][nrad-k] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for (j=0; j < 2; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for (j = 0; j < 3; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
					wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=", icoords,3,"\n");
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	    break;
	case 1:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] + i;
		index = d_index(ic,top_gmax,dim);
		if (!needBufferFromIntfc(comp,cell_center[index].comp))
		{
		    vst->dens[n+nrad+i-1] = m_vst->dens[index];
                    if(eqn_params->multi_comp_non_reactive == YES)
                    {
                        int ii;
                        for(ii = 0; ii < eqn_params->n_comps; ii++)
                        {
                            vst->pdens[ii][n+nrad+i-1] = m_vst->pdens[ii][index];
                        }
                    }
		    vst->engy[n+nrad+i-1] = m_vst->engy[index];
		    vst->pres[n+nrad+i-1] = m_vst->pres[index];
		    
		    for (j = 0; j < 3; j++)
			vst->momn[j][n+nrad+i-1] = 0.0;
		    if (dim == 1)
			vst->momn[0][n+nrad+i-1] = 
			         	m_vst->momn[0][index];
		    else if (dim == 2)
			for(j = 0; j < 2; j++)
			    	vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3)
			for (j = 0; j < 3; j++)
			    vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind3[idir][j]][index];
		}
		else
		{
		    for (k = 0; k < dim; ++k)
			ic_next[k] = ic[k];
		    ic_next[idir]--;
		    if (!FT_StateStructAtGridCrossing(front,ic_next,
			rdir[idir],comp,(POINTER*)&state,&hs,crx_coords))
		    {
		    	(void) printf("In appendGhostBuffer()\n");
		    	(void) printf("ERROR: No crossing found!\n");
		    	(void) print_int_vector("icoords=",ic_next,3,"\n");
		    	(void) printf("direction: %s side %d\n",
		            	grid_direction_name(rdir[idir]), nb);
		    	clean_up(ERROR);
		    }
                    state->eos = &(eqn_params->eos[comp]);
		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,n,i,comp);
		    	break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,idir,nb,
						n,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	GFMGhostState(ic_next,ic,comp,&ghost_st,1,m_vst,idir);
		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[n+nrad+k-1] = ghost_st.dens;
                            if(eqn_params->multi_comp_non_reactive == YES)
                            {
                                int ii;
                                for(ii = 0; ii < eqn_params->n_comps; ii++)
                                {
                                    vst->pdens[ii][n+nrad+k-1] = ghost_st.pdens[ii];
                                }
                            }
		    	    vst->engy[n+nrad+k-1] = ghost_st.engy;
		    	    vst->pres[n+nrad+k-1] = ghost_st.pres;
			
			    for(j=0; j<3; j++)
			    	vst->momn[j][n+nrad+k-1] = 0.0;
			    if (dim == 1)
				vst->momn[0][n+nrad+k-1] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for(j = 0; j < 2; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for(j = 0; j < 3; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
				wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=",icoords,3,"\n");
		    	(void) printf("nb = %d\n",nb);
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	}
}	/* end appendGhostBuffer */

//ghost fluid method.

void G_CARTESIAN::solve_exp_value()
{
	int		i, j, k, n;
	int		index;
	double		**gnor = eqn_params->gnor;

	fflush(NULL);

	get_normal_from_front();

	if (dim == 1)
	{
	    for(k=0; k<dim; k++)
	    {
		for (i = imin[0]; i <= imax[0]; ++i)
		{
	 	    index = d_index1d(i,top_gmax);
		    array[index] = gnor[k][index];
		}
		scatMeshArray();
        	for (i = 0; i <= top_gmax[0]; i++)
        	{
	    	    index  = d_index1d(i,top_gmax);
	    	    gnor[k][index] = array[index];
		}
	    }
	}
	else if (dim == 2)
	{
	    for(k=0; k<dim; k++)
	    {
		for (j = imin[1]; j <= imax[1]; ++j)
		for (i = imin[0]; i <= imax[0]; ++i)
		{
	 	    index = d_index2d(i,j,top_gmax);
		    array[index] = gnor[k][index];
		}
		scatMeshArray();
        	for (j = 0; j <= top_gmax[1]; j++)
        	for (i = 0; i <= top_gmax[0]; i++)
        	{
	    	    index  = d_index2d(i,j,top_gmax);
	    	    gnor[k][index] = array[index];
		}
	    }
	}
	else
	{
	    for(k=0; k<dim; k++)
	    {
		for (n = imin[2]; n <= imax[2]; ++n)
		for (j = imin[1]; j <= imax[1]; ++j)
		for (i = imin[0]; i <= imax[0]; ++i)
		{
	    	    index = d_index3d(i,j,n,top_gmax);
	    	    array[index] = gnor[k][index];
		}
		scatMeshArray();
        	for (n = 0; n <= top_gmax[2]; n++)
        	for (j = 0; j <= top_gmax[1]; j++)
        	for (i = 0; i <= top_gmax[0]; i++)
        	{
	    	    index  = d_index3d(i,j,n,top_gmax);
	    	    gnor[k][index] = array[index];
		}
	    }
	}
}

void G_CARTESIAN::scatMeshGhost()
{
	int		i, j, k, n, index;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
        double          ***Gpdens = eqn_params->Gpdens;
	double		**Gpres = eqn_params->Gpres;

	if(dim == 2)
	{
	for(k=0; k<2; k++)
	{
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gdens[k][index];
	}
	scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gdens[k][index] = array[index];
	}
	
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gpres[k][index];
	}
	scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gpres[k][index] = array[index];
	}

	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gvel[k][0][index];
	}
	scatMeshArray();
	state_reflect(0,array);	//Dan	FIXME
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gvel[k][0][index] = array[index];
	}

	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gvel[k][1][index];
	}
	scatMeshArray();
	state_reflect(1,array);	//Dan	FIXME
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gvel[k][1][index] = array[index];
	}

        if(eqn_params->multi_comp_non_reactive == YES)
        {
	    int ii;
	    for(ii = 0; ii < eqn_params->n_comps; ii++)
	    {
		for (j = imin[1]; j <= imax[1]; ++j)
		for (i = imin[0]; i <= imax[0]; ++i)
		{
		    index = d_index2d(i,j,top_gmax);
		    array[index] = Gpdens[k][ii][index];
		}
		scatMeshArray();
		for (j = 0; j <= top_gmax[1]; j++)
		for (i = 0; i <= top_gmax[0]; i++)
		{
		    index  = d_index2d(i,j,top_gmax);
		    Gpdens[k][ii][index] = array[index];
		}
	    } //for ii 
        } //if

	}    //for k
	}    //if dim == 2
	else if(dim == 3)
	{
	for(k=0; k<2; k++)
	{
	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gdens[k][index];
	}
	scatMeshArray();
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gdens[k][index] = array[index];
	}

	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gpres[k][index];
	}
	scatMeshArray();
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gpres[k][index] = array[index];
	}
	
	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gvel[k][0][index];
	}
	scatMeshArray();
	state_reflect(0,array);	//Dan	FIXME
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gvel[k][0][index] = array[index];
	}

	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gvel[k][1][index];
	}
	scatMeshArray();
	state_reflect(1,array);	//Dan	FIXME
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gvel[k][1][index] = array[index];
	}
	
	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gvel[k][2][index];
	}
	scatMeshArray();
	state_reflect(2,array);	//Dan	FIXME
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gvel[k][2][index] = array[index];
	}

        if(eqn_params->multi_comp_non_reactive == YES)
        {
	    int ii;
	    for(ii = 0; ii < eqn_params->n_comps; ii++)
	    {
		for (n = imin[2]; n <= imax[2]; ++n)
		for (j = imin[1]; j <= imax[1]; ++j)
		for (i = imin[0]; i <= imax[0]; ++i)
		{
		    index = d_index3d(i,j,n,top_gmax);
		    array[index] = Gpdens[k][ii][index];
		}
		scatMeshArray();
		for (n = 0; n <= top_gmax[2]; n++)
		for (j = 0; j <= top_gmax[1]; j++)
		for (i = 0; i <= top_gmax[0]; i++)
		{
		    index  = d_index3d(i,j,n,top_gmax);
		    Gpdens[k][ii][index] = array[index];
		}
	    } //for ii 
        } //if

	}    //for k
	}    //for dim == 3
}

#define	corner_index(p,i,gr)	irint(floor(((p)-(gr)->L[i])/(gr)->h[i]-0.5))

boolean	find_block(
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

boolean G_CARTESIAN::get_ave_normal(
	int		*ic,
	int		***norset)
{
	double		f;
	int		i, j, k, n, ic1[3], ic2[3], dir, num;
	boolean		found;
	int		index0, index;
	double		**gnor = eqn_params->gnor;
	
	found = NO;

	for(i=0; i<dim; i++)
	for(j=0; j<2; j++)
	{
		dir = j == 0 ? -1 : 1;
		ft_assign(ic1, ic, 3*INT);
		ic1[i] = ic[i] + dir;

		if(ic1[i] < 0 || ic1[i] > top_gmax[i])
		    continue;
		if(norset[ic1[0]][ic1[1]][ic1[2]] == 1)
		    found = YES;
	}
	if(!found)
	    return NO;
	
	index0  = d_index(ic,top_gmax,dim);

	gnor[0][index0] = 0.0;
	if (dim > 1)
	    gnor[1][index0] = 0.0;
	if (dim > 2)
	    gnor[2][index0] = 0.0;

	num = 0;
	for(i=ic[0]-1; i<=ic[0]+1; i++)
	for(j=ic[1]-1; j<=ic[1]+1; j++)
	for(k=ic[2]-1; k<=ic[2]+1; k++)
	{
	    if(i < 0 || i > top_gmax[0] || 
	       j < 0 || j > top_gmax[1] || 
	       k < 0 || k > top_gmax[2]) 
		continue;
	    if(norset[i][j][k] != 1)
		continue;

	    ic2[0] = i;
	    ic2[1] = j;
	    ic2[2] = k;
	    index  = d_index(ic2,top_gmax,dim);
		    
	    //do not use length weighted normal direction
	    gnor[0][index0] += gnor[0][index];
	    if(dim > 1)
	    	gnor[1][index0] += gnor[1][index];
	    if(dim > 2)
		gnor[2][index0] += gnor[2][index];
	    num++;
	}
	
	f = 0.0;
	for(n=0; n<dim; n++)
	    f += sqr(gnor[n][index0]);
	f = sqrt(f);

	if(f < 1.0e-6)
	{
	    gnor[0][index0] = 0.0;
	    if (dim > 1) gnor[1][index0] = 0.0;
	    if (dim > 2) gnor[2][index0] = 0.0;
	}
	else
	{
	    gnor[0][index0] /= f;
	    if (dim > 1) gnor[1][index0] /= f;
	    if (dim > 2) gnor[2][index0] /= f;
	}

	return YES;
}

boolean	find_block(double*,int*,double*,RECT_GRID*);

void	get_normal_from_front();

//it will fill gnor field by interface normals
void G_CARTESIAN::get_normal_from_front()
{
	INTERFACE               *intfc = front->interf;
	RECT_GRID		*rgr = front->rect_grid;
	HYPER_SURF              *hs;
	HYPER_SURF_ELEMENT      *hse;
	POINT                   *p;
	int			i,j,k,n, num;
	int			ic[3];
	double			curv,nor[3],d[3],f,d1,d2,d3,*pt,tol;
	int			ix, iy, iz, index;
	boolean			found;
	double			**gnor = eqn_params->gnor;
	static	int		***norset;
	int ict[3];
	double ptt[3];
	boolean status;

	if (norset == NULL)
	    FT_TriArrayMemoryAlloc((POINTER*)&norset,top_gmax[0]+1,
				top_gmax[1]+1,top_gmax[2]+1,INT);

	tol = hmin*1.0e-6;

	for (i = 0; i <= top_gmax[0]; i++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (k = 0; k <= top_gmax[2]; k++)
	{
	    ic[0] = i;
	    ic[1] = j;
	    ic[2] = k;
	    index = d_index(ic,top_gmax,dim);
	    gnor[0][index] = 0.0;
	    if (dim > 1)
	    	gnor[1][index] = 0.0;
	    if (dim > 2)
		gnor[2][index] = 0.0;
	    norset[i][j][k] = 0;
	}

	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (Boundary_point(p))
	    {
		p->_nor[0] = 0.0;
		p->_nor[1] = 0.0;
		p->_nor[2] = 0.0;
		continue;
	    }

	    normal(p,hse,hs,nor,front);
	    curv = p->curvature;

	    pt = Coords(p);
	   
	    status = rect_in_which(pt,ict,top_grid);
	    if (!status) continue;
	    for(i = 0; i < dim; i++)
		ptt[i] = top_grid->L[i] + ict[i]*top_grid->h[i];

	    for(i = 0; i < dim; i++)
	    {
		d[i] = fabs(pt[i]-ptt[i])/rgr->h[i];
	        if(d[i] < -tol || d[i] > 1.0 + tol)
		{
		    status = NO;
		}
	    }
	    if (status == NO) continue;

	    if (dim == 1)
	    {
	    	for(i = 0; i < 2; i++)
		{
		    ix = ict[0] + i;
		    d1 = (i == 0) ? fabs(1.0-d[0]) : d[0];
		    f = d1;

		    index = d_index1d(ix,top_gmax);
		    gnor[0][index] += nor[0]*f;
		    norset[ix][0][0] = 1;
		}
	    }
	    else if (dim == 2)
	    {
	    	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
		    ix = ict[0] + i;
		    iy = ict[1] + j;
		    d1 = (i == 0) ? fabs(1.0-d[0]) : d[0];
		    d2 = (j == 0) ? fabs(1.0-d[1]) : d[1];
		    f = d1*d2;
		    index = d_index2d(ix,iy,top_gmax);
		    gnor[0][index] += nor[0]*f;
		    gnor[1][index] += nor[1]*f;
		    norset[ix][iy][0] = 1;
		}
	    }
	    else if (dim == 3)
	    {
	    	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		for (k = 0; k < 2; k++)
		{
		    ix = ict[0] + i;
		    iy = ict[1] + j;
		    iz = ict[2] + k;
		    d1 = (i == 0) ? fabs(1.0-d[0]) : d[0];
		    d2 = (j == 0) ? fabs(1.0-d[1]) : d[1];
		    d3 = (k == 0) ? fabs(1.0-d[2]) : d[2];
		    f = d1*d2*d3;
		    index = d_index3d(ix,iy,iz,top_gmax);
		    gnor[0][index] += nor[0]*f;
		    gnor[1][index] += nor[1]*f;
		    gnor[2][index] += nor[2]*f;
		    norset[ix][iy][iz] = 1;
		}
	    }
	}

	for (i = 0; i <= top_gmax[0]; i++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (k = 0; k <= top_gmax[2]; k++)
	{
	    //make sure Vel(st) is assigned 
	    if (norset[i][j][k] != 1)
		continue;

	    ic[0] = i;
	    ic[1] = j;
	    ic[2] = k;
	    index = d_index(ic,top_gmax,dim);
	    f = 0.0;
	    for(n=0; n<dim; n++)
		f += sqr(gnor[n][index]);
	    f = sqrt(f);

	    if (f < 1.0e-10)
	    {
		gnor[0][index] = 0.0;
		if (dim > 1) gnor[1][index] = 0.0;
		if (dim > 2) gnor[2][index] = 0.0;
	    }
	    else
	    {
		gnor[0][index] /= f;
		if (dim > 1) gnor[1][index] /= f;
		if (dim > 2) gnor[2][index] /= f;
	    }
	}

	found = YES;
	num = 1;
	while (found && num > 0)
	{
	    found = NO;
	    num = 0;

	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]; k++)
	    {
		if(norset[i][j][k] != 0)
		    continue;

		found = YES;
		ic[0] = i;
		ic[1] = j;
		ic[2] = k;

		if(get_ave_normal(ic,norset))
		{
		    num++;
		    norset[i][j][k] = 2;
		}
	    }

	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]; k++)
		if(norset[i][j][k] == 2)
		    norset[i][j][k] = 1;
	}
}

void G_CARTESIAN::tecplot_interior_states(
			char	*bname)
{
	char		s[1000];
	double		coords[3];
	int		ix, iy, iz, comp, i, imin[3], imax[3];
	double		**vel = eqn_params->vel;
	double		*dens = eqn_params->dens;
	double		*pres = eqn_params->pres;
	double		**gnor = eqn_params->gnor;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	FILE		*fp;
	int		index;

	sprintf(s,"%s-%d.plt", bname,pp_mynode());
	printf("tecplot_interior_states  file name %s \n",s);

	fp = fopen(s, "w");
	if(fp == NULL)
	{
	    printf("WARNING tecplot_interior_states, can not open file %s\n", s);
	    return; 
	}
	
	for(i=0; i<3; i++)
	{
	    imin[i] = 0;
	    imax[i] = top_gmax[i];
	}

	fprintf(fp, "TITLE = \"inner states\" ");
	if(dim == 2)
	{
	    fprintf(fp, "VARIABLES = \"x\", \"y\", \"comp\",  ");
	    fprintf(fp, "\"dens\", \"press\", \"u\", \"v\",  " );
	    fprintf(fp, "\"nx\", \"ny\",  " );
	    fprintf(fp, "\"dens1\", \"press1\", \"u1\", \"v1\",  " );
	    fprintf(fp, "\"dens2\", \"press2\", \"u2\", \"v2\"  \n" );
	}
	else
	{
	    fprintf(fp, "VARIABLES = \"x\", \"y\", \"z\", \"comp\",  ");
	    fprintf(fp, "\"dens\", \"press\", \"u\", \"v\", \"w\", " );
	    fprintf(fp, "\"nx\", \"ny\", \"nz\", " );
	    fprintf(fp, "\"dens1\", \"press1\", \"u1\", \"v1\", \"w1\" " );
	    fprintf(fp, "\"dens2\", \"press2\", \"u2\", \"v2\"  \"w2\"\n" );
	}

	if(dim == 2)
	    fprintf(fp, "ZONE i=%d, j=%d \n", imax[0]-imin[0]+1, imax[1]-imin[1]+1);
	else
	    fprintf(fp, "ZONE i=%d, j=%d, k=%d \n", imax[0]-imin[0]+1, imax[1]-imin[1]+1, imax[2]-imin[2]+1);

	if(dim == 2)
	{
	    for(iy=imin[1]; iy <= imax[1]; iy++)
		  for(ix=imin[0]; ix <= imax[0]; ix++)
		  {
			index = d_index2d(ix,iy,top_gmax);
			
			getRectangleCenter(index, coords);
			comp = cell_center[index].comp;

			fprintf(fp, "%f ", coords[0]);
			fprintf(fp, "%f ", coords[1]);
			fprintf(fp, "%d ", comp);
		
			if(!gas_comp(comp))
			{
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e   %12.5e %12.5e  ", 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ", 
			    		0.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ", 
			    		0.0, 0.0, 0.0, 0.0);
			}
			else
			{
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e   %12.5e %12.5e  ", 
					dens[index], pres[index], 
					vel[0][index],
					vel[1][index],
					gnor[0][index],
					gnor[1][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ",
			    		Gdens[0][index], Gpres[0][index],
					Gvel[0][0][index], Gvel[0][1][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ",
			    		Gdens[1][index], Gpres[1][index],
					Gvel[1][0][index], Gvel[1][1][index]);
			}

			fprintf(fp, "\n");
		  }
	}
	else
	{
	for(iz=imin[2]; iz <= imax[2]; iz++)
	    for(iy=imin[1]; iy <= imax[1]; iy++)
		  for(ix=imin[0]; ix <= imax[0]; ix++)
		  {
			index = d_index3d(ix,iy,iz,top_gmax);
			
			getRectangleCenter(index, coords);
			comp = cell_center[index].comp;

			fprintf(fp, "%f %f %f ", coords[0], coords[1], coords[2]);
			fprintf(fp, "%d ", comp);
		
			if(!gas_comp(comp))
			{
			    fprintf(fp, "%12.5e %12.5e  %12.5e %12.5e %12.5e  %12.5e %12.5e %12.5e ", 
					0.0,0.0, 0.0,0.0,0.0,  0.0,0.0,0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ", 
			    		0.0, 0.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ", 
			    		0.0, 0.0, 0.0, 0.0, 0.0);
			}
			else
			{
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e  %12.5e %12.5e %12.5e ", 
					dens[index], pres[index], 
					vel[0][index],
					vel[1][index],
					vel[2][index],
					gnor[0][index],
					gnor[1][index],
					gnor[2][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ",
			    		Gdens[0][index], Gpres[0][index],
					Gvel[0][0][index], Gvel[0][1][index], Gvel[0][2][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ",
			    		Gdens[1][index], Gpres[1][index],
					Gvel[1][0][index], Gvel[1][1][index], Gvel[1][2][index]);
			}
			fprintf(fp, "\n");
		  }
	}

	fclose(fp);

}

EXPORT  void    tecplot_surface_states(
	const char	*bname,
	FILE		*file,
	SURFACE		*s)
{
	TRI	*tri;
	POINT	*p;
	int	i,npts,ntri;
	Locstate  sl, sr;

	if (bname != NULL)//direct call
	{
	    if ((file = fopen(bname,"w")) == NULL)
	    {
		screen("WARNING in tecplot_surface(), "
		       "can't open %s\n",bname);
		return;
	    }
	    (void) fprintf(file,"TITLE = \"tecplot surface\"\n"
		   	    "VARIABLES = \"x\", \"y\", \"z\", \"PL\", \"PR\", \"DL\", \"DR\" "
			    "\"u\", \"v\", \"w\", \"u1\", \"v1\", \"w1\" \n");
	}
	
	//called from tecplot_interface
	if (file == NULL)
	{
	    screen("ERROR, in tecplot_surface, file is NULL\n");
	    clean_up(ERROR);
	}
	if (!(first_tri(s)))
	{
	    screen("WARNING, first bond of the curve is NULL\n");
	    return;
	}

	//count number of points(npts) and number of tris(ntri)
	for (tri=first_tri(s),ntri=0; !at_end_of_tri_list(tri,s); tri=tri->next,ntri++)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri=first_tri(s),npts=0; !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		}
	    }
	}
	//end counting
	
	fprint_wave_type(file, "ZONE T=\"", wave_type(s), "\"", s->interface);
    	fprintf(file, " N=%d E=%d\nF=FEPOINT, ET=TRIANGLE\n",npts,ntri);

	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri=first_tri(s),npts=0; !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		    
		    FT_GetStatesAtPoint(p,Hyper_surf_element(tri),Hyper_surf(s),&sl,&sr);
		    fprintf(file,"%15.8e %15.8e %15.8e  %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",Coords(p)[0],
		    	 Coords(p)[1],Coords(p)[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		}
	    }
	}
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		fprintf(file,"%d ",Index_of_point(Point_of_tri(tri)[i]));
	    }
	    fprintf(file,"\n");
	}

	if (ntri != s->num_tri)
	{
	    printf("WARNING, num of tri in surface is wrong\n"); 
	}
	if (bname != NULL)
	    fclose(file);

}	/* end tecplot_surface_states */


EXPORT  void    tecplot_interface_states(const char*, INTERFACE	*);

EXPORT  void    tecplot_interface_states(
	const char	*bname,
	INTERFACE	*intfc)
{
	SURFACE	**s;
	char    bname1[200];
	FILE	*file;

	sprintf(bname1, "%s.plt", bname);
	if ((file = fopen(bname1,"w")) == NULL)
	{
	    screen("WARNING in tecplot_interface_states(), "
	           "can't open %s\n",bname1);
	    return;
	}
	(void) fprintf(file,"TITLE = \"tecplot interface\"\n"
		   	    "VARIABLES = \"x\", \"y\", \"z\", \"PL\", \"PR\", \"DL\", \"DR\" "
			    "\"u\", \"v\", \"w\", \"u1\", \"v1\", \"w1\" \n");

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    tecplot_surface_states(NULL,file,*s);
	}
	fclose(file);
}	/* end tecplot_interface */

boolean G_CARTESIAN::get_ave_state(
	SWEEP 		m_vst,
	int		*ic,
	int		***norset,
	int		comp,
	int		ind)
{
	int		i, j, k, l, num, ic1[3], dir;
	double		gd, gp, gvel[3], gpd[2];
	boolean		found;
	double		**momn = m_vst.momn;
	double		*dens = m_vst.dens;
        double          **pdens = m_vst.pdens;
	double		*pres = m_vst.pres;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
        double          ***Gpdens = eqn_params->Gpdens;
	double		**Gpres = eqn_params->Gpres;
	int		index, index0;
	int		icoords[MAXD];

	found = NO;

	for(i=0; i<dim; i++)
	    for(j=0; j<2; j++)
	    {
		dir = j == 0 ? -1 : 1;
		ft_assign(ic1, ic, 3*INT);
		ic1[i] = ic[i] + dir;

		if(ic1[i] < 0 || ic1[i] > top_gmax[i])
		    continue;
		if(norset[ic1[0]][ic1[1]][ic1[2]] == 1)
		    found = YES;
	    }

	if(!found)
	    return NO;

	index0 = d_index(ic,top_gmax,dim);

	num = 0;
	gd = 0.0;
	gp = 0.0;
	gvel[0] = 0.0;
	gvel[1] = 0.0;
	gvel[2] = 0.0;

	for (i = ic[0]-1; i <= ic[0]+1; i++)
	for (j = ic[1]-1; j <= ic[1]+1; j++)
	for (k = ic[2]-1; k <= ic[2]+1; k++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    if(i < 0 || i > top_gmax[0] || 
	       j < 0 || j > top_gmax[1] ||
	       k < 0 || k > top_gmax[2])
		continue;
	    if(norset[i][j][k] != 1)
		continue;

	    index = d_index(icoords,top_gmax,dim);

	    if(cell_center[index].comp == comp)
	    {
		gd += dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        gpd[ii] += pdens[ii][index];
                    }
                }
		gp += pres[index];
		for (l = 0; l < dim; ++l)
		    gvel[l] += momn[l][index]/dens[index];
	    }
	    else
	    {
		gd += Gdens[ind][index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        gpd[ii] += Gpdens[ind][ii][index];
                    }
                }
		gp += Gpres[ind][index];
		for (l = 0; l < dim; ++l)
		    gvel[l] += Gvel[ind][l][index];
	    }
	    num++;
	}

	Gdens[ind][index0] = gd/num;
        if(eqn_params->multi_comp_non_reactive == YES)
        {
            int ii;
            for(ii = 0; ii < eqn_params->n_comps; ii++)
            {
                Gpdens[ind][ii][index0] = gpd[ii]/num;
            }
        }
	Gpres[ind][index0] = gp/num;
	for (l = 0; l < dim; ++l)
	    Gvel[ind][l][index0] = gvel[l]/num;

	return YES;
}

void G_CARTESIAN::get_ghost_state(
	SWEEP 		m_vst,
	int		comp,
	int		ind)
{
	int			i,j,k;
	int			ic[3],index;
	int			c, num;
	boolean			found;
	double			**momn = m_vst.momn;
	double			*dens = m_vst.dens;
        double                  **pdens = m_vst.pdens;
	double			*pres = m_vst.pres;
	double			***Gvel = eqn_params->Gvel;
	double			**Gdens = eqn_params->Gdens;
        double                  ***Gpdens = eqn_params->Gpdens;
	double			**Gpres = eqn_params->Gpres;
	static	int		***norset;
	static 	int 		loop_count = 0;
	std::list<ToFill> resetThese;
	std::list<ToFill> fillThese;


	if (norset == NULL)
	{
	    int ft_vec_size = 0;
	    for (i = 0; i < dim; ++i)
	    {
		if (top_gmax[i]+8 > ft_vec_size)
		    ft_vec_size = top_gmax[i]+8;
	    }
	    if(dim == 1)
	    	FT_TriArrayMemoryAlloc((POINTER*)&norset,ft_vec_size,
				   1,1,INT);

	    if(dim == 2)
	    	FT_TriArrayMemoryAlloc((POINTER*)&norset,ft_vec_size,
				   ft_vec_size,1,INT);
	    if(dim == 3)
	    	FT_TriArrayMemoryAlloc((POINTER*)&norset,ft_vec_size,
				   ft_vec_size,ft_vec_size,INT);
	}

	ToFill aghst;
	
	for (i=0; i<=top_gmax[0]; i++)
	for (j=0; j<=top_gmax[1]; j++)
	for (k=0; k<=top_gmax[2]; k++)
	{
	    if (dim == 1)
	    	index = d_index1d(i,top_gmax);
	    else if (dim == 2)
	    	index = d_index2d(i,j,top_gmax);
	    else if (dim == 3)
	    	index = d_index3d(i,j,k,top_gmax);
	    c = cell_center[index].comp;

	    // for each cell that has component "comp" we
	    // set G* values and mark norset  for that cell to 1
	    if(c == comp)
	    {
		norset[i][j][k] = 1;
		Gdens[ind][index] = dens[index];
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        Gpdens[ind][ii][index] = pdens[ii][index];
                    }
                }
		Gpres[ind][index] = pres[index];
		Gvel[ind][0][index] = momn[0][index]/dens[index];
		if(dim > 1)
		    Gvel[ind][1][index] = momn[1][index]/dens[index];
		if(dim > 2)
		    Gvel[ind][2][index] = momn[2][index]/dens[index];
	    }
	    else
	    {
		aghst.icoords[0] = i;
		aghst.icoords[1] = j;
		aghst.icoords[2] = k;
		// hardcoded the stencil size to 4.... 
		if(withinStencilLen(aghst.icoords, 3) )
		{
		    fillThese.push_back(aghst);
		}
		norset[i][j][k] = 0;
	    }
	}

	found = YES;
	num = 1;
	while(found && (num > 0))
	{
	    std::list<ToFill>::iterator it;

	    found = NO;
	    loop_count++;
	    num = 0;

	    resetThese.clear();	
	    for (it=fillThese.begin() ; it != fillThese.end(); )
	    {
		found = YES;
		ic[0] = it->icoords[0]; 
		ic[1] = it->icoords[1]; 
		ic[2] = it->icoords[2]; 

		// if no neighbors are 1, return 0.
		if(get_ave_state(m_vst, ic,norset,comp,ind))
		{
		    num++;
		    norset[ ic[0] ][ ic[1] ][ ic[2] ] = 2;
 		    aghst.icoords[0] = ic[0];
		    aghst.icoords[1] = ic[1];
		    aghst.icoords[2] = ic[2]; 
		    resetThese.push_back(aghst);
		    it=fillThese.erase(it);// erase returns the next valid entery after the one we just erased.
		}
		else
		{
		    ++it;
		}
	    }
	    for (it=resetThese.begin(); it != resetThese.end(); it++)
		 norset[it->icoords[0]][it->icoords[1]][it->icoords[2]] = 1;
	}
	fillThese.clear();
	resetThese.clear();	
	loop_count = 0;
}

void G_CARTESIAN::GFMGhostState_old(
	int	*ic,
	int	comp,
	STATE	*ghost_st)
{
	int		i, index;
	double		ncor;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	double		**Gnor = eqn_params->gnor;
	EOS_PARAMS	*eos = eqn_params->eos;

	index = d_index(ic,top_gmax,dim);
	ghost_st->eos = &(eos[comp]);
	ghost_st->dim = dim;

	ncor = 0.0;
	for(i=0; i<dim; i++)
	    ncor += (Gvel[1][i][index] - Gvel[0][i][index])*Gnor[i][index];
		    
	if(comp == 2)
	{
	    ghost_st->pres = Gpres[1][index];
	    ghost_st->dens = Gdens[0][index];
	    for(i=0; i<dim; i++)
		ghost_st->vel[i] = Gvel[0][i][index] + ncor*Gnor[i][index];
	}
	else
	{
	    ghost_st->pres = Gpres[0][index];
	    ghost_st->dens = Gdens[1][index];
	    for(i=0; i<dim; i++)
		ghost_st->vel[i] = Gvel[1][i][index] - ncor*Gnor[i][index];
	}
	for(i=0; i<dim; i++)
	    ghost_st->momn[i] = ghost_st->dens*ghost_st->vel[i];
	
	ghost_st->engy = EosEnergy(ghost_st);
}

void G_CARTESIAN::find_cross2d(double *c1, double *c2, double *cc, int idir, double *nor, int dir)
{
    	//double tol_new = top_h[idir]*1e-6;
        double tol = fabs(c1[idir]-c2[idir])*1e-6;

        //if(fabs(tol_new-tol)>1.e-12)
        //printf("HK tol_old %e tol_new %e top_h[idir] %e \n",fabs(c1[idir]-c2[idir])*1e-6,top_h[idir]*1e-6,top_h[idir]);

	if (idir!=0 && idir!=1)
	{
		printf("Wrong direction in find_cross().\n");
		clean_up(ERROR);
	}

	if ((idir==0 && fabs(c1[1]-c2[1])>1.e-12)
		|| (idir==1 && fabs(c1[0]-c2[0])>1.e-12))
	{
		printf("Error in find_cross(). Check coordinates.\n");
		clean_up(ERROR);
	}

	INTERFACE *intfc;
	POINT *ps, *pe;
	POINT *p1, *p2;
	CURVE **c;
	BOND *b;
	int cross = 0;
	int found = 0;
	double tmpc[2];
	double dist, min_dist;

	min_dist = HUGE_VAL;
	intfc = front->interf;

	FT_ScalarMemoryAlloc((POINTER*)&p1,sizeof(POINT));
	FT_ScalarMemoryAlloc((POINTER*)&p2,sizeof(POINT));

	for (c=intfc->curves; (!found)&&c&&*c; ++c)
		for (b=(*c)->first; b; b=b->next)
		{
			ps = b->start;
			pe = b->end;
			if (idir==0 && (ps->_coords[1]-c1[1])*(pe->_coords[1]-c1[1]) < tol)
			{
				tmpc[1] = c1[1];
				tmpc[0] = (pe->_coords[0]*(tmpc[1]-ps->_coords[1])+ps->_coords[0]*(pe->_coords[1]-tmpc[1]))/(pe->_coords[1]-ps->_coords[1]);
				if ((c1[0]-tmpc[0])*(c2[0]-tmpc[0]) < tol)
				    cross = 1;
				if (cross)
				{
				    if (dir == 1)
					dist = fabs(tmpc[0] - c1[0]);
				    else
					dist = fabs(c2[0] - tmpc[0]);
				    if (dist < min_dist)
				    {
					min_dist = dist;
					cc[0] = tmpc[0];
					cc[1] = tmpc[1];
					*p1 = *ps;
					*p2 = *pe;
					found = 1;
				    }
				    cross = 0;
				}
			}
			else if (idir==1 && (ps->_coords[0]-c1[0])*(pe->_coords[0]-c1[0]) < tol)
			{
				tmpc[0] = c1[0];
				tmpc[1] = (pe->_coords[1]*(tmpc[0]-ps->_coords[0])+ps->_coords[1]*(pe->_coords[0]-tmpc[0]))/(pe->_coords[0]-ps->_coords[0]);
				if ((c1[1]-tmpc[1])*(c2[1]-tmpc[1]) < tol)
					cross = 1;
				if (cross)
				{
				    if (dir == 1)
					dist = fabs(tmpc[1] - c1[1]);
				    else
					dist = fabs(c2[1] - tmpc[1]);
				    if (dist < min_dist)
				    {
					min_dist = dist;
					cc[0] = tmpc[0];
					cc[1] = tmpc[1];
					*p1 = *ps;
					*p2 = *pe;
					found = 1;
				    }
				    cross = 0;
				}
			}
		}

	if (!found)
	{
		printf("Can't find cross point.\n");
		clean_up(ERROR);
	}

	int			i;
	double			nor1[MAXD], nor2[MAXD], tmp;
	POINT			*p;
	HYPER_SURF		*hs;
	HYPER_SURF_ELEMENT	*hse;

	p = p1;
	hs = p1->hs;
	hse = p1->hse;
	normal(p,hse,hs,nor1,front);

	p = p2;
	hs = p2->hs;
	hse = p2->hse;
	normal(p,hse,hs,nor2,front);

	tmp = 0;
	for (i=0; i<dim; i++)
	{
		nor[i] = (p2->_coords[0]-cc[0])/(p2->_coords[0]-p1->_coords[0])*nor1[i]+(cc[0]-p1->_coords[0])/(p2->_coords[0]-p1->_coords[0])*nor2[i];
		tmp += nor[i]*nor[i];
	}
	for (i=0; i<dim; i++)
		nor[i] /= tmp;

	FT_FreeThese(2, p1, p2);
}

void G_CARTESIAN::find_cross3d(double *c1, double *c2, double *cc, int idir, double *nor, int dir)
{

	struct Table *T = table_of_interface(front->grid_intfc);
	TRI **t, *tri;
	int num_tris;
	int i, j, k0, k1, k2, n;
	int ic1[3], ic2[3];
	POINT **pt, *pp;
	double *pt0, *pt1, *pt2;
	double ps[3], pe[3], p[3];
	int found = 0;
	int in = 0;
	double dist, min_dist = 1e10;
	double f[3], pnor[3], len;

	if (!rect_in_which(c1,ic1,top_grid) ||
		!rect_in_which(c2,ic2,top_grid))
	{
		printf("Error in find_cross3d().\n");
		clean_up(ERROR);
	}
	if (c1[0]>c2[0] || c1[1]>c2[1] || c1[2]>c2[2])
	{
		printf("Wrong c1 or c2 in find_cross3d().\n");
		clean_up(ERROR);
	}

	for (k2=ic1[2]; k2<=ic2[2]; k2++)
	for (k1=ic1[1]; k1<=ic2[1]; k1++)
	for (k0=ic1[0]; k0<=ic2[0]; k0++)
	{
		t = T->tris[k2][k1][k0];
		num_tris = T->num_of_tris[k2][k1][k0];
		for (i=0; i<num_tris; i++)
		{
			pt = Point_of_tri(t[i]);
			pt0 = Coords(pt[0]);
			pt1 = Coords(pt[1]);
			pt2 = Coords(pt[2]);
			if (cross_line_segment_triangle(pt0,pt1,pt2,c1,c2,p))
			{
				found = 1;
				dist = 0;
				if (dir == 1)
					for (j=0; j<3; j++)
						dist += sqr(p[j]-c2[j]);
				else
					for (j=0; j<3; j++)
						dist += sqr(p[j]-c1[j]);
				if (dist < min_dist)
				{
					cc = p;
					tri = t[i];
					in = 1;
					min_dist = dist;
				}
			}
			else if (point_on_segment(pt0,c1,c2))
			{
				found = 1;
				dist = 0;
				if (dir == 1)
					for (j=0; j<3; j++)
						dist += sqr(pt0[j]-c2[j]);
				else
					for (j=0; j<3; j++)
						dist += sqr(pt0[j]-c1[j]);
				if (dist < min_dist)
				{
					cc = pt0;
					tri = t[i];
					n = 0;
					in = 0;
					min_dist = dist;
				}
			}
			else if (point_on_segment(pt1,c1,c2))
			{
				found = 1;
				dist = 0;
				if (dir == 1)
					for (j=0; j<3; j++)
						dist += sqr(pt1[j]-c2[j]);
				else
					for (j=0; j<3; j++)
						dist += sqr(pt1[j]-c1[j]);
				if (dist < min_dist)
				{
					cc = pt1;
					tri = t[i];
					n = 1;
					in = 0;
					min_dist = dist;
				}
			}
			else if (point_on_segment(pt2,c1,c2))
			{
				found = 1;
				dist = 0;
				if (dir == 1)
					for (j=0; j<3; j++)
						dist += sqr(pt2[j]-c2[j]);
				else
					for (j=0; j<3; j++)
						dist += sqr(pt2[j]-c1[j]);
				if (dist < min_dist)
				{
					cc = pt2;
					tri = t[i];
					n = 2;
					in = 0;
					min_dist = dist;
				}
			}
		}
	}

	if (found == 1 && in == 1)
	{
		pt = Point_of_tri(tri);
		pt0 = Coords(pt[0]);
		pt1 = Coords(pt[1]);
		pt2 = Coords(pt[2]);

		linear_interp_coefs_three_pts(f,cc,pt0,pt1,pt2);

		for(j=0; j<3; j++)
	    		nor[j] = 0.0;
	
		for(i=0; i<3; i++)
		{
			if(f[i] < 0.0)
				f[i] = 0.0;
			if(f[i] > 1.0)
				f[i] = 1.0;

	    		pp = Point_of_tri(tri)[i];
			normal(pp, Hyper_surf_element(tri), pp->hs, pnor, front);
	    
			for(j=0; j<3; j++)
				nor[j] += pnor[j]*f[i];
	    
		}

		len = Mag3d(nor);
		for(j=0; j<3; j++)
	    		nor[j] /= len;
	}
	else if (found == 1 && in == 0)
	{
		pp = Point_of_tri(tri)[n];
		normal(pp, Hyper_surf_element(tri), pp->hs, nor, front);
	}
	else
	{
		printf("Can't find crossing point.\n");
		clean_up(ERROR);
	}
}

boolean G_CARTESIAN::point_on_segment(
	double		*p,
	double		*ps,
	double		*pe)
{
	int		i;
	double		v1[3], v2[3], nv[3];
	double		tol = 1e-5;

	for (i=0; i<3; i++)
	{
		v1[i] = p[i] - ps[i];
		v2[i] = pe[i] - p[i];
	}

	Cross3d(v1,v2,nv);

	if (fabs(nv[0]) < tol && fabs(nv[1]) < tol && fabs(nv[2]) < tol)
		return YES;

	else
		return NO;
}


boolean G_CARTESIAN::cross_line_segment_triangle(
        double           *p1t1,
        double           *p1t2,
        double           *p1t3,
        double           *p2s,
        double           *p2e,
        double           *p)
{
        int             i, max_index;
        bool            intersect_with_plane;
        double          v1[3],v2[3],v3[3],v4[2],v5[2],v6[2];
        double          det1,det2;
        double          tmp1[3],tmp2[3],tmp3[3];
        double          p1t1_tmp[2],p1t2_tmp[2],p1t3_tmp[2],p_tmp[2];
        double          n1,n2,n3,n4,t,nv1[3],nv2[3],norm[3];


        for(i = 0; i < 3; i++)
        {
            v1[i] = p2s[i] - p1t1[i];
            v2[i] = p1t2[i] - p1t1[i];
            v3[i] = p1t3[i] - p1t1[i];
        }

        det1 = Det3d(v1,v2,v3);

        for(i = 0; i < 3; i++)
            v1[i] = p2e[i] - p1t1[i];

        det2 = Det3d(v1,v2,v3);

        if(det1 == 0 && det2 == 0)
            intersect_with_plane = NO;
        else if(det1*det2 > 0)
            intersect_with_plane = NO;
        else
            intersect_with_plane = YES;

        if(intersect_with_plane != YES)
            return NO;

        tmp1[0] = p1t1[1];  tmp1[1] = p1t1[2];  tmp1[2] = 1.0;
        tmp2[0] = p1t2[1];  tmp2[1] = p1t2[2];  tmp2[2] = 1.0;
        tmp3[0] = p1t3[1];  tmp3[1] = p1t3[2];  tmp3[2] = 1.0;
        n1 = Det3d(tmp1,tmp2,tmp3);

        tmp1[0] = p1t1[0];  tmp1[1] = p1t1[2];  tmp1[2] = 1.0;
        tmp2[0] = p1t2[0];  tmp2[1] = p1t2[2];  tmp2[2] = 1.0;
        tmp3[0] = p1t3[0];  tmp3[1] = p1t3[2];  tmp3[2] = 1.0;
        n2 = Det3d(tmp1,tmp2,tmp3);

        tmp1[0] = p1t1[0];  tmp1[1] = p1t1[1];  tmp1[2] = 1.0;
        tmp2[0] = p1t2[0];  tmp2[1] = p1t2[1];  tmp2[2] = 1.0;
        tmp3[0] = p1t3[0];  tmp3[1] = p1t3[1];  tmp3[2] = 1.0;
        n3 = Det3d(tmp1,tmp2,tmp3);

        n4 = Det3d(p1t1,p1t2,p1t3);

        t = (p2s[1]*n2-p2s[0]*n1-p2s[2]*n3+n4)
            /((p2e[0]-p2s[0])*n1-(p2e[1]-p2s[1])*n2+(p2e[2]-p2s[2])*n3);

        p[0] = p2s[0] + (p2e[0] - p2s[0])*t;
        p[1] = p2s[1] + (p2e[1] - p2s[1])*t;
        p[2] = p2s[2] + (p2e[2] - p2s[2])*t;

        for(i = 0; i < 3; i++)
        {
            nv1[i] = p1t2[i] - p1t1[i];
            nv2[i] = p1t3[i] - p1t1[i];
        }

        Cross3d(nv1,nv2,norm);

        if(within_tri(p,p1t1,p1t2,p1t3,norm,1.0e-6) != YES)
            return NO;

        return YES;
}

void G_CARTESIAN::GFMGhostState(
	int	*real_ic,
	int	*ic,
	int	comp,
	STATE	*ghost_st,
	int	dir,	//dir == 1 or -1
        SWEEP	*m_vst,
	int	idir)	//idir == 0 or 1
{
	int		i, index_r, index_g;
        double          **Rmomn = m_vst->momn;
        double          *Rdens = m_vst->dens;
        double          **Rpdens = m_vst->pdens;
        double          *Rpres = m_vst->pres;
	EOS_PARAMS	*eos = eqn_params->eos;
	double		lvnor, rvnor;
	double		lvtang[MAXD], rvtang[MAXD];
	int		li, ri, lcomp, rcomp;
	double		nor[MAXD];
	double		c1[MAXD], c2[MAXD], cc[MAXD];

	index_r = d_index(real_ic,top_gmax,dim);
	index_g = d_index(ic,top_gmax,dim);
        /*
	if (dim==2)
	{
	    if (dir==1)
		for (i=0; i<dim; i++)
		{
		    c1[i] = cell_center[index_r].m_coords[i];
		    c2[i] = cell_center[index_g].m_coords[i];
		}
	    else
		for (i=0; i<dim; i++)
		{
		    c1[i] = cell_center[index_g].m_coords[i];
		    c2[i] = cell_center[index_r].m_coords[i];
		}
	    find_cross2d(c1,c2,cc,idir,nor,dir);
	}
	else if (dim==3)
	{
	    for (i=0; i<dim; i++)
		nor[i] = eqn_params->gnor[i][index_g];
	    if (dir==1)
		for (i=0; i<dim; i++)
		{
		    c1[i] = cell_center[index_r].m_coords[i];
		    c2[i] = cell_center[index_g].m_coords[i];
		}
	    else
		for (i=0; i<dim; i++)
		{
		    c1[i] = cell_center[index_g].m_coords[i];
		    c2[i] = cell_center[index_r].m_coords[i];
		}
		find_cross3d(c1,c2,cc,idir,nor,dir);
	}
        */

        for (i = 0; i < dim; i++)
            nor[i] = 0.0;
        nor[idir] = 1.0;

	ghost_st->eos = &(eos[comp]);
	ghost_st->dim = dim;
	/*
	if (nor[idir]<0)
	    for (i=0; i<dim; i++)
		nor[i] = -nor[i];
        */
	const double epsilon = 1.e-10;
	const double delta = 1.e-10;

	if (dir==1)
	{
	    li = index_r;
	    ri = index_g;
	}
	else
	{
	    li = index_g;
	    ri = index_r;
	}

	double pml, pmr, uml, umr, ml, mr;
	RIEMANN_SOLVER_WAVE_TYPE l_wave,r_wave;
	STATE *stl2, *str2;
	STATE *ansl, *ansr;
	int tmpc;

	FT_ScalarMemoryAlloc((POINTER*)&stl2,sizeof(STATE));
	FT_ScalarMemoryAlloc((POINTER*)&str2,sizeof(STATE));

	if (dir == 1)
	{
	    str2->eos = &(eos[comp]);
	    tmpc = (comp == GAS_COMP1) ? GAS_COMP2 : GAS_COMP1;
	    stl2->eos = &(eos[tmpc]);
	}
	else if (dir == -1)
	{
	    stl2->eos = &(eos[comp]);
	    tmpc = (comp == GAS_COMP1) ? GAS_COMP2 : GAS_COMP1;
	    str2->eos = &(eos[tmpc]);
	}
	else
	{
	    printf("wrong dir in GFMGhostState().\n");
	    clean_up(ERROR);
	}
	stl2->dim = dim;
	str2->dim = dim;
        stl2->dens = Rdens[li];
        str2->dens = Rdens[ri];
        stl2->pres = Rpres[li];
        str2->pres = Rpres[ri];
        stl2->pdens[0] = Rpdens[0][li];
        stl2->pdens[1] = Rpdens[1][li];
        str2->pdens[0] = Rpdens[0][ri];
        str2->pdens[1] = Rpdens[1][ri];

	lvnor = 0;
	rvnor = 0;
	for (i=0; i<dim; i++)
	{
	    lvnor += (Rmomn[i][li]/Rdens[li])*nor[i];
	    rvnor += (Rmomn[i][ri]/Rdens[ri])*nor[i];
	}
	for (i=0; i<dim; i++)
	{
	    lvtang[i] = (Rmomn[i][li]/Rdens[li]) - lvnor*nor[i];
	    rvtang[i] = (Rmomn[i][ri]/Rdens[ri]) - rvnor*nor[i];
	}
	stl2->vel[0] = lvnor;
	str2->vel[0] = rvnor;
	for (i = 1; i < dim; i++)
	    stl2->vel[i] = str2->vel[i] = 0.0;

	FT_ScalarMemoryAlloc((POINTER*)&ansl,sizeof(STATE));
	FT_ScalarMemoryAlloc((POINTER*)&ansr,sizeof(STATE));

	find_mid_state(stl2,str2,0.0/*pjump*/,&pml,&pmr,&uml,&umr,&ml,&mr,&l_wave,&r_wave);
	midstate(stl2,ansl,ml,uml,pml,l_wave,1);
	midstate(str2,ansr,mr,umr,pmr,r_wave,-1);

	if (dir==1)
	{
	    ghost_st->pres = ansl->pres;
	    ghost_st->dens = ansl->dens;
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                int ii;
                for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    //ghost_st->pdens[ii] = (Rpdens[ii][li]/Rdens[li])*ghost_st->dens;
                    ghost_st->pdens[ii] = ansl->pdens[ii];
                }
            }
	    for (i=0; i<dim; i++)
	    {
		ghost_st->vel[i] = ansl->vel[0]*nor[i] + rvtang[i];
		ghost_st->momn[i] = ghost_st->dens*ghost_st->vel[i];
	    }
	}
	else
	{
	    ghost_st->pres = ansr->pres;
	    ghost_st->dens = ansr->dens;
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                int ii;
                for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    //ghost_st->pdens[ii] = (Rpdens[ii][ri]/Rdens[ri])*ghost_st->dens;
                    ghost_st->pdens[ii] = ansr->pdens[ii];
                }
            }
	    for (i=0; i<dim; i++)
	    {
		ghost_st->vel[i] = ansr->vel[0]*nor[i] + lvtang[i];
		ghost_st->momn[i] = ghost_st->dens*ghost_st->vel[i];
	    }
	}
	ghost_st->engy = EosEnergy(ghost_st);

	FT_FreeThese(2,stl2,str2);
	FT_FreeThese(2,ansl,ansr);
}

void G_CARTESIAN::setNeumannStates(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
	int 		i,j,index;
	int             ind2[2][2] = {{0,1},{1,0}};
        int             ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int 		ic[MAXD];
	double		*vel_ref = state->vel;
	double		coords[MAXD],coords_ref[MAXD],crx_coords[MAXD];
	double		nor[MAXD],vn,v[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	GRID_DIRECTION  dir;
	STATE		st_tmp;

	st_tmp.eos = state->eos;
	st_tmp.dim = dim;
	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic[i] = icoords[i];
	}
	dir = (nb == 0) ? ldir[idir] : rdir[idir];
	FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);

	if (debugging("neumann_buffer"))
	{
	    (void) printf("Entering setNeumannStates()\n");
	    (void) printf("comp = %d\n",comp);
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],
				icoords[2]);
	    (void) printf("idir = %d nb = %d\n",idir,nb);
	    (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
	    (void) print_general_vector("coords = ",coords,dim,"\n");
	    (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
	    (void) print_general_vector("nor = ",nor,dim,"\n");
	    (void) print_general_vector("vel_ref = ",vel_ref,dim,"\n");
	}

	for (i = istart; i <= nrad; ++i)
	{
	    /* Find ghost point */
	    ic[idir] = (nb == 0) ? icoords[idir] - i : icoords[idir] + i;
	    for (j = 0; j < dim; ++j)
		coords_ref[j] = top_L[j] + ic[j]*top_h[j];

	    /* Reflect ghost point through intfc-mirror at crossing */
	    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ref[idir];
	    vn = 0.0;
	    for (j = 0; j < dim; ++j)
	    {
		v[j] = coords_ref[j] - crx_coords[j];
		vn += v[j]*nor[j];
	    }
	    for (j = 0; j < dim; ++j)
		v[j] = 2.0*vn*nor[j] - v[j];
	    for (j = 0; j < dim; ++j)
		coords_ref[j] = crx_coords[j] + v[j];
			
	    /* Interpolate the state at the reflected point */
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->dens,getStateDens,&st_tmp.dens,&m_vst->dens[index]);

            if(eqn_params->multi_comp_non_reactive == YES)
            {
                //int ii;
                //for(ii = 0; ii < eqn_params->n_comps; ii++)
                {
                    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                        m_vst->pdens[0],getStatePdens0,&st_tmp.pdens[0],&m_vst->pdens[0][index]);
                    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                        m_vst->pdens[1],getStatePdens1,&st_tmp.pdens[1],&m_vst->pdens[1][index]);
                }
            }

	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->engy,getStateEngy,&st_tmp.engy,&m_vst->engy[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->pres,getStatePres,&st_tmp.pres,&m_vst->pres[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
			m_vst->momn[0],getStateXmom,&st_tmp.momn[0],
			&m_vst->momn[0][index]);
	    if (dim > 1)
		FT_IntrpStateVarAtCoords(front,comp,coords_ref,
			m_vst->momn[1],getStateYmom,&st_tmp.momn[1],
			&m_vst->momn[1][index]);
	    if (dim > 2)
		FT_IntrpStateVarAtCoords(front,comp,coords_ref,
			m_vst->momn[2],getStateZmom,&st_tmp.momn[2],
			&m_vst->momn[2][index]);
		/* Galileo Transformation */
	    vn = 0.0;
	    for (j = 0; j < dim; j++)
	    {
		v[j] = st_tmp.momn[j]/st_tmp.dens - vel_ref[j];
		vn += v[j]*nor[j];
	    }
	    for (j = 0; j < dim; j++)
	    {
		v[j] += vel_ref[j] - 2.0*vn*nor[j];
		st_tmp.momn[j] = v[j]*st_tmp.dens;
	    }
	    st_tmp.pres = EosPressure(&st_tmp);
	    if (st_tmp.pres < min_pres) st_tmp.pres = min_pres;
	    st_tmp.engy = EosEnergy(&st_tmp);

	    if (nb == 0)
	    {
		vst->dens[nrad-i] = st_tmp.dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        vst->pdens[ii][nrad-i] = st_tmp.pdens[ii];
                    }
                }
		vst->engy[nrad-i] = st_tmp.engy;
		vst->pres[nrad-i] = st_tmp.pres;
	    	for (j = 0; j < 3; j++)
		    vst->momn[j][nrad-i] = 0.0;
		if (dim == 1)
		   vst->momn[0][nrad-i] = st_tmp.momn[0];
	    	else if (dim == 2)
		    for (j = 0; j < 2; j++)
		    	vst->momn[j][nrad-i] = 
				st_tmp.momn[ind2[idir][j]];
	    	else if (dim == 3)
		    for (j = 0; j < 3; j++)
		    	vst->momn[j][nrad-i] = 
				st_tmp.momn[ind3[idir][j]];
	    }
	    else
	    {
		vst->dens[n+nrad+i-1] = st_tmp.dens;
                if(eqn_params->multi_comp_non_reactive == YES)
                {
                    int ii;
                    for(ii = 0; ii < eqn_params->n_comps; ii++)
                    {
                        vst->pdens[ii][n+nrad+i-1] = st_tmp.pdens[ii];
                    }
                }
		vst->engy[n+nrad+i-1] = st_tmp.engy;
		vst->pres[n+nrad+i-1] = st_tmp.pres;
	    	for (j = 0; j < 3; j++)
		    vst->momn[j][n+nrad+i-1] = 0.0;
		if (dim == 1)
		   vst->momn[0][n+nrad+i-1] = st_tmp.momn[0];
	    	else if (dim == 2)
		    for (j = 0; j < 2; j++)
		    	vst->momn[j][n+nrad+i-1] = 
				st_tmp.momn[ind2[idir][j]];
	    	else if (dim == 3)
		    for (j = 0; j < 3; j++)
		    	vst->momn[j][n+nrad+i-1] = 
				st_tmp.momn[ind3[idir][j]];
	    }
	}
	if (debugging("neumann_buffer"))
	    (void) printf("Leaving setNeumannStates()\n");
}	/* end setNeumannStates */

void G_CARTESIAN::setDirichletStates(
	STATE		*crx_state,
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	int		*icoords,
	int		dir,
	int		nb,
	int		n,
	int		istart)
{
	int		j, k, index;
	STATE 		*state;
	int		ind2[2][2] = {{0,1},{1,0}};
	int		ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};

	if (nb == 0)
	{
	  if (boundary_state(hs) != NULL)
	  {
	    //preset state bdry
	    state = (STATE*)boundary_state(hs);
	    for (k = istart; k <= nrad; ++k)
	    {
		vst->dens[nrad-k] = state->dens;
		vst->engy[nrad-k] = state->engy;
		vst->pres[nrad-k] = state->pres;
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-k] = 0.0;
		if (dim == 1)
		    vst->momn[0][nrad-k] = state->momn[0];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][nrad-k] = state->momn[ind2[dir][j]];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][nrad-k] = state->momn[ind3[dir][j]];
	    }
	  }
	  else if (boundary_state_function(hs) &&
              strcmp(boundary_state_function_name(hs),
	      "cF_flowThroughBoundaryState") == 0)
	  {
	    //flow through bdry
	    for (k = istart; k <= nrad; ++k)
	    {
		index = d_index(icoords,top_gmax, dim);
		vst->dens[nrad-k] = m_vst->dens[index];
		vst->engy[nrad-k] = m_vst->engy[index];
		vst->pres[nrad-k] = m_vst->pres[index];
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-k] = 0.0;
		if (dim == 1)
		    vst->momn[0][nrad-k] = m_vst->momn[0][index];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][nrad-k] = m_vst->momn[ind2[dir][j]][index];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][nrad-k] = m_vst->momn[ind3[dir][j]][index];
	    }
	  }
	  else
	  {
	    (void) printf("Unimplemented Dirichlet boundary type!\n");
	    clean_up(ERROR);
	  }
	}
	else
	{
	  if (boundary_state(hs) != NULL)
	  {
	    state = (STATE*)boundary_state(hs);
	    for (k = istart; k <= nrad; ++k)
	    {
		vst->dens[n+nrad+k-1] = state->dens;
		vst->engy[n+nrad+k-1] = state->engy;
		vst->pres[n+nrad+k-1] = state->pres;
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+k-1] = 0.0;
		if (dim == 1)
		    vst->momn[0][n+nrad+k-1] = state->momn[0];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][n+nrad+k-1] = state->momn[ind2[dir][j]];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][n+nrad+k-1] = state->momn[ind3[dir][j]];
	    }
	  }
	  else if (boundary_state_function(hs) &&
              strcmp(boundary_state_function_name(hs),
	      "cF_flowThroughBoundaryState") == 0)
	  {
	    for (k = istart; k <= nrad; ++k)
	    {
		index = d_index(icoords,top_gmax, dim);
		vst->dens[n+nrad+k-1] = m_vst->dens[index];
		vst->engy[n+nrad+k-1] = m_vst->engy[index];
		vst->pres[n+nrad+k-1] = m_vst->pres[index];
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+k-1] = 0.0;
		if (dim == 1)
		    vst->momn[0][n+nrad+k-1] = m_vst->momn[0][index];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][n+nrad+k-1] = 
					m_vst->momn[ind2[dir][j]][index];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][n+nrad+k-1] = 
					m_vst->momn[ind3[dir][j]][index];
	    }
	  }
	  else
	  {
	    (void) printf("Unimplemented Dirichlet boundary type!\n");
	    clean_up(ERROR);
	  }
	}
}

void G_CARTESIAN::initSampleVelocity(char *in_name)
{
        FILE *infile;
	static SAMPLE *sample;
	char *sample_type;
	double *sample_line;

	infile = fopen(in_name,"r");
	FT_ScalarMemoryAlloc((POINTER*)&sample,sizeof(SAMPLE));
	sample_type = sample->sample_type;
	sample_line = sample->sample_coords;
	dim = front->rect_grid->dim;

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
            (void) printf(" %f %f\n",sample_line[0],sample_line[1]);
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

void G_CARTESIAN::checkCorrectForTolerance(STATE *state)
{
	if (state->dens < min_dens)
	    state->dens = min_dens;
	if (state->pres < min_pres)
	    state->pres = min_pres;
	state->engy = EosEnergy(state);
}	/* end checkCorrectForTolerance */

boolean G_CARTESIAN::needBufferFromIntfc(
	COMPONENT domain_comp,
	COMPONENT comp)
{
	if (eqn_params->tracked)
	    return (domain_comp != comp) ? YES : NO;
	else
	    return (gas_comp(comp)) ? NO : YES;
}	/* needBufferFromIntfc */


bool G_CARTESIAN::withinStencilLen( int *icrds, int stencil )
{
        int istart = std::max(0, icrds[0] - stencil);
        int jstart = std::max(0, icrds[1] - stencil);
        int kstart = std::max(0, icrds[2] - stencil);
        int iend = std::min(top_gmax[0],icrds[0]+stencil);
        int jend = std::min(top_gmax[1],icrds[1]+stencil);
        int kend = std::min(top_gmax[2],icrds[2]+stencil);

        int index  =  d_index(icrds,top_gmax,dim);
        int mycomp =  top_comp[index];

        int i,j,k;
        if (dim == 2)
        {
            for (i = istart; i <= iend; i++)
            for (j = jstart; j <= jend; j++)
            {
                int ic[3];
                ic[0] = i; ic[1] = j; ic[2] = 0;
                index  =  d_index(ic,top_gmax,dim);
                if(mycomp != top_comp[index])
                    return true;
            }
            return NO;
        }
        else if (dim == 3)
	{
            for (i = istart; i <= iend; i++)
            for (j = jstart; j <= jend; j++)
            for (k = kstart; k <= kend; k++)
	    {
                int ic[3];
                ic[0] = i; ic[1] = j; ic[2] = k;
                index  =  d_index(ic,top_gmax,dim);
                if(mycomp != top_comp[index])
                    return true;
	    }
            return NO;
	}
	return YES;
}

void G_CARTESIAN::checkIntfc(char *out_name)
{
	INTERFACE *intfc = front->interf;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double h_min, h_max;
	double tol = 1e-5;

	if (dim == 1)
		return;

    char filename[500];
    FILE *outfile;


    if(pp_mynode() == 0)
    {
        sprintf(filename, "%s/growth_rate.dat",out_name);
        if(front->step == 1)
            outfile = fopen(filename,"w");
        else
            outfile = fopen(filename,"a");
    }

	h_min = front->rect_grid->U[dim-1];
	h_max = front->rect_grid->L[dim-1];

       	next_point(intfc,NULL,NULL,NULL);
       	while (next_point(intfc,&p,&hse,&hs))
	{
		if (Coords(p)[dim-1] > front->rect_grid->U[dim-1]-tol ||
		    Coords(p)[dim-1] < front->rect_grid->L[dim-1]+tol)
			continue;
		if (Coords(p)[dim-1] > h_max)
			h_max = Coords(p)[dim-1];
		if (Coords(p)[dim-1] < h_min)
			h_min = Coords(p)[dim-1];
	}

    pp_global_max(&h_max,1);
    pp_global_min(&h_min,1);

    if(pp_mynode() == 0)
    {
        if(front->step == 1)
        {
            (void) fprintf(outfile,"#time    h_max    h_min\n");
        }
        (void) fprintf(outfile,"%.5E     ",front->time);
        (void) fprintf(outfile,"%.5E     ",h_max);
        (void) fprintf(outfile,"%.5E     ",h_min);
        (void) fprintf(outfile,"\n");
        fclose(outfile);
    }

}

void G_CARTESIAN::record_intfc_extrema()
{
	INTERFACE		*intfc = front->interf;
	EQN_PARAMS		*eqn_params = (EQN_PARAMS*)front->extra1;
        POINT			*p;
        HYPER_SURF		*hs;
        HYPER_SURF_ELEMENT	*hse;
	COMPONENT		comp_at_min, comp_at_max;
	int			i;
	double			h_at_min, h_at_max;
	double			g_h_at_min, g_h_at_max;
	double			crds_at_min[MAXD], crds_at_max[MAXD];
	double			tol = 1e-5;
	Big_State		tmp_bst[2];
	STATE			tmp_st;

	iext = eqn_params->iext;

	comp_at_max = GAS_COMP1;
	comp_at_min = GAS_COMP2;

	/*intfc extrema*/
	if (dim == 1)
		return;
	h_at_min = front->rect_grid->GU[dim-1];
	h_at_max = front->rect_grid->GL[dim-1];

       	next_point(intfc,NULL,NULL,NULL);
       	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (dim == 2 &&
		Coords(p)[0] < front->rect_grid->GL[0]+tol ||
		Coords(p)[0] > front->rect_grid->GU[0]-tol)
		continue;
	    if (Coords(p)[dim-1] > front->rect_grid->GU[dim-1]-tol ||
		Coords(p)[dim-1] < front->rect_grid->GL[dim-1]+tol)
		continue;
	    if (Coords(p)[dim-1] > h_at_max)
	    {
		h_at_max = Coords(p)[dim-1];
		for (i = 0; i < dim; i++)
		    crds_at_max[i] = Coords(p)[i];
	    }
	    if (Coords(p)[dim-1] < h_at_min)
	    {
		h_at_min = Coords(p)[dim-1];
		for (i = 0; i < dim; i++)
		    crds_at_min[i] = Coords(p)[i];
	    }
	}

	g_h_at_max = h_at_max;
	g_h_at_min = h_at_min;

	pp_global_max(&g_h_at_max,1L);
	pp_global_min(&g_h_at_min,1L);

	iext->h_max = g_h_at_max;
	iext->h_min = g_h_at_min;

	zero_scalar(&bst_00[0],sizeof(Big_State));
	zero_scalar(&bst_00[1],sizeof(Big_State));

	double crds_at_g_min[MAXD], crds_at_g_max[MAXD];
	for (i = 0; i < dim; i++)
	{
	    crds_at_g_min[i] = 0.0;
	    crds_at_g_max[i] = 0.0;
	}

	if (fabs(h_at_min-g_h_at_min) < 1e-12)
	    for (i = 0; i < dim; i++)
		crds_at_g_min[i] = crds_at_min[i];
	for (i = 0; i < dim; i++)
	    pp_global_max(&(crds_at_g_min[i]),1L);

	hyp_sol(crds_at_g_min,&tmp_st,comp_at_min);

	bst_00[0].d = tmp_st.dens;
	bst_00[0].p = tmp_st.pres;
	bst_00[0].e = tmp_st.engy;
	pp_global_max(&(bst_00[0].d),1L);
	pp_global_max(&(bst_00[0].p),1L);
	pp_global_max(&(bst_00[0].e),1L);
	bst_00[0].v[dim-1] = tmp_st.momn[dim-1] / bst_00[0].d;
	double tmp_momn[MAXD];
	for (i = 0; i < dim; i++)
	    tmp_momn[i] = tmp_st.momn[i];
	for (i = 0; i < dim; i++)
	    pp_global_max(&(tmp_momn[i]),1L);
	if (tmp_momn[dim-1] < 1e-12)
	    for (i = 0; i < dim; i++)
		pp_global_min(&(bst_00[0].v[i]),1L);
	else
	    for (i = 0; i < dim; i++)
		pp_global_max(&(bst_00[0].v[i]),1L);

	if (fabs(h_at_max-g_h_at_max) < 1e-12)
	    for (i = 0; i < dim; i++)
		crds_at_g_max[i] = crds_at_max[i];
	for (i = 0; i < dim; i++)
	    pp_global_max(&(crds_at_g_max[i]),1L);

	hyp_sol(crds_at_g_max,&tmp_st,comp_at_max);

	bst_00[1].d = tmp_st.dens;
	bst_00[1].p = tmp_st.pres;
	bst_00[1].e = tmp_st.engy;
	pp_global_max(&(bst_00[1].d),1L);
	pp_global_max(&(bst_00[1].p),1L);
	pp_global_max(&(bst_00[1].e),1L);
	bst_00[1].v[dim-1] = tmp_st.momn[dim-1] / bst_00[1].d;
	for (i = 0; i < dim; i++)
	    tmp_momn[i] = tmp_st.momn[i];
	for (i = 0; i < dim; i++)
	    pp_global_max(&(tmp_momn[i]),1L);
	if (tmp_momn[dim-1] < 1e-12)
	    for (i = 0; i < dim; i++)
		pp_global_min(&(bst_00[1].v[i]),1L);
	else
	    for (i = 0; i < dim; i++)
		pp_global_max(&(bst_00[1].v[i]),1L);

	accumulate_state_in_layer(g_h_at_max,&amb_bst_00[1],YES);
	accumulate_state_in_layer(g_h_at_min,&amb_bst_00[0],YES);

	/*volume fraction levels*/
	double dh = (g_h_at_max - g_h_at_min)/100.0;

	if (iext->do_01 == YES)
	{
	    iext->h_max_01 = height_at_fraction(g_h_at_max,dh,0.01,-1,comp_at_max);
	    iext->h_min_01 = height_at_fraction(g_h_at_min,dh,0.01,1,comp_at_min);

	    accumulate_state_in_layer(iext->h_max_01,tmp_bst,NO);
	    int index_at_max = 0;
	    copy_Big_State(&tmp_bst[index_at_max],&bst_01[1]);
	    copy_Big_State(&tmp_bst[1-index_at_max],&amb_bst_01[1]);

	    accumulate_state_in_layer(iext->h_min_01,tmp_bst,NO);
	    int index_at_min = 1;
	    copy_Big_State(&tmp_bst[index_at_min],&bst_01[0]);
	    copy_Big_State(&tmp_bst[1-index_at_min],&amb_bst_01[0]);

	}

	if (iext->do_05 == YES)
	{
	    iext->h_max_05 = height_at_fraction(g_h_at_max,dh,0.05,-1,comp_at_max);
	    iext->h_min_05 = height_at_fraction(g_h_at_min,dh,0.05,1,comp_at_min);

	    accumulate_state_in_layer(iext->h_max_05,tmp_bst,NO);
	    int index_at_max = 0;
	    copy_Big_State(&tmp_bst[index_at_max],&bst_05[1]);
	    copy_Big_State(&tmp_bst[1-index_at_max],&amb_bst_05[1]);

	    accumulate_state_in_layer(iext->h_min_05,tmp_bst,NO);
	    int index_at_min = 1;
	    copy_Big_State(&tmp_bst[index_at_min],&bst_05[0]);
	    copy_Big_State(&tmp_bst[1-index_at_min],&amb_bst_05[0]);

	}
}

double G_CARTESIAN::height_at_fraction(
	double		h0,
	double		dh,
	double		frac,
	int		dir,
	COMPONENT	vol_comp)
{
	bool stop;
	double height = h0;
	double old_frac, new_frac;

	accumulate_fractions_in_layer(height, &old_frac, vol_comp);

	if(old_frac >= frac)
	{
	    height += dir*dh*frac/old_frac;
	    return height;
	}

	stop = NO;
	do
	{
	    accumulate_fractions_in_layer(height+dir*dh, &new_frac, vol_comp);
	    if (((old_frac <= frac) && (frac <= new_frac)) ||
		((new_frac <= frac) && (frac <= old_frac)))
	    {
		height += dir*dh*(frac-old_frac) / (new_frac-old_frac);
		stop = YES;
	    }
	    else
	    {
		old_frac = new_frac;
		height += dir*dh;
	    }
	}
	while (stop == NO &&
		height >= front->rect_grid->GL[dim-1] &&
		height <= front->rect_grid->GU[dim-1]);

	return height;
}

void G_CARTESIAN::accumulate_fractions_in_layer(
	double		height,
	double		*frac,
	COMPONENT	vol_comp)
{
	RECT_GRID *gr = computational_grid(front->interf);
	int zdir = dim - 1;
	int Nx, Ny;
	int i, j;
	double dh[MAXD-1];
	double coords[MAXD];
	double light_fluid_volume = 0.0;

	coords[zdir] = height;
	Nx = irint(gr->gmax[0]*iext->rfactor);
	dh[0] = (gr->U[0]-gr->L[0])/Nx;

	switch(dim)
	{
	    case 1:
		printf("1D not supported for accumulate_fractions_in_layer().\n");
		clean_up(ERROR);
		break;
	    case 2:
		for (i = 0; i < Nx; i++)
		{
		    coords[0] = gr->L[0] + i*dh[0];
		    light_fluid_volume +=
			find_particular_fluid_cell_volume(coords,dh,vol_comp);
		}
		if (pp_numnodes() > 1)
		    pp_global_sum(&light_fluid_volume, 1);
		*frac = light_fluid_volume / (gr->GU[0]-gr->GL[0]);
		break;
	    case 3:
		Ny = irint(gr->gmax[1]*iext->rfactor);
		dh[1] = (gr->U[1]-gr->L[1])/Nx;
		for (j = 0; j < Ny; j++)
		{
		    coords[1] = gr->L[1] + j*dh[1];
		    for (i = 0; i < Nx; i++)
		    {
			coords[0] = gr->L[0] + i*dh[0];
			light_fluid_volume +=
			    find_particular_fluid_cell_volume(coords,dh,vol_comp);
		    }
		}
		if (height < gr->L[zdir] || height >= gr->U[zdir])
		    light_fluid_volume = 0.0;

		if (pp_numnodes() > 1)
		{
		    pp_global_sum(&light_fluid_volume, 1);
		    pp_gsync();
		}

		*frac = light_fluid_volume / ((gr->GU[0]-gr->GL[0])*(gr->GU[1]-gr->GL[1]));

		break;
	}
}

double G_CARTESIAN::find_particular_fluid_cell_volume(
	double		*crds,
	double		*dh,
	COMPONENT	vol_comp)
{
	/*This part is copied from gas*/
        INTERFACE	*intfc = front->interf;
        double		ans = 0.0;
	double		previous_coords[MAXD], coords[MAXD];
	double		*intersection_coord;
	/* The following, intersection_length, is by default computed
	   as the length of the line segment from the LOWER x or y bound 
	   (respectively) of a grid cell line to where the interface cuts 
	   that x or y (resp.) grid cell line. The values are 
	   adjusted appropriately in case it is the length on 
	   the other side of the intersection which is needed.
	*/
	double		*intersection_length;
	COMPONENT	*comp;
	int		vol_comp_count, num_points;
	int		i;


	uni_array(&intersection_coord, (int)pow(2,MAXD-1), sizeof(double));
	uni_array(&intersection_length, (int)pow(2,MAXD-1), sizeof(double));
	uni_array(&comp, (int)pow(2,MAXD-1), sizeof(int));

	for (i = 0; i < MAXD; i++)
	    previous_coords[i] = coords[i] = crds[i];
	/* Note in 2d, we really only need intersection_coord and 
	   intersection_length to be scalars.
	*/
	for (i = 0; i < (int)pow(2,MAXD-1); i++)
	    intersection_coord[i] = intersection_length[i] = -HUGE;


	/* If the dimension of the interface is not 2 or 3, return 
	   the initial value of ans  (0.0).
	*/
	if (front->rect_grid->dim == 2)
	{
	    coords[0] += dh[0]; 
	    comp[0] = component(previous_coords, intfc);
	    comp[1] = component(coords, intfc);
		  
	    if (comp[0] != comp[1])
	    {
		intersection_coord[0] = find_intersection_by_bisection(
					     front,previous_coords,coords,0);
		intersection_length[0] = intersection_coord[0] - 
		                         previous_coords[0];
	    }

	    if (comp[0] == vol_comp && comp[1] == vol_comp)
	        ans = dh[0];
	    else if (comp[0] == vol_comp && comp[1] != vol_comp)
	        ans = intersection_length[0];
	    else if (comp[0] != vol_comp && comp[1] == vol_comp)
	        ans = dh[0] - intersection_length[0];
	    /* else, return 0.0 */
	}



	else if (front->rect_grid->dim == 3)
	{
	    comp[0] = component(coords,intfc);


	    coords[0] += dh[0];
	    comp[1] = component(coords,intfc);
	    if (comp[1] != comp[0])
	    {
		intersection_coord[0] = find_intersection_by_bisection(
					 front,previous_coords, coords, 0);
		intersection_length[0] = intersection_coord[0] - 
		                         previous_coords[0];
	    }
	    previous_coords[0] = coords[0];


	    coords[1] += dh[1]; 
	    comp[2] = component(coords,intfc);
	    if (comp[2] != comp[1])
	    {
		intersection_coord[1] = find_intersection_by_bisection(
					  front,previous_coords, coords, 1);
		intersection_length[1] = intersection_coord[1] - 
		                         previous_coords[1];
	    }
	    previous_coords[1] = coords[1];
		    

	    coords[0] -= dh[0];
	    comp[3] = component(coords,intfc);
	    if (comp[3] != comp[2])
	    {
		intersection_coord[2] = find_intersection_by_bisection(
					  front,previous_coords, coords, 0);
		intersection_length[2] = intersection_coord[2] - coords[0];
	    }
	    previous_coords[0] = coords[0];


	    coords[1] -= dh[1];
	    if (comp[0] != comp[3])
	    {
		intersection_coord[3] = find_intersection_by_bisection(
					  front,previous_coords, coords, 1);
		intersection_length[3] = intersection_coord[3] - coords[1];
	    }




	    vol_comp_count = 0;
	    num_points = 4;
	    for (i = 0; i < num_points; i++)
	    {
		if (comp[i] == vol_comp)
		    vol_comp_count++;
	    }
	    
	    /* There are sixteen basic cases to consider, which 
	       I have amalgamated into 3 main cases.  
	    */
	    
	    /* CASE 1: No cell edges are crossed by the interface. 
	       This accounts for 2 of the basic cases.
	    */
	    if (vol_comp_count == 0 || vol_comp_count == 4)
	    {
		if (comp[0] == vol_comp)
		    ans = dh[0]*dh[1];
		else
		    ans = 0.0;
	    }
	    
	    
	    /* CASE 2: Two adjacent edges are crossed by the interface
	       and no other edges are crossed. This contains four 
	       sub-cases. In total, this main case takes care of 
	       all eight of the basic cases in which only one of the 
	       cell nodes differs in component value from the other three.
	    */
	    else if (vol_comp_count == 1 || vol_comp_count == 3)
	    {
		if ((comp[0] != comp[1]) && (comp[1] != comp[2]))
		{
		    ans = 0.5 * (dh[0] - intersection_length[0]) * 
		      (intersection_length[1]);
		    if (comp[1] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
		
		
		else if ((comp[1] != comp[2]) && (comp[2] != comp[3]))
		{
		    ans = 0.5 * (dh[1]-intersection_length[1]) * 
		      (dh[0]-intersection_length[2]);
		    if (comp[2] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
		
		
		else if ((comp[2] != comp[3]) && (comp[3] != comp[0]))
		{
		    ans = 0.5 * (intersection_length[2]) * 
		      (dh[1]-intersection_length[3]);
		    if (comp[3] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
		
		else
		{
		    ans = 0.5 * (intersection_length[3]) * 
		      (intersection_length[0]);
		    if (comp[0] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
	    }
	    


	    /* CASE 3: Two non-adjacent edges are crossed by the interface
	       and no other edges are crossed or all four edges are crossed 
	       by the interface. This contains the remaining six basic cases.
	    */
	    else if (vol_comp_count == 2)
	    {
		/* First the four cases in which two adjacent cell nodes 
		   have the same component value
		*/
		if (comp[0] == comp[1])
		{
		    ans = dh[0] * 0.5 * (intersection_length[1] + 
					 intersection_length[3]); 
		    if (comp[0] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
		
		else if (comp[1] == comp[2])
		{
		    ans = dh[1] * 0.5* (intersection_length[0] + 
					intersection_length[2]); 
		    if (comp[1] == vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }  
		}
		
		/* Finally, there are the two tricky cases in which 
		   the interface cuts each of the four cell edges. 
		   There are two possible solutions to each case.
		   I break ties here by determining the component 
		   in the cell centre.
		*/
		else 
		{
		    coords[0] = crds[0] + 0.5*dh[0];
		    coords[1] = crds[1] + 0.5*dh[1];
		    
		    if (comp[0] == vol_comp)
		    {
			if (component(coords,intfc) != vol_comp)
			{
			    ans = 0.5 * (intersection_length[0]) * 
			          intersection_length[3];
			    ans += 0.5 * (dh[0] - intersection_length[2]) *
			          (dh[1]-intersection_length[1]);
			}
			else
			{
			    ans = 0.5 * (dh[0] - intersection_length[0]) *
			          intersection_length[1];
			    ans += 0.5 * intersection_length[2] * 
			           (dh[1] - intersection_length[3]);
			    ans *= -1.0;
			    ans += dh[0]*dh[1];
			}
		    }
		    else /* if comp[1] == vol_comp */
		    {
			if (component(coords,intfc) != vol_comp)
			{
			    ans = 0.5 * (dh[0]-intersection_length[0]) * 
			          intersection_length[1];
			    ans += 0.5 * (intersection_length[2]) * 
			           (dh[1]-intersection_length[3]);
			}
			else
			{
			    ans = 0.5 * (dh[1]-intersection_length[1]) * 
			          (dh[0]-intersection_length[2]);
			    ans += 0.5 * (intersection_length[0] * 
					  intersection_length[3]);
			    ans *= -1.0;
			    ans += dh[0]*dh[1];
			}
		    }
		}
	    }
	  
	    
	    /* Finally, if none of the four cell nodes has the same 
	       component value as vol_comp. This may seem to be a sub-case 
	       of CASE 1 above, but it in fact also takes care of the 
	       situation in which at least one of the cell nodes has 
	       NULL component value - which seems to happen within about 
	       3 cell blocks of the top and bottom of the global domain 
	       at time = 0. 4 Apr. 2002 - egeorge.
	    */
	    if ((comp[0] != vol_comp) && (comp[1] != vol_comp))
	    {
		if ((comp[2] != vol_comp) && (comp[3] != vol_comp))
		{
		    ans = 0.0;
		}
	    }
	} /* end of else if (front->rect_grid->dim == 3) */


	free(intersection_coord);
	free(intersection_length);
	free(comp);
  
	return ans;
}

/*This function is copied from gas*/
double G_CARTESIAN::find_intersection_by_bisection(
	Front	*front,
	double	*coords_a,
	double	*coords_b,
	int	dir)
{
        double		ans;
	double		crds_a[MAXD], crds_b[MAXD], cross_crds[MAXD];
	INTERFACE	*intfc = front->interf;
	COMPONENT	comp, comp_a, comp_b;
	int		i;
	int		max_num_bisections = 10;



	for (i = 0; i < MAXD; i++)
	{
	    crds_a[i] = cross_crds[i] = coords_a[i];
	    crds_b[i] = coords_b[i];
	}

	/* Now arrange it so that the crds_a[dir] < crds_b[dir] */

	if (coords_a[dir] > coords_b[dir])
	{
	    for (i = 0; i < MAXD; i++)
	    {
		crds_a[i] = crds_b[i];
		crds_b[i] = cross_crds[i];
	    }
	}
	  

	/* First check if coords_a or coords_b happen to be ONFRONT */
	if (component(crds_a, intfc) == ONFRONT)
	{
	    ans = crds_a[dir];
	}
	else if (component(crds_b, intfc) == ONFRONT)
	{
	    ans = crds_b[dir];
	}
	else
	{
	    for (i = 0; i < max_num_bisections; i++)
	    {      	
		cross_crds[dir] = 0.5*(crds_a[dir] + crds_b[dir]);
		comp = component(cross_crds,intfc);

		if (comp == ONFRONT)
		    break;		
		comp_a = component(crds_a,intfc);
		comp_b = component(crds_b,intfc);
	
		if (comp == comp_a)
		    crds_a[dir] = cross_crds[dir];
		else
		    crds_b[dir] = cross_crds[dir];
	    }
	    ans = cross_crds[dir];
	}
   
	return ans;
}		/*end find_intersection_by_bisection*/ 
/*This function is copied from gas*/

void G_CARTESIAN::print_intfc_extrema(char *out_name)
{
	char filename[500];
	FILE *outfile;

	if(pp_mynode() == 0)
	{
	    /*for max_00*/
	    sprintf(filename, "%s/intfc_extrema.max00",out_name);

	    if(front->step == 0)
	    {
	        outfile = fopen(filename,"w");
	        (void) fprintf(outfile," %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n",
	    	    "#time","h_max","vel","amb_vel","dens","amb_dens","pres","amb_pres","engy","amb_engy");
	    }
	    else
	        outfile = fopen(filename,"a");

	    (void) fprintf(outfile,"%15.5e", front->time);
	    (void) fprintf(outfile,"%15.5e", iext->h_max);
	    (void) fprintf(outfile,"%15.5e", bst_00[1].v[dim-1]);
	    (void) fprintf(outfile,"%15.5e", amb_bst_00[1].v[dim-1]);
	    (void) fprintf(outfile,"%15.5e", bst_00[1].d);
	    (void) fprintf(outfile,"%15.5e", amb_bst_00[1].d);
	    (void) fprintf(outfile,"%15.5e", bst_00[1].p);
	    (void) fprintf(outfile,"%15.5e", amb_bst_00[1].p);
	    (void) fprintf(outfile,"%15.5e", bst_00[1].e);
	    (void) fprintf(outfile,"%15.5e", amb_bst_00[1].e);
	    (void) fprintf(outfile,"\n");

	    fclose(outfile);

	    /*for min_00*/
	    sprintf(filename, "%s/intfc_extrema.min00",out_name);

	    if(front->step == 0)
	    {
	        outfile = fopen(filename,"w");
	        (void) fprintf(outfile," %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n",
	    	    "#time","h_min","vel","amb_vel","dens","amb_dens","pres","amb_pres","engy","amb_engy");
	    }
	    else
	        outfile = fopen(filename,"a");

	    (void) fprintf(outfile,"%15.5e", front->time);
	    (void) fprintf(outfile,"%15.5e", iext->h_min);
	    (void) fprintf(outfile,"%15.5e", bst_00[0].v[dim-1]);
	    (void) fprintf(outfile,"%15.5e", amb_bst_00[0].v[dim-1]);
	    (void) fprintf(outfile,"%15.5e", bst_00[0].d);
	    (void) fprintf(outfile,"%15.5e", amb_bst_00[0].d);
	    (void) fprintf(outfile,"%15.5e", bst_00[0].p);
	    (void) fprintf(outfile,"%15.5e", amb_bst_00[0].p);
	    (void) fprintf(outfile,"%15.5e", bst_00[0].e);
	    (void) fprintf(outfile,"%15.5e", amb_bst_00[0].e);
	    (void) fprintf(outfile,"\n");

	    fclose(outfile);

	    if (iext->do_01 == YES)
	    {
		/*for max_01*/
		sprintf(filename, "%s/intfc_extrema.max01",out_name);

		if(front->step == 0)
		{
		    outfile = fopen(filename,"w");
		    (void) fprintf(outfile," %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n",
			    "#time","h_max","vel","amb_vel","dens","amb_dens","pres","amb_pres","engy","amb_engy");
		}
		else
		    outfile = fopen(filename,"a");

		(void) fprintf(outfile,"%15.5e", front->time);
		(void) fprintf(outfile,"%15.5e", iext->h_max_01);
		(void) fprintf(outfile,"%15.5e", bst_01[1].v[dim-1]);
		(void) fprintf(outfile,"%15.5e", amb_bst_01[1].v[dim-1]);
		(void) fprintf(outfile,"%15.5e", bst_01[1].d);
		(void) fprintf(outfile,"%15.5e", amb_bst_01[1].d);
		(void) fprintf(outfile,"%15.5e", bst_01[1].p);
		(void) fprintf(outfile,"%15.5e", amb_bst_01[1].p);
		(void) fprintf(outfile,"%15.5e", bst_01[1].e);
		(void) fprintf(outfile,"%15.5e", amb_bst_01[1].e);
		(void) fprintf(outfile,"\n");

		fclose(outfile);

		/*for min_01*/
		sprintf(filename, "%s/intfc_extrema.min01",out_name);

		if(front->step == 0)
		{
		    outfile = fopen(filename,"w");
		    (void) fprintf(outfile," %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n",
			    "#time","h_min","vel","amb_vel","dens","amb_dens","pres","amb_pres","engy","amb_engy");
		}
		else
		    outfile = fopen(filename,"a");

		(void) fprintf(outfile,"%15.5e", front->time);
		(void) fprintf(outfile,"%15.5e", iext->h_min_01);
		(void) fprintf(outfile,"%15.5e", bst_01[0].v[dim-1]);
		(void) fprintf(outfile,"%15.5e", amb_bst_01[0].v[dim-1]);
		(void) fprintf(outfile,"%15.5e", bst_01[0].d);
		(void) fprintf(outfile,"%15.5e", amb_bst_01[0].d);
		(void) fprintf(outfile,"%15.5e", bst_01[0].p);
		(void) fprintf(outfile,"%15.5e", amb_bst_01[0].p);
		(void) fprintf(outfile,"%15.5e", bst_01[0].e);
		(void) fprintf(outfile,"%15.5e", amb_bst_01[0].e);
		(void) fprintf(outfile,"\n");

		fclose(outfile);
	    }
	    if (iext->do_05 == YES)
	    {
		/*for max_05*/
		sprintf(filename, "%s/intfc_extrema.max05",out_name);

		if(front->step == 0)
		{
		    outfile = fopen(filename,"w");
		    (void) fprintf(outfile," %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n",
			    "#time","h_max","vel","amb_vel","dens","amb_dens","pres","amb_pres","engy","amb_engy");
		}
		else
		    outfile = fopen(filename,"a");

		(void) fprintf(outfile,"%15.5e", front->time);
		(void) fprintf(outfile,"%15.5e", iext->h_max_05);
		(void) fprintf(outfile,"%15.5e", bst_05[1].v[dim-1]);
		(void) fprintf(outfile,"%15.5e", amb_bst_05[1].v[dim-1]);
		(void) fprintf(outfile,"%15.5e", bst_05[1].d);
		(void) fprintf(outfile,"%15.5e", amb_bst_05[1].d);
		(void) fprintf(outfile,"%15.5e", bst_05[1].p);
		(void) fprintf(outfile,"%15.5e", amb_bst_05[1].p);
		(void) fprintf(outfile,"%15.5e", bst_05[1].e);
		(void) fprintf(outfile,"%15.5e", amb_bst_05[1].e);
		(void) fprintf(outfile,"\n");

		fclose(outfile);

		/*for min_05*/
		sprintf(filename, "%s/intfc_extrema.min05",out_name);

		if(front->step == 0)
		{
		    outfile = fopen(filename,"w");
		    (void) fprintf(outfile," %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n",
			    "#time","h_min","vel","amb_vel","dens","amb_dens","pres","amb_pres","engy","amb_engy");
		}
		else
		    outfile = fopen(filename,"a");

		(void) fprintf(outfile,"%15.5e", front->time);
		(void) fprintf(outfile,"%15.5e", iext->h_min_05);
		(void) fprintf(outfile,"%15.5e", bst_05[0].v[dim-1]);
		(void) fprintf(outfile,"%15.5e", amb_bst_05[0].v[dim-1]);
		(void) fprintf(outfile,"%15.5e", bst_05[0].d);
		(void) fprintf(outfile,"%15.5e", amb_bst_05[0].d);
		(void) fprintf(outfile,"%15.5e", bst_05[0].p);
		(void) fprintf(outfile,"%15.5e", amb_bst_05[0].e);
		(void) fprintf(outfile,"%15.5e", amb_bst_05[0].e);
		(void) fprintf(outfile,"\n");
		(void) fprintf(outfile,"\n");

		fclose(outfile);
	    }
	}
}

void G_CARTESIAN::accumulate_state_in_layer(
	double			height,
	Big_State		*bst,
	bool			ignore_params)
{
	INTERFACE		*intfc = front->interf;
	const RECT_GRID		*rgrid = front->rect_grid;
	const int		dim = rgrid->dim;
	const int		zdir = dim - 1;
	const int		n_params = 2;

	int		i, j, p;
	int		Nx, Ny;
	int		n_points;
	int		index;
	double		dx, dy;
	double		x_left, y_left;
	double		coords[3];
	COMPONENT	comp;

	STATE		*st;
	Big_State	*tmpbst;

	FT_ScalarMemoryAlloc((POINTER*)&st,sizeof(STATE));
	uni_array(&tmpbst,n_params,sizeof(Big_State));

	for (i = 0; i < n_params; i++)
	    zero_scalar(&tmpbst[i],sizeof(Big_State));

	n_points = 0;

	coords[zdir] = height;

	Nx = irint(rgrid->gmax[0]*iext->rfactor);
	dx = (rgrid->U[0] - rgrid->L[0])/Nx;
	x_left = rgrid->L[0] + 0.5*dx;

	switch(dim)
	{
	    case 1:
	        screen("ERROR in accumulate_state_in_layer(), "
		       "1D not supported\n");
	        clean_up(ERROR);
	        break;
	    case 2:
		if (height < rgrid->L[zdir] || height >= rgrid->U[zdir])
		    break;
		for (i = 0; i < Nx; i++)
		{
		    n_points++;
		    coords[0] = x_left + i*dx;
		    comp = component(coords,intfc);
		    hyp_sol(coords,st,comp);
		    index = comp - 2;
		    if (index != 0 && index != 1)
		    {
			printf("Not supported by accumulate_state_in_layer().\n");
			clean_up(ERROR);
		    }
		    add_state_to_totals(st,&tmpbst[index]);
		}
		break;
	    case 3:
		if (height < rgrid->L[zdir] || height >= rgrid->U[zdir])
		    break;
		Ny = irint(rgrid->gmax[1]*iext->rfactor);
		dy = (rgrid->U[1] - rgrid->L[1])/Ny;
		y_left = rgrid->L[1] + 0.5*dy;

		for (j = 0; j < Ny; j++)
		{
		    coords[1] = y_left + j*dy;
		    for (i = 0; i < Nx; i++)
		    {
		        n_points++;
			coords[0] = x_left + i*dx;
			comp = component(coords,intfc);
			hyp_sol(coords,st,comp);
			index = comp - 2;
			if (index != 0 && index != 1)
			{
			    printf("Not supported by accumulate_state_in_layer().\n");
			    clean_up(ERROR);
			}
			add_state_to_totals(st,&tmpbst[index]);
		    }
		}
		break;
	}

	const int	nn = pp_numnodes();
	if (nn > 1)
	{
	    for (p = 0; p < n_params; p++)
	    {
		pp_global_isum(&(tmpbst[p].count),1L);
		pp_global_sum(&(tmpbst[p].d),1L);
		for (i = 0; i < dim; i++)
		    pp_global_sum(&(tmpbst[p].v[i]),1L);
		pp_global_sum(&(tmpbst[p].p),1L);
		pp_global_sum(&(tmpbst[p].e),1L);
	    }
	}


	if (ignore_params == YES)	//for 0%
	{
	    zero_scalar(&bst[0],sizeof(Big_State));
	    for (p = 0; p < n_params; p++)
		accumulate_state_totals(&tmpbst[p],bst);
	    normalize_state_totals(bst);
	}
	else				//for 1% and 5%
	{
	    zero_scalar(&bst[0],sizeof(Big_State));
	    zero_scalar(&bst[1],sizeof(Big_State));
	    for (p = 0; p < n_params; p++)
	    {
		accumulate_state_totals(&tmpbst[p],&bst[p]);
		normalize_state_totals(&bst[p]);
	    }
	}

	free(st);
	free(tmpbst);
}

void G_CARTESIAN::hyp_sol(double *coords, STATE *st, COMPONENT comp)
{
	double		default_var;
	double		**m_mom = eqn_params->mom;
	double		*m_dens = eqn_params->dens;
	double		*m_engy = eqn_params->engy;
	double		*m_pres = eqn_params->pres;
	int		i;

	FT_NearestRectGridVarInRange(front,comp,coords,m_dens,2,&default_var);
	FT_IntrpStateVarAtCoords(front,comp,coords,m_dens,getStateDens,&st->dens,&default_var);
	FT_NearestRectGridVarInRange(front,comp,coords,m_engy,2,&default_var);
	FT_IntrpStateVarAtCoords(front,comp,coords,m_engy,getStateEngy,&st->engy,&default_var);
	FT_NearestRectGridVarInRange(front,comp,coords,m_pres,2,&default_var);
	FT_IntrpStateVarAtCoords(front,comp,coords,m_pres,getStatePres,&st->pres,&default_var);
	for (i = 0; i < dim; i++)
	{
	    FT_NearestRectGridVarInRange(front,comp,coords,m_mom[i],2,&default_var);
	    FT_IntrpStateVarAtCoords(front,comp,coords,m_mom[i],getStateMom[i],&st->momn[i],&default_var);
	    st->vel[i] = st->momn[i]/st->dens;
	}
	st->dim = dim;
}

void G_CARTESIAN::add_state_to_totals(STATE* st, Big_State *bst)
{
    	int i;

	bst->count++;
	bst->d += st->dens;
	for (i = 0; i < dim; i++)
	    bst->v[i] = st->vel[i];
	bst->p += st->pres;
	bst->e += st->engy;
}

void G_CARTESIAN::copy_Big_State(const Big_State *bst0, Big_State *bst)
{
    	int i;

	bst->count = bst0->count;
	bst->d = bst0->d;
	for (i = 0; i < dim; i++)
	    bst->v[i] = bst0->v[i];
	bst->p = bst0->p;
	bst->e = bst0->e;
}

void G_CARTESIAN::accumulate_state_totals(const Big_State *tmpbst, Big_State *bst)
{
    	int i;

    	bst->count += tmpbst->count;
	bst->d += tmpbst->d;
	for (i = 0; i < dim; i++)
	    bst->v[i] += tmpbst->v[i];
	bst->p += tmpbst->p;
	bst->e += tmpbst->e;
}

void G_CARTESIAN::normalize_state_totals(Big_State *bst)
{
    	int 	i;
	double	nf = (double)bst->count;

	bst->d /= nf;
	for (i = 0; i < dim; i++)
	    bst->v[i] /= nf;
	bst->p /= nf;
	bst->e /= nf;
}

//Dan	FIXME
void G_CARTESIAN::state_reflect(
	int	dir,
	double	*solute)
{
    	int		side, i;
    	INTERFACE	*intfc = front->grid_intfc;
	RECT_GRID	*comp_grid, *top_grid;
	int		lbuf[MAXD],ubuf[MAXD],*gmax;

	comp_grid = computational_grid(intfc);
	top_grid = &topological_grid(intfc);
	gmax = top_grid->gmax;
	for (i = 0; i < dim; ++i)
	{
	    lbuf[i] = comp_grid->lbuf[i];
	    ubuf[i] = comp_grid->ubuf[i];
	}

	for (side = 0; side < 2; side++)
	    if (rect_boundary_type(intfc,dir,side) == REFLECTION_BOUNDARY)
		reflect_buffer(dim,dir,side,gmax,lbuf,ubuf,solute);

	return;
}

//Dan	FIXME
void G_CARTESIAN::reflect_buffer(
        int	dim,
        int	dir,
        int	side,
        int	*gmax,
        int	*lbuf,
        int	*ubuf,
        double	*solute)
{
        int i, j, k;
        int index;
        switch (dim)
        {
        case 1:
            if (side == 0)
            {
                for (i = 0; i < lbuf[0]; ++i)
                {
                    index = d_index1d(lbuf[0]-1-i,gmax);
                    solute[index] = -solute[index];
                }
            }
            else
            {
                for (i = 0; i < ubuf[0]; ++i)
                {
                    index = d_index1d(gmax[0]-ubuf[0]+1+i,gmax);
                    solute[index] = -solute[index];
                }
            }
            break;
        case 2:
            if (side == 0)
            {
		switch (dir)
		{
		case 0:
		    for (j = 0; j <= gmax[1]; ++j)
                    for (i = 0; i < lbuf[0]; ++i)
		    {
                    	index = d_index2d(lbuf[0]-1-i,j,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		case 1:
		    for (i = 0; i <= gmax[0]; ++i)
                    for (j = 0; j < lbuf[1]; ++j)
		    {
                    	index = d_index2d(i,lbuf[1]-1-j,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		}
	    }
            else
            {
		switch (dir)
		{
		case 0:
		    for (j = 0; j <= gmax[1]; ++j)
                    for (i = 0; i < ubuf[0]; ++i)
		    {
                    	index = d_index2d(gmax[0]-ubuf[0]+1+i,j,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		case 1:
		    for (i = 0; i <= gmax[0]; ++i)
                    for (j = 0; j < ubuf[1]; ++j)
		    {
                    	index = d_index2d(i,gmax[1]-ubuf[1]+1+j,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		}
            }
            break;
        case 3:
            if (side == 0)
            {
		switch (dir)
		{
		case 0:
		    for (k = 0; k <= gmax[2]; ++k)
		    for (j = 0; j <= gmax[1]; ++j)
                    for (i = 0; i < lbuf[0]; ++i)
		    {
                    	index = d_index3d(lbuf[0]-1-i,j,k,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		case 1:
		    for (k = 0; k <= gmax[2]; ++k)
		    for (i = 0; i <= gmax[0]; ++i)
                    for (j = 0; j < lbuf[1]; ++j)
		    {
                    	index = d_index3d(i,lbuf[1]-1-j,k,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		case 2:
		    for (i = 0; i <= gmax[0]; ++i)
		    for (j = 0; j <= gmax[1]; ++j)
                    for (k = 0; k < lbuf[2]; ++k)
		    {
                    	index = d_index3d(i,j,lbuf[2]-1-k,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		}
	    }
            else
            {
		switch (dir)
		{
		case 0:
		    for (k = 0; k <= gmax[2]; ++k)
		    for (j = 0; j <= gmax[1]; ++j)
                    for (i = 0; i < ubuf[0]; ++i)
		    {
                    	index = d_index3d(gmax[0]-ubuf[0]+1+i,j,k,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		case 1:
		    for (k = 0; k <= gmax[2]; ++k)
		    for (i = 0; i <= gmax[0]; ++i)
                    for (j = 0; j < ubuf[1]; ++j)
		    {
                    	index = d_index3d(i,gmax[1]-ubuf[1]+1+j,k,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		case 2:
		    for (i = 0; i <= gmax[0]; ++i)
		    for (j = 0; j <= gmax[1]; ++j)
                    for (k = 0; k < ubuf[2]; ++k)
		    {
                    	index = d_index3d(i,j,gmax[2]-ubuf[2]+1+k,gmax);
                    	solute[index] = -solute[index];
		    }
		    break;
		}
            }
            break;
        }
}
