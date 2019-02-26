#include <cFluid.h>

void getPolyhCent(CPOLYHEDRON*,double*);
void getRMState(STATE*,EQN_PARAMS*,double*,COMPONENT);
void getTESTState(STATE*,EQN_PARAMS*,double*,COMPONENT);
void behind_state(int,double,double*,int,STATE*,STATE*);
void set_dir_and_side(int*,int*,CPOLYGON*,CELL*);
bool need_exterior_pam(CPAM*,int*,int*,int*);
void print_cpoint(CPOINT*);
void print_cpolygon(CPOLYGON*);
void print_cpolyh(CPOLYHEDRON*);
bool corr_polyh_in_cell(CPOLYHEDRON*,CELL*,int*);

//currently for TWO_FLUID_RM only
void G_CARTESIAN::cft_set_init_polyh_states()
{
	int i, j, k, index;
	double coords[3];
	STATE state;
	CPOLYHEDRON *polyh;
	bool update;

	if (eqn_params->prob_type != TWO_FLUID_RM)
	{
	    printf("ERROR: the initialization for CFT has not been implemented yet"
		   "for problems other than TWO_FLUID_RM.\n");
	    clean_up(ERROR);
	}

	cells = cells_old;	//FIXME

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    update = false;
	    index = d_index3d(i,j,k,top_gmax);
	    polyh = cells[index].polyhs;
	    //if (polyh->cut == TRUE)	TODO
	    while (polyh)
	    {
		getPolyhCent(polyh,coords);
		//set component
		//this part should be removed	TODO
		/*
		if (polyh->scp_dir == INW)
		{
		    polyh->comp = GAS_COMP2;
		}
		else if (polyh->scp_dir == OUTW)
		{
		    polyh->comp = GAS_COMP1;
		}
		else
		{
		    printf("ERROR in cft_set_init_polyh_states(): unset scp_dir.\n");
		    clean_up(ERROR);
		}
		*/
		//set state
		getRMState(&(polyh->state),eqn_params,coords,polyh->comp);
		//set test state	FIXME
		//getTESTState(&(polyh->state),eqn_params,coords,polyh->comp);
		update = true;
		polyh = polyh->next;
	    }

	    //update cell states also
	    if (update)
	    {
		cft_update_cell_states(&cells[index]);
	    }
	}

	return;
}

void G_CARTESIAN::cft_update_cell_states(CELL *c)
{
	int i, index;
	double *dens = field.dens;
	double **pdens = field.pdens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
	double mass, cvol;
	double vel[3];
	CPOLYHEDRON *polyh;
	STATE state;
	COMPONENT comp;

	//states need to be kept (vel, pres)
	index = d_index3d(c->icrds[0],c->icrds[1],c->icrds[2],top_gmax);
	cvol = top_h[0] * top_h[1] * top_h[2];
	for (int i = 0; i < 3; i++)
	    vel[i] = momn[i][index]/(dens[index]*cvol);

	//dens
	mass = 0;
	polyh = c->polyhs;
	while (polyh)
	{
	    mass += polyh->state.dens * polyh->vol;
	    polyh = polyh->next;
	}
	dens[index] = mass/cvol;

	//momn
	for (i = 0; i < 3; i++)
	    momn[i][index] = vel[i] * mass;

	//engy
	comp = top_comp[index];
	state.dim = eqn_params->dim;
	state.pres = pres[index];
	state.dens = dens[index];
	for (i = 0; i < 3; i++)
	    state.momn[i] = momn[i][index];
	state.eos = &(eqn_params->eos[comp]);
	engy[index] = EosEnergy(&state);

	//pdens?	TODO

	return;
}

//TODO This coords is not as good as that of the centroid.
//TODO However, results will not be changed a lot.
//TODO To be improved.
void getPolyhCent(
	CPOLYHEDRON	*polyh,
	double		*coords)
{
	int i, count;
	double fcrds[3];
	CPOINT *v;
	CPOLYGON *face;

	if (polyh->iscell)
	{
	    for (i = 0; i < 3; i++)
		coords[i] = (polyh->cell->celll[i] + polyh->cell->cellu[i]) / 2;
	    return;
	}

	count = 0;
	for (i = 0; i < 3; i++)
	    coords[i] = 0.0;
	face = polyh->faces;
	while (face)
	{
	    for (i = 0; i < 3; i++)
		fcrds[i] = 0.0;
	    v = face->vertices;
	    while (v)
	    {
		count++;
		for (i = 0; i < 3; i++)
		    fcrds[i] += v->crds[i];
		v = v->next;
	    }
	    for (i = 0; i < 3; i++)
		coords[i] += fcrds[i];
	    face = face->next;
	}

	for (i = 0; i < 3; i++)
	    coords[i] = coords[i]/count;

	return;
}

void getTESTState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	FOURIER_POLY *wave_params;
	EOS_PARAMS	*eos;
	double rho1 = eqn_params->rho1;
	double rho2 = eqn_params->rho2;
	double p0 = eqn_params->p0;
	double shock_position = eqn_params->shock_position;
	double Mach_number = eqn_params->Mach_number;
	double shock_speed;
	double csp = eqn_params->contact_vel;
	int shock_dir = eqn_params->shock_dir;
	int i,dim;
 
	wave_params = (FOURIER_POLY*)eqn_params->level_func_params;
	dim = wave_params->dim;

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->vel[i] = state->momn[i] = 0.0;
	state->dim = dim;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;

	
	switch (comp)
	{
	case GAS_COMP1:
	    state->dens = rho1;
	    state->pres = p0;
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                {
                    state->pdens[0] = state->dens;
                    state->pdens[1] = 0.0;
                }
            }
	    state->engy = EosInternalEnergy(state);
	    break;
	case GAS_COMP2:
	    state->dens = rho2;
	    state->pres = p0;
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                {
                    state->pdens[0] = 0.0;
                    state->pdens[1] = state->dens;
                }
            }
	    state->engy = EosInternalEnergy(state);
	    break;
	case EXT_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                state->pdens[0] = 0.0;
                state->pdens[1] = 0.0;
            }
	    state->engy = 0.0;
	    return;	//Dan	FIXME
//	    break;
	default:
	    printf("ERROR: Unknown component %d in getRMState()!\n",comp);
	    clean_up(ERROR);
	}

	//initial diffusion layer
	/*
	double idl = eqn_params->thickness_idl;
	double z0 = wave_params->z0;
	double z_intfc_pert = level_wave_func(wave_params,coords);
	if (eqn_params->multi_comp_non_reactive == YES)
	{
	    state->pdens[0] = 0.5*rho1 + 0.5*(0.0-rho1)*erf((coords[dim-1]-z_intfc_pert)/idl*4);
	    state->pdens[1] = 0.5*rho2 + 0.5*(rho2-0.0)*erf((coords[dim-1]-z_intfc_pert)/idl*4);
	}
	*/
	//state->dens = 0.5*(rho1+rho2) + 0.5*(rho1-rho2)*erf((coords[dim-1]-z_intfc_pert)/idl*4);

	if (debugging("rm_state"))
	{
	    printf("Before calling behind_state()\n");
	    printf("state = %f %f %f\n",state->dens,state->pres,
					state->vel[0]);
	}
	/*
	if ((shock_dir ==  1 && coords[dim-1] < shock_position) ||
	    (shock_dir == -1 && coords[dim-1] > shock_position))
	{
	    behind_state(SHOCK_MACH_NUMBER,Mach_number,
			&shock_speed,shock_dir,state,state);	
	    state->engy = EosEnergy(state);
	    if (debugging("rm_state"))
	    {
	    	printf("After calling behind_state()\n");
	    	printf("state = %f %f %f\n",state->dens,state->pres,
			state->vel[0]);
	    }
	}
	*/
	//state->vel[dim-1] -= csp;
	state->vel[dim-1] = csp;
	state->momn[dim-1] = state->vel[dim-1]*state->dens;
	state->engy = EosEnergy(state);

	return;
}

void getRMState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	FOURIER_POLY *wave_params;
	EOS_PARAMS	*eos;
	double rho1 = eqn_params->rho1;
	double rho2 = eqn_params->rho2;
	double p0 = eqn_params->p0;
	double shock_position = eqn_params->shock_position;
	double Mach_number = eqn_params->Mach_number;
	double shock_speed;
	double csp = eqn_params->contact_vel;
	int shock_dir = eqn_params->shock_dir;
	int i,dim;
 
	wave_params = (FOURIER_POLY*)eqn_params->level_func_params;
	dim = wave_params->dim;

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->vel[i] = state->momn[i] = 0.0;
	state->dim = dim;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;

	
	switch (comp)
	{
	case GAS_COMP1:
	    state->dens = rho1;
	    state->pres = p0;
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                {
                    state->pdens[0] = state->dens;
                    state->pdens[1] = 0.0;
                }
            }
	    state->engy = EosInternalEnergy(state);
	    break;
	case GAS_COMP2:
	    state->dens = rho2;
	    state->pres = p0;
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                {
                    state->pdens[0] = 0.0;
                    state->pdens[1] = state->dens;
                }
            }
	    state->engy = EosInternalEnergy(state);
	    break;
	case EXT_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
            if(eqn_params->multi_comp_non_reactive == YES)
            {
                state->pdens[0] = 0.0;
                state->pdens[1] = 0.0;
            }
	    state->engy = 0.0;
	    return;	//Dan	FIXME
//	    break;
	default:
	    printf("ERROR: Unknown component %d in getRMState()!\n",comp);
	    clean_up(ERROR);
	}

	//initial diffusion layer
	/*
	double idl = eqn_params->thickness_idl;
	double z0 = wave_params->z0;
	double z_intfc_pert = level_wave_func(wave_params,coords);
	if (eqn_params->multi_comp_non_reactive == YES)
	{
	    state->pdens[0] = 0.5*rho1 + 0.5*(0.0-rho1)*erf((coords[dim-1]-z_intfc_pert)/idl*4);
	    state->pdens[1] = 0.5*rho2 + 0.5*(rho2-0.0)*erf((coords[dim-1]-z_intfc_pert)/idl*4);
	}
	*/
	//state->dens = 0.5*(rho1+rho2) + 0.5*(rho1-rho2)*erf((coords[dim-1]-z_intfc_pert)/idl*4);

	if (debugging("rm_state"))
	{
	    printf("Before calling behind_state()\n");
	    printf("state = %f %f %f\n",state->dens,state->pres,
					state->vel[0]);
	}
	if ((shock_dir ==  1 && coords[dim-1] < shock_position) ||
	    (shock_dir == -1 && coords[dim-1] > shock_position))
	{
	    behind_state(SHOCK_MACH_NUMBER,Mach_number,
			&shock_speed,shock_dir,state,state);	
	    state->engy = EosEnergy(state);
	    if (debugging("rm_state"))
	    {
	    	printf("After calling behind_state()\n");
	    	printf("state = %f %f %f\n",state->dens,state->pres,
			state->vel[0]);
	    }
	}
	state->vel[dim-1] -= csp;
	state->momn[dim-1] = state->vel[dim-1]*state->dens;
	state->engy = EosEnergy(state);

	return;
}	/* end getRMState */

void behind_state(
	int		which_parameter,
	double		parameter,
	double		*shock_speed,
	int		shock_dir,
	STATE		*ahead_state,
	STATE		*behind_state)
{
	double		r0, p0, u0;		/* ahead state */
	double		r1, p1, u1;		/* behind state */
	double		U;			/* shock speed */
	double		M0n;			/* shock mack number,
						   relative to ahead flow */
	double		M0nsq;			/* steady normal ahead Mach
						   number squared */
	int		dim;

	dim = ahead_state->dim;
	r0  = ahead_state->dens;
	p0  = ahead_state->pres;
	u0  = ahead_state->vel[dim-1]*shock_dir;

	switch(which_parameter)
	{
	case SHOCK_MACH_NUMBER:
	    M0n = parameter;
	    *shock_speed = U = u0 + M0n*EosSoundSpeed(ahead_state);
	    M0nsq = sqr(M0n);
	    p1 = EosMaxBehindShockPres(M0nsq,ahead_state);
	    u1 =  u0 + (p0 - p1) / (r0*(u0 - U)); 
	    r1 = r0*((u0 - U)/(u1 - U));
	    if (debugging("rm_state"))
	    {
		printf("M0n = %f  shock_speed = %f\n",M0n,*shock_speed);
		printf("p1 = %f  u1 = %f  r1 = %f\n",p1,u1,r1);
	    }
	    break;
	default:
	    screen("ERROR in behind_state(), "
	           "unknown parameter %d\n",which_parameter);
	    clean_up(ERROR);
	}	
	behind_state->dens = r1;
	//Dan	FIXME
	if (behind_state->pdens[0] < 1e-12)
	    behind_state->pdens[1] = r1;
	else
	    behind_state->pdens[0] = r1;
	//Dan	FIXME
	behind_state->pres = p1;
	behind_state->vel[dim-1] = u1*shock_dir;
	behind_state->momn[dim-1] = r1*u1*shock_dir;
}		/*end behind_state */

void G_CARTESIAN::cft_set_face_flux()
{
	int i, j, k, index, ii, jj, kk, indexx, ic;
	int dir, side, sign;
	double cf_area, f, deltax;
	double min_crds[3], max_crds[3], cflux_gmax[3];
	CELL *c, *cells;
	CPAM *pam;
	CPOLYHEDRON *polyh;
	CPOLYGON *face;
	COMPONENT comp;
	bool debugcft = false;

	cells = cells_halft;
	for (i = 0; i < 3; i++)
	    cflux_gmax[i] = top_gmax[i]+1;
	cf_area = top_h[0]*top_h[1];	//requires top_h[i]s are equivalent to each other	FIXME
	deltax = top_h[0];

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);

	    //debugdan	FIXME
	    /*
	    if (i == 4 && j == 4)
		debugcft = true;
	    else
		debugcft = false;
	    */

	    if (c->merged == FALSE)
		continue;

	    //check if this comp is the same as comp[index] or not
	    /*
	    comp = c->pams->polyh->comp;
	    if (comp == GAS_COMP1)
		ic = 0;
	    else if (comp == GAS_COMP2)
		ic = 1;
	    else
	    {
		printf("ERROR in cft_set_face_flux(): incorrect comp.\n");
		clean_up(ERROR);
	    }
	    */

	    //pam = c->pams;
	    polyh = c->polyhs;
	    //while (pam)
	    while (polyh)
	    {
		//polyh = pam->polyh;
		polyh->mflux = 0;
		comp = polyh->comp;
		if (comp == GAS_COMP1)
		    ic = 0;
		else if (comp == GAS_COMP2)
		    ic = 1;
		else
		{
		    printf("ERROR in cft_set_face_flux(): incorrect comp.\n");
		    clean_up(ERROR);
		}

		//debugdan	FIXME
		if (debugcft)
		    printf("Polyh:\n");
		if (polyh->iscell)
		{
		    //do we need to reset flux?	TODO
		    //set flux for polyh
		    //polyh->mflux = 0;
		    for (dir = 0; dir < 3; dir++)
		    {
			for (side = 0; side < 2; side++)
			{
			    ii = i;
			    jj = j;
			    kk = k;
			    if (dir == 0)
				ii = ii + side;
			    else if (dir == 1)
				jj = jj + side;
			    else
				kk = kk + side;
			    if (side == 0)
				sign = -1;
			    else
				sign = 1;
			    indexx = d_index3d(ii,jj,kk,cflux_gmax);
			    polyh->mflux += sign*cflux[ic][dir].dens_flux[indexx];
			    //debugdan	FIXME
			    if (debugcft)
				printf("%d %d %d: dir = %d, side = %d, "
					"dens_flux[%d] = %e, mflux = %e.\n",
					i, j, k, dir, side, indexx, 
					cflux[ic][dir].dens_flux[indexx], 
					polyh->mflux);
			}
		    }
		    polyh->mflux *= polyh->vol;
		}
		else
		{
		    if (debugcft)
		    {
			printf("polyh is not a cell.\n");
		    }
		    face = polyh->faces;
		    while (face)
		    {
			if (face->oncf)
			{
			    //decide direction and side
			    set_dir_and_side(&dir,&side,face,c);

			    //set index for cflux
			    ii = i;
			    jj = j;
			    kk = k;
			    if (dir == 0)
				ii = (side==0) ? ii : ii+1;
			    else if (dir == 1)
				jj = (side==0) ? jj : jj+1;
			    else if (dir == 2)
				kk = (side==0) ? kk : kk+1;
			    else
			    {
				printf("ERROR in cft_set_face_flux(): incorrect dir.\n");
				clean_up(ERROR);
			    }
			    if (side == 0)
				sign = -1;
			    else
				sign = 1;
			    indexx = d_index3d(ii,jj,kk,cflux_gmax);

			    //set dens flux
			    //face->dens_flux = sign*face->area/cf_area*cflux[ic][dir].dens_flux[indexx];
			    //set mass flux
			    //face->mass_flux = face->dens_flux*polyh->vol;
			    face->mass_flux = sign*deltax*face->area*cflux[ic][dir].dens_flux[indexx];
			    polyh->mflux += face->mass_flux;

			    //debugdan	FIXME
			    if (debugcft)
			    {
				printf("%d %d %d: comp = %d, cflux[%d][%d] = %e, mflux = %e.\n",
					i, j, k, polyh->comp, ic, dir, 
					cflux[ic][dir].dens_flux[indexx],
					polyh->mflux);
			    }
			    /*
			    if (i == 4 && j == 4 && k == 21 && polyh->comp ==2)
			    {
				printf("4 4 21 vol = %e, dir = %d, side = %d, dens_flux[%d] = %e (%e).\n",
					polyh->vol, dir, side, indexx, face->dens_flux,
					cflux[ic][dir].dens_flux[indexx]);
			    }
			    */
			    
			}
			face = face->next;
		    }
		}

		//pam = pam->next;
		polyh = polyh->next;
	    }
	}

	return;
}

void set_dir_and_side(
	int		*dir,
	int		*side,
	CPOLYGON	*face,
	CELL		*c)
{
	int i;
	double min_crds[3], max_crds[3];
	CPOINT *p;
	double tol = 1e-12;

	p = face->vertices;
	for (i = 0; i < 3; i++)
	{
	    min_crds[i] = p->crds[i];
	    max_crds[i] = p->crds[i];
	}
	p = p->next;
	while (p)
	{
	    for (i = 0; i < 3; i++)
	    {
		if (p->crds[i] < min_crds[i])
		    min_crds[i] = p->crds[i];
		if (p->crds[i] > max_crds[i])
		    max_crds[i] = p->crds[i];
	    }
	    p = p->next;
	}

	for (i = 0; i < 3; i++)
	{
	    if (max_crds[i]-min_crds[i] < tol)
	    {
		*dir = i;
		break;
	    }
	}

	if (fabs(min_crds[i]-c->celll[i]) < tol)
	    *side = 0;
	else if (fabs(min_crds[i]-c->cellu[i]) < tol)
	    *side = 1;
	else
	{
	    printf("ERROR in set_dir_and_side().\n");
	    clean_up(ERROR);
	}

	return;
}

bool debugdan = false;

void G_CARTESIAN::cft_update_states_new()
{
	int i, j, k, index, ii, jj, kk, indexx, count, dd;
	int iii, jjj, kkk, indexxx;
	CELL *c;
	CFTCELL *cftc;
	CPAM *pam;
	CPOLYHEDRON *polyh;

	//set CFTCELLs
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells_new[index]);

	    //debugdan	FIXME
	    /*
	    if (i == 7 && j == 4 && k == 19)
		debugdan = true;
	    else
		debugdan = false;
	    if (debugdan)
	    {
		printf("%d %d %d index = %d, comp = %d.\n",
			i, j, k, index, c->pams->polyh->comp);
		
		polyh = c->polyhs;
		while (polyh)
		{
		    printf("comp = %d:\n", polyh->comp);
		    //print_cpolyh(polyh);
		    polyh = polyh->next;
		}
		
	    }
	    */
	    //debugdan	FIXME

	    cft_init_cftcell(c);
	    c->cftcell->pams = c->pams;

	    //set indices based on the direction of velocity
	    cft_set_indices_cftcell(c);

	    //set maps
	    cft_set_maps_for_cftcell(c->cftcell);

	    //set vol_new and vol_old
	    cft_set_vols_for_cftcell(c->cftcell);

	    //set mass at old time step
	    cft_set_mass_at_oldt_cftcell(c->cftcell);

	    //set flux at half time step
	    cft_set_flux_at_halft_cftcell(c->cftcell);
	}

	//update polyhs' states
	cft_update_polyhs_states();

	//update cells' states
	cft_update_cells_states();

	//printf("End of cft_update_states_new().\n");
	//exit(0);
}

void G_CARTESIAN::cft_update_cells_states()
{
	int i, j, k, index, ii;
	double *dens = field.dens;
	double **pdens = field.pdens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
	double vol_cell = top_h[0]*top_h[1]*top_h[2];
	STATE state;
	CELL *c;
	CPOLYHEDRON *polyh;

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells_new[index]);

	    if (c->merged == false)
		continue;

	    state.dens = 0.0;
	    for (int ii = 0; ii < 2; ii++)
		state.pdens[ii] = 0.0;
	    for (int ii = 0; ii < 3; ii++)
		state.momn[ii] = 0.0;
	    polyh = c->polyhs;
	    while (polyh)
	    {
		state.dens += polyh->state.dens * polyh->vol;
		for (int ii = 0; ii < 2; ii++)
		{
		    if (polyh->comp == GAS_COMP1)
			state.pdens[0] += polyh->state.dens * polyh->vol;
		    else
			state.pdens[1] += polyh->state.dens * polyh->vol;
		}
		for (int ii = 0; ii < 3; ii++)
		    state.momn[ii] += polyh->state.momn[ii] * polyh->vol;
		polyh = polyh->next;
	    }
	    state.dens /= vol_cell;
	    for (int ii = 0; ii < 2; ii++)
		state.pdens[ii] /= vol_cell;
	    for (int ii = 0; ii < 3; ii++)
		state.momn[ii] /= vol_cell;

	    dens[index] = state.dens;
	    for (int ii = 0; ii < 2; ii++)
		pdens[ii][index] = state.pdens[ii];
	    for (int ii = 0; ii < 3; ii++)
		momn[ii][index] = state.momn[ii];

	}

	return;
}

void G_CARTESIAN::cft_update_polyhs_states()
{
	int i, j, k, index, ii, jj, kk, indexx, d;
	double dens_new;
	double vel[3];
	double **momn = field.momn;
	double *dens = field.dens;
	double **pdens = field.pdens;
	double *engy = field.engy;
	double *pres = field.pres;
	double vol_cell = top_h[0]*top_h[1]*top_h[2];
	CELL *c;
	CFTCELL *cftc;
	CPAM *pam;
	CPOLYHEDRON *polyh;
	COMPONENT comp;

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells_new[index]);

	    if (c->merged == false)
	    {
		polyh = c->polyhs;
		polyh->state.dens = dens[index];
		polyh->state.engy = engy[index];
		polyh->state.pres = pres[index];
		for (int d = 0; d < 3; d++)
		{
		    polyh->state.momn[d] = momn[d][index];
		}
		if (comp == GAS_COMP1)
		{
		    polyh->state.pdens[0] = polyh->state.dens;
		    polyh->state.pdens[1] = 0.0;
		}
		else
		{
		    polyh->state.pdens[0] = 0.0;
		    polyh->state.pdens[1] = polyh->state.dens;
		}
	    }
	    else
	    {
		cftc = c->cftcell;
		cftc->mass_new = cftc->mass_old + cftc->mflux;
		dens_new = cftc->mass_new / cftc->vol_new;
		comp = cftc->pams->polyh->comp;
		pam = cftc->pams;
		while (pam)
		{
		    polyh = pam->polyh;

		    //cell index
		    ii = pam->polyh->cell->icrds[0];
		    jj = pam->polyh->cell->icrds[1];
		    kk = pam->polyh->cell->icrds[2];
		    indexx = d_index3d(ii,jj,kk,top_gmax);

		    //update dens
		    polyh->state.dens = dens_new;

		    //kepp vel and pres
		    for (int d = 0; d < 3; d++)
			vel[d] = momn[d][indexx] / dens[indexx];
		    polyh->state.pres = pres[indexx];

		    //update momn
		    for (int d = 0; d < 3; d++)
			polyh->state.momn[d] = vel[d] * polyh->state.dens;

		    //update engy
		    polyh->state.dim = eqn_params->dim;
		    if (comp == GAS_COMP1)
		    {
			polyh->state.pdens[0] = polyh->state.dens;
			polyh->state.pdens[1] = 0.0;
		    }
		    else
		    {
			polyh->state.pdens[0] = 0.0;
			polyh->state.pdens[1] = polyh->state.dens;
		    }
		    polyh->state.eos = &(eqn_params->eos[comp]);
		    polyh->state.engy = EosEnergy(&(polyh->state));

		    pam = pam->next;
		}
	    }
	}

	return;
}

void G_CARTESIAN::cft_set_indices_cftcell(CELL *c)
{
	int ii, jj, kk, indexx, iii, jjj, kkk, indexxx, dd;
	int vdir[3], c_icrds[3], c_imin[3], c_imax[3];
	CPAM *pam = c->pams;
	CPOLYHEDRON *polyh;
	double **momn = field.momn;

	while (pam)
	{
	    polyh = pam->polyh;
	    ii = polyh->cell->icrds[0];
	    jj = polyh->cell->icrds[1];
	    kk = polyh->cell->icrds[2];
	    indexx = d_index3d(ii,jj,kk,top_gmax);
	    c_icrds[0] = ii;
	    c_icrds[1] = jj;
	    c_icrds[2] = kk;
	    //debugdan	FIXME
	    if (debugdan)
	    {
		printf("merged %d %d %d:\n", ii, jj, kk);
		//print_cpolyh(polyh);
	    }
	    //debugdan	FIXME
	    for (dd = 0; dd < 3; dd++)
	    {
		if (momn[dd][indexx] > 0)
		{
		    c_imin[dd] = c_icrds[dd] - 1;
		    c_imax[dd] = c_icrds[dd];
		}
		else
		{
		    c_imin[dd] = c_icrds[dd];
		    c_imax[dd] = c_icrds[dd] + 1;
		}
	    }
	    for (kkk = c_imin[2]; kkk <= c_imax[2]; kkk++)
	    for (jjj = c_imin[1]; jjj <= c_imax[1]; jjj++)
	    for (iii = c_imin[0]; iii <= c_imax[0]; iii++)
	    {
		indexxx = d_index3d(iii,jjj,kkk,top_gmax);
		if (cft_effective_index(indexxx,c))
		{
		    //debugdan	FIXME
		    if (debugdan)
		    {
			printf("add (%d, %d, %d), index = %d.\n",
				iii, jjj, kkk, indexxx);
		    }
		    //debugdan	FIXME
		    cft_add_cftcell_index(c->cftcell,indexxx);
		}
	    }
	    pam = pam->next;
	}

	return;
}

void G_CARTESIAN::cft_set_flux_at_halft_cftcell(CFTCELL *cftc)
{
	CPAM *pam = cftc->halfts_pams;

	cftc->mflux = 0.0;
	while (pam)
	{
	    cftc->mflux += pam->polyh->mflux;
	    pam = pam->next;
	}

	return;
}

void G_CARTESIAN::cft_set_mass_at_oldt_cftcell(CFTCELL *cftc)
{
	CPAM *pam = cftc->oldts_pams;

	cftc->mass_old = 0.0;
	while (pam)
	{
	    cftc->mass_old += pam->polyh->state.dens * pam->polyh->vol;
	    pam = pam->next;
	}

	return;
}

void G_CARTESIAN::cft_set_vols_for_cftcell(CFTCELL *cftc)
{
	CPAM *pam;

	//vol_new
	cftc->vol_new = 0.0;
	pam = cftc->pams;
	while (pam)
	{
	    cftc->vol_new += pam->polyh->vol;
	    pam = pam->next;
	}

	//vol_old
	cftc->vol_old = 0.0;
	pam = cftc->oldts_pams;
	while (pam)
	{
	    cftc->vol_old += pam->polyh->vol;
	    pam = pam->next;
	}

	return;
}

bool G_CARTESIAN::cft_effective_index(
	int	index,
	CELL	*c)
{
	int i, j, k, indexx;
	COMPONENT comp;
	CPAM *pam;
	CPOLYHEDRON *polyh;

	//index is in pams
	pam = c->pams;
	comp = pam->polyh->comp;
	while (pam)
	{
	    i = pam->polyh->cell->icrds[0];
	    j = pam->polyh->cell->icrds[1];
	    k = pam->polyh->cell->icrds[2];
	    indexx = d_index3d(i,j,k,top_gmax);
	    if (index == indexx)
		return true;
	    pam = pam->next;
	}

	//index is not in pams
	//and cells_new[index] does not have polyh with comp
	polyh = cells_new[index].polyhs;
	while (polyh)
	{
	    if (polyh->comp == comp)
		return false;
	    polyh = polyh->next;
	}
	return true;
}

void G_CARTESIAN::cft_add_cftcell_index(
	CFTCELL	*cc,
	int	index)
{
	int i;
	for (i = 0; i < cc->icount; i++)
	{
	    if (cc->indices[i] == index)
		return;
	}

	if (i == 125)
	{
	    printf("ERROR in cft_add_cftcell_index(): too many indices.\n");
	    clean_up(ERROR);
	}

	cc->indices[i] = index;
	cc->icount++;

	return;
}

void G_CARTESIAN::cft_set_maps_for_cftcell(CFTCELL *cftc)
{
	int i, index;
	COMPONENT comp;
	CELL *c;
	CPOLYHEDRON *polyh;

	//debugdan	FIXME
	if (debugdan)
	{
	    printf("set maps for cftcell.\n");
	}
	//debugdan	FIXME

	comp = cftc->pams->polyh->comp;
	for (i = 0; i < cftc->icount; i++)
	{
	    index = cftc->indices[i];

	    //half time step
	    c = &(cells_halft[index]);
	    polyh = c->polyhs;
	    while (polyh)
	    {
		if (polyh->comp == comp)
		{
		    cft_add_pam(polyh,&(cftc->halfts_pams));
		    //debugdan	FIXME
		    if (debugdan)
		    {
			printf("half time step, add polyh in (%d %d %d).\n",
				polyh->cell->icrds[0],
				polyh->cell->icrds[1],
				polyh->cell->icrds[2]);
		    }
		    //debugdan	FIXME
		}
		polyh = polyh->next;
	    }

	    //old time step
	    c = &(cells_old[index]);
	    polyh = c->polyhs;
	    while (polyh)
	    {
		if (polyh->comp == comp)
		{
		    cft_add_pam(polyh,&(cftc->oldts_pams));
		    //debugdan	FIXME
		    if (debugdan)
		    {
			printf("old time step, add polyh in (%d %d %d).\n",
				polyh->cell->icrds[0],
				polyh->cell->icrds[1],
				polyh->cell->icrds[2]);
		    }
		    //debugdan	FIXME
		}
		polyh = polyh->next;
	    }
	}

	return;
}

void G_CARTESIAN::cft_add_pam(
	CPOLYHEDRON	*polyh,
	CPAM		**pamlist)
{
	CPAM *pam;

	FT_ScalarMemoryAlloc((POINTER*)&pam,sizeof(CPAM));
	pam->polyh = polyh;
	pam->targetc = NULL;
	pam->next = NULL;

	pam->next = *pamlist;
	*pamlist = pam;

	return;
}

void G_CARTESIAN::cft_init_cftcell(CELL *c)
{
	CFTCELL *cftcell;

	FT_ScalarMemoryAlloc((POINTER*)&cftcell,sizeof(CFTCELL));
	cftcell->pams = NULL;
	cftcell->halfts_pams = NULL;
	cftcell->oldts_pams = NULL;
	cftcell->icount = 0;
	cftcell->vol_old = 0.0;
	cftcell->vol_new = 0.0;
	cftcell->mass_old = 0.0;
	cftcell->mass_new = 0.0;

	c->cftcell = cftcell;

	return;
}

void G_CARTESIAN::cft_update_states()
{
	int i, j, k, index, ii, jj, kk, indexx;
	int c_imin[3], c_imax[3], vdir[3];
	double vol_old, vol_new, vol_cell;
	double mass_oldt, total_mflux, dens_new;
	double vel[3];
	CELL *c, *cc;
	CPOLYHEDRON *polyh;
	CPOLYGON *face;
	CPAM *pam;
	COMPONENT comp;
	bool debugcft = false;
	double *dens = field.dens;
	double **pdens = field.pdens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

	vol_cell = top_h[0]*top_h[1]*top_h[2];

	cells = cells_new;

	//go through all new cells
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells_new[index]);

	    if (c->merged == FALSE)
	    {
		c->polyhs->state.dens = dens[index];
		continue;
	    }

	    for (ii = 0; ii < 3; ii++)
	    {
		c_imin[ii] = imax[ii];
		c_imax[ii] = imin[ii];
		//should use old ts momn?	FIXME
		if (momn[ii][index] > 0)
		    vdir[ii] = 1;
		else
		    vdir[ii] = -1;
	    }

	    comp = c->pams->polyh->comp;
	    //comp = top_comp[index];	FIXME

	    //debugdan	FIXME
	    /*
	    if (i == 4 && j == 4)
		debugcft = true;
	    else
		debugcft = false;
	    */
	    if (debugcft)
	    {
		printf("\ni = %d, j = %d, k = %d.\n", i, j, k);
	    }
	    //debugdan	FIXME

	    //volume at new time step
	    vol_new = 0.0;
	    pam = c->pams;
	    while (pam)
	    {
		polyh = pam->polyh;
		/*
		if (polyh->comp != comp)
		{
		    pam = pam->next;
		    continue;
		}
		*/
		for (ii = 0; ii < 3; ii++)
		{
		    //IMPORTANT
		    //FIXME
		    //FIXME
		    //FIXME
		    //FIXME
		    //FIXME
		    //instead of using c_imin and c_imax, we should record exact cell index!!!	FIXME
		    //FIXME
		    //FIXME
		    //FIXME
		    //FIXME
		    //FIXME
		    //IMPORTANT
		    if (polyh->cell->icrds[ii] < c_imin[ii])
			c_imin[ii] = polyh->cell->icrds[ii];
		    if (polyh->cell->icrds[ii] > c_imax[ii])
			c_imax[ii] = polyh->cell->icrds[ii];
		}

		vol_new += polyh->vol;

		//debugdan	FIXME
		if (debugcft)
		    printf("%d %d %d newt comp %d, add vol %e.\n",
			    polyh->cell->icrds[0], polyh->cell->icrds[1], 
			    polyh->cell->icrds[2], comp, polyh->vol);

		pam = pam->next;
	    }

	    //debugdan	FIXME
	    if (debugcft)
	    {
		printf("vol_new = %e.\n", vol_new);
		printf("c_imin = (%d, %d, %d), c_imax = (%d, %d, %d).\n",
			c_imin[0], c_imin[1], c_imin[2],
			c_imax[0], c_imax[1], c_imax[2]);
	    }


	    //volume and mass at old time step
	    vol_old = 0.0;
	    mass_oldt = 0.0;
	    for (kk = c_imin[2]; kk <= c_imax[2]; kk++)
	    for (jj = c_imin[1]; jj <= c_imax[1]; jj++)
	    for (ii = c_imin[0]; ii <= c_imax[0]; ii++)
	    {
		indexx = d_index3d(ii,jj,kk,top_gmax);
		cc = &(cells_old[indexx]);

		polyh = cc->polyhs;
		while (polyh)
		{
		    //what is the polyh is merged to another cell/in another direction?
		    //check pam->mdir.	TODO
		    if (polyh->comp == comp)
		    {
			vol_old += polyh->vol;
			mass_oldt += polyh->state.dens * polyh->vol;
		    }
		    polyh = polyh->next;
		}
		pam = cc->pams;
		while (pam != NULL)
		{
		    polyh = pam->polyh;
		    if (polyh->comp != comp)
			break;
		    if (need_exterior_pam(pam,c_imin,c_imax,vdir))
		    //if (polyh->cell->icrds[0] < c_imin[0] || polyh->cell->icrds[0] > c_imax[0] ||
			//polyh->cell->icrds[1] < c_imin[1] || polyh->cell->icrds[1] > c_imax[1] ||
			//polyh->cell->icrds[2] < c_imin[2] || polyh->cell->icrds[2] > c_imax[2])
		    {
			//if we have corresponding polyh in new ts, 
			//do NOT add it.
			if (corr_polyh_in_cell(polyh,cells_new,top_gmax))
			{
			    pam = pam->next;
			    continue;
			}
			vol_old += polyh->vol;
			mass_oldt += polyh->state.dens * polyh->vol;
		    }
		    pam = pam->next;
		}
	    }

	    //debugdan	FIXME
	    if (debugcft)
		printf("%d %d %d total: vol_old = %lf, mass_oldt = %lf, old dens = %lf.\n",
			i, j, k, vol_old, mass_oldt, mass_oldt/vol_old);

	    //flux at mid time step
	    total_mflux = 0;
	    for (kk = c_imin[2]; kk <= c_imax[2]; kk++)
	    for (jj = c_imin[1]; jj <= c_imax[1]; jj++)
	    for (ii = c_imin[0]; ii <= c_imax[0]; ii++)
	    {
		indexx = d_index3d(ii,jj,kk,top_gmax);
		cc = &(cells_halft[indexx]);

		polyh = cc->polyhs;
		while (polyh)
		{
		    if (polyh->comp != comp)
		    {
			polyh = polyh->next;
			continue;
		    }

		    //debugdan	FIXME
		    if (debugcft)
		    {
			printf("%d %d %d halft polyh mflux = %e, dens flux = %e.\n",
				ii, jj, kk, polyh->mflux, polyh->mflux/polyh->vol);
		    }

		    total_mflux += polyh->mflux;
		    polyh = polyh->next;
		}
		//find polyhs which are out of range but merged to this cell
		pam = cc->pams;
		while (pam != NULL)
		{
		    polyh = pam->polyh;
		    //debugdan	FIXME
		    if (debugcft)
		    {
			printf("%d %d %d halft polyh mflux = %e, dens flux = %e.\n",
				ii, jj, kk, polyh->mflux, polyh->mflux/polyh->vol);
		    }
		    if (polyh->comp != comp)
			break;
		    if (need_exterior_pam(pam,c_imin,c_imax,vdir))
		    //if (polyh->cell->icrds[0] < c_imin[0] || polyh->cell->icrds[0] > c_imax[0] ||
			//polyh->cell->icrds[1] < c_imin[1] || polyh->cell->icrds[1] > c_imax[1] ||
			//polyh->cell->icrds[2] < c_imin[2] || polyh->cell->icrds[2] > c_imax[2])
		    {
			if (corr_polyh_in_cell(polyh,cells_new,top_gmax))
			{
			    pam = pam->next;
			    continue;
			}
			total_mflux += polyh->mflux;
		    }
		    pam = pam->next;
		}
	    }

	    dens_new = (mass_oldt + total_mflux) / vol_new;

	    //debugdan	FIXME
	    if (debugcft)
		printf("%d %d %d dens_new = %e.\n", i, j, k, dens_new);
	    if (dens_new > 3.3)
	    {
		printf("%d %d %d dens_new = %e.\n", i, j, k, dens_new);
		printf("vol_old = %e, old dens = %e, vol_new = %e.\n",
			vol_old, mass_oldt/vol_old, vol_new);
		printf("c_imin: %d %d %d. c_imax: %d %d %d.\n",
			c_imin[0], c_imin[1], c_imin[2],
			c_imax[0], c_imax[1], c_imax[2]);
		//vol_old
		printf("%d %d %d vol_old:\n", c->icrds[0], c->icrds[1], c->icrds[2]);
		for (kk = c_imin[2]; kk <= c_imax[2]; kk++)
		for (jj = c_imin[1]; jj <= c_imax[1]; jj++)
		for (ii = c_imin[0]; ii <= c_imax[0]; ii++)
		{
		    indexx = d_index3d(ii,jj,kk,top_gmax);
		    cc = &(cells_old[indexx]);
		    polyh = cc->polyhs;
		    while (polyh)
		    {
			if (polyh->comp == comp)
			{
			    printf("%d %d %d + %e.\n", ii, jj, kk, polyh->vol);
			    print_cpolyh(polyh);
			}
			polyh = polyh->next;
		    }
		    pam = cc->pams;
		    while (pam != NULL)
		    {
			polyh = pam->polyh;
			if (polyh->comp != comp)
			    break;
			if (need_exterior_pam(pam,c_imin,c_imax,vdir))
			{
			    printf("%d %d %d pam + %e. vdir = %d %d %d.\n",
				    polyh->cell->icrds[0], polyh->cell->icrds[1], polyh->cell->icrds[2],
				    polyh->vol, vdir[0], vdir[1], vdir[2]);
			    print_cpolyh(polyh);
			}
			pam = pam->next;
		    }
		}
		//vol_new
		printf("%d %d %d vol_new:\n", c->icrds[0], c->icrds[1], c->icrds[2]);
		pam = c->pams;
		while (pam)
		{
		    printf(" + %e\n", pam->polyh->vol);
		    print_cpolyh(pam->polyh);
		    pam = pam->next;
		}
	    }
	    /*
	    if (j == 4 && k == 20)
	    {
		debugvol[i] = vol_new;
	    }
	    if (j == 4 && k == 20)
		printf("i = %d, dens_new = %lf, total_mflux = %e, Dvol = %e.\n",
			i, dens_new, total_mflux, debugvol[i]-debugvol[i-1]);
	    */

	    pam = c->pams;
	    while (pam)
	    {
		pam->polyh->state.dens = dens_new;
		pam = pam->next;
	    }
	}

	double cmass;
	STATE state;
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells_new[index]);

	    if (c->merged == FALSE)
	    {
		continue;
	    }

	    //debugdan	FIXME
	    /*
	    if (i == 4 && j == 4)
		debugcft = true;
	    else
		debugcft = false;
	    */
	    if (debugcft)
	    {
		printf("\ni = %d, j = %d, k = %d.\n", i, j, k);
	    }
	    //debugdan	FIXME

	    cmass = 0;
	    polyh = c->polyhs;
	    while (polyh)
	    {
		//debugdan	FIXME
		if (debugcft)
		{
		    printf("new ts polyh comp = %d, vol = %lf, dens = %lf.\n",
			    polyh->comp, polyh->vol, polyh->state.dens);
		}
		cmass += polyh->state.dens*polyh->vol;
		polyh = polyh->next;
	    }

	    //KEY
	    //keep vel
	    for (ii = 0; ii < 3; ii ++)
		vel[ii] = momn[ii][index]/dens[index];
		//vel[ii] = momn[ii][index]/(dens[index]*vol_cell);
	    //update dens
	    dens[index] = cmass/vol_cell;
	    //update pdens???	TODO
	    //update momn
	    for (ii = 0; ii < 3; ii++)
		momn[ii][index] = vel[ii]*dens[index];
		//momn[ii][index] = vel[ii]*cmass;
	    //update engy
	    comp = top_comp[index];
	    state.dim = eqn_params->dim;
	    state.pres = pres[index];
	    state.dens = dens[index];
	    for (ii = 0; ii < eqn_params->n_comps; ii++)
		state.pdens[ii] = pdens[ii][index];	//FIXME
	    for (ii = 0; ii < 3; ii++)
		state.momn[ii] = momn[ii][index];
	    state.eos = &(eqn_params->eos[comp]);
	    engy[index] = EosEnergy(&state);

	    //update polyhs states of cells_new
	    polyh = c->polyhs;
	    while (polyh)
	    {
		polyh->state.dim = eqn_params->dim;
		polyh->state.pres = pres[index];
		for (ii = 0; ii < 3; ii++)
		    polyh->state.momn[ii] = vel[ii]*polyh->state.dens;
		polyh->state.eos = &(eqn_params->eos[polyh->comp]);
		//pdens	FIXME
		//for current tests, we don't need to consider a changing gamma,
		//so pdens is not important.
		//if we have a changing gamma, weno scheme needs to be modified.
		if (eqn_params->multi_comp_non_reactive == YES)
		{
		    for (ii = 0; ii < eqn_params->n_comps; ii++)
			polyh->state.pdens[ii] = pdens[ii][index];
		}
		polyh->state.engy = EosEnergy(&(polyh->state));
		//debugdan	FIXME
		/*
		if (i == 4 && j == 4)
		{
		    printf("k = %d, polyh states: %e, %e, %e, %e.\n",
			    k, polyh->state.dens, polyh->state.pres, polyh->state.momn[2],
			    polyh->state.engy);
		}
		*/
		//debugdan	FIXME
		polyh = polyh->next;
	    }

	    //debugdan	FIXME
	    //debug dens
	    if (dens[index] > 3.3)
	    {
		printf("%d %d %d incorrect dens %e.\n",
			i, j, k, dens[index]);
		polyh = c->polyhs;
		while (polyh)
		{
		    printf("polyh dens = %e. vol = %e.\n",
			    polyh->state.dens, polyh->vol);
		    polyh = polyh->next;
		}
	    }
	    /*
	    if (debugcft)
		printf("newdens = %lf.\n\n", field.dens[index]);
	    */
	}

	return;
}

bool corr_polyh_in_cell(
	CPOLYHEDRON	*polyh,
	CELL		*cells,
	int		*gmax)
{
	int i, j, k, index;
	COMPONENT comp = polyh->comp;
	CELL *c;
	CPOLYHEDRON *plh;

	i = polyh->cell->icrds[0];
	j = polyh->cell->icrds[1];
	k = polyh->cell->icrds[2];
	index = d_index3d(i,j,k,gmax);
	c = &(cells[index]);

	plh = c->polyhs;
	while (plh)
	{
	    if (plh->comp == comp)
		return true;
	    plh = plh->next;
	}

	return false;
}

bool need_exterior_pam(
	CPAM	*pam,
	int	*c_imin,
	int	*c_imax,
	int	*vdir)
{
	int i;
	int *icrds = pam->polyh->cell->icrds;

	for (i = 0; i < 3; i++)
	{
	    if (icrds[i] < c_imin[i])
	    {
		if (vdir[i] == 1)
		    return true;
		else
		    return false;
	    }
	    else if (icrds[i] > c_imax[i])
	    {
		if (vdir[i] == -1)
		    return true;
		else
		    return false;
	    }
	}

	return false;
}

void G_CARTESIAN::cft_newts_cut_cells()
{
	CELL *tmp;

	tmp = cells_old;
	cells_old = cells_new;
	cells_new = tmp;

	cft_free_cells(HALFTS);
	cft_free_cells(NEWTS);
	//cft_reset_halft_and_new_cells();

	return;
}

void G_CARTESIAN::cft_free_cells(TS_LEVEL ts)
{
	switch (ts)
	{
	    case OLDTS:
		cells = cells_old;
		break;
	    case HALFTS:
		cells = cells_halft;
		break;
	    case NEWTS:
		cells = cells_new;
		break;
	    default:
		printf("ERROR in cft_free_cells: unknown ts level.\n");
		clean_up(ERROR);
	}

	int i, j, k, index;
	CELL *c;

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);
	    cft_reset_cell(c);
	}

	free(cells);
	cells = NULL;
	return;
}

/*
void G_CARTESIAN::cft_reset_halft_and_new_cells()
{
	int i, j, k, index;
	CELL *c;

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells_halft[index]);
	    cft_reset_cell(c);

	    c = &(cells_new[index]);
	    cft_reset_cell(c);
	}

	return;
}
*/

void G_CARTESIAN::cft_reset_cell(CELL *c)
{
	cft_clean_edges(c);
	cft_clean_faces(c);

	cft_free_ctri_list(&c->ctris);
	cft_free_polyg_list(&c->ctri_polygs);
	cft_free_polyg_list(&c->cf_polygs);
	cft_free_scp_list(&c->scpocs);
	cft_free_scp_list(&c->scpics);
	cft_free_polyh_list(&c->polyhs);
	cft_free_pam_list(&c->pams);

	return;
}

void G_CARTESIAN::cft_clean_edges(CELL *c)
{
	for (int i = 0; i < 12; i++)
	    cft_free_p_list(&c->edges[i].crxps);
}

void G_CARTESIAN::cft_clean_faces(CELL *c)
{
}

void G_CARTESIAN::cft_free_ctri_list(CTRI **ctri)
{
	if (!*ctri)
	    return;

	CTRI *node = *ctri;
	while (node)
	{
	    *ctri = node->next;
	    cft_free_ctri(&node);
	    node = *ctri;
	}

	*ctri = NULL;
	return;
}

void G_CARTESIAN::cft_free_polyg_list(CPOLYGON **polyg)
{
	if (!*polyg)
	    return;

	CPOLYGON *node = *polyg;
	while (node)
	{
	    *polyg = node->next;
	    cft_free_polyg(&node);
	    free(node);
	    node = *polyg;
	}

	*polyg = NULL;
	return;
}

void G_CARTESIAN::cft_free_p_list(CPOINT **p)
{
	if (!*p)
	    return;

	CPOINT *node = *p;
	while (node)
	{
	    *p = node->next;
	    free(node);
	    node = *p;
	}

	*p = NULL;
	return;
}

void G_CARTESIAN::cft_free_edge_list(CEDGE **edge)
{
	if (!*edge)
	    return;

	CEDGE *node = *edge;
	while (node)
	{
	    *edge = node->next;
	    cft_free_edge(&node);
	    node = *edge;
	}

	*edge = NULL;
	return;
}

void G_CARTESIAN::cft_free_scp_list(SETOFCPOLYGS **scp)
{
	if (!*scp)
	    return;

	SETOFCPOLYGS *node = *scp;
	while (node)
	{
	    *scp = node->next;
	    cft_free_scp(&node);
	    node = *scp;
	}

	*scp = NULL;
	return;
}

void G_CARTESIAN::cft_free_bdry_list(CBOUNDARY **bdry)
{
	if (!*bdry)
	    return;

	CBOUNDARY *node = *bdry;
	while (node)
	{
	    *bdry = node->next;
	    cft_free_bdry(&node);
	    node = *bdry;
	}

	*bdry = NULL;
	return;
}

void G_CARTESIAN::cft_free_polyh_list(CPOLYHEDRON **polyh)
{
	if (!*polyh)
	    return;

	CPOLYHEDRON *node = *polyh;
	while (node)
	{
	    *polyh = node->next;
	    cft_free_polyh(&node);
	    node = *polyh;
	}

	*polyh = NULL;
	return;
}

void G_CARTESIAN::cft_free_pam_list(CPAM **pam)
{
	if (!*pam)
	    return;

	CPAM *node = *pam;
	while (node)
	{
	    *pam = node->next;
	    cft_free_pam(&node);
	    node = *pam;
	}

	*pam = NULL;
	return;
}

void G_CARTESIAN::cft_free_pam(CPAM **pam)
{
	if (!*pam)
	    return;

	cft_free_polyh(&(*pam)->polyh);

	*pam = NULL;
	return;
}

void G_CARTESIAN::cft_free_polyh(CPOLYHEDRON **polyh)
{
	if (!*polyh)
	    return;

	cft_free_polyg_list(&(*polyh)->faces);
	cft_free_bdry_list(&(*polyh)->boundaries);

	*polyh = NULL;
	return;
}

void G_CARTESIAN::cft_free_bdry(CBOUNDARY **bdry)
{
	if (!*bdry)
	    return;

	cft_free_edge_list(&(*bdry)->edges);
	free(*bdry);

	*bdry = NULL;
	return;
}

void G_CARTESIAN::cft_free_scp(SETOFCPOLYGS **scp)
{
	if (!*scp)
	    return;

	cft_free_polyg_list(&(*scp)->polygs);
	cft_free_bdry_list(&(*scp)->boundaries);
	free(*scp);

	*scp = NULL;
	return;
}

void G_CARTESIAN::cft_free_edge(CEDGE **edge)
{
	if (!*edge)
	    return;

	cft_free_p_list(&(*edge)->crxps);
	free(*edge);

	*edge = NULL;
	return;
}

void G_CARTESIAN::cft_free_ctri(CTRI **ctri)
{
	if (!*ctri)
	    return;

	cft_free_polyg(&(*ctri)->polyg);
	free(*ctri);

	*ctri = NULL;
	return;
}

void G_CARTESIAN::cft_free_polyg(CPOLYGON **polyg)
{
	if (!*polyg)
	    return;

	cft_free_p_list(&(*polyg)->vertices);
	cft_free_edge_list(&(*polyg)->edges);
	cft_free_edge_list(&(*polyg)->undir_edges);
	free(*polyg);

	*polyg = NULL;
	return;
}

void G_CARTESIAN::cft_check_mass()
{
	int i, j, k, index;
	CELL *c;
	double cvol, tm1, tm2;

	printf("This function is for test ONLY.\n");
	printf("cft_check_mass():\n");

	cells = cells_new;
	tm1 = 0.0;
	tm2 = 0.0;
	cvol = top_h[0]*top_h[1]*top_h[2];

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);

	    if (c->merged == FALSE)
	    {
		if (top_comp[index] == GAS_COMP1)
		{
		    tm1 += cvol * field.dens[index];
		}
		else if (top_comp[index] == GAS_COMP2)
		{
		    tm2 += cvol * field.dens[index];
		}
		continue;
	    }
	    else
	    {
		CPOLYHEDRON *polyh = c->polyhs;
		while (polyh)
		{
		    if (polyh->comp == GAS_COMP1)
		    {
			tm1 += polyh->state.dens * polyh->vol;
		    }
		    else if (polyh->comp == GAS_COMP2)
		    {
			tm2 += polyh->state.dens * polyh->vol;
		    }
		    polyh = polyh->next;
		}
	    }
	}

	printf("total mass for comp 1 = %e.\n", tm1);
	printf("total mass for comp 2 = %e.\n", tm2);

	return;
}

void G_CARTESIAN::ncft_check_mass()
{
	int i, j, k, index;
	CELL *c;
	double cvol, tm1, tm2;

	printf("This function is for test ONLY.\n");
	printf("ncft_check_mass():\n");

	tm1 = 0.0;
	tm2 = 0.0;
	cvol = top_h[0]*top_h[1]*top_h[2];

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);

	    if (top_comp[index] == GAS_COMP1)
	    {
		tm1 += cvol * field.dens[index];
	    }
	    else if (top_comp[index] == GAS_COMP2)
	    {
		tm2 += cvol * field.dens[index];
	    }
	}

	printf("total mass for comp 1 = %e.\n", tm1);
	printf("total mass for comp 2 = %e.\n", tm2);

	return;
}

void print_cpolyh(CPOLYHEDRON *polyh)
{
	CPOLYGON *face;
	printf("Polyhedron:\n");
	if (polyh->iscell)
	{
	    printf("is a cell.\n");
	    return;
	}
	face = polyh->faces;
	while (face)
	{
	    print_cpolygon(face);
	    face = face->next;
	}
}

void print_cpolygon(CPOLYGON *pg)
{
	//format it for Mathematica
	CPOINT *p;
	printf("Polygon[{\n");
	p = pg->vertices;
	while(p)
	{
	    print_cpoint(p);
	    if (p->next)
		printf(",");
	    printf("\n");
	    p = p->next;
	}
	printf("}],\n");
	return;
}

void print_cpoint(CPOINT *p)
{
	//printf("(%lf, %lf, %lf)", p->crds[0], p->crds[1], p->crds[2]);
	printf("{%.12g, %.12g, %.12g}", p->crds[0], p->crds[1], p->crds[2]);
}

