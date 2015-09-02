/*******************************************************************
 * 			iFcartsn3d.cpp	
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

//--------------------------------------------------------------------------
// 		   Incompress_Solver_Smooth_3D_Cartesian	
//--------------------------------------------------------------------------

void Incompress_Solver_Smooth_3D_Cartesian::printInteriorVelocity(char *out_name, bool binary)
{
        int   i,j,k,l,index,totalpoints;
        double coord_x,coord_y,coord_z;
        double vel_x, vel_y, vel_z;
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
        if(iFparams->movie_option->plot_pres || iFparams->movie_option->plot_velo || iFparams->movie_option->plot_vort)
        {
	    outfile = fopen(filename,"w");
            fprintf(outfile,"# vtk DataFile Version 3.0\n");
            fprintf(outfile,"States of the whole computational domain\n");
            fprintf(outfile,"ASCII\n");
            fprintf(outfile,"DATASET STRUCTURED_GRID\n");
            int pointsx, pointsy, pointsz;
        
	    pointsz = top_gmax[2] + 1;
            pointsy = top_gmax[1] + 1;
            pointsx = top_gmax[0] + 1;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile, "POINTS %d double\n", totalpoints);
            for(i = 0; i <= top_gmax[0]; ++i)
            for(j = 0; j <= top_gmax[1]; ++j)
            for(k = 0; k <= top_gmax[2]; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                coord_x = cell_center[index].m_coords[0];
                coord_y = cell_center[index].m_coords[1];
                coord_z = cell_center[index].m_coords[2];
                fprintf(outfile, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }
            fprintf(outfile, "POINT_DATA %i\n", totalpoints);
            if(iFparams->movie_option->plot_velo)
            {
                fprintf(outfile, "VECTORS velocity double\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    vel_x = cell_center[index].m_state.m_U[0];
                    vel_y = cell_center[index].m_state.m_U[1];
                    vel_z = cell_center[index].m_state.m_U[2];
                    fprintf(outfile, "%.16g %.16g %16g\n",vel_x, vel_y, vel_z);
                }
            }
            if(iFparams->movie_option->plot_pres)
            {
                fprintf(outfile, "SCALARS pressure double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_P);
                } 
            }
            if(iFparams->movie_option->plot_vort)
            {
            }
            fclose(outfile);
        }
}

void Incompress_Solver_Smooth_3D_Cartesian::printExpandedMesh(char *out_name, bool binary)
{
        int   i,j,k,l,index,totalpoints;
        double coord_x,coord_y,coord_z;
        int first;
        char filename[200];
        FILE *outfile;
        sprintf(filename,"%s/mesh-ts%s",out_name,
                right_flush(front->step,7));
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif
        sprintf(filename,"%s-liquid.vtk",filename);
        i = 0;

	    outfile = fopen(filename,"w");
            fprintf(outfile,"# vtk DataFile Version 3.0\n");
            fprintf(outfile,"Mesh of the topological grid\n");
            fprintf(outfile,"ASCII\n");
            fprintf(outfile,"DATASET STRUCTURED_GRID\n");
            int pointsx, pointsy, pointsz;
	    double xmin, xmax, ymin, ymax, zmin, zmax;
	    index = d_index3d(0,0,0,top_gmax);
	    xmin = cell_center[index].m_coords[0] - top_h[0]/2.0;
	    ymin = cell_center[index].m_coords[1] - top_h[1]/2.0;
	    zmin = cell_center[index].m_coords[2] - top_h[2]/2.0;
	    index = d_index3d(top_gmax[0],top_gmax[1],top_gmax[2],top_gmax);
	    xmax = cell_center[index].m_coords[0] + top_h[0]/2.0;
	    ymax = cell_center[index].m_coords[1] + top_h[1]/2.0;
	    zmax = cell_center[index].m_coords[2] + top_h[2]/2.0;
        
	    pointsz = top_gmax[2] + 2;
            pointsy = top_gmax[1] + 2;
            pointsx = top_gmax[0] + 2;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile, "POINTS %d double\n", totalpoints);
            for(i = 0; i <= top_gmax[0]+1; ++i)
            for(j = 0; j <= top_gmax[1]+1; ++j)
            for(k = 0; k <= top_gmax[2]+1; ++k)
            {
                coord_x = xmin + i*top_h[0];
                coord_y = ymin + j*top_h[1];
                coord_z = zmin + k*top_h[2];
                fprintf(outfile, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }

            fclose(outfile);
}

void Incompress_Solver_Smooth_3D_Cartesian::computeAdvection(void)
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
	    index00 = d_index3d(i-1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
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
			        &intfc_state,&hs,crx_coords,m_t_old) &&
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
			        &intfc_state,&hs,crx_coords,m_t_old) &&
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
			        &intfc_state,&hs,crx_coords,m_t_old) &&
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
			        &intfc_state,&hs,crx_coords,m_t_old) &&
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
			        &intfc_state,&hs,crx_coords,m_t_old) &&
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



void Incompress_Solver_Smooth_3D_Cartesian::
	compDiffWithSmoothProperty_1st_coupled(void)
{
        COMPONENT comp;
        int index,index_nb[18],size;
        int I,I_nb[18];
        double coords[MAXD],crx_coords[MAXD];
	double coeff[18],mu[6],mu_edge[6],mu0,rho,rhs,U0_nb[18],U1_nb[18],U2_nb[18];
	int flag[6]; //denote whether this is dirichlet or neumann boundary
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

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 15, 15);
	// 7u + 4v + 4w for the first equation
	// 7v + 4u + 4w for the second equation
	// 7w + 4u + 4v for the third equation

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
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	//xy cut neighbours
	    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
	    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
	    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
	    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	//yz cut neighbours
	    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
	    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
	    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
	    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);
	//xz cut neighbours
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
	
	//xy cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
	    I_nb[7] = ijk_to_I[i+1][j-1][k];
	    I_nb[8] = ijk_to_I[i+1][j+1][k];
	    I_nb[9] = ijk_to_I[i-1][j+1][k];
	//yz cut neighbours
	    I_nb[10] = ijk_to_I[i][j-1][k-1];
	    I_nb[11] = ijk_to_I[i][j+1][k-1];
	    I_nb[12] = ijk_to_I[i][j+1][k+1];
	    I_nb[13] = ijk_to_I[i][j-1][k+1];
	//xz cut neighbours
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
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    flag[nb] = 1;
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);

		    if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			wave_type(hs) == NEUMANN_BOUNDARY)
		    {
			mu[nb] = mu0;
			mu_edge[nb] = mu0;
		    }
			
		    else
		    {
			    mu[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
			    mu_edge[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
		}
                    else
		    {
			flag[nb] = 0;
                    	U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
			U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
			U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
			mu_edge[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
			mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
	    }

	    //xy cut neighbour cells
	    
	    if (flag[0] == 1 && flag[2] == 1) //neighbourcell 6
	    {
		U0_nb[6] = 0.5*(U0_nb[0] + U0_nb[2]);
		U1_nb[6] = 0.5*(U1_nb[0] + U1_nb[2]);
		U2_nb[6] = 0.5*(U2_nb[0] + U2_nb[2]);
	    }
	    else
	    {
		U0_nb[6] = cell_center[index_nb[6]].m_state.m_U[0];
		U1_nb[6] = cell_center[index_nb[6]].m_state.m_U[1];
		U2_nb[6] = cell_center[index_nb[6]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[2] == 1) //neighbourcell 7 
	    {
		U0_nb[7] = 0.5*(U0_nb[1] + U0_nb[2]);
		U1_nb[7] = 0.5*(U1_nb[1] + U1_nb[2]);
		U2_nb[7] = 0.5*(U2_nb[1] + U2_nb[2]);
	    }
	    else
	    {
		U0_nb[7] = cell_center[index_nb[7]].m_state.m_U[0];
		U1_nb[7] = cell_center[index_nb[7]].m_state.m_U[1];
		U2_nb[7] = cell_center[index_nb[7]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[3] == 1) //neighbourcell 8 
	    {
		U0_nb[8] = 0.5*(U0_nb[1] + U0_nb[3]);
		U1_nb[8] = 0.5*(U1_nb[1] + U1_nb[3]);
		U2_nb[8] = 0.5*(U2_nb[1] + U2_nb[3]);
	    }
	    else
	    {
		U0_nb[8] = cell_center[index_nb[8]].m_state.m_U[0];
		U1_nb[8] = cell_center[index_nb[8]].m_state.m_U[1];
		U2_nb[8] = cell_center[index_nb[8]].m_state.m_U[2];
	    }
	    if (flag[0] == 1 && flag[3] == 1) //neighbourcell 9 
	    {
		U0_nb[9] = 0.5*(U0_nb[0] + U0_nb[3]);
		U1_nb[9] = 0.5*(U1_nb[0] + U1_nb[3]);
		U2_nb[9] = 0.5*(U2_nb[0] + U2_nb[3]);
	    }
	    else
	    {
		U0_nb[9] = cell_center[index_nb[9]].m_state.m_U[0];
		U1_nb[9] = cell_center[index_nb[9]].m_state.m_U[1];
		U2_nb[9] = cell_center[index_nb[9]].m_state.m_U[2];
	    }
	    //yz cut neighbours
	    if (flag[2] == 1 && flag[4] == 1) //neighbourcell 10 
	    {
		U0_nb[10] = 0.5*(U0_nb[2] + U0_nb[4]);
		U1_nb[10] = 0.5*(U1_nb[2] + U1_nb[4]);
		U2_nb[10] = 0.5*(U2_nb[2] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[10] = cell_center[index_nb[10]].m_state.m_U[0];
		U1_nb[10] = cell_center[index_nb[10]].m_state.m_U[1];
		U2_nb[10] = cell_center[index_nb[10]].m_state.m_U[2];
	    }
	    if (flag[3] == 1 && flag[4] == 1) //neighbourcell 11 
	    {
		U0_nb[11] = 0.5*(U0_nb[3] + U0_nb[4]);
		U1_nb[11] = 0.5*(U1_nb[3] + U1_nb[4]);
		U2_nb[11] = 0.5*(U2_nb[3] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[11] = cell_center[index_nb[11]].m_state.m_U[0];
		U1_nb[11] = cell_center[index_nb[11]].m_state.m_U[1];
		U2_nb[11] = cell_center[index_nb[11]].m_state.m_U[2];
	    }
	    if (flag[3] == 1 && flag[5] == 1) //neighbourcell 12 
	    {
		U0_nb[12] = 0.5*(U0_nb[3] + U0_nb[5]);
		U1_nb[12] = 0.5*(U1_nb[3] + U1_nb[5]);
		U2_nb[12] = 0.5*(U2_nb[3] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[12] = cell_center[index_nb[12]].m_state.m_U[0];
		U1_nb[12] = cell_center[index_nb[12]].m_state.m_U[1];
		U2_nb[12] = cell_center[index_nb[12]].m_state.m_U[2];
	    }
	    if (flag[2] == 1 && flag[5] == 1) //neighbourcell 13
	    {
		U0_nb[13] = 0.5*(U0_nb[2] + U0_nb[5]);
		U1_nb[13] = 0.5*(U1_nb[2] + U1_nb[5]);
		U2_nb[13] = 0.5*(U2_nb[2] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[13] = cell_center[index_nb[13]].m_state.m_U[0];
		U1_nb[13] = cell_center[index_nb[13]].m_state.m_U[1];
		U2_nb[13] = cell_center[index_nb[13]].m_state.m_U[2];
	    }
	    //xz cut neighbours
	    if (flag[0] == 1 && flag[4] == 1) //neighbourcell 14
	    {
		U0_nb[14] = 0.5*(U0_nb[0] + U0_nb[4]);
		U1_nb[14] = 0.5*(U1_nb[0] + U1_nb[4]);
		U2_nb[14] = 0.5*(U2_nb[0] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[14] = cell_center[index_nb[14]].m_state.m_U[0];
		U1_nb[14] = cell_center[index_nb[14]].m_state.m_U[1];
		U2_nb[14] = cell_center[index_nb[14]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[4] == 1) //neighbourcell 15
	    {
		U0_nb[15] = 0.5*(U0_nb[1] + U0_nb[4]);
		U1_nb[15] = 0.5*(U1_nb[1] + U1_nb[4]);
		U2_nb[15] = 0.5*(U2_nb[1] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[15] = cell_center[index_nb[15]].m_state.m_U[0];
		U1_nb[15] = cell_center[index_nb[15]].m_state.m_U[1];
		U2_nb[15] = cell_center[index_nb[15]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[5] == 1) //neighbourcell 16
	    {
		U0_nb[16] = 0.5*(U0_nb[1] + U0_nb[5]);
		U1_nb[16] = 0.5*(U1_nb[1] + U1_nb[5]);
		U2_nb[16] = 0.5*(U2_nb[1] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[16] = cell_center[index_nb[16]].m_state.m_U[0];
		U1_nb[16] = cell_center[index_nb[16]].m_state.m_U[1];
		U2_nb[16] = cell_center[index_nb[16]].m_state.m_U[2];
	    }
	    if (flag[0] == 1 && flag[5] == 1) //neighbourcell 17
	    {
		U0_nb[17] = 0.5*(U0_nb[0] + U0_nb[5]);
		U1_nb[17] = 0.5*(U1_nb[0] + U1_nb[5]);
		U2_nb[17] = 0.5*(U2_nb[0] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[17] = cell_center[index_nb[17]].m_state.m_U[0];
		U1_nb[17] = cell_center[index_nb[17]].m_state.m_U[1];
		U2_nb[17] = cell_center[index_nb[17]].m_state.m_U[2];
	    }



            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


    //Setting the coeffecients for the first equation


	    coeff[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    coeff[6] =  0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[8] =  0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);
	    coeff[9] = -0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);

	    coeff[14] =  0.5*m_dt/rho * mu[4]/(4*top_h[0]*top_h[2]);
	    coeff[15] = -0.5*m_dt/rho * mu[4]/(4*top_h[0]*top_h[2]);
	    coeff[16] =  0.5*m_dt/rho * mu[5]/(4*top_h[0]*top_h[2]);
	    coeff[17] = -0.5*m_dt/rho * mu[5]/(4*top_h[0]*top_h[2]);

	    solver.Set_A(I*3,I*3,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[0];

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
	    for (nb = 6; nb < 10; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }
	    for (nb = 14; nb < 18; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;

	    solver.Set_b(I*3, rhs);

	    /************************************************************************/


    //Setting the coeffecients for the second equation


	    coeff[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    coeff[6] =  0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[8] =  0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[9] = -0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);

	    coeff[10] =  0.5*m_dt/rho * mu[4]/(4*top_h[1]*top_h[2]);
	    coeff[11] = -0.5*m_dt/rho * mu[4]/(4*top_h[1]*top_h[2]);
	    coeff[12] =  0.5*m_dt/rho * mu[5]/(4*top_h[1]*top_h[2]);
	    coeff[13] = -0.5*m_dt/rho * mu[5]/(4*top_h[1]*top_h[2]);

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
	    for (nb = 6; nb < 10; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 10; nb < 14; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;

	    solver.Set_b(I*3+1, rhs);

	    /************************************************************************/

    //Setting the coeffecients for the third equation


	    coeff[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff[4] = m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff[5] = m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    coeff[14] =  0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[2]);
	    coeff[15] = -0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[2]);
	    coeff[16] =  0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[2]);
	    coeff[17] = -0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[2]);

	    coeff[10] =  0.5*m_dt/rho * mu[2]/(4*top_h[1]*top_h[2]);
	    coeff[11] = -0.5*m_dt/rho * mu[3]/(4*top_h[1]*top_h[2]);
	    coeff[12] =  0.5*m_dt/rho * mu[3]/(4*top_h[1]*top_h[2]);
	    coeff[13] = -0.5*m_dt/rho * mu[2]/(4*top_h[1]*top_h[2]);

	    solver.Set_A(I*3+2,I*3+2,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[2];

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
	    for (nb = 14; nb < 18; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 10; nb < 14; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

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

	//if (rel_residual > 1)
	//{
	//    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	//    solver.Reset_x();
	//    solver.Solve_GMRES();
	//    solver.GetNumIterations(&num_iter);
	//    solver.GetFinalRelativeResidualNorm(&rel_residual);
	//}

	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_1st_coupled: "
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
}       /* end compDiffWithSmoothProperty3d_1st_coupled */

void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_2nd_coupled(void)
{
        COMPONENT comp;
        int index,index_nb[18],size;
        int I,I_nb[18];
        double coords[MAXD],crx_coords[MAXD];
	double coeff0[6],coeff1[6],coeff2[6],coeff_temp;
	double mu[6],mu_edge[6],mu0,rho,rhs;
	double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
	int flag[6]; //denote whether this is dirichlet or neumann boundary
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

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 30, 30);
	// 7u + 9v + 9w for the first equation
	// 7v + 9u + 9w for the second equation
	// 7w + 9u + 9v for the third equation

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
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	//xy cut neighbours
	    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
	    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
	    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
	    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	//yz cut neighbours
	    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
	    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
	    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
	    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);
	//xz cut neighbours
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
	
	//xy cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
	    I_nb[7] = ijk_to_I[i+1][j-1][k];
	    I_nb[8] = ijk_to_I[i+1][j+1][k];
	    I_nb[9] = ijk_to_I[i-1][j+1][k];
	//yz cut neighbours
	    I_nb[10] = ijk_to_I[i][j-1][k-1];
	    I_nb[11] = ijk_to_I[i][j+1][k-1];
	    I_nb[12] = ijk_to_I[i][j+1][k+1];
	    I_nb[13] = ijk_to_I[i][j-1][k+1];
	//xz cut neighbours
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
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];
	    U2_center = cell_center[index].m_state.m_U[2];

            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    flag[nb] = 1;
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);

		    if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			wave_type(hs) == NEUMANN_BOUNDARY)
		    {
			mu[nb] = mu0;
			mu_edge[nb] = mu0;
		    }
		    else
		    {
			mu[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
			mu_edge[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
		}
		else
		{
		    flag[nb] = 0;
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		    mu_edge[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
		    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		}
	    }

	    for (nb = 6; nb < 18; nb++) //cut corner values, interior
	    {
		if (I_nb[nb] != -1)
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }


	    // source term
            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


	    //Setting the coeffecients for U0 in the first equation


	    coeff0[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff0[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff0[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff0[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff0[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff0[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    solver.Add_A(I*3,I*3,1.0);
	    rhs = U0_center;

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Add_A(I*3,I_nb[nb]*3,-coeff0[nb]);
		    rhs += coeff0[nb]*U0_nb[nb];
		}
		else
		{
		    coeff0[nb] = 2.0*coeff0[nb];
		    rhs += 2.0*coeff0[nb]*U0_nb[nb];
		}
		solver.Add_A(I*3,I*3,coeff0[nb]);
		rhs -= coeff0[nb]*U0_center;
	    }

	    //set the coefficients for U1 in the first equation (mu*v_x)_y
	    //traverse the four corners
	    //(i-1/2,j-1/2,k),(i+1/2,j-1/2,k),(i+1/2,j+1/2,k),(i-1/2,j+1/2,k)
	    
	    //corner (i-1/2,j-1/2,k)

	    coeff_temp = mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[0] == 1 && flag[2] == 0)
		rhs += coeff_temp*U1_nb[0];
	    else if(flag[0] == 0 && flag[2] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else if(flag[0] == 1 && flag[2] == 1)
		rhs += coeff_temp*U1_nb[0];
	    else {
		rhs += coeff_temp*(U1_nb[0]+U1_nb[2]+U1_nb[6]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[6]*3+1, -coeff_temp/8.0);
	    }
		    
	    //corner (i+1/2,j-1/2,k)

	    coeff_temp = -mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[2] == 1 && flag[1] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if(flag[2] == 0 && flag[1] == 1)
		rhs += coeff_temp*U1_nb[1];
	    else if(flag[2] == 1 && flag[1] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[1]+U1_nb[7]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[7]*3+1, -coeff_temp/8.0);
	    }

	    
	    //corner (i+1/2,j+1/2,k)

	    coeff_temp = mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[1] == 1 && flag[3] == 0)
		rhs += coeff_temp*U1_nb[1];
	    else if(flag[1] == 0 && flag[3] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else if(flag[1] == 1 && flag[3] == 1)
		rhs += coeff_temp*U1_nb[1];
	    else {
		rhs += coeff_temp*(U1_nb[1]+U1_nb[3]+U1_nb[8]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[3]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[8]*3+1, -coeff_temp/8.0);
	    }

    
	    //corner (i-1/2,j+1/2,k)

	    coeff_temp = -mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[3] == 1 && flag[0] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if(flag[3] == 0 && flag[0] == 1)
		rhs += coeff_temp*U1_nb[0];
	    else if(flag[3] == 1 && flag[0] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[0]+U1_nb[9]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[3]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[9]*3+1, -coeff_temp/8.0);
	    }


	    //set the coefficients for U2 in the first equation (mu*w_x)_z
	    //traverse the four corners
	    //(i-1/2,j,k-1/2),(i+1/2,j,k-1/2),(i+1/2,j,k+1/2),(i-1/2,j,k+1/2)
	    
	    //corner (i-1/2,j,k-1/2)

	    coeff_temp = mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[0] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[0];
	    else if(flag[0] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if(flag[0] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[0];
	    else {
		rhs += coeff_temp*(U2_nb[0]+U2_nb[4]+U2_nb[14]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[14]*3+2, -coeff_temp/8.0);
	    }
		    
	    //corner (i+1/2,j,k-1/2)

	    coeff_temp = -mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[1] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[1];
	    else if(flag[1] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if(flag[1] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[1];
	    else {
		rhs += coeff_temp*(U2_nb[1]+U2_nb[4]+U2_nb[15]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[15]*3+2, -coeff_temp/8.0);
	    }

	    
	    //corner (i+1/2,j,k+1/2)

	    coeff_temp = mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[1] == 1 && flag[5] == 0)
		rhs += coeff_temp*U2_nb[1];
	    else if(flag[1] == 0 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else if(flag[1] == 1 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[1];
	    else {
		rhs += coeff_temp*(U2_nb[1]+U2_nb[5]+U2_nb[16]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[5]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[16]*3+2, -coeff_temp/8.0);
	    }

    
	    //corner (i-1/2,j,k+1/2)

	    coeff_temp = -mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[5] == 1 && flag[0] == 0)
		rhs += coeff_temp*U2_nb[5];
	    else if(flag[5] == 0 && flag[0] == 1)
		rhs += coeff_temp*U2_nb[0];
	    else if(flag[5] == 1 && flag[0] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else {
		rhs += coeff_temp*(U2_nb[5]+U2_nb[0]+U2_nb[17]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[5]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[17]*3+2, -coeff_temp/8.0);
	    }

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[0];

	    solver.Add_b(I*3, rhs);

	    /************************************************************************/


	    //Setting the coeffecients for U1 in the second equation


	    coeff1[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff1[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff1[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff1[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff1[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff1[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    solver.Add_A(I*3+1,I*3+1,1.0);
	    rhs = U1_center;

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Add_A(I*3+1,I_nb[nb]*3+1,-coeff1[nb]);
		    rhs += coeff1[nb]*U1_nb[nb];
		}
		else
		{
		    coeff1[nb] = 2.0*coeff1[nb];
		    rhs += 2.0*coeff1[nb]*U1_nb[nb];
		}
		solver.Add_A(I*3+1,I*3+1,coeff1[nb]);
		rhs -= coeff1[nb]*U1_center;
	    }

	    //set the coefficients for U0 in the second equation (mu*u_y)_x
	    //traverse the four corners
	    //(i-1/2,j-1/2,k),(i+1/2,j-1/2,k),(i+1/2,j+1/2,k),(i-1/2,j+1/2,k)

	    //corner (i-1/2,j-1/2,k)
	    
	    coeff_temp = mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[0] == 1 && flag[2] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[2];
	    else if (flag[0] == 1 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[2]+U0_nb[6]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[6]*3,  -coeff_temp/8.0);
	    }

	    //corner (i+1/2,j-1/2,k)
	    
	    coeff_temp = -mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[1] == 1 && flag[2] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[2];
	    else if (flag[1] == 1 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[2]+U0_nb[7]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[7]*3,  -coeff_temp/8.0);
	    }
	    //corner (i+1/2,j+1/2,k)
	    
	    coeff_temp = mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[1] == 1 && flag[3] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[3];
	    else if (flag[1] == 1 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[3]+U0_nb[8]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[8]*3,  -coeff_temp/8.0);
	    }
	    //corner (i-1/2,j+1/2,k)
	    
	    coeff_temp = -mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[0] == 1 && flag[3] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[3];
	    else if (flag[0] == 1 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[3]+U0_nb[9]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[9]*3,  -coeff_temp/8.0);
	    }

	    //set the coefficients for U2 in the second equation (mu*w_y)_z
	    //traverse the four corners
	    //(i,j-1/2,k-1/2),(i,j+1/2,k-1/2),(i,j+1/2,k+1/2),(i,j-1/2,k+1/2)

	    //corner (i,j-1/2,k-1/2)
	    
	    coeff_temp = mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[2] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[2];
	    else if (flag[2] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if (flag[2] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[2];
	    else {
		rhs += coeff_temp*(U2_nb[2]+U2_nb[4]+U2_nb[10]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[4]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[10]*3+2,  -coeff_temp/8.0);
	    }
	
	    //corner (i,j+1/2,k-1/2)
	    
	    coeff_temp = -mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[3] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[3];
	    else if (flag[3] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if (flag[3] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[3];
	    else {
		rhs += coeff_temp*(U2_nb[3]+U2_nb[4]+U2_nb[11]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[4]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[11]*3+2,  -coeff_temp/8.0);
	    }

	    
	    //corner (i,j+1/2,k+1/2)
	    
	    coeff_temp = mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[3] == 1 && flag[5] == 0)
		rhs += coeff_temp*U2_nb[3];
	    else if (flag[3] == 0 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else if (flag[3] == 1 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[3];
	    else {
		rhs += coeff_temp*(U2_nb[3]+U2_nb[5]+U2_nb[12]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[5]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[12]*3+2,  -coeff_temp/8.0);
	    }

	    //corner (i,j-1/2,k+1/2)
	    
	    coeff_temp = -mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[2] == 1 && flag[5] == 0)
		rhs += coeff_temp*U2_nb[2];
	    else if (flag[2] == 0 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else if (flag[2] == 1 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[2];
	    else {
		rhs += coeff_temp*(U2_nb[2]+U2_nb[5]+U2_nb[13]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[5]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[13]*3+2,  -coeff_temp/8.0);
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[1];

	    solver.Add_b(I*3+1, rhs);

	    /************************************************************************/

	    //Setting the coeffecients of U2 for the third equation


	    coeff2[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff2[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff2[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff2[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff2[4] = m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff2[5] = m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    solver.Add_A(I*3+2,I*3+2,1.0);
	    rhs = U2_center;

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Add_A(I*3+2,I_nb[nb]*3+2,-coeff2[nb]);
		    rhs += coeff2[nb]*U2_nb[nb];
		}
		else
		{
		    coeff2[nb] = 2.0*coeff2[nb];
		    rhs += 2.0*coeff2[nb]*U2_nb[nb];
		}
		solver.Add_A(I*3+2,I*3+2,coeff2[nb]);
		rhs -= coeff2[nb]*U2_center;
	    }

	    //set the coefficients for U0 in the thrid equation (mu*u_z)_x
	    //traverse the four corners
	    //(i-1/2,j,k-1/2),(i+1/2,j,k-1/2),(i+1/2,j,k+1/2),(i-1/2,j,k+1/2)

	    //corner (i-1/2,j,k-1/2)

	    coeff_temp = mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[0] == 1 && flag[4] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[4];
	    else if (flag[0] == 1 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[4]+U0_nb[14]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[14]*3, -coeff_temp/8.0);
	    }

	    //corner (i+1/2,j,k-1/2)

	    coeff_temp = -mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[1] == 1 && flag[4] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[4];
	    else if (flag[1] == 1 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[4]+U0_nb[15]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[15]*3, -coeff_temp/8.0);
	    }


	    //corner (i+1/2,j,k+1/2)

	    coeff_temp = mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[1] == 1 && flag[5] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[5];
	    else if (flag[1] == 1 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[5]+U0_nb[16]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[16]*3, -coeff_temp/8.0);
	    }

	    //corner (i-1/2,j,k+1/2)

	    coeff_temp = -mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[0] == 1 && flag[5] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[5];
	    else if (flag[0] == 1 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[5]+U0_nb[17]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[17]*3, -coeff_temp/8.0);
	    }

	    //set the coefficients for U1 in the thrid equation (mu*v_z)_y
	    //traverse the four corners
	    //(i,j-1/2,k-1/2),(i,j+1/2,k-1/2),(i,j+1/2,k+1/2),(i,j-1/2,k+1/2)

	    //corner (i,j-1/2,k-1/2)

	    coeff_temp = mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[2] == 1 && flag[4] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if (flag[2] == 0 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[4];
	    else if (flag[2] == 1 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[4]+U1_nb[10]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[2]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[10]*3+1, -coeff_temp/8.0);
	    }

	    //corner (i,j+1/2,k-1/2)

	    coeff_temp = -mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[3] == 1 && flag[4] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if (flag[3] == 0 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[4];
	    else if (flag[3] == 1 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[4]+U1_nb[11]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[3]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[11]*3+1, -coeff_temp/8.0);
	    }

	    //corner (i,j+1/2,k+1/2)

	    coeff_temp = mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[3] == 1 && flag[5] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if (flag[3] == 0 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[5];
	    else if (flag[3] == 1 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[5]+U1_nb[12]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[3]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[12]*3+1, -coeff_temp/8.0);
	    }

	    //corner (i,j-1/2,k+1/2)

	    coeff_temp = -mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[2] == 1 && flag[5] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if (flag[2] == 0 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[5];
	    else if (flag[2] == 1 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[5]+U1_nb[13]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[2]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[13]*3+1, -coeff_temp/8.0);
	    }

	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[2];

	    solver.Add_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	//if (rel_residual > 1)
	//{
	//    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	//    solver.Reset_x();
	//    solver.Solve_GMRES();
	//    solver.GetNumIterations(&num_iter);
	//    solver.GetFinalRelativeResidualNorm(&rel_residual);
	//}

	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_2nd_coupled: "
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
}       /* end compDiffWithSmoothProperty3d_2nd_coupled */


void Incompress_Solver_Smooth_3D_Cartesian::computeProjection(void)
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
	int num_nb;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	boolean use_neumann_solver = YES;
	max_value = 0.0;
	double value;
	double sum_div;
	double sum_phi;
	sum_div = 0.0;
	sum_phi = 0.0;
	PetscInt num_iter = 0;
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
		sum_div = sum_div + cell_center[index].m_state.div_U;
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
	    num_nb = 0;
	    for (l = 0; l < 6; ++l)
	    {
		if (I_nb[l] == -1)
		    index_nb[l] = index;
		else num_nb++;
		rho[l] = 1.0/2*(rho0 + cell_center[index_nb[l]].m_state.m_rho);
		coeff[l] = 1.0/rho[l]/(top_h[l/2]*top_h[l/2]);
	    }

	    rhs = cell_center[index].m_state.div_U/accum_dt;

	    aII = 0.0;
	    for (l = 0; l < 6; ++l)
	    {
		if (num_nb <= 1) break;
	    	if (I_nb[l] != -1)
		{
		    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
		}
		else
                {
                    if (isDirichletPresetBdry(front,icoords,dir[l],comp))
                    {
                        /* The boundary condition of phi at preset
                           Dirichlet boundary is set to zero, it then
                           uses non-Neumann poisson solver
                        */
                        use_neumann_solver = NO;
                        aII += -coeff[l];
                    }
                }
	    }
            if (num_nb > 1)
	    {
                solver.Set_A(I,I,aII);
	    }
            else
            {
		printf("\nUnsolvable psudo-pressure, using original pressure!\n");
                solver.Set_A(I,I,1.0);
                rhs = cell_center[index].m_state.m_P;
            }
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Before Petsc Solver in Projection step");
	if (use_neumann_solver)
	{
	    printf("\nUsing Neumann Solver!\n");
	    solver.Solve_withPureNeumann();
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

	}
	else
	{
	    printf("\nUsing non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if(rel_residual > 1)
	    {
		printf("\n The solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }
	}
	stop_clock("After Petsc Solver in Projection step");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("Incompress_Solver_Smooth_3D_Cartesian::"
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
		sum_phi = sum_phi + cell_center[index].m_state.m_phi*top_h[0]*top_h[1]*top_h[2];
		if (value > max_value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_phi,1);
	    printf("\nThe integration of phi is %.16g\n", sum_phi);
	    pp_global_max(&max_value,1);
	    printf("\nThe max value of phi is %.16g\n",max_value);
	}

	FT_FreeThese(1,x);
}	/* end computeProjection3d */

void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocity(void)
{
	int i, j, k, l, index;
	double grad_phi[3], rho;
	COMPONENT comp;
	double speed;
	int icoords[MAXD];

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

void Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm_Adv(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
}


void Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
}
void Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm(double *coords, double t, L_STATE &state) 
{
	computeSourceTerm(coords, state);
}

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void Incompress_Solver_Smooth_3D_Cartesian::solve(double dt)
{

        printf("\nEntering Incompress Solve 3D!! The dt for this step is %.16g\n",dt);

	m_t_old = front->time;
	m_t_int = front->time + front->dt/2.0;
	m_t_new = front->time + front->dt;

	static boolean first = YES;
	if (first)
	{
	    accum_dt = 0.0;
	    first = NO;
	}
	m_dt = dt;
	max_speed = 0.0;

	start_clock("solve");
	setDomain();

	setComponent();
	if (debugging("trace"))
	    printf("Passed setComponent()\n");
	setGlobalIndex();
	if (debugging("trace"))
	    printf("Passed setGlobalIndex()\n");
	start_clock("setSmoothedProperties");
	setSmoothedProperties();
	stop_clock("setSmoothedProperties");
	if (debugging("trace"))
	    printf("Passed setSmoothedProperties()\n");
	
	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	computeAdvection();
	stop_clock("computeAdvection");
	if (debugging("trace"))
	    printf("max_speed after computeAdvection(): %20.14f\n",
				max_speed);
	if (debugging("sample_velocity"))
	    sampleVelocity();
	
	start_clock("compDiffWithSmoothProperty");
	compDiffWithSmoothProperty_1st_decoupled();
	//compDiffWithSmoothProperty();
	stop_clock("compDiffWithSmoothProperty");
	if (debugging("sample_velocity"))
	    sampleVelocity();

        start_clock("compSGS");
        //compSGS();	//Subgrid model by Hyunkyun Lim
        stop_clock("compSGS");

	if (debugging("trace"))
	    printf("max_speed after compDiffWithSmoothProperty(): %20.14f\n",
				max_speed);

	// 2) projection step
	accum_dt += m_dt;
	if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    computeProjection();
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
	    sampleVelocity();

	if (debugging("trace"))
	    printf("max_speed after computeNewVelocity(): %20.14f\n",
				max_speed);

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
	stop_clock("solve");
}	/* end solve */


double Incompress_Solver_Smooth_3D_Cartesian::getVorticityX(int i, int j, int k)
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

double Incompress_Solver_Smooth_3D_Cartesian::getVorticityY(int i, int j, int k)
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

double Incompress_Solver_Smooth_3D_Cartesian::getVorticityZ(int i, int j, int k)
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


void Incompress_Solver_Smooth_3D_Cartesian::copyMeshStates()
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


void Incompress_Solver_Smooth_3D_Cartesian::
	compDiffWithSmoothProperty_1st_decoupled(void)
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
	PetscInt num_iter;
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
                                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
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

            	coeff[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = 0.5*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = 0.5*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, state);

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

	    //if(rel_residual > 1)
	    //{
	    //	printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	    //	solver.Reset_x();
	    //	solver.Solve_GMRES();
	    //	solver.GetNumIterations(&num_iter);
            //	solver.GetFinalRelativeResidualNorm(&rel_residual);
	    //}
	    stop_clock("After Petsc solve");

            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
			"compDiffWithSmoothProperty_1st_decoupled: "
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
}       /* end compDiffWithSmoothProperty3d_1st_decoupled */


void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmI(void)
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

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmII(void)
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
	for(k = 0; k <= top_gmax[2]; k++)
	for(j = 0; j <= top_gmax[1]; j++)
	for(i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}
}        /* end computePressurePmII3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmIII(void)
{
        int i,j,k,index;
        double mu0;

	for(k = 0; k <= top_gmax[2]; k++)
	for(j = 0; j <= top_gmax[1]; j++)
	for(i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    mu0 = 0.5*cell_center[index].m_state.m_mu;
	    cell_center[index].m_state.m_P = 
		cell_center[index].m_state.m_phi -
		accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = 0.0;
	}
}        /* end computePressurePmIII3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressure(void)
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

void Incompress_Solver_Smooth_3D_Cartesian::computeGradientQ(void)
{
	int i,j,k,l,index;
	double *grad_q;
	int icoords[MAXD];

	for (k = kmin; k <= kmax; ++k)
	for (j = jmin; j <= jmax; ++j)
	for (i = imin; i <= imax; ++i)
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

	    computeFieldPointGrad(icoords,array,grad_q);
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
void Incompress_Solver_Smooth_3D_Cartesian::surfaceTension_Fedkiw(
	double delta,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
	double sigma,
	double rho,
	double *weight)
{
    	/*
    	double cellVolume = top_h[0]*top_h[1]*top_h[2];
	int i,j,k,num_tris;
	TRI *tri,*tri_list[MAX_TRI_FOR_INTEGRAL];
	double kappa_tmp,kappa,mag_nor,area,delta;
	double median[MAXD],nor[MAXD];
	POINT *p;

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
	    mag_nor = mag_vector(nor,3);
	    area = 0.5*mag_nor;
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
	    force[j] /= cellVolume;
	*/
}	/* end surfaceTension3d */

void Incompress_Solver_Smooth_3D_Cartesian::setInitialCondition()
{
	int i;
	COMPONENT comp;
	double coords[MAXD];

	FT_MakeGridIntfc(front);
	setDomain();

        m_rho[0] = iFparams->rho1;
        m_rho[1] = iFparams->rho2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	m_sigma = iFparams->surf_tension;
	mu_min = rho_min = HUGE;
	for (i = 0; i < 2; ++i)
	{
	    if (ifluid_comp(m_comp[i]))
	    {
        	mu_min = std::min(mu_min,m_mu[i]);
        	rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

	// Initialize state at cell_center
        for (i = 0; i < cell_center.size(); i++)
        {
            getRectangleCenter(i, coords);
	    cell_center[i].m_state.setZero();
	    comp = top_comp[i];
	    if (getInitialState != NULL)
	    	(*getInitialState)(comp,coords,cell_center[i].m_state,dim,
						iFparams);
        }
	computeGradientQ();
        copyMeshStates();
	setAdvectionDt();
}       /* end setInitialCondition */

double Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointDiv(
        int *icoords,
        double **field)
{
        int index;
        COMPONENT comp;
        int i, j,k,nb;
	int index_nb[6];
        double div,u_nb[2],v_nb[2],w_nb[2];
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



	if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    u_nb[0] = getStateXvel(intfc_state);
	else
	    u_nb[0] = (field[0][index] + field[0][index_nb[0]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,EAST,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    u_nb[1] = getStateXvel(intfc_state);
	else
	    u_nb[1] = (field[0][index] + field[0][index_nb[1]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,SOUTH,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    v_nb[0] = getStateYvel(intfc_state);
	else
	    v_nb[0] = (field[1][index] + field[1][index_nb[2]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,NORTH,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    v_nb[1] = getStateYvel(intfc_state);
	else
	    v_nb[1] = (field[1][index] + field[1][index_nb[3]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    w_nb[0] = getStateZvel(intfc_state);
	else
	    w_nb[0] = (field[2][index] + field[2][index_nb[4]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    w_nb[1] = getStateZvel(intfc_state);
	else
	    w_nb[1] = (field[2][index] + field[2][index_nb[5]])/2.0;

        div = (u_nb[1] - u_nb[0])/top_h[0] + (v_nb[1] - v_nb[0])/top_h[1] + (w_nb[1] - w_nb[0])/top_h[2];
        return div;
}       /* end computeFieldPointDiv */

void Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointGradPhi(
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

	for (nb = 0; nb < 6; nb++)
	{
	    if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state(hs) == NULL)
		    p_nbedge[nb] = 0.0;
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
		p_nbedge[nb] = (p0 + field[index_nb[nb]])/2.0;
	}
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/top_h[0];
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
	grad_field[2] = (p_nbedge[5] - p_nbedge[4])/top_h[2];
}


void Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointGrad(
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

	for (nb = 0; nb < 6; nb++)
	{
	    if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state(hs) == NULL)
		    p_nbedge[nb] = 0.0;
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
		p_nbedge[nb] = (p0 + field[index_nb[nb]])/2.0;
	}
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/top_h[0];
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
	grad_field[2] = (p_nbedge[5] - p_nbedge[4])/top_h[2];
}

void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_2nd_decoupled(void)
{

    printf("\nEntering the Incompress_debug and the m_dt for diffusion step is %.16g\n",m_dt);

    COMPONENT comp;
    int index,index_nb[6],size;
    int I,I_nb[6];
    int i,j,k,l,nb,icoords[MAXD];
    L_STATE source_term;
    INTERFACE *intfc = front->interf;
    double coords[MAXD], crx_coords[MAXD];
    double coeff[6],mu[6],mu0,rho,corner[6],rhs;

    // U_nb contains states at neighbor cell or states on the boundary.
    double U_nb[6], U_nb_new[6], U_center;

    double speed;
    double *x;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    double (*getStateVel[3])(POINTER) =
    {getStateXvel,getStateYvel,getStateZvel};
    POINTER intfc_state;
    HYPER_SURF *hs;

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
		    U_center =  cell_center[index].m_state.m_U[l];

		    for (nb = 0; nb < 6; nb++)
		    {
			if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
				comp,&intfc_state,&hs,crx_coords,m_t_old) &&
				wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
			{
			    // old boundary condition
			    U_nb[nb] = getStateVel[l](intfc_state);
			    // new boundary condition
			    FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
				    comp,&intfc_state,&hs,crx_coords,m_t_new);
			    U_nb_new[nb] = getStateVel[l](intfc_state);

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


		    getRectangleCenter(index, coords);
		    computeSourceTerm(coords, source_term);


		    rhs = 0;
		    double dh[3] = {top_h[0], top_h[1], top_h[2]};
		    for(nb = 0; nb<6; nb++)
		    {
			// use dh[1]*dh[2] as the face area
			if(nb<2)
			{
			    dh[0] = top_h[0];
			    dh[1] = top_h[1];
			    dh[2] = top_h[2];
			}
			else if(nb<4 && nb>=2)
			{
			    dh[0] = top_h[1];
			    dh[1] = top_h[2];
			    dh[2] = top_h[0];
			}
			else if(nb<6 && nb>=4)
			{
			    dh[0] = top_h[2];
			    dh[1] = top_h[0];
			    dh[2] = top_h[1];
			}

			if(I_nb[nb]>=0)
			{
			    // u^{*}
			    solver.Add_A(I, I,        (-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
			    solver.Add_A(I, I_nb[nb], (-1) *  0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
			    // u^{n}
			    rhs += 0.5*m_dt/rho*mu[nb] *
				    (U_nb[nb]-U_center) /(dh[0]*dh[0]);
			}
			else		// boundary
			{
			    // u^{*}
			    solver.Add_A(I, I,       (-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2);
			    rhs += 0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2 * U_nb_new[nb];
			    // u^{n}
			    rhs += 0.5*m_dt/rho*mu[nb] *
				    (U_nb[nb]-U_center) /(dh[0]*dh[0]) * 2;
			}
		    }

		    // interior
		    rhs += m_dt*source_term.m_U[l];
		    rhs += m_dt*cell_center[index].m_state.f_surf[l];
		    rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
		    //printf("gradq = %.16g\n",cell_center[index].m_state.grad_q[l]);

		    rhs -= m_dt * cell_center[index].m_state.m_adv[l]; //advection source term

		    solver.Add_A(I, I, 1.0);
		    rhs += U_center;

		    solver.Add_b(I, rhs);
		}

	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);
	solver.Solve_GMRES();

	// get back the solution
	solver.Get_x(x);

	PetscInt num_iter;
	double rel_residual;
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	if (debugging("PETSc"))
	    (void) printf("L_CARTESIAN::"
		    "compDiffWithSmoothProperty_2nd_decoupled: "
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

}

void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_decoupled(void)
{
    int index;
    int I;
    int i,j,k,icoords[MAXD];
 
    setIndexMap();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I  = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index  = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    double convectionTerm[3];
	    getAdvectionTerm_decoupled(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	    cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}
	
}

void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_coupled(void)
{
    int index;
    int I;
    int i,j,k,icoords[MAXD];
 
    setIndexMap();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I  = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index  = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    double convectionTerm[3];
	    getAdvectionTerm_coupled(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	    cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}
	
}


void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_decoupled_upgraded(void)
{
    int index;
    int I;
    int i,j,k,icoords[MAXD];
 
    setIndexMap();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I  = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index  = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    double convectionTerm[3];
	    getAdvectionTerm_decoupled_upgraded(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	    cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}
	
}

void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_coupled_upgraded(void)
{
    int index;
    int I;
    int i,j,k,icoords[MAXD];
 
    setIndexMap();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I  = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index  = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    double convectionTerm[3];
	    getAdvectionTerm_coupled_upgraded(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	    cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}
	
}


void Incompress_Solver_Smooth_3D_Cartesian::getAdvectionTerm_decoupled(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west, state_east, state_south, state_north, state_lower, state_upper;

    double dx = top_h[0];
    double dy = top_h[1];
    double dz = top_h[2];

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
    }
    getFaceVelocity_middleStep(icoords,WEST,sr);
    getRiemannSolution(COORD_X,sl,sr,state_west);

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,WEST,sr);
    }
    getFaceVelocity_middleStep(icoords,EAST,sl);
    getRiemannSolution(COORD_X,sl,sr,state_east);

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,NORTH,sl);
    }
    getFaceVelocity_middleStep(icoords,SOUTH,sr);
    getRiemannSolution(COORD_Y,sl,sr,state_south);

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,SOUTH,sr);
    }
    getFaceVelocity_middleStep(icoords,NORTH,sl);
    getRiemannSolution(COORD_Y,sl,sr,state_north);

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep(ICoords,UPPER,sl);
    }
    getFaceVelocity_middleStep(icoords,LOWER,sr);
    getRiemannSolution(COORD_Z,sl,sr,state_lower);

    // UPPER
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep(ICoords,LOWER,sr);
    }
    getFaceVelocity_middleStep(icoords,UPPER,sl);
    getRiemannSolution(COORD_Z,sl,sr,state_upper);

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[0]-state_west.m_U[0])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[0]-state_south.m_U[0])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[0]-state_lower.m_U[0])/dz;
    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[1]-state_west.m_U[1])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[1]-state_south.m_U[1])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[1]-state_lower.m_U[1])/dz;
    convectionTerm[2] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[2]-state_west.m_U[2])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[2]-state_south.m_U[2])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[2]-state_lower.m_U[2])/dz;
}


void Incompress_Solver_Smooth_3D_Cartesian::getAdvectionTerm_coupled(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west, state_east, state_south, state_north, state_lower, state_upper;

    double dx = top_h[0];
    double dy = top_h[1];
    double dz = top_h[2];

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
    }
    getFaceVelocity_middleStep_coupled(icoords,WEST,sr);
    getRiemannSolution(COORD_X,sl,sr,state_west);

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,WEST,sr);
    }
    getFaceVelocity_middleStep_coupled(icoords,EAST,sl);
    getRiemannSolution(COORD_X,sl,sr,state_east);

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,NORTH,sl);
    }
    getFaceVelocity_middleStep_coupled(icoords,SOUTH,sr);
    getRiemannSolution(COORD_Y,sl,sr,state_south);

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,SOUTH,sr);
    }
    getFaceVelocity_middleStep_coupled(icoords,NORTH,sl);
    getRiemannSolution(COORD_Y,sl,sr,state_north);

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_coupled(ICoords,UPPER,sl);
    }
    getFaceVelocity_middleStep_coupled(icoords,LOWER,sr);
    getRiemannSolution(COORD_Z,sl,sr,state_lower);

    // UPPER
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_coupled(ICoords,LOWER,sr);
    }
    getFaceVelocity_middleStep_coupled(icoords,UPPER,sl);
    getRiemannSolution(COORD_Z,sl,sr,state_upper);

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[0]-state_west.m_U[0])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[0]-state_south.m_U[0])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[0]-state_lower.m_U[0])/dz;
    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[1]-state_west.m_U[1])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[1]-state_south.m_U[1])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[1]-state_lower.m_U[1])/dz;
    convectionTerm[2] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[2]-state_west.m_U[2])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[2]-state_south.m_U[2])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[2]-state_lower.m_U[2])/dz;
}


void Incompress_Solver_Smooth_3D_Cartesian::getAdvectionTerm_decoupled_upgraded(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west_bar, state_east_bar, state_south_bar, state_north_bar, state_lower_bar, state_upper_bar;

    L_STATE state_west_hat, state_east_hat, state_south_hat, state_north_hat, state_lower_hat, state_upper_hat;

    L_STATE state_west_hat_l, state_west_hat_r;
    L_STATE state_east_hat_l, state_east_hat_r;
    L_STATE state_south_hat_l, state_south_hat_r;
    L_STATE state_north_hat_l, state_north_hat_r;
    L_STATE state_lower_hat_l, state_lower_hat_r;
    L_STATE state_upper_hat_l, state_upper_hat_r;

    double transverseD[3];

    double dx = top_h[0];
    double dy = top_h[1];
    double dz = top_h[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;


    ///////////////////////////// Get the state_hat on six faces first //////////////////////

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,state_west_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_west_hat = state_west_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,WEST,state_west_hat_r);
	//getRiemannSolution(COORD_X,state_west_hat_l,state_west_hat_r,state_west_hat);
    }
    else
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,EAST,state_west_hat_l);
	getFaceVelocity_middleStep_hat(icoords,WEST,state_west_hat_r);
	getRiemannSolution(COORD_X,state_west_hat_l,state_west_hat_r,state_west_hat);
    }

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,state_east_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_east_hat = state_east_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,EAST,state_east_hat_l);
	//getRiemannSolution(COORD_X,state_east_hat_l,state_east_hat_r,state_east_hat);
    }
    else
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,WEST,state_east_hat_r);
	getFaceVelocity_middleStep_hat(icoords,EAST,state_east_hat_l);
	getRiemannSolution(COORD_X,state_east_hat_l,state_east_hat_r,state_east_hat);
    }

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,state_south_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_south_hat = state_south_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	//getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,NORTH,state_south_hat_l);
	getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,state_north_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_north_hat = state_north_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	//getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,SOUTH,state_north_hat_r);
	getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);

    }

    // LOWER 
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,state_lower_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_lower_hat = state_lower_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	//getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_hat(ICoords,UPPER,state_lower_hat_l);
	getFaceVelocity_middleStep_hat(icoords,LOWER,state_lower_hat_r);
	getRiemannSolution(COORD_Z,state_lower_hat_l,state_lower_hat_r,state_lower_hat);
    }

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,state_upper_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_upper_hat = state_upper_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	//getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_hat(ICoords,LOWER,state_upper_hat_r);
	getFaceVelocity_middleStep_hat(icoords,UPPER,state_upper_hat_l);
	getRiemannSolution(COORD_Z,state_upper_hat_l,state_upper_hat_r,state_upper_hat);

    }

    ///////////////////////////// get the state_bar on four edges ////////////////////////////////

    // WEST

    transverseD[0] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_west_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	//getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }
    else
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_bar(ICoords,EAST,sl,transverseD,state_west_hat_l);
	getFaceVelocity_middleStep_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }

    // EAST
    

    transverseD[0] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_east_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	//getRiemannSolution(COORD_X,sl,sr,state_east_bar);

    }
    else
    {

	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_bar(ICoords,WEST,sr,transverseD,state_east_hat_r);
	getFaceVelocity_middleStep_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	getRiemannSolution(COORD_X,sl,sr,state_east_bar);
    }

    // SOUTH
    //

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_south_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	//getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_bar(ICoords,NORTH,sl,transverseD,state_south_hat_l);
	getFaceVelocity_middleStep_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }

    // NORTH

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_north_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_bar(ICoords,SOUTH,sr,transverseD,state_north_hat_r);
	getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	getRiemannSolution(COORD_Y,sl,sr,state_north_bar);

    }

    // LOWER 

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);

    
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_lower_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_bar(ICoords,UPPER,sl,transverseD,state_lower_hat_l);
	getFaceVelocity_middleStep_bar(icoords,LOWER,sr,transverseD,state_lower_hat_r);
	getRiemannSolution(COORD_Z,sl,sr,state_lower_bar);

    }

    // UPPER 

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);

    
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_upper_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_bar(ICoords,LOWER,sr,transverseD,state_upper_hat_r);
	getFaceVelocity_middleStep_bar(icoords,UPPER,sl,transverseD,state_upper_hat_l);
	getRiemannSolution(COORD_Z,sl,sr,state_upper_bar);

    }

    convectionTerm[0] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[0]-state_west_bar.m_U[0])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[0]-state_south_bar.m_U[0])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[0]-state_lower_bar.m_U[0])/dz;
    convectionTerm[1] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[1]-state_west_bar.m_U[1])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[1]-state_south_bar.m_U[1])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[1]-state_lower_bar.m_U[1])/dz;
    convectionTerm[2] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[2]-state_west_bar.m_U[2])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[2]-state_south_bar.m_U[2])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[2]-state_lower_bar.m_U[2])/dz;
}


void Incompress_Solver_Smooth_3D_Cartesian::getAdvectionTerm_coupled_upgraded(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west_bar, state_east_bar, state_south_bar, state_north_bar, state_lower_bar, state_upper_bar;

    L_STATE state_west_hat, state_east_hat, state_south_hat, state_north_hat, state_lower_hat, state_upper_hat;

    L_STATE state_west_hat_l, state_west_hat_r;
    L_STATE state_east_hat_l, state_east_hat_r;
    L_STATE state_south_hat_l, state_south_hat_r;
    L_STATE state_north_hat_l, state_north_hat_r;
    L_STATE state_lower_hat_l, state_lower_hat_r;
    L_STATE state_upper_hat_l, state_upper_hat_r;

    double transverseD[3];

    double dx = top_h[0];
    double dy = top_h[1];
    double dz = top_h[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;


    ///////////////////////////// Get the state_hat on six faces first //////////////////////

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,state_west_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_west_hat = state_west_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,WEST,state_west_hat_r);
	//getRiemannSolution(COORD_X,state_west_hat_l,state_west_hat_r,state_west_hat);
    }
    else
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,EAST,state_west_hat_l);
	getFaceVelocity_middleStep_hat(icoords,WEST,state_west_hat_r);
	getRiemannSolution(COORD_X,state_west_hat_l,state_west_hat_r,state_west_hat);
    }

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,state_east_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_east_hat = state_east_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,EAST,state_east_hat_l);
	//getRiemannSolution(COORD_X,state_east_hat_l,state_east_hat_r,state_east_hat);
    }
    else
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,WEST,state_east_hat_r);
	getFaceVelocity_middleStep_hat(icoords,EAST,state_east_hat_l);
	getRiemannSolution(COORD_X,state_east_hat_l,state_east_hat_r,state_east_hat);
    }

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,state_south_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_south_hat = state_south_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	//getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,NORTH,state_south_hat_l);
	getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,state_north_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_north_hat = state_north_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	//getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,SOUTH,state_north_hat_r);
	getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);

    }

    // LOWER 
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,state_lower_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_lower_hat = state_lower_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	//getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_hat(ICoords,UPPER,state_lower_hat_l);
	getFaceVelocity_middleStep_hat(icoords,LOWER,state_lower_hat_r);
	getRiemannSolution(COORD_Z,state_lower_hat_l,state_lower_hat_r,state_lower_hat);
    }

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,state_upper_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_upper_hat = state_upper_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	//getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_hat(ICoords,LOWER,state_upper_hat_r);
	getFaceVelocity_middleStep_hat(icoords,UPPER,state_upper_hat_l);
	getRiemannSolution(COORD_Z,state_upper_hat_l,state_upper_hat_r,state_upper_hat);

    }

    ///////////////////////////// get the state_bar on four edges ////////////////////////////////

    // WEST

    transverseD[0] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_west_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	//getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }
    else
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled_bar(ICoords,EAST,sl,transverseD,state_west_hat_l);
	getFaceVelocity_middleStep_coupled_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }

    // EAST
    

    transverseD[0] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_east_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	//getRiemannSolution(COORD_X,sl,sr,state_east_bar);

    }
    else
    {

	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled_bar(ICoords,WEST,sr,transverseD,state_east_hat_r);
	getFaceVelocity_middleStep_coupled_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	getRiemannSolution(COORD_X,sl,sr,state_east_bar);
    }

    // SOUTH
    //

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_south_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	//getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled_bar(ICoords,NORTH,sl,transverseD,state_south_hat_l);
	getFaceVelocity_middleStep_coupled_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }

    // NORTH

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_north_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled_bar(ICoords,SOUTH,sr,transverseD,state_north_hat_r);
	getFaceVelocity_middleStep_coupled_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	getRiemannSolution(COORD_Y,sl,sr,state_north_bar);

    }

    // LOWER 

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);

    
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_lower_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_coupled_bar(ICoords,UPPER,sl,transverseD,state_lower_hat_l);
	getFaceVelocity_middleStep_coupled_bar(icoords,LOWER,sr,transverseD,state_lower_hat_r);
	getRiemannSolution(COORD_Z,sl,sr,state_lower_bar);

    }

    // UPPER 

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);

    
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_upper_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_coupled_bar(ICoords,LOWER,sr,transverseD,state_upper_hat_r);
	getFaceVelocity_middleStep_coupled_bar(icoords,UPPER,sl,transverseD,state_upper_hat_l);
	getRiemannSolution(COORD_Z,sl,sr,state_upper_bar);

    }

    convectionTerm[0] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[0]-state_west_bar.m_U[0])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[0]-state_south_bar.m_U[0])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[0]-state_lower_bar.m_U[0])/dz;
    convectionTerm[1] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[1]-state_west_bar.m_U[1])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[1]-state_south_bar.m_U[1])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[1]-state_lower_bar.m_U[1])/dz;
    convectionTerm[2] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[2]-state_west_bar.m_U[2])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[2]-state_south_bar.m_U[2])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[2]-state_lower_bar.m_U[2])/dz;
}



void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_hat(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_hat)
{
    L_STATE  state_orig;
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state_orig = cell_center[index].m_state;
    double sL;

    double dx = 0, dy = 0, dz = 0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};



    switch(dir)
    {
    case WEST:
	getLimitedSlope(icoords,COORD_X,slope_x_limited);
	dx = top_h[0];
	if (state_orig.m_U[0] <= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[2];
	break;


    case EAST:
	getLimitedSlope(icoords,COORD_X,slope_x_limited);
	dx = top_h[0];
	if (state_orig.m_U[0] >= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[2];
	break;


    case SOUTH:
	getLimitedSlope(icoords,COORD_Y,slope_y_limited);
	dy = top_h[1];
	if (state_orig.m_U[1] <= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[2];
	break;

    case NORTH:
	getLimitedSlope(icoords,COORD_Y,slope_y_limited);
	dy = top_h[1];
	if (state_orig.m_U[1] >= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[2];
	break;

    case LOWER:
	getLimitedSlope(icoords,COORD_Z,slope_z_limited);
	dz = top_h[2];
	if (state_orig.m_U[2] <= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[2];
	break;

    case UPPER:
	getLimitedSlope(icoords,COORD_Z,slope_z_limited);
	dz = top_h[2];
	if (state_orig.m_U[2] >= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[2];
	break;
    default:
	assert(false);
    }

}


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_bar(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_bar,
	double transverseD[3],
	L_STATE state_hat)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[3], gradP[3];
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];
    gradP[2] = cell_center[index].m_state.grad_q[2];

    getDifffusion(icoords,diffusion);
    state_bar.m_U[0] = state_hat.m_U[0] + m_dt/2.0 * (-transverseD[0] + diffusion[0]/rho - gradP[0]/rho);
    state_bar.m_U[1] = state_hat.m_U[1] + m_dt/2.0 * (-transverseD[1] + diffusion[1]/rho - gradP[1]/rho);
    state_bar.m_U[2] = state_hat.m_U[2] + m_dt/2.0 * (-transverseD[2] + diffusion[2]/rho - gradP[2]/rho);

    double coords[3];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);

    state_bar.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_bar.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_bar.m_U[2] += m_dt/2 * source_term.m_U[2];
}


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_coupled_bar(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_bar,
	double transverseD[3],
	L_STATE state_hat)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[3], gradP[3];
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];
    gradP[2] = cell_center[index].m_state.grad_q[2];

    getDiffusion_coupled(icoords,diffusion);
    state_bar.m_U[0] = state_hat.m_U[0] + m_dt/2.0 * (-transverseD[0] + diffusion[0]/rho - gradP[0]/rho);
    state_bar.m_U[1] = state_hat.m_U[1] + m_dt/2.0 * (-transverseD[1] + diffusion[1]/rho - gradP[1]/rho);
    state_bar.m_U[2] = state_hat.m_U[2] + m_dt/2.0 * (-transverseD[2] + diffusion[2]/rho - gradP[2]/rho);

    double coords[3];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);

    state_bar.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_bar.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_bar.m_U[2] += m_dt/2 * source_term.m_U[2];
}


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    L_STATE state;
    state = cell_center[index].m_state;
    double rho = cell_center[index].m_state.m_rho;


    double dx = 0, dy = 0,dz = 0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};

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
    double diffusion[3], gradP[3];
    getDifffusion(icoords,diffusion);
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];
    gradP[2] = cell_center[index].m_state.grad_q[2];

    state_face.m_U[0] += m_dt/2 *
	    (diffusion[0]/rho - state.m_U[0]*slope_x_limited[0] - state.m_U[1]*slope_y_limited[0] - state.m_U[2]*slope_z_limited[0] - gradP[0]/rho);
    state_face.m_U[1] += m_dt/2 *
	    (diffusion[1]/rho - state.m_U[0]*slope_x_limited[1] - state.m_U[1]*slope_y_limited[1] - state.m_U[2]*slope_z_limited[1] - gradP[1]/rho);
    state_face.m_U[2] += m_dt/2 *
    	    (diffusion[2]/rho - state.m_U[0]*slope_x_limited[2] - state.m_U[1]*slope_y_limited[2] - state.m_U[2]*slope_z_limited[2] - gradP[2]/rho);

    //    return;
    // rhs
    double coords[3];
    L_STATE source_term;

    getRectangleCenter(d_index3d(icoords[0],icoords[1],icoords[2],top_gmax), coords);
    computeSourceTerm_Adv(coords, source_term);

    state_face.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_face.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_face.m_U[2] += m_dt/2 * source_term.m_U[2];
}

void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_coupled(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    L_STATE state;
    state = cell_center[index].m_state;
    double rho = cell_center[index].m_state.m_rho;


    double dx = 0, dy = 0,dz = 0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};

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
    double diffusion[3], gradP[3];
    getDiffusion_coupled(icoords,diffusion);
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];
    gradP[2] = cell_center[index].m_state.grad_q[2];

    state_face.m_U[0] += m_dt/2 *
	    (diffusion[0]/rho - state.m_U[0]*slope_x_limited[0] - state.m_U[1]*slope_y_limited[0] - state.m_U[2]*slope_z_limited[0] - gradP[0]/rho);
    state_face.m_U[1] += m_dt/2 *
	    (diffusion[1]/rho - state.m_U[0]*slope_x_limited[1] - state.m_U[1]*slope_y_limited[1] - state.m_U[2]*slope_z_limited[1] - gradP[1]/rho);
    state_face.m_U[2] += m_dt/2 *
    	    (diffusion[2]/rho - state.m_U[0]*slope_x_limited[2] - state.m_U[1]*slope_y_limited[2] - state.m_U[2]*slope_z_limited[2] - gradP[2]/rho);

    //    return;
    // rhs
    double coords[3];
    L_STATE source_term;

    getRectangleCenter(d_index3d(icoords[0],icoords[1],icoords[2],top_gmax), coords);
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


void Incompress_Solver_Smooth_3D_Cartesian::getDifffusion(
	int *icoords,
	double diffusion[3])
{
    double Uxx[3], Uyy[3], Uzz[3];
    getDU2(icoords,COORD_X,Uxx);
    getDU2(icoords,COORD_Y,Uyy);
    getDU2(icoords,COORD_Z,Uzz);

    double mu = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state.m_mu;

    diffusion[0] = mu * (Uxx[0] + Uyy[0] + Uzz[0]);
    diffusion[1] = mu * (Uxx[1] + Uyy[1] + Uzz[1]);
    diffusion[2] = mu * (Uxx[2] + Uyy[2] + Uzz[2]);
}


void Incompress_Solver_Smooth_3D_Cartesian::getDiffusion_coupled(
	int *icoords,
	double diffusion[3])
{
    int index,index_nb[18];
    double mu[6],mu_edge[6],mu0;
    L_STATE Unb;
    double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
    // U0[6] -- U0[17] U1[6] -- U1[17]  U2[6] -- U2[17] are corner values
    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    bool bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];

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

    //xy cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	
    //yz cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);
	
    //xz cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    mu0 = cell_center[index].m_state.m_mu;

    U0_center = cell_center[index].m_state.m_U[0];
    U1_center = cell_center[index].m_state.m_U[1];
    U2_center = cell_center[index].m_state.m_U[2];

    for (nb = 0; nb < 6; nb++)
    {
	bNoBoundary[nb] = getNeighborOrBoundaryState(icoords,dir[nb],Unb,m_t_old);
	U0_nb[nb] = Unb.m_U[0];
	U1_nb[nb] = Unb.m_U[1];
	U2_nb[nb] = Unb.m_U[2];
	if(!bNoBoundary[nb])
	{
	    mu[nb] = mu0;
	    mu_edge[nb] = mu0;
	}
	else
	{
	    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
	    mu_edge[nb] = 0.5*(mu[nb] + mu0);
	}
    }

    // non-cross derivative terms

    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    if (bNoBoundary[0])
	dh0[0] = top_h[0];
    else
	dh0[0] = top_h[0]/2.0;
    if (bNoBoundary[1])
	dh1[0] = top_h[0];
    else
	dh1[0] = top_h[0]/2.0;

    if (bNoBoundary[2])
	dh0[1] = top_h[1];
    else
	dh0[1] = top_h[1]/2.0;
    if (bNoBoundary[3])
	dh1[1] = top_h[1];
    else
	dh1[1] = top_h[1]/2.0;

    if (bNoBoundary[4])
	dh0[2] = top_h[2];
    else
	dh0[2] = top_h[2]/2.0;
    if (bNoBoundary[5])
	dh1[2] = top_h[2];
    else
	dh1[2] = top_h[2]/2.0;

    diffusion[0] += 2.0*(mu_edge[1]*(U0_nb[1]-U0_center)/dh1[0] - mu_edge[0]*(U0_center-U0_nb[0])/dh0[0])/dh[0];// (2*mu*u_x)_x
    diffusion[1] +=     (mu_edge[1]*(U1_nb[1]-U1_center)/dh1[0] - mu_edge[0]*(U1_center-U1_nb[0])/dh0[0])/dh[0];// (mu*v_x)_x
    diffusion[2] +=     (mu_edge[1]*(U2_nb[1]-U2_center)/dh1[0] - mu_edge[0]*(U2_center-U2_nb[0])/dh0[0])/dh[0];// (mu*w_x)_x

    diffusion[0] +=     (mu_edge[3]*(U0_nb[3]-U0_center)/dh1[1] - mu_edge[2]*(U0_center-U0_nb[2])/dh0[1])/dh[1];// (mu*u_y)_y
    diffusion[1] += 2.0*(mu_edge[3]*(U1_nb[3]-U1_center)/dh1[1] - mu_edge[2]*(U1_center-U1_nb[2])/dh0[1])/dh[1];// (2*mu*v_y)_y
    diffusion[2] +=     (mu_edge[3]*(U2_nb[3]-U2_center)/dh1[1] - mu_edge[2]*(U2_center-U2_nb[2])/dh0[1])/dh[1];// (mu*w_y)_y
 
    diffusion[0] +=     (mu_edge[5]*(U0_nb[5]-U0_center)/dh1[2] - mu_edge[4]*(U0_center-U0_nb[4])/dh0[2])/dh[2];// (mu*u_z)_z
    diffusion[1] +=     (mu_edge[5]*(U1_nb[5]-U1_center)/dh1[2] - mu_edge[4]*(U1_center-U1_nb[4])/dh0[2])/dh[2];// (mu*v_z)_z
    diffusion[2] += 2.0*(mu_edge[5]*(U2_nb[5]-U2_center)/dh1[2] - mu_edge[4]*(U2_center-U2_nb[4])/dh0[2])/dh[2];// (2*mu*w_z)_z

    //cross derivative terms

    //traverse the corners on 3 cut planes

    //corner (i-1/2,j-1/2,k)

    if (!bNoBoundary[0] && bNoBoundary[2])
    {
	U0_nb[6] = U0_nb[0];
	U1_nb[6] = U1_nb[0];
	U2_nb[6] = U2_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[2])
    {
	U0_nb[6] = U0_nb[2];
	U1_nb[6] = U1_nb[2];
	U2_nb[6] = U2_nb[2];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[2])
    {
	U0_nb[6] = U0_nb[0];
	U1_nb[6] = U1_nb[0];
	U2_nb[6] = U2_nb[0];
    }
    else
    {
	U0_nb[6] = (U0_nb[0]+U0_nb[2]+cell_center[index_nb[6]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[6] = (U1_nb[0]+U1_nb[2]+cell_center[index_nb[6]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[6] = (U2_nb[0]+U2_nb[2]+cell_center[index_nb[6]].m_state.m_U[2]+U2_center)/4.0;
    }

    //corner (i+1/2,j-1/2,k)

    if (!bNoBoundary[1] && bNoBoundary[2])
    {
	U0_nb[7] = U0_nb[1];
	U1_nb[7] = U1_nb[1];
	U2_nb[7] = U2_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[2])
    {
	U0_nb[7] = U0_nb[2];
	U1_nb[7] = U1_nb[2];
	U2_nb[7] = U2_nb[2];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[2])
    {
	U0_nb[7] = U0_nb[1];
	U1_nb[7] = U1_nb[1];
	U2_nb[7] = U2_nb[1];
    }
    else
    {
	U0_nb[7] = (U0_nb[1]+U0_nb[2]+cell_center[index_nb[7]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[7] = (U1_nb[1]+U1_nb[2]+cell_center[index_nb[7]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[7] = (U2_nb[1]+U2_nb[2]+cell_center[index_nb[7]].m_state.m_U[2]+U2_center)/4.0;
    }


    //corner (i+1/2,j+1/2,k)

    if (!bNoBoundary[1] && bNoBoundary[3])
    {
	U0_nb[8] = U0_nb[1];
	U1_nb[8] = U1_nb[1];
	U2_nb[8] = U2_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[3])
    {
	U0_nb[8] = U0_nb[3];
	U1_nb[8] = U1_nb[3];
	U2_nb[8] = U2_nb[3];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[3])
    {
	U0_nb[8] = U0_nb[1];
	U1_nb[8] = U1_nb[1];
	U2_nb[8] = U2_nb[1];
    }
    else
    {
	U0_nb[8] = (U0_nb[1]+U0_nb[3]+cell_center[index_nb[8]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[8] = (U1_nb[1]+U1_nb[3]+cell_center[index_nb[8]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[8] = (U2_nb[1]+U2_nb[3]+cell_center[index_nb[8]].m_state.m_U[2]+U2_center)/4.0;
    }


    //corner (i-1/2,j+1/2,k)

    if (!bNoBoundary[0] && bNoBoundary[3])
    {
	U0_nb[9] = U0_nb[0];
	U1_nb[9] = U1_nb[0];
	U2_nb[9] = U2_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[3])
    {
	U0_nb[9] = U0_nb[3];
	U1_nb[9] = U1_nb[3];
	U2_nb[9] = U2_nb[3];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[3])
    {
	U0_nb[9] = U0_nb[0];
	U1_nb[9] = U1_nb[0];
	U2_nb[9] = U2_nb[0];
    }
    else
    {
	U0_nb[9] = (U0_nb[0]+U0_nb[3]+cell_center[index_nb[9]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[9] = (U1_nb[0]+U1_nb[3]+cell_center[index_nb[9]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[9] = (U2_nb[0]+U2_nb[3]+cell_center[index_nb[9]].m_state.m_U[2]+U2_center)/4.0;
    }


    //corner (i,j-1/2,k-1/2)

    if (!bNoBoundary[2] && bNoBoundary[4])
    {
	U0_nb[10] = U0_nb[2];
	U1_nb[10] = U1_nb[2];
	U2_nb[10] = U2_nb[2];
    }
    else if(bNoBoundary[2] && !bNoBoundary[4])
    {
	U0_nb[10] = U0_nb[4];
	U1_nb[10] = U1_nb[4];
	U2_nb[10] = U2_nb[4];
    }
    else if(!bNoBoundary[2] && !bNoBoundary[4])
    {
	U0_nb[10] = U0_nb[2];
	U1_nb[10] = U1_nb[2];
	U2_nb[10] = U2_nb[2];
    }
    else
    {
	U0_nb[10] = (U0_nb[2]+U0_nb[4]+cell_center[index_nb[10]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[10] = (U1_nb[2]+U1_nb[4]+cell_center[index_nb[10]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[10] = (U2_nb[2]+U2_nb[4]+cell_center[index_nb[10]].m_state.m_U[2]+U2_center)/4.0;
    }


    //corner (i,j+1/2,k-1/2)

    if (!bNoBoundary[3] && bNoBoundary[4])
    {
	U0_nb[11] = U0_nb[3];
	U1_nb[11] = U1_nb[3];
	U2_nb[11] = U2_nb[3];
    }
    else if(bNoBoundary[3] && !bNoBoundary[4])
    {
	U0_nb[11] = U0_nb[4];
	U1_nb[11] = U1_nb[4];
	U2_nb[11] = U2_nb[4];
    }
    else if(!bNoBoundary[3] && !bNoBoundary[4])
    {
	U0_nb[11] = U0_nb[3];
	U1_nb[11] = U1_nb[3];
	U2_nb[11] = U2_nb[3];
    }
    else
    {
	U0_nb[11] = (U0_nb[3]+U0_nb[4]+cell_center[index_nb[11]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[11] = (U1_nb[3]+U1_nb[4]+cell_center[index_nb[11]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[11] = (U2_nb[3]+U2_nb[4]+cell_center[index_nb[11]].m_state.m_U[2]+U2_center)/4.0;
    }

    //corner (i,j+1/2,k+1/2)

    if (!bNoBoundary[3] && bNoBoundary[5])
    {
	U0_nb[12] = U0_nb[3];
	U1_nb[12] = U1_nb[3];
	U2_nb[12] = U2_nb[3];
    }
    else if(bNoBoundary[3] && !bNoBoundary[5])
    {
	U0_nb[12] = U0_nb[5];
	U1_nb[12] = U1_nb[5];
	U2_nb[12] = U2_nb[5];
    }
    else if(!bNoBoundary[3] && !bNoBoundary[5])
    {
	U0_nb[12] = U0_nb[3];
	U1_nb[12] = U1_nb[3];
	U2_nb[12] = U2_nb[3];
    }
    else
    {
	U0_nb[12] = (U0_nb[3]+U0_nb[5]+cell_center[index_nb[12]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[12] = (U1_nb[3]+U1_nb[5]+cell_center[index_nb[12]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[12] = (U2_nb[3]+U2_nb[5]+cell_center[index_nb[12]].m_state.m_U[2]+U2_center)/4.0;
    }


    //corner (i,j-1/2,k+1/2)

    if (!bNoBoundary[2] && bNoBoundary[5])
    {
	U0_nb[13] = U0_nb[2];
	U1_nb[13] = U1_nb[2];
	U2_nb[13] = U2_nb[2];
    }
    else if(bNoBoundary[2] && !bNoBoundary[5])
    {
	U0_nb[13] = U0_nb[5];
	U1_nb[13] = U1_nb[5];
	U2_nb[13] = U2_nb[5];
    }
    else if(!bNoBoundary[2] && !bNoBoundary[5])
    {
	U0_nb[13] = U0_nb[2];
	U1_nb[13] = U1_nb[2];
	U2_nb[13] = U2_nb[2];
    }
    else
    {
	U0_nb[13] = (U0_nb[2]+U0_nb[5]+cell_center[index_nb[13]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[13] = (U1_nb[2]+U1_nb[5]+cell_center[index_nb[13]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[13] = (U2_nb[2]+U2_nb[5]+cell_center[index_nb[13]].m_state.m_U[2]+U2_center)/4.0;
    }


    //corner (i-1/2,j,k-1/2)

    if (!bNoBoundary[0] && bNoBoundary[4])
    {
	U0_nb[14] = U0_nb[0];
	U1_nb[14] = U1_nb[0];
	U2_nb[14] = U2_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[4])
    {
	U0_nb[14] = U0_nb[4];
	U1_nb[14] = U1_nb[4];
	U2_nb[14] = U2_nb[4];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[4])
    {
	U0_nb[14] = U0_nb[0];
	U1_nb[14] = U1_nb[0];
	U2_nb[14] = U2_nb[0];
    }
    else
    {
	U0_nb[14] = (U0_nb[0]+U0_nb[4]+cell_center[index_nb[14]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[14] = (U1_nb[0]+U1_nb[4]+cell_center[index_nb[14]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[14] = (U2_nb[0]+U2_nb[4]+cell_center[index_nb[14]].m_state.m_U[2]+U2_center)/4.0;
    }


    //corner (i+1/2,j,k-1/2)

    if (!bNoBoundary[1] && bNoBoundary[4])
    {
	U0_nb[15] = U0_nb[1];
	U1_nb[15] = U1_nb[1];
	U2_nb[15] = U2_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[4])
    {
	U0_nb[15] = U0_nb[4];
	U1_nb[15] = U1_nb[4];
	U2_nb[15] = U2_nb[4];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[4])
    {
	U0_nb[15] = U0_nb[1];
	U1_nb[15] = U1_nb[1];
	U2_nb[15] = U2_nb[1];
    }
    else
    {
	U0_nb[15] = (U0_nb[1]+U0_nb[4]+cell_center[index_nb[15]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[15] = (U1_nb[1]+U1_nb[4]+cell_center[index_nb[15]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[15] = (U2_nb[1]+U2_nb[4]+cell_center[index_nb[15]].m_state.m_U[2]+U2_center)/4.0;
    }

    //corner (i+1/2,j,k+1/2)

    if (!bNoBoundary[1] && bNoBoundary[5])
    {
	U0_nb[16] = U0_nb[1];
	U1_nb[16] = U1_nb[1];
	U2_nb[16] = U2_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[5])
    {
	U0_nb[16] = U0_nb[5];
	U1_nb[16] = U1_nb[5];
	U2_nb[16] = U2_nb[5];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[5])
    {
	U0_nb[16] = U0_nb[1];
	U1_nb[16] = U1_nb[1];
	U2_nb[16] = U2_nb[1];
    }
    else
    {
	U0_nb[16] = (U0_nb[1]+U0_nb[5]+cell_center[index_nb[16]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[16] = (U1_nb[1]+U1_nb[5]+cell_center[index_nb[16]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[16] = (U2_nb[1]+U2_nb[5]+cell_center[index_nb[16]].m_state.m_U[2]+U2_center)/4.0;
    }


    //corner (i-1/2,j,k+1/2)

    if (!bNoBoundary[0] && bNoBoundary[5])
    {
	U0_nb[17] = U0_nb[0];
	U1_nb[17] = U1_nb[0];
	U2_nb[17] = U2_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[5])
    {
	U0_nb[17] = U0_nb[5];
	U1_nb[17] = U1_nb[5];
	U2_nb[17] = U2_nb[5];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[5])
    {
	U0_nb[17] = U0_nb[0];
	U1_nb[17] = U1_nb[0];
	U2_nb[17] = U2_nb[0];
    }
    else
    {
	U0_nb[17] = (U0_nb[0]+U0_nb[5]+cell_center[index_nb[17]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[17] = (U1_nb[0]+U1_nb[5]+cell_center[index_nb[17]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[17] = (U2_nb[0]+U2_nb[5]+cell_center[index_nb[17]].m_state.m_U[2]+U2_center)/4.0;
    }

    diffusion[0] += (mu_edge[2]*U1_nb[6]- mu_edge[2]*U1_nb[7]+ mu_edge[3]*U1_nb[8]- mu_edge[3]*U1_nb[9])/ (top_h[0]*top_h[1]);//(mu*v_x)_y
    diffusion[0] += (mu_edge[4]*U2_nb[14]-mu_edge[4]*U2_nb[15]+mu_edge[5]*U2_nb[16]-mu_edge[5]*U2_nb[17])/(top_h[0]*top_h[2]);//(mu*w_x)_z

    diffusion[1] += (mu_edge[0]*U0_nb[6]- mu_edge[1]*U0_nb[7]+ mu_edge[1]*U0_nb[8]- mu_edge[0]*U0_nb[9])/ (top_h[0]*top_h[1]);//(mu*u_y)_x
    diffusion[1] += (mu_edge[4]*U2_nb[10]-mu_edge[4]*U2_nb[11]+mu_edge[5]*U2_nb[12]-mu_edge[5]*U2_nb[13])/(top_h[1]*top_h[2]);//(mu*w_y)_z


    diffusion[2] += (mu_edge[0]*U0_nb[14]-mu_edge[1]*U0_nb[15]+mu_edge[1]*U0_nb[16]-mu_edge[0]*U0_nb[17])/(top_h[0]*top_h[2]);//(mu*u_z)_x
    diffusion[2] += (mu_edge[2]*U1_nb[10]-mu_edge[3]*U1_nb[11]+mu_edge[3]*U1_nb[12]-mu_edge[2]*U1_nb[13])/(top_h[1]*top_h[2]);//(mu*v_z)_y

}


/**
* calculate Uxx,Uyy,Uzz and Px,Py,Pz.
* @param dir
* @param icoords
* @param dU2
* @param dP
*/
void Incompress_Solver_Smooth_3D_Cartesian::getDU2(
	int *icoords,
	EBM_COORD xyz,
	double dU2[3])
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

    dU2[0] = ((U2.m_U[0] - U1.m_U[0])/dh1 - (U1.m_U[0] - U0.m_U[0])/dh0) / dh;
    dU2[1] = ((U2.m_U[1] - U1.m_U[1])/dh1 - (U1.m_U[1] - U0.m_U[1])/dh0) / dh;
    dU2[2] = ((U2.m_U[2] - U1.m_U[2])/dh1 - (U1.m_U[2] - U0.m_U[2])/dh0) / dh;

}

// Minmod slope limiter
void Incompress_Solver_Smooth_3D_Cartesian::getLimitedSlope(
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
    else  //xyz == COORD_Z
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

// Van Leer Slope Limiter


void Incompress_Solver_Smooth_3D_Cartesian::getLimitedSlope_Vanleer(
	int *icoords,
	EBM_COORD xyz,
	double slope[3])
{
    double dx,dy,dz;
    L_STATE U0, U1, U2;
    double u_lim, v_lim, w_lim;
    double u_slope, v_slope, w_slope;

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);

    U1 = cell_center[index].m_state;
    dx = top_h[0];
    dy = top_h[1];
    dz = top_h[2];

    bool bNoBoundary[2];


    if(xyz==COORD_X)
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,WEST,U0,m_t_old);
	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,EAST,U2,m_t_old);
	if(bNoBoundary[0] || bNoBoundary[1])
	{

	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dx, 2.0*(U1.m_U[0]-U0.m_U[0])/dx);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dx, 2.0*(U1.m_U[1]-U0.m_U[1])/dx);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dx, 2.0*(U1.m_U[2]-U0.m_U[2])/dx);

	    u_slope = (U2.m_U[0] - U0.m_U[0])/(2*dx);
	    v_slope = (U2.m_U[1] - U0.m_U[1])/(2*dx);
	    w_slope = (U2.m_U[2] - U0.m_U[2])/(2*dx);
	    

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else if ( (!bNoBoundary[0]) || bNoBoundary[1])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dx, 2.0*(U1.m_U[0]-U0.m_U[0])/(dx/2.0));
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dx, 2.0*(U1.m_U[1]-U0.m_U[1])/(dx/2.0));
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dx, 2.0*(U1.m_U[2]-U0.m_U[2])/(dx/2.0));

	    u_slope = (U2.m_U[0] + U1.m_U[0] + 2.0*U0.m_U[0])/(2.0*dx);
	    v_slope = (U2.m_U[1] + U1.m_U[1] + 2.0*U0.m_U[1])/(2.0*dx);
	    w_slope = (U2.m_U[2] + U1.m_U[2] + 2.0*U0.m_U[2])/(2.0*dx);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);
	}
	else if ( (!bNoBoundary[1]) || bNoBoundary[0])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/(dx/2.0), 2.0*(U1.m_U[0]-U0.m_U[0])/dx);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/(dx/2.0), 2.0*(U1.m_U[1]-U0.m_U[1])/dx);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/(dx/2.0), 2.0*(U1.m_U[2]-U0.m_U[2])/dx);

	    u_slope = (2.0*U2.m_U[0] - U1.m_U[0] - U0.m_U[0])/(2.0*dx);
	    v_slope = (2.0*U2.m_U[1] - U1.m_U[1] - U0.m_U[1])/(2.0*dx);
	    w_slope = (2.0*U2.m_U[2] - U1.m_U[2] - U0.m_U[2])/(2.0*dx);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else
	{
	    printf("\nThe number of cell in x direction is less than 1!!\n");
	    slope[0] = slope[1] = slope[2] = 0.0;
	}

    }
    else if (xyz == COORD_Y)	//
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,SOUTH,U0,m_t_old);
	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,NORTH,U2,m_t_old);

	if(bNoBoundary[0] || bNoBoundary[1])
	{

	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dy, 2.0*(U1.m_U[0]-U0.m_U[0])/dy);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dy, 2.0*(U1.m_U[1]-U0.m_U[1])/dy);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dy, 2.0*(U1.m_U[2]-U0.m_U[2])/dy);

	    u_slope = (U2.m_U[0] - U0.m_U[0])/(2*dy);
	    v_slope = (U2.m_U[1] - U0.m_U[1])/(2*dy);
	    w_slope = (U2.m_U[2] - U0.m_U[2])/(2*dy);
	    

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else if ( (!bNoBoundary[0]) || bNoBoundary[1])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dy, 2.0*(U1.m_U[0]-U0.m_U[0])/(dy/2.0));
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dy, 2.0*(U1.m_U[1]-U0.m_U[1])/(dy/2.0));
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dy, 2.0*(U1.m_U[2]-U0.m_U[2])/(dy/2.0));

	    u_slope = (U2.m_U[0] + U1.m_U[0] + 2.0*U0.m_U[0])/(2.0*dy);
	    v_slope = (U2.m_U[1] + U1.m_U[1] + 2.0*U0.m_U[1])/(2.0*dy);
	    w_slope = (U2.m_U[2] + U1.m_U[2] + 2.0*U0.m_U[2])/(2.0*dy);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);
	}
	else if ( (!bNoBoundary[1]) || bNoBoundary[0])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/(dy/2.0), 2.0*(U1.m_U[0]-U0.m_U[0])/dy);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/(dy/2.0), 2.0*(U1.m_U[1]-U0.m_U[1])/dy);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/(dy/2.0), 2.0*(U1.m_U[2]-U0.m_U[2])/dy);

	    u_slope = (2.0*U2.m_U[0] - U1.m_U[0] - U0.m_U[0])/(2.0*dy);
	    v_slope = (2.0*U2.m_U[1] - U1.m_U[1] - U0.m_U[1])/(2.0*dy);
	    w_slope = (2.0*U2.m_U[2] - U1.m_U[2] - U0.m_U[2])/(2.0*dy);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else
	{
	    printf("\nThe number of cell in y direction is less than 1!!\n");
	    slope[0] = slope[1] = slope[2] = 0.0;
	}
    }
    else // xyz == COORD_Z
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,LOWER,U0,m_t_old);
	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,UPPER,U2,m_t_old);

	if(bNoBoundary[0] || bNoBoundary[1])
	{

	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dz, 2.0*(U1.m_U[0]-U0.m_U[0])/dz);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dz, 2.0*(U1.m_U[1]-U0.m_U[1])/dz);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dz, 2.0*(U1.m_U[2]-U0.m_U[2])/dz);

	    u_slope = (U2.m_U[0] - U0.m_U[0])/(2*dz);
	    v_slope = (U2.m_U[1] - U0.m_U[1])/(2*dz);
	    w_slope = (U2.m_U[2] - U0.m_U[2])/(2*dz);
	    

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else if ( (!bNoBoundary[0]) || bNoBoundary[1])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dz, 2.0*(U1.m_U[0]-U0.m_U[0])/(dz/2.0));
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dz, 2.0*(U1.m_U[1]-U0.m_U[1])/(dz/2.0));
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dz, 2.0*(U1.m_U[2]-U0.m_U[2])/(dz/2.0));

	    u_slope = (U2.m_U[0] + U1.m_U[0] + 2.0*U0.m_U[0])/(2.0*dz);
	    v_slope = (U2.m_U[1] + U1.m_U[1] + 2.0*U0.m_U[1])/(2.0*dz);
	    w_slope = (U2.m_U[2] + U1.m_U[2] + 2.0*U0.m_U[2])/(2.0*dz);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);
	}
	else if ( (!bNoBoundary[1]) || bNoBoundary[0])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/(dz/2.0), 2.0*(U1.m_U[0]-U0.m_U[0])/dz);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/(dz/2.0), 2.0*(U1.m_U[1]-U0.m_U[1])/dz);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/(dz/2.0), 2.0*(U1.m_U[2]-U0.m_U[2])/dz);

	    u_slope = (2.0*U2.m_U[0] - U1.m_U[0] - U0.m_U[0])/(2.0*dz);
	    v_slope = (2.0*U2.m_U[1] - U1.m_U[1] - U0.m_U[1])/(2.0*dz);
	    w_slope = (2.0*U2.m_U[2] - U1.m_U[2] - U0.m_U[2])/(2.0*dz);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else
	{
	    printf("\nThe number of cell in z direction is less than 1!!\n");
	    slope[0] = slope[1] = slope[2] = 0.0;
	}
    }
}


double Incompress_Solver_Smooth_3D_Cartesian::EBM_minmod(
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
bool Incompress_Solver_Smooth_3D_Cartesian::getNeighborOrBoundaryState(
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


/**
*
* @param state_left
* @param state_right
* @param ans
*/
void Incompress_Solver_Smooth_3D_Cartesian::getRiemannSolution(
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

