#include "solver.h"

void poisson_solver2d(
	Front *front,	    /* front structure containing interface geometry */
	int ilower,
	int iupper,
	int **ij_to_I,
	double *source,	    /* source term in poisson's equation */
	double *k,	    /* coefficient function in poisson's equation */
	double *soln,	    /* solution of the poisson's equation */
	double *max_soln,   /* maximum value of the solution */
	double *min_soln)   /* minimum value of the solution */
{
	int index,index_nb[4],size;
	double k0,k_nb[4];
	double rhs,coeff[4];
	int I,I_nb[4];
	int i,j,l,icoords[MAXD];
	int imin,imax,jmin,jmax;
	COMPONENT comp;
	double aII;
	int num_nb;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;
	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
	struct Table *T = table_of_interface(front->grid_intfc);
	COMPONENT *top_comp = T->components;
	int *lbuf = front->rect_grid->lbuf;
	int *ubuf = front->rect_grid->ubuf;
	int *top_gmax = rgr->gmax;
	double *top_h = rgr->h;
	HYPER_SURF *hs;

	PETSc solver;
	solver.Create(ilower, iupper-1, 5, 0);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;
	*max_soln = -HUGE;
	*min_soln = HUGE;

	imin = (lbuf[0] == 0) ? 1 : lbuf[0];
        jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
        jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    I = ij_to_I[i][j];
	    if (I == -1) continue;

	    index_nb[0] = d_index2d(i-1,j,top_gmax);
	    index_nb[1] = d_index2d(i+1,j,top_gmax);
	    index_nb[2] = d_index2d(i,j-1,top_gmax);
	    index_nb[3] = d_index2d(i,j+1,top_gmax);
	    I_nb[0] = ij_to_I[i-1][j];
	    I_nb[1] = ij_to_I[i+1][j];
	    I_nb[2] = ij_to_I[i][j-1];
	    I_nb[3] = ij_to_I[i][j+1];
	    icoords[0] = i;
	    icoords[1] = j;
	
	    k0 = k[index];
	    num_nb = 0;
	    for (l = 0; l < 4; ++l)
	    {
		if (I_nb[l] == -1)
		    index_nb[l] = index;
		else num_nb++;
		k_nb[l] = 0.5*(k0 + k[index_nb[l]]);
	    	coeff[l] = 1.0/k_nb[l]/(top_h[l/2]*top_h[l/2]); 
	    }

	    rhs = source[index];

	    aII = 0.0;
	    for (l = 0; l < 4; ++l)
	    {
		if (num_nb < 2) break;
	    	if (I_nb[l] != -1)
		{
		    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
		}
		else
		{
		    hs = FT_HyperSurfAtGridCrossing(front,icoords,dir[l]);
		    if (hs != NULL && wave_type(hs) == DIRICHLET_BOUNDARY)
		    {
			if (!boundary_state(hs))
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
	    }
	    /*
	     * This change reflects the need to treat point with only one
	     * interior neighbor (a convex point). Not sure why PETSc cannot
	     * handle such case. If we have better understanding, this should
	     * be changed back.
	     */
	    if(num_nb >= 2)
	    {
                solver.Set_A(I,I,aII);
	    }
            else
            {
	
                solver.Set_A(I,I,1.0);
		rhs = soln[index];
            }
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-14);

	start_clock("Before Petsc Solver");
	if (use_neumann_solver)
	{
	    printf("\nUsing Neumann Solver!\n");
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if(rel_residual > 1)
	    {
		printf("\n The solution diverges! The residual "
		       "is %le. Solve again using GMRES!\n",rel_residual);
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
		printf("\n The solution diverges! The residual "
		       "is %le. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	stop_clock("After Petsc Solver");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %le \n", 
			num_iter, rel_residual);
	
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    I = ij_to_I[i][j];
	    soln[index] = x[I-ilower];
	    if (*max_soln < soln[index]) *max_soln = soln[index];
	    if (*min_soln > soln[index]) *min_soln = soln[index];
	}
	FT_ParallelExchGridArrayBuffer(soln,front);
	pp_global_max(max_soln,1);
	pp_global_min(min_soln,1);

	if(debugging("step_size"))
	{
	    printf("\nThe max solution value is %.16g\n",*max_soln);
	    printf("\nThe min solution value is %.16g\n",*min_soln);
	}

	FT_FreeThese(1,x);
}	/* end computeProjection2d */
