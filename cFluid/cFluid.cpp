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


/*
*				cFluid.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a double vortex field. It demonstrates
*	the resolution of the front tracking method.
*
*/

#include "cFluid.h"

	/*  Function Declarations */
static void gas_driver(Front*,G_CARTESIAN&);
static void cft_driver(Front*,G_CARTESIAN&);
static int g_cartesian_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReadFromInput;
int RestartStep;
boolean binary = NO;


/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	static EQN_PARAMS eqn_params;

	G_CARTESIAN	g_cartesian(front);


	/* Initialize basic computational data */

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);
	
        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_state_name,"%s/state.ts%s",restart_name,
                        right_flush(RestartStep,7));
        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
			right_flush(RestartStep,7));
	if (pp_numnodes() > 1)
	{
            sprintf(restart_name,"%s-nd%s",restart_name,
			right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-nd%s",restart_state_name,
                        right_flush(pp_mynode(),4));
	}

	eqn_params.dim = f_basic.dim;
	front.extra1 = (POINTER)&eqn_params;
	FT_ReadSpaceDomain(in_name,&f_basic);
	read_cFluid_params(in_name,&eqn_params);
	if (debugging("trace")) printf("Passed read_cFluid_params()\n");

	g_cartesian.setProbParams(&eqn_params, f_basic);
	if (debugging("trace"))
	    printf("Passed g_cartesian.setProbParams()\n");

	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);
	if (debugging("sample_velocity"))
            g_cartesian.initSampleVelocity(in_name);

	if (debugging("trace")) printf("Passed FT_StartUp()\n");

	read_movie_options(in_name,&eqn_params);
	if (eqn_params.dim != 1)
	    read_statistics_options(in_name,&eqn_params);

	/* Initialize interface through level function */

	if (!RestartRun)
	{
	    g_cartesian.setInitialIntfc(&level_func_pack,in_name);
	    if (f_basic.dim == 3) level_func_pack.set_3d_bdry = YES;
	    FT_InitIntfc(&front,&level_func_pack);

	    FT_PromptSetMixedTypeBoundary2d(in_name,&front);
	    read_dirichlet_bdry_data(in_name,&front);
	    if (f_basic.dim < 3)
	    	FT_ClipIntfcToSubdomain(&front);
	    FT_RedistMesh(&front);
	    if (debugging("trace"))
	    {
                char test_name[100];
                (void) printf("Passed FT_InitIntfc()\n");
                switch (f_basic.dim)
                {
                case 2:
                    sprintf(test_name,"init_intfc-%d.xg",pp_mynode());
                    xgraph_2d_intfc(test_name,front.interf);
                    break;
                case 3:
                    gview_plot_interface("gv-init",front.interf);
                    break;
                }
	    }
	}
	else
	    restart_set_dirichlet_bdry_function(&front);

	/* Initialize velocity field function */

	front._compute_force_and_torque = cfluid_compute_force_and_torque;
	velo_func_pack.func_params = (POINTER)&g_cartesian;
	velo_func_pack.func = g_cartesian_vel;
	velo_func_pack.point_propagate = cFluid_point_propagate;
	FT_InitVeloFunc(&front,&velo_func_pack);
	if (debugging("trace"))
	    printf("Passed FT_InitVeloFunc()\n");

	FT_MakeGridIntfc(&front);
	g_cartesian.initMesh(
		topological_grid(front.grid_intfc).L,
		topological_grid(front.grid_intfc).U,
		topological_grid(front.grid_intfc).gmax,
		front.rect_grid->lbuf,
		front.rect_grid->ubuf
		);
	FT_FreeGridIntfc(&front);

	if (debugging("trace"))
	    printf("Passed g_cartesian.initMesh()\n");

	FT_AddMovieFrame(&front,out_name,binary);

	if (RestartRun)
	{
	    readFrontStates(&front,restart_state_name);
	    g_cartesian.readInteriorStates(restart_state_name);
	}
	else
	{
	    g_cartesian.setInitialStates();
	}
	if (debugging("trace"))
	    printf("Passed state initialization()\n");

	/* Propagate the front */

	//gas_driver(&front, g_cartesian);
	cft_driver(&front, g_cartesian);

	//PetscFinalize();
	clean_up(0);
}

static  void cft_driver(
        Front *front,
	G_CARTESIAN &g_cartesian)
{
	printf("\nCalling cft_driver.\n");
	printf("Conservative front tracking is now for RM only.\n");
	printf("Delta h in 3 directions should be equivalent to each other.\n");

	if (g_cartesian.dim != 3)
	{
	    printf("Conservative front tracking is implemented for 3D simulations only.\n");
	    clean_up(ERROR);
	}

        double CFL;
	double dt;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

	if (!RestartRun)
	{
	    FT_ResetTime(front);

	    //debugdan	FIXME
	    //g_cartesian.printInteriorVtk(out_name);
	    //vtk_interface_plot("debug-intfc/intfc0",front->interf,FALSE,0,0);
	    //debugdan	FIXME

	    //CFT
	    dt = 0.5*front->dt;
	    front->dt = dt;

	    //set old t cells
	    g_cartesian.cft_init_cut_cells(OLDTS);
	    g_cartesian.cft_set_cut_cells(OLDTS);
	    g_cartesian.cft_merge_polyhs(OLDTS);

	    //initialization
	    g_cartesian.cft_set_init_polyh_states();

	    //scatter states
	    g_cartesian.cft_scatMeshStates();
	    //g_cartesian.copyMeshStates();	???	FIXME

	    //first propagation
	    FrontPreAdvance(front);
	    FT_Propagate(front);
	    g_cartesian.cft_init_cut_cells(HALFTS);
	    g_cartesian.cft_set_cut_cells(HALFTS);
	    g_cartesian.cft_merge_polyhs(HALFTS);

	    //second propagation
	    FrontPreAdvance(front);
	    FT_Propagate(front);
	    g_cartesian.cft_init_cut_cells(NEWTS);
	    g_cartesian.cft_set_cut_cells(NEWTS);
	    g_cartesian.cft_merge_polyhs(NEWTS);

	    //solve
	    g_cartesian.cft_solve(2*dt);
	    //g_cartesian.cft_set_face_flux();
	    //g_cartesian.cft_update_states();
	    g_cartesian.cft_update_states_new();
	    //CFT

	    //g_cartesian.printInteriorVtk(out_name);
	    //vtk_interface_plot("debug-intfc/intfc1",front->interf,FALSE,0,0);

	    FT_SetTimeStep(front);
	    dt = front->dt = std::min(front->dt,CFL*g_cartesian.max_dt);
	    FT_SetOutputCounter(front);
        }
        else
	{
	    printf("Restart has not been setup for conservative front tracking.\n");
	    clean_up(ERROR);
	    /*
	    g_cartesian.time = front->time;
	    g_cartesian.step = front->step;
	    dt = front->dt = std::min(front->dt,CFL*g_cartesian.max_dt);
	    FT_SetOutputCounter(front);
	    */
	}

	FT_TimeControlFilter(front);
	printf("\ntime = %20.14f   step = %5d   next dt = %20.14f\n",
		g_cartesian.time, g_cartesian.step, dt);

	if (!RestartRun && g_cartesian.dim != 1)
	{
	    g_cartesian.record_intfc_extrema();
	    g_cartesian.print_intfc_extrema(out_name);
	}

        while (true)
        {
	    if (g_cartesian.dim != 1)
		g_cartesian.record_intfc_extrema();

	    if (front->step == 0)
            {
                g_cartesian.initMovieVariables();
                FT_AddMovieFrame(front,out_name,binary);
		g_cartesian.printInteriorVtk(out_name);
            }

	    //reset cells
	    //new cells to old cells
	    g_cartesian.cft_newts_cut_cells();

	    //debugdan	FIXME
	    //g_cartesian.cft_check_mass_oldts();
	    //debugdan	FIXME

	    dt = 0.5*dt;
	    front->dt = dt;

	    //scatter states
	    g_cartesian.cft_scatMeshStates();
	    //g_cartesian.copyMeshStates();	???	FIXME

	    //first propagation
	    //printf("First propagate.\n");
	    FrontPreAdvance(front);
	    FT_Propagate(front);
	    g_cartesian.cft_init_cut_cells(HALFTS);
	    g_cartesian.cft_set_cut_cells(HALFTS);
	    g_cartesian.cft_merge_polyhs(HALFTS);

	    //second propagation
	    //printf("Second propagate.\n");
	    FrontPreAdvance(front);
	    FT_Propagate(front);
	    g_cartesian.cft_init_cut_cells(NEWTS);
	    g_cartesian.cft_set_cut_cells(NEWTS);
	    g_cartesian.cft_merge_polyhs(NEWTS);

	    //solve
	    g_cartesian.cft_solve(2*dt);
	    //g_cartesian.cft_set_face_flux();
	    g_cartesian.cft_update_states_new();

	    dt = 2*dt;
	    front->dt = dt;
	    //CFT

	    /*
	    FrontPreAdvance(front);
	    FT_Propagate(front);

	    g_cartesian.solve(dt);
	    */

	    FT_AddTimeStepToCounter(front);
	    ++g_cartesian.step;
	    g_cartesian.time += dt;

	    //debugdan	FIXME
	    g_cartesian.cft_check_mass();
	    //debugdan	FIXME
				
	    FT_SetTimeStep(front);
            dt = front->dt = std::min(front->dt,CFL*g_cartesian.max_dt);
	
            /* Output section */
            if (FT_IsSaveTime(front))
	    {
            	FT_Save(front,out_name);
		g_cartesian.printFrontInteriorStates(out_name);
	    }
            if (FT_IsMovieFrameTime(front))
	    {
		if (g_cartesian.dim != 1)
		    g_cartesian.print_intfc_extrema(out_name);
	        g_cartesian.initMovieVariables();
            	FT_AddMovieFrame(front,out_name,binary);
		if (g_cartesian.dim == 1)
		    g_cartesian.printOneDimStates(out_name);
		else
		    g_cartesian.printInteriorVtk(out_name);
	    }

            if (FT_TimeLimitReached(front))
	    {
		if (!FT_IsSaveTime(front))
		{
            	    FT_Save(front,out_name);
		    g_cartesian.printFrontInteriorStates(out_name);
		}
		if (!FT_IsMovieFrameTime(front))
		{
                    g_cartesian.initMovieVariables();
                    FT_AddMovieFrame(front,out_name,binary);
		}
		(void) printf("\ntime = %20.14f   step = %5d   ",
                                front->time,front->step);
                (void) printf("next dt = %20.14f\n",front->dt);
                break;
	    }
	    FT_TimeControlFilter(front);
	    dt = front->dt;
	    (void) printf("\ntime = %20.14f   step = %5d   next dt = %20.14f\n",
                        front->time,front->step,front->dt);

            fflush(stdout);
        }

	printf("End of cft_driver.\n");
}

static  void gas_driver(
        Front *front,
	G_CARTESIAN &g_cartesian)
{
        double CFL;
	double dt;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

	if (!RestartRun)
	{
	    FT_ResetTime(front);

	    if (debugging("trace"))
		printf("Calling initial FT_Propagate()\n");

	    FrontPreAdvance(front);
	    FT_Propagate(front);
	    g_cartesian.solve(front->dt);

	    //debugdan	FIXME
	    //g_cartesian.cft_init_cut_cells(NEWTS);
	    //g_cartesian.cft_set_cut_cells(NEWTS);
	    //g_cartesian.cft_free_cells(NEWTS);
	    //debugdan	FIXME

	    FT_SetTimeStep(front);
	    dt = front->dt = std::min(front->dt,CFL*g_cartesian.max_dt);
	    FT_SetOutputCounter(front);
        }
        else
	{
	    g_cartesian.time = front->time;
	    g_cartesian.step = front->step;
	    dt = front->dt = std::min(front->dt,CFL*g_cartesian.max_dt);
	    FT_SetOutputCounter(front);
	}

	FT_TimeControlFilter(front);
	assert(g_cartesian.time == front->time);
	assert(g_cartesian.step == front->step);
	(void) printf("\ntime = %20.14f   step = %5d   next dt = %20.14f\n",
                        g_cartesian.time, g_cartesian.step, dt);

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));
	}

	if (!RestartRun && g_cartesian.dim != 1)
	{
	    g_cartesian.record_intfc_extrema();
	    g_cartesian.print_intfc_extrema(out_name);
	}

	//Dan
	/*
	if (true && !RestartRun && g_cartesian.dim == 3)
	{
	    g_cartesian.cvol();
	}
	*/
	//Dan

	if (debugging("trace")) printf("Before time loop\n");
        for (;;)
        {
	    if (g_cartesian.dim != 1)
		g_cartesian.record_intfc_extrema();

	    if (front->step == 0 && g_cartesian.dim != 1)
            {
                g_cartesian.initMovieVariables();
                FT_AddMovieFrame(front,out_name,binary);
		g_cartesian.printInteriorVtk(out_name);
            }

	    print_storage("Storage at start of time step","step_storage");
	    if (debugging("trace")) printf("Begin a time step\n");
	    FrontPreAdvance(front);
	    FT_Propagate(front);

	    if (debugging("trace")) printf("Begin calling solve()\n");
	    g_cartesian.solve(dt);
	    if (debugging("trace")) 
	    {
		printf("Passed solve()\n");
		print_storage("Storage after time step","trace");
	    }

	    //debugdan	FIXME
	    //g_cartesian.cft_init_cut_cells(NEWTS);
	    //g_cartesian.cft_set_cut_cells(NEWTS);
	    //debugdan	FIXME

	    FT_AddTimeStepToCounter(front);
	    ++g_cartesian.step;
	    g_cartesian.time += dt;
	    assert(g_cartesian.time==front->time);
	    assert(g_cartesian.step==front->step);

	    //debugdan	FIXME
	    g_cartesian.ncft_check_mass();
	    //g_cartesian.cft_free_cells(NEWTS);
	    //debugdan	FIXME
				
            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
	    if (debugging("step_size"))
	    {
		(void) printf("Step size from front:    %20.14f\n",front->dt);
		(void) printf("Step size from interior: %20.14f\n",
					CFL*g_cartesian.max_dt);
	    }
            dt = front->dt = std::min(front->dt,CFL*g_cartesian.max_dt);

            /* Output section */

            if (FT_IsSaveTime(front))
	    {
            	FT_Save(front,out_name);
		g_cartesian.printFrontInteriorStates(out_name);
	    }
            if (FT_IsMovieFrameTime(front))
	    {
		if (g_cartesian.dim != 1)
		    g_cartesian.print_intfc_extrema(out_name);
	    	if (debugging("trace")) 
		    printf("Calling initMovieVariable()\n");
	        g_cartesian.initMovieVariables();
	    	if (debugging("trace")) 
		    printf("Calling FT_AddMovieFrame()\n");
            	FT_AddMovieFrame(front,out_name,binary);
		if (g_cartesian.dim == 1)
		    g_cartesian.printOneDimStates(out_name);
		else
		    g_cartesian.printInteriorVtk(out_name);
	    }

            if (FT_TimeLimitReached(front))
	    {
		if (!FT_IsSaveTime(front))
		{
            	    FT_Save(front,out_name);
		    g_cartesian.printFrontInteriorStates(out_name);
		}
		if (!FT_IsMovieFrameTime(front))
		{
		    if (debugging("trace")) 
                    	printf("Calling initMovieVariable()\n");
                    g_cartesian.initMovieVariables();
                    if (debugging("trace"))
                    	printf("Calling FT_AddMovieFrame()\n");
                    FT_AddMovieFrame(front,out_name,binary);
		}
		(void) printf("\ntime = %20.14f   step = %5d   ",
                                front->time,front->step);
                (void) printf("next dt = %20.14f\n",front->dt);
                break;
	    }
	    FT_TimeControlFilter(front);
	    dt = front->dt;
	    (void) printf("\ntime = %20.14f   step = %5d   next dt = %20.14f\n",
                        front->time,front->step,front->dt);

	    //Dan
	    /*
	    if (front->step > 0 && true && !RestartRun && g_cartesian.dim == 3)
	    {
		g_cartesian.cvol();
//		if (front->step > 0)
//		    exit(0);	//FIXME
	    }
	    */
	    //Dan

            fflush(stdout);
        }
	if (debugging("trace")) printf("After time loop\n");
}       /* end gas_driver */

static int g_cartesian_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	double *coords = Coords(p);
	((G_CARTESIAN*)params)->getVelocity(coords, vel);
	return YES;
}	/* end g_cartesian_vel */

