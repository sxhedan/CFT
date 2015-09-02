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
*				ifluid.c
* This program is modified from example0.c for solving incompressible flow.
* 
* The solver is define in lcartsn.h/c.
*
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a double vortex field. It demonstrates
*	the resolution of the front tracking method.
*
*/

#include "iFluid.h"
#include "ifluid_basic.h"
#include "iFluid_debug.h"
#include "Graph.h"
#include "iFDroplet.h"

	/*  Function Declarations */
static void init_io( int,char**);
static void ifluid_driver(Front*,Incompress_Solver_Smooth_Basis*,F_BASIC_DATA*);
static void ifluid_driver_debug(Front*,Incompress_Solver_Smooth_Basis*);
//Accuracy testing ifluid_driver

static int l_cartesian_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);

/************************************************************
 * 	init the whole domain as a single component.
 * This is used as testing.
 */
void initIntfc_oneComponent(Front *front,LEVEL_FUNC_PACK *level_func_pack);
double level_oneComponent_2D(POINTER func_params,double *coords);
double level_oneComponent_3D(POINTER func_params,double *coords);
double level_oneComponent_3D_Cylindrical(POINTER func_params,double *coords);
/******************************************************/

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean RegridRun;
boolean RegridRestart;
boolean ReadFromInput;
int RestartStep;
boolean binary = YES;
boolean isTesting = NO; //The switch for the testing mode
//vector<int>  restart_step;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	IF_PARAMS iFparams;
	IF_PROB_TYPE prob_type;

	/* Initialize basic computational data */

	FT_Init(argc,argv,&f_basic);//Read parameters from command line
	f_basic.size_of_intfc_state = sizeof(STATE);
	
	//Initialize Petsc before the FrontStartUp
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	//if (debugging("trace")) printf("Passed PetscInitialize()\n");

	/*Construct Incompress Solver l_cartesian*/

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;


	Incompress_Solver_Smooth_Basis *l_cartesian;//Construct Incompress Solver



	if(f_basic.dim == 2)
	{
	    if (isTesting == YES)
	    {
		printf("\nUsing the debugging class!\n");
		l_cartesian = new Incompress_Solver_Smooth_2D_Cartesian_Debug(front);
	    }
	    else
		l_cartesian = new Incompress_Solver_Smooth_2D_Cartesian(front);
	}
	else if(f_basic.dim == 3)
	{
	    if (f_basic.coord_system == CYLINDRICAL_REMAP) 
	    //Cylindrical coordinate in 3D
	    {
		if (isTesting == YES)//Testing mode using exact solution
		{
		    printf("\nUsing the debugging class!\n");
		    l_cartesian = new Incompress_Solver_Smooth_3D_Cylindrical_Debug(front);
		}
		else
		    l_cartesian = new Incompress_Solver_Smooth_3D_Cylindrical(front);
	    }
	    else if(f_basic.coord_system == SPHERICAL_REMAP)
	    //Sperical coordinate in 3D, not implemented yet
	    {
	    }
	    else 
	    //Default: Rectangular coordinate in 3D
	    {
		if (isTesting == YES)
		{
		    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian_Debug(front);
		    printf("\nUsing the debugging class!\n");
		}
		else
		    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);
	    }
	}


        sprintf(restart_state_name,"%s/restart/state.t%s",restart_name,
                        right_flush(RestartStep,7));
        sprintf(restart_name,"%s/restart/intfc-t%s",restart_name,
                        right_flush(RestartStep,7));

	if (pp_numnodes() > 1)
	{
            sprintf(restart_name,"%s-p%s",restart_name,
			right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-p%s",restart_state_name,
                        right_flush(pp_mynode(),4));
	}

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

        RegridRun = f_basic.RegridRun;
        RegridRestart = f_basic.RegridRestart;

	if (debugging("trace")) printf("Passed FT_StartUp()\n");
	iFparams.dim = f_basic.dim;
	front.extra1 = (POINTER)&iFparams;

	if (isTesting != YES)
	    read_iF_prob_type(in_name,&prob_type);

	read_iFparams(in_name,&iFparams); 
	//read the projection type (current B.C. support BELL_COLELLA)
	
	read_iF_movie_options(in_name,&iFparams);
	if (debugging("trace")) printf("Passed read_iFparams()\n");

	/* Initialize interface through level function */

	if (isTesting == YES)
	    initIntfc_oneComponent(&front, &level_func_pack);
	else
	    setInitialIntfc(&front,&level_func_pack,in_name,prob_type);
	    //set iFparams according to PROB_TYPE

	if (debugging("trace")) printf("Passed setInitialIntfc()\n");

	if (!RestartRun)
	{
	    if (f_basic.dim == 3) level_func_pack.set_3d_bdry = YES;
	    FT_InitIntfc(&front,&level_func_pack);
	    if (debugging("trace"))
	    {
		char test_name[100];

		printf("Passed FT_InitIntfc()\n");
		switch (f_basic.dim)
		{
		case 2:
		    sprintf(test_name,"init_intfc-%d.xg",pp_mynode());
		    xgraph_2d_intfc(test_name,front.interf);
		    break;
		case 3:
		    sprintf(test_name,"init_intfc-%d.xg",pp_mynode());
		    //gview_plot_interface("gv-init",front.interf);
		    break;
		}
	    }
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (f_basic.dim < 3)
	    	FT_ClipIntfcToSubdomain(&front);
	    if (debugging("trace")) 
		printf("Passed read_iF_dirichlet_bdry_data()\n");
	}

	/* Initialize velocity field function */
	velo_func_pack.func_params = (POINTER)l_cartesian;
	velo_func_pack.func = l_cartesian_vel;
	velo_func_pack.point_propagate = ifluid_point_propagate;
	FT_InitVeloFunc(&front,&velo_func_pack);
	if (debugging("trace"))
	    printf("Passed FT_InitVeloFunc()\n");

	l_cartesian->initMesh();



	if (debugging("sample_velocity"))
	    l_cartesian->initSampleVelocity(in_name);

	init_fluid_state_func(l_cartesian,prob_type);
	if (debugging("trace"))
	    printf("Passed l_cartesian.initMesh()\n");

	if (RestartRun)
	    l_cartesian->readFrontInteriorStates(restart_state_name, binary, RegridRestart);
	else
	    l_cartesian->setInitialCondition();

	if (debugging("trace"))
            printf("Passed state initialization()\n");


	/* Enter the iFluid Driver */

	if (isTesting == YES)
	{
	    printf("\nEntering ifluid driver for testing problem!\n");
	    ifluid_driver_debug(&front, l_cartesian);
	}
	else
	    ifluid_driver(&front, l_cartesian, &f_basic);

	PetscFinalize();
	clean_up(0);
}


static  void ifluid_driver(
        Front *front,
	Incompress_Solver_Smooth_Basis *l_cartesian,
        F_BASIC_DATA *f_basic)
{
    	printf("\nnode = %d, PID = %d\n",pp_mynode(), getpid());

	CellCornerGraph majority(l_cartesian);
	DropletAnalysis droplet(l_cartesian);

        double CFL;
	int i,dim = front->rect_grid->dim;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

	if (!RestartRun || RegridRestart)
	{
	    FT_RedistMesh(front);
	}

        if (!RestartRun)
        {
	    FT_ResetTime(front);
	    if (debugging("trace"))
		printf("Zeroth step: Before FT_Propagate() front->dt = %f,"
			"l_cartesian->max_dt = %f\n",front->dt, l_cartesian->max_dt);

	    ((Incompress_Solver_Smooth_3D_Cylindrical*)l_cartesian)->adjustReflectionBuffer();
            FT_Propagate(front);
	    
	    if (debugging("trace")) printf("Zeroth step: Calling ifluid solve()\n");
            l_cartesian->solve(front->dt);
	    if (debugging("trace")) printf("Zeroth step: Passed ifluid solve()\n");

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
            FT_SetOutputCounter(front);

        }
        else
        {
	    FT_SetOutputCounter(front);
        }

	FT_TimeControlFilter(front);


	/**** The zero time step output  ****/

       if(!RestartRun)
       {	   
	   l_cartesian->multiDiagnosisInterval(out_name);
	   l_cartesian->printInteriorVelocity(out_name, binary);
	   
	   /*
	   l_cartesian->printInteriorVelocity_vtu(out_name, binary);
	   if( (pp_numnodes()) > 1 && (pp_mynode() == 1))
	       l_cartesian->outputParallelVisitFile(out_name, binary);
	   */
	   l_cartesian->printExpandedMesh(out_name, binary);
	    
	    if (debugging("trace")) 
		(void) printf("Calling initMovieVariables()\n");
	    l_cartesian->initMovieVariables();
	    if (debugging("trace")) 
		(void) printf("Calling FT_AddMovieFrame()\n");
	    FT_AddMovieFrame(front,out_name,binary);
	    //FT_AddMovieFrame_vtp(front,out_name,binary);

	    if (majority.doCalculation())
	    {
		majority.Init();
		start_clock("Getting Raw Data");
		majority.getRawData();
		stop_clock("Getting Raw Data");
		start_clock("Getting Histogram");
		majority.getPhyCompHistogram(out_name);
		stop_clock("Getting Histogram");
		majority.ClearAllCorner();
	    }

	    if (droplet.doCalculation())
	    {

		//droplet.InitFromInputFile();
		//droplet.InitFromInputFile_binary();
		droplet.preInit();
		droplet.Init();
		//droplet.InitFromInputFile_noMalloc();

		droplet.GetCompIndexForTris();
		printf("\n First step finished in droplet analysis\n");


		droplet.GetLocalVolume();
		printf("\n Second step finished in droplet analysis\n");
	/*	
		droplet.outputAllTris();
		droplet.outputAllCells();
		droplet.outputAllTris_binary();
		droplet.outputAllCells_binary();
*/
		//droplet.printAllTris();
		//droplet.printAllCells();

		droplet.FillVolumeForBlankCells();
		printf("\n Third step finished in droplet analysis\n");

		//droplet.printAllCells();

		droplet.outputHistogram(out_name);
		printf("\n After output histogram\n");
		droplet.ClearAll();
		droplet.clearTempFront();
		printf("\n After droplet clear \n");
	    }
	    
       }

       //clean_up(0);

	/******************************************************************/

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));
	}

	if (debugging("step_size"))
                printf("Time step from start: %f\n",front->dt);


        front->restart_step = RestartStep;
        front->start_small_dt = NO;
        front->end_small_dt = NO;
        for (;;)
        {
            /* Propagating interface for time step dt */

	    start_clock("Main_Loop");
	    if (debugging("trace"))
                printf("Before FT_Propagate()\n");

	    ((Incompress_Solver_Smooth_3D_Cylindrical*)l_cartesian)->adjustReflectionBuffer();

	    start_clock("FT_Propagate");
            FT_Propagate(front);
	    stop_clock("FT_Propagate");
            if(front->restart_small_dt == YES)
            {
                fprintf(stdout, "use of small CFL number from the time step: %d\n",front->last_restart_step);
                front->start_small_dt = YES;
                front->end_small_dt = NO;
           
                restart_state_name      = f_basic->restart_state_name;
                restart_name            = f_basic->restart_name;
                sprintf(restart_state_name,"%s/restart/state.t%s",f_basic->out_name,
                                right_flush(front->last_restart_step,7));
                sprintf(restart_name,"%s/restart/intfc-t%s",f_basic->out_name,
                                right_flush(front->last_restart_step,7));
                if (pp_numnodes() > 1)
                {
                    sprintf(restart_name,"%s-p%s",restart_name,
                                right_flush(pp_mynode(),4));
                    sprintf(restart_state_name,"%s-p%s",restart_state_name,
                                right_flush(pp_mynode(),4));
                }
                FT_ReadSpaceDomain(in_name,f_basic);
                FT_StartUp(front,f_basic);
                FT_InitDebug(in_name);
                l_cartesian->readFrontInteriorStates(restart_state_name, binary,RegridRestart);
                start_clock("FT_Propagate");
                FT_Propagate(front);
                stop_clock("FT_Propagate");
                front->restart_small_dt = NO;
            }

	    if (debugging("trace")) printf("Passed FT_Propagate()\n");

            if (debugging("trace")) printf("Calling ifluid solve()\n");
	    start_clock("Fluid_Solve");
	    l_cartesian->solve(front->dt);
	    stop_clock("Fluid_Solve");
	    if (debugging("trace")) printf("Passed ifluid solve()\n");
	    if (debugging("trace"))
            {
                (void) printf("After solve()\n");
                (void) print_storage("at end of time step","trace");
            }

	    FT_AddTimeStepToCounter(front);
            
            if(front->step == 0)
                front->time = 0.0;
 	
            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
	    //if (debugging("step_size"))
            //    (void) printf("Time step from FT_SetTimeStep(): %20.14f\n",
	    //				front->dt);
            front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
            if(front->start_small_dt == YES && front->end_small_dt == NO)
                front->dt = front->dt/9.0;

	    if (debugging("step_size"))
                (void) printf("Time step from l_cartesian->max_dt(): %20.14f\n",
					front->dt);
	
	    start_clock("Output_Section");
            /* Output section */

            (void) printf("\ntime = %20.14f   step = %5d   next dt = %20.14f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            printf("\n Entering the Step Diagnosis part \n");

            l_cartesian->multiDiagnosisStep(out_name);

	    if (majority.doCalculation())
	    {
		majority.Init();
		start_clock("Getting Raw Data");
		majority.getRawData();
		stop_clock("Getting Raw Data");

		start_clock("Getting Histogram");
		majority.getPhyCompHistogram(out_name);
		stop_clock("Getting Histogram");
		majority.ClearAllCorner();
	    }

	    //if (droplet.doCalculation() && ( (front->step == 131770) ) )
	    if (droplet.doCalculation())
	    {
		droplet.preInit();
		droplet.Init();
		droplet.GetCompIndexForTris();
		printf("\n First step finished in droplet analysis\n");



		droplet.GetLocalVolume();
		printf("\n Second step finished in droplet analysis\n");

	/*	
		droplet.outputAllTris();
		droplet.outputAllCells();

		droplet.outputAllTris_binary();
		droplet.outputAllCells_binary();
	*/
		
/*
		droplet.printAllTris();
		droplet.printAllCells();
*/
		droplet.FillVolumeForBlankCells();
		printf("\n Third step finished in droplet analysis\n");

		//droplet.printAllCells();
		//droplet.outputHistogram(out_name);
		printf("\n After output histogram\n");
		droplet.ClearAll();
		droplet.clearTempFront();
		printf("\n After droplet clear \n");
	    }

            printf("\n Passed the Step Diagnosis part \n");

            if (RegridRun || front->repartition == YES)
            {
                FT_Save(front,out_name);
                l_cartesian->printFrontInteriorStatesRegridRep(out_name, binary, RegridRun);
            }
            front->regrid_restart == NO;
            if (RegridRun || front->repartition == YES)
                clean_up(0);

            if (FT_IsSaveTime(front) || front->step == RestartStep+1)
	    {
                front->end_small_dt = YES;
                front->last_restart_step = front->step;
                fprintf(stdout, "last_restart_step %d\n",front->last_restart_step);
            	FT_Save(front,out_name);
		l_cartesian->printFrontInteriorStates(out_name, binary);
	    }

	    //if (debugging("trace"))
            //    printf("After print output()\n");
            if (FT_IsMovieFrameTime(front))
	    {
		l_cartesian->multiDiagnosisInterval(out_name);
		l_cartesian->printInteriorVelocity(out_name, binary);
		/*
		l_cartesian->printInteriorVelocity_vtu(out_name, binary);
		if( (pp_numnodes()) > 1 && (pp_mynode() == 1))
		    l_cartesian->outputParallelVisitFile(out_name, binary);
		*/
	    	if (debugging("trace")) 
		    (void) printf("Calling initMovieVariables()\n");
	        l_cartesian->initMovieVariables();
	    	if (debugging("trace")) 
		    (void) printf("Calling FT_AddMovieFrame()\n");
            	FT_AddMovieFrame(front,out_name,binary);
            	//FT_AddMovieFrame_vtp(front,out_name,binary);
	    }

            if (FT_TimeLimitReached(front))
                    break;

	    if (debugging("storage"))
	    {
		char s[100];
		sprintf(s,"Storage at end of time step %d",front->step);
		print_storage(s,"storage");
	    }
	    FT_TimeControlFilter(front);
	    if (debugging("step_size"))
                (void) printf("Time step from FT_TimeControlFilter(): %f\n",
                                        front->dt);
	stop_clock("Output_Section");
	stop_clock("Main_Loop");
        }
	if (debugging("trace")) printf("After time loop\n");
}       /* end ifluid_driver */


// The testing mode using exact solution

static  void ifluid_driver_debug(
        Front *front,
	Incompress_Solver_Smooth_Basis *l_cartesian)
{
        double CFL;
	int i,dim = front->rect_grid->dim;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

	if (RestartRun)
	{
	    FT_ParallelExchIntfcBuffer(front);
	}
	else
	{
	    FT_RedistMesh(front);
	}

        if (!RestartRun)
        {
	    FT_ResetTime(front);
            FT_SetOutputCounter(front);
	    if (debugging("trace"))
		printf("Before FrontProp() front->dt = %f\n",front->dt);
            FT_Propagate(front);
	    if (debugging("trace")) printf("Calling ifluid solve()\n");
            l_cartesian->solve(front->dt);
	    if (debugging("trace")) printf("Passed ifluid solve()\n");
	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
        }
        else
        {
	    FT_SetOutputCounter(front);
        }

	FT_TimeControlFilter(front);

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));
	}

	if (debugging("step_size"))
                printf("Time step from start: %f\n",front->dt);

	double dh;
	dh = std::min(front->rect_grid->h[0], front->rect_grid->h[1]);
	if (dim == 3)
	    dh = std::min(dh, front->rect_grid->h[2]);

	double dt = CFL * dh;
	int nStep = int (front->max_time/dt) + 1;
	dt = front->max_time/ nStep;
	front->dt = dt;

	if(dim==2)
	{
	    ((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,false);
	    //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,false);
	    //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
	}
	else if(dim==3)
	{
	    if (front->coordinate == 'R' || front->coordinate == 'r')
	    {
	    	((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,false);
	    	//((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,false);
	    	//((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
	    }
	    else if (front->coordinate == 'C' || front->coordinate == 'c')
	    {
	    	((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,true);
	    	//((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,true);
	    	//((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,true);
	    }
	}

	for(int i=1; i<=nStep; i++)
	{
	    front->step = i;
	    front->dt = dt;

	    l_cartesian->solve(front->dt);

	    front->time += front->dt;
	    printf("\ntime = %f   step = %5d   next dt = %f\n",
		    front->time,front->step,front->dt);
	    if(dim==2)
	    {
		((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,false);
		//((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,false);
		//((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
	    }
	    else if(dim==3)
	    {
		if (front->coordinate == 'R' || front->coordinate == 'r')
		{
		    ((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,false);
		    //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,false);
		    //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
		}
		else
		{
		((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,true);
		}
	    }
	}

	if (debugging("trace")) printf("After time loop\n");
}       /* end ifluid_driver */

static int l_cartesian_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	double *coords = Coords(p);
	((Incompress_Solver_Smooth_Basis*)params)->getVelocity(coords, vel);
}	/* end l_cartesian_vel */

void initIntfc_oneComponent(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack)
{
    level_func_pack->neg_component = LIQUID_COMP1;
    level_func_pack->pos_component = LIQUID_COMP2;
    if(front->rect_grid->dim==2)
	level_func_pack->func = level_oneComponent_2D;
    else if(front->rect_grid->dim==3)
    {
	if (front->coordinate == 'R' || front->coordinate == 'r')
	    level_func_pack->func = level_oneComponent_3D;
	else if(front->coordinate == 'C' || front->coordinate == 'c')
	{
	    printf("\nSetting the level_func_pack to cylindrical case!!\n");
	    level_func_pack->func = level_oneComponent_3D_Cylindrical;
	}
    }
    level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
}	/* end initCirclePlaneIntfc */

double level_oneComponent_2D(
        POINTER func_params,
        double *coords)
{
    double dist = sqrt(
	    sqr(coords[0]-0.5) +
	    sqr(coords[1]-0.5)) - 0.2;
    return dist;
}       /* end level_circle_func */

double level_oneComponent_3D(
        POINTER func_params,
        double *coords)
{
    double dist = sqrt(
	    sqr(coords[0]-0.5) +
	    sqr(coords[1]-0.5) +
	    sqr(coords[2]-0.5)  ) - 0.2;
    return dist;
}       /* end level_circle_func */

double level_oneComponent_3D_Cylindrical(
        POINTER func_params,
        double *coords)
{
    double dist = sqrt(
	    sqr(coords[0]-3) +
	    sqr(coords[1]-0) +
	    sqr(coords[2]-1.5)  ) - 0.3;
    return dist;
//    return 1;
}       /* end level_circle_func */
