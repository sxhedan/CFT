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
*				example.c:
*
*	3D surface mesh redistribution example for Front Package
*	The triangle size and quality control is only approximate
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static double sphere_func(POINTER,double*);
static int test_curvature_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);


char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
        double center[3];
	double radius;
} TEST_SPHERE_PARAMS;


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	TEST_SPHERE_PARAMS s_params;
	static LEVEL_FUNC_PACK level_func_pack;
	char dname[100];
	int i,count;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0; 	f_basic.L[2] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0; 	f_basic.U[2] = 1.0;
	f_basic.gmax[0] = 32;	f_basic.gmax[1] = 32; f_basic.gmax[2] = 32;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s.ts%s",restart_name,right_flush(RestartStep,7));
#if defined(__MPI__)
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

	FT_StartUp(&front,&f_basic);

	if (!RestartRun)
	{
	    /* Initialize interface through level function */
	    s_params.center[0] = 0.5;
	    s_params.center[1] = 0.5;
	    s_params.center[2] = 0.5;
	    s_params.radius = 0.4;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&s_params;
	    level_func_pack.func = sphere_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    FT_InitIntfc(&front,&level_func_pack);
	}

	/* Original mesh produced my marching cube method */
	printf("Original mesh produced my marching cube method\n");
	printf("NumOfIntfcTris = %d\n",NumOfIntfcTris(front.interf));
	sprintf(dname,"intfc-%d",count++);
	gview_plot_interface(dname,front.interf);

	/* Default redistribute */
	/* Default values
		2.0,	 Max tri area factor of standard area
		0.5,	 Min tri area factor of standard area
		15,	 Min angle at tri vertex
		1.3.	 Max scaled tri side length
	*/
	sprintf(dname,"intfc-%d",count++);
	FT_RedistMesh(&front);
	printf("After redistribution using default parameters\n");
	printf("NumOfIntfcTris = %d\n",NumOfIntfcTris(front.interf));
	printf("Max_tri_sqr_area = %24.18g\n",Max_tri_sqr_area(&front,0));
	printf("Min_tri_sqr_area = %24.18g\n",Min_tri_sqr_area(&front,0));
	printf("Max_scaled_tri_side_sqr_length = %24.18g\n",
			Max_scaled_tri_side_sqr_length(&front));
	printf("Aspect_ratio_tolerance = %24.18g\n",
			Aspect_ratio_tolerance(&front,0));
	gview_plot_interface(dname,front.interf);
	if (pp_numnodes() > 1)
	{
	    clean_up(0);
	}

	FrontSetTriParams(&front,
		1.0,	// Max tri area factor of standard area
		0.25,	// Min tri area factor of standard area
		30,	// Min angle at tri vertex
		1.2);	// Max scaled tri side length
	sprintf(dname,"intfc-%d",count++);
	FT_RedistMesh(&front);
	printf("After redistribution\n");
	printf("NumOfIntfcTris = %d\n",NumOfIntfcTris(front.interf));
	printf("Max_tri_sqr_area = %24.18g\n",Max_tri_sqr_area(&front,0));
	printf("Min_tri_sqr_area = %24.18g\n",Min_tri_sqr_area(&front,0));
	printf("Max_scaled_tri_side_sqr_length = %24.18g\n",
			Max_scaled_tri_side_sqr_length(&front));
	printf("Aspect_ratio_tolerance = %24.18g\n",
			Aspect_ratio_tolerance(&front,0));
	gview_plot_interface(dname,front.interf);

	FrontSetTriParams(&front,4.0,1.0,30,1.2);
	sprintf(dname,"intfc-%d",count++);
	for (i = 0; i < 3; ++i)
	    FT_RedistMesh(&front);
	printf("After redistribution\n");
	printf("NumOfIntfcTris = %d\n",NumOfIntfcTris(front.interf));
	printf("Max_tri_sqr_area = %24.18g\n",Max_tri_sqr_area(&front,0));
	printf("Min_tri_sqr_area = %24.18g\n",Min_tri_sqr_area(&front,0));
	printf("Max_scaled_tri_side_sqr_length = %24.18g\n",
			Max_scaled_tri_side_sqr_length(&front));
	printf("Aspect_ratio_tolerance = %24.18g\n",
			Aspect_ratio_tolerance(&front,0));
	gview_plot_interface(dname,front.interf);

	clean_up(0);
}

/********************************************************************
 *	Sample (sphere 3D) level function for the initial interface    *
 ********************************************************************/

static double sphere_func(
        POINTER func_params,
        double *coords)
{
        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
	double x0,y0,z0,R;
	double distance;

        x0 = s_params->center[0];
        y0 = s_params->center[1];
        z0 = s_params->center[2];
	R = s_params->radius;

	distance = sqrt(sqr(coords[0] - x0) + sqr(coords[1] - y0) +
			sqr(coords[2] - z0)) - R;

        return distance;

}       /* end sphere_func */

