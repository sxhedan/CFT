#include <cFluid.h>

extern void read_cFluid_params(
	char *inname,
	EQN_PARAMS *eqn_params)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	int i,dim = eqn_params->dim;

	eqn_params->prob_type = ERROR_TYPE;
	eqn_params->tracked = YES;		// Default
	CursorAfterString(infile,"Enter problem type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'T' || string[0] == 't')
	{
	    if (string[10] == 'B' || string[10] == 'b')
	    	eqn_params->prob_type = TWO_FLUID_BUBBLE;
	    else if (string[10] == 'R' || string[10] == 'r')
	    {
		if (string[11] == 'T' || string[11] == 't')
	    	    eqn_params->prob_type = TWO_FLUID_RT;
		else if (string[11] == 'M' || string[11] == 'm')
	    	    eqn_params->prob_type = TWO_FLUID_RM;
	    }
	    else if (string[10] == 'V' || string[10] == 'v')
	    	eqn_params->prob_type = TWO_FLUID_VST_RM;
	} 
	else if (string[0] == 'F' || string[0] == 'f')
	    eqn_params->prob_type = FLUID_SOLID_CIRCLE;
	else if (string[0] == 'B' || string[0] == 'b')
	    eqn_params->prob_type = BUBBLE_SURFACE;
	else if (string[0] == 'I' || string[0] == 'i')
	    eqn_params->prob_type = IMPLOSION;
	else if (string[0] == 'P' || string[0] == 'p')
	    eqn_params->prob_type = PROJECTILE;
	else if (string[0] == 'R' || string[0] == 'r')
	    eqn_params->prob_type = RIEMANN_PROB;
	else if (string[0] == 'M' || string[0] == 'm')
	    eqn_params->prob_type = MT_FUSION;
	else if (string[0] == 'O' || string[0] == 'o')
	{
	    if (string[5] == 'S' || string[5] == 's')
	    	eqn_params->prob_type = ONED_SSINE;
	    else if (string[5] == 'B' || string[5] == 'b')
	    	eqn_params->prob_type = ONED_BLAST;
	}
	else if (string[0] == 'S' || string[0] == 's')
	{
	    if (string[4] == 'O' || string[4] == 'o')
		eqn_params->prob_type = SOD_OBLIQ;
	    else if (string[4] == '3')
		eqn_params->prob_type = SOD_3D;
	}
	else if (string[0] == 'C' || string[0] == 'c')
	    eqn_params->prob_type = CFT_TEST;
	CursorAfterString(infile,"Enter numerical scheme for interior solver:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'T':
	case 't':
	    switch (string[4])
	    {
	    case '1':
		eqn_params->num_scheme = TVD_FIRST_ORDER;
		break;
	    case '2':
		eqn_params->num_scheme = TVD_SECOND_ORDER;
		break;
	    case '4':
		eqn_params->num_scheme = TVD_FOURTH_ORDER;
		break;
	    default:
		printf("Numerical scheme %s not implemented!\n",string);
		clean_up(ERROR);
	    }
	    break;
	case 'W':
	case 'w':
	    switch (string[5])
	    {
        case 'S':
        eqn_params->num_scheme = WENO_STRANG_SPLIT;
        break;
	    case '1':
		eqn_params->num_scheme = WENO_FIRST_ORDER;
		break;
	    case '2':
		eqn_params->num_scheme = WENO_SECOND_ORDER;
		break;
	    case '4':
		eqn_params->num_scheme = WENO_FOURTH_ORDER;
		break;
	    default:
		printf("Numerical scheme %s not implemented!\n",string);
		clean_up(ERROR);
	    }
	    break;
	default:
	    printf("Numerical scheme %s not implemented!\n",string);
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter order of point propagator:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
    case 'S':
        eqn_params->point_prop_scheme = STRANG_SPLIT;
        break;
	case '1':
	    eqn_params->point_prop_scheme = FIRST_ORDER;
	    break;
	case '2':
	    eqn_params->point_prop_scheme = SECOND_ORDER;
	    break;
	case '4':
	    eqn_params->point_prop_scheme = FOURTH_ORDER;
	    break;
	default:
	    printf("Point propagator order %s not implemented!\n",string);
	    clean_up(ERROR);
	}

	assert(eqn_params->prob_type != ERROR_TYPE);
	fclose(infile);

}	/* end read_cFluid_params */

extern void read_movie_options(
	char *inname,
	EQN_PARAMS *eqn_params)
{
	static MOVIE_OPTION *movie_option;
	FILE *infile = fopen(inname,"r");
	char string[100];

	FT_ScalarMemoryAlloc((POINTER*)&movie_option,sizeof(MOVIE_OPTION));
	eqn_params->movie_option = movie_option;
	movie_option->set_bounds = NO;	// default
	if (CursorAfterStringOpt(infile,"Type y to set movie bounds:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->set_bounds = YES;
	}

	CursorAfterString(infile,"Type y to make movie of density:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	{
	    movie_option->plot_dens = YES;
	    if (movie_option->set_bounds)
	    {
		CursorAfterString(infile,"Enter min and max density:");
		fscanf(infile,"%lf %lf",&movie_option->min_dens,
				&movie_option->max_dens);
		(void) printf("%f %f\n",movie_option->min_dens,
				movie_option->max_dens);
	    }
	}
	CursorAfterString(infile,"Type y to make movie of pressure:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	{
	    movie_option->plot_pres = YES;
	    if (movie_option->set_bounds)
	    {
		CursorAfterString(infile,"Enter min and max pressure:");
		fscanf(infile,"%lf %lf",&movie_option->min_pres,
				&movie_option->max_pres);
		(void) printf("%f %f\n",movie_option->min_pres,
				movie_option->max_pres);
	    }
	}
	if (eqn_params->dim != 1)
	{
	    CursorAfterString(infile,"Type y to make movie of vorticity:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	movie_option->plot_vort = YES;
	    }
	}
	CursorAfterString(infile,"Type y to make movie of velocity:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	{
	    movie_option->plot_velo = YES;
	    if (movie_option->set_bounds)
	    {
		CursorAfterString(infile,"Enter min and max velocity:");
		fscanf(infile,"%lf %lf",&movie_option->min_velo,
				&movie_option->max_velo);
		(void) printf("%f %f\n",movie_option->min_velo,
				movie_option->max_velo);
	    }
	}

	if (eqn_params->dim == 3)
	{
	    CursorAfterString(infile,"Type y to make yz cross section movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_cross_section[0] = YES;
	    CursorAfterString(infile,"Type y to make xz cross section movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_cross_section[1] = YES;
	    CursorAfterString(infile,"Type y to make xy cross section movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_cross_section[2] = YES;
	}
	fclose(infile);
}	/* end read_movie_options */

extern void read_statistics_options(
	char		*inname,
	EQN_PARAMS	*eqn_params)
{
	static Intfc_extrema *iext;
	FILE *infile = fopen(inname,"r");
	char string[100];

	FT_ScalarMemoryAlloc((POINTER*)&iext,sizeof(Intfc_extrema));
	eqn_params->iext = iext;

	CursorAfterString(infile,"Enter a sub-grid refinement factor:");
	fscanf(infile,"%lf",&iext->rfactor);
	(void) printf("%lf\n",iext->rfactor);

	CursorAfterString(infile,"Type y to get data for 1%-99% levels:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	    iext->do_01 = YES;
	else
	    iext->do_01 = NO;
	CursorAfterString(infile,"Type y to get data for 5%-95% levels:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	    iext->do_05 = YES;
	else
	    iext->do_05 = NO;

	fclose(infile);
}
