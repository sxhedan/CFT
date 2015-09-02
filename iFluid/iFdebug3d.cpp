/**
* iFdebug3d.cpp
*
*  Created on: Jun 23, 2010
*      Author: shuqiangwang & YijieZhou
*
*  Testin the algorithm for solving incompressible N-S equation using exact solution in 3D
*/
#include "iFluid_debug.h"
#include "ifluid_basic.h"
#include "solver.h"

/*****************************************************************************
*			Incompress_Solver_Smooth_3D_Cartesian_Debug
* 3D debug.
*****************************************************************************/


Incompress_Solver_Smooth_3D_Cartesian_Debug::Incompress_Solver_Smooth_3D_Cartesian_Debug(Front &front)
: Incompress_Solver_Smooth_3D_Cartesian(front)
{
}

/**
* see also:
* Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition
*/
void Incompress_Solver_Smooth_3D_Cartesian_Debug::setInitialCondition(void)
{
    //	Incompress_Solver_Smooth_3D_Cartesian_Debug::setInitialCondition();
    printf("\nEnter 3D_Debug::setInitialCondition()\n");
    int i,j,k,l,index;
    COMPONENT comp;
    double coords[MAXD];

    FT_MakeGridIntfc(front);
    setDomain();

    m_rho[0] = 1;
    m_rho[1] = 1;
    m_mu[0] = 1;
    m_mu[1] = 1;
    m_comp[0] = iFparams->m_comp1;
    m_comp[1] = iFparams->m_comp2;
    m_smoothing_radius = 0;
    m_sigma = 0;
    mu_min = rho_min = HUGE;
    for (i = 0; i < 3; ++i)
    {
	if (ifluid_comp(m_comp[i]))
	{
	    mu_min = std::min(mu_min,m_mu[i]);
	    rho_min = std::min(rho_min,m_rho[i]);
	}
    }

    double CFL = 0.1;
    front->max_time = 0.2;
    double dh;
    dh = std::min(front->rect_grid->h[0], front->rect_grid->h[1]);
    dh = std::min(dh, front->rect_grid->h[2]);
    double dt = CFL * dh;
    int nStep = int (front->max_time/dt) + 1;
    dt = front->max_time/nStep;
    front->dt = dt;
    // Initialize state at cell_center
    for (i = 0; i < cell_center.size(); i++)
    {
	getRectangleCenter(i, coords);
	cell_center[i].m_state.setZero();
	getExactSolution(coords,front->time,cell_center[i].m_state);
	cell_center[i].m_state.m_q = cell_center[i].m_state.m_P;
    }
    front->dt = 0.0;

    for (l = 0; l < dim; l++)
    {
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_U[l];
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_U[l] = array[index];
	}
    }

    // Using initial pressure as p - paverage
    double averagep, averagep_test;
    double cellVolume = top_h[0]*top_h[1]*top_h[2];
    double totalVolume;
    averagep = 0.0;
    averagep_test = 0.0;
    totalVolume = 0.0;
    printf("\nThe cell volume is : %.16g\n",cellVolume);

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);
	totalVolume += cellVolume;
	averagep += cell_center[index].m_state.m_P * cellVolume;
    }

    pp_global_sum(&totalVolume, 1);
    pp_global_sum(&averagep, 1);


    printf("\nThe total volume is : %.16g\n",totalVolume);

    averagep = averagep/totalVolume;

    printf("\nThe integral of pressure is %.16g\n",averagep);

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);
	array[index] = cell_center[index].m_state.m_P - averagep;
	averagep_test += array[index] * cellVolume;
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
    //averagep_test = averagep_test / totalVolume;

    printf("\nAfter adjusting, the integral of pressure right now is %.16g\n", averagep_test);

    computeGradientQ();
    copyMeshStates();
    setAdvectionDt();
}

void Incompress_Solver_Smooth_3D_Cartesian_Debug::solve(double dt)
{

      printf("\nEntering Incompress Solve 3D!! The dt for this step is %.16g\n",dt);
//    Incompress_Solver_Smooth_3D_Cartesian::solve(dt);

//    solve_diffusionOnly(dt);
//    return;


    m_t_old = front->time;
    m_t_int = front->time + front->dt/2;
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
    setSmoProOnePhase();
    stop_clock("setSmoothedProperties");
    if (debugging("trace"))
	printf("Passed setSmoothedProperties()\n");

//    // 1) solve for intermediate velocity
//    start_clock("computeAdvection");
//    computeAdvection();
//    stop_clock("computeAdvection");
//    if (debugging("trace"))
//	printf("max_speed after computeAdvection(): %20.14f\n",
//		max_speed);

    start_clock("compAdvectionTerm");
    compAdvectionTerm_decoupled();
    stop_clock("compAdvectionTerm");

    start_clock("compDiffWithSmoothProperty");
    compDiffWithSmoothProperty_2nd_decoupled();
//    compDiffWithSmoothProperty();
    stop_clock("compDiffWithSmoothProperty");

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
    {
	sampleVelocity();
    }

    if (debugging("trace"))
	printf("max_speed after computeNewVelocity(): %20.14f\n",
		max_speed);

    start_clock("copyMeshStates");
    copyMeshStates();
    stop_clock("copyMeshStates");

    setAdvectionDt();
    stop_clock("solve");
}

void Incompress_Solver_Smooth_3D_Cartesian_Debug::solve_diffusionOnly(double dt)
{
    m_t_old = front->time;
    m_t_int = front->time + front->dt/2;
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


    start_clock("compDiffWithSmoothProperty");
    compDiffWithSmoothProperty_2nd_decoupled();
//    compDiffWithSmoothProperty();
    stop_clock("compDiffWithSmoothProperty");


    start_clock("copyMeshStates");
    copyMeshStates();
    stop_clock("copyMeshStates");

    setAdvectionDt();
    stop_clock("solve");
}

void Incompress_Solver_Smooth_3D_Cartesian_Debug::getVelocity(double *p, double *U)
{
    int i;
    for (i = 0; i < dim; ++i)
	U[i] = 0.0;
}

/**
* This function is written here so that it can be overwritten by a derived
* class.
*/
boolean Incompress_Solver_Smooth_3D_Cartesian_Debug::FT_StateStructAtGridCrossing_tmp(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	COMPONENT comp,
	Locstate *state,
	HYPER_SURF **hs,
	double *crx_coords,
	double t)
{
    boolean bRet = FT_StateStructAtGridCrossing(front,icoords,dir,comp,state,hs,crx_coords);

    if(bRet)
    {
	static int debug_index = 0;
	L_STATE exact_state;

	getExactSolution(crx_coords,t,exact_state);

	STATE *fstate = (STATE*)(*state);
	fstate->vel[0] = exact_state.m_U[0];
	fstate->vel[1] = exact_state.m_U[1];
	fstate->vel[2] = exact_state.m_U[2];
    }
    return bRet;
}


/**
* see also:
* EBM2D_Liquid_Solution::getExactSolution, using TEST_BROWN.
*/

void Incompress_Solver_Smooth_3D_Cartesian_Debug::computeSourceTerm_Adv(
	double *coords,
	L_STATE &state)
{
    //    Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm(coords,state);

    double t = m_t_old;
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    double mu = m_mu[0];
    double a = 1;
    double d = 1;

    // incompressible
    state.m_U[0] = -a * d * d * (-exp(-d * d * t + a * x) * sin(a * y + d * z) - exp(-d * d * t + a * z) * cos(a * x + d * y) + mu * exp(-d * d * t + a * x) * sin(a * y + d * z) + mu * exp(-d * d * t + a * z) * cos(a * x + d * y));
    state.m_U[1] = -a * d * d * (-exp(-d * d * t + a * y) * sin(a * z + d * x) - exp(-d * d * t + a * x) * cos(a * y + d * z) + mu * exp(-d * d * t + a * y) * sin(a * z + d * x) + mu * exp(-d * d * t + a * x) * cos(a * y + d * z));
    state.m_U[2] = -a * d * d * (-exp(-d * d * t + a * z) * sin(a * x + d * y) - exp(-d * d * t + a * y) * cos(a * z + d * x) + mu * exp(-d * d * t + a * z) * sin(a * x + d * y) + mu * exp(-d * d * t + a * y) * cos(a * z + d * x));

//    // diffusion
//        state.m_U[0] = 0.4e1 / 0.125e3 * exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.5e1) * (double) (5 * x * x + 5 * y * y + 5 * z * z - 25 * t - 25 + 30 * mu * t + 30 * mu - 4 * mu * x * x - 4 * mu * y * y - 4 * mu * z * z) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1;
//        state.m_U[1] = exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.2e1) * (double) (x * x + y * y + z * z - 2 * t - 2 + 6 * mu * t + 6 * mu - 2 * mu * x * x - 2 * mu * y * y - 2 * mu * z * z) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 / 0.4e1;
//    //    state.m_U[2] = exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.3e1) * (double) (3 * x * x + 3 * y * y + 3 * z * z - 9 * t - 9 + 18 * mu * t + 18 * mu - 4 * mu * x * x - 4 * mu * y * y - 4 * mu * z * z) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 / 0.27e2;
//        state.m_U[2] = state.m_U[1];
}

void Incompress_Solver_Smooth_3D_Cartesian_Debug::computeSourceTerm(
	double *coords,
	L_STATE &state)
{
    //    Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm(coords,state);

//    double t = front->time + front->dt/2;
    double t = m_t_int;
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    double mu = m_mu[0];
    double a = 1;
    double d = 1;

    // incompressible
    state.m_U[0] = -a * d * d * (-exp(-d * d * t + a * x) * sin(a * y + d * z) - exp(-d * d * t + a * z) * cos(a * x + d * y) + mu * exp(-d * d * t + a * x) * sin(a * y + d * z) + mu * exp(-d * d * t + a * z) * cos(a * x + d * y));
    state.m_U[1] = -a * d * d * (-exp(-d * d * t + a * y) * sin(a * z + d * x) - exp(-d * d * t + a * x) * cos(a * y + d * z) + mu * exp(-d * d * t + a * y) * sin(a * z + d * x) + mu * exp(-d * d * t + a * x) * cos(a * y + d * z));
    state.m_U[2] = -a * d * d * (-exp(-d * d * t + a * z) * sin(a * x + d * y) - exp(-d * d * t + a * y) * cos(a * z + d * x) + mu * exp(-d * d * t + a * z) * sin(a * x + d * y) + mu * exp(-d * d * t + a * y) * cos(a * z + d * x));

//    // diffusion
//        state.m_U[0] = 0.4e1 / 0.125e3 * exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.5e1) * (double) (5 * x * x + 5 * y * y + 5 * z * z - 25 * t - 25 + 30 * mu * t + 30 * mu - 4 * mu * x * x - 4 * mu * y * y - 4 * mu * z * z) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1;
//        state.m_U[1] = exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.2e1) * (double) (x * x + y * y + z * z - 2 * t - 2 + 6 * mu * t + 6 * mu - 2 * mu * x * x - 2 * mu * y * y - 2 * mu * z * z) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 / 0.4e1;
//    //    state.m_U[2] = exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.3e1) * (double) (3 * x * x + 3 * y * y + 3 * z * z - 9 * t - 9 + 18 * mu * t + 18 * mu - 4 * mu * x * x - 4 * mu * y * y - 4 * mu * z * z) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 / 0.27e2;
//        state.m_U[2] = state.m_U[1];
}

/**
* TEST_EthierSteinman in EBM3D_Liquid_Solution.
*/
void Incompress_Solver_Smooth_3D_Cartesian_Debug::getExactSolution(
	double *coords,
	double t,
	L_STATE &state)
{
    //    double pi = 0.3141592654e1;
    //    double omega = 1+sin(2*pi*t*t);
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    double mu = m_mu[0];
    double a = 1;
    double d = 1;
    double t_int = t - (front->dt)/2.0;

    // incompressible
    state.m_U[0] = -a * (exp(a * x) * sin(a * y + d * z) + exp(a * z) * cos(a * x + d * y)) * exp(-d * d * t);
    state.m_U[1] = -a * (exp(a * y) * sin(a * z + d * x) + exp(a * x) * cos(a * y + d * z)) * exp(-d * d * t);
    state.m_U[2] = -a * (exp(a * z) * sin(a * x + d * y) + exp(a * y) * cos(a * z + d * x)) * exp(-d * d * t);
    state.m_P = -(double) (a * a) * (exp((double) (2 * a * x)) + exp((double) (2 * a * y)) + exp((double) (2 * a * z)) + 0.2e1 * sin((double) (a * x + d * y)) * cos((double) (a * z + d * x)) * exp((double) (a * (y + z))) + 0.2e1 * sin((double) (a * y + d * z)) * cos((double) (a * x + d * y)) * exp((double) (a * (z + x))) + 0.2e1 * sin((double) (a * z + d * x)) * cos((double) (a * y + d * z)) * exp((double) (a * (x + y)))) * exp(-(double) (2 * d * d * t_int)) / 0.2e1;

    // diffusion
//    state.m_U[0] = 0.4e1 / 0.5e1 * exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.5e1) / 0.3141592654e1 / (double) (t + 1);
//    state.m_U[1] = exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.2e1) / 0.3141592654e1 / (double) (t + 1) / 0.2e1;
////    state.m_U[2] = exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.3e1) / 0.3141592654e1 / (double) (t + 1) / 0.3e1;
//    state.m_U[2] = state.m_U[1];

}

void Incompress_Solver_Smooth_3D_Cartesian_Debug::saveStates_Tecplot(
	const char*out_name,
	double t,
	bool bPrintCoords,
	bool bPrintExact,
	bool bPrintError)
{
    bPrintExact = true;
    /*
    char filename[200];
    sprintf(filename,"%s/states_%s.plt",out_name,
	    right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, Z, U, V, W, P");
    if(bPrintExact)
	fprintf(hfile, ", U_exact, V_exact, W_exact, P_exact");
    if(bPrintError)
	fprintf(hfile, ", U_error, V_error, W_error, P_error");
    fprintf(hfile, "\n");

    fprintf(hfile, "ZONE I=%d, J=%d, K=%d, F=POINT\n",
	    imax-imin+1, jmax-jmin+1, kmax-kmin+1);
*/
    int c, icoords[3];
    double *coords, cellVolume = top_h[0]*top_h[1]*top_h[2];

    int max_ij[4][4], index;
    max_ij[0][0] = max_ij[0][1] = max_ij[0][2] =-1;
    max_ij[1][0] = max_ij[1][1] = max_ij[1][2] =-1;
    max_ij[2][0] = max_ij[2][1] = max_ij[2][2] =-1;
    max_ij[3][0] = max_ij[3][1] = max_ij[3][2] =-1;

    L_STATE error, L1, L2, LInf;

    L_STATE state, state_exact;

    double average_exactp = 0.0;
    double totalVolume = 0.0;

    for (int k=kmin; k<=kmax; k++)
	for (int j=jmin; j<=jmax; j++)
	    for (int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;
		state = cell_center[index].m_state;

		getExactSolution(coords,t,state_exact);
		average_exactp = average_exactp + state_exact.m_P*cellVolume;
		totalVolume += cellVolume;
	    }
    average_exactp = average_exactp / totalVolume;
    printf("\nIn error calculation, time is %.16g, Total Volume is %.16g\n",front->time, totalVolume);
    printf("\nIn error calculation, time is %.16g, average_exactp is %.16g\n",front->time, average_exactp);

    for(int k=kmin; k<=kmax; k++)
	for(int j=jmin; j<=jmax; j++)
	    for(int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;
		state = cell_center[index].m_state;

		if(bPrintExact)
		    getExactSolution(coords, t, state_exact);
		else
		{
		    state_exact.m_U[0] = state_exact.m_U[1] =
			    state_exact.m_U[2] = state_exact.m_P = 0;
		}

		state_exact.m_P = state_exact.m_P - average_exactp;

		error.m_U[0] = fabs(state.m_U[0] - state_exact.m_U[0]);
		error.m_U[1] = fabs(state.m_U[1] - state_exact.m_U[1]);
		error.m_U[2] = fabs(state.m_U[2] - state_exact.m_U[2]);
		error.m_P    = fabs(state.m_P - state_exact.m_P);

		// LInfinity
		for(unsigned int kk=0; kk<3; kk++)
		    if(error.m_U[kk] > LInf.m_U[kk])
		    {
			LInf.m_U[kk] = error.m_U[kk];
			max_ij[kk][0] = i;
			max_ij[kk][1] = j;
			max_ij[kk][2] = k;
		    }
		if(error.m_P > LInf.m_P)
		{
		    LInf.m_P = error.m_P;
		    max_ij[3][0] = i;
		    max_ij[3][1] = j;
		    max_ij[3][2] = k;
		}

		// L1
		L1.m_U[0] += error.m_U[0]*cellVolume;
		L1.m_U[1] += error.m_U[1]*cellVolume;
		L1.m_U[2] += error.m_U[2]*cellVolume;
		L1.m_P += error.m_P*cellVolume;

		// L2
		L2.m_U[0] += error.m_U[0]*error.m_U[0]*cellVolume;
		L2.m_U[1] += error.m_U[1]*error.m_U[1]*cellVolume;
		L2.m_U[2] += error.m_U[2]*error.m_U[2]*cellVolume;
		L2.m_P += error.m_P*error.m_P*cellVolume;

		//			state_print = state;
/*
		if(bPrintCoords)
		    fprintf(hfile, "% 12.8e % 12.8e % 12.8e",
			    coords[0], coords[1], coords[2]);
		else
		    fprintf(hfile, "%3d %3d %3d",
			    i, j, k);

		fprintf(hfile, " % 12.8e % 12.8e % 12.8e % 12.8e",
			state.m_U[0], state.m_U[1], state.m_U[2], state.m_P);
		if(bPrintExact)
		    fprintf(hfile, " % 12.8e % 12.8e % 12.8e % 12.8e",
			    state_exact.m_U[0], state_exact.m_U[1], state_exact.m_U[2], state_exact.m_P);
		if(bPrintError)
		    fprintf(hfile, " % 12.8e % 12.8e % 12.8e % 12.8e",
			    state_exact.m_U[0] - state.m_U[0],
			    state_exact.m_U[1] - state.m_U[1],
			    state_exact.m_U[2] - state.m_U[2],
			    state_exact.m_P - state.m_P);
		fprintf(hfile, "\n");
		*/
	    }
    pp_global_sum(&L1.m_U[0],1);
    pp_global_sum(&L1.m_U[1],1);
    pp_global_sum(&L1.m_U[2],1);
    pp_global_sum(&L1.m_P,1);

    pp_global_sum(&L2.m_U[0],1);
    pp_global_sum(&L2.m_U[1],1);
    pp_global_sum(&L2.m_U[2],1);
    pp_global_sum(&L2.m_P,1);

    pp_global_max(&LInf.m_U[0],1);
    pp_global_max(&LInf.m_U[1],1);
    pp_global_max(&LInf.m_U[2],1);
    pp_global_max(&LInf.m_P,1);

    // LInfinity
    // L1
    // L2
    L2.m_U[0] = sqrt(L2.m_U[0]);
    L2.m_U[1] = sqrt(L2.m_U[1]);
    L2.m_U[2] = sqrt(L2.m_U[2]);
    L2.m_P = sqrt(L2.m_P);

    //fclose(hfile);

    printf("Incompress_Solver_Smooth_3D_Cartesian_Debug::saveStates_Tecplot: \n");
    printf("\t#### L1  ={%12.8f, %12.8f, %12.8f, %12.8f}\n", L1.m_U[0], L1.m_U[1], L1.m_U[2], L1.m_P);
    printf("\t#### L2  ={%12.8f, %12.8f, %12.8f, %12.8f}\n", L2.m_U[0], L2.m_U[1], L2.m_U[2], L2.m_P);
    printf("\t#### LInf={%12.8f, %12.8f, %12.8f, %12.8f}\n", LInf.m_U[0], LInf.m_U[1], LInf.m_U[2], LInf.m_P);
}

void Incompress_Solver_Smooth_3D_Cartesian_Debug::saveDivUPhi_Tecplot(
	const char*out_name,
	double t,
	bool bPrintCoords)
{
    char filename[200];
    sprintf(filename,"%s/DivUPhi_%s.plt",out_name,
	    right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, Z, DivU, Phi\n");

    fprintf(hfile, "ZONE I=%d, J=%d, K=%d, F=POINT\n",
	    imax-imin+1, jmax-jmin+1, kmax-kmin+1);

    double *coords;
    int index;
    L_STATE state;

    double max_div_U = 0, max_phi = 0;

    for(int k=kmin; k<=kmax; k++)
	for(int j=jmin; j<=jmax; j++)
	    for(int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;

		state = cell_center[index].m_state;
		if(fabs(state.div_U)>max_div_U)
		    max_div_U = fabs(state.div_U);
		if(fabs(state.m_phi)>max_phi)
		    max_phi = fabs(state.m_phi);

		if(bPrintCoords)
		    fprintf(hfile, "% 12.8e % 12.8e % 12.8e",
			    coords[0], coords[1], coords[2]);
		else
		    fprintf(hfile, "%3d %3d %3d", i, j, k);

		fprintf(hfile, " % 12.8e % 12.8e",
			state.div_U, state.m_phi);

		fprintf(hfile, "\n");
	    }

    fclose(hfile);

    printf("Incompress_Solver_Smooth_3D_Cartesian_Debug::saveDivUPhi_Tecplot: max_div_U = %12.8f, max_phi = %12.8f\n",
	    max_div_U, max_phi);
}

void Incompress_Solver_Smooth_3D_Cartesian_Debug::saveParameters_Tecplot(
	const char*out_name,
	double t,
	bool bPrintCoords)
{
    char filename[200];
    sprintf(filename,"%s/parameters_%s.plt",out_name,
	    right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, Z, rho, mu\n");

    fprintf(hfile, "ZONE I=%d, J=%d, K=%d, F=POINT\n",
	    imax-imin+1, jmax-jmin+1, kmax-kmin+1);

    double *coords;
    int index;
    L_STATE state;

    for(int k=kmin; k<=kmax; k++)
	for(int j=jmin; j<=jmax; j++)
	    for(int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;
		state = cell_center[index].m_state;
		if(bPrintCoords)
		    fprintf(hfile, "% 12.8e % 12.8e % 12.8e",
			    coords[0], coords[1], coords[2]);
		else
		    fprintf(hfile, "%3d %3d %3d", i, j, k);

		fprintf(hfile, " % 12.8e % 12.8e",
			state.m_mu, state.m_rho);
		fprintf(hfile, "\n");
	    }

    fclose(hfile);
}

/*****************************************************************************
*			Incompress_Solver_Smooth_3D_Cylindrical_Debug
* 3D debug.
*****************************************************************************/
Incompress_Solver_Smooth_3D_Cylindrical_Debug::Incompress_Solver_Smooth_3D_Cylindrical_Debug(Front &front)
: Incompress_Solver_Smooth_3D_Cylindrical(front)
{
}

/**
* see also:
* Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition
*/
void Incompress_Solver_Smooth_3D_Cylindrical_Debug::setInitialCondition(void)
{
    printf("\nUsing setInitialCondition for Cylindrical Debug class!\n");
    //	Incompress_Solver_Smooth_3D_Cylindrical_Debug::setInitialCondition();
    int i,j,k,index,l;
    COMPONENT comp;
    double coords[MAXD];

    FT_MakeGridIntfc(front);
    setDomain();

    m_rho[0] = 1;
    m_rho[1] = 1;
    m_mu[0] = 1;
    m_mu[1] = 1;
    m_comp[0] = iFparams->m_comp1;
    m_comp[1] = iFparams->m_comp2;
    m_smoothing_radius = 0;
    m_sigma = 0;
    mu_min = rho_min = HUGE;
    for (i = 0; i < 2; ++i)
    {
	if (ifluid_comp(m_comp[i]))
	{
	    mu_min = std::min(mu_min,m_mu[i]);
	    rho_min = std::min(rho_min,m_rho[i]);
	}
    }

    double CFL = 0.1;
    front->max_time = 0.2;
    double dh;
    dh = std::min(front->rect_grid->h[0], front->rect_grid->h[1]);
    dh = std::min(dh, front->rect_grid->h[2]);
    double dt = CFL * dh;
    int nStep = int (front->max_time/dt) + 1;
    dt = front->max_time/nStep;
    front->dt = dt;

    // Initialize state at cell_center
    for (i = 0; i < cell_center.size(); i++)
    {
	getRectangleCenter(i, coords);
	cell_center[i].m_state.setZero();
	getExactSolution(coords,front->time,cell_center[i].m_state);
	cell_center[i].m_state.m_q = cell_center[i].m_state.m_P;
    }
    front->dt = 0.0;

    for (l = 0; l < dim; l++)
    {
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_U[l];
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_U[l] = array[index];
	}
    }

    double averagep, averagep_test;
    double totalVolume, cellVolume;
    totalVolume = 0.0;
    averagep =  0.0;
    averagep_test = 0.0;
    double r;

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);
	r = cell_center[index].m_coords[2];
	cellVolume = r*top_h[0]*top_h[1]*top_h[2];
	totalVolume += cellVolume;
	averagep += cell_center[index].m_state.m_P * cellVolume;
    }
    
    pp_global_sum(&totalVolume, 1);
    pp_global_sum(&averagep, 1);

    printf("\nThe total volume is : %20.16g\n", totalVolume);
    
    averagep = averagep/totalVolume;

    printf("\nThe integral of pressure is : %20.16g\n", averagep);

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);
	r = cell_center[index].m_coords[2];
	cellVolume = r*top_h[0]*top_h[1]*top_h[2];
	array[index] = cell_center[index].m_state.m_P - averagep;
	averagep_test += array[index] * cellVolume;
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

    printf("\nAfter adjusting, the integral of pressure right now is %20.16g\n", averagep_test);

    computeGradientQ();
    copyMeshStates();
    setAdvectionDt();

}

void Incompress_Solver_Smooth_3D_Cylindrical_Debug::solve(double dt)
{
    //    Incompress_Solver_Smooth_3D_Cartesian::solve(dt);

    printf("\nIn this solve(), dt = %.10g, old time = %.10g\n",dt, front->time);
    if(m_testingCase==TEST_DIFFUSION)
    {
	solve_diffusionOnly(dt);
	return;
    }

    m_t_old = front->time;
    m_t_int = front->time + front->dt/2;
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
    setSmoProOnePhase();
    stop_clock("setSmoothedProperties");
    if (debugging("trace"))
	printf("Passed setSmoothedProperties()\n");

    //    // 1) solve for intermediate velocity
        start_clock("computeAdvection");
    //    computeAdvection();
 
	compAdvectionTerm_decoupled();
        stop_clock("computeAdvection");
    //    if (debugging("trace"))
    //	printf("max_speed after computeAdvection(): %20.14f\n",
    //		max_speed);

    start_clock("compDiffWithSmoothProperty");
    compDiffWithSmoothProperty_2nd_decoupled_Shuqiang();
    stop_clock("compDiffWithSmoothProperty");

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
	computeProjection_Shuqiang();
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
    {
	sampleVelocity();
    }

    if (debugging("trace"))
	printf("max_speed after computeNewVelocity(): %20.14f\n",
		max_speed);

    start_clock("copyMeshStates");
    copyMeshStates();
    stop_clock("copyMeshStates");

    setAdvectionDt();
    stop_clock("solve");
}

void Incompress_Solver_Smooth_3D_Cylindrical_Debug::solve_diffusionOnly(double dt)
{
    m_t_old = front->time;
    m_t_int = front->time + front->dt/2;
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



    start_clock("compDiffWithSmoothProperty");
//    compDiffWithSmoothProperty_decoupled();
    compDiffWithSmoothProperty_2nd_decoupled_Shuqiang();
    stop_clock("compDiffWithSmoothProperty");

    start_clock("compSGS");
    //compSGS();	//Subgrid model by Hyunkyun Lim
    stop_clock("compSGS");

    if (debugging("trace"))
	printf("max_speed after compDiffWithSmoothProperty(): %20.14f\n",
		max_speed);


    if (debugging("trace"))
	printf("max_speed after computeNewVelocity(): %20.14f\n",
		max_speed);

    start_clock("copyMeshStates");
    copyMeshStates();
    stop_clock("copyMeshStates");

    setAdvectionDt();
    stop_clock("solve");
}

void Incompress_Solver_Smooth_3D_Cylindrical_Debug::getVelocity(double *p, double *U)
{
    int i;
    for (i = 0; i < dim; ++i)
	U[i] = 0.0;
}

/**
* This function is written here so that it can be overwritten by a derived
* class.
*/
boolean Incompress_Solver_Smooth_3D_Cylindrical_Debug::FT_StateStructAtGridCrossing_tmp(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	COMPONENT comp,
	Locstate *state,
	HYPER_SURF **hs,
	double *crx_coords,
	double t)
{
    boolean bRet = FT_StateStructAtGridCrossing(front,icoords,dir,comp,state,hs,crx_coords);

    if(bRet)
    {
	static int debug_index = 0;
	L_STATE exact_state;
	//	getExactSolution(crx_coords,front->time,exact_state);
	getExactSolution(crx_coords,t,exact_state);
	//	if(m_bStartDebugging)
	//	{
	//	    printf("crx_coords={% 12.8f,% 12.8f}, U={% 12.8f,% 12.8f}, debug_index=%5d\n",
	//		    crx_coords[0], crx_coords[1],
	//		    exact_state.m_U[0], exact_state.m_U[1],
	//		    debug_index++);
	//
	//	}

	STATE *fstate = (STATE*)(*state);
	fstate->vel[0] = exact_state.m_U[0];
	fstate->vel[1] = exact_state.m_U[1];
	fstate->vel[2] = exact_state.m_U[2];
    }
    return bRet;
}

/**
 * the source term does not contain density.
* see also:
* EBM2D_Liquid_Solution::getExactSolution, using TEST_BROWN.
*/

void Incompress_Solver_Smooth_3D_Cylindrical_Debug::computeSourceTerm_Adv(
	double *coords,
	L_STATE &state)
{
    //    Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm(coords,state);

    double t = front->time;
    double theta = coords[0];
    double z     = coords[1];
    double r     = coords[2];

//    double x = r;
//    double y = theta;

    double mu = m_mu[0];
    double a = 1;
    double d = 1;

    switch(m_testingCase)
    {
    case TEST_DIFFUSION:
	// diffusion
//	state.m_U[0] = a * exp(-d * d * t) * d * d * (-sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) - sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) + cos(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) + cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) + mu * sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) + mu * sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) - mu * cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) - mu * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * cos(theta));
//	state.m_U[1] = -a * exp(-d * d * t) * d * d * (-exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) - exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + mu * exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) + mu * exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta)));
//	state.m_U[2] = -a * exp(-d * d * t) * d * d * (-cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) - cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) - sin(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) - sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) + mu * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * sin(theta) + mu * cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) + mu * sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) + mu * cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))));
	state.m_U[0] = (0.25e2 * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.25e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * theta * theta + 0.25e2 * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r - 0.50e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t - 0.50e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.150e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.150e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.150e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.100e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * theta * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r + 0.64e2 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t + 0.64e2 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.50e2 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) (t * t)) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * pow(r, -0.2e1) / 0.100e3;
//	state.m_U[0] = -exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) * (double) (-3 *  pow((double) r, (double) 4) - 3 * r * r * theta * theta - 3 * r * r * z * z + 9 * r * r * t + 9 * r * r - 18 * mu * r * r * t - 18 * mu * r * r + 4 * mu *  pow((double) r, (double) 4) - 6 * mu * t - 6 * mu + 4 * mu * theta * theta + 4 * mu * z * z * r * r) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * (double)  pow((double) r, (double) (-2)) / 0.27e2;
	state.m_U[1] = -exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) * (double) (-3 *  pow((double) r, (double) 4) - 3 * r * r * theta * theta - 3 * r * r * z * z + 9 * r * r * t + 9 * r * r - 18 * mu * r * r * t - 18 * mu * r * r + 4 * mu *  pow((double) r, (double) 4) - 6 * mu * t - 6 * mu + 4 * mu * theta * theta + 4 * mu * z * z * r * r) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * (double)  pow((double) r, (double) (-2)) / 0.27e2;
//	state.m_U[2] = -exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) * (double) (-3 *  pow((double) r, (double) 4) - 3 * r * r * theta * theta - 3 * r * r * z * z + 9 * r * r * t + 9 * r * r - 18 * mu * r * r * t - 18 * mu * r * r + 4 * mu *  pow((double) r, (double) 4) - 6 * mu * t - 6 * mu + 4 * mu * theta * theta + 4 * mu * z * z * r * r) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * (double)  pow((double) r, (double) (-2)) / 0.27e2;
	state.m_U[2] = -(-0.20e2 * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.20e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * theta * theta - 0.20e2 * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * r * r + 0.100e3 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t + 0.100e3 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.120e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t - 0.120e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.240e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t - 0.140e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * theta * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * r * r - 0.100e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) (t * t) + 0.125e3 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.125e3 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1)) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * pow(r, -0.2e1) / 0.125e3;

//	state.m_U[0] = (0.25e2 * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.25e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * theta * theta + 0.25e2 * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r - 0.50e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t - 0.50e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.50e2 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) (t * t) + 0.150e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.100e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.150e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r * (double) t + 0.150e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r - 0.50e2 * mu * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * theta * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r + 0.64e2 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t + 0.64e2 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1)) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * pow(r, -0.2e1) / 0.100e3;
//	state.m_U[1] = -exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) * (double) (-3 *  pow((double) r, (double) 4) - 3 * r * r * theta * theta - 3 * r * r * z * z + 9 * r * r * t + 9 * r * r - 18 * mu * r * r * t - 18 * mu * r * r + 4 * mu *  pow((double) r, (double) 4) - 6 * mu * t - 6 * mu + 4 * mu * theta * theta + 4 * mu * z * z * r * r) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * (double)  pow((double) r, (double) (-2)) / 0.27e2;
//	state.m_U[2] = -(-0.20e2 * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.20e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * theta * theta - 0.20e2 * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * r * r + 0.100e3 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t + 0.100e3 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.120e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t - 0.120e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.240e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t - 0.140e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * theta * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * r * r - 0.100e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) (t * t) + 0.125e3 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.125e3 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1)) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * pow(r, -0.2e1) / 0.125e3;
	break;

    case TEST_EthierSteinman:
	// incompressible without gradP
//	state.m_U[0] = a * exp(-d * d * t) * (mu * sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * d * d - mu * cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * d * d - d * d * sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) + d * d * cos(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) - sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * d * d + cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * d * d + a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * cos(theta) * cos(a * z + d * r * cos(theta)) - a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * cos(theta) * sin(a * r * sin(theta) + d * z) * d - a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * cos(a * z + d * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z) * d - a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - a * a * exp(-d * d * t + 0.2e1 * a * r * cos(theta)) * sin(theta) + a * a * exp(-d * d * t + 0.2e1 * a * r * sin(theta)) * cos(theta) - mu * cos(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * d * d - a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * sin(theta) - a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * cos(theta) + a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * cos(theta) + mu * sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) * d * d + a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) * sin(r * (a * cos(theta) + d * sin(theta))) + a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(theta) * sin(a * z + d * r * cos(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * d + a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * cos(a * r * sin(theta) + d * z) + a * exp(-d * d * t + a * z + a * r * sin(theta)) * cos(r * (a * cos(theta) + d * sin(theta))) * cos(theta) * cos(a * z + d * r * cos(theta)) * d - a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * sin(a * r * sin(theta) + d * z));
//	state.m_U[1] = a * exp(-d * d * t) * (-a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * sin(a * z + d * r * cos(theta)) + a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) + d * d * exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) + d * d * exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + a * a * exp(-d * d * t + 0.2e1 * a * z) + a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * cos(a * z + d * r * cos(theta)) - mu * exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) * d * d - a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * sin(a * z + d * r * cos(theta)) * d + a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * cos(r * (a * cos(theta) + d * sin(theta))) - mu * exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) * d * d + a * exp(-d * d * t + a * z + a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * cos(r * (a * cos(theta) + d * sin(theta))) * d);
//	state.m_U[2] = -a * exp(-d * d * t) * (a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * sin(theta) - a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * sin(theta) - d * d * cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) - d * d * sin(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) - a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * sin(a * r * sin(theta) + d * z) - a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * cos(a * z + d * r * cos(theta)) * d + a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * cos(theta) * sin(r * (a * cos(theta) + d * sin(theta))) * d + mu * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * d * d * sin(theta) - cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * d * d - sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * d * d + a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * cos(theta) * sin(r * (a * cos(theta) + d * sin(theta))) - a * a * exp(-d * d * t + 0.2e1 * a * r * cos(theta)) * cos(theta) - a * a * exp(-d * d * t + 0.2e1 * a * r * sin(theta)) * sin(theta) + mu * cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * d * d + mu * sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * d * d - a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * sin(theta) * cos(a * z + d * r * cos(theta)) + a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * sin(theta) * sin(a * r * sin(theta) + d * z) * d - a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * cos(a * z + d * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z) * d - a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * cos(theta) + mu * cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) * d * d - a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * cos(a * r * sin(theta) + d * z));
//	state.m_U[0] = a * d * d * (mu * exp(-d * d * t + a * z) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - mu * exp(-d * d * t + a * r * sin(theta)) * cos(theta) * sin(a * z + d * r * cos(theta)) - mu * exp(-d * d * t + a * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) - exp(-d * d * t + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) - exp(-d * d * t + a * z) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) + exp(-d * d * t + a * r * sin(theta)) * cos(theta) * sin(a * z + d * r * cos(theta)) + exp(-d * d * t + a * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z));
	state.m_U[0] = a * d * d * (mu * exp(-d * d * t + a * z) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - mu * exp(-d * d * t + a * r * sin(theta)) * cos(theta) * sin(a * z + d * r * cos(theta)) - mu * exp(-d * d * t + a * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) + exp(-d * d * t + a * r * sin(theta)) * cos(theta) * sin(a * z + d * r * cos(theta)) - exp(-d * d * t + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) - exp(-d * d * t + a * z) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) + exp(-d * d * t + a * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z));
//	state.m_U[1] = -a * d * d * (-exp(-d * d * t + a * z) * sin(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + mu * exp(-d * d * t + a * z) * sin(r * (a * cos(theta) + d * sin(theta))) + mu * exp(-d * d * t + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)));
	state.m_U[1] = -a * d * d * (-exp(-d * d * t + a * z) * sin(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + mu * exp(-d * d * t + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + mu * exp(-d * d * t + a * z) * sin(r * (a * cos(theta) + d * sin(theta))));
//	state.m_U[2] = -a * d * d * (mu * exp(-d * d * t + a * r * cos(theta)) * cos(theta) * sin(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * sin(theta) + mu * exp(-d * d * t + a * z) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * z) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * sin(theta)) * sin(theta) * sin(a * z + d * r * cos(theta)) - exp(-d * d * t + a * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z) - exp(-d * d * t + a * r * cos(theta)) * cos(theta) * sin(a * r * sin(theta) + d * z));
	state.m_U[2] = -a * d * d * (mu * exp(-d * d * t + a * r * cos(theta)) * cos(theta) * sin(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * sin(theta) + mu * exp(-d * d * t + a * z) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * cos(theta)) * cos(theta) * sin(a * r * sin(theta) + d * z) - exp(-d * d * t + a * z) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * sin(theta)) * sin(theta) * sin(a * z + d * r * cos(theta)) - exp(-d * d * t + a * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z));
	break;
    case TEST_EthierSteinman_Simple:
	state.m_U[0] = exp(-t) * (-0.2e1 * r * z * z * exp(-t) * pow(cos(theta), 0.2e1) - 0.6e1 * r * r * cos(theta) + 0.2e1 * exp(-t) * pow(r, 0.3e1) * pow(cos(theta), 0.2e1) + 0.6e1 * r * r * sin(theta) - exp(-t) * pow(r, 0.3e1) + 0.6e1 * z * z * sin(theta) - 0.6e1 * z * z * cos(theta) + 0.36e2 * mu * sin(theta) - 0.36e2 * mu * cos(theta) + r * z * z * exp(-t)) / 0.72e2;
	state.m_U[1] = z * exp(-t) * (0.6e1 * r * cos(theta) + 0.6e1 * r * sin(theta) - exp(-t) * z * z + 0.2e1 * exp(-t) * r * r * cos(theta) * sin(theta) + 0.6e1) / 0.36e2;
	state.m_U[2] = exp(-t) * (-0.2e1 * r * z * z * exp(-t) * cos(theta) * sin(theta) - 0.6e1 * r * r * cos(theta) + 0.12e2 * r - 0.6e1 * r * r * sin(theta) + 0.2e1 * exp(-t) * pow(r, 0.3e1) * sin(theta) * cos(theta) + exp(-t) * pow(r, 0.3e1) - 0.6e1 * z * z * sin(theta) - 0.6e1 * z * z * cos(theta) - 0.36e2 * mu * sin(theta) - 0.36e2 * mu * cos(theta) - r * z * z * exp(-t)) / 0.72e2;
	break;

    case TEST_EthierSteinman_Simple2:
    	state.m_U[0] = exp(-t) * (-0.6e1 * theta * exp(-t) * r * z * z - 0.180e3 * pow(r, 0.4e1) + 0.360e3 * theta * r + 0.12e2 * exp(-t) * pow(r, 0.3e1) * theta * theta - 0.6e1 * exp(-t) * pow(r, 0.3e1) * z * z + 0.6e1 * theta * exp(-t) * pow(r, 0.3e1) + 0.6e1 * pow(theta, 0.3e1) * exp(-t) * r + pow(z, 0.4e1) * exp(-t) * r + 0.3e1 * exp(-t) * r * pow(theta, 0.4e1) - (double) (360 * mu) - 0.180e3 * r * r * theta * theta - 0.180e3 * r * r * z * z + 0.9e1 * exp(-t) * pow(r, 0.5e1) - 0.900e3 * (double) mu * r * r - 0.720e3 * (double) mu * theta + 0.180e3 * (double) mu * theta * theta + 0.180e3 * (double) mu * z * z) * pow(r, -0.2e1) / 0.10800e5;
    	state.m_U[1] = z * exp(-t) * (0.4e1 * theta * exp(-t) * r * z * z - 0.3e1 * z * z * exp(-t) * r + 0.4e1 * z * z * exp(-t) * r * theta * theta + 0.270e3 * pow(r, 0.4e1) + 0.180e3 * pow(r, 0.3e1) + 0.3e1 * exp(-t) * r * theta * theta + 0.180e3 * r * r * theta - 0.3e1 * exp(-t) * pow(r, 0.3e1) + 0.6e1 * exp(-t) * pow(r, 0.3e1) * theta * theta + 0.2e1 * exp(-t) * pow(r, 0.3e1) * z * z + 0.18e2 * theta * exp(-t) * pow(r, 0.3e1) + 0.6e1 * pow(theta, 0.3e1) * exp(-t) * r + pow(z, 0.4e1) * exp(-t) * r + 0.3e1 * exp(-t) * r * pow(theta, 0.4e1) + (double) (180 * mu) + 0.90e2 * r * r * theta * theta + 0.30e2 * r * r * z * z + 0.9e1 * exp(-t) * pow(r, 0.5e1) + 0.450e3 * (double) mu * r * r + 0.180e3 * (double) mu * theta + 0.90e2 * (double) mu * theta * theta + 0.30e2 * (double) mu * z * z) * pow(r, -0.3e1) / 0.5400e4;
    	state.m_U[2] = exp(-t) * (-0.6e1 * theta * exp(-t) * r * z * z - 0.12e2 * z * z * exp(-t) * r * theta * theta - 0.180e3 * pow(r, 0.4e1) + 0.360e3 * pow(r, 0.3e1) + 0.3e1 * exp(-t) * pow(r, 0.5e1) + 0.6e1 * theta * exp(-t) * pow(r, 0.3e1) - 0.18e2 * exp(-t) * pow(r, 0.3e1) * z * z - 0.5e1 * pow(z, 0.4e1) * exp(-t) * r + 0.6e1 * pow(theta, 0.3e1) * exp(-t) * r - 0.3e1 * exp(-t) * r * pow(theta, 0.4e1) - (double) (360 * mu) - 0.180e3 * r * r * theta * theta - 0.180e3 * r * r * z * z - 0.900e3 * (double) mu * r * r + 0.720e3 * (double) mu * theta + 0.180e3 * (double) mu * theta * theta + 0.180e3 * (double) mu * z * z) * pow(r, -0.2e1) / 0.10800e5;
    	break;
    default:
	assert(false);
    }
}


void Incompress_Solver_Smooth_3D_Cylindrical_Debug::computeSourceTerm(
	double *coords,
	L_STATE &state)
{
    //    Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm(coords,state);

    //    double t = front->time + front->dt/2;
    double t = m_t_int;
    double theta = coords[0];
    double z     = coords[1];
    double r     = coords[2];

//    double x = r;
//    double y = theta;

    double mu = m_mu[0];
    double a = 1;
    double d = 1;

    switch(m_testingCase)
    {
    case TEST_DIFFUSION:
	// diffusion
//	state.m_U[0] = a * exp(-d * d * t) * d * d * (-sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) - sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) + cos(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) + cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) + mu * sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) + mu * sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) - mu * cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) - mu * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * cos(theta));
//	state.m_U[1] = -a * exp(-d * d * t) * d * d * (-exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) - exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + mu * exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) + mu * exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta)));
//	state.m_U[2] = -a * exp(-d * d * t) * d * d * (-cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) - cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) - sin(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) - sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) + mu * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * sin(theta) + mu * cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) + mu * sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) + mu * cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))));
	state.m_U[0] = (0.25e2 * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.25e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * theta * theta + 0.25e2 * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r - 0.50e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t - 0.50e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.150e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.150e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.150e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.100e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * theta * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r + 0.64e2 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t + 0.64e2 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.50e2 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) (t * t)) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * pow(r, -0.2e1) / 0.100e3;
//	state.m_U[0] = -exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) * (double) (-3 *  pow((double) r, (double) 4) - 3 * r * r * theta * theta - 3 * r * r * z * z + 9 * r * r * t + 9 * r * r - 18 * mu * r * r * t - 18 * mu * r * r + 4 * mu *  pow((double) r, (double) 4) - 6 * mu * t - 6 * mu + 4 * mu * theta * theta + 4 * mu * z * z * r * r) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * (double)  pow((double) r, (double) (-2)) / 0.27e2;
	state.m_U[1] = -exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) * (double) (-3 *  pow((double) r, (double) 4) - 3 * r * r * theta * theta - 3 * r * r * z * z + 9 * r * r * t + 9 * r * r - 18 * mu * r * r * t - 18 * mu * r * r + 4 * mu *  pow((double) r, (double) 4) - 6 * mu * t - 6 * mu + 4 * mu * theta * theta + 4 * mu * z * z * r * r) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * (double)  pow((double) r, (double) (-2)) / 0.27e2;
//	state.m_U[2] = -exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) * (double) (-3 *  pow((double) r, (double) 4) - 3 * r * r * theta * theta - 3 * r * r * z * z + 9 * r * r * t + 9 * r * r - 18 * mu * r * r * t - 18 * mu * r * r + 4 * mu *  pow((double) r, (double) 4) - 6 * mu * t - 6 * mu + 4 * mu * theta * theta + 4 * mu * z * z * r * r) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * (double)  pow((double) r, (double) (-2)) / 0.27e2;
	state.m_U[2] = -(-0.20e2 * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.20e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * theta * theta - 0.20e2 * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * r * r + 0.100e3 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t + 0.100e3 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.120e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t - 0.120e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.240e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t - 0.140e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * theta * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * r * r - 0.100e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) (t * t) + 0.125e3 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.125e3 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1)) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * pow(r, -0.2e1) / 0.125e3;

//	state.m_U[0] = (0.25e2 * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.25e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * theta * theta + 0.25e2 * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r - 0.50e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t - 0.50e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.50e2 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) (t * t) + 0.150e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.100e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) + 0.150e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r * (double) t + 0.150e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r - 0.50e2 * mu * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * theta * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) - 0.50e2 * mu * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * r * r + 0.64e2 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t + 0.64e2 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1)) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * pow(r, -0.2e1) / 0.100e3;
//	state.m_U[1] = -exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) * (double) (-3 *  pow((double) r, (double) 4) - 3 * r * r * theta * theta - 3 * r * r * z * z + 9 * r * r * t + 9 * r * r - 18 * mu * r * r * t - 18 * mu * r * r + 4 * mu *  pow((double) r, (double) 4) - 6 * mu * t - 6 * mu + 4 * mu * theta * theta + 4 * mu * z * z * r * r) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * (double)  pow((double) r, (double) (-2)) / 0.27e2;
//	state.m_U[2] = -(-0.20e2 * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.20e2 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * theta * theta - 0.20e2 * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * r * r + 0.100e3 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t + 0.100e3 * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.120e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t - 0.120e3 * mu * r * r * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * pow(r, 0.4e1) * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) - 0.240e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) t - 0.140e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * theta * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) + 0.16e2 * mu * z * z * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * r * r - 0.100e3 * mu * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.5e1) * (double) (t * t) + 0.125e3 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1) * (double) t + 0.125e3 * mu * theta * exp(-(r * r + theta * theta + z * z) / (double) (t + 1) / 0.2e1)) * (double)  pow((double) (t + 1), (double) (-3)) / 0.3141592654e1 * pow(r, -0.2e1) / 0.125e3;
	break;

    case TEST_EthierSteinman:
	// incompressible without gradP
//	state.m_U[0] = a * exp(-d * d * t) * (mu * sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * d * d - mu * cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * d * d - d * d * sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) + d * d * cos(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) - sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * d * d + cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * d * d + a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * cos(theta) * cos(a * z + d * r * cos(theta)) - a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * cos(theta) * sin(a * r * sin(theta) + d * z) * d - a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * cos(a * z + d * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z) * d - a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - a * a * exp(-d * d * t + 0.2e1 * a * r * cos(theta)) * sin(theta) + a * a * exp(-d * d * t + 0.2e1 * a * r * sin(theta)) * cos(theta) - mu * cos(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * d * d - a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * sin(theta) - a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * cos(theta) + a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * cos(theta) + mu * sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) * d * d + a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) * sin(r * (a * cos(theta) + d * sin(theta))) + a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(theta) * sin(a * z + d * r * cos(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * d + a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * cos(a * r * sin(theta) + d * z) + a * exp(-d * d * t + a * z + a * r * sin(theta)) * cos(r * (a * cos(theta) + d * sin(theta))) * cos(theta) * cos(a * z + d * r * cos(theta)) * d - a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * sin(a * r * sin(theta) + d * z));
//	state.m_U[1] = a * exp(-d * d * t) * (-a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * sin(a * z + d * r * cos(theta)) + a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) + d * d * exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) + d * d * exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + a * a * exp(-d * d * t + 0.2e1 * a * z) + a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * cos(a * z + d * r * cos(theta)) - mu * exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) * d * d - a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * sin(a * z + d * r * cos(theta)) * d + a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * cos(r * (a * cos(theta) + d * sin(theta))) - mu * exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) * d * d + a * exp(-d * d * t + a * z + a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * cos(r * (a * cos(theta) + d * sin(theta))) * d);
//	state.m_U[2] = -a * exp(-d * d * t) * (a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * sin(theta) - a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * sin(theta) - d * d * cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) - d * d * sin(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) - a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * sin(a * r * sin(theta) + d * z) - a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * cos(a * z + d * r * cos(theta)) * d + a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * cos(theta) * sin(r * (a * cos(theta) + d * sin(theta))) * d + mu * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * d * d * sin(theta) - cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * d * d - sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * d * d + a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * cos(theta) * sin(r * (a * cos(theta) + d * sin(theta))) - a * a * exp(-d * d * t + 0.2e1 * a * r * cos(theta)) * cos(theta) - a * a * exp(-d * d * t + 0.2e1 * a * r * sin(theta)) * sin(theta) + mu * cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) * d * d + mu * sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * d * d - a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * sin(theta) * cos(a * z + d * r * cos(theta)) + a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(r * (a * cos(theta) + d * sin(theta))) * sin(theta) * sin(a * r * sin(theta) + d * z) * d - a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * cos(a * z + d * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z) * d - a * a * exp(-d * d * t + a * z + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - a * a * exp(-d * d * t + a * r * sin(theta) + a * r * cos(theta)) * sin(a * z + d * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * cos(theta) + mu * cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) * d * d - a * a * exp(-d * d * t + a * z + a * r * cos(theta)) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) * cos(a * r * sin(theta) + d * z));
//	state.m_U[0] = a * d * d * (mu * exp(-d * d * t + a * z) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - mu * exp(-d * d * t + a * r * sin(theta)) * cos(theta) * sin(a * z + d * r * cos(theta)) - mu * exp(-d * d * t + a * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) - exp(-d * d * t + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) - exp(-d * d * t + a * z) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) + exp(-d * d * t + a * r * sin(theta)) * cos(theta) * sin(a * z + d * r * cos(theta)) + exp(-d * d * t + a * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z));
	state.m_U[0] = a * d * d * (mu * exp(-d * d * t + a * z) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - mu * exp(-d * d * t + a * r * sin(theta)) * cos(theta) * sin(a * z + d * r * cos(theta)) - mu * exp(-d * d * t + a * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) + exp(-d * d * t + a * r * sin(theta)) * cos(theta) * sin(a * z + d * r * cos(theta)) - exp(-d * d * t + a * r * cos(theta)) * sin(theta) * sin(a * r * sin(theta) + d * z) - exp(-d * d * t + a * z) * sin(theta) * cos(r * (a * cos(theta) + d * sin(theta))) + exp(-d * d * t + a * r * cos(theta)) * cos(theta) * cos(a * r * sin(theta) + d * z));
//	state.m_U[1] = -a * d * d * (-exp(-d * d * t + a * z) * sin(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + mu * exp(-d * d * t + a * z) * sin(r * (a * cos(theta) + d * sin(theta))) + mu * exp(-d * d * t + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)));
	state.m_U[1] = -a * d * d * (-exp(-d * d * t + a * z) * sin(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + mu * exp(-d * d * t + a * r * sin(theta)) * cos(a * z + d * r * cos(theta)) + mu * exp(-d * d * t + a * z) * sin(r * (a * cos(theta) + d * sin(theta))));
//	state.m_U[2] = -a * d * d * (mu * exp(-d * d * t + a * r * cos(theta)) * cos(theta) * sin(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * sin(theta) + mu * exp(-d * d * t + a * z) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * z) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * sin(theta)) * sin(theta) * sin(a * z + d * r * cos(theta)) - exp(-d * d * t + a * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z) - exp(-d * d * t + a * r * cos(theta)) * cos(theta) * sin(a * r * sin(theta) + d * z));
	state.m_U[2] = -a * d * d * (mu * exp(-d * d * t + a * r * cos(theta)) * cos(theta) * sin(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z) + mu * exp(-d * d * t + a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) * sin(theta) + mu * exp(-d * d * t + a * z) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * cos(theta)) * cos(theta) * sin(a * r * sin(theta) + d * z) - exp(-d * d * t + a * z) * cos(theta) * cos(r * (a * cos(theta) + d * sin(theta))) - exp(-d * d * t + a * r * sin(theta)) * sin(theta) * sin(a * z + d * r * cos(theta)) - exp(-d * d * t + a * r * cos(theta)) * sin(theta) * cos(a * r * sin(theta) + d * z));
	break;
    case TEST_EthierSteinman_Simple:
	state.m_U[0] = exp(-t) * (-0.2e1 * r * z * z * exp(-t) * pow(cos(theta), 0.2e1) - 0.6e1 * r * r * cos(theta) + 0.2e1 * exp(-t) * pow(r, 0.3e1) * pow(cos(theta), 0.2e1) + 0.6e1 * r * r * sin(theta) - exp(-t) * pow(r, 0.3e1) + 0.6e1 * z * z * sin(theta) - 0.6e1 * z * z * cos(theta) + 0.36e2 * mu * sin(theta) - 0.36e2 * mu * cos(theta) + r * z * z * exp(-t)) / 0.72e2;
	state.m_U[1] = z * exp(-t) * (0.6e1 * r * cos(theta) + 0.6e1 * r * sin(theta) - exp(-t) * z * z + 0.2e1 * exp(-t) * r * r * cos(theta) * sin(theta) + 0.6e1) / 0.36e2;
	state.m_U[2] = exp(-t) * (-0.2e1 * r * z * z * exp(-t) * cos(theta) * sin(theta) - 0.6e1 * r * r * cos(theta) + 0.12e2 * r - 0.6e1 * r * r * sin(theta) + 0.2e1 * exp(-t) * pow(r, 0.3e1) * sin(theta) * cos(theta) + exp(-t) * pow(r, 0.3e1) - 0.6e1 * z * z * sin(theta) - 0.6e1 * z * z * cos(theta) - 0.36e2 * mu * sin(theta) - 0.36e2 * mu * cos(theta) - r * z * z * exp(-t)) / 0.72e2;
	break;

    case TEST_EthierSteinman_Simple2:
    	state.m_U[0] = exp(-t) * (-0.6e1 * theta * exp(-t) * r * z * z - 0.180e3 * pow(r, 0.4e1) + 0.360e3 * theta * r + 0.12e2 * exp(-t) * pow(r, 0.3e1) * theta * theta - 0.6e1 * exp(-t) * pow(r, 0.3e1) * z * z + 0.6e1 * theta * exp(-t) * pow(r, 0.3e1) + 0.6e1 * pow(theta, 0.3e1) * exp(-t) * r + pow(z, 0.4e1) * exp(-t) * r + 0.3e1 * exp(-t) * r * pow(theta, 0.4e1) - (double) (360 * mu) - 0.180e3 * r * r * theta * theta - 0.180e3 * r * r * z * z + 0.9e1 * exp(-t) * pow(r, 0.5e1) - 0.900e3 * (double) mu * r * r - 0.720e3 * (double) mu * theta + 0.180e3 * (double) mu * theta * theta + 0.180e3 * (double) mu * z * z) * pow(r, -0.2e1) / 0.10800e5;
    	state.m_U[1] = z * exp(-t) * (0.4e1 * theta * exp(-t) * r * z * z - 0.3e1 * z * z * exp(-t) * r + 0.4e1 * z * z * exp(-t) * r * theta * theta + 0.270e3 * pow(r, 0.4e1) + 0.180e3 * pow(r, 0.3e1) + 0.3e1 * exp(-t) * r * theta * theta + 0.180e3 * r * r * theta - 0.3e1 * exp(-t) * pow(r, 0.3e1) + 0.6e1 * exp(-t) * pow(r, 0.3e1) * theta * theta + 0.2e1 * exp(-t) * pow(r, 0.3e1) * z * z + 0.18e2 * theta * exp(-t) * pow(r, 0.3e1) + 0.6e1 * pow(theta, 0.3e1) * exp(-t) * r + pow(z, 0.4e1) * exp(-t) * r + 0.3e1 * exp(-t) * r * pow(theta, 0.4e1) + (double) (180 * mu) + 0.90e2 * r * r * theta * theta + 0.30e2 * r * r * z * z + 0.9e1 * exp(-t) * pow(r, 0.5e1) + 0.450e3 * (double) mu * r * r + 0.180e3 * (double) mu * theta + 0.90e2 * (double) mu * theta * theta + 0.30e2 * (double) mu * z * z) * pow(r, -0.3e1) / 0.5400e4;
    	state.m_U[2] = exp(-t) * (-0.6e1 * theta * exp(-t) * r * z * z - 0.12e2 * z * z * exp(-t) * r * theta * theta - 0.180e3 * pow(r, 0.4e1) + 0.360e3 * pow(r, 0.3e1) + 0.3e1 * exp(-t) * pow(r, 0.5e1) + 0.6e1 * theta * exp(-t) * pow(r, 0.3e1) - 0.18e2 * exp(-t) * pow(r, 0.3e1) * z * z - 0.5e1 * pow(z, 0.4e1) * exp(-t) * r + 0.6e1 * pow(theta, 0.3e1) * exp(-t) * r - 0.3e1 * exp(-t) * r * pow(theta, 0.4e1) - (double) (360 * mu) - 0.180e3 * r * r * theta * theta - 0.180e3 * r * r * z * z - 0.900e3 * (double) mu * r * r + 0.720e3 * (double) mu * theta + 0.180e3 * (double) mu * theta * theta + 0.180e3 * (double) mu * z * z) * pow(r, -0.2e1) / 0.10800e5;
    	break;
    default:
	assert(false);
    }
}

/**
* TEST_EthierSteinman in EBM3D_Liquid_Solution.
*/
void Incompress_Solver_Smooth_3D_Cylindrical_Debug::getExactSolution(
	double *coords,
	double t,
	L_STATE &state)
{
    //    double pi = 0.3141592654e1;
    //    double omega = 1+sin(2*pi*t*t);
    double theta = coords[0];
    double z     = coords[1];
    double r     = coords[2];

//    double x = r;
//    double y = theta;

    double mu = m_mu[0];
    double a = 1;
    double d = 1;
    double t_int = t - (front->dt)/2.0;

    switch(m_testingCase)
    {
    case TEST_DIFFUSION:
	// diffusion
//	state.m_U[0] = -a * exp(-d * d * t) * (-sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) - sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) + cos(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) + cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z));
//	state.m_U[1] = -a * (exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) + exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta))) * exp(-d * d * t);
//	state.m_U[2] = -a * exp(-d * d * t) * (cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) + cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) + sin(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) + sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z));
	state.m_U[0] = exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.2e1) / 0.3141592654e1 / (double) (t + 1) / 0.2e1;
//	state.m_U[0] = exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) / 0.3141592654e1 / (double) (t + 1) / 0.3e1;
	state.m_U[1] = exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) / 0.3141592654e1 / (double) (t + 1) / 0.3e1;
//	state.m_U[2] = exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.3e1) / 0.3141592654e1 / (double) (t + 1) / 0.3e1;
	state.m_U[2] = 0.4e1 / 0.5e1 * exp(-(double) ((r * r + theta * theta + z * z) / (t + 1)) / 0.5e1) / 0.3141592654e1 / (double) (t + 1);

//	state.m_U[0] = exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.2e1) / 0.3141592654e1 / (double) (t + 1) / 0.2e1;
//	state.m_U[1] = exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.3e1) / 0.3141592654e1 / (double) (t + 1) / 0.3e1;
//	state.m_U[2] = 0.4e1 / 0.5e1 * exp(-(double) ((x * x + y * y + z * z) / (t + 1)) / 0.5e1) / 0.3141592654e1 / (double) (t + 1);

	break;

    case TEST_EthierSteinman:
	// incompressible
	state.m_U[0] = a * exp(-d * d * t) * (sin(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) + sin(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) - cos(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) - cos(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z));
	state.m_U[1] = -a * (exp(a * z) * sin(r * (a * cos(theta) + d * sin(theta))) + exp(a * r * sin(theta)) * cos(a * z + d * r * cos(theta))) * exp(-d * d * t);
	state.m_U[2] = -a * exp(-d * d * t) * (cos(theta) * exp(a * r * cos(theta)) * sin(a * r * sin(theta) + d * z) + cos(theta) * exp(a * z) * cos(r * (a * cos(theta) + d * sin(theta))) + sin(theta) * exp(a * r * sin(theta)) * sin(a * z + d * r * cos(theta)) + sin(theta) * exp(a * r * cos(theta)) * cos(a * r * sin(theta) + d * z));
	state.m_P    = -a * a * (exp(0.2e1 * a * r * cos(theta)) + exp(0.2e1 * a * r * sin(theta)) + exp(0.2e1 * a * z) + 0.2e1 * sin(r * (a * cos(theta) + d * sin(theta))) * cos(a * z + d * r * cos(theta)) * exp(a * (r * sin(theta) + z)) + 0.2e1 * sin(a * r * sin(theta) + d * z) * cos(r * (a * cos(theta) + d * sin(theta))) * exp(a * (z + r * cos(theta))) + 0.2e1 * sin(a * z + d * r * cos(theta)) * cos(a * r * sin(theta) + d * z) * exp(a * r * (cos(theta) + sin(theta)))) * exp(-0.2e1 * d * d * t_int) / 0.2e1;
	break;
    case TEST_EthierSteinman_Simple:
	state.m_U[0] = exp(-t) * (-r * r * sin(theta) + r * r * cos(theta) - z * z * sin(theta) + z * z * cos(theta)) / 0.12e2;
	state.m_U[1] = -z * exp(-t) * r * (cos(theta) + sin(theta)) / 0.6e1;
	state.m_U[2] = exp(-t) * (r * r * cos(theta) + r * r * sin(theta) + z * z * cos(theta) + z * z * sin(theta)) / 0.12e2;
	state.m_P = exp(-t_int) * (r * r + z * z) / 0.12e2;
	 break;
    case TEST_EthierSteinman_Simple2:
	state.m_U[0] = (r * r + theta * theta + z * z) * exp(-t) / 0.60e2;
	state.m_U[1] = -(double) z * (double) (9 * r * r + 3 * theta * theta + 6 * theta + z * z) / (double) r * exp(-t) / 0.180e3;
	state.m_U[2] = (r * r + theta * theta + z * z) * exp(-t) / 0.60e2;
	state.m_P = (r * r + theta * theta + z * z) * exp(-t_int) / 0.60e2;
	break;
    default:
	assert(false);
    }
}


double Incompress_Solver_Smooth_3D_Cylindrical_Debug::getAveragePressure(
	double t)
{
    int icoords[3], index;
    double *coords, cellVolume; // = top_h[0]*top_h[1]*top_h[2];

    L_STATE state_exact;
    // calculate an average of the pressure
    double P_average = 0, totalVolume = 0;

    for(int k=kmin; k<=kmax; k++)
	for(int j=jmin; j<=jmax; j++)
	    for(int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;
		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		cellVolume = getCellVolume(icoords);

		totalVolume += cellVolume;

//		if(bPrintExact)
		{
		    getExactSolution(coords, t, state_exact);
		    P_average += state_exact.m_P * cellVolume;
		}
	    }
    P_average = P_average/totalVolume;
    return P_average;
}
void Incompress_Solver_Smooth_3D_Cylindrical_Debug::saveStates_Tecplot(
	const char*out_name,
	double t,
	bool bPrintCoords,
	bool bPrintExact,
	bool bPrintError)
{
    bPrintExact = true;
    /*
    char filename[200];
    sprintf(filename,"%s/states_%s.plt",out_name,
	    right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, Z, U, V, W, P");
    if(bPrintExact)
	fprintf(hfile, ", U_exact, V_exact, W_exact, P_exact");
    if(bPrintError)
	fprintf(hfile, ", U_error, V_error, W_error, P_error");
    fprintf(hfile, "\n");

    fprintf(hfile, "ZONE I=%d, J=%d, K=%d, F=POINT\n",
	    imax-imin+1, jmax-jmin+1, kmax-kmin+1);
*/
    int c, icoords[3];
    double *coords, cellVolume; // = top_h[0]*top_h[1]*top_h[2];

    int max_ij[4][4], index;
    max_ij[0][0] = max_ij[0][1] = max_ij[0][2] =-1;
    max_ij[1][0] = max_ij[1][1] = max_ij[1][2] =-1;
    max_ij[2][0] = max_ij[2][1] = max_ij[2][2] =-1;
    max_ij[3][0] = max_ij[3][1] = max_ij[3][2] =-1;

    L_STATE error, L1, L2, LInf;

    L_STATE state, state_exact;

    double average_exactp = 0.0;
    double totalVolume = 0.0;
    double r;

    for(int k=kmin; k<=kmax; k++)
	for(int j=jmin; j<=jmax; j++)
	    for(int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;
		r = coords[2];
		cellVolume = r*top_h[0]*top_h[1]*top_h[2];
		
		getExactSolution(coords, t, state_exact);
		average_exactp = average_exactp + state_exact.m_P*cellVolume;
		totalVolume += cellVolume;
	    }
    average_exactp = average_exactp / totalVolume;

    printf("\nThe total Volume is : %20.16g\n", totalVolume);

    for(int k=kmin; k<=kmax; k++)
	for(int j=jmin; j<=jmax; j++)
	    for(int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;
		state = cell_center[index].m_state;
		cellVolume = coords[2]*top_h[0]*top_h[1]*top_h[2];

		if(bPrintExact)
		    getExactSolution(coords, t, state_exact);
		else
		{
		    state_exact.m_U[0] = state_exact.m_U[1] =
			    state_exact.m_U[2] = state_exact.m_P = 0;
		}

		state_exact.m_P = state_exact.m_P - average_exactp;

		error.m_U[0] = fabs(state.m_U[0] - state_exact.m_U[0]);
		error.m_U[1] = fabs(state.m_U[1] - state_exact.m_U[1]);
		error.m_U[2] = fabs(state.m_U[2] - state_exact.m_U[2]);
		error.m_P    = fabs(state.m_P - state_exact.m_P);

		// LInfinity
		for(unsigned int kk=0; kk<3; kk++)
		    if(error.m_U[kk] > LInf.m_U[kk])
		    {
			LInf.m_U[kk] = error.m_U[kk];
			max_ij[kk][0] = i;
			max_ij[kk][1] = j;
			max_ij[kk][2] = k;
		    }
		if(error.m_P > LInf.m_P)
		{
		    LInf.m_P = error.m_P;
		    max_ij[3][0] = i;
		    max_ij[3][1] = j;
		    max_ij[3][2] = k;
		}

		// L1
		L1.m_U[0] += error.m_U[0]*cellVolume;
		L1.m_U[1] += error.m_U[1]*cellVolume;
		L1.m_U[2] += error.m_U[2]*cellVolume;
		L1.m_P += error.m_P*cellVolume;

		// L2
		L2.m_U[0] += error.m_U[0]*error.m_U[0]*cellVolume;
		L2.m_U[1] += error.m_U[1]*error.m_U[1]*cellVolume;
		L2.m_U[2] += error.m_U[2]*error.m_U[2]*cellVolume;
		L2.m_P += error.m_P*error.m_P*cellVolume;

		//			state_print = state;
/*
		if(bPrintCoords)
		    fprintf(hfile, "% 12.8e % 12.8e % 12.8e",
			    coords[0], coords[1], coords[2]);
		else
		    fprintf(hfile, "%3d %3d %3d",
			    i, j, k);

		fprintf(hfile, " % 12.8e % 12.8e % 12.8e % 12.8e",
			state.m_U[0], state.m_U[1], state.m_U[2], state.m_P);
		if(bPrintExact)
		    fprintf(hfile, " % 12.8e % 12.8e % 12.8e % 12.8e",
			    state_exact.m_U[0], state_exact.m_U[1], state_exact.m_U[2], state_exact.m_P);
		if(bPrintError)
		    fprintf(hfile, " % 12.8e % 12.8e % 12.8e % 12.8e",
			    state_exact.m_U[0] - state.m_U[0],
			    state_exact.m_U[1] - state.m_U[1],
			    state_exact.m_U[2] - state.m_U[2],
			    state_exact.m_P - state.m_P);
		fprintf(hfile, "\n");
		*/
	    }

    pp_global_sum(&L1.m_U[0],1);
    pp_global_sum(&L1.m_U[1],1);
    pp_global_sum(&L1.m_U[2],1);
    pp_global_sum(&L1.m_P,1);

    pp_global_sum(&L2.m_U[0],1);
    pp_global_sum(&L2.m_U[1],1);
    pp_global_sum(&L2.m_U[2],1);
    pp_global_sum(&L2.m_P,1);

    pp_global_max(&LInf.m_U[0],1);
    pp_global_max(&LInf.m_U[1],1);
    pp_global_max(&LInf.m_U[2],1);
    pp_global_max(&LInf.m_P,1);

    L2.m_U[0] = sqrt(L2.m_U[0]);
    L2.m_U[1] = sqrt(L2.m_U[1]);
    L2.m_U[2] = sqrt(L2.m_U[2]);
    L2.m_P = sqrt(L2.m_P);

    //fclose(hfile);

    printf("Incompress_Solver_Smooth_3D_Cylindrical_Debug::saveStates_Tecplot: \n");
    printf("\t#### L1  ={%12.8f, %12.8f, %12.8f, %12.8f}\n", L1.m_U[0], L1.m_U[1], L1.m_U[2], L1.m_P);
    printf("\t#### L2  ={%12.8f, %12.8f, %12.8f, %12.8f}\n", L2.m_U[0], L2.m_U[1], L2.m_U[2], L2.m_P);
    printf("\t#### LInf={%12.8f, %12.8f, %12.8f, %12.8f}\n", LInf.m_U[0], LInf.m_U[1], LInf.m_U[2], LInf.m_P);
//    printf("\t###- debug_totalVolume = %12.8f \n", debug_totalVolume);
}

void Incompress_Solver_Smooth_3D_Cylindrical_Debug::saveDivUPhi_Tecplot(
	const char*out_name,
	double t,
	bool bPrintCoords)
{
    char filename[200];
    sprintf(filename,"%s/DivUPhi_%s.plt",out_name,
	    right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, Z, DivU, Phi\n");

    fprintf(hfile, "ZONE I=%d, J=%d, K=%d, F=POINT\n",
	    imax-imin+1, jmax-jmin+1, kmax-kmin+1);

    double *coords;
    int index;
    L_STATE state;

    double max_div_U = 0, max_phi = 0;

    for(int k=kmin; k<=kmax; k++)
	for(int j=jmin; j<=jmax; j++)
	    for(int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;

		state = cell_center[index].m_state;
		if(fabs(state.div_U)>max_div_U)
		    max_div_U = fabs(state.div_U);
		if(fabs(state.m_phi)>max_phi)
		    max_phi = fabs(state.m_phi);

		if(bPrintCoords)
		    fprintf(hfile, "% 12.8e % 12.8e % 12.8e",
			    coords[0], coords[1], coords[2]);
		else
		    fprintf(hfile, "%3d %3d %3d", i, j, k);

		fprintf(hfile, " % 12.8e % 12.8e",
			state.div_U, state.m_phi);

		fprintf(hfile, "\n");
	    }

    fclose(hfile);

    printf("Incompress_Solver_Smooth_3D_Cylindrical_Debug::saveDivUPhi_Tecplot: max_div_U = %12.8f, max_phi = %12.8f\n",
	    max_div_U, max_phi);
}

void Incompress_Solver_Smooth_3D_Cylindrical_Debug::saveParameters_Tecplot(
	const char*out_name,
	double t,
	bool bPrintCoords)
{
    char filename[200];
    sprintf(filename,"%s/parameters_%s.plt",out_name,
	    right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, Z, rho, mu\n");

    fprintf(hfile, "ZONE I=%d, J=%d, K=%d, F=POINT\n",
	    imax-imin+1, jmax-jmin+1, kmax-kmin+1);

    double *coords;
    int index;
    L_STATE state;

    for(int k=kmin; k<=kmax; k++)
	for(int j=jmin; j<=jmax; j++)
	    for(int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;
		state = cell_center[index].m_state;
		if(bPrintCoords)
		    fprintf(hfile, "% 12.8e % 12.8e % 12.8e",
			    coords[0], coords[1], coords[2]);
		else
		    fprintf(hfile, "%3d %3d %3d", i, j, k);

		fprintf(hfile, " % 12.8e % 12.8e",
			state.m_mu, state.m_rho);
		fprintf(hfile, "\n");
	    }

    fclose(hfile);
}

void Incompress_Solver_Smooth_3D_Cylindrical_Debug::getMaximumSourceTerm(
	double t)
{
    double t_bak = m_t_int;

    m_t_int = t;
    double *coords;
    int index;


    L_STATE sourceTerm, sourceTerm_max;

    for(int k=kmin; k<=kmax; k++)
	for(int j=jmin; j<=jmax; j++)
	    for(int i=imin; i<=imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		coords = cell_center[index].m_coords;
		computeSourceTerm(coords,sourceTerm);
		for(int d=0; d<3; d++)
		    if(fabs(sourceTerm.m_U[d]) > sourceTerm_max.m_U[d])
			sourceTerm_max.m_U[d] = fabs(sourceTerm.m_U[d]);
	    }
    printf("Incompress_Solver_Smooth_3D_Cylindrical_Debug::getMaximumSourceTerm: "
	    "sourceTerm_max={%12.8f, %12.8f, %12.8f}\n",
	    sourceTerm_max.m_U[0], sourceTerm_max.m_U[1], sourceTerm_max.m_U[2]);

    m_t_int = t_bak;
}

Incompress_Solver_Smooth_3D_Cylindrical_Debug::TEST_CASE
//Incompress_Solver_Smooth_3D_Cylindrical_Debug::m_testingCase = TEST_EthierSteinman;
Incompress_Solver_Smooth_3D_Cylindrical_Debug::m_testingCase = TEST_EthierSteinman_Simple;
//Incompress_Solver_Smooth_3D_Cylindrical_Debug::m_testingCase = TEST_EthierSteinman_Simple2;
//Incompress_Solver_Smooth_3D_Cylindrical_Debug::m_testingCase = TEST_DIFFUSION;
//Incompress_Solver_Smooth_3D_Cylindrical_Debug::m_testingCase = TEST_ELLIPTIC;
