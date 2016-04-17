#include <cFluid.h>

static double weno5_scal(double *f);
static void matmvec(double *b, double L[5][5], double *x);
static void f2is(double *f, double *s);
static void u2f(double *u, double *f);
static void weno5_get_flux(POINTER,int,int,double**,double**);

extern void WENO_flux(
        POINTER params,
        SWEEP *vst,
        FSWEEP *vflux,
        int n)
{
	static double **u_old;
	static double **flux;
	int i,extend_size,ghost_size;
	SCHEME_PARAMS *scheme_params = (SCHEME_PARAMS*)params;
	double lambda = scheme_params->lambda;

	if (u_old == NULL)
	{
            //NEED TO FIX for general case. This is only for 2 components
	    FT_VectorMemoryAlloc((POINTER*)&u_old,8,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&flux,7,sizeof(double*));
	}
	u_old[0] = vst->dens;
	u_old[1] = vst->momn[0];
	u_old[2] = vst->momn[1];
	u_old[3] = vst->momn[2];
	u_old[4] = vst->engy;
	u_old[5] = vst->pres;
        u_old[6] = vst->pdens[0];
        u_old[7] = vst->pdens[1];
        u_old[8] = vst->gamma;
	flux[0] = vflux->dens_flux;
	flux[1] = vflux->momn_flux[0];
	flux[2] = vflux->momn_flux[1];
	flux[3] = vflux->momn_flux[2];
	flux[4] = vflux->engy_flux;
        flux[5] = vflux->pdens_flux[0];
        flux[6] = vflux->pdens_flux[1];

	ghost_size = 3;
	extend_size = n + 2*ghost_size;
	weno5_get_flux(params,extend_size,ghost_size,u_old,flux);

	for (i = ghost_size; i < n+ghost_size; ++i)
	{
	    vflux->dens_flux[i] = -lambda*(flux[0][i+1] - flux[0][i]);
	    vflux->momn_flux[0][i] = -lambda*(flux[1][i+1] - flux[1][i]);
	    vflux->momn_flux[1][i] = -lambda*(flux[2][i+1] - flux[2][i]);
	    vflux->momn_flux[2][i] = -lambda*(flux[3][i+1] - flux[3][i]);
	    vflux->engy_flux[i] = -lambda*(flux[4][i+1] - flux[4][i]);
            vflux->pdens_flux[0][i] = -lambda*(flux[5][i+1] - flux[5][i]);
            vflux->pdens_flux[1][i] = -lambda*(flux[6][i+1] - flux[6][i]);
	    if (isnan(vflux->dens_flux[i]))
	    {
		int j;
		for (j = 0; j < extend_size; ++j)
		    printf("u = %f %f %f %f %f\n",
			u_old[0][j],u_old[1][j],u_old[2][j],
			u_old[3][j],u_old[4][j]);
		clean_up(ERROR);
	    }
	}
}	/* end weno5_flux */

static void weno5_get_flux(
	POINTER params,
	int extend_size, 
	int ghost_size, 
	double **u_old,
        double **flux)
{
    	int i, j, k;
    	double u_mid[12];	// rho,mx,my,mz,e,u,v,w,p,u^2+v^2+w^2,a
        double rho0[2];         // pdens0,pdens1
    	double h; 	//total enthalpy h=0.5*(u^2+v^2+w^2) + (a^2)/(gma - 1)
    	double R[5][5], L[5][5];
    	double sten_u[6][5], sten_f[6][5]; //f_-2, f_-1, f_0, f_1, f_2, f_3
    	double gfluxp[5][5], gfluxm[5][5]; 
    	//double maxeig[5] = {0 ,0 ,0, 0, 0};
        double **maxeig;
	double v,u[9],f_tmp[5],ff[5]; //u[6]? u[5]?
	static double **f;
	static int max_n = 0;
	SCHEME_PARAMS *scheme_params = (SCHEME_PARAMS*)params;
	double a;
//	double *gamma;
	double *gamma = u_old[8];
        //double gamma = scheme_params->gamma;
    	//double a,gm = gamma - 1.0;

	for (i = 0; i < extend_size; i++)
	{
	    if (fabs(gamma[i]) < 1e-12)
	    {
		printf("Wrong gamma in WENO solver.\n");
		clean_up(ERROR);
	    }
//	    printf("gamma[%d] = %lf.\n", i, gamma[i]);
	}

	if (max_n < extend_size)
        {
            if (max_n != 0)
                FT_FreeThese(1,f);
            max_n = extend_size;
            FT_MatrixMemoryAlloc((POINTER*)&f,extend_size,5,sizeof(double));
        }
//        FT_VectorMemoryAlloc((POINTER*)&gamma,extend_size,sizeof(double));
        FT_MatrixMemoryAlloc((POINTER*)&maxeig,5,extend_size,sizeof(double));        

        //if(eqn_params->multi_comp_non_reactive == YES)
/*
        {
        int ii;
        double          w;
        double          cv, mgam[2];
        double          M[2];
        double          mole_dens;

        mgam[0] = 1.67;
        mgam[1] = 1.67;
        M[0] = 18.02;
        M[1] = 18.02;

        //NEED TO FIX. This is only for M EOS test
        for(i = 0; i < extend_size; ++i)
        {
            mole_dens = 0.0;
            cv = 0.0;
            for (ii = 0; ii < 2; ++ii)
            {
                w = u_old[6+ii][ii]/(M[ii]*u_old[0][ii]);
                mole_dens += w;
                cv += w/(mgam[ii] - 1.0);
            }
            gamma[i] = 1.0 + mole_dens/cv;
            //fprintf(stdout, "PRINT-HK gamma in WENO %e\n",gamma[i]);
        }
        }
*/
        /*
        for(i = 0; i < extend_size; ++i)
            gamma[i] = scheme_params->gamma;
        */

    	for(i = 0; i < extend_size; ++i)
    	{
	    //a = sqrt(gamma * u_old[5][i]/u_old[0][i]);
	    //v = u_old[1][i]/u_old[0][i];
	    //maxeig[0] = std::max(maxeig[0], fabs(v - a));
	    //maxeig[1] = std::max(maxeig[1], fabs(v));
	    //maxeig[4] = std::max(maxeig[4], fabs(v + a));
	    for (j = 0; j < 6; ++j)
                u[j] = u_old[j][i];
	    u2f(u,f[i]);
    	}
//	maxeig[2] = maxeig[1];
//	maxeig[3] = maxeig[1];

	for(i = ghost_size; i < extend_size - ghost_size +1; i++)	//Dan
	{
	    for (j = 0; j < 5; j++)
		maxeig[j][i] = 0;
	    for (j = i-ghost_size; j < i+ghost_size; j++)
	    {
		a = sqrt(gamma[i] * u_old[5][j]/u_old[0][j]);
		v = u_old[1][j]/u_old[0][j];
		maxeig[0][i] = std::max(maxeig[0][i], fabs(v - a));
		maxeig[1][i] = std::max(maxeig[1][i], fabs(v));
		maxeig[4][i] = std::max(maxeig[4][i], fabs(v + a));
	    }
	    maxeig[2][i] = maxeig[1][i];
	    maxeig[3][i] = maxeig[1][i];
	}

    	for(i = ghost_size; i < extend_size - ghost_size + 1; ++i)
    	{
	    /*** Get u_1/2 ***/

	    for(j = 0; j < 5; ++j)
	    {
	    	u_mid[j] = 0.5*(u_old[j][i-1] + u_old[j][i]);
	    }
            u_mid[11] = 0.5*(gamma[i-1] + gamma[i]);
            //u_mid[11] = 0.5*(u_old[8][i-1] + u_old[8][i]);
	    u_mid[5] = u_mid[1]/u_mid[0];
	    u_mid[6] = u_mid[2]/u_mid[0];
	    u_mid[7] = u_mid[3]/u_mid[0];
	    u_mid[9] = sqr(u_mid[5]) + sqr(u_mid[6]) + sqr(u_mid[7]);
            u_mid[8] = (u_mid[4] - 0.5*u_mid[0]*u_mid[9])*(u_mid[11] - 1.0); //u_mid[11] = gamma and u_old[8] = gamma
	    u_mid[10] = sqrt(u_mid[11]*u_mid[8]/u_mid[0]); //u_mid[11] = gamma and u_old[8] = gamma
 
            rho0[0] = (0.5*(u_old[6][i-1] + u_old[6][i]))/u_mid[0]; //mass fraction for pdens0
            rho0[1] = (0.5*(u_old[7][i-1] + u_old[7][i]))/u_mid[0]; //mass fraction for pdens1

	    /*** R(u_1/2) & R^-1(u_1/2) ***/

	    h = 0.5*u_mid[9] + u_mid[11]*u_mid[8]/(u_mid[11] - 1.0)/u_mid[0]; //u_mid[11] = gamma and u_old[8] = gamma

	    R[0][0] = 1.0;
	    R[0][1] = 1.0;
	    R[0][2] = 0.0;
	    R[0][3] = 0.0;
	    R[0][4] = 1.0;

	    R[1][0] = u_mid[5] - u_mid[10];
	    R[1][1] = u_mid[5];
	    R[1][2] = 0.0;
	    R[1][3] = 0.0;
	    R[1][4] = u_mid[5] + u_mid[10];

	    R[2][0] = u_mid[6];
	    R[2][1] = u_mid[6];
	    R[2][2] = 1.0;
	    R[2][3] = 0.0;
	    R[2][4] = u_mid[6];

	    R[3][0] = u_mid[7];
	    R[3][1] = u_mid[7];
	    R[3][2] = 0.0;
	    R[3][3] = 1.0;
	    R[3][4] = u_mid[7];

	    R[4][0] = h - u_mid[5]*u_mid[10];
	    R[4][1] = 0.5*u_mid[9];
	    R[4][2] = u_mid[6];
	    R[4][3] = u_mid[7];
	    R[4][4] = h + u_mid[5]*u_mid[10];;

	    L[0][0] = h + u_mid[10]*(u_mid[5] - u_mid[10]) / (u_mid[11] - 1.0); //gm = gamma -1.0 = u_mid[11] - 1.0;
	    L[0][1] = -1.0*(u_mid[5] + u_mid[10] / (u_mid[11] - 1.0)); //gm = gamma -1.0 = u_mid[11] - 1.0;
	    L[0][2] = -1.0*u_mid[6];
	    L[0][3] = -1.0*u_mid[7];
	    L[0][4] = 1.0;

	    L[1][0] = -2.0*h + 4.0*sqr(u_mid[10])/(u_mid[11] - 1.0); //gm = gamma -1.0 = u_mid[11] - 1.0;
	    L[1][1] = 2.0*u_mid[5];
	    L[1][2] = 2.0*u_mid[6];
	    L[1][3] = 2.0*u_mid[7];
	    L[1][4] = -2.0;

	    L[2][0] = -2.0*u_mid[6]*sqr(u_mid[10])/ (u_mid[11] - 1.0); //gm = gamma -1.0 = u_mid[11] - 1.0;
	    L[2][1] = 0;
	    L[2][2] = 2.0*sqr(u_mid[10])/(u_mid[11] - 1.0); //gm = gamma -1.0 = u_mid[11] - 1.0;
	    L[2][3] = 0;
	    L[2][4] = 0;

	    L[3][0] = -2.0*u_mid[7]*sqr(u_mid[10])/ (u_mid[11] - 1.0); //gm = gamma -1.0 = u_mid[11] - 1.0;
	    L[3][1] = 0;
	    L[3][2] = 0;
	    L[3][3] = 2.0*sqr(u_mid[10])/ (u_mid[11] - 1.0); //gm = gamma -1.0 = u_mid[11] - 1.0;
	    L[3][4] = 0;

	    L[4][0] = h - u_mid[10]*(u_mid[5] + u_mid[10])/ (u_mid[11] - 1.0); //gm = gamma -1.0 = u_mid[11] - 1.0;
	    L[4][1] = -1.0*u_mid[5] + u_mid[10]/ (u_mid[11] - 1.0); //gm = gamma -1.0 = u_mid[11] - 1.0;
	    L[4][2] = -1.0*u_mid[6];
	    L[4][3] = -1.0*u_mid[7];
	    L[4][4] = 1.0;

	    for(j = 0; j < 5; ++j)
	    for(k = 0; k < 5; ++k)
		L[j][k] *= (u_mid[11] - 1.0)/(2.0*sqr(u_mid[10]));  //gm = gamma -1.0 = u_mid[11] - 1.0;

	    /*** Get R^-1 * u and R^-1 * F ***/	    

	    for(j = 0; j < 6; ++j)
	    {
		for (k = 0; k < 5; ++k)
                    u[k] = u_old[k][i - ghost_size + j];
	    	matmvec(sten_u[j],L,u);
	    	matmvec(sten_f[j],L,f[i - ghost_size + j]);
	    }

	    for(j = 0; j < 5; ++j)
	    for(k = 0; k < 5; ++k)
	    {
		gfluxp[j][k] = 0.5*(sten_f[k][j] + maxeig[j][i]*sten_u[k][j]);
		gfluxm[j][k] = 0.5*(sten_f[5 - k][j] - maxeig[j][i]*
				sten_u[5 - k][j]);
	    }

	    for(j = 0; j < 5; ++j)
	    {
	    	f_tmp[j] = weno5_scal(gfluxp[j]);
	    	f_tmp[j] += weno5_scal(gfluxm[j]);
	    }

	    matmvec(ff, R, f_tmp);
	    for(j = 0; j < 5; ++j)
	    {
                flux[j][i] = ff[j];
		if (isnan(flux[j][i]))
		{
		    (void) printf("In weno5_flux(): f[%d][%d] = %f\n",j,i,
					f[j][i]);
		    for (i = 0; i < extend_size; ++i)
                        printf("u[%d] = %f %f %f %f %f\n",k,u_old[0][i],
				u_old[1][i],u_old[2][i],u_old[3][i],
				u_old[4][i]);
                    clean_up(ERROR);
		}
	    }
            flux[5][i] = rho0[0] * flux[0][i];
            flux[6][i] = rho0[1] * flux[0][i];
    	}
        FT_FreeThese(1,maxeig);
}

static double weno5_scal(double *f)
{
    	int i, j;
    	const double eps = 1.e-8;
    	const int p = 2;

    	double c[3] = {0.1, 0.6, 0.3}; //*** Optimal weights C_k 
    	double is[3]; //*** a smoothness measurement of the flux function
    	double alpha[3];
    	double omega[3]; // weights for stencils
    	double sum;
    	double a[3][5] = {{1.0/3, -7.0/6, 11.0/6, 0, 0}, 
			  {0, -1.0/6, 5.0/6, 1.0/3, 0}, 
			  {0, 0, 1.0/3, 5.0/6, -1.0/6}}; 
		//*** coefficients for 2nd-order ENO interpolation stencil
    	double w[5]; //weight for every point
    	double flux = 0.0;

    	f2is(f,is);

    	sum = 0.0;

    	for(i = 0; i < 3; ++i)
    	{
	    alpha[i] = c[i]/pow(eps + is[i],p);
	    sum += alpha[i];
    	}
    	for(i = 0; i < 3; ++i)
    	{
	    omega[i] = alpha[i] / sum;
    	}

    	for(i = 0; i < 5; ++i)
    	{
	    w[i] = 0.0;
	    for(j = 0; j < 3; ++j)
	    {
	    	w[i] += omega[j] * a[j][i];
	    }
    	}
    	for(i = 0; i < 5; ++i)
    	{
	    flux += w[i] * f[i];
    	}  	  
    	return flux;
}

static void f2is(
	double *f, 
	double *s)
{
	s[0] = 13.0/12*sqr((f[0] - 2.0*f[1] + f[2])) +
                0.25*sqr((f[0] - 4.0*f[1] + 3.0*f[2]));
        s[1] = 13.0/12*sqr((f[1] - 2.0*f[2] + f[3])) +
                0.25*sqr((f[1] - f[3]));
        s[2] = 13.0/12*sqr((f[2] - 2.0*f[3] + f[4])) +
                0.25*sqr((3.0*f[2] - 4.0*f[3] + f[4]));
}

static void matmvec(
	double *b, 
	double L[5][5], 
	double *x)
{
    	int i, j;

    	for(i = 0; i < 5; ++i)
    	{
	    b[i] = 0.0;
	    for(j = 0; j < 5; ++j)
	    {
	    	b[i] += L[i][j] * x[j]; 
	    }
    	}
}

static void u2f(
	double *u,
        double *f)
{
	double v = u[1]/u[0];

    	f[0] = u[1];
    	f[1] = v*u[1] + u[5];
    	f[2] = v*u[2];
    	f[3] = v*u[3];
    	f[4] = v*(u[4] + u[5]);
}
