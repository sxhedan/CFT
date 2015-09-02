#include <cFluid.h>

#if !defined(max)
#  define     max(a,b)     (((a) > (b)) ? (a) : (b))
#endif /* !defined(max) */
#define floor_pressure(p,pmin)	max((p),(pmin))
#define Godunov_pressure(pl,vxl,ml,pr,vxr,mr)                              \
    ((((vxl) - (vxr))*(mr)*(ml) + (pr)*(ml) + (pl)*(mr))/ ((ml) + (mr)))
#define vacuum_produced(pstar,ul_mid,ur_mid)                            \
    ((pstar <= 2.0*p_min) /*TOLERANCE*/ && (ul_mid < ur_mid))

LOCAL bool secant_find_mid_state(STATE*,STATE*,double,double*,double*,double*,double*,double*,
				  double*,RIEMANN_SOLVER_WAVE_TYPE*,RIEMANN_SOLVER_WAVE_TYPE*,
				  double,double,double,double,double,double);
LOCAL bool secant_press_conv(double,double,double,double);
LOCAL bool secant_vel_conv(double,double,double,double,double,double);
LOCAL bool invalid_shock(double,double,double,double);

extern bool find_mid_state(
	STATE				*stl,
	STATE				*str,
	double				pjump,
	double				*ppstarl,
	double				*ppstarr,
	double				*pustarl,
	double				*pustarr,
	double				*pml,
	double				*pmr,
	RIEMANN_SOLVER_WAVE_TYPE	*pl_wave,
	RIEMANN_SOLVER_WAVE_TYPE	*pr_wave)
{
	double		vlmax, vrmin;
	double		pr, pl;
	double		p_start, p_min;
	double		eps_u,eps_p;		/* convergence limits for ustar, pstar */
	double		pmax, pmin;
	bool		status;
	const double	SECANT_EPS = 1e-10;	/* TOLERANCE */

	//Initial guess using the linearized Godunov equation
	initialize_riemann_solver(stl,str,&p_start,&p_min,
	                          SECANT_EPS,&eps_u,&eps_p);

	//Check for vacuum in solution
	pr = p_min;
	pl = floor_pressure(p_min+pjump,p_min);
	vrmin = str->vel[0] + riemann_wave_curve(str,pr);
	vlmax = stl->vel[0] - riemann_wave_curve(stl,pl);

	if (vlmax < vrmin)
	{
	    //wave curves do not intersect,  a vacuum is formed
	    *ppstarl = pl;
	    *pustarl = vlmax;
 	    *pml = mass_flux(pl,stl);
	    *pl_wave = (pl >= stl->pres) ? SHOCK : RAREFACTION;

	    *ppstarr = pr;
	    *pustarr = vrmin;
 	    *pmr = mass_flux(pr,str);
	    *pr_wave = (pr >= str->pres) ? SHOCK : RAREFACTION;

	    return FUNCTION_SUCCEEDED;
	}

   	pmax = HUGE_VAL;
	pmin = p_min;

	status = secant_find_mid_state(stl,str,pjump,ppstarl,ppstarr,
				       pustarl,pustarr,pml,pmr,pl_wave,
				       pr_wave,p_start,p_min,eps_u,eps_p,
				       pmin,pmax);

	return status;
}

extern	void midstate(
	STATE				*ahead,
	STATE				*ans,
	double				m,
	double				um,
	double				pm,
	RIEMANN_SOLVER_WAVE_TYPE	wave,
	int				sgn)
{
	double pa;

	switch (wave) 
	{
	case SHOCK:
	    pa = ahead->pres;
	    ans->dim = ahead->dim;
	    ans->eos = ahead->eos;
	    ans->dens = sgn*m/(um - ahead->vel[0] + sgn*m/ahead->dens);
	    if (invalid_shock(pa,ahead->dens,pm,ans->dens))
		ans->dens = dens_Hugoniot(pm,ahead);
	    ans->vel[0] = um;
	    ans->pres = pm;
            ans->pdens[0] = (ahead->pdens[0]/ahead->dens)*ans->dens;
            ans->pdens[1] = (ahead->pdens[1]/ahead->dens)*ans->dens;
	    break;

	case RAREFACTION:
	    state_on_adiabat_with_pr(ahead,pm,ans);
	    ans->vel[0] = um;
	    break;

	default:
	    screen("ERROR: unknown wave %d in midstate()\n",wave);
	    clean_up(ERROR);
	    break;
	}

}	/*end midstate*/

LOCAL bool secant_find_mid_state(
	STATE				*stl,	/* left TGas or VGas states */
	STATE				*str,	/* right TGas or VGas states */
	double				pjump,	/* Pressure jump across contact, pl-pr*/
	double				*ppstarl,	/* left middle pressure */
	double				*ppstarr,	/* right middle pressure */
	double				*pustarl,	/* left middle velocity */
	double				*pustarr,	/* right middle velocity */
	double				*pml,
	double				*pmr,
	RIEMANN_SOLVER_WAVE_TYPE	*pl_wave,
	RIEMANN_SOLVER_WAVE_TYPE	*pr_wave, /* wave families */
	double				p_start,
	double				p_min,
	double				eps_u,	/*convergence limits for ustar, pstar*/
	double				eps_p,	/*convergence limits for ustar, pstar*/
	double				pmin,  /*solution lies in range */
	double				pmax)  /*pmin <= p <= pmax      */
{
	const double meps = 10.0*MACH_EPS;/*TOLERANCE*/
	int	    i;
	double	    pstarr, pstarl;
	double	    ml,	mr;
	double	    ul_mid, ur_mid, ul_mid_old, ur_mid_old, p_old, p_old2;
	double	    gul_mid, gur_mid;
	double	    vul, vur;
	double	    numer, denom;
	double	    pl, pr, vxl, vxr;
	double	    poj, plj;
	double	    du_p_start, du_pstar, du_ppstar;
	bool	    status;
	const double MIN_DENOM = 1.0e-10;	//TOLERANCE
	const int   NUM_SEC_ITER = 160;		//Max # secant iterations

	p_old = p_start;

	poj = floor_pressure(p_old+pjump,p_min);
 	ml = mass_flux(poj,stl);
	pl = stl->pres;	
	vxl = stl->vel[0];		

 	mr = mass_flux(p_old,str);
	pr = str->pres;
	vxr = str->vel[0];

	plj = floor_pressure(pl-pjump,p_min);
	pstarr = Godunov_pressure(plj,vxl,ml,pr,vxr,mr);
	if (pstarr <= p_min)
	    pstarr = p_min;
	pstarl = floor_pressure(pstarr+pjump,p_min);
	ul_mid = vxl;
	if (fabs(ml) > MIN_DENOM)
	    ul_mid += (plj - pstarl)/(ml);
	ur_mid = vxr;
	if (fabs(mr) > MIN_DENOM)
	    ur_mid += (pstarr - pr)/(mr);

	//   secant method

	for (i = 0; i < NUM_SEC_ITER; ++i)
	{
	    ur_mid_old = ur_mid;
	    ul_mid_old = ul_mid;
	    p_old2 = p_old;
	    p_old = pstarr;
	    poj = floor_pressure(p_old+pjump,p_min);
	    ur_mid = vxr + riemann_wave_curve(str,p_old);
	    ul_mid = vxl - riemann_wave_curve(stl,poj);
	    if (ur_mid <= ul_mid)
	    {
		if (pmin < p_old)
		    pmin = p_old;
	    }
	    else
	    {
		if (p_old < pmax)
		    pmax = p_old;
	    }
	    denom = (ur_mid - ur_mid_old) - (ul_mid - ul_mid_old);
	    numer = (ur_mid - ul_mid)*(p_old - p_old2);
	    if (fabs(denom) >= MIN_DENOM &&
	        fabs(numer) > 0.5*eps_p*p_old*fabs(denom))
	    {
	        pstarr = p_old - numer/denom;
	        if (pstarr <= p_min)
	            pstarr = p_min;
	        pstarl = floor_pressure(pstarr+pjump,p_min);
	    }

	    if (secant_press_conv(pstarr,p_old,eps_p,meps))
	    {
	        if (secant_vel_conv(ul_mid,ul_mid_old,ur_mid,
	                            ur_mid_old,eps_u,meps))
	        {
	            *ppstarl = pstarl;
	            *ppstarr = pstarr;
	            *pustarl = *pustarr = 0.5*(ur_mid + ul_mid);
 		    *pml = mass_flux(pstarl,stl);
 		    *pmr = mass_flux(pstarr,str);
	            *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	            *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
	            return FUNCTION_SUCCEEDED;
	        }
	        else if (vacuum_produced(pstarr,ul_mid,ur_mid))
	        {
	            // Vacuum produced
	            *ppstarr = pstarr = p_min;
	            *ppstarl = pstarl = floor_pressure(p_min+pjump,p_min);
	            *pustarl = ul_mid;
	            *pustarr = ur_mid;
 		    *pml = mass_flux(pstarl,stl);
 		    *pmr = mass_flux(pstarr,str);
	            *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	            *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
	            return FUNCTION_SUCCEEDED;
	        }
	    }
	}

	printf("Does not converge in Riemann solver!.\n");
	return FUNCTION_FAILED;
}

LOCAL	bool	secant_press_conv(
	double	pstar,
	double	p_old,
	double	eps_p,
	double	meps)
{
	double	eps = eps_p*pstar;

	eps = max(eps,meps);
	return (fabs(pstar-p_old) <= eps) ? YES : NO;
}

LOCAL	bool secant_vel_conv(
	double ulm,
	double ulm_old,
	double urm,
	double urm_old,
	double eps_u,
	double meps)
{
	double	eps, eps_l, eps_r;

	eps = 0.5*(fabs(ulm)+fabs(urm))*eps_u;	eps = max(eps,meps);
	if (fabs(ulm-urm) <= eps)
	    return YES;

	eps_l = fabs(ulm)*eps_u;	eps_l = max(eps_l,meps);
	eps_r = fabs(urm)*eps_u;	eps_r = max(eps_r,meps);

	return ((fabs(ulm-ulm_old) <= eps_l) && (fabs(urm-urm_old) <= eps_r)) ?
	    YES : NO;
}

LOCAL	bool invalid_shock(
	double pa,
	double rhoa,
	double pb,
	double rhob)
{
	return (((pb - pa)*(rhob - rhoa) < 0.0) || (rhob < 0.0)) ? YES : NO;
}
