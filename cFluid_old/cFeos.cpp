#include <cFluid.h>

#define R 1.25
#define EPS 1e-6

#if !defined(max)
#  define     max(a,b)     (((a) > (b)) ? (a) : (b))
#endif /* !defined(max) */

#if !defined(min)
#define     min(a,b)     (((a) < (b)) ? (a) : (b))
#endif /* !defined(min) */

extern double EosPressure(
	STATE *state)
{
	double 		dens = state->dens;
	double 		engy = state->engy;
	double 		*momn = state->momn;
	double 		ke,pres;
	int		i;
	int		dim = state->dim;
	EOS_PARAMS	*eos = state->eos;
	double		gamma = eos->gamma;

	if (dens <= 0.0)
	    return 0.0;
	ke = 0.0;
	for (i = 0; i < dim; ++i)
	    ke += sqr(momn[i]);
	ke *= 0.5/dens;
	pres = (gamma - 1.0)*(engy - ke + dens*eos->einf) - gamma*eos->pinf;
	
	return pres;
}	/* end EosPressure */

extern double EosCV(STATE *state)
{
        double          dens = state->dens;
        double          pres = state->pres;
        double          cv;
        double          gamma = state->eos->gamma;

        cv = R/(gamma-1);

        return cv;
}       /* end EosCV */

extern double  EosCP(STATE *state)
{
        double          cp;
        double          x;
        double          dens = state->dens;
        double          pres = state->pres;
        double          gamma = state->eos->gamma;
        EOS_PARAMS      *eos = state->eos;
        double          cv;


        cv = R/(gamma-1);
        cp = gamma*cv;
        if (eos->et!=0.0)
        {
           x = eos->et*(gamma-1)*(gamma-1)*pow(dens,gamma-1);
           cp = 1.0 + x/( (pres + eos->pinf)/dens -x );
        }

        return cp;
}       /* end EosCP */

extern double EosTemperature(STATE *state)
{
        double          dens = state->dens;
        double          pres = state->pres;
        double          temp;
        double          gamma = state->eos->gamma;
        EOS_PARAMS      *eos = state->eos;
        double          cv;

        cv = R/(gamma-1);
        if (dens <= 0.0)
            return 0.0;
        temp = (pres + eos->pinf)/(R*dens);
        if (eos->et!=0)
            temp = eos->et*pow(dens, gamma-1)/cv;

        return temp;
}       /* end EosTemperature */

extern double EosSoundSpeedSqr(
	STATE *state)
{
	double		pres = state->pres;
	double		dens = state->dens;
	EOS_PARAMS	*eos = state->eos;
	
	return eos->gamma*(pres + eos->pinf)/dens;
}

extern double EosSoundSpeed(
	STATE *state)
{
	return sqrt(EosSoundSpeedSqr(state));
}

extern double EosInternalEnergy(
	STATE *state)
{
	double		pres = state->pres;
	double		dens = state->dens;
	EOS_PARAMS	*eos = state->eos;
	double		gamma = eos->gamma;

	return (pres+gamma*eos->pinf)/(gamma-1) - dens*eos->einf;
}

extern double EosEnergy(
	STATE *state)
{
	int	i,dim = state->dim;
	double	dens = state->dens;
	double	*momn = state->momn;
	double	e;
	
	e = 0.0;
	for (i = 0; i < dim; ++i)
	    e += 0.5*sqr(momn[i])/dens;
	e += EosInternalEnergy(state);

	return e;
}

extern double EosMaxBehindShockPres(
	double M2,		// Mach number sqaure
	STATE *state)
{
	double im2;		// acoustic impedance squared
	double gamma = state->eos->gamma;
	double dens = state->dens;
	double pres = state->pres;
	double c4,c5,p1;

	im2 = M2*gamma*pres*dens;
	c4 = (gamma -1.0)/(gamma + 1.0);
	c5 = 2.0/(gamma + 1.0);
	p1 = c5*im2/dens - c4*pres;
	return p1;
}	/* end EosMaxBehindShockPres */

extern void CovertVstToState(
	STATE		*state,
	SWEEP		*vst,
	EOS_PARAMS	*eos,
	int		ind,
	int		dim)
{
	int	i;

	state->dim = dim;
	state->eos = eos;
	state->dens = vst->dens[ind];
	state->engy = vst->engy[ind];
	for (i = 0; i < dim; ++i)
	    state->momn[i] = vst->momn[i][ind];
	state->pres = EosPressure(state);
}

extern void EosSetTVDParams(
	SCHEME_PARAMS	*scheme_params,
	EOS_PARAMS	*eos)
{
	scheme_params->gamma = eos->gamma;
	scheme_params->einf = eos->einf;
}

extern double GAM(STATE *st)
{
	double gamma = st->eos->gamma;
	return gamma - 1.0;
}

extern double Coef1(STATE *st)
{
	double	gamma = st->eos->gamma;
	return 0.5*(gamma + 1.0);
}

extern double Coef2(STATE *st)
{
	double gamma = st->eos->gamma;
	return 0.5*(gamma - 1.0);
}

extern double Coef3(STATE *st)
{
	double gamma = st->eos->gamma;
	return 0.5*(gamma - 1.0)/gamma;
}

extern double Coef4(STATE *st)
{
	double gamma = st->eos->gamma;
	return (gamma - 1.0)/(gamma+1);
}

extern void initialize_riemann_solver(
	STATE	*stl,
	STATE	*str,
	double	*pstar,
	double	*p_min,
	double	eps,
	double	*eps_u,
	double	*eps_p)
{
	double pl, pr;
	double cl, cr, ul_tdl, ur_tdl, z;
	double c2l, c2r, c3l, c3r;
	double dutdl;
	double vl, vr;

	vl = stl->vel[0];
	vr = str->vel[0];
	pl = stl->pres;
	pr = str->pres;

	*eps_u = *eps_p = eps;
	*p_min = 1e-12;

	    /*
	     * Setting pstar = 0.5*(pl + pr) at an SPOLY--POLY contact can
	     * yield pstar < 0 for situations in which cavitation DOES NOT
	     * occur (that is, vacuum DOES NOT form at contact).
	     * In this case, the (negative) initial guess for pstar may be
	     * passed to POLY_mass_flux(), from secant_find_mid_state() or
	     * godunov_find_mid_state(), which should never happen.
	     */

	*pstar = 0.5*(pl + pr);
	*pstar = max(*pstar,*p_min);
	return;
}

extern double riemann_wave_curve(
	STATE	*st0,
	double	pstar)
{
	double rho0, p0;
	double c1, c2, c3;

	rho0 = st0->dens;
	p0 = st0->pres;

	if (pstar < 1e-12)
	    pstar = 1e-12;

	c1 = Coef1(st0);
	c2 = Coef2(st0);
	c3 = Coef3(st0);

	return (pstar < p0) ?
	        EosSoundSpeed(st0)*(pow(pstar/p0,c3) - 1.0)/ c2 :
	        (pstar-p0)/sqrt(rho0*(c1*pstar+c2*p0));
}

extern double acoustic_impedance_squared(STATE *st)
{
	double gamma, i2;

	gamma = st->eos->gamma;
	i2 =  gamma * st->pres * st->dens;

	return i2;
}

extern double acoustic_impedance(STATE *st)
{
	return sqrt(acoustic_impedance_squared(st));
}

extern double mass_flux(
	double	p,
	STATE	*st0)
{
	double p0, rho0;
	double xi, m, i0;

	p0 = st0->pres;
	rho0 = st0->dens;
	if (p < p0)
	{
	    i0 = acoustic_impedance(st0);
	    xi = (p > 0.0) ? p/p0 : 0.0;
	    if ( (1.0 - xi) < EPS)
	    	return i0;
	    m = Coef3(st0)*(1.0-xi)/(1.0-pow(xi,Coef3(st0)));
	    return i0*m;
	}
	else
	    return sqrt(rho0*(Coef1(st0)*p + Coef2(st0)*p0));
}

extern double dens_Hugoniot(
	double	p1,
	STATE	*st0)
{
	double p0, c4;

	p0 = st0->pres;
	c4 = Coef4(st0);
	return st0->dens*(p1 + p0*c4)/(p0 + p1*c4);
}

extern void state_on_adiabat_with_pr(
	STATE	*st0,
	double	p1,
	STATE	*st1)
{
    	int i;

	st1->dim = st0->dim;
	st1->eos = st0->eos;
    	for (i = 0; i < st0->dim; i++)
	    st1->vel[i] = 0.0;
	st1->dens = st0->dens*pow(p1/st0->pres,1.0/st0->eos->gamma);
	st1->pres = p1;
}
