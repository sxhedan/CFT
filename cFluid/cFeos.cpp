#include <cFluid.h>

LOCAL double Gamma(STATE*);
LOCAL double R(STATE*);

//#define R 1.25
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
        double          dens = state->dens;
        double          engy = state->engy;
        double          *momn = state->momn;
        double          ke,pres;
        int             i;
        int             dim = state->dim;
        EOS_PARAMS      *eos = state->eos;
        double          gamma = eos->gamma;
        int             eos_type = eos->eos_type;
        double          re;

        switch (eos_type)
        {
        case    POLYTROPIC:
        case    MULTI_COMP_POLYTROPIC:
            if (dens <= 0.0)
                pres = 0.0;
            else
            {
                ke = 0.0;
                for (i = 0; i < dim; ++i)
                    ke += sqr(momn[i]);
                ke *= 0.5/dens;
                pres = (Gamma(state) - 1.0)*(engy - ke);
            }
            break;
        case    STIFFENED_POLYTROPIC:
            if (dens <= 0.0)
                pres = 0.0;
            else
            {
                ke = 0.0;
                for (i = 0; i < dim; ++i)
                    ke += sqr(momn[i]);
                ke *= 0.5/dens;
                pres = (Gamma(state) - 1.0)*(engy - ke + dens*eos->einf) - Gamma(state)*eos->pinf;
            }
            break;
        default:
            screen("ERROR in EosPressure() "
                   "Unknown equation of state\n");
            clean_up(ERROR);
            break;
        }

        return pres;
}       /* end EosPressure */

/* Terms added by PRAO */
extern double EosCV(STATE *state)
{
        double          dens = state->dens;
        double          pres = state->pres;
        double          cv;
        double          gamma = state->eos->gamma;

        cv = R(state)/(Gamma(state)-1.0);

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


        cv = R(state)/(Gamma(state)-1);
        cp = Gamma(state)*cv;
        if (eos->et!=0.0)
        {
           x = eos->et*(Gamma(state)-1)*(Gamma(state)-1)*pow(dens,Gamma(state)-1);
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

        cv = R(state)/(Gamma(state)-1);
        if (dens <= 0.0)
            return 0.0;
        temp = (pres + eos->pinf)/(R(state)*dens);
        if (eos->et!=0)
            temp = eos->et*pow(dens, Gamma(state)-1)/cv;

        return temp;
}       /* end EosTemperature */

/* end of terms added by PRAO */

extern double EosSoundSpeedSqr(
	STATE *state)
{
	double		pres = state->pres;
	double		dens = state->dens;
	EOS_PARAMS	*eos = state->eos;
        int             eos_type = eos->eos_type;
        double          soundspeedsqr;

        switch (eos_type)
        {
        case    POLYTROPIC:
        case    MULTI_COMP_POLYTROPIC:
            //return eos->gamma*(pres + eos->pinf)/dens;
            soundspeedsqr = Gamma(state)*pres/dens;
            break;
        case    STIFFENED_POLYTROPIC:
            soundspeedsqr = Gamma(state)*EosPressure(state)/dens; //same as Gamma(state)*pres/dens NEED TO CHECK
            break;
        default:
	    screen("ERROR in EosSoundSpeedSqr() "
	    "Unknown equation of state\n");
	    clean_up(ERROR);
            break;
        }
      
        return soundspeedsqr;	
}

extern double EosSoundSpeed(
	STATE *state)
{
	return sqrt(EosSoundSpeedSqr(state));
}

extern double EosInternalEnergy(
	STATE *state)
{
        double          pres = state->pres;
        double          dens = state->dens;
        EOS_PARAMS      *eos = state->eos;
        double          gamma = eos->gamma;
        int             eos_type = eos->eos_type;
        double          internalenergy;

        switch (eos_type)
        {
        case    POLYTROPIC:
        case    MULTI_COMP_POLYTROPIC:
            internalenergy = pres/(Gamma(state)-1.0);  //TGAS NEED TO CHECK
            break;
        case    STIFFENED_POLYTROPIC:
            internalenergy = (pres+Gamma(state)*eos->pinf)/(Gamma(state)-1) - dens*eos->einf;
            break;
        default:
            screen("ERROR in EosInternalEnergy() "
                   "Unknown equation of state\n");
            clean_up(ERROR);
            break;
        }

	return internalenergy;
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

extern double EosEnthalpyDifference(
        STATE *state)
{
        /* H_h = e_h + p/rho, where H_h = partial specific enthalpy, e_h = specific internal energy, p = pressure and rho= density. In 2008 Subgrid paper by Hyunkyung, in eqn (7), the parabolic term contains (H_h-H_l), which is equal to (e_h-e_l). This is what is being calculated here */

        double se1, se2, se;
        double dens = state->dens;
        double mgamma0 = state->eos->mgamma[0];
        double mgamma1 = state->eos->mgamma[1];

        int    eos_type = state->eos->eos_type;

        switch (eos_type)
        {
        case MULTI_COMP_POLYTROPIC:
             se1 = EosPressure(state)/(dens*(mgamma0-1.0));
             se2 = EosPressure(state)/(dens*(mgamma1-1.0));
             se = se1 - se2;
             break;
        default:
            screen("ERROR in EosEnthalpyDifference() "
                   "Unknown equation of state\n");
            clean_up(ERROR);
            break;
        }
        return se;

} /* added by PRAO */


extern double EosMaxBehindShockPres(
        double M2,              // Mach number sqaure
        STATE *state)
{
        double im2;             // acoustic impedance squared
        double gamma = state->eos->gamma;
        double dens = state->dens;
        double pres = state->pres;
        double c4,c5,p1;
        int    eos_type = state->eos->eos_type;

        switch (eos_type)
        {
        case    POLYTROPIC:
        case    MULTI_COMP_POLYTROPIC:
        case    STIFFENED_POLYTROPIC:
            im2 = M2*Gamma(state)*pres*dens;
            c4 = (Gamma(state) -1.0)/(Gamma(state) + 1.0);
            c5 = 2.0/(Gamma(state) + 1.0);
            p1 = c5*im2/dens - c4*pres;
            break;
        default:
	    p1 = 0;
//            screen("ERROR in EosMaxBehindShockPres() "
//                   "Unknown equation of state\n");
//            clean_up(ERROR);
            break;
        }

        return p1;
}       /* end EosMaxBehindShockPres */

extern double EosGruneisenGamma(
        STATE *state)
{
        return Gamma(state) - 1.0;
}       /* end EosGruneisenGamma */

LOCAL double Gamma(
        STATE *state)
{
        EOS_PARAMS      *eos = state->eos;
        double          gamma = eos->gamma;
        int             eos_type = eos->eos_type;
        double          dens = state->dens;
        double          *pdens = state->pdens;
        double          gam;
        int             i, nc = eos->n_comps;
        double          w;
        double          cv, *mgam = eos->mgamma;
        double          *M = eos->M;
        double          mole_dens;

        switch (eos_type)
        {
        case    POLYTROPIC:
        case    STIFFENED_POLYTROPIC:
            gam = gamma;
            break;
        case    MULTI_COMP_POLYTROPIC:
            mole_dens = 0.0;
            cv = 0.0;
            for (i = 0; i < nc; ++i)
            {
                w = pdens[i]/(M[i]*dens);
                mole_dens += w;
                cv += w/(mgam[i] - 1.0);
            }
            gam = 1.0 + mole_dens/cv;
            break;
        default:
	    gam = 0;
//            screen("ERROR in Gamma() "
//                   "Unknown equation of state\n");
//            clean_up(ERROR);
            break;
        }
        return gam;
}       /* end Gamma */

extern double R(
        STATE *state)
{
        EOS_PARAMS      *eos = state->eos;
        int             eos_type = eos->eos_type;
        double          dens = state->dens;
        double          *pdens = state->pdens;
        double          r;
        int             i, nc = eos->n_comps;
        double          n;
        double          *M = eos->M;

        switch (eos_type)
        {
        case    POLYTROPIC:
        case    STIFFENED_POLYTROPIC:
            //r = 1.0; // FIX ME?
            r = eos->R;
            break;
        case MULTI_COMP_POLYTROPIC:
            n = 0.0;
            for (i = 0; i < nc; ++i)
                n += pdens[i]/M[i];
            r = eos->R/(dens/n);
            break;
        default:
	    r = 0;
//            screen("ERROR in R() "
//                   "Unknown equation of state\n");
//            clean_up(ERROR);
            break;
        }

        return r;
}

extern double EosGamma(
        STATE *state)
{
        return Gamma(state);
}


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
        int nc = eos->n_comps;
        //if(eqn_params->multi_comp_non_reactive == YES)
        {
            int ii;
            for(ii = 0; ii < nc; ii++)
            //for(ii = 0; ii < 2; ii++)
            {
                state->pdens[ii] = vst->pdens[ii][ind];
            }
        }
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

extern double GAM(STATE *state)
{
	//double gamma = st->eos->gamma;
	return Gamma(state) - 1.0;
}

extern double Coef1(STATE *state)
{
	//double	gamma = st->eos->gamma;
	return 0.5*(Gamma(state) + 1.0);
}

extern double Coef2(STATE *state)
{
	//double gamma = st->eos->gamma;
	return 0.5*(Gamma(state) - 1.0);
}

extern double Coef3(STATE *state)
{
	//double gamma = st->eos->gamma;
	return 0.5*(Gamma(state) - 1.0)/Gamma(state);
}

extern double Coef4(STATE *state)
{
	//double gamma = st->eos->gamma;
	return (Gamma(state) - 1.0)/(Gamma(state)+1);
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
	/*this part is working well now for RT simulation, but is not tested yet for other problems*/
	/*FIXME*/
//	if (stl->eos != str->eos)
//	{
//#if defined(UNRESTRICTED_THERMODYNAMICS)
//	    *p_min = -HUGE_VAL;
//#else //defined(UNRESTRICTED_THERMODYNAMICS) 
//	    *p_min = max(Min_pressure(Tsl),Min_pressure(Tsr));
//	    *p_min = max(stl->min_pres, str->min_pres);
	    *p_min = 1e-12;
//#endif //defined(UNRESTRICTED_THERMODYNAMICS) 

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
//	}
/*
	c2l = Coef2(stl);
	c3l = Coef3(stl);
	c2r = Coef2(str);
	c3r = Coef3(str);
//	cl = sound_speed(Tsl);
//	cr = sound_speed(Tsr);
	cl = EosSoundSpeed(stl);
	cr = EosSoundSpeed(str);
	ul_tdl = vl + cl/c2l;
	ur_tdl = vr - cr/c2r;
	dutdl = ul_tdl - ur_tdl;
	if (pl >= pr)
	{
	    z = (c2l*cr/(c2r*cl))*pow(pl/pr,c3l);
	    *pstar = (dutdl > 0.0) ? pl*pow(c2l*dutdl/((1.0 + z)*cl),1.0/c3l) :
				     0.5*min(pl,pr);
	}
	else
	{
	    z = (c2r*cl/(c2l*cr))*pow(pr/pl,c3r);
	    *pstar = (dutdl > 0.0) ? pr*pow(c2r*dutdl/((1.0 + z)*cr),1.0/c3r) :
				     0.5*min(pl,pr);
	}

//#if defined(UNRESTRICTED_THERMODYNAMICS)	//what is this?	FIXME
	*p_min = -HUGE_VAL;
//#else //defined(UNRESTRICTED_THERMODYNAMICS)
//	*pstar = max(*pstar,Min_pressure(Tsl));
//	*pstar = max(*pstar,Min_pressure(Tsr));
//	*p_min = Min_pressure(Tsl);
//#endif //defined(UNRESTRICTED_THERMODYNAMICS)
	*pstar = max(*pstar,*p_min);
 */
}

extern double riemann_wave_curve(
	STATE	*st0,
	double	pstar)
{
	double rho0, p0;
	double c1, c2, c3;

	rho0 = st0->dens;
	p0 = st0->pres;

//#if !defined(UNRESTRICTED_THERMODYNAMICS)
//	if (pstar < Min_pressure(st0))
//	    pstar = Min_pressure(st0);
	if (pstar < 1e-12)	//should be min_pressure	FIXME
	    pstar = 1e-12;
//#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	c1 = Coef1(st0);
	c2 = Coef2(st0);
	c3 = Coef3(st0);

	return (pstar < p0) ?
	        EosSoundSpeed(st0)*(pow(pstar/p0,c3) - 1.0)/ c2 :
	        (pstar-p0)/sqrt(rho0*(c1*pstar+c2*p0));
}

extern double acoustic_impedance_squared(STATE *state)
{
	double gamma, i2;

	//gamma = st->eos->gamma;
	i2 =  Gamma(state) * state->pres * state->dens;

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
	//st1->dens = st0->dens*pow(p1/st0->pres,1.0/st0->eos->gamma);
        st1->dens = st0->dens*pow(p1/st0->pres,1.0/Gamma(st0));
	st1->pres = p1;

        //HK NEED TO CHECK for pdens
        st1->pdens[0] = (st0->pdens[0]/st0->dens)*st1->dens;
        st1->pdens[1] = (st0->pdens[1]/st0->dens)*st1->dens;
}

extern void EosViscTherm(
	STATE	*st,
	double	*visc,
	double	*therm)
{
    	*visc = st->eos->tbl_visc;
	*therm = st->eos->tbl_therm;
}
