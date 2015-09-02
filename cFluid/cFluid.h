/**********************************************************************
 * 		cFluid.h
 **********************************************************************/

#include <FronTier.h>
#include <vector>
#include <list>
#include <assert.h>

#define         EXT_COMP		0
#define         SOLID_COMP		0
#define         GAS_COMP1		2
#define         GAS_COMP2		3
#define         MAX_COMP                10

#define		gas_comp(comp)   (((comp) == GAS_COMP1 || 	\
		comp == GAS_COMP2) ? YES : NO)

enum _PROB_TYPE {
        ERROR_TYPE = -1,
        TWO_FLUID_BUBBLE = 1,
        TWO_FLUID_RT,
        FLUID_SOLID_CIRCLE,
        BUBBLE_SURFACE,
        FLUID_RIGID_BODY,
	ROTOR_ONE_FLUID,
	ROTOR_TWO_FLUID,
        TWO_FLUID_RM,
        IMPLOSION,
        MT_FUSION,
        PROJECTILE,
        RIEMANN_PROB,
	FLUID_CRYSTAL,
	ONED_BLAST,
	ONED_SSINE,
	TWO_FLUID_VST_RM
};
typedef enum _PROB_TYPE PROB_TYPE;

enum _EOS_TYPE {
        ERROR_EOSTYPE = -1,
        POLYTROPIC = 1,
        STIFFENED_POLYTROPIC,
        MULTI_COMP_POLYTROPIC
};
typedef enum _EOS_TYPE EOS_TYPE;

struct _EOS_PARAMS {
        EOS_TYPE eos_type;
        int n_comps;
        double  mgamma[2];
        double  M[2];  /* Molecular weights of components */
        double  R;  /* Ideal gas constant PV = RT */
        double  gamma;
        double  pinf;
        double  einf;
        double  et;
	double	tbl_visc;
	double	tbl_therm;
};
typedef struct _EOS_PARAMS EOS_PARAMS;

struct _STATE {
	double dens;			/* density */
        double pdens[2];         /* partial density */
	double engy;			/* energy density */
	double momn[MAXD];		/* momentum deisnty */
	double pres;			/* Pressure */
	double vel[MAXD];		/* Velocities */
	double vort;			/* Vorticity */
	double vort3d[MAXD];            /* Vorticity in 3D */
	EOS_PARAMS      *eos;
	int dim;
	double gamma;	//Dan
};
typedef struct _STATE STATE;

struct _MOVIE_OPTION {
        boolean plot_dens;
        boolean plot_pres;
        boolean plot_vort;
        boolean plot_velo;
        boolean plot_cross_section[MAXD];  /* 3D 0: yz; 1: zx; 2: xy */
	boolean set_bounds;
	double min_dens,max_dens;
	double min_pres,max_pres;
	double min_velo,max_velo;
	double min_vort,max_vort;
};
typedef struct _MOVIE_OPTION MOVIE_OPTION;

enum _NUM_SCHEME {
	TVD_FIRST_ORDER		=	1,
        TVD_SECOND_ORDER,
        TVD_FOURTH_ORDER,
        WENO_FIRST_ORDER,
        WENO_SECOND_ORDER,
        WENO_FOURTH_ORDER
};
typedef enum _NUM_SCHEME NUM_SCHEME;

enum _POINT_PROP_SCHEME {
	FIRST_ORDER		=	1,
        SECOND_ORDER,
        FOURTH_ORDER
};
typedef enum _POINT_PROP_SCHEME POINT_PROP_SCHEME;

enum _SHOCK_PARAMETER {
	BEHIND_PRESSURE 	=	1,
	BEHIND_VELOCITY,
	SHOCK_SPEED,
	SHOCK_MACH_NUMBER
};
typedef enum _SHOCK_PARAMETER SHOCK_PARAMETER;

enum _RIEMANN_SOLVER_WAVE_TYPE {
	UNSET_RIEMANN_SOLVER_WAVE_TYPE = -1,
	RAREFACTION = 1,
	SHOCK,
	COMPOSITE,
	STRONG_DET,
	CJ_DET,
	LAST_GAS_RIEMANN_SOLVER_WAVE_TYPE
};
typedef enum _RIEMANN_SOLVER_WAVE_TYPE RIEMANN_SOLVER_WAVE_TYPE;

struct _Intfc_extrema{
	bool	do_01, do_05;
	double	rfactor;
	double	h_min, h_min_01, h_min_05;
	double	h_max, h_max_01, h_max_05;
};
typedef struct _Intfc_extrema Intfc_extrema;

struct _Big_State{
    	int	count;		/*num of pts summed within layer*/
    	double	d;		/*density*/
	double	v[MAXD];	/*velocity*/
	double	p;		/*pressure*/
	double	e;		/*energy*/
};
typedef struct _Big_State Big_State;

typedef struct {
        int dim;
        PROB_TYPE prob_type;
        POINTER level_func_params;
	MOVIE_OPTION *movie_option;
	Intfc_extrema *iext;
	NUM_SCHEME num_scheme;
        POINT_PROP_SCHEME point_prop_scheme;
	EOS_PARAMS      eos[MAX_COMP];
	boolean tracked;
        boolean multi_comp_non_reactive;
        int n_comps;
	int shock_dir;
	double p0;
        double p1;
        double p2;
	double rho0;
        double rho1;
        double rho2;
	double v0[MAXD];
        double v1[MAXD];
        double v2[MAXD];
	double mu1;
	double mu2;
        double mu; 
        bool parabolic_step; 
        bool viscosity; 
        bool thermal_conduction;
        bool mass_diffusion; 
        double kappa;  
        double D;  
        bool subgrid; 
        bool subgrid_vis; 
        bool subgrid_con; 
        bool subgrid_md;  
	double gamma;
	double gravity[MAXD];
	double Mach_number;
	double shock_position;
	double contact_vel;
	double **vel;
	double **mom;
	double *dens;
        double **pdens;
	double *engy;
	double *pres;
	double *vort;
	double **vort3d;
	//GFM
	double **gnor;
	double **Gdens;
        double ***Gpdens;
	double ***Gvel;
	double **Gpres;

	//initial diffusion layer
	double thickness_idl;

	//TBL
	bool TBL;
	double tbl_visc;
	double tbl_therm;

	// Base front for comparison
	boolean use_base_soln;
        char base_dir_name[200];
        int num_step;
        int *steps;
        F_BASIC_DATA *f_basic;
} EQN_PARAMS;

struct _SCHEME_PARAMS
{
        double lambda;
        double beta;
	double gamma, einf;
};
typedef struct _SCHEME_PARAMS SCHEME_PARAMS;

struct _FLOW_THROUGH_PARAMS {
        POINT *oldp;
        COMPONENT comp;
        EQN_PARAMS *eqn_params;
};
typedef struct _FLOW_THROUGH_PARAMS FLOW_THROUGH_PARAMS;

struct _RG_PARAMS {
        int dim;
        double  total_mass;             /* Total mass */
        double  moment_of_inertial;     /* Moment of inertial about the axis */
        double  center_of_mass[MAXD];   /* Center of mass */
        double  rotation_dir[MAXD];     /* Direction of rotation */
        double  rotation_cen[MAXD];     /* Center of rotation */
        double  cen_of_mass_velo[MAXD]; /* Center of mass velocity */
        double  angular_velo;           /* Angular velocity of rotation */
        MOTION_TYPE motion_type;
};
typedef struct _RG_PARAMS RG_PARAMS;

typedef class G_CARTESIAN G_CARTESIAN_EB;


struct _VAR_BDRY_PARAMS {
	int dim;
        double center[MAXD];        /* Center of disk/sphere */
        double *angles_pistons;     /* Angles to the pistons' centers */
        double half_angular_width; /* Half-angle of the piston's surface */
        double bdry_vel;            /* Boundary velocity */
        double bdry_dens;           /* Boundary density */
        double bdry_pres;           /* Boundary pressure */
        int number_pistons;         /* Number of pistons */
	double jet_duration_time;   /* Time duration for the jet */
};
typedef struct _VAR_BDRY_PARAMS VAR_BDRY_PARAMS;


/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 * 	G_CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class G_CARTESIAN;

//typedef class G_CARTESIAN G_CARTESIAN_EB;
//enum VISITED_TYPE {UNVISITED, VISITED, PARTIAL_VISITED, FULL_VISITED};

//------------------------------------------------------
//		MESH
//------------------------------------------------------
// note that the following VERTEX2D/RECT_EDGE are different 
// from those defined in MESH2D.h

class L_RECTANGLE {
public:
	int m_index;			// rectangle index
	int comp;			 
	double m_coords[MAXD];	
	int icoords[MAXD];

	L_RECTANGLE();

	void setCoords(double*,int);
};

struct _FIELD
{
	double **vel;
	double **momn;
	double *dens;
	double *engy;
	double *pres;
	double *vort;
	double **vort3d;
	double **pdens;
        double *qt;
        double ***tau;
	double *gamma;	//Dan
};

struct _SWEEP
{
        double *dens;            /* density vector */
        double **pdens;            /* partial density vector */
        double **momn;      /* momentum vector */
        double *engy;            /* internal energy vector */
        double *pres;        /* used for EOS */
        double *GAM;
        double *c;
        COMPONENT *comp;
	double *gamma;	//Dan
};

struct _FSWEEP
{
        double *dens_flux;       /* density flux */
        double **pdens_flux; /* partial density flux */
        double **momn_flux; /* momentum flux */
        double *engy_flux;       /* internal energy flux */
	double *gamma;	//Dan
};

typedef struct _FIELD FIELD;

typedef struct _FSWEEP FSWEEP;

typedef struct _SWEEP SWEEP;

class G_CARTESIAN{
	Front *front;
public:
	//G_CARTESIAN();
	G_CARTESIAN(Front &front);
	int dim;
	double m_dt;			// time increment
	double max_dt;			// max_dt from cartesian
	double hmin;			// smallest spacing
	Intfc_extrema	*iext;	
	Big_State	bst_00[2], amb_bst_00[2];
	Big_State	bst_01[2], amb_bst_01[2];
	Big_State	bst_05[2], amb_bst_05[2];

	void setInitialIntfc(LEVEL_FUNC_PACK*,char*);// setup initial geometry
	void setInitialStates(); 	// setup initial state
	void setProbParams(char*); 	// setup initial state
	void initMesh(void);		// setup the cartesian grid
	void readInteriorStates(char*);
	void printFrontInteriorStates(char*);
	void initMovieVariables();
	void getVelocity(double*,double*);
	void initSampleVelocity(char *in_name);
	void compareWithBaseData(char *out_name);
	void freeBaseFront();
	void printOneDimStates(char*);
	void printInteriorVtk(char*);
	void checkIntfc(char*);
	void record_intfc_extrema();
	void print_intfc_extrema(char*);
	double height_at_fraction(double,double,double,int,COMPONENT);
	void accumulate_fractions_in_layer(double,double*,COMPONENT);
	double find_particular_fluid_cell_volume(double*,double*,COMPONENT);
	double find_intersection_by_bisection(Front*,double*,double*,int);
	void accumulate_state_in_layer(double,Big_State*,bool);
	void hyp_sol(double*,STATE*,COMPONENT);
	void add_state_to_totals(STATE*,Big_State*);
	void copy_Big_State(const Big_State*,Big_State*);
	void accumulate_state_totals(const Big_State*,Big_State*);
	void normalize_state_totals(Big_State*);

	// main step function
	void solve(double dt);		

	// constructor
	~G_CARTESIAN();

private:
	// On topological grid
	RECT_GRID *top_grid;
	double *array;		// for scatter states;
	COMPONENT *top_comp;
	EQN_PARAMS *eqn_params;
	FIELD field;
	FIELD *base_field;
	Front *base_front;

	int top_gmax[MAXD];
	int lbuf[MAXD],ubuf[MAXD];
	double top_L[MAXD],top_U[MAXD],top_h[MAXD];
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D
	int nrad;			// Buffer size for a given solver

	// Sweeping limites
	int imin[MAXD];
	int imax[MAXD];

	// member data: mesh storage
	std::vector<L_RECTANGLE> 	cell_center;

	// member data: 
	int m_comp[2];
	double m_mu[2];
	double m_dens[2];		// two component at most
	double m_smoothing_radius;	// used by getSmoothingFunction()

	double m_t;                     // time
	double max_speed;		// for stability of convection
	double min_dens,min_pres;	// minimum physical variables


	// for parallel partition
	int             NLblocks,ilower,iupper;
        int             *n_dist;

	// mesh: full cells mesh
	void setComponent_old(void);	// init components	
        void setComponent(void);        
        void setRPGhost(int,STATE*,int);
	void setDomain();
	void augmentMovieVariables(void);
	void copyMeshStates();
	void sampleVelocity();
	void sampleVelocity2d();
	void sampleVelocity3d();

	/*TMP*/
	void checkVst(SWEEP*);
	void checkFlux(FSWEEP*);

	// parallelization related functions
	//
	void scatMeshArray();
	void scatMeshStates();
	void scatMeshVst(SWEEP*);
	void scatMeshFlux(FSWEEP*);
	void state_reflect(int,double*);
	void reflect_buffer(int,int,int,int*,int*,int*,double*);

	// -------------------------------------------------------
	// 		compressible solver functions
	// -------------------------------------------------------
	void setAdvectionDt(void);
	void computeAdvection(void);
	void computeParab(void);

	void computeTBL();	//Dan	FIXME
	void oned_turbulence_boundary_layer(SWEEP*,int,int);	//Dan
	void TBL2D(SWEEP*,int,int);
	void TBL3D(SWEEP*,int,int);
	void TBLsolver(SWEEP*,int*,int,int);	//Dan
	void computeTBLflux(STATE,double*,double,double,int);	//Dan

	/* Mesh memory management */
	bool withinStencilLen(int*,int);
	void allocMeshVst(SWEEP*);
	void allocMeshFlux(FSWEEP*);
	void allocDirVstFlux(SWEEP*,FSWEEP*);
	void freeVst(SWEEP*);
	void freeFlux(FSWEEP*);

	/* Mesh operations */
	void solveRungeKutta(int);
	void addMeshFluxToVst(SWEEP*,FSWEEP,double);
	void computeMeshFlux(SWEEP,FSWEEP*,double);
	void copyMeshVst(SWEEP,SWEEP*);
	void copyFromMeshVst(SWEEP);
	void copyToMeshVst(SWEEP*);
	void addSourceTerm(SWEEP*,FSWEEP*,double);

	/* Directional flux solver */
	void resetFlux(FSWEEP*);
	void addFluxInDirection(int,SWEEP*,FSWEEP*,double);
	void addFluxInDirection1d(int,SWEEP*,FSWEEP*,double);
	void addFluxInDirection2d(int,SWEEP*,FSWEEP*,double);
	void addFluxInDirection3d(int,SWEEP*,FSWEEP*,double);
	void augmentOneDimBuffer(int,int);
	void numericalFlux(POINTER,SWEEP*,FSWEEP*,int);
	void appendStencilBuffer2d(SWEEP*,SWEEP*,int,int);
	void appendStencilBuffer3d(SWEEP*,SWEEP*,int,int,int);
	void appendGhostBuffer(SWEEP*,SWEEP*,int,int*,int,int);
	// -------------------------------------------------------
	// 		initialization functions
	// -------------------------------------------------------
	void initSinePertIntfc(LEVEL_FUNC_PACK*,char*);
	void initVSTRMIntfc(LEVEL_FUNC_PACK*,char*);
	void initCirclePlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initImplosionIntfc(LEVEL_FUNC_PACK*,char*);
	void initProjectileIntfc(LEVEL_FUNC_PACK*,char*);
	void initMTFusionIntfc(LEVEL_FUNC_PACK*,char*);
	void initRiemannProb(LEVEL_FUNC_PACK*,char*);
	void initRayleiTaylorStates();
	void initRichtmyerMeshkovStates();
	void initVSTRMStates();
	void initBubbleStates();
	void initImplosionStates();
	void initMTFusionStates();
	void initProjectileStates();
	void initRiemProbStates();
	void initBlastWaveStates();
	void initShockSineWaveStates();
	void setRayleiTaylorParams(char*);
	void setRichtmyerMeshkovParams(char*);
	void setVSTRMParams(char*);
	void setBubbleParams(char*);
	void setImplosionParams(char*);
	void setMTFusionParams(char*);
	void setProjectileParams(char*);
	void setRiemProbParams(char*);
	void setOnedParams(char*);
	void readBaseFront(int i);
	void readBaseStates(char *restart_name);
	void readFrontInteriorStates(char *restart_state_name);
        void compSGS2D(SWEEP*);       
        void compSGS3D(SWEEP*);
        void parab_step2D(void);
        void parab2D(SWEEP*);  
        void parab3D(SWEEP*);   
        void appendGhostBufferforParab3D(SWEEP*);       
        void appendGhostBufferforParab2D(SWEEP*);       
        void setNeumannStatesforParab(SWEEP*,HYPER_SURF*,STATE*,int*,int,int,int);      
        void getPressureJumpParameter(double *coords0, double *coords1,
                        double &theta, double &jumpPressure,
                        double &jumpDerivative);

	// velocity field query
	void getVelocityGradient(double *p,double *gradU,double *gradV);

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	void getRectangleIndex(int indexRectangle, int &i, int &j);
	void getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	int getRectangleComponent(int index);	// the center component
	void getRectangleCenter(int index, double *coords);
	void getRectangleCenter(int index0, int index1, double *coords);
	
	int getInteger(double i);
	boolean isInteger(double i);

	double getVorticity(int i, int j);
        double getVorticityX(int i, int j, int k);
        double getVorticityY(int i, int j, int k);
        double getVorticityZ(int i, int j, int k);
	double getDistance(double *coords0, double *coords1);

			// incompletely implemented
	void getNearestInterfacePoint(double *q,double *p); 
        
	int  getComponent(int index);		
	int  getComponent(int *icoords);	
	int  getComponent(double *coords);	
	
	void save(char *filename);

	void setDirichletStates(STATE*,SWEEP*,SWEEP*,HYPER_SURF*,int*,int,
					int,int,int);
	void setNeumannStates(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,int*,int,int,int,
					int,int);

	//GFM
	void solve_exp_value();
	boolean get_ave_normal(int*,int***);
	boolean get_ave_state(SWEEP, int*,int***,int,int);
	boolean needBufferFromIntfc(COMPONENT,COMPONENT);
	void get_normal_from_front();
	void get_ghost_state(SWEEP, int,int);
	void tecplot_interior_states(char*);
	void scatMeshGhost();
	void GFMGhostState_old(int*,int,STATE*);
	void GFMGhostState(int*,int*,int,STATE*,int,SWEEP*,int);
	void find_cross2d(double*,double*,double*,int,double*,int);
	void find_cross3d(double*,double*,double*,int,double*,int);
	boolean point_on_segment(double*,double*,double*);
	boolean cross_line_segment_triangle(double*,double*,double*,double*,double*,double*);
	void checkCorrectForTolerance(STATE*);
};

extern double getStateDens(POINTER);
extern double getStatePdens0(POINTER);
extern double getStatePdens1(POINTER);
extern double getStateEngy(POINTER);
extern double getStateXmom(POINTER);
extern double getStateYmom(POINTER);
extern double getStateZmom(POINTER);
extern double getStatePres(POINTER);
extern double getStateVort(POINTER);
extern double getStateXvel(POINTER);
extern double getStateYvel(POINTER);
extern double getStateZvel(POINTER);
extern double getStateXvort(POINTER);
extern double getStateYvort(POINTER);
extern double getStateZvort(POINTER);
extern double burger_flux(double,double,double);
extern double linear_flux(double,double,double,double);
extern void read_dirichlet_bdry_data(char*,Front*);
extern void restart_set_dirichlet_bdry_function(Front*);
extern void cF_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void cF_variableBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void cFluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
/*	rgbody.c functions */
extern void init_moving_bodies(EQN_PARAMS,char*,Front*);
extern void cfluid_compute_force_and_torque(Front*,HYPER_SURF*,double,double*,
			double*);
extern void record_moving_body_data(char*,Front*);
extern void read_cFluid_params(char*,EQN_PARAMS*);
extern void read_movie_options(char*,EQN_PARAMS*);
extern void read_statistics_options(char*,EQN_PARAMS*);
extern void readFrontStates(Front*,char*);
extern void reflectVectorThroughPlane(double*,double*,double*,int);
extern boolean reflectNeumannState(Front*,HYPER_SURF*,double*,COMPONENT,SWEEP*,
			STATE*);
extern double EosPressure(STATE*);
extern double EosCV(STATE *state);
extern double EosCP(STATE *state);
extern double EosTemperature(STATE *state);
extern double EosInternalEnergy(STATE*);
extern double EosEnergy(STATE*);
extern double EosEnthalpyDifference(STATE*);
extern double EosSoundSpeed(STATE*);
extern double EosSoundSpeedSqr(STATE*);
extern double EosMaxBehindShockPres(double,STATE*);
extern void   EosSetTVDParams(SCHEME_PARAMS*,EOS_PARAMS*);
extern double EosGamma(STATE*);
extern void   CovertVstToState(STATE*,SWEEP*,EOS_PARAMS*,int,int);
extern void findGhostState(STATE,STATE,STATE*);
extern double GAM(STATE*);
extern double Coef1(STATE*);
extern double Coef2(STATE*);
extern double Coef2(STATE*);
extern double Coef3(STATE*);
extern double Coef3(STATE*);
extern void EosViscTherm(STATE*,double*,double*);
extern void initialize_riemann_solver(STATE*,STATE*,double*,double*,double,double*,double*);
extern double riemann_wave_curve(STATE*,double);
extern double acoustic_impedance_squared(STATE*);
extern double acoustic_impedance(STATE*);
extern double mass_flux(double,STATE*);	
extern double dens_Hugoniot(double,STATE*);
extern void state_on_adiabat_with_pr(STATE*,double,STATE*);
extern void state_on_adiabat_with_pr(STATE,double,STATE*);
extern bool find_mid_state(STATE*,STATE*,double,double*,double*,double*,double*,double*,double*,
			   RIEMANN_SOLVER_WAVE_TYPE*,RIEMANN_SOLVER_WAVE_TYPE*);
extern void midstate(STATE*,STATE*,double,double,double,RIEMANN_SOLVER_WAVE_TYPE,int);

	/* Structures and functions for TVD scheme */

extern	void TVD_flux(POINTER,SWEEP*,FSWEEP*,int);
extern	void WENO_flux(POINTER,SWEEP*,FSWEEP*,int);

class EOS{
	EOS_PARAMS *params;
public:
	EOS(EOS_PARAMS &params);
	double Pressure(STATE);
	double Energy(STATE);
};

class COMPRESSIBLE_GAS_SOLVER{
	Front *front;
public: 
	COMPRESSIBLE_GAS_SOLVER(Front &front);

	void set_solver_domain(void);
        void runge_kutta(void);
        void solve(void);
        void solve1d(void);
        void solve2d(void);
        void solve3d(void);
private:
	int dim;
        COMPONENT *top_comp;
        double *top_h;
        int *top_gmax;
        int imin,jmin,kmin;
        int imax,jmax,kmax;
};
