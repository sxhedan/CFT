/**********************************************************************
 * 		iFluid.h
 **********************************************************************/

#ifndef _FT_IFLUID_H_
#define _FT_IFLUID_H_

#include <vector>
#include <petscksp.h>
#include <assert.h>

#include "FronTier.h"
#include "solver.h"

#define         SOLID_COMP		1
#define         LIQUID_COMP1		2
#define         LIQUID_COMP2		3

#define		ifluid_comp(comp)   (((comp) == LIQUID_COMP1 || 	\
		comp == LIQUID_COMP2) ? YES : NO)

enum EBM_COORD
{
    COORD_X = 0,  COORD_Y = 1,  COORD_Z = 2
};

struct _IF_FIELD {
	double *pres;			/* Pressure */
	double **vel;			/* Velocities */
	double *vort;			/* Vorticity in 2D */
	double **vort3d;		/* Vorticity in 3D */
	double *div;			/* Divergence */
};
typedef struct _IF_FIELD IF_FIELD;

struct _IF_MOVIE_OPTION {
	boolean plot_comp;
	boolean plot_pres;
	boolean plot_vort;
	boolean plot_velo;
	boolean plot_cross_section[MAXD]; /* 3D 0: yz; 1: zx; 2: xy */
	boolean plot_vel_vector;	  /* Plot velocity vector field */
	boolean plot_phase_1;
	boolean plot_phase_2;
	boolean plot_dens;
	boolean plot_pres_error;
	boolean plot_divU;
	boolean plot_rhot;
	boolean plot_MB_residual;
	boolean output_droplet_analysis;
	boolean output_connectivity;

	boolean output_tridia;
	boolean output_error;
};
typedef struct _IF_MOVIE_OPTION IF_MOVIE_OPTION;

enum _NS_SCHEME {
	ERROR_SCHEME		= -1,
        SIMPLE			=  1,
        BELL_COLELLA,
        KIM_MOIN,
        PEROT_BOTELLA
};
typedef enum _NS_SCHEME NS_SCHEME;

typedef struct {
        int dim;
        POINTER level_func_params;
	IF_MOVIE_OPTION *movie_option;
	NS_SCHEME num_scheme;
        double rho1;
        double rho2;
	double mu1;
	double mu2;
	double U1[MAXD];
	double U2[MAXD];
        double bvel[2][MAXD];
	double gravity[MAXD];
	double surf_tension;
	double smoothing_radius;
	COMPONENT m_comp1;
	COMPONENT m_comp2;
	IF_FIELD *field;
        double init_vel[MAXD];
        boolean use_couette_init_vel;
} IF_PARAMS;

struct _FLOW_THROUGH_PARAMS {
        POINT *oldp;
        COMPONENT comp;
};
typedef struct _FLOW_THROUGH_PARAMS FLOW_THROUGH_PARAMS;

enum _TIME_FUNC_TYPE {
	CONSTANT		=  1,
	PULSE_FUNC,
	SINE_FUNC	
};
typedef enum _TIME_FUNC_TYPE TIME_FUNC_TYPE;

struct _TIME_DEPENDENT_PARAMS {
	TIME_FUNC_TYPE td_type;
	double v_base[MAXD],p_base;
	double v_peak[MAXD],p_peak;
	double v_tail[MAXD],p_tail;
	double v_amp[MAXD],p_amp;
	double omega,phase;
	double T[10];
};
typedef struct _TIME_DEPENDENT_PARAMS TIME_DEPENDENT_PARAMS;

/******************************************************************************
 * 		iFluid.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 *      Incompress_Solver_Basis::setInitialCondition() //set initial states for
 *      different problems
 * 	Incompress_Solver_Basis::solve() //using projection method to solve the
 * 	problem
 *
 * References:
 ******************************************************************************/

class SOLVER;
class Incompress_Solver_Basis;

//enum VISITED_TYPE {UNVISITED, VISITED, PARTIAL_VISITED, FULL_VISITED};

//-------------------------------------------------
//		STATES
// NOTE:
//      L_STATE/L_STATE_RECT_EDGEshould be put into 
// lcartsn.h. However, there are some trouble
// to compile in that way.
//-------------------------------------------------
// states inside a cell

class L_STATE{
public:
	double m_U[MAXD];		// velocity vector
	double m_exactU[MAXD];		// exact velocity solution
	double m_P;			// pressure
	double m_exactP;		// exact pressure solution

	double m_phi;			// (Brown's) phi
	double m_q;			// lagged pressure
	double m_adv[MAXD];		// advection source term for Diffusion solver

	double m_mu;		        // smoothed in IBM
	double m_mu_old;		// viscosity in previous time step
	double m_rho;		        // smoothed in IBM
	double m_rho_old;		// density in previous time step
	double div_U;			// Velocity divergence
	double grad_q[MAXD];		// Gradient of q
	double f_surf[MAXD];		// Surface force (such as tension)

	L_STATE();			// Default constructor function
	void setZero(void);		// Set back to zero
};
// states on edge

//------------------------------------------------------
//		MESH
//------------------------------------------------------
// note that the following VERTEX2D/RECT_EDGE are different 
// from those defined in MESH2D.h
class Graph;
class CellCornerGraph;
class CellCenterGraph;

class L_POINT {
public:
    int pg_index; //global index
    int pl_index; //local index
    int comp;     //component

    double p_coords[MAXD];

    L_POINT():pg_index(-1),pl_index(-1),comp(-1) {};

    void setCoords(double*, int);
};

class L_RECTANGLE {
public:
	int m_index;			// rectangle index
	int comp;			// component
	L_STATE m_state;		// store the states
	double m_coords[MAXD];		//cell center coordinate
	int icoords[MAXD]; 		//rectangle index i,j,k
	L_POINT *corner[8];

	L_RECTANGLE();

	void setCoords(double*,int);	//set coordinate of cell center
};

class Incompress_Solver_Basis{
public:
       	Incompress_Solver_Basis() {}; 			// constructor
	virtual ~Incompress_Solver_Basis() {};

	virtual void setInitialCondition(void) = 0;	//Initialization
	virtual void solve(double dt) = 0; 		//main step function
};

class Incompress_Solver_Smooth_Basis:public Incompress_Solver_Basis{
public:
    	friend class CellCornerGraph;
    	friend class DropletAnalysis;
	Incompress_Solver_Smooth_Basis(Front &front);   //constructor
	virtual ~Incompress_Solver_Smooth_Basis() {};

	double m_dt;
	double accum_dt;
	double max_speed;
	double max_value; 				//for debugging
	double max_dt;
	double min_dt;
	double *top_h;
	int dim;

	int has_small_edge;

	void PrintMesh(void); //TMP

	void initMesh(void);
	virtual void setAdvectionDt(void); 
	//using max speed and hmin to determine max_dt, min_dt
	void readFrontInteriorStates(char *state_name, bool binary, bool RegridRestart);
	void printFrontInteriorStates(char *state_name, bool binary);
        void printFrontInteriorStatesRegridRep(char *state_name, bool binary, bool RegridRun);
	void initMovieVariables(void);
	void augmentMovieVariables(void);
	void getVelocity(double *p, double *U);
	void initSampleVelocity(char *in_name);

	//Initialization of States
	double (*getInitialState) (COMPONENT,double*,L_STATE&,int,IF_PARAMS*);
	virtual boolean FT_StateStructAtGridCrossing_tmp(Front*, int*, GRID_DIRECTION, COMPONENT, Locstate*, HYPER_SURF**, double*, double t=0);


	//User interface
	virtual void printInteriorVelocity(char *outname, bool binary) = 0;
	virtual void printInteriorVelocity_vtu(char *outname, bool binary) = 0;
	virtual void outputParallelVisitFile(char *outname, bool binary) = 0;
	virtual void printExpandedMesh(char *outname, bool binary) = 0;
	virtual void setInitialCondition(void) = 0;
	virtual void solve(double dt) = 0; // main step function
	virtual void multiDiagnosisStep(char *outname) = 0;
	//Diagnositc code for each step
	virtual void multiDiagnosisInterval(char *outname) = 0;
	//Diagnositc code for movie interval

protected:
	Front *front;
	Front *tempfront;
	// On topological grid
	RECT_GRID *top_grid;
	RECT_GRID *comp_grid; //Computational grid
	PP_GRID *pp_grid;
	double *array;
	double *source;
	double *diff_coeff;
	COMPONENT *top_comp;
	IF_PARAMS *iFparams;
	IF_FIELD  *field;

	int *top_gmax;
	int *lbuf, *ubuf;
	double *top_L, *top_U;
	int **ij_to_I, **I_to_ij;
	int ***ijk_to_I, ***I_to_ijk;

	// Sweeping limites
	int imin, jmin, kmin;
	int imax, jmax, kmax;

	//member data: mesh storage
	std::vector<L_RECTANGLE>   cell_center;
	std::vector<L_POINT> cell_corner;

	//member data:
	int    m_comp[2];
	double m_mu[2];
	double m_rho[2];// two components at most
	double m_sigma; //surface tension
	double m_smoothing_radius;// used by smoothing function

	double hmin; //smallest spacing
	double mu_min; //smallest viscocity
	double rho_min;// smallest density
	double m_t;
	double m_t_old, m_t_int, m_t_new;

	// for parallel partition
	int NLblocks, ilower, iupper;
	int Ntris, Ntris_all, Npoints;
        int tri_ilower, tri_iupper;
	int *n_dist;
	int *n_dist_tri;

	double *x_pp_bdry, *y_pp_bdry, *z_pp_bdry;
	int x_pp_index, y_pp_index, z_pp_index;

protected:
	void setCornerMesh(void);

	void setCornerComp(void);//Set cell corner components
	void setComponent(void); //Set cell center components
	void setDomain();

	// parallelization related functions
	void scatMeshArray(void);
	void setGlobalIndex(void);
	void setGlobalTriIndex(void);
	void commGlobalTriIndex(void);
	void commGlobalTriIndex2(void);
	void setTempFrontAndRedist(void);
	void checkTempFront(void);
	void checkUnsetGlobalIndex(void);

	void makeTrisForSend(int dim, int dir, int *&num_tris, double *&tris0, double *&tris1, double *&tris2, int *&tris_index);
	void makeTrisForSend2(int dim, int dir, int *&num_tris, double *&tris0, double *&tris1, double *&tris2, int *&tris_index);

	void sendTrisTo(int dim, int dir, int dst_id, int *&num_tris, double *&tris0, double *&tris1, double *&tris2, int *&tris_index);

	void recvTrisFrom(int dim, int dir, int rec_id);
	void recvTrisFrom2(int dim, int dir, int rec_id);

	void setIndexMap(void);

/*  These functions should be rewritten in 2D basis and 3D basis classes */
	virtual double getSmoothingFunction(double r) = 0; //Heaviside function
	virtual double getSmoothingFunctionD(double*, double*) = 0; 
		//Heaviside function
	virtual double smoothedDeltaFunction(double*, double*) = 0;
	virtual double smoothedStepFunction(double*, double*, int) = 0;
	virtual void   sampleVelocity() = 0;
	virtual void   setSmoothedProperties(void) = 0;
	virtual void   setSmoProOnePhase(void) = 0;
	virtual void   setSmoProExact(void) = 0;
		//smooth discontinuous properties
	
/****************  Functions related to solve() *********/

	virtual void copyMeshStates(void) = 0;
	virtual void computeAdvection(void) = 0; 
	//compute advection step in first order scheme

	virtual void compAdvectionTerm_coupled(void) = 0; 
	//get 2nd order advection term for coupled system
	virtual void compAdvectionTerm_coupled_upgraded(void) = 0; 
	//Upgraded algorithm for getting 2nd order advection term for coupled system
	
	//The following 2 functions are for decoupled system
	virtual void compAdvectionTerm_decoupled(void) = 0; 
	//get 2nd order advection term
	virtual void compAdvectionTerm_decoupled_upgraded(void) = 0; 
	//Upgraded algorithm for getting 2nd order advection term


	virtual void compDiffWithSmoothProperty_1st_coupled(void) = 0; //1st order coupled diffusion solver
	virtual void compDiffWithSmoothProperty_1st_decoupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_coupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_decoupled(void) = 0;//2nd order decoupled diffusion solver

	virtual void computeProjection(void) = 0;
	virtual void computePressure(void) = 0;
	virtual void computePressurePmI(void) = 0;
	virtual void computePressurePmII(void) = 0;
	virtual void computePressurePmIII(void) = 0;
	virtual void computeGradientQ(void) = 0;
	virtual void computeNewVelocity(void) = 0;
	virtual void computeSourceTerm(double *coords, L_STATE &state) = 0;
	virtual void computeSourceTerm(double *coords, double t, L_STATE &state) = 0;
	virtual void computeSourceTerm_Adv(double *coords, L_STATE &state) = 0;


	virtual void surfaceTension_Fedkiw(double, HYPER_SURF_ELEMENT*, 
		HYPER_SURF*, double*, double, double, double*) = 0;
	virtual void surfaceTension_Peskin(void) = 0;

	virtual void computeSubgridModel(void) = 0;    // subgrid model 

	
/***********************  Utility functions  *******************/

	void   computeExactSolution(double *coords, L_STATE &state);
	void   getRectangleIndex(int indexRectangle, int &i, int &j);
	void   getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	int    getRectangleComponent(int index);	// the center component
	void   getRectangleCenter(int index, double *coords);
	double getDistance(double *coords0, double *coords1);
	int    getComponent(int *icoords);	
	int    getComponent(double *coords);	
	void   save(char *filename);

	bool   interiorTri(TRI *tri);

/************* TMP Functions which are not implemented or used ***********/

	void getNearestInterfacePoint(COMPONENT,double*,double*,double*,
					double*); 
};

///////////////Interface for Embedded Boundary Method////////////////////

class Incompress_Solver_EBM:public Incompress_Solver_Basis{
public:
        Incompress_Solver_EBM(Front &front) {};//constructor
	~Incompress_Solver_EBM() {};

	virtual void setInitialCondition() {};
	virtual void solve(double dt) {};
};




/////////////////////////////////////////////////////////////////////////////////

class Incompress_Solver_Smooth_2D_Basis:
public 	Incompress_Solver_Smooth_Basis{


public:
        Incompress_Solver_Smooth_2D_Basis(Front &front):
	Incompress_Solver_Smooth_Basis(front) {};
	virtual ~Incompress_Solver_Smooth_2D_Basis() {};

	virtual void printExpandedMesh(char* outname, bool binary) = 0;
	virtual void printInteriorVelocity(char *outname, bool binary) = 0;
	virtual void printInteriorVelocity_vtu(char *outname, bool binary) {}
	void outputParallelVisitFile(char *outname, bool binary) {};
	virtual void setInitialCondition(void) = 0;
	virtual void solve(double dt) = 0;
	void multiDiagnosisStep(char* outname) {};
	void multiDiagnosisInterval(char* outname) {};
protected:
	double getSmoothingFunction(double r);
	double getSmoothingFunctionD(double*, double*);
	double smoothedDeltaFunction(double*, double*);
	double smoothedStepFunction(double*, double*, int);
	void sampleVelocity();
	void setSmoothedProperties(void);
	void setSmoProOnePhase(void);
	void setSmoProExact(void) {};

	virtual void copyMeshStates(void) = 0;
	virtual void computeAdvection(void) = 0;
	virtual void compAdvectionTerm_coupled(void) = 0;
	virtual void compAdvectionTerm_coupled_upgraded(void) = 0;
	virtual void compAdvectionTerm_decoupled(void) = 0;
	virtual void compAdvectionTerm_decoupled_upgraded(void) = 0;

	virtual void compDiffWithSmoothProperty_1st_coupled(void) = 0; 
	//1st order coupled diffusion solver
	virtual void compDiffWithSmoothProperty_1st_decoupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_coupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_decoupled(void) = 0;
	//2nd order decoupled diffusion solver

	virtual void computeProjection(void) = 0;
	virtual void computePressure(void) = 0;
	virtual void computePressurePmI(void) = 0;
	virtual void computePressurePmII(void) = 0;
	virtual void computePressurePmIII(void) = 0;
	virtual void computeGradientQ(void) = 0;
	virtual void computeNewVelocity(void) = 0;
	virtual void computeSourceTerm(double *coords, L_STATE &state) = 0;
	virtual void computeSourceTerm(double *coords, double t, L_STATE &state) = 0;
	virtual void computeSourceTerm_Adv(double *coords, L_STATE &state) = 0;

	virtual void surfaceTension_Fedkiw(double, HYPER_SURF_ELEMENT*,
		HYPER_SURF*, double*, double, double, double*) = 0;
	virtual void surfaceTension_Peskin(void) = 0;

	virtual void computeSubgridModel(void) = 0;    // subgrid model 
};


class Incompress_Solver_Smooth_3D_Basis:
public 	Incompress_Solver_Smooth_Basis{

public:
        Incompress_Solver_Smooth_3D_Basis(Front &front):
	Incompress_Solver_Smooth_Basis(front) {};
	virtual ~Incompress_Solver_Smooth_3D_Basis() {};

	virtual void printExpandedMesh(char* outname, bool binary) = 0;
	virtual void printInteriorVelocity(char *outname, bool binary) = 0;
	virtual void printInteriorVelocity_vtu(char *outname, bool binary) = 0;
	virtual void outputParallelVisitFile(char *outname, bool binary) = 0;
	virtual void setInitialCondition(void) = 0;
	virtual void solve(double dt) = 0;
	virtual void multiDiagnosisStep(char* outname) = 0;
	virtual void multiDiagnosisInterval(char* outname) = 0;

protected:
	double getSmoothingFunction(double r);
	double getSmoothingFunctionD(double*, double*);
	double smoothedDeltaFunction(double*, double*);
	double smoothedStepFunction(double*, double*, int);
	void sampleVelocity();
	void setSmoothedProperties(void);
	void setSmoProOnePhase(void);
	void setSmoProExact(void);

	virtual void copyMeshStates(void) = 0;
	virtual void computeAdvection(void) = 0;
	virtual void compAdvectionTerm_coupled(void) = 0;
	virtual void compAdvectionTerm_coupled_upgraded(void) = 0;
	virtual void compAdvectionTerm_decoupled(void) = 0;
	virtual void compAdvectionTerm_decoupled_upgraded(void) = 0;


	virtual void compDiffWithSmoothProperty_1st_coupled(void) = 0; 
	//1st order coupled diffusion solver
	virtual void compDiffWithSmoothProperty_1st_decoupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_coupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_decoupled(void) = 0;
	//2nd order decoupled diffusion solver

	virtual void computeProjection(void) = 0;
	virtual void computePressure(void) = 0;
	virtual void computePressurePmI(void) = 0;
	virtual void computePressurePmII(void) = 0;
	virtual void computePressurePmIII(void) = 0;
	virtual void computeGradientQ(void) = 0;
	virtual void computeNewVelocity(void) = 0;
	virtual void computeSourceTerm(double *coords, L_STATE &state) = 0;
	virtual void computeSourceTerm(double *coords, double t, L_STATE &state) = 0;
	virtual void computeSourceTerm_Adv(double *coords, L_STATE &state) = 0;


	virtual void surfaceTension_Fedkiw(double, HYPER_SURF_ELEMENT*, 
		HYPER_SURF*, double*, double, double, double*) = 0;
	virtual void surfaceTension_Peskin(void) = 0;

	virtual void computeSubgridModel(void) = 0;    // subgrid model 
};

class Incompress_Solver_Smooth_2D_Cartesian:
public 	Incompress_Solver_Smooth_2D_Basis{

public:
        Incompress_Solver_Smooth_2D_Cartesian(Front &front):
	Incompress_Solver_Smooth_2D_Basis(front) {}
	~Incompress_Solver_Smooth_2D_Cartesian() {}

	void printExpandedMesh(char* outname, bool binary) {}
	void printInteriorVelocity(char* outname, bool binary) {}
	void setInitialCondition(void);
	void solve(double dt);
protected:
	void copyMeshStates(void);
	void computeAdvection(void); //first order advection

	void compAdvectionTerm_coupled(void); //not implemented yet
	void compAdvectionTerm_coupled_upgraded(void); //not implemented yet
	void compAdvectionTerm_decoupled(void);
	void compAdvectionTerm_decoupled_upgraded(void);


	void compDiffWithSmoothProperty_1st_coupled(void);
	void compDiffWithSmoothProperty_1st_decoupled(void);
	void compDiffWithSmoothProperty_2nd_coupled(void); //Not implemented yet
	void compDiffWithSmoothProperty_2nd_decoupled(void);

	void computeProjection(void);
	void computeProjection_new(void);
	void computePressure(void);
	void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computeGradientQ(void);
	void computeNewVelocity(void);

	void surfaceTension_Fedkiw(double, HYPER_SURF_ELEMENT*, HYPER_SURF*, double*, double, double, double*);
	//Need to be fixed
	void surfaceTension_Peskin(void) {};
	//Not implemented

	void computeSubgridModel(void);    // subgrid model 

	/***************   Low level computation functions  *************/
	virtual void computeSourceTerm(double *coords, L_STATE &state);
	virtual void computeSourceTerm(double *coords, double t, L_STATE &state);
	virtual void computeSourceTerm_Adv(double *coords, L_STATE &state);
	double computeFieldPointDiv(int*, double**);
	void   computeFieldPointGrad(int*, double*, double*);
	void   computeFieldPointGradPhi(int*, double*, double*);
	double computeFieldPointCurl(int*, double**, double*);
	double getVorticity(int i, int j);
	//---------------------------------------------------------------
	//         utility functions for the advection step
	//---------------------------------------------------------------
	
	void getAdvectionTerm_decoupled(int *icoords, double convectionTerm[2]);         
	//get second-order advection term at each cell
	void getAdvectionTerm_decoupled_upgraded(int *icoords, double convectionTerm[2]); 
	//Upgraded version of getAdvectionTerm
	void getAdvectionTerm_coupled(int *icoords, double convectionTerm[2]);
	//get second-order advection term at each cell with cross derivative terms
	void getAdvectionTerm_coupled_upgraded(int *icoords, double convectionTerm[2]);
	//Upgraded version of getAdvectionTerm_coupled


	void getFaceVelocity_middleStep(int *icoords,GRID_DIRECTION dir, L_STATE &state_face);
	void getFaceVelocity_middleStep_hat(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat); 
	//For the upgraded algorithm
	void getFaceVelocity_middleStep_bar(int *icoords,GRID_DIRECTION dir,L_STATE &state_bar, double transverseD[2], L_STATE state_hat); 
	//For the upgraded algorithm

	void getFaceVelocity_middleStep_coupled(int *icoords,GRID_DIRECTION dir, L_STATE &state_face); //for variable mu
	void getFaceVelocity_middleStep_coupled_hat(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat) {}; 
	//For the upgraded algorithm with variable mu
	void getFaceVelocity_middleStep_coupled_bar(int *icoords,GRID_DIRECTION dir,L_STATE &state_bar, double transverseD[2], L_STATE state_hat); 
	//For the upgraded algorithm with variable mu

	void getDifffusion(int *icoords,double diffusion[2]); //Get the diffusion terms
	void getDiffusion_coupled(int *icoords, double diffusion[2]); //Get the diffusion terms for variable mu

	void getDU2(int *icoords,EBM_COORD xyz,double dU2[2]); //Get the U_xx, U_yy, U_zz
	void getLimitedSlope(int *icoords,EBM_COORD xzy,double slope[2]); //mimmod slope limiter
	void getLimitedSlope_Vanleer(int *icoords,EBM_COORD xyz, double slope[2]); //Van Leer slope limiter

	double EBM_minmod(double x, double y); //minmod function

	bool getNeighborOrBoundaryState(int icoords[2],GRID_DIRECTION dir,L_STATE &state,double t); //get the neighbor state or boundary state
	void getRiemannSolution(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans); //Compute Riemann solution using left and right state


	//void   computeVelDivergence(void);
	//void getVelocityGradient(double *p,double *gradU,double *gradV);

};

class Incompress_Solver_Smooth_3D_Cartesian:
public 	Incompress_Solver_Smooth_3D_Basis{

public:
        Incompress_Solver_Smooth_3D_Cartesian(Front &front):
	Incompress_Solver_Smooth_3D_Basis(front) {};
	~Incompress_Solver_Smooth_3D_Cartesian() {};

	void printExpandedMesh(char* outname,bool binary);
	void printInteriorVelocity(char* outname,bool binary);
	void printInteriorVelocity_vtu(char* outname,bool binary) {}
	void outputParallelVisitFile(char *outname, bool binary) {}
	void setInitialCondition(void);
	void solve(double dt);
	void multiDiagnosisStep(char* outname) {};
	void multiDiagnosisInterval(char* outname) {};
protected:
	void copyMeshStates(void);
	void computeAdvection(void);
	void compAdvectionTerm_coupled(void); 
	void compAdvectionTerm_coupled_upgraded(void); 
	void compAdvectionTerm_decoupled(void);
	void compAdvectionTerm_decoupled_upgraded(void);

	void compDiffWithSmoothProperty_1st_coupled(void);
	void compDiffWithSmoothProperty_1st_decoupled(void);
	void compDiffWithSmoothProperty_2nd_coupled(void); 
	void compDiffWithSmoothProperty_2nd_decoupled(void);


	void computeProjection(void);
	void computePressure(void);
	void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computeGradientQ(void);
	void computeNewVelocity(void);

	void surfaceTension_Fedkiw(double, HYPER_SURF_ELEMENT*, HYPER_SURF*, double*, double, double, double*);

	//Need to be fixed
	void surfaceTension_Peskin(void) {};
	//Not implemented
	void computeSubgridModel(void) {};    // subgrid model 
	//Not implemented

	virtual void computeSourceTerm(double *coords, L_STATE &state);
	virtual void computeSourceTerm(double *coords, double t, L_STATE &state);
	virtual void computeSourceTerm_Adv(double *coords, L_STATE &state);
	double computeFieldPointDiv(int*, double**);
	void   computeFieldPointGrad(int*, double*, double*);
	void   computeFieldPointGradPhi(int*, double*, double*);
	double computeFieldPointCurl(int*, double**, double*);
	double getVorticityX(int i, int j, int k);
	double getVorticityY(int i, int j, int k);
	double getVorticityZ(int i, int j, int k);

	//void   computeVelDivergence(void);
	//void getVelocityGradient(double *p,double *gradU,double *gradV);
	//-------------------------------------------------
	//	utility function for the advection step
	//-------------------------------------------------
	void getAdvectionTerm_decoupled(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_decoupled_upgraded(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_coupled(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_coupled_upgraded(int *icoords, double convectionTerm[3]);
	//not implemented yet 
	void getFaceVelocity_middleStep(int *icoords,GRID_DIRECTION dir,L_STATE &state_face);
	void getFaceVelocity_middleStep_hat(int *icoords,GRID_DIRECTION dir, L_STATE &state_hat);
	void getFaceVelocity_middleStep_bar(int *icoords,GRID_DIRECTION dir, L_STATE &state_bar, double transverseD[3], L_STATE state_hat);
 
	void getFaceVelocity_middleStep_coupled(int *icoords,GRID_DIRECTION dir,L_STATE &state_face);
	void getFaceVelocity_middleStep_coupled_hat(int *icoords,GRID_DIRECTION dir, L_STATE &state_hat) {};
	void getFaceVelocity_middleStep_coupled_bar(int *icoords,GRID_DIRECTION dir, L_STATE &state_bar, double transverseD[3], L_STATE state_hat);


	void getDifffusion(int *icoords,double diffusion[3]);
	void getDiffusion_coupled(int *icoords, double diffusion[3]);

	void getDU2(int *icoords,EBM_COORD xyz,double dU2[3]);
	void getLimitedSlope(int *icoords,EBM_COORD xzy,double slope[3]);
	void getLimitedSlope_Vanleer(int *icoords, EBM_COORD xyz, double slope[3]);

	double EBM_minmod(double x, double y);

	bool getNeighborOrBoundaryState(int icoords[3],GRID_DIRECTION dir,L_STATE &state,double t);
	void getRiemannSolution(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans);
};




////////////////////////////////////////////////////////////////////////////////////////////////////////////////





class Incompress_Solver_Smooth_3D_Cylindrical:public Incompress_Solver_Smooth_3D_Basis{
public:
        Incompress_Solver_Smooth_3D_Cylindrical(Front &front):Incompress_Solver_Smooth_3D_Basis(front) {};
	~Incompress_Solver_Smooth_3D_Cylindrical() {};

	void printExpandedMesh(char *out_name,bool binary);
	void printInteriorVelocity(char *out_name,bool binary);
	void printInteriorVelocity_vtu(char *out_name,bool binary);
	void outputParallelVisitFile(char *outname, bool binary);

	void setInitialCondition(void);
	void solve(double dt);
	void multiDiagnosisStep(char *outname);
	void multiDiagnosisInterval(char *outname);
	void setAdvectionDt(void);

protected:
	void copyMeshStates(void);
	void computeAdvection(void); //First order advection, operator splitting
	void computeAdvection_test(void); //Paper version, first order advection, operator splitting

	void compAdvectionTerm_coupled(void); //not implemented
	void compAdvectionTerm_coupled_upgraded(void) {}; //not implemented

	void compAdvectionTerm_decoupled(void);
	void compAdvectionTerm_decoupled_upgraded(void) {}; //not implemented

	void compDiffWithSmoothProperty_1st_coupled(void) {}; //not implemented
	void compDiffWithSmoothProperty_1st_decoupled(void);
	void compDiffWithSmoothProperty_1st_decoupled_test(void); //Paper version
	void compDiffWithSmoothProperty_1st_decoupled_source(void); //half Crank-Nicholson

	void compDiffWithSmoothProperty_2nd_coupled(void); 
	void compDiffWithSmoothProperty_2nd_decoupled(void); 
	//just adjusting coefficients, do not support complext B.C.

	//---------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------
	//   utility functions for diffusion and projection 
	//----------------------------------------------------------------------------

	void compDiff_CellFace( //Set the coefficient for cell corners
		PETSc *pSolver,
		int I,
		int I_nb[18],
		double U_center[3],
		double U_nb[3][18],
		int flag[6],
		int equation_index,
		int vel_comp,
		int face_index,
		double coeff);


	void compDiff_CellCorner( //Set the coefficient for cell corners
		PETSc *pSolver,
		int I,
		int I_nb[18],
		double U_center[3],
		double U_nb[3][18],
		int flag[6],
		int equation_index,
		int vel_comp,
		int corner_index,
		double coeff);

	double getFaceArea(int *icoords, GRID_DIRECTION dir);
	double getCellVolume(int *icoords);
	void getFaceCenter(int *icoords, GRID_DIRECTION dir, double faceCenter[3]);

	void compDiffWithSmoothProperty_2nd_decoupled_Shuqiang(void);
	//used in accuracy test for cylindrical coordinate
	void compDiffWithSmoothProperty_cellFace(
		PETSc *pSolver,
		int *icoords,
		int I, int I_nb,
		GRID_DIRECTION dir,
		double dh,
		double faceArea,
		double cellVolume,
		double r,
		double mu,
		double rho,
		L_STATE &U_nb,
		L_STATE &U_nb_new,
		L_STATE &U_center);
	
	void compDiffWithSmoothProperty_Dirichlet(
		PETSc *pSolver,
		int *icoords,
		int I, int I_nb,
		GRID_DIRECTION dir,
		double dh,
		double faceArea,
		double cellVolume,
		double r,
		double mu,
		double rho,
		L_STATE &U_nb,
		L_STATE &U_nb_new,
		L_STATE &U_center);

	void compDiffWithSmoothProperty_cellInterior(
		PETSc *pSolver,
		int *icoords,
		int I,
		double cellVolume,
		double r,
		double mu,
		double rho,
		L_STATE &U_center);


	void computeProjection_Shuqiang(void);
	//used in accuracy test for cylindrical coordinate

	//--------------------------------------------------------

	void computeProjection(void);
	void computePressure(void);
	void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computeGradientQ(void);
	void computeNewVelocity(void);

	void surfaceTension_Fedkiw(double, HYPER_SURF_ELEMENT*, HYPER_SURF*, double*, double, double, double*);
	void surfaceTension_Peskin(void);
	void compSurfaceTension_Tri(TRI *triangle);
	double getDiscreteDelta(const double*, const double*);
	double getDiscreteDelta_Phi(double);


	void computeSubgridModel(void);    // subgrid model 

	virtual void computeSourceTerm(double *coords, L_STATE &state);
	virtual void computeSourceTerm(double *coords, double t, L_STATE &state);
	virtual void computeSourceTerm_Adv(double *coords, L_STATE &state);
	double computeFieldPointDiv(int*, double**);
	void   computeFieldPointGrad(int*, double*, double*);
	void   computeFieldPointGradPhi(int*, double*, double*);
	double computeFieldPointCurl(int*, double**, double*);
	double getVorticityX(int i, int j, int k);
	double getVorticityY(int i, int j, int k);
	double getVorticityZ(int i, int j, int k);

	//void   computeVelDivergence(void);
	//void getVelocityGradient(double *p,double *gradU,double *gradV);
	//-------------------------------------------------
	//	utility function for the advection step
	//-------------------------------------------------
	void getAdvectionTerm_decoupled(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_decoupled_upgraded(int *icoords, double convectionTerm[3]) {};
	//not implemented yet

	void getAdvectionTerm_coupled(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_coupled_upgraded(int *icoords, double convectionTerm[3]) {};
	//not implemented yet 
	
	void getFaceVelocity_middleStep(int *icoords,GRID_DIRECTION dir,L_STATE &state_face);
	void getFaceVelocity_middleStep_hat(int *icoords,GRID_DIRECTION dir, L_STATE &state_hat) {}; 
	//not implemented
	void getFaceVelocity_middleStep_bar(int *icoords,GRID_DIRECTION dir, L_STATE &state_bar, double transverseD[3], L_STATE state_hat) {}; 
	//not implemented
 
	void getFaceVelocity_middleStep_coupled(int *icoords,GRID_DIRECTION dir,L_STATE &state_face);
	void getFaceVelocity_middleStep_coupled_hat(int *icoords,GRID_DIRECTION dir, L_STATE &state_hat) {}; 
	//not implemented
	void getFaceVelocity_middleStep_coupled_bar(int *icoords,GRID_DIRECTION dir, L_STATE &state_bar, double transverseD[3], L_STATE state_hat) {}; 
	//not implemented


	void getDifffusion(int *icoords,double diffusion[3]);
	void getDiffusion_coupled(int *icoords, double diffusinon[3]);

	void getDU2(int *icoords,EBM_COORD xyz,double dU2[3], double dU[3]);
	void getLimitedSlope(int *icoords,EBM_COORD xzy,double slope[3]);
	void getLimitedSlope_Vanleer(int *icoords, EBM_COORD xyz, double slope[3]);

	double EBM_minmod(double x, double y);

	bool getNeighborOrBoundaryState(int icoords[3],GRID_DIRECTION dir,L_STATE &state,double t);
	void getRiemannSolution(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans);

	//////////// Diagnosis Code ///////////////////////
	
	void computeError(FILE* outfile);
	void recordMeshLength(void);
	void printNormalAndCurvature(void);
	void computeError_part(FILE* outfile);
	void setRhoMuOld(void);
	void getMassConservation(double &mass_conserv_max, double &mass_conserv_ave, double &divu_max, double &divu_ave);

	void getInterfaceArea(double &area);
	double getTriArea(TRI* triangle);
        void getMaxMinDiameterTri(double &min, double &max);

	void getPhaseVolume(double &vol_phone, double &vol_phtwo);

	void printTridia(char* out_name);

	////////// Special Code for contactor project ///////
	//
	void getIcoords(int index, int *icoords, int *gmax);

	/////////  Utility Code for binary output ///////////
	void printExpandedMesh_little_endian(char *out_name);
	void printExpandedMesh_big_endian(char *out_name);
	void printExpandedMesh_ascii(char* out_name);

	void printInteriorVelocity_little_endian(char *out_name);
	void printInteriorVelocity_big_endian(char *out_name);
	void printInteriorVelocity_ascii(char *out_name);


	void printInteriorVelocity_vtu_ascii(char *out_name);
	void printInteriorVelocity_vtu_binary(char *out_name);

public:
	void adjustReflectionBuffer(void);

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////




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
extern void fluid_print_front_states(FILE*,Front*,bool binary);
extern void fluid_read_front_states(FILE*,Front*,bool binary);
extern void read_iF_movie_options(char*,IF_PARAMS*);
extern void read_iF_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern boolean isDirichletPresetBdry(Front*,int*,GRID_DIRECTION,COMPONENT);

#endif
