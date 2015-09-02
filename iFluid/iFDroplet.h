#ifndef _IFDROPLET_H_
#define _IFDROPLET_H_

#include "iFluid.h"
#include <map>

using namespace std;

typedef int L3ARRAY [3];

class DropletAnalysis
{
public:
    DropletAnalysis(Incompress_Solver_Smooth_Basis *solver):pSolver(solver),
    num_tris(0),num_points(0),num_cells(0),num_corners(0),num_tris_interior(0),num_comps(0),
    grid_spacing(0),ps(0),tris(0),tri_nb(0),comp_index(0),global_index(0),
    pos_map(0),neg_map(0),pos_vol(0),neg_vol(0),comp_list(0),
    tri_interior(0),c_centers(0),c_corners(0),face_type(0),face_sign(0),blank_proc(false),first_non_blank_cell(0),
    num_comps_global(0),gmax_cell(0)
    {}
    ~DropletAnalysis() { ClearAll(); }

    boolean doCalculation(void)
    {
	return pSolver->iFparams->movie_option->output_droplet_analysis;
    }
    void InitFromInputFile(void);
    void InitFromInputFile_noMalloc(void); 
 
    void InitFromInputFile_binary(void);
 
    void Init(void);
    void preInit(void);

    void clearTempFront(void)
    {
	if (pSolver->tempfront != NULL)
	    free_front(pSolver->tempfront);
    }

    void ClearAll(void)
    {
	DeleteMappings();
	DeleteCells();
	DeleteTris();
	DeletePoints();
	
    }

    void printAllTris(void);

    void printAllCells(void);

    void outputAll(void)
    {
	outputAllTris();
	outputAllCells();
    }

    void outputAll_binary(void)
    {
	outputAllTris_binary();
	outputAllCells_binary();
    }

    void outputAllTris_binary(void);
    void outputAllCells_binary(void);


    void outputAllTris(void);
    void outputAllCells(void);

    void outputAllPhyComps(void); //new visulaization 

    void outputHistogram(char* out_name);

    void GetCompIndexForTris(void); //Step 1

    void GetLocalVolume(void); //Step 2

    void adjustSurface(void); //prestep in Step2

    void FillVolumeForBlankCells(void); //Step 3

public:
    int getNumTris(void) {return num_tris;}
    int* getGlobalIndex(void) {return global_index;}
    L3ARRAY* getTriNb(void) {return tri_nb;}

private:

    void CreatePoints(void);
    void CreateTris(void);
    void CreateCells(void);
    void CreateMappings(void);

    void DeletePoints(void);
    void DeleteTris(void);
    void DeleteCells(void);
    void DeleteMappings(void);

private:

    void getLocalCompIndexForTris(void); 
    //Use connected component algorithm to fill the tri comp index for each
    //processor, substep 1.1

    void getGlobalCompIndexForTris(void);
    //Use communication for tri index and mappings to fill the tri comp index
    //for each processor, substep 1.2

    void adjustCompIndexForTris(void);
    //Adjust the comp index for all tris in a continuous way.
    //count the number of global connected surface, substep 1.3

    void getLocalCompIndexForCells(void);
    //Use a DFS through all the cells to fill the index for blank cells 
    //and constuct the mappings between index. Substep 3.1
    //
    void getGlobalCompIndexForCells(void);
    //Communicate the face index and mappings through processors and final
    //adjust. Substep 3.2
    //
    void getFinalVolume(void);
    //Use mmapping and index information to update volume and merge to processor
    //0 for output. Substep 3.3

private:

    /******************Step 1 Functions *****************/

    void commCompIndexForTris(void);
    //Communicate component index for triangles across processors
    //directly change the index, needs verifying
    void sendTriCompTo(int dim, int dir, int dst_id, MPI_Request* request);
    void mergeTriCompFrom(int dim, int dir, int rec_id);

    void commCompIndexForTris2(void);
    //Communicate component index for triangles across processors
    //get index mapping from communication, do not change the comp index for tris
    void makeTriCompForSend2(int dim, int dir, int *&num_tris_send, int *&global_index_send, int *&comp_index_send);
    //void sendTriCompTo2(int dim, int dir, int dst_id, int *&num_tris_send, int *&global_index_send, int *&comp_index_send);
    void sendTriCompTo2(int dim, int dir, int dst_id, MPI_Request* request);
    void mergeTriCompFrom2(int dim, int dir, int rec_id);

    /*********************Step 2 Functions ****************/


    /********************Step 3 Functions ****************/

    void initMapping(void); //Step 3 Initialization

    void adjustFaceIndexFromMapping(void);//Used in 3.1 & 3.2
    void adjustVolumeFromMapping(void);//Used in 3.1 & 3.3


    void visitCell(int from_cell, int cur_cell, int from_dir, int **cell_nb_list, bool *visit_flag);
    //DFS visit sub-function, used in 3.1

    void fillBlankProcFaces(void);
    //Use some abitrary element to initialize blank processors, in 3.2
    void sendIndexTo(int dim, int dir, int dst_id, int *&index_for_send, unsigned char *&sign_for_send);
    void mergeIndexFrom(int dim, int dir, int rec_id);
    void makeIndexForSend(int dim, int dir, int *&index_for_send, unsigned char *&sing_for_send);

    void commCellFace(void);
    //commumication cell faces to construct mapping, in 3.2
    void sendCellFaceTo(int dim, int dir, int dst_id, int *&index_for_send, unsigned char *&sign_for_send);
    void mergeCellFaceFrom(int dim, int dir, int rec_id);
    void makeCellFaceForSend(int dim, int dir, int *&index_for_send, unsigned char *&sign_for_send);

    void mergeVolumeInfo(map<int,double> &vol_for_comm);
    //Merge volume information to processor 0, in 3.3
    void sendVolumeTo(int dst_id, int dim, map<int,double> &vol_for_send);
    void mergeVolumeFrom(int rec_id, int dim, map<int,double> &vol_for_recv);

    void getBlankProcVolume(void);//Used in 3.3
    void adjustVolumeForBlankCells(void);//Used in 3.3

    /******************** Common Functions ******************/
    void commMapping(map<int,int> &map_for_comm); 
    //communicate the mapping across processors
    //first combine the mapping to processor 0, then broadcast to each processor
    void sendMappingTo(int dst_id, int dim, map<int,int> &map_for_send);
    void mergeMappingFrom(int rec_id, int dim, map<int,int> &map_for_recv);

    void adjustMapping(map<int,int> &map_for_comm);
    //Adjust the mapping so that each index would map to the smallest index
    //that it coul map default lessd map to. 
    //E.g, 19--->4, 4--->0 would generate 19--->0 and 4--->0

    void commSet(set<int> &set_for_comm); 
    //communicate the set across processors
    //first combine the mapping to processor 0, then broadcast to each processor
    void sendSetTo(int dst_id, int dim, set<int> &set_for_send);
    void mergeSetFrom(int rec_id, int dim, set<int> &set_for_recv);


    /************* Utility Functions ****************/

    bool allBlankCells(void);
    bool blankCell(int cell_index);
    int findNbCellIndex(int current_id, int dir);
    int getCellIndex(int i, int j, int k)
    {
	return (i + j*gmax_cell[0] + k*gmax_cell[0]*gmax_cell[1]);
    }
    void getCellCorners(int index, double *corners);


    void addPairToMap(int big, int small, map<int,int> &map_for_modify);
    void addPairToPosMap(int big, int small);
    void addPairToNegMap(int big, int small);


private:
    Incompress_Solver_Smooth_Basis *pSolver;

    int num_tris; // m tris
    int num_points; //n points
    int num_comps; // K comps
    int num_cells; // Q cells;
    int num_corners; // Q cells;

    int num_tris_interior;// m' interior tris

    bool *tri_interior; //flag for interior triangles
    double (*ps)[3]; //points n*3 matrix
    int (*tris)[3]; //tris, index to ps m*3 matrix
    int (*tri_nb)[3]; //tris neighbours m*3 matrix
    int *comp_index; //tris component index m*1 vector initial value -2
    int *global_index; //tris global index m*1 vector
    int *comp_list; //comp index list K*1 vector
    int *pos_map; //positive mapping K*1 vector initial to itself
    int *neg_map; //negative mapping K*1 vector initial to itself
    double *pos_vol; //positive volume K*1 vector initial zero
    double *neg_vol; //negative volume K*1 vector initial zero

    double (*c_centers)[3]; //cell centers Q cells Q*3 matrix
    int (*face_type)[6]; //faces type for cells Q*6 matrix
    // initial value = blank cell value = -2
    // crossing face value = -1
    unsigned char (*face_sign)[6]; //true/1 is +, false/0 is - initial value to true

    double (*c_corners)[3]; //cell corners which are consistent globally

    double *grid_spacing; //uniform grid spacing 1*3 vector

    map<int, int> map_comm; //mapping for communication

    map<int, int> map_comm_pos; //mapping for communication
    map<int, int> map_comm_neg; //mapping for communication

    map<int, double> vol_comm_pos;
    map<int, double> vol_comm_neg;

    bool blank_proc;
    int first_non_blank_cell;
    int *gmax_cell;

    int num_comps_global;

    Front *tempfront;

};





































#endif
