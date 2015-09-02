#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <list>
#include <vector>
#include <map>
#include <set>

#include "iFluid.h"
#include "iFDroplet.h"

enum FACE_DIR
{
    XY = 0, YZ = 1, XZ = 2
};

using namespace std;

class Incompress_Solver_Smooth_Basis;

class Graph_Point
{
public:
    int comp;
    int set_index;
    int global_index;
    int local_index;
    std::list<std::list<Graph_Point*> >::iterator pSet; 
    //point to the corresponding connected component
    std::list<Graph_Point*> edge; //adjacency list

    Graph_Point():comp(-1),global_index(-1),local_index(-1),set_index(-1) {}

};



class Graph
{
public:
    virtual ~Graph() {ClearAll();}

protected:
    Graph(void):num_v(0) {}
    virtual void ConstructVertex(void) = 0;
    virtual void ConstructEdge(void) = 0;
    void MakeSet(void);
    void ClearAll(void) 
    {
	ClearEdge(); ClearSet();
	connected_sets.clear(); vertices.clear();
	num_v = 0;
    }
    
    void FindConnectedComp(void);//For each processor
    void NumberOfComponents();//For each processor

    void ClearEdge(void);
    void ClearSet(void);

    void UnionSet1(vector<Graph_Point>::iterator to, list<Graph_Point*>::iterator from);
    void UnionSet2(list<Graph_Point*>::iterator to, vector<Graph_Point>::iterator from);
    

    /**************************************************************************/
protected:
    std::list<std::list<Graph_Point*> > connected_sets;
    std::vector<Graph_Point> vertices;

    int num_v;
};

class CellCornerGraph:public Graph
{


public:
    CellCornerGraph(Incompress_Solver_Smooth_Basis *solver):Graph(),pSolver(solver),
    xy_choose(false), yz_choose(false), xz_choose(false) {}

    ~CellCornerGraph() {ClearAllCorner();}


    void Init(void);
    void getRawData(void);
    void getPhyCompHistogram(char* out_name);
    void PrintOut(void); //For debugging
    boolean doCalculation(void)
    {
	return pSolver->iFparams->movie_option->output_connectivity;
    }

    void ClearAllCorner(void) 
    {
	phy_comp.clear(); bdry_comp.clear();
	ClearAll();
    }


private:

    void ChooseMarchingCube(bool xy, bool yz, bool xz)
    {
	xy_choose = xy;
	yz_choose = yz;
	xz_choose = xz;
    }

    void ConstructVertex(void);
    void ConstructEdge(void);

    void MergeSet(void); //For parallel


    
    void MergeBdryFrom(int dim, int dir, int rec_id);
    void SendBdryTo(int dim, int dir, int dst_id, int *&bdry);
    void MakeBdryForSend(int dim, int dir, int *&bdry);

    void AddFaceEdge(int n1, int n2, int n3, int n4, FACE_DIR dir);
    inline void AddEdge(int m, int n); //local index m,n

    void SetComp(void);

    /* These functions are realted to combine the information to processor 0 */

    void SendCompTo(int dst_id, int dim);
    void MergeCompFrom(int rec_id, int dim);

    void SendBdryCompTo(int dst_id, int dim);
    void MergeBdryCompFrom(int rec_id, int dim);

    void MergeToSingleProc(void);
    void getPhyComp(void);
    void getPhyComp2(void);
    void getBdryComp(void);
    void adjustPhyCompForMerge(void);
    void outputHistogram(char* out_name);

    void outputCoordList(char* out_name);



private:
    Incompress_Solver_Smooth_Basis *pSolver;
    bool xy_choose, xz_choose, yz_choose;

    std::map<int, pair<int, int> > phy_comp; 
    //set_index as key, pair<number of points, fluid comp> as value
    std::map<int, pair<int, int> > bdry_comp;
    //shared boundary points, 
    //set_index as key, pair<number of points, fluid comp> as value
};

class TriGraph:public Graph
{
public:
    TriGraph(DropletAnalysis *drop):Graph(),pDrop(drop) {}

    ~TriGraph() {ClearAll();}

    Graph_Point* getVertexPointer(void)
    {
	if (vertices.empty())
	    return 0;
	else
	    return (&(vertices[0]));
    }


    void Init(void)
    {
	ConstructVertex();
	ConstructEdge();
    }
    void GetTriComp(void)
    {
	MakeSet();
	FindConnectedComp();
	//NumberOfComponents();
    }

private:

    void ConstructVertex(void) { conTriVertices(pDrop->getNumTris(), pDrop->getGlobalIndex()); }
    void ConstructEdge(void) { conTriEdges(pDrop->getTriNb()); }

    void conTriVertices(int n_tris, int *g_index);
    void conTriEdges(int (*t_nb)[3]);


private:
    DropletAnalysis *pDrop;
};



#endif
