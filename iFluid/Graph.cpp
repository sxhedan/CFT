/*
 * Graph.cpp
 *
 *  Created on: May 4, 2011
 *      Author: Yijie
 */

#include "Graph.h"
#include "iFluid.h"



void Graph::ClearEdge(void)
{
    vector<Graph_Point>::iterator itp;

    for (itp = vertices.begin(); itp != vertices.end(); itp++)
	itp->edge.clear();
}

void Graph::ClearSet(void)
{
    vector<Graph_Point>::iterator itp;
    for (itp = vertices.begin(); itp != vertices.end(); itp++)
    {
	itp->set_index = -1;
	itp->pSet = connected_sets.begin();
    }

    list<list<Graph_Point*> >::iterator itset;

    for (itset = connected_sets.begin(); itset != connected_sets.end(); itset++)
	itset->clear();
    connected_sets.clear();
}



void Graph::MakeSet(void)
{
    ClearSet();
    int i;

    for (i = 0; i < num_v; i++)
    {
        vertices[i].set_index = vertices[i].global_index;
        list<Graph_Point*> self(1, &(vertices[i]));
        connected_sets.push_back(self);
        vertices[i].pSet = --(connected_sets.end());
    }
}

void Graph::UnionSet1(vector<Graph_Point>::iterator to, list<Graph_Point*>::iterator from)
{

    list<list<Graph_Point*> >::iterator iterase = (*from)->pSet;
    list<Graph_Point*>::iterator itp;

    int num_merged = (int) (*from)->pSet->size();

    for (itp = (*from)->pSet->begin(); itp != (*from)->pSet->end(); itp++)
    {
        (*itp)->set_index = to->set_index;
        to->pSet->push_back((*itp));
    }

    for (itp = to->pSet->end(); num_merged > 0; num_merged--)
    {
        itp--;
        (*itp)->pSet = to->pSet;
    }

    connected_sets.erase(iterase);

}

void Graph::UnionSet2(list<Graph_Point*>::iterator to, vector<Graph_Point>::iterator from)
{

    list<list<Graph_Point*> >::iterator iterase = from->pSet;
    list<Graph_Point*>::iterator itp;

    int num_merged = (int) from->pSet->size();

    for (itp = from->pSet->begin(); itp != from->pSet->end(); itp++)
    {
        (*itp)->set_index = (*to)->set_index;
        (*to)->pSet->push_back((*itp));
    }

    for (itp = (*to)->pSet->end(); num_merged > 0; num_merged--)
    {
        itp--;
        (*itp)->pSet = (*to)->pSet;
    }

    connected_sets.erase(iterase);

}

void Graph::FindConnectedComp()
{
    vector<Graph_Point>::iterator itp;
    list<Graph_Point*>::iterator itedge;
    for (itp = vertices.begin(); itp != vertices.end(); itp++)
    {
        for (itedge = (itp->edge).begin(); itedge != (itp->edge).end(); itedge++)
        {
            if (itp->set_index == (*itedge)->set_index)
                continue;

            if (itp->set_index < (*itedge)->set_index)
                UnionSet1(itp, itedge);
            else if (itp->set_index > (*itedge)->set_index)
                UnionSet2(itedge, itp);
        }
    }
}


void Graph::NumberOfComponents()
{
    int size = (int) connected_sets.size();
    printf("The number of components is: %d\n", size);

    list<list<Graph_Point*> >::iterator itsets;

    int i = 1;

    for(itsets = connected_sets.begin(); itsets != connected_sets.end(); itsets++)
    {
	printf("\nThe %d-th set has length %d:"
		,i, (int) itsets->size());
	printf("\nSet index = %d", (*itsets->begin())->set_index);
	printf("\nSet comp = %d\n", (*itsets->begin())->comp);
	i++;
    }

}

/***************************************************************************/

void CellCornerGraph::PrintOut()
{
    int size = (int) vertices.size();
    printf("There are %d points here.\n",size);
    int i = 0;

    vector<Graph_Point>::iterator it;
    list<Graph_Point*>::iterator ite;

    for(it = vertices.begin(); it != vertices.end(); it++)
    {
	printf("Point %d: Global index: %d Comp: %d\n",i, it->global_index, it->comp);
	printf("Edges connected to it:\n");
    
	for(ite = it->edge.begin(); ite != it->edge.end(); ite++)
	    printf("(%d,%d)\n",it->local_index, (*ite)->local_index);
	i++;
    }
    if (connected_sets.empty())
	printf("The current set is empty\n");
}

void CellCornerGraph::ConstructVertex()
{
    vertices.clear();
    num_v = (int) (pSolver->cell_corner.size());
    printf("\nnumber of vertices: %d\n", num_v);
    
    Graph_Point m_point;
    vertices.insert(vertices.end(), num_v, m_point);

    for (int i = 0; i < num_v; i++)
    {
	vertices[i].local_index = pSolver->cell_corner[i].pl_index;
	vertices[i].global_index = pSolver->cell_corner[i].pg_index;
    }
}
void CellCornerGraph::SetComp(void)
{
    int i;

    for (i = 0; i < num_v; i++)
	vertices[i].comp = pSolver->cell_corner[i].comp;
}

void CellCornerGraph::AddFaceEdge(int n1, int n2, int n3, int n4, FACE_DIR dir)
{
    if (vertices[n1].comp == vertices[n2].comp)
	AddEdge(n1, n2);
    if (vertices[n2].comp == vertices[n3].comp)
	AddEdge(n2, n3);
    if (vertices[n3].comp == vertices[n4].comp)
	AddEdge(n3, n4);
    if (vertices[n4].comp == vertices[n1].comp)
	AddEdge(n4, n1);

    if (    (vertices[n1].comp == LIQUID_COMP1) && (vertices[n3].comp == LIQUID_COMP1)
	 && (vertices[n2].comp == LIQUID_COMP2) && (vertices[n4].comp == LIQUID_COMP2)  )
    {
	switch(dir)
	{
	case XY:
	    if (xy_choose) //If choose == true, then comp1 are connected
		AddEdge(n1,n3);
	    else	   //If choose == false, then comp2 are connected
		AddEdge(n2,n4);
	    break;
	case YZ:
	    if (yz_choose)
		AddEdge(n1,n3);
	    else
		AddEdge(n2,n4);
	    break;
	case XZ:
	    if (xz_choose)
		AddEdge(n1,n3);
	    else
		AddEdge(n2,n4);
	    break;
	default:
	    printf("\n Error in choosing face direction! \n");
	    clean_up(ERROR);
	}
    }
    else if (    (vertices[n1].comp == LIQUID_COMP2) && (vertices[n3].comp == LIQUID_COMP2)
	      && (vertices[n2].comp == LIQUID_COMP1) && (vertices[n4].comp == LIQUID_COMP1)  )
    {
	switch(dir)
	{
	case XY:
	    if (xy_choose) //If choose == true, then comp1 are connected
		AddEdge(n2,n4);
	    else	   //If choose == false, then comp2 are connected
		AddEdge(n1,n3);
	    break;
	case YZ:
	    if (yz_choose)
		AddEdge(n2,n4);
	    else
		AddEdge(n1,n3);
	    break;
	case XZ:
	    if (xz_choose)
		AddEdge(n2,n4);
	    else
		AddEdge(n1,n3);
	    break;
	default:
	    printf("\n Error in choosing face direction! \n");
	    clean_up(ERROR);
	}
    }
}

void CellCornerGraph::AddEdge(int m, int n)
{
    vertices[m].edge.push_back(&(vertices[n]));
    vertices[n].edge.push_back(&(vertices[m]));
}

void CellCornerGraph::ConstructEdge()
{
    ClearEdge();
    SetComp();

    int i,j,k,index,l;

    int plindex[8];

    int imin = pSolver->imin;
    int imax = pSolver->imax;
    int jmin = pSolver->jmin;
    int jmax = pSolver->jmax;
    int kmin = pSolver->kmin;
    int kmax = pSolver->kmax;

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,pSolver->top_gmax);
	for (l = 0; l < 8; l++)
	    plindex[l] = pSolver->cell_center[index].corner[l]->pl_index;

	AddFaceEdge(plindex[0],plindex[1],plindex[2],plindex[3], XY);
	AddFaceEdge(plindex[0],plindex[4],plindex[5],plindex[1], XZ);
	AddFaceEdge(plindex[0],plindex[3],plindex[7],plindex[4], YZ);
    }

    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	k = kmax;
	index = d_index3d(i,j,k,pSolver->top_gmax);
	for (l = 0; l < 8; l++)
	    plindex[l] = pSolver->cell_center[index].corner[l]->pl_index;

	AddFaceEdge(plindex[4],plindex[5],plindex[6],plindex[7], XY);
    }

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    {
	i = imax;
	index = d_index3d(i,j,k,pSolver->top_gmax);
	for (l = 0; l < 8; l++)
	    plindex[l] = pSolver->cell_center[index].corner[l]->pl_index;

	AddFaceEdge(plindex[1],plindex[2],plindex[6],plindex[5], YZ);
    }

    for (k = kmin; k <= kmax; k++)
    for (i = imin; i <= imax; i++)
    {
	j = jmax;
	index = d_index3d(i,j,k,pSolver->top_gmax);
	for (l = 0; l < 8; l++)
	    plindex[l] = pSolver->cell_center[index].corner[l]->pl_index;

	AddFaceEdge(plindex[2],plindex[3],plindex[7],plindex[6], XY);
    }

    vector<Graph_Point>::iterator itp;

    for (itp = vertices.begin(); itp != vertices.end(); itp++)
    {
	(itp->edge).sort();
	(itp->edge).unique();
    }
    
}


void CellCornerGraph::MergeBdryFrom(int dim, int dir, int rec_id)
{
    int tag = 1;
    int gmax[MAXD];
    int num_x, num_y, num_z;
    num_x = 0; num_y = 0; num_z = 0;
    int num_points;
    int *bdry;
    MPI_Status stat;
    int idim,i,j,k,l,index;

    for (idim = 0; idim < 3; idim++)
	gmax[idim] = pSolver->comp_grid->gmax[idim];

    switch(dim)
    {
	case 0: //In theta-direction
	    num_y = pSolver->comp_grid->gmax[1] + 1;
	    num_z = pSolver->comp_grid->gmax[2] + 1;
	    num_points = num_y * num_z;

	    bdry = new int[num_points];

	    MPI_Recv(bdry, num_points, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

	    switch(dir)
	    {
		case 0:
		    i = 0;
		    l = 0;
		    for (k = 0; k <= gmax[2]; k++)
		    for (j = 0; j <= gmax[1]; j++)
		    {
			index = d_index3d(i,j,k,gmax);
			if (bdry[l] >= vertices[index].set_index)
			{
			    l++;
			    continue;
			}
			else
			{
			    std::list<Graph_Point*>::iterator piter;
			    std::list<std::list<Graph_Point*> >::iterator pset = vertices[index].pSet;
			    for (piter = pset->begin(); piter != pset->end(); piter++)
			    {
				(*piter)->set_index = bdry[l];
			    }
			    l++;
			}
		    }
		    break;
		case 1:
		    i = gmax[0];
		    l = 0;
		    for (k = 0; k <= gmax[2]; k++)
		    for (j = 0; j <= gmax[1]; j++)
		    {
			index = d_index3d(i,j,k,gmax);
			if (bdry[l] >= vertices[index].set_index)
			{
			    l++;
			    continue;
			}
			else
			{
			    std::list<Graph_Point*>::iterator piter;
			    std::list<std::list<Graph_Point*> >::iterator pset = vertices[index].pSet;
			    for (piter = pset->begin(); piter != pset->end(); piter++)
			    {
				(*piter)->set_index = bdry[l];
			    }
			    l++;
			}
		    }
		    break;
	    }
	    delete[] bdry;
	    bdry = 0;
	    break;
	case 1: //In z-direction
	    num_x = pSolver->comp_grid->gmax[0] + 1;
	    num_z = pSolver->comp_grid->gmax[2] + 1;
	    num_points = num_x * num_z;

	    bdry = new int[num_points];


	    MPI_Recv(bdry, num_points, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

	    switch(dir)
	    {
		case 0:
		    j = 0;
		    l = 0;
		    for (k = 0; k <= gmax[2]; k++)
		    for (i = 0; i <= gmax[0]; i++)
		    {
			index = d_index3d(i,j,k,gmax);
			if (bdry[l] >= vertices[index].set_index)
			{
			    l++;
			    continue;
			}
			else
			{
			    std::list<Graph_Point*>::iterator piter;
			    std::list<std::list<Graph_Point*> >::iterator pset = vertices[index].pSet;
			    for (piter = pset->begin(); piter != pset->end(); piter++)
			    {
				(*piter)->set_index = bdry[l];
			    }
			    l++;
			}
		    }
		    break;
		case 1:
		    j = gmax[1];
		    l = 0;
		    for (k = 0; k <= gmax[2]; k++)
		    for (i = 0; i <= gmax[0]; i++)
		    {
			index = d_index3d(i,j,k,gmax);
			if (bdry[l] >= vertices[index].set_index)
			{
			    l++;
			    continue;
			}
			else
			{
			    std::list<Graph_Point*>::iterator piter;
			    std::list<std::list<Graph_Point*> >::iterator pset = vertices[index].pSet;
			    for (piter = pset->begin(); piter != pset->end(); piter++)
			    {
				(*piter)->set_index = bdry[l];
			    }
			    l++;
			}
		    }
		    break;
	    }

	    delete[] bdry;
	    bdry = 0;
	    break;
	case 2: //In r-direction
	    num_x = pSolver->comp_grid->gmax[0] + 1;
	    num_y = pSolver->comp_grid->gmax[1] + 1;
	    num_points = num_x * num_y;

	    bdry = new int[num_points];
	    
	    MPI_Recv(bdry, num_points, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

	    switch(dir)
	    {
		case 0:
		    k = 0;
		    l = 0;
		    for (j = 0; j <= gmax[1]; j++)
		    for (i = 0; i <= gmax[0]; i++)
		    {
			index = d_index3d(i,j,k,gmax);
			if (bdry[l] >= vertices[index].set_index)
			{
			    l++;
			    continue;
			}
			else
			{
			    std::list<Graph_Point*>::iterator piter;
			    std::list<std::list<Graph_Point*> >::iterator pset = vertices[index].pSet;
			    for (piter = pset->begin(); piter != pset->end(); piter++)
			    {
				(*piter)->set_index = bdry[l];
			    }
			    l++;
			}
		    }
		    break;
		case 1:
		    k = gmax[2];
		    l = 0;
		    for (j = 0; j <= gmax[1]; j++)
		    for (i = 0; i <= gmax[0]; i++)
		    {
			index = d_index3d(i,j,k,gmax);
			if (bdry[l] >= vertices[index].set_index)
			{
			    l++;
			    continue;
			}
			else
			{
			    std::list<Graph_Point*>::iterator piter;
			    std::list<std::list<Graph_Point*> >::iterator pset = vertices[index].pSet;
			    for (piter = pset->begin(); piter != pset->end(); piter++)
			    {
				(*piter)->set_index = bdry[l];
			    }
			    l++;
			}
		    }
		    break;
	    }

	    delete[] bdry;
	    bdry = 0;
	    break;
	default:
	    printf("\n Wrong dimension! \n");
	    clean_up(ERROR);
    }

}


void CellCornerGraph::MakeBdryForSend(int dim, int dir, int *&bdry)
{
    int i,j,k,index,l,idim;
    int gmax[MAXD];
    int num_x, num_y, num_z;
    num_x = 0; num_y = 0; num_z = 0;
    int num_points;

    for (idim = 0; idim < 3; idim++)
	gmax[idim] = pSolver->comp_grid->gmax[idim];

    switch(dim)
    {
	case 0:
	    num_y = pSolver->comp_grid->gmax[1] + 1;
	    num_z = pSolver->comp_grid->gmax[2] + 1;
	    num_points = num_y * num_z;

	    bdry = new int[num_points];

	    switch(dir)
	    {
		case 0:
		    i = 0;
		    l = 0;
		    for (k = 0; k <= gmax[2]; k++)
		    for (j = 0; j <= gmax[1]; j++)
		    {
			index = d_index3d(i,j,k,gmax);
			bdry[l] = vertices[index].set_index;
			l++;
		    }
		    break;
		case 1:
		    i = gmax[0];
		    l = 0;
		    for (k = 0; k <= gmax[2]; k++)
		    for (j = 0; j <= gmax[1]; j++)
		    {
			index = d_index3d(i,j,k,gmax);
			bdry[l] = vertices[index].set_index;
			l++;
		    }
		    break;
		default:
		    printf("\nWrong Direction!\n");
		    clean_up(ERROR);
		    break;
	    }
	    break;
	case 1:
	    num_x = pSolver->comp_grid->gmax[0] + 1;
	    num_z = pSolver->comp_grid->gmax[2] + 1;
	    num_points = num_x * num_z;

	    bdry = new int[num_points];

	    switch(dir)
	    {
		case 0:
		    j = 0;
		    l = 0;
		    for (k = 0; k <= gmax[2]; k++)
		    for (i = 0; i <= gmax[0]; i++)
		    {
			index = d_index3d(i,j,k,gmax);
			bdry[l] = vertices[index].set_index;
			l++;
		    }
		    break;
		case 1:
		    j = gmax[1];
		    l = 0;
		    for (k = 0; k <= gmax[2]; k++)
		    for (i = 0; i <= gmax[0]; i++)
		    {
			index = d_index3d(i,j,k,gmax);
			bdry[l] = vertices[index].set_index;
			l++;
		    }
		    break;
		default:
		    printf("\nWrong Direction!\n");
		    clean_up(ERROR);
		    break;
	    }
	    break;
	case 2:
	    num_x = pSolver->comp_grid->gmax[0] + 1;
	    num_y = pSolver->comp_grid->gmax[1] + 1;
	    num_points = num_x * num_y;

	    bdry = new int[num_points];

	    switch(dir)
	    {
		case 0:
		    k = 0;
		    l = 0;
		    for (j = 0; j <= gmax[1]; j++)
		    for (i = 0; i <= gmax[0]; i++)
		    {
			index = d_index3d(i,j,k,gmax);
			bdry[l] = vertices[index].set_index;
			l++;
		    }
		    break;
		case 1:
		    k = gmax[2];
		    l = 0;
		    for (j = 0; j <= gmax[1]; j++)
		    for (i = 0; i <= gmax[0]; i++)
		    {
			index = d_index3d(i,j,k,gmax);
			bdry[l] = vertices[index].set_index;
			l++;
		    }
		    break;
		default:
		    printf("\nWrong Direction!\n");
		    clean_up(ERROR);
		    break;
	    }
	    break;
	default:
	    printf("\n Wrong dimension! \n");
	    clean_up(ERROR);
	    break;
    }

}
void CellCornerGraph::SendBdryTo(int dim, int dir, int dst_id, int *&bdry)
{
    MPI_Request request;
    int idim;
    int tag = 1;
    int gmax[MAXD];
    int num_points;
    int num_x,num_y,num_z;

    for (idim = 0; idim < 3; idim++)
	gmax[idim] = pSolver->comp_grid->gmax[idim];

    switch(dim)
    {
	case 0:
	    num_y = pSolver->comp_grid->gmax[1] + 1;
	    num_z = pSolver->comp_grid->gmax[2] + 1;
	    num_points = num_y * num_z;
	    break;
	case 1:
	    num_x = pSolver->comp_grid->gmax[0] + 1;
	    num_z = pSolver->comp_grid->gmax[2] + 1;
	    num_points = num_x * num_z;
	    break;
	case 2:
	    num_x = pSolver->comp_grid->gmax[0] + 1;
	    num_y = pSolver->comp_grid->gmax[1] + 1;
	    num_points = num_x * num_y;
	    break;
	default:
	    printf("\n Wrong dimension! \n");
	    clean_up(ERROR);
	    break;
    }
    MPI_Isend(bdry, num_points, MPI_INT, dst_id, tag, MPI_COMM_WORLD, &request);
}

void CellCornerGraph::MergeSet(void)
{
    int i,j,k;
    int bdry_type[MAXD][2];
    int dim = pSolver->dim;
    MPI_Request request;
    MPI_Status status;

    int *G = pSolver->pp_grid->gmax;

    int me[MAXD], him[MAXD];

    int my_id, dst_id, rcv_id;

    INTERFACE *intfc = pSolver->front->interf;

    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < 2; j++)
	{
	    bdry_type[i][j] = rect_boundary_type(intfc,i,j);

	}
    }

    my_id = pp_mynode();
    find_Cartesian_coordinates(my_id, pSolver->pp_grid, me);


    for (i = 0; i < dim; i++)
    {
	for (int l = 0; l < G[i] - 1; l++)
	{
	    j = 1;

	    int *bdry_send = 0;
	    for (k = 0; k < dim; k++)
		him[k] = me[k];

	    pp_gsync();
	    if (bdry_type[i][j] == SUBDOMAIN_BOUNDARY)
	    {
		MakeBdryForSend(i, j, bdry_send);

		him[i] = me[i] + 2*j - 1;
		him[i] = (him[i] + G[i]) % G[i];
		dst_id = domain_id(him,G,dim);
		SendBdryTo(i, j, dst_id, bdry_send);
	    }

	    if (bdry_type[i][((j+1)%2)] == SUBDOMAIN_BOUNDARY)
	    {
		him[i] = me[i] - 2*j + 1;
		him[i] = (him[i] + G[i]) % G[i];
		rcv_id = domain_id(him,G,dim);
		MergeBdryFrom(i, ((j+1)%2), rcv_id);
	    }
	    pp_gsync();

	    delete[] bdry_send;

	}
    }

}

void CellCornerGraph::getPhyComp2(void)
{

    phy_comp.clear();
    int i,j,k,index,idim,sindex;
    int gmax[MAXD];
    map<int, pair<int, int> >::iterator itfind;
    pair<int, pair<int, int> > item;

    for (idim = 0; idim < 3; idim++)
	gmax[idim] = pSolver->comp_grid->gmax[idim];

    for (k = 0; k <= gmax[2]; k++)
    for (j = 0; j <= gmax[1]; j++)
    for (i = 0; i <= gmax[0]; i++)
    {
	index = d_index3d(i,j,k,gmax);
	sindex = vertices[index].set_index;

	itfind = phy_comp.find(sindex);
	if (itfind == phy_comp.end())
	{
	    item.first = sindex;
	    item.second.first = 1;
	    item.second.second = vertices[index].comp;
	    phy_comp.insert(item);
	}
	else
	{
	    itfind->second.first = itfind->second.first + 1;
	    itfind->second.second = vertices[index].comp;
	}

    }

}

void CellCornerGraph::getPhyComp(void)
{
    phy_comp.clear();
    std::list<std::list<Graph_Point*> >::iterator siter;
    std::list<Graph_Point*>::iterator piter;
    std::pair<int, pair<int, int> > elem;

    for (siter = connected_sets.begin(); siter != connected_sets.end(); siter++)
    {
	piter = siter->begin();
	elem.first = (*piter)->set_index;
	elem.second.first = (int) siter->size();
	elem.second.second = (*piter)->comp;
	phy_comp.insert(elem);
    }
}

void CellCornerGraph::getBdryComp(void)
{

    bdry_comp.clear();

    int my_id;
    int me[MAXD];
    map<int, pair<int, int> >::iterator itfind;
    pair<int, pair<int, int> > item;

    my_id = pp_mynode();
    find_Cartesian_coordinates(my_id, pSolver->pp_grid, me);


    int i,j,k,index,idim;
    int sindex;
    int gmax[MAXD];

    for (idim = 0; idim < 3; idim++)
	gmax[idim] = pSolver->comp_grid->gmax[idim];


    if (me[0] != 0 && me[1] != 0 && me[2] !=0)
    {
	i = 0;
	for (k = 0; k <= gmax[2]; k++)
	for (j = 0; j <= gmax[1]; j++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }

	}
	j = 0;
	for (i = 0; i <= gmax[0]; i++)
	for (k = 0; k <= gmax[2]; k++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }

	}
	k = 0;
	for (j = 0; j <= gmax[1]; j++)
	for (i = 0; i <= gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }

	}
	i = 0;
	j = 0;
	for (k = 0; k <= gmax[2]; k++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;
	    ((bdry_comp[sindex]).first)--;
	}
	i = 0;
	k = 0;
	for (j = 0; j <= gmax[1]; j++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;
	    ((bdry_comp[sindex]).first)--;
	}
	j = 0;
	k = 0;
	for (i = 0; i <= gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;
	    ((bdry_comp[sindex]).first)--;
	}
	i = 0; j = 0; k = 0;
	index = d_index3d(i,j,k,gmax);
	sindex = vertices[index].set_index;
	((bdry_comp[sindex]).first)++;
    }
    else if (me[0] == 0 && me[1] != 0 && me[2] != 0)
    {
	j = 0;
	for (i = 0; i <= gmax[0]; i++)
	for (k = 0; k <= gmax[2]; k++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }

	}
	k = 0;
	for (j = 0; j <= gmax[1]; j++)
	for (i = 0; i <= gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }

	}
	j = 0;
	k = 0;
	for (i = 0; i <= gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;
	    ((bdry_comp[sindex]).first)--;
	}
    }
    else if (me[0] != 0 && me[1] == 0 && me[2] != 0)
    {
	i = 0;
	for (k = 0; k <= gmax[2]; k++)
	for (j = 0; j <= gmax[1]; j++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }

	}
	k = 0;
	for (j = 0; j <= gmax[1]; j++)
	for (i = 0; i <= gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }
	}
	i = 0;
	k = 0;
	for (j = 0; j <= gmax[1]; j++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;
	    ((bdry_comp[sindex]).first)--;
	}
    }
    else if (me[0] == 0 && me[1] == 0 && me[2] != 0)
    {
	k = 0;
	for (j = 0; j <= gmax[1]; j++)
	for (i = 0; i <= gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }
	}
    }
    else if (me[0] != 0 && me[1] != 0 && me[2] == 0)
    {
	i = 0;
	for (k = 0; k <= gmax[2]; k++)
	for (j = 0; j <= gmax[1]; j++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }
	}
	j = 0;
	for (i = 0; i <= gmax[0]; i++)
	for (k = 0; k <= gmax[2]; k++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }
	}
	i = 0;
	j = 0;
	for (k = 0; k <= gmax[2]; k++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;
	    ((bdry_comp[sindex]).first)--;
	}
    }
    else if (me[0] == 0 && me[1] != 0 && me[2] == 0)
    {
	j = 0;
	for (i = 0; i <= gmax[0]; i++)
	for (k = 0; k <= gmax[2]; k++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }
	}
    }
    else if (me[0] != 0 && me[1] == 0 && me[2] == 0)
    {
	i = 0;
	for (k = 0; k <= gmax[2]; k++)
	for (j = 0; j <= gmax[1]; j++)
	{
	    index = d_index3d(i,j,k,gmax);
	    sindex = vertices[index].set_index;

	    itfind = bdry_comp.find(sindex);
	    if (itfind == bdry_comp.end())
	    {
		item.first = sindex;
		item.second.first = 1;
		item.second.second = vertices[index].comp;
		bdry_comp.insert(item);
	    }
	    else
	    {
		itfind->second.first = itfind->second.first + 1;
		itfind->second.second = vertices[index].comp;
	    }
	}
    }
}

void CellCornerGraph::SendCompTo(int dst_id, int dim)
{

    printf("\nI'm node %d, I'm sending comp to %d\n", pp_mynode(), dst_id);

    int tag = 1;
    int *num_comps_send;
    int ncomps = phy_comp.size();
    num_comps_send = new int (ncomps);
    printf("\nBefore sending numbers\n");

    MPI_Send(num_comps_send, 1, MPI_INT, dst_id, tag, MPI_COMM_WORLD);

    printf("\n Finished sending the number %d \n", (*num_comps_send));
    int *l_set_index_send;
    int *l_set_size_send;
    int *l_set_comp_send;

    printf ("\nBefore making the array\n");

    l_set_index_send = new int [(*num_comps_send)];
    l_set_size_send = new int [(*num_comps_send)];
    l_set_comp_send = new int [(*num_comps_send)];

    printf ("\nFinished making the array\n");

    int i = 0;
    map<int, pair<int, int> >::iterator iphy;

    for (iphy = phy_comp.begin(); iphy != phy_comp.end(); iphy++)
    {
	l_set_index_send[i] = iphy->first;
	l_set_size_send[i] = iphy->second.first;
	l_set_comp_send[i] = iphy->second.second;
	i++;
    }

    printf("\nFinishing filling the array to send\n");
    MPI_Send(l_set_index_send, (*num_comps_send), MPI_INT, dst_id, tag, MPI_COMM_WORLD);
    printf("\n Finished sending the index \n");
    MPI_Send(l_set_size_send, (*num_comps_send), MPI_INT, dst_id, tag, MPI_COMM_WORLD);
    printf("\n Finished sending the size \n");
    MPI_Send(l_set_comp_send, (*num_comps_send), MPI_INT, dst_id, tag, MPI_COMM_WORLD);
    printf("\n Finished sending the fluid comp \n");
    delete num_comps_send;
    delete [] l_set_index_send;
    delete [] l_set_size_send;
    delete [] l_set_comp_send;

    printf("\nI'm node %d, I finished sending comp to %d\n", pp_mynode(), dst_id);
}

void CellCornerGraph::MergeCompFrom(int rec_id, int dim)
{
    printf("\nI'm node %d, I'm receiving comp from %d\n", pp_mynode(), rec_id);

    
    int tag = 1;
    int *num_comps_recv = new int (2);
    MPI_Status stat;
    MPI_Recv(num_comps_recv, 1, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

    printf("\nFinished receiving the number %d\n", (*num_comps_recv));
    int *l_set_index_recv;
    int *l_set_size_recv;
    int *l_set_comp_recv;

    l_set_index_recv = new int [(*num_comps_recv)];
    l_set_size_recv = new int [(*num_comps_recv)];
    l_set_comp_recv = new int [(*num_comps_recv)];

    printf("\nfinished making the space\n");
    MPI_Recv(l_set_index_recv, (*num_comps_recv), MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    printf("\nFinished receiving the index\n");
    MPI_Recv(l_set_size_recv, (*num_comps_recv), MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    printf("\nFinished receiving the size\n");
    MPI_Recv(l_set_comp_recv, (*num_comps_recv), MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    printf("\nFinished receiving the comp\n");

    for (int i = 0; i < (*num_comps_recv); i++)
    {
	int set_index = l_set_index_recv[i];
	map<int, pair<int, int> >::iterator itfind = phy_comp.find(set_index);
	pair<int, pair<int, int> > item;
	if (itfind == phy_comp.end())
	{
	    item.first = set_index;
	    item.second.first = l_set_size_recv[i];
	    item.second.second = l_set_comp_recv[i];
	    phy_comp.insert(item);
	}
	else
	{
	    itfind->second.first += l_set_size_recv[i];
	    itfind->second.second = l_set_comp_recv[i];
	}
    }
    delete num_comps_recv;
    delete [] l_set_index_recv;
    delete [] l_set_size_recv;
    delete [] l_set_comp_recv;

    printf("\nI'm node %d, I finished receiving comp from %d\n", pp_mynode(), rec_id);
}

void CellCornerGraph::outputHistogram(char* out_name)
{
    char filename[256];
    FILE *outfile;
    map<int, pair<int, int> >::iterator itphy;
    int my_id = pp_mynode();

    int i = 1;
    if (my_id == 0)
    {
	sprintf(filename, "%s/connect-log", out_name);
	if(pSolver->front->step == 0)
	    outfile = fopen(filename, "w");
	else
	{
	    outfile = fopen(filename, "a");
	}

	(void) fprintf(outfile, "\n\n\ntime = %20.16g   step = %5d\n",
		pSolver->front->time, pSolver->front->step);

	fprintf(outfile, "\nNumber of connected physical components: %d", (int) phy_comp.size());
	fprintf(outfile, "\nSize and fluid component list:");
	for (itphy = phy_comp.begin(); itphy != phy_comp.end(); itphy++)
	{
	    fprintf(outfile, "\nThe %d-th physical component:", i);
	    fprintf(outfile, "\nSize = %d, Fluid comp = %d", itphy->second.first, itphy->second.second);
	    i++;
	}
    }
    printf("\n");
}

void CellCornerGraph::MergeToSingleProc(void)
{
    pp_gsync();
    int i,j,k;
    int *G = pSolver->pp_grid->gmax;
    int dim = pSolver->dim;

    int my_id, dst_id, rcv_id, act_id;

    int me[MAXD], him[MAXD], prev[MAXD];

    my_id = pp_mynode();

    // First Merge in X-direction

    if (G[0] > 1)
    {
	for (j = 0; j < G[1]; j++)
	for (k = 0; k < G[2]; k++)
	{
	    me[0] = G[0] - 1;
	    me[1] = j;
	    me[2] = k;
	    him[0] = me[0] - 1;
	    him[1] = me[1];
	    him[2] = me[2];
	    act_id = domain_id(me,  G, dim);
	    dst_id = domain_id(him, G, dim);

	    if(my_id == act_id)
		SendCompTo(dst_id, 0);

	    for (i = G[0] - 2; i > 0; i--)
	    {
		me[0] = i;
		me[1] = j;
		me[2] = k;
		him[0] = i - 1;
		him[1] = j;
		him[2] = k;
		prev[0] = i + 1;
		prev[1] = j;
		prev[2] = k;

		act_id = domain_id(me,  G, dim);
		dst_id = domain_id(him, G, dim);
		rcv_id = domain_id(prev,G, dim);

		if (my_id == act_id)
		{
		    MergeCompFrom(rcv_id, 0);
		    SendCompTo(dst_id, 0);
		}
	    }

	    me[0] = 0;
	    me[1] = j;
	    me[2] = k;
	    prev[0] = 1;
	    prev[1] = j;
	    prev[2] = k;

	    act_id = domain_id(me,  G, dim);
	    rcv_id = domain_id(prev,G, dim);

	    if (my_id == act_id)
		MergeCompFrom(rcv_id, 0);
	}
    }

    pp_gsync();

    printf("\nFinished merging in x-direction\n");

    // Then merge in Y-direction

    if (G[1] > 1)
    {
	for (k = 0; k < G[2]; k++)
	{
	    me[0] = 0;
	    me[1] = G[1] - 1;
	    me[2] = k;
	    him[0] = me[0];
	    him[1] = me[1] - 1;
	    him[2] = me[2];
	    act_id = domain_id(me,  G, dim);
	    dst_id = domain_id(him, G, dim);

	    if (my_id == act_id)
		SendCompTo(dst_id, 1);

	    for (j = G[1] - 2; j > 0; j--)
	    {
		me[0] = 0;
		me[1] = j;
		me[2] = k;
		him[0] = me[0];
		him[1] = me[1] - 1;
		him[2] = me[2];
		prev[0] = me[0];
		prev[1] = me[1] + 1;
		prev[2] = me[2];

		act_id = domain_id(me,  G, dim);
		dst_id = domain_id(him, G, dim);
		rcv_id = domain_id(prev,G, dim);

		if (my_id == act_id)
		{
		    MergeCompFrom(rcv_id, 1);
		    SendCompTo(dst_id, 1);
		}
	    }

	    me[0] = 0;
	    me[1] = 0;
	    me[2] = k;
	    prev[0] = 0;
	    prev[1] = 1;
	    prev[2] = k;

	    act_id = domain_id(me,  G, dim);
	    rcv_id = domain_id(prev,G, dim);

	    if (my_id == act_id)
		MergeCompFrom(rcv_id, 1);
	}
    }

    pp_gsync();

    printf("\nFinished merging in y-direction\n");
    // Last Merge in Z-direction

    if (G[2] > 1)
    {
	me[0] = 0;
	me[1] = 0;
	me[2] = G[2] - 1;
	him[0] = me[0];
	him[1] = me[1];
	him[2] = me[2] - 1;
	act_id = domain_id(me,  G, dim);
	dst_id = domain_id(him, G, dim);

	if (my_id == act_id)
	    SendCompTo(dst_id, 2);

	for (k = G[2] - 2; k > 0; k--)
	{
	    me[0] = 0;
	    me[1] = 0;
	    me[2] = k;
	    him[0] = me[0];
	    him[1] = me[1];
	    him[2] = me[2] - 1;
	    prev[0] = me[0];
	    prev[1] = me[1];
	    prev[2] = me[2] + 1;

	    act_id = domain_id(me,  G, dim);
	    dst_id = domain_id(him, G, dim);
	    rcv_id = domain_id(prev,G, dim);

	    if (my_id == act_id)
	    {
		MergeCompFrom(rcv_id, 2);
		SendCompTo(dst_id, 2);
	    }
	}

	me[0] = 0;
	me[1] = 0;
	me[2] = 0;
	prev[0] = 0;
	prev[1] = 0;
	prev[2] = 1;

	act_id = domain_id(me,  G, dim);
	rcv_id = domain_id(prev,G, dim);

	if (my_id == act_id)
	    MergeCompFrom(rcv_id, 2);
    }


    pp_gsync();

    printf("\nFinished merging in z-direction\n");
}

void CellCornerGraph::Init(void)
{
    ChooseMarchingCube(false, true, false);
    ConstructVertex();
}

void CellCornerGraph::getRawData(void)
{
    ConstructEdge();
    MakeSet();
    FindConnectedComp();

    MergeSet();

    NumberOfComponents();
}

void CellCornerGraph::getPhyCompHistogram(char* out_name)
{
    start_clock("getPhyCom2");
    getPhyComp2();
    stop_clock("getPhyCom2");

    start_clock("adjustPhyCompForMerge");
    adjustPhyCompForMerge();
    stop_clock("adjustPhyCompForMerge");

    start_clock("MergeToSingleProc");
    MergeToSingleProc();
    stop_clock("MergeToSingleProc");

    start_clock("outputHistogram");
    outputHistogram(out_name);
    stop_clock("outputHistogram");
}

void CellCornerGraph::adjustPhyCompForMerge(void)
{
    int my_id;
    int me[MAXD];
    int set_index;
    map<int, pair<int, int> >::iterator itbdry;

    getBdryComp();


    for(itbdry = bdry_comp.begin(); itbdry != bdry_comp.end(); itbdry++)
    {
	set_index = itbdry->first;
	(phy_comp[set_index]).first -= itbdry->second.first;
    }

    int i = 1;
    printf("\nAfter adjusting the comps:\n");
   for (map<int, pair<int, int> >::iterator itcomp = phy_comp.begin(); itcomp != phy_comp.end(); itcomp++)
   {
       printf("\nelement %d in phy_comp, key = %d, number = %d, fluid comp = %d\n", i, itcomp->first, itcomp->second.first, itcomp->second.second);
   }
}

void CellCornerGraph::outputCoordList(char* out_name)
{
   int my_id;
   char dirname[256];
   char filename[256];

   if (pp_numnodes() > 1)
       sprintf(dirname,"%s/P-%s",out_name,right_flush(pp_mynode(),4));
   else
       sprintf(dirname,"%s",out_name);
   sprintf(filename, "%s/connect-log", dirname);

   if (pp_numnodes() > 1)
       sprintf(filename, "%s-p%s", filename, right_flush(pp_mynode(), 4));

   FILE* outfile;
   outfile = fopen(filename, "w");

   list<list<Graph_Point*> >::iterator itsets;
   list<Graph_Point*>::iterator itp;

   for (itsets = connected_sets.begin(); itsets != connected_sets.end(); itsets++)
   {
       
       printf("Droplet index: %d, length = %d\n", (*(itsets->begin()))->set_index, (int) itsets->size());
       for (itp = itsets->begin(); itp != itsets->end(); itp++)
       {
	   int lindex = (*itp)->local_index;
	   double coord_theta = (pSolver->cell_corner[lindex]).p_coords[0];
	   double coord_z = (pSolver->cell_corner[lindex]).p_coords[1];
	   double coord_r = (pSolver->cell_corner[lindex]).p_coords[2];
	   printf("(%10.8g, %10.8g, %10.8g)\n", coord_theta, coord_z, coord_r);
       }
   }
}


void TriGraph::conTriVertices(int n_tris, int *g_index)
{
    num_v = n_tris;

    Graph_Point m_point;
    vertices.insert(vertices.end(), num_v, m_point);

    int i;
    for (i = 0; i < n_tris; i++)
    {
	vertices[i].local_index = i;
	vertices[i].global_index = g_index[i];
    }
}
void TriGraph::conTriEdges(int (*t_nb)[3])
{
    ClearEdge();

    int i,k;
    int i_nb0, i_nb1, i_nb2;

    for (i = 0; i < num_v; i++)
    {
	i_nb0 = t_nb[i][0];
	i_nb1 = t_nb[i][1];
	i_nb2 = t_nb[i][2];
	
	if(i_nb0 != -1)
	    vertices[i].edge.push_back(&(vertices[i_nb0]));
	if(i_nb1 != -1)
	    vertices[i].edge.push_back(&(vertices[i_nb1]));
	if(i_nb2 != -1)
	    vertices[i].edge.push_back(&(vertices[i_nb2]));
    }
}
