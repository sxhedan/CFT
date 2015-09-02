#include "Graph.h"
#include "iFDroplet.h"
#include "./volcomp/include/volumecompfnc.h"


void DropletAnalysis::InitFromInputFile_noMalloc(void)
{
    int i;
    char dirname[256];
    char *out_name = pSolver->front->out_name;

    sprintf(dirname, "%s", out_name);

    char filename_intfc[256];
    char filename_cell[256];

    sprintf(filename_intfc, "%s/intfc-nd%s", dirname, right_flush(pp_mynode(),4));
    sprintf(filename_cell, "%s/cell-nd%s", dirname, right_flush(pp_mynode(),4));

    FILE *infile_intfc = fopen(filename_intfc, "r");
    FILE *infile_cell = fopen(filename_cell, "r");

    /*
    fscanf(infile_intfc, "%d", &num_points);
    ps = new double [num_points][3];
  
    for (i = 0; i < num_points; i++) 
    {
      fscanf(infile_intfc, "%lf", &ps[i][0]);
      fscanf(infile_intfc, "%lf", &ps[i][1]);
      fscanf(infile_intfc, "%lf", &ps[i][2]);
    }
    */

    int num_points_read;
    fscanf(infile_intfc, "%d", &num_points_read);
    if (num_points_read != num_points)
    {
	printf("\nWrong num of points!\n");
	clean_up(ERROR);
    }

    double pt0, pt1, pt2;
    for (i = 0; i < num_points; i++) 
    {
      fscanf(infile_intfc, "%lf", &pt0);
      fscanf(infile_intfc, "%lf", &pt1);
      fscanf(infile_intfc, "%lf", &pt2);

      if ( (fabs(pt0 - ps[i][0]) > 1e-15) || (fabs(pt1 - ps[i][1]) > 1e-15) || (fabs(pt2 - ps[i][2]) > 1e-15) )
      {
	  printf("\nWrong Point\n");
	  clean_up(ERROR);
      }
    }
/*
    fscanf(infile_intfc, "%d", &num_tris);
    tris = new int [num_tris][3];
    tri_nb = new int [num_tris][3];
    comp_index = new int [num_tris];
    global_index = new int [num_tris];
    tri_interior = new bool [num_tris];
    for (i = 0; i < num_tris; i++) 
    {
	fscanf(infile_intfc, "%d", &tris[i][0]);
	fscanf(infile_intfc, "%d", &tris[i][1]);
	fscanf(infile_intfc, "%d", &tris[i][2]);
	fscanf(infile_intfc, "%d", &comp_index[i]);
	tri_nb[i][0] = tri_nb[i][1] = tri_nb[i][2] = -1;
	global_index[i] = -1;
	tri_interior[i] = true;
    }
    num_tris_interior = pSolver->Ntris;

    */

    int tri0,tri1,tri2,comp_indexi;
    int num_tris_read;
    fscanf(infile_intfc, "%d", &num_tris_read);

    if (num_tris != num_tris_read)
    {
	printf("\nWrong num of tris!\n");
	clean_up(ERROR);
    }


    for (i = 0; i < num_tris; i++) 
    {
	fscanf(infile_intfc, "%d", &tri0);
	fscanf(infile_intfc, "%d", &tri1);
	fscanf(infile_intfc, "%d", &tri2);
	fscanf(infile_intfc, "%d", &comp_indexi);

	if ( (tri0 != tris[i][0]) || (tri1 != tris[i][1]) || (tri2 != tris[i][2]) )
	{
	    printf("\nWrong tris\n");
	    clean_up(ERROR);
	}

    }
   /* 
    fscanf(infile_cell, "%d", &num_cells);
    c_centers = new double [num_cells][3];
    face_type = new int [num_cells][6];
    face_sign = new unsigned char [num_cells][6];
    gmax_cell = new int [3];

    for (i = 0; i < num_cells; i++)
    {
	fscanf(infile_cell, "%lf", &c_centers[i][0]);
	fscanf(infile_cell, "%lf", &c_centers[i][1]);
	fscanf(infile_cell, "%lf", &c_centers[i][2]);
	for (int j = 0; j < 6; j++) {
	    face_type[i][j] = -2;
	    face_sign[i][j] = 1;
	}
    }
    grid_spacing = new double [3];
    fscanf(infile_cell, "%lf", &(grid_spacing[0]));
    fscanf(infile_cell, "%lf", &(grid_spacing[1]));
    fscanf(infile_cell, "%lf", &(grid_spacing[2]));
    */
    int num_cells_read;
    fscanf(infile_cell, "%d", &num_cells_read);

    printf("\n Number of elements read from file = %d\n", num_cells_read);
    printf("\n Number of elements set from Init = %d\n", num_cells);
    if (num_cells != num_cells_read)
    {
	printf("\nWrong num of cells!\n");
	clean_up(ERROR);
    }

    for (i = 0; i < num_cells; i++)
    {
	
	double c0,c1,c2;
	fscanf(infile_cell, "%lf", &c0);
	fscanf(infile_cell, "%lf", &c1);
	fscanf(infile_cell, "%lf", &c2);
	if ( (fabs(c0 - c_centers[i][0]) > 1e-15) || (fabs(c1 - c_centers[i][1]) > 1e-15) || (fabs(c2 - c_centers[i][2]) > 1e-15) )
	{
	    printf("\nWrong Cell\n");
	    clean_up(ERROR);
	}
/*	
	fscanf(infile_cell, "%lf", &c_centers[i][0]);
	fscanf(infile_cell, "%lf", &c_centers[i][1]);
	fscanf(infile_cell, "%lf", &c_centers[i][2]);
	
	for (int j = 0; j < 6; j++) {
	    face_type[i][j] = -2;
	    face_sign[i][j] = 1;
	}
*/	
    }

    double gr0, gr1, gr2;
    fscanf(infile_cell, "%lf", &gr0);
    fscanf(infile_cell, "%lf", &gr1);
    fscanf(infile_cell, "%lf", &gr2);

    if ( (fabs(gr0 - grid_spacing[0]) > 1e-15) || (fabs(gr1 - grid_spacing[1]) > 1e-15) || (fabs(gr2 - grid_spacing[2]) > 1e-15) )
    {
	printf("\n Wrong spacing!\n");
	clean_up(ERROR);
    }
/*
    fscanf(infile_cell, "%lf", &(grid_spacing[0]));
    fscanf(infile_cell, "%lf", &(grid_spacing[1]));
    fscanf(infile_cell, "%lf", &(grid_spacing[2]));
*/
    fclose(infile_intfc);
    fclose(infile_cell);

    	map_comm.clear();
	map_comm_pos.clear();
	map_comm_neg.clear();
	vol_comm_pos.clear();
	vol_comm_neg.clear();
	blank_proc = false;
	first_non_blank_cell = 0;
	num_comps_global = 0;

    printf ("Finished Reading for no malloc !\n");

}

void DropletAnalysis::InitFromInputFile_binary(void)
{
    int i;
    char dirname[256];
    char *out_name = pSolver->front->out_name;

    sprintf(dirname, "%s", out_name);

    double val[3];
    int ival[4];

    char filename_intfc[256];
    char filename_cell[256];

    sprintf(filename_intfc, "%s/bintfc-nd%s", dirname, right_flush(pp_mynode(),4));
    sprintf(filename_cell, "%s/bcell-nd%s", dirname, right_flush(pp_mynode(),4));

    FILE *infile_intfc = fopen(filename_intfc, "rb");
    FILE *infile_cell = fopen(filename_cell, "rb");

    if (hardware_is_little_endian())
    {
	fread(ival, sizeof(int), 1, infile_intfc);
	num_points = endian_int_swap(ival[0]);
	ps = new double [num_points][3];

	for (i = 0; i < num_points; i++)
	{
	    fread(val, sizeof(double), 3, infile_intfc);
	    ps[i][0] = endian_double_swap(val[0]);
	    ps[i][1] = endian_double_swap(val[1]);
	    ps[i][2] = endian_double_swap(val[2]);
	}

	fread(ival, sizeof(int), 1, infile_intfc);
	num_tris = endian_int_swap(ival[0]);
	tris = new int [num_tris][3];
	tri_nb = new int [num_tris][3];
	comp_index = new int [num_tris];
	global_index = new int [num_tris];
	tri_interior = new bool [num_tris];

	for (i = 0; i < num_tris; i++)
	{
	    fread(ival, sizeof(int), 4, infile_intfc);
	    tris[i][0] = endian_int_swap(ival[0]);
	    tris[i][1] = endian_int_swap(ival[1]);
	    tris[i][2] = endian_int_swap(ival[2]);
	    comp_index[i] = endian_int_swap(ival[3]);
	    tri_nb[i][0] = tri_nb[i][1] = tri_nb[i][2] = -1;
	    global_index[i] = -1;
	    tri_interior[i] = true;
	}
	num_tris_interior = pSolver->Ntris;

	fread(ival, sizeof(int), 1, infile_cell);
	num_cells = endian_int_swap(ival[0]);
	c_centers = new double [num_cells][3];
	face_type = new int [num_cells][6];
	face_sign = new unsigned char [num_cells][6];
	for (i = 0; i < num_cells; i++)
	{
	    fread(val, sizeof(double), 3, infile_cell);
	    c_centers[i][0] = endian_double_swap(val[0]);
	    c_centers[i][1] = endian_double_swap(val[1]);
	    c_centers[i][2] = endian_double_swap(val[2]);
	    for (int j = 0; j < 6; j++) {
		face_type[i][j] = -2;
		face_sign[i][j] = 1;
	    }
	}

	grid_spacing = new double [3];
	fread(val, sizeof(double), 3, infile_cell);
	grid_spacing[0] = endian_double_swap(val[0]);
	grid_spacing[1] = endian_double_swap(val[1]);
	grid_spacing[2] = endian_double_swap(val[2]);

	fread(ival, sizeof(int), 1, infile_cell);
	num_corners = endian_int_swap(ival[0]);
	c_corners = new double [num_corners][3];
	for (i = 0; i < num_corners; i++)
	{
	    fread(val, sizeof(double), 3, infile_cell);
	    c_corners[i][0] = endian_double_swap(val[0]);
	    c_corners[i][1] = endian_double_swap(val[1]);
	    c_corners[i][2] = endian_double_swap(val[2]);
	}
	gmax_cell = new int [3];
	fread(ival, sizeof(int), 3, infile_cell);
	gmax_cell[0] = endian_int_swap(ival[0]);
	gmax_cell[1] = endian_int_swap(ival[1]);
	gmax_cell[2] = endian_int_swap(ival[2]);
    }	   
    else
    {
	fread(ival, sizeof(int), 1, infile_intfc);
	num_points = ival[0];
	ps = new double [num_points][3];

	for (i = 0; i < num_points; i++)
	{
	    fread(val, sizeof(double), 3, infile_intfc);
	    ps[i][0] = val[0];
	    ps[i][1] = val[1];
	    ps[i][2] = val[2];
	}

	fread(ival, sizeof(int), 1, infile_intfc);
	num_tris = ival[0];
	tris = new int [num_tris][3];
	tri_nb = new int [num_tris][3];
	comp_index = new int [num_tris];
	global_index = new int [num_tris];
	tri_interior = new bool [num_tris];

	for (i = 0; i < num_tris; i++)
	{
	    fread(ival, sizeof(int), 4, infile_intfc);
	    tris[i][0] = ival[0];
	    tris[i][1] = ival[1];
	    tris[i][2] = ival[2];
	    comp_index[i] = ival[3];
	    tri_nb[i][0] = tri_nb[i][1] = tri_nb[i][2] = -1;
	    global_index[i] = -1;
	    tri_interior[i] = true;
	}
	num_tris_interior = pSolver->Ntris;

	fread(ival, sizeof(int), 1, infile_cell);
	num_cells = ival[0];
	c_centers = new double [num_cells][3];
	face_type = new int [num_cells][6];
	face_sign = new unsigned char [num_cells][6];
	for (i = 0; i < num_cells; i++)
	{
	    fread(val, sizeof(double), 3, infile_cell);
	    c_centers[i][0] = val[0];
	    c_centers[i][1] = val[1];
	    c_centers[i][2] = val[2];
	    for (int j = 0; j < 6; j++) {
		face_type[i][j] = -2;
		face_sign[i][j] = 1;
	    }
	}

	grid_spacing = new double [3];
	fread(val, sizeof(double), 3, infile_cell);
	grid_spacing[0] = val[0];
	grid_spacing[1] = val[1];
	grid_spacing[2] = val[2];

	fread(ival, sizeof(int), 1, infile_cell);
	num_corners = ival[0];
	c_corners = new double [num_corners][3];
	for (i = 0; i < num_corners; i++)
	{
	    fread(val, sizeof(double), 3, infile_cell);
	    c_corners[i][0] = val[0];
	    c_corners[i][1] = val[1];
	    c_corners[i][2] = val[2];
	}
	gmax_cell = new int [3];
	fread(ival, sizeof(int), 3, infile_cell);
	gmax_cell[0] = ival[0];
	gmax_cell[1] = ival[1];
	gmax_cell[2] = ival[2];
    }	   

    

    fclose(infile_intfc);
    fclose(infile_cell);

    	map_comm.clear();
	map_comm_pos.clear();
	map_comm_neg.clear();
	vol_comm_pos.clear();
	vol_comm_neg.clear();
	blank_proc = false;
	first_non_blank_cell = 0;
	num_comps_global = 0;

    printf ("Finished Reading !\n");

}

void DropletAnalysis::InitFromInputFile(void)
{
    int i;
    char dirname[256];
    char *out_name = pSolver->front->out_name;

    sprintf(dirname, "%s", out_name);

    char filename_intfc[256];
    char filename_cell[256];

    sprintf(filename_intfc, "%s/intfc-nd%s", dirname, right_flush(pp_mynode(),4));
    sprintf(filename_cell, "%s/cell-nd%s", dirname, right_flush(pp_mynode(),4));

    FILE *infile_intfc = fopen(filename_intfc, "r");
    FILE *infile_cell = fopen(filename_cell, "r");

    fscanf(infile_intfc, "%d", &num_points);
    ps = new double [num_points][3];
  
    for (i = 0; i < num_points; i++) 
    {
      fscanf(infile_intfc, "%lf", &ps[i][0]);
      fscanf(infile_intfc, "%lf", &ps[i][1]);
      fscanf(infile_intfc, "%lf", &ps[i][2]);
    }
  
    fscanf(infile_intfc, "%d", &num_tris);
    tris = new int [num_tris][3];
    tri_nb = new int [num_tris][3];
    comp_index = new int [num_tris];
    global_index = new int [num_tris];
    tri_interior = new bool [num_tris];
    for (i = 0; i < num_tris; i++) 
    {
	fscanf(infile_intfc, "%d", &tris[i][0]);
	fscanf(infile_intfc, "%d", &tris[i][1]);
	fscanf(infile_intfc, "%d", &tris[i][2]);
	fscanf(infile_intfc, "%d", &comp_index[i]);
	tri_nb[i][0] = tri_nb[i][1] = tri_nb[i][2] = -1;
	global_index[i] = -1;
	tri_interior[i] = true;
    }
    num_tris_interior = pSolver->Ntris;

    fscanf(infile_cell, "%d", &num_cells);
    c_centers = new double [num_cells][3];
    face_type = new int [num_cells][6];
    face_sign = new unsigned char [num_cells][6]; 

    for (i = 0; i < num_cells; i++)
    {
	fscanf(infile_cell, "%lf", &c_centers[i][0]);
	fscanf(infile_cell, "%lf", &c_centers[i][1]);
	fscanf(infile_cell, "%lf", &c_centers[i][2]);

	for (int j = 0; j < 6; j++) {
	    face_type[i][j] = -2;
	    face_sign[i][j] = 1;
	}
    }

    grid_spacing = new double [3];
    fscanf(infile_cell, "%lf", &(grid_spacing[0]));
    fscanf(infile_cell, "%lf", &(grid_spacing[1]));
    fscanf(infile_cell, "%lf", &(grid_spacing[2]));

    fscanf(infile_cell, "%d", &num_corners);
    c_corners = new double [num_corners][3];
    for (i = 0; i < num_corners; i++)
    {
	fscanf(infile_cell, "%lf", &c_corners[i][0]);
	fscanf(infile_cell, "%lf", &c_corners[i][1]);
	fscanf(infile_cell, "%lf", &c_corners[i][2]);
    }
    gmax_cell = new int [3];

    fscanf(infile_cell, "%d", &(gmax_cell[0]));
    fscanf(infile_cell, "%d", &(gmax_cell[1]));
    fscanf(infile_cell, "%d", &(gmax_cell[2]));

    fclose(infile_intfc);
    fclose(infile_cell);

    	map_comm.clear();
	map_comm_pos.clear();
	map_comm_neg.clear();
	vol_comm_pos.clear();
	vol_comm_neg.clear();
	blank_proc = false;
	first_non_blank_cell = 0;
	num_comps_global = 0;

    printf ("Finished Reading !\n");

}
void DropletAnalysis::CreatePoints(void)
{

    INTERFACE *intfc = pSolver->tempfront->interf;

    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
    POINT *p;

    SURFACE **s;
    TRI *tri;

    int i = 0;

    num_points = pSolver->Npoints;

    ps = new double [num_points][3];

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    for (int k = 0; k < 3; k++)
	    {
		p = Point_of_tri(tri)[k];
		p->local_index = -1;
	    }
	}
    }

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    for (int k = 0; k < 3; k++)
	    {
		p = Point_of_tri(tri)[k];
		if (p->local_index == -1)
		{
		    ps[i][0] = p->_coords[0];
		    ps[i][1] = p->_coords[1];
		    ps[i][2] = p->_coords[2];
		    p->local_index = i;
		    i++;
		}
	    }
	}
    }
    
    
    /*
    next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
	ps[i][0] = p->_coords[0];
	ps[i][1] = p->_coords[1];
	ps[i][2] = p->_coords[2];
	i++;
    }
    
    */

}

void DropletAnalysis::adjustCompIndexForTris(void)
{
    int i;
    set<int> comps;
    for (i = 0; i < num_tris; i++)
	comps.insert(comp_index[i]);

    commSet(comps); //now comps is the global set

    num_comps_global = (int) comps.size();

    int k = 1;
    map<int,int> map_comm_global;

    set<int>::iterator its;

    for (its = comps.begin(); its != comps.end(); its++)
    {
	map_comm_global.insert( pair<int,int>((*its), k) );
	k++;
    }

    for (i = 0; i < num_tris; i++)
    {
	map<int,int>::iterator itm = map_comm_global.find(comp_index[i]);
	if (itm != map_comm_global.end())
	    comp_index[i] = itm->second;
    }

}

void DropletAnalysis::CreateMappings(void)
{
    int i;
    set<int> comps_proc;
    for (i = 0; i < num_tris; i++)
	comps_proc.insert(comp_index[i]);
    num_comps = comps_proc.size();

    //printf("\nnum_comps = %d\n", num_comps);

    comp_list = new int [num_comps];
    pos_map = new int [num_comps];
    neg_map = new int [num_comps];
    pos_vol = new double [num_comps];
    neg_vol = new double [num_comps];

    int p = 0;
    for (set<int>::iterator its = comps_proc.begin(); its != comps_proc.end(); its++)
    {
	comp_list[p] = (*its);
	p++;
    }

    if (p != num_comps)
    {
	printf("\nWrong number of comps!\n");
	clean_up(ERROR);
    }

    for (i = 0; i < num_comps; i++)
    {
	pos_map[i] = comp_list[i]; 
	neg_map[i] = comp_list[i];
	pos_vol[i] = 0.0;
	neg_vol[i] = 0.0;
    }
/*
    for (i = 0; i < num_comps; i++)
    {
	printf("comp_list[%d] = %d\n", i, comp_list[i]);
	printf("pos_map[%d] = %d\n", i, pos_map[i]);
	printf("neg_map[%d] = %d\n", i, neg_map[i]);
	printf("pos_vol[%d] = %10.8g\n", i, pos_vol[i]);
	printf("neg_vol[%d] = %10.8g\n", i, neg_vol[i]);
    }
*/
}

void DropletAnalysis::CreateTris(void)
{
    num_tris = pSolver->Ntris_all;
    num_tris_interior = pSolver->Ntris;

    INTERFACE *intfc = pSolver->tempfront->interf;

    SURFACE **s;
    TRI *tri;
    tris = new int [num_tris][3];
    tri_nb = new int [num_tris][3];
    comp_index = new int [num_tris];
    global_index = new int [num_tris];
    tri_interior = new bool [num_tris];

    int i = 0;

    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    global_index[i] = tri->global_index;
	    comp_index[i] = -2;

	    tris[i][0] = (Point_of_tri(tri)[0])->local_index;
	    tris[i][1] = (Point_of_tri(tri)[1])->local_index;
	    tris[i][2] = (Point_of_tri(tri)[2])->local_index;

	    if(is_side01_a_bond(tri))
		tri_nb[i][0] = -1;
	    else
		tri_nb[i][0] = Tri_on_side01(tri)->local_index;

	    if(is_side12_a_bond(tri))
		tri_nb[i][1] = -1;
	    else
		tri_nb[i][1] = Tri_on_side12(tri)->local_index;

	    if(is_side20_a_bond(tri))
		tri_nb[i][2] = -1;
	    else
		tri_nb[i][2] = Tri_on_side20(tri)->local_index;


	    if (tri->interior == TRUE)
		tri_interior[i] = true;
	    else
		tri_interior[i] = false;
	    
	    i++; 
	}
    }
}

void DropletAnalysis::CreateCells(void)
{
    int i,j,k,l,index;
    int imin = pSolver->imin;
    int imax = pSolver->imax;
    int jmin = pSolver->jmin;
    int jmax = pSolver->jmax;
    int kmin = pSolver->kmin;
    int kmax = pSolver->kmax;
    int t_max[3];
    t_max[0] = pSolver->top_gmax[0];
    t_max[1] = pSolver->top_gmax[1];
    t_max[2] = pSolver->top_gmax[2];

    gmax_cell = new int [3];

    gmax_cell[0] = imax - imin + 1;
    gmax_cell[1] = jmax - jmin + 1;
    gmax_cell[2] = kmax - kmin + 1;


    num_cells = (imax - imin + 1)*(jmax - jmin + 1)*(kmax - kmin + 1);
    num_corners = (imax - imin + 2)*(jmax - jmin + 2)*(kmax - kmin + 2);

    c_centers = new double [num_cells][3];
    face_type = new int [num_cells][6];
    face_sign = new unsigned char [num_cells][6];
    grid_spacing = new double [3];

    c_corners = new double [num_corners][3];

    l = 0;

    double proc_xmin = (pSolver->x_pp_bdry)[(pSolver->x_pp_index)];
    double proc_xmax = (pSolver->x_pp_bdry)[(pSolver->x_pp_index)+1];
    double proc_ymin = (pSolver->y_pp_bdry)[(pSolver->y_pp_index)];
    double proc_ymax = (pSolver->y_pp_bdry)[(pSolver->y_pp_index)+1];
    double proc_zmin = (pSolver->z_pp_bdry)[(pSolver->z_pp_index)];
    double proc_zmax = (pSolver->z_pp_bdry)[(pSolver->z_pp_index)+1];

    for (k = 0; k <= gmax_cell[2]; k++)
    for (j = 0; j <= gmax_cell[1]; j++)
    for (i = 0; i <= gmax_cell[0]; i++)
    {
	double cornerx, cornery, cornerz;

	if (i == 0)
	    cornerx = proc_xmin;
	else if (i == gmax_cell[0])
	    cornerx = proc_xmax;
	else
	    cornerx = proc_xmin + i*(proc_xmax - proc_xmin)/gmax_cell[0];

	if (j == 0)
	    cornery = proc_ymin;
	else if (j == gmax_cell[1])
	    cornery = proc_ymax;
	else
	    cornery = proc_ymin + j*(proc_ymax - proc_ymin)/gmax_cell[1];

	if (k == 0)
	    cornerz = proc_zmin;
	else if (k == gmax_cell[2])
	    cornerz = proc_zmax;
	else
	    cornerz = proc_zmin + k*(proc_zmax - proc_zmin)/gmax_cell[2];

	c_corners[l][0] = cornerx;
	c_corners[l][1] = cornery;
	c_corners[l][2] = cornerz;

	l++;
    }

    l = 0;
    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,t_max);

	c_centers[l][0] = (pSolver->cell_center)[index].m_coords[0];
	c_centers[l][1] = (pSolver->cell_center)[index].m_coords[1];
	c_centers[l][2] = (pSolver->cell_center)[index].m_coords[2];

	for (int p = 0; p < 6; p++)
	{
	    face_type[l][p] = -2;
	    face_sign[l][p] = 1;
	}
	l++;
    }

    grid_spacing[0] = pSolver->top_h[0];
    grid_spacing[1] = pSolver->top_h[1];
    grid_spacing[2] = pSolver->top_h[2];
}


void DropletAnalysis::DeleteMappings(void)
{
    delete[] comp_list;
    delete[] pos_map;
    delete[] neg_map;
    delete[] pos_vol;
    delete[] neg_vol;

    comp_list = 0; 
    pos_map = 0;
    neg_map = 0;
    pos_vol = 0;
    neg_vol = 0;
}

void DropletAnalysis::DeleteCells(void)
{

    delete[] gmax_cell;
    delete[] c_centers;
    delete[] face_type;
    delete[] face_sign;
    delete[] grid_spacing;
    delete[] c_corners;

    c_centers = 0;
    c_corners = 0;
    face_type = 0;
    face_sign = 0;
    grid_spacing = 0;
}

void DropletAnalysis::DeleteTris(void)
{
    delete[] comp_index;
    delete[] global_index;

    delete[] tris;
    delete[] tri_nb;
    delete[] tri_interior;

    comp_index = 0;
    global_index = 0;
    tris = 0;
    tri_nb = 0;
    tri_interior = 0;
}

void DropletAnalysis::DeletePoints(void)
{
    delete[] ps;
    ps = 0;
}

void DropletAnalysis::printAllTris(void)
{
    printf("\nInside DropletAnalysis, before finding connected components\n");
    printf("The number of triangles: %d\n", num_tris);
    printf("The number of interior triangles: %d\n", num_tris_interior);
    printf("The number of points: %d\n", num_points);

    printf("\n");
    for (int i = 0; i < num_tris; i++)
    {
	printf("%d-th triangle, global_index = %d, comp_index = %d, three indices:(%d, %d, %d), neighbour tris: (%d, %d, %d)\n",
		i,global_index[i],comp_index[i],tris[i][0],tris[i][1],tris[i][2],tri_nb[i][0],tri_nb[i][1],tri_nb[i][2]);
    }
}



void DropletAnalysis::printAllCells(void)
{
    printf("\nAll Cells:\n");

    printf("\ngmax_cell = (%d, %d, %d)\n", gmax_cell[0], gmax_cell[1], gmax_cell[2]);

    for (int i = 0; i < num_cells; i++)
    {
	double cell_corner[6];
	printf("%d-th cell, face_type = %d, %d, %d, %d, %d, %d\n", i, face_type[i][0], face_type[i][1], face_type[i][2], face_type[i][3], face_type[i][4], face_type[i][5]); 
	printf("%d-th cell, face_sign = %d, %d, %d, %d, %d, %d\n", i, face_sign[i][0], face_sign[i][1], face_sign[i][2], face_sign[i][3], face_sign[i][4], face_sign[i][5]); 
	printf("cell_center = (%10.8g, %10.8g, %10.8g)\n", c_centers[i][0], c_centers[i][1], c_centers[i][2]);
	//getCellCorners(i, cell_corner);
	//printf("xdir = (%10.8g, %10.8g), ydir = (%10.8g,  %10.8g), zdir = (%10.8g, %10.8g)\n", cell_corner[0], cell_corner[1], cell_corner[2], cell_corner[3], cell_corner[4], cell_corner[5]);
    }
}

void DropletAnalysis::GetCompIndexForTris(void)
{
    start_clock("getLocalComp");
    getLocalCompIndexForTris();
    stop_clock("getLocalComp");


    start_clock("commTriComp");
    getGlobalCompIndexForTris();
    stop_clock("commTriComp");


    start_clock("adjustTriComp");
    adjustCompIndexForTris();
    stop_clock("adjustTriComp");

}

void DropletAnalysis::adjustSurface(void)
{
  int i,j;
  double px,py,pz,dxmin,dxmax,dymin,dymax,dzmin,dzmax;
  double tol = 1e-02;

  for (i = 0; i < num_points; i++)
  {
    px = ps[i][0];
    py = ps[i][1];
    pz = ps[i][2];
    int count = 0;
  
    for (j = 0; j< num_cells; j++)
      {
        double cell_bnds[6];
        getCellCorners(j, cell_bnds);

        dxmin = fabs(px-cell_bnds[0])/(cell_bnds[1]-cell_bnds[0]);
        dxmax = fabs(px-cell_bnds[1])/(cell_bnds[1]-cell_bnds[0]);
        dymin = fabs(py-cell_bnds[2])/(cell_bnds[3]-cell_bnds[2]);
        dymax = fabs(py-cell_bnds[3])/(cell_bnds[3]-cell_bnds[2]);
        dzmin = fabs(pz-cell_bnds[4])/(cell_bnds[5]-cell_bnds[4]);
        dzmax = fabs(pz-cell_bnds[5])/(cell_bnds[5]-cell_bnds[4]);

        if (dxmin < tol)
          {
           ps[i][0] = cell_bnds[0];
          }
       else if (dxmax < tol)
          {
           ps[i][0] = cell_bnds[1];
	  }
        if (dymin < tol)
          {
	    ps[i][1] = cell_bnds[2];
	  }
        else if (dymax < tol)
          {
	    ps[i][1] = cell_bnds[3];
	  }
        if (dzmin < tol)
          {
	    ps[i][2] = cell_bnds[4];
	  }
        else if (dzmax < tol)
          {
	    ps[i][2] = cell_bnds[5];
	  }
      }
  }

}

void DropletAnalysis::GetLocalVolume(void)
{ 
    CreateMappings();

    adjustSurface();

    
    start_clock("getLocalVolume");

    for (int i = 0; i < num_cells; i++)
    {
	//printf ("%d-th cell \n", i);
	double cell_bnds[6];
	getCellCorners(i, cell_bnds);
	volumecompfnc(cell_bnds, face_type[i], face_sign[i], 
		ps[0], num_points, tris[0], num_tris, comp_index, num_comps, comp_list, 
		pos_map, neg_map, pos_vol, neg_vol);
    }
    
    stop_clock("getLocalVolume");
    
    printf("\n Number of comps = %d\n", num_comps);
}



void DropletAnalysis::fillBlankProcFaces(void)
{
    int i,j,k;
    int bdry_type[MAXD][2];
    int dim = pSolver->dim;
    MPI_Request request;
    MPI_Status status;

    int *G = pSolver->pp_grid->gmax;

    int me[MAXD], him[MAXD];

    int my_id, dst_id, rcv_id;

    INTERFACE *intfc = pSolver->tempfront->interf;

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

	    int *index_send = 0;
	    unsigned char *sign_send = 0;

	    for (k = 0; k < dim; k++)
		him[k] = me[k];

	    pp_gsync();
	    if (bdry_type[i][j] == SUBDOMAIN_BOUNDARY)
	    {
		makeIndexForSend(i, j, index_send, sign_send);

		him[i] = me[i] + 2*j - 1;
		him[i] = (him[i] + G[i]) % G[i];
		dst_id = domain_id(him,G,dim);
		sendIndexTo(i, j, dst_id, index_send, sign_send);
	    }

	    if (bdry_type[i][((j+1)%2)] == SUBDOMAIN_BOUNDARY)
	    {
		him[i] = me[i] - 2*j + 1;
		him[i] = (him[i] + G[i]) % G[i];
		rcv_id = domain_id(him,G,dim);
		mergeIndexFrom(i, ((j+1)%2), rcv_id);
	    }
	    pp_gsync();

	    delete index_send;
	    delete sign_send;

	}
    }// First direction

    for (i = 0; i < dim; i++)
    {
	for (int l = 0; l < G[i] - 1; l++)
	{
	    j = 0;

	    int *index_send = 0;
	    unsigned char *sign_send = 0;

	    for (k = 0; k < dim; k++)
		him[k] = me[k];

	    pp_gsync();
	    if (bdry_type[i][j] == SUBDOMAIN_BOUNDARY)
	    {
		makeIndexForSend(i, j, index_send, sign_send);

		him[i] = me[i] + 2*j - 1;
		him[i] = (him[i] + G[i]) % G[i];
		dst_id = domain_id(him,G,dim);
		sendIndexTo(i, j, dst_id, index_send, sign_send);
	    }

	    if (bdry_type[i][((j+1)%2)] == SUBDOMAIN_BOUNDARY)
	    {
		him[i] = me[i] - 2*j + 1;
		him[i] = (him[i] + G[i]) % G[i];
		rcv_id = domain_id(him,G,dim);
		mergeIndexFrom(i, ((j+1)%2), rcv_id);
	    }
	    pp_gsync();

	    delete index_send;
	    delete sign_send;

	}
    }// Second direction

}


void DropletAnalysis::sendIndexTo(int dim, int dir, int dst_id, int *&index_for_send, unsigned char *&sign_for_send)
{
    int tag = 1;

    MPI_Request request_index;
    MPI_Request request_sign;

    MPI_Isend(index_for_send, 1, MPI_INT, dst_id, tag, MPI_COMM_WORLD, &request_index);
    MPI_Isend(sign_for_send, 1, MPI_UNSIGNED_CHAR, dst_id, tag, MPI_COMM_WORLD, &request_sign);
}

void DropletAnalysis::mergeIndexFrom(int dim, int dir, int rec_id)
{
    int tag = 1;
    MPI_Status stat;

    int *index_recv = new int;
    unsigned char *sign_recv = new unsigned char;

    MPI_Recv(index_recv, 1, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(sign_recv, 1, MPI_UNSIGNED_CHAR, rec_id, tag, MPI_COMM_WORLD, &stat);

    if (blank_proc == true)
    {
	for (int i = 0; i < num_cells; i++)
	{
	    for (int k = 0; k < 6; k++)
	    {
		face_type[i][k] = (*index_recv);
		face_sign[i][k] = (*sign_recv);
	    }
	}
    }

    delete index_recv;
    delete sign_recv;
}

void DropletAnalysis::makeIndexForSend(int dim, int dir, int *&index_send, unsigned char *&sign_send)
{
    index_send = new int;
    sign_send = new unsigned char;

    int i,j,k;
    int index;

    switch(dim)
    {
	case 0:
	    if (dir == 0)
	    {
		i = 0; j = 0; k = 0;
		index = getCellIndex(i,j,k);
		(*index_send) = face_type[index][3];
		(*sign_send) = face_sign[index][3];
	    }
	    else if (dir == 1)
	    {
		i = gmax_cell[0] - 1; j = 0; k = 0;
		index = getCellIndex(i,j,k);
		(*index_send) = face_type[index][1];
		(*sign_send) = face_sign[index][1];
	    }
	    break;
	case 1:
	    if (dir == 0)
	    {
		i = 0; j = 0; k = 0;
		index = getCellIndex(i,j,k);
		(*index_send) = face_type[index][4];
		(*sign_send) = face_sign[index][4];
	    }
	    else if (dir == 1)
	    {
		i = 0; j = gmax_cell[1] - 1; k = 0;
		index = getCellIndex(i,j,k);
		(*index_send) = face_type[index][2];
		(*sign_send) = face_sign[index][2];
	    }
	    break;
	case 2:
	    if (dir == 0)
	    {
		i = 0; j = 0; k = 0;
		index = getCellIndex(i,j,k);
		(*index_send) = face_type[index][0];
		(*sign_send) = face_sign[index][0];
	    }
	    else if (dir == 1)
	    {
		i = 0; j = 0; k = gmax_cell[2] - 1;
		index = getCellIndex(i,j,k);
		(*index_send) = face_type[index][5];
		(*sign_send) = face_sign[index][5];
	    }
	    break;
	default:
	    printf("\nWrong dimension!\n");
	    clean_up(ERROR);
	    break;
    }
}

void DropletAnalysis::getBlankProcVolume(void)
{
    vol_comm_pos.clear();
    vol_comm_neg.clear();

    double proc_vol = 0.0;

    for (int i = 0; i < num_cells; i++)
    {
	double cell_vol = c_centers[i][2]*grid_spacing[0]*grid_spacing[1]*grid_spacing[2];
	proc_vol = proc_vol + cell_vol;
    }

    if (face_sign[0][0] != 0)
	vol_comm_pos.insert(pair<int,double>(face_type[0][0], proc_vol));
    else
	vol_comm_neg.insert(pair<int,double>(face_type[0][0], proc_vol));
}

void DropletAnalysis::mergeVolumeInfo(map<int,double> &vol_for_comm)
{
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
		sendVolumeTo(dst_id, 0, vol_for_comm);

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
		    mergeVolumeFrom(rcv_id, 0, vol_for_comm);
		    sendVolumeTo(dst_id, 0, vol_for_comm);
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
		mergeVolumeFrom(rcv_id, 0, vol_for_comm);
	}
    }

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
		sendVolumeTo(dst_id, 1, vol_for_comm);

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
		    mergeVolumeFrom(rcv_id, 1, vol_for_comm);
		    sendVolumeTo(dst_id, 1, vol_for_comm);
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
		mergeVolumeFrom(rcv_id, 1, vol_for_comm);
	}
    }

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
	    sendVolumeTo(dst_id, 2, vol_for_comm);

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
		mergeVolumeFrom(rcv_id, 2, vol_for_comm);
		sendVolumeTo(dst_id, 2, vol_for_comm);
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
	    mergeVolumeFrom(rcv_id, 2, vol_for_comm);
    }

}

void DropletAnalysis::adjustVolumeForBlankCells(void)
{
    for (int i = 0; i < num_cells; i++)
    {
	if (blankCell(i))
	{
	    double cell_volume = c_centers[i][2]*grid_spacing[0]*grid_spacing[1]*grid_spacing[2];

	    unsigned char cell_sign = face_sign[i][0];
	    int cell_comp = face_type[i][0];

	    if (cell_sign != 0)
	    {
		map<int,double>::iterator itm = vol_comm_pos.find(cell_comp);

		if (itm != vol_comm_pos.end())
		    itm->second += cell_volume;
		else
		{
		    printf("Cell %d:\n", i);
		    printf("ERROR!!!\n");
		}
	    }
	    else
	    {
		map<int,double>::iterator itm = vol_comm_neg.find(cell_comp);
		if(itm != vol_comm_neg.end())
		    itm->second += cell_volume;
		else
		{
		    printf("Cell %d:\n", i);
		    printf("ERROR!!!\n");
		}
	    }

	}
    }
}

void DropletAnalysis::adjustVolumeFromMapping(void)
{
    int from_comp, to_comp;
    int from_index, to_index;

    for (map<int,int>::iterator itpos = map_comm_pos.begin(); itpos != map_comm_pos.end(); itpos++)
    {
	from_comp = itpos->first;
	to_comp = itpos->second;

	map<int,double>::iterator itvol_from = vol_comm_pos.find(from_comp);
	map<int,double>::iterator itvol_to = vol_comm_pos.find(to_comp);

	if (itvol_from != vol_comm_pos.end())
	{
	    if (itvol_to != vol_comm_pos.end())
	    {
		itvol_to->second += itvol_from->second;
		itvol_from->second = 0.0;
	    }
	    else
	    {
		printf("\npos volume mapped to not included in this processor\n");
		vol_comm_pos.insert(pair<int,double>(to_comp, itvol_from->second));
		itvol_from->second = 0.0;
	    }
	}

    }

    for (map<int,int>::iterator itneg = map_comm_neg.begin(); itneg != map_comm_neg.end(); itneg++)
    {
	from_comp = itneg->first;
	to_comp = itneg->second;

	map<int,double>::iterator itvol_from = vol_comm_neg.find(from_comp);
	map<int,double>::iterator itvol_to = vol_comm_neg.find(to_comp);

	if (itvol_from != vol_comm_neg.end())
	{
	    if (itvol_to != vol_comm_neg.end())
	    {
		itvol_to->second += itvol_from->second;
		itvol_from->second = 0.0;
	    }
	    else
	    {
		printf("\nneg volume mapped to not included in this processor\n");
		vol_comm_neg.insert(pair<int,double>(to_comp, itvol_from->second));
		itvol_from->second = 0.0;
	    }
	}
	itvol_to->second += itvol_from->second;
	itvol_from->second = 0.0;
    }
}


void DropletAnalysis::addPairToMap(int big, int small, map<int,int> &map_for_modify)
{
    map<int,int>::iterator itm = map_for_modify.find(big);
    if(itm == map_for_modify.end())
        map_for_modify.insert(pair<int,int> (big, small));
    else
    {
        if (small < itm->second)
            addPairToMap(itm->second, small, map_for_modify);
        else if (small > itm->second)
            addPairToMap(small, itm->second, map_for_modify);
    }

}
void DropletAnalysis::addPairToPosMap(int big, int small)
{
    map<int,int>::iterator itm = map_comm_pos.find(big);
    if(itm == map_comm_pos.end())
        map_comm_pos.insert(pair<int,int> (big, small));
    else
    {
        if (small < itm->second)
            addPairToPosMap(itm->second, small);
        else if (small > itm->second)
            addPairToPosMap(small, itm->second);
    }
}

void DropletAnalysis::addPairToNegMap(int big, int small)
{
    map<int,int>::iterator itm = map_comm_neg.find(big);
    if(itm == map_comm_neg.end())
        map_comm_neg.insert(pair<int,int> (big, small));
    else
    {
        if (small < itm->second)
            addPairToNegMap(itm->second, small);
        else if (small > itm->second)
            addPairToNegMap(small, itm->second);
    }
}

void DropletAnalysis::initMapping(void)
{
    map_comm_pos.clear();
    map_comm_neg.clear();

    vol_comm_pos.clear();
    vol_comm_neg.clear();

    for (int i = 0; i < num_comps; i++)
    {
	int orig_index = comp_list[i];
	int pos_index = pos_map[i];
	int neg_index = neg_map[i];

	if (orig_index != pos_index)
	    map_comm_pos.insert(pair<int,int> (orig_index, pos_index));
	if (orig_index != neg_index)
	    map_comm_neg.insert(pair<int,int> (orig_index, neg_index));

	vol_comm_pos.insert(pair<int,double> (orig_index, pos_vol[i]));
	vol_comm_neg.insert(pair<int,double> (orig_index, neg_vol[i]));
    }
}

void DropletAnalysis::adjustFaceIndexFromMapping(void)
{
    for (int i = 0; i < num_cells; i++)
    {
	for (int j = 0; j < 6; j++)
	{
	    if (face_sign[i][j] != 0)
	    {
		map<int,int>::iterator itpos = map_comm_pos.find(face_type[i][j]);
		if (itpos != map_comm_pos.end())
		    face_type[i][j] = itpos->second;
	    }
	    else
	    {
		map<int,int>::iterator itneg = map_comm_neg.find(face_type[i][j]);
		if (itneg != map_comm_neg.end())
		    face_type[i][j] = itneg->second;
	    }
	}
    }
}

void DropletAnalysis::FillVolumeForBlankCells(void)
{
    initMapping();

/*
    double vol_all = 0.0;
    for (int i = 0; i < num_cells; i++)
    {
	if ( (!blankCell(i)) )
	{
	    double vol_cell = c_centers[i][2]*grid_spacing[0]*grid_spacing[1]*grid_spacing[2];
	    vol_all = vol_all + vol_cell;
	}
    }

    printf("\nTotal non blank volume is %20.16g\n", vol_all); 
 

    double total_vol_pos = 0.0;
    double total_vol_neg = 0.0;

    printf("\nPositive volume:\n");
    for (int i = 0; i < num_comps; i++)
    {
	printf("Comp Index = %d, volume = %20.16g\n", comp_list[i], pos_vol[i]);
	total_vol_pos = total_vol_pos + pos_vol[i];

    }

    printf("\nNegative volume:\n");
    for (int i = 0; i < num_comps; i++)
    {
	printf("Comp Index = %d, volume = %20.16g\n", comp_list[i], neg_vol[i]);
	total_vol_neg = total_vol_neg + neg_vol[i];

    }
    double total_vol_compare = total_vol_pos + total_vol_neg;

    printf("\nTotal volume from computation is %20.16g\n", total_vol_compare);
*/

    if ( allBlankCells() )
	blank_proc = true;
    else
	blank_proc = false;

    if (blank_proc == false)
	printf("\nFirst non blank cell = %d\n", first_non_blank_cell);
    if (blank_proc == false)
	getLocalCompIndexForCells();

    printf("\nFinished Step 3.1\n");


    getGlobalCompIndexForCells();
    printf("\nFinished Step 3.2\n");
    getFinalVolume();

    
    //

/*
    double total_vol_pos2 = 0.0;
    double total_vol_neg2 = 0.0;

    printf("\nPositive volume:\n");

    for (map<int,double>::iterator itp = vol_comm_pos.begin(); itp != vol_comm_pos.end(); itp++)
    {
	printf("Comp Index = %d, volume = %20.16g\n", itp->first, itp->second);
	total_vol_pos2 = total_vol_pos2 + itp->second;
    }
    printf("\nNegative volume:\n");
    for (map<int,double>::iterator itn = vol_comm_neg.begin(); itn != vol_comm_neg.end(); itn++)
    {
	printf("Comp Index = %d, volume = %20.16g\n", itn->first, itn->second);
	total_vol_neg2 = total_vol_neg2 + itn->second;
    }
    double total_vol_compare2 = total_vol_pos2 + total_vol_neg2;

    printf("\nTotal volume from computation is %20.16g\n", total_vol_compare2);
*/
}


void DropletAnalysis::getLocalCompIndexForCells(void)
{
    adjustFaceIndexFromMapping();

    int i;

    int **cell_nbs;
    bool *cell_visited;

    cell_nbs = new int* [num_cells];
    cell_visited = new bool [num_cells];

    for (i = 0; i < num_cells; i++)
    {
	int *cell_nb_for_insert = new int [6];

	for (int k = 0; k < 6; k++)
	{
	    if (face_type[i][k] != -1)
		cell_nb_for_insert[k] = findNbCellIndex(i, k);
	    else
		cell_nb_for_insert[k] = -1;
	}
	cell_nbs[i] = cell_nb_for_insert;
	cell_visited[i] = false;
    }

    for (i = first_non_blank_cell; i < num_cells; i++)//DFS visit
    {
	if ( cell_visited[i] == false )
	    visitCell(-1, i, 0, cell_nbs, cell_visited);
    }
    

    for (i = 0; i < num_cells; i++)
	delete[] cell_nbs[i];
    delete[] cell_nbs;
    delete[] cell_visited;

    adjustMapping(map_comm_pos);
    adjustMapping(map_comm_neg);

    adjustFaceIndexFromMapping();

    adjustVolumeFromMapping();
    adjustVolumeForBlankCells();

}

void DropletAnalysis::visitCell(int from_cell, int cur_cell, int from_dir, int **cell_nb_list, bool *visit_flag)
{
    visit_flag[cur_cell] = true;

    //printf("\n From cell %d, cur cell %d \n", from_cell, cur_cell);

    if (from_cell != -1)
    {
	unsigned char from_sign, to_sign;
	int from_type, to_type;
	int big_type, small_type;


	switch (from_dir)
	{
	    case 0:
		from_sign = face_sign[from_cell][0];
		to_sign = face_sign[cur_cell][5];
		from_type = face_type[from_cell][0];
		to_type = face_type[cur_cell][5];
		break;
	    case 1:
		from_sign = face_sign[from_cell][1];
		to_sign = face_sign[cur_cell][3];
		from_type = face_type[from_cell][1];
		to_type = face_type[cur_cell][3];
		break;
	    case 2:
		from_sign = face_sign[from_cell][2];
		to_sign = face_sign[cur_cell][4];
		from_type = face_type[from_cell][2];
		to_type = face_type[cur_cell][4];
		break;
	    case 3:
		from_sign = face_sign[from_cell][3];
		to_sign = face_sign[cur_cell][1];
		from_type = face_type[from_cell][3];
		to_type = face_type[cur_cell][1];
		break;
	    case 4:
		from_sign = face_sign[from_cell][4];
		to_sign = face_sign[cur_cell][2];
		from_type = face_type[from_cell][4];
		to_type = face_type[cur_cell][2];
		break;
	    case 5:
		from_sign = face_sign[from_cell][5];
		to_sign = face_sign[cur_cell][0];
		from_type = face_type[from_cell][5];
		to_type = face_type[cur_cell][0];
		break;
	    default:
		printf("\n Wrong direction !\n");
		clean_up(ERROR);
	}


	if (to_type == -2) //if blank cell, fill the cell with previous face type and sign
	{
	    for (int k = 0; k < 6; k++)
	    {
		face_type[cur_cell][k] = from_type;
		face_sign[cur_cell][k] = from_sign;
	    }
	}
	else //check whether we need a mapping
	{
	    if (from_sign != to_sign)
	    {
		printf("\n Face sign not matching between two nb cells! \n");
		printf("\n Error! from cell = %d, current cell = %d, from dir = %d \n", from_cell, cur_cell, from_dir);
		clean_up(ERROR);
	    }

	    if (from_type != to_type) 
	    {
		big_type = std::max(from_type, to_type);
		small_type = std::min(from_type, to_type);
    		if (from_sign != 0)
	    	    addPairToPosMap(big_type, small_type);
		else
		    addPairToNegMap(big_type, small_type);
	    }
	}
    }

    for (int i = 0; i < 6; i++)
    {
	int next_cell = cell_nb_list[cur_cell][i];

	if ( (next_cell != -1) && (visit_flag[next_cell] == false) )
	    visitCell(cur_cell, next_cell, i, cell_nb_list, visit_flag);
    }	

}

int DropletAnalysis::findNbCellIndex(int current_id, int dir)
{
    int icoords[3];
    int icoords_nb[3];
    int id = current_id;

    for (int d = 0; d < 3; d++)
    {
	int G = gmax_cell[d];
	icoords[d] = id % G;
	id = (id - icoords[d])/G;
    }

    switch (dir)
    {
	case 0://z-
	    if (icoords[2] == 0)
		return -1;
	    else
		return getCellIndex(icoords[0], icoords[1], icoords[2]-1);
	    break;
	case 1://x+
	    if (icoords[0] == (gmax_cell[0] -1))
		return -1;
	    else
		return getCellIndex(icoords[0]+1, icoords[1], icoords[2]);
	    break;
	case 2://y+
	    if (icoords[1] == (gmax_cell[1] -1))
		return -1;
	    else
		return getCellIndex(icoords[0], icoords[1]+1, icoords[2]);
	    break;
	case 3://x-
	    if (icoords[0] == 0)
		return -1;
	    else
		return getCellIndex(icoords[0]-1, icoords[1], icoords[2]);
	    break;
	case 4://y-
	    if (icoords[1] == 0)
		return -1;
	    else
		return getCellIndex(icoords[0], icoords[1]-1, icoords[2]);
	    break;
	case 5://z+
	    if (icoords[2] == (gmax_cell[2] - 1))
		return -1;
	    else
		return getCellIndex(icoords[0], icoords[1], icoords[2]+1);
	    break;
	default:
	    printf("\n Wrong dir in findNbCellIndex! \n");
	    clean_up(ERROR);
	    break;

    }
}


void DropletAnalysis::getFinalVolume(void)
{
    if (blank_proc == false)
	adjustVolumeFromMapping();

    if (blank_proc == true)
	getBlankProcVolume();
    mergeVolumeInfo(vol_comm_pos);
    mergeVolumeInfo(vol_comm_neg);
    
}

void DropletAnalysis::getGlobalCompIndexForCells(void)
{
    fillBlankProcFaces();

    commCellFace(); //communiate cell faces (both sign and index)
    //change the mapping_pos and mapping_neg again
    commMapping(map_comm_pos);
    commMapping(map_comm_neg);

    adjustMapping(map_comm_pos);
    adjustMapping(map_comm_neg);

    adjustFaceIndexFromMapping();
}

void DropletAnalysis::commCellFace(void)
{
    int i,j,k,l;
    int bdry_type[MAXD][2];
    MPI_Request request;
    MPI_Status status;

    int dim = 3;

    int *G = pSolver->pp_grid->gmax;

    int me[MAXD], him[MAXD];

    int my_id, dst_id, rcv_id;

    INTERFACE *intfc = pSolver->tempfront->interf;

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
	for(j = 1; j >= 0; j--)
	{
	    int *index_send = 0;
	    unsigned char *sign_send = 0;

	    for (k = 0; k < dim; k++)
		him[k] = me[k];

	    pp_gsync();
	    if (bdry_type[i][j] == SUBDOMAIN_BOUNDARY)
	    {
		makeCellFaceForSend(i, j, index_send, sign_send);
		him[i] = me[i] + 2*j - 1;
		him[i] = (him[i] + G[i]) % G[i];
		dst_id = domain_id(him,G,dim);
		sendCellFaceTo(i, j, dst_id, index_send, sign_send);
	    }

	    if (bdry_type[i][((j+1)%2)] == SUBDOMAIN_BOUNDARY)
	    {
		him[i] = me[i] - 2*j + 1;
		him[i] = (him[i] + G[i]) % G[i];
		rcv_id = domain_id(him,G,dim);
		mergeCellFaceFrom(i, ((j+1)%2), rcv_id);
	    }
	    pp_gsync();

	    delete[] index_send;
	    delete[] sign_send;
	}
    }
}


bool DropletAnalysis::allBlankCells(void)
{
    for (int i = 0; i < num_cells; i++)
    {
	if ( !(blankCell(i)) )
	{
	    first_non_blank_cell = i;
	    return false;
	}
    }
    return true;
}

bool DropletAnalysis::blankCell(int cell_index)
{
    for (int i = 0; i < 5; i++)
    {
	if (face_type[cell_index][i] != face_type[cell_index][i+1])
	    return false;
    }

    if (face_type[cell_index][0] == -1)
	return false;
    else
	return true;
}

void DropletAnalysis::getLocalCompIndexForTris(void)
{

    TriGraph tri_connect(this);
    tri_connect.Init();
    tri_connect.GetTriComp();

    Graph_Point* g_tri = tri_connect.getVertexPointer();

    for (int i = 0; i < num_tris; i++)
	comp_index[i] = (g_tri[i]).set_index;
}



void DropletAnalysis::commCompIndexForTris(void)
{
    int i,j,k,l;
    int bdry_type[MAXD][2];
    MPI_Request request;
    MPI_Status status;

    int dim = 3;

    int *G = pSolver->pp_grid->gmax;

    int me[MAXD], him[MAXD];

    int my_id, dst_id, rcv_id;

    INTERFACE *intfc = pSolver->tempfront->interf;

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
	for (l = 0; l < G[i] - 1; l++)
	{
	    for(j = 1; j >= 0; j--)
	    {
		pp_gsync();
	    
		for (k = 0; k < dim; k++)
		    him[k] = me[k];

		if (bdry_type[i][j] == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = me[i] + 2*j - 1;
		    him[i] = (him[i] + G[i]) % G[i];
		    dst_id = domain_id(him,G,dim);
		    sendTriCompTo(i, j, dst_id, &request);
		}

		if (bdry_type[i][((j+1)%2)] == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = me[i] - 2*j + 1;
		    him[i] = (him[i] + G[i]) % G[i];
		    rcv_id = domain_id(him,G,dim);
		    mergeTriCompFrom(i, ((j+1)%2), rcv_id);
		}
		pp_gsync();
	    }
	}
    }
}


void DropletAnalysis::makeCellFaceForSend(int dim, int dir, int *&index_for_send, unsigned char *&sign_for_send)
{
    int i,j,k,index,l;
    int num_faces;

    switch(dim)
    {
	case 0:
	    num_faces = gmax_cell[1] * gmax_cell[2];
	    index_for_send = new int [num_faces];
	    sign_for_send = new unsigned char [num_faces];
	    switch(dir)
	    {
		case 0:
		    i = 0;
		    l = 0;
		    for (k = 0; k < gmax_cell[2]; k++)
		    for (j = 0; j < gmax_cell[1]; j++)
		    {
			index = getCellIndex(i,j,k);
			index_for_send[l] = face_type[index][3];
			sign_for_send[l] = face_sign[index][3];
			l++;
		    }
		    break;
		case 1:
		    i = gmax_cell[0] - 1;
		    l = 0;
		    for (k = 0; k < gmax_cell[2]; k++)
		    for (j = 0; j < gmax_cell[1]; j++)
		    {
			index = getCellIndex(i,j,k);
			index_for_send[l] = face_type[index][1];
			sign_for_send[l] = face_sign[index][1];
			l++;
		    }
		    break;
		default:
		    printf("\nWrong direction!\n");
		    clean_up(ERROR);
		    break;
	    }
	    break;
	case 1:
	    num_faces = gmax_cell[0] * gmax_cell[2];
	    index_for_send = new int [num_faces];
	    sign_for_send = new unsigned char [num_faces];
	    switch(dir)
	    {
		case 0:
		    j = 0;
		    l = 0;
		    for (k = 0; k < gmax_cell[2]; k++)
		    for (i = 0; i < gmax_cell[0]; i++)
		    {
			index = getCellIndex(i,j,k);
			index_for_send[l] = face_type[index][4];
			sign_for_send[l] = face_sign[index][4];
			l++;
		    }
		    break;
		case 1:
		    j = gmax_cell[1] - 1;
		    l = 0;
		    for (k = 0; k < gmax_cell[2]; k++)
		    for (i = 0; i < gmax_cell[0]; i++)
		    {
			index = getCellIndex(i,j,k);
			index_for_send[l] = face_type[index][2];
			sign_for_send[l] = face_sign[index][2];
			l++;
		    }
		    break;
		default:
		    printf("\nWrong direction!\n");
		    clean_up(ERROR);
		    break;
	    }
	    break;
	case 2:
	    num_faces = gmax_cell[0] * gmax_cell[1];
	    index_for_send = new int [num_faces];
	    sign_for_send = new unsigned char [num_faces];
	    switch(dir)
	    {
		case 0:
		    k = 0;
		    l = 0;
		    for (j = 0; j < gmax_cell[1]; j++)
		    for (i = 0; i < gmax_cell[0]; i++)
		    {
			index = getCellIndex(i,j,k);
			index_for_send[l] = face_type[index][0];
			sign_for_send[l] = face_sign[index][0];
			l++;
		    }
		    break;
		case 1:
		    k = gmax_cell[2] - 1;
		    l = 0;
		    for (j = 0; j < gmax_cell[1]; j++)
		    for (i = 0; i < gmax_cell[0]; i++)
		    {
			index = getCellIndex(i,j,k);
			index_for_send[l] = face_type[index][5];
			sign_for_send[l] = face_sign[index][5];
			l++;
		    }
		    break;
		default:
		    printf("\nWrong direction!\n");
		    clean_up(ERROR);
		    break;
	    }
	    break;
	default:
	    printf("\nWrong dimension!\n");
	    clean_up(ERROR);
	    break;
    }
    
}
void DropletAnalysis::sendCellFaceTo(int dim, int dir, int dst_id, int *&index_for_send, unsigned char *&sign_for_send)
{
    int tag = 1;

    MPI_Request index_request;
    MPI_Request sign_request;

    int num_faces;

    switch(dim)
    {
	case 0:
	    num_faces = gmax_cell[1]*gmax_cell[2];
	    break;
	case 1:
	    num_faces = gmax_cell[0]*gmax_cell[2];
	    break;
	case 2:
	    num_faces = gmax_cell[0]*gmax_cell[1];
	    break;
	default:
	    printf("\nWrong dimension\n");
	    clean_up(ERROR);
	    break;
    }


    MPI_Isend(index_for_send, num_faces, MPI_INT, dst_id, tag, MPI_COMM_WORLD, &index_request);
    MPI_Isend(sign_for_send, num_faces, MPI_UNSIGNED_CHAR, dst_id, tag, MPI_COMM_WORLD, &sign_request);
}

void DropletAnalysis::mergeCellFaceFrom(int dim, int dir, int rec_id)
{
    int tag = 1;
    int num_faces;

    int *index_recv;
    unsigned char *sign_recv;

    MPI_Status stat;

    int i,j,k,l,index;
    switch(dim)
    {
	case 0: //In theta-direction
	    num_faces = gmax_cell[1] * gmax_cell[2];

	    index_recv = new int[num_faces];
	    sign_recv = new unsigned char[num_faces];

	    MPI_Recv(index_recv, num_faces, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
	    MPI_Recv(sign_recv, num_faces, MPI_UNSIGNED_CHAR, rec_id, tag, MPI_COMM_WORLD, &stat);

	    switch(dir)
	    {
		case 0:
		    i = 0;
		    l = 0;
		    for (k = 0; k < gmax_cell[2]; k++)
		    for (j = 0; j < gmax_cell[1]; j++)
		    {
			index = getCellIndex(i,j,k);
			if (sign_recv[l] == face_sign[index][3])
			{
			    if (index_recv[l] > face_type[index][3])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(index_recv[l], face_type[index][3]);
				else
				    addPairToNegMap(index_recv[l], face_type[index][3]);
			    }
			    else if (index_recv[l] < face_type[index][3])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(face_type[index][3], index_recv[l]);
				else
				    addPairToNegMap(face_type[index][3], index_recv[l]);
			    }
			}
			else
			{
			    printf("\nNot matching sign between two processors!\n");
			    printf("\nDim = %d, Dir = %d, Current proc = %d, Rec id = %d\n",dim, dir, pp_mynode(), rec_id);
			    printf("\nIndex = %d, i = %d, j = %d, k = %d, face_sign = %d, face_type = %d\n", index, i, j, k, face_sign[index][3],face_type[index][3]);
			    printf("\nFrom face_sign = %d, face_type = %d\n", sign_recv[l], index_recv[l]);
			    clean_up(ERROR);
			}
			l++;
		    }
		    break;
		case 1:
		    i = gmax_cell[0]-1;
		    l = 0;
		    for (k = 0; k < gmax_cell[2]; k++)
		    for (j = 0; j < gmax_cell[1]; j++)
		    {
			index = getCellIndex(i,j,k);
			if (sign_recv[l] == face_sign[index][1])
			{
			    if (index_recv[l] > face_type[index][1])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(index_recv[l], face_type[index][1]);
				else
				    addPairToNegMap(index_recv[l], face_type[index][1]);
			    }
			    else if (index_recv[l] < face_type[index][1])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(face_type[index][1], index_recv[l]);
				else
				    addPairToNegMap(face_type[index][1], index_recv[l]);
			    }
			}
			else
			{	
			    printf("\nNot matching sign between two processors!\n");
			    printf("\nDim = %d, Dir = %d, Current proc = %d, Rec id = %d\n", dim, dir, pp_mynode(), rec_id);
			    printf("\nIndex = %d, i = %d, j = %d, k = %d, face_sign = %d, face_type = %d\n", index, i, j, k, face_sign[index][1],face_type[index][1]);
			    printf("\nFrom face_sign = %d, face_type = %d\n", sign_recv[l], index_recv[l]);
			    clean_up(ERROR);
			}
			l++;
		    }
		    break;
	    }
	    delete[] index_recv;
	    delete[] sign_recv;
	    index_recv = 0;
	    sign_recv = 0;
	    break;
	case 1: 
	    num_faces = gmax_cell[0] * gmax_cell[2];

	    index_recv = new int[num_faces];
	    sign_recv = new unsigned char[num_faces];

	    MPI_Recv(index_recv, num_faces, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
	    MPI_Recv(sign_recv, num_faces, MPI_UNSIGNED_CHAR, rec_id, tag, MPI_COMM_WORLD, &stat);

	    switch(dir)
	    {
		case 0:
		    j = 0;
		    l = 0;
		    for (k = 0; k < gmax_cell[2]; k++)
		    for (i = 0; i < gmax_cell[0]; i++)
		    {
			index = getCellIndex(i,j,k);
			if (sign_recv[l] == face_sign[index][4])
			{
			    if (index_recv[l] > face_type[index][4])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(index_recv[l], face_type[index][4]);
				else
				    addPairToNegMap(index_recv[l], face_type[index][4]);
			    }
			    else if (index_recv[l] < face_type[index][4])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(face_type[index][4], index_recv[l]);
				else
				    addPairToNegMap(face_type[index][4], index_recv[l]);
			    }
			}
			else
			{	
			    printf("\nNot matching sign between two processors!\n");
			    printf("\nDim = %d, Dir = %d, Current proc = %d, Rec id = %d\n", dim, dir, pp_mynode(), rec_id);
			    printf("\nIndex = %d, i = %d, j = %d, k = %d, face_sign = %d, face_type = %d\n", index, i, j, k, face_sign[index][4],face_type[index][4]);
			    printf("\nFrom face_sign = %d, face_type = %d\n", sign_recv[l], index_recv[l]);
			    clean_up(ERROR);
			}
			l++;
		    }
		    break;
		case 1:
		    j = gmax_cell[1]-1;
		    l = 0;
		    for (k = 0; k < gmax_cell[2]; k++)
		    for (i = 0; i < gmax_cell[0]; i++)
		    {
			index = getCellIndex(i,j,k);
			if (sign_recv[l] == face_sign[index][2])
			{
			    if (index_recv[l] > face_type[index][2])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(index_recv[l], face_type[index][2]);
				else
				    addPairToNegMap(index_recv[l], face_type[index][2]);
			    }
			    else if (index_recv[l] < face_type[index][2])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(face_type[index][2], index_recv[l]);
				else
				    addPairToNegMap(face_type[index][2], index_recv[l]);
			    }
			}
			else
			{	
			    printf("\nNot matching sign between two processors!\n");
			    printf("\nDim = %d, Dir = %d, Current proc = %d, Rec id = %d\n", dim, dir, pp_mynode(), rec_id);
			    printf("\nIndex = %d, i = %d, j = %d, k = %d, face_sign = %d, face_type = %d\n", index, i, j, k, face_sign[index][2],face_type[index][2]);
			    printf("\nFrom face_sign = %d, face_type = %d\n", sign_recv[l], index_recv[l]);
			    clean_up(ERROR);
			}
			l++;
		    }
		    break;
	    }
	    delete[] index_recv;
	    delete[] sign_recv;
	    index_recv = 0;
	    sign_recv = 0;
	    break;
	case 2: //In r-direction
	    num_faces = gmax_cell[0] * gmax_cell[1];

	    index_recv = new int[num_faces];
	    sign_recv = new unsigned char[num_faces];

	    MPI_Recv(index_recv, num_faces, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
	    MPI_Recv(sign_recv, num_faces, MPI_UNSIGNED_CHAR, rec_id, tag, MPI_COMM_WORLD, &stat);

	    switch(dir)
	    {
		case 0:
		    k = 0;
		    l = 0;
		    for (j = 0; j < gmax_cell[1]; j++)
		    for (i = 0; i < gmax_cell[0]; i++)
		    {
			index = getCellIndex(i,j,k);
			if (sign_recv[l] == face_sign[index][0])
			{
			    if (index_recv[l] > face_type[index][0])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(index_recv[l], face_type[index][0]);
				else
				    addPairToNegMap(index_recv[l], face_type[index][0]);
			    }
			    else if (index_recv[l] < face_type[index][0])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(face_type[index][0], index_recv[l]);
				else
				    addPairToNegMap(face_type[index][0], index_recv[l]);
			    }
			}
			else
			{	
			    printf("\nNot matching sign between two processors!\n");
			    printf("\nDim = %d, Dir = %d, Current proc = %d, Rec id = %d\n", dim, dir, pp_mynode(), rec_id);
			    printf("\nIndex = %d, i = %d, j = %d, k = %d, face_sign = %d, face_type = %d\n", index, i, j, k, face_sign[index][0],face_type[index][0]);
			    printf("\nFrom face_sign = %d, face_type = %d\n", sign_recv[l], index_recv[l]);
			    clean_up(ERROR);
			}
			l++;
		    }
		    break;
		case 1:
		    k = gmax_cell[2]-1;
		    l = 0;
		    for (j = 0; j < gmax_cell[1]; j++)
		    for (i = 0; i < gmax_cell[0]; i++)
		    {
			index = getCellIndex(i,j,k);
			if (sign_recv[l] == face_sign[index][5])
			{
			    if (index_recv[l] > face_type[index][5])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(index_recv[l], face_type[index][5]);
				else
				    addPairToNegMap(index_recv[l], face_type[index][5]);
			    }
			    else if (index_recv[l] < face_type[index][5])
			    {
				if(sign_recv[l] != 0)
				    addPairToPosMap(face_type[index][5], index_recv[l]);
				else
				    addPairToNegMap(face_type[index][5], index_recv[l]);
			    }
			}
			else
			{	
			    printf("\nNot matching sign between two processors!\n");
			    printf("\nDim = %d, Dir = %d, Current proc = %d, Rec id = %d\n", dim, dir, pp_mynode(), rec_id);
			    printf("\nIndex = %d, i = %d, j = %d, k = %d, face_sign = %d, face_type = %d\n", index, i, j, k, face_sign[index][5],face_type[index][5]);
			    printf("\nFrom face_sign = %d, face_type = %d\n", sign_recv[l], index_recv[l]);
			    clean_up(ERROR);
			}
			l++;
		    }
		    break;
	    }
	    delete[] index_recv;
	    delete[] sign_recv;
	    index_recv = 0;
	    sign_recv = 0;
	    break;
	default:
	    printf("\n Wrong dimension! \n");
	    clean_up(ERROR);
	    break;
    }

}

void DropletAnalysis::sendTriCompTo(int dim, int dir, int dst_id, MPI_Request* request)
{
    int tag = 1;
    MPI_Isend(&num_tris, 1, MPI_INT, dst_id, tag, MPI_COMM_WORLD, request);
    MPI_Isend(global_index, num_tris, MPI_INT, dst_id, tag, MPI_COMM_WORLD, request);
    MPI_Isend(comp_index, num_tris, MPI_INT, dst_id, tag, MPI_COMM_WORLD, request);
}

void DropletAnalysis::mergeTriCompFrom(int dim, int dir, int rec_id)
{
    int tag = 1;
    MPI_Status stat;
    int *num_tris_rec = new int;
    int *global_index_recv;
    int *comp_index_recv;

    MPI_Recv(num_tris_rec, 1, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    global_index_recv = new int [(*num_tris_rec)];
    comp_index_recv = new int [(*num_tris_rec)];

    MPI_Recv(global_index_recv, (*num_tris_rec), MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(comp_index_recv, (*num_tris_rec), MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

    for (int i = 0; i < (*num_tris_rec); i++)
    {
	for (int j = 0; j < num_tris; j++)
	{
	    if(global_index[j] == global_index_recv[i])
	    {
		if(comp_index[j] > comp_index_recv[i])
		{
		    int origin_comp_index = comp_index[j];
		    for (int k = 0; k < num_tris; k++)
			if(comp_index[k] == origin_comp_index)
			    comp_index[k] = comp_index_recv[i];
		}
	    }
	}
    }
   
    delete[] global_index_recv;
    delete[] comp_index_recv;

    delete num_tris_rec;
}

void DropletAnalysis::commSet(set<int> &set_for_comm)
{
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
		sendSetTo(dst_id, 0, set_for_comm);

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
		    mergeSetFrom(rcv_id, 0, set_for_comm);
		    sendSetTo(dst_id, 0, set_for_comm);
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
		mergeSetFrom(rcv_id, 0, set_for_comm);
	}
    }

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
		sendSetTo(dst_id, 1, set_for_comm);

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
		    mergeSetFrom(rcv_id, 1, set_for_comm);
		    sendSetTo(dst_id, 1, set_for_comm);
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
		mergeSetFrom(rcv_id, 1, set_for_comm);
	}
    }

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
	    sendSetTo(dst_id, 2, set_for_comm);

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
		mergeSetFrom(rcv_id, 2, set_for_comm);
		sendSetTo(dst_id, 2, set_for_comm);
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
	    mergeSetFrom(rcv_id, 2, set_for_comm);
    }

    pp_gsync();

    int num_set_bcast = (int) set_for_comm.size();

    MPI_Bcast(&num_set_bcast, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *set_index_bcast = new int [num_set_bcast];

    int p = 0;
    if (pp_mynode() == 0)
    {
	for (set<int>::iterator it = set_for_comm.begin(); it != set_for_comm.end(); it++)
	{
	    set_index_bcast[p] = (*it);
	    p++;
	}
    }

    MPI_Bcast(set_index_bcast, num_set_bcast, MPI_INT, 0, MPI_COMM_WORLD);

    set_for_comm.clear();

    for (p = 0; p < num_set_bcast; p++)
	set_for_comm.insert(set_index_bcast[p]);

    delete[] set_index_bcast;
}


void DropletAnalysis::commMapping(map<int,int> &map_for_comm)
{
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
		sendMappingTo(dst_id, 0, map_for_comm);

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
		    mergeMappingFrom(rcv_id, 0, map_for_comm);
		    sendMappingTo(dst_id, 0, map_for_comm);
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
		mergeMappingFrom(rcv_id, 0, map_for_comm);
	}
    }

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
		sendMappingTo(dst_id, 1, map_for_comm);

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
		    mergeMappingFrom(rcv_id, 1, map_for_comm);
		    sendMappingTo(dst_id, 1, map_for_comm);
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
		mergeMappingFrom(rcv_id, 1, map_for_comm);
	}
    }

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
	    sendMappingTo(dst_id, 2, map_for_comm);

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
		mergeMappingFrom(rcv_id, 2, map_for_comm);
		sendMappingTo(dst_id, 2, map_for_comm);
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
	    mergeMappingFrom(rcv_id, 2, map_for_comm);
    }

    pp_gsync();

    int num_map_bcast = (int) map_for_comm.size();

    MPI_Bcast(&num_map_bcast, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *orig_index_bcast = new int [num_map_bcast];
    int *map_index_bcast = new int [num_map_bcast];

    int p = 0;
    if (pp_mynode() == 0)
    {
	for (map<int,int >::iterator it = map_for_comm.begin(); it != map_for_comm.end(); it++)
	{
	    orig_index_bcast[p] = it->first;
	    map_index_bcast[p] = it->second;
	    p++;
	}
    }

    MPI_Bcast(orig_index_bcast, num_map_bcast, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(map_index_bcast, num_map_bcast, MPI_INT, 0, MPI_COMM_WORLD);

    map_for_comm.clear();

    for (p = 0; p < num_map_bcast; p++)
	map_for_comm.insert(pair<int,int>(orig_index_bcast[p], map_index_bcast[p]));

    delete[] orig_index_bcast;
    delete[] map_index_bcast;
}

void DropletAnalysis::commCompIndexForTris2(void)
{
    map_comm.clear();
    int i,j,k,l;
    int bdry_type[MAXD][2];
    MPI_Request request;
    MPI_Status status;

    int dim = 3;

    int *G = pSolver->pp_grid->gmax;

    int me[MAXD], him[MAXD];

    int my_id, dst_id, rcv_id;

    INTERFACE *intfc = pSolver->tempfront->interf;

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
	    for(j = 0; j < 2; j++)
	    {
		int *tris_num_for_send = 0;

		int *global_index_for_send = 0;
		int *comp_index_for_send = 0;
	    
		for (k = 0; k < dim; k++)
		    him[k] = me[k];

		pp_gsync();
		if (bdry_type[i][j] == SUBDOMAIN_BOUNDARY)
		{
		    //makeTriCompForSend2(i, j, tris_num_for_send, global_index_for_send, comp_index_for_send);
		    him[i] = me[i] + 2*j - 1;
		    him[i] = (him[i] + G[i]) % G[i];
		    dst_id = domain_id(him,G,dim);
		    //sendTriCompTo2(i, j, dst_id, tris_num_for_send, global_index_for_send, comp_index_for_send);
		    sendTriCompTo2(i, j, dst_id, &request);
		}

		if (bdry_type[i][((j+1)%2)] == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = me[i] - 2*j + 1;
		    him[i] = (him[i] + G[i]) % G[i];
		    rcv_id = domain_id(him,G,dim);
		    mergeTriCompFrom2(i, ((j+1)%2), rcv_id);
		}
		pp_gsync();

		delete tris_num_for_send;
		delete[] global_index_for_send;
		delete[] comp_index_for_send;
	    }
    }
}


void DropletAnalysis::makeTriCompForSend2(int dim, int dir, int *&num_tris_send, int *&global_index_send, int *&comp_index_send)
{
    vector<int> global_index_bdry;
    vector<int> comp_index_bdry;

    global_index_bdry.clear(); comp_index_bdry.clear();

    int *G = pSolver->pp_grid->gmax;
    int i;

    double lBound, uBound;

    double *x_pp_bdry = pSolver->x_pp_bdry;
    double *y_pp_bdry = pSolver->y_pp_bdry;
    double *z_pp_bdry = pSolver->z_pp_bdry;

    double *top_h = pSolver->top_h;

    int x_pp_index = pSolver->x_pp_index;
    int y_pp_index = pSolver->y_pp_index;
    int z_pp_index = pSolver->z_pp_index;

    switch(dim)
    {
	case 0:
	    if (dir == 0)
	    {
		lBound = 3.0*(x_pp_bdry[x_pp_index] - 3.0*top_h[0]);
		uBound = 3.0*(x_pp_bdry[x_pp_index] + 6.0*top_h[0]);

	    }
	    else if (dir == 1)
	    {
		lBound = 3.0*(x_pp_bdry[x_pp_index+1] - 6.0*top_h[0]);
		uBound = 3.0*(x_pp_bdry[x_pp_index+1] + 3.0*top_h[0]);
	    }
	    else
	    {
		printf("\nWrong Direction!\n");
		clean_up(ERROR);
	    }
	    for (i = 0; i < num_tris; i++)
	    {
		    double center0;
		    center0 = ps[(tris[i][0])][0] + ps[(tris[i][1])][0] + ps[(tris[i][2])][0];
		    if( (center0 >= lBound) && (center0 <= uBound) )
		    {
			global_index_bdry.push_back(global_index[i]);
			comp_index_bdry.push_back(comp_index[i]);
		    }
	    }
	    break;
	case 1:
	    if (dir == 0)
	    {
		lBound = 3.0*(y_pp_bdry[y_pp_index] - 3.0*top_h[1]);
		uBound = 3.0*(y_pp_bdry[y_pp_index] + 6.0*top_h[1]);

	    }
	    else if (dir == 1)
	    {
		lBound = 3.0*(y_pp_bdry[y_pp_index+1] - 6.0*top_h[1]);
		uBound = 3.0*(y_pp_bdry[y_pp_index+1] + 3.0*top_h[1]);
	    }
	    else
	    {
		printf("\nWrong Direction!\n");
		clean_up(ERROR);
	    }
	    for (i = 0; i < num_tris; i++)
	    {
		    double center1;
		    center1 = ps[(tris[i][0])][1] + ps[(tris[i][1])][1] + ps[(tris[i][2])][1];
		    if( (center1 >= lBound) && (center1 <= uBound) )
		    {
			global_index_bdry.push_back(global_index[i]);
			comp_index_bdry.push_back(comp_index[i]);
		    }
	    }	    
	    break;
	case 2:
	    if (dir == 0)
	    {
		lBound = 3.0*(z_pp_bdry[z_pp_index] - 3.0*top_h[2]);
		uBound = 3.0*(z_pp_bdry[z_pp_index] + 6.0*top_h[2]);

	    }
	    else if (dir == 1)
	    {
		lBound = 3.0*(z_pp_bdry[z_pp_index+1] - 6.0*top_h[2]);
		uBound = 3.0*(z_pp_bdry[z_pp_index+1] + 3.0*top_h[2]);
	    }
	    else
	    {
		printf("\nWrong Direction!\n");
		clean_up(ERROR);
	    }
	    for (i = 0; i < num_tris; i++)
	    {
		    double center2;
		    center2 = ps[(tris[i][0])][2] + ps[(tris[i][1])][2] + ps[(tris[i][2])][2];
		    if( (center2 >= lBound) && (center2 <= uBound) )
		    {
			global_index_bdry.push_back(global_index[i]);
			comp_index_bdry.push_back(comp_index[i]);
		    }
	    }	 	   
	    break;
	default:
	    printf("\nWrong Dimension!\n");
	    clean_up(ERROR);
	    break;
    }

    int num_tris_bdry = (int) global_index_bdry.size();

    num_tris_send = new int (num_tris_bdry);

    global_index_send = new int [num_tris_bdry];
    comp_index_send = new int [num_tris_bdry];


    for (i = 0; i < num_tris_bdry; i++)
    {
	global_index_send[i] = global_index_bdry[i];
	comp_index_send[i] = comp_index_bdry[i];
    }

    global_index_bdry.clear();
    comp_index_bdry.clear();


}

/*
void DropletAnalysis::sendTriCompTo2(int dim, int dir, int dst_id, int *&num_tris_send, int *&global_index_send, int *&comp_index_send)
{
    int tag = 1;
    MPI_Request request_num;
    MPI_Request request_globali;
    MPI_Request request_compi;

    MPI_Isend(num_tris_send, 1, MPI_INT, dst_id, tag, MPI_COMM_WORLD, request_num);
    MPI_Isend(global_index_send, (*num_tris_send), MPI_INT, dst_id, tag, MPI_COMM_WORLD, request_globali);
    MPI_Isend(comp_index_send, (*num_tris_send), MPI_INT, dst_id, tag, MPI_COMM_WORLD, request_compi);




}
  */  
void DropletAnalysis::sendTriCompTo2(int dim, int dir, int dst_id, MPI_Request* request)
{
    int tag = 1;
    MPI_Isend(&num_tris, 1, MPI_INT, dst_id, tag, MPI_COMM_WORLD, request);
    MPI_Isend(global_index, num_tris, MPI_INT, dst_id, tag, MPI_COMM_WORLD, request);
    MPI_Isend(comp_index, num_tris, MPI_INT, dst_id, tag, MPI_COMM_WORLD, request);

}

void DropletAnalysis::mergeTriCompFrom2(int dim, int dir, int rec_id)
{
    int tag = 1;
    MPI_Status stat;
    int *num_tris_rec = new int;
    int *global_index_recv;
    int *comp_index_recv;

    MPI_Recv(num_tris_rec, 1, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    global_index_recv = new int [(*num_tris_rec)];
    comp_index_recv = new int [(*num_tris_rec)];

    MPI_Recv(global_index_recv, (*num_tris_rec), MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(comp_index_recv, (*num_tris_rec), MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

    for (int i = 0; i < (*num_tris_rec); i++)
    {
	for (int j = 0; j < num_tris; j++)
	{
	    if(global_index[j] == global_index_recv[i])
	    {
		if(comp_index[j] > comp_index_recv[i])
		    addPairToMap(comp_index[j], comp_index_recv[i], map_comm);
		else if(comp_index[j] < comp_index_recv[i])
		    addPairToMap(comp_index_recv[i], comp_index[j], map_comm);
	    }
	}
    }
   
    delete[] global_index_recv;
    delete[] comp_index_recv;

    delete num_tris_rec;
}

void DropletAnalysis::sendSetTo(int dst_id, int dim, set<int> &set_for_send)
{
    int tag = 1;
    int num_set = (int) set_for_send.size();
    printf("\nnum_set = %d\n", num_set);

    MPI_Send(&num_set, 1, MPI_INT, dst_id, tag, MPI_COMM_WORLD);

    if (num_set == 0)
	return;

    int *set_index;

    set_index = new int [num_set];

    int i = 0;

    set<int>::iterator its;

    for (its = set_for_send.begin(); its != set_for_send.end(); its++)
    {
	set_index[i] = (*its);
	i++;
    }

    MPI_Send(set_index, num_set, MPI_INT, dst_id, tag, MPI_COMM_WORLD);

    delete[] set_index;
}

void DropletAnalysis::mergeSetFrom(int rec_id, int dim, set<int> &set_for_recv)
{
    int tag = 1;
    int num_set_recv;
    MPI_Status stat;
    MPI_Recv(&num_set_recv, 1, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

    if (num_set_recv == 0)
	return;

    int *set_index_recv = new int [num_set_recv];


    MPI_Recv(set_index_recv, num_set_recv, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

    for (int i = 0; i < num_set_recv; i++)
    {
	set_for_recv.insert(set_index_recv[i]);
    }

    delete[] set_index_recv;
}

void DropletAnalysis::sendVolumeTo(int dst_id, int dim, map<int,double> &vol_for_send)
{
    int tag = 1;
    int num_map = (int) vol_for_send.size();

    MPI_Send(&num_map, 1, MPI_INT, dst_id, tag, MPI_COMM_WORLD);

    if (num_map == 0)
	return;

    int *phy_index;
    double *vol_value;

    phy_index = new int [num_map];
    vol_value = new double [num_map];

    int i = 0;

    map<int, double>::iterator itm;

    for (itm = vol_for_send.begin(); itm != vol_for_send.end(); itm++)
    {
	phy_index[i] = itm->first;
	vol_value[i] = itm->second;
	i++;
    }

    MPI_Send(phy_index, num_map, MPI_INT, dst_id, tag, MPI_COMM_WORLD);
    MPI_Send(vol_value, num_map, MPI_DOUBLE, dst_id, tag, MPI_COMM_WORLD);

    delete[] phy_index;
    delete[] vol_value;
}

void DropletAnalysis::mergeVolumeFrom(int rec_id, int dim, map<int,double> &vol_for_recv)
{
    int tag = 1;
    int num_map_recv;
    MPI_Status stat;
    MPI_Recv(&num_map_recv, 1, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

    if (num_map_recv == 0)
	return;

    int *phy_index_recv = new int [num_map_recv];
    double *vol_value_recv = new double [num_map_recv];


    MPI_Recv(phy_index_recv, num_map_recv, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(vol_value_recv, num_map_recv, MPI_DOUBLE, rec_id, tag, MPI_COMM_WORLD, &stat);

    for (int i = 0; i < num_map_recv; i++)
    {
	map<int,double>::iterator itm = vol_for_recv.find(phy_index_recv[i]);

	if (itm != vol_for_recv.end())
	    itm->second += vol_value_recv[i];
	else
	    vol_for_recv.insert(pair<int,double>(phy_index_recv[i], vol_value_recv[i]));
    }


    delete[] phy_index_recv;
    delete[] vol_value_recv;
}



void DropletAnalysis::sendMappingTo(int dst_id, int dim, map<int,int> &map_for_send)
{
    int tag = 1;
    int num_map = (int) map_for_send.size();

    MPI_Send(&num_map, 1, MPI_INT, dst_id, tag, MPI_COMM_WORLD);

    if (num_map == 0)
	return;

    int *orig_index;
    int *map_index;

    orig_index = new int [num_map];
    map_index = new int [num_map];

    int i = 0;

    map<int, int>::iterator itm;

    for (itm = map_for_send.begin(); itm != map_for_send.end(); itm++)
    {
	orig_index[i] = itm->first;
	map_index[i] = itm->second;
	i++;
    }

    MPI_Send(orig_index, num_map, MPI_INT, dst_id, tag, MPI_COMM_WORLD);
    MPI_Send(map_index, num_map, MPI_INT, dst_id, tag, MPI_COMM_WORLD);

    delete[] orig_index;
    delete[] map_index;
}

void DropletAnalysis::mergeMappingFrom(int rec_id, int dim, map<int,int> &map_for_recv)
{
    int tag = 1;
    int num_map_recv;
    MPI_Status stat;
    MPI_Recv(&num_map_recv, 1, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

    if (num_map_recv == 0)
	return;

    int *orig_index_recv = new int [num_map_recv];
    int *map_index_recv = new int [num_map_recv];


    MPI_Recv(orig_index_recv, num_map_recv, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);
    MPI_Recv(map_index_recv, num_map_recv, MPI_INT, rec_id, tag, MPI_COMM_WORLD, &stat);

    for (int i = 0; i < num_map_recv; i++)
	addPairToMap(orig_index_recv[i], map_index_recv[i], map_for_recv);

    delete[] orig_index_recv;
    delete[] map_index_recv;
}



void DropletAnalysis::getGlobalCompIndexForTris(void)
{
    start_clock("commCompIndexForTris2");
    commCompIndexForTris2();
    stop_clock("commCompIndexForTris2");

    start_clock("commMapping");
    commMapping(map_comm);
    stop_clock("commMapping");

    start_clock("adjustMapping");
    adjustMapping(map_comm);
    stop_clock("adjustMapping");

    map<int, int>::iterator itm;

    for (int i = 0; i < num_tris; i++)
    {
	itm = map_comm.find(comp_index[i]);

	if(itm != map_comm.end())
	    comp_index[i] = itm->second;
    }
    
}

void DropletAnalysis::adjustMapping(map<int,int> &map_for_adjust)
{
    map<int, int>::iterator it, itm;

    for (it = map_for_adjust.begin(); it != map_for_adjust.end(); it++)
    {
	itm = map_for_adjust.find(it->second);
	if (itm != map_for_adjust.end())
	    it->second = itm->second;
    }

}


void DropletAnalysis::outputAllPhyComps(void)
{
    char dirname[256];
    char *out_name = pSolver->front->out_name;

    if (pp_numnodes() > 1)
	sprintf(dirname, "%s/P-%s", out_name, right_flush(pp_mynode(), 4));
    else
	sprintf(dirname, "%s", out_name);

    sprintf(dirname, "%s/phycomps", dirname);

    if (!create_directory(dirname, YES))
    {
	screen ("Cannot create directory %s\n", dirname);
	clean_up(ERROR);
    }

    for (int i = 0; i < num_comps; i++)
    {
	char filename[256];
	int comp_index_cur = comp_list[i];

	sprintf(filename, "%s/intfc-nd%s-Comp%s", dirname, right_flush(pp_mynode(),4), right_flush(comp_index_cur, 4));

	FILE *outfile = fopen(filename, "w");

	fprintf(outfile, "# vtk DataFile Version 3.0\n");
	fprintf(outfile, "Interface number %d\n", comp_index);
	fprintf(outfile, "ASCII\n");
	fprintf(outfile, "DATASET POLYDATA\n");

	fprintf(outfile, "POINTS %d double\n", num_points);

	for (int j = 0; j < num_points; j++)
	    fprintf(outfile, "%-9g %-9g %-9g\n", ps[i][2]*cos(ps[i][0]), ps[i][2]*sin(ps[i][0]), ps[i][1]);

	int num_tris_for_this_comp = 0;

	for (int k = 0; k < num_tris; k++)
	{
	    if (comp_index[k] == comp_index_cur)
		num_tris_for_this_comp++;
	}

	fprintf(outfile, "POLYGONS %d %d \n", num_tris_for_this_comp, num_tris_for_this_comp*4);

	for (int k = 0; k < num_tris; k++)
	{
	    if (comp_index[k] == comp_index_cur)
		fprintf(outfile, "%d %d %d %d \n", 3, tris[k][0], tris[k][1], tris[k][2]);
	}

    }


}

void DropletAnalysis::outputAllTris_binary(void)
{
    char dirname[256];
    char *out_name = pSolver->front->out_name;

    sprintf(dirname, "%s", out_name);

    char filename[256];

    sprintf(filename, "%s/bintfc-nd%s", dirname, right_flush(pp_mynode(),4));

    int ival[4];
    double val[3];

    FILE *outfile = fopen(filename, "wb");

    if (hardware_is_little_endian())
    {
	ival[0] = endian_int_swap(num_points);
	fwrite(ival, sizeof(int), 1, outfile);

	for (int i = 0; i < num_points; i++)
	{
	    val[0] = endian_double_swap(ps[i][0]);
	    val[1] = endian_double_swap(ps[i][1]);
	    val[2] = endian_double_swap(ps[i][2]);
	    fwrite(val, sizeof(double), 3, outfile);
	}
	ival[0] = endian_int_swap(num_tris);
	fwrite(ival, sizeof(int), 1, outfile);

	for(int i = 0; i < num_tris; i++)
	{
	    ival[0] = endian_int_swap(tris[i][0]);
	    ival[1] = endian_int_swap(tris[i][1]);
	    ival[2] = endian_int_swap(tris[i][2]);
	    ival[3] = endian_int_swap(comp_index[i]);
	    fwrite(ival, sizeof(int), 4, outfile);
	}
	ival[0] = endian_int_swap(num_comps);
	fwrite(ival, sizeof(int), 1, outfile);

	for (int i = 0; i < num_comps; i++)
	{
	    ival[0] = endian_int_swap(comp_list[i]);
	    fwrite(ival, sizeof(int), 1, outfile);
	}
    }
    else
    {
	ival[0] = num_points;
	fwrite(ival, sizeof(int), 1, outfile);

	for (int i = 0; i < num_points; i++)
	{
	    val[0] = ps[i][0];
	    val[1] = ps[i][1];
	    val[2] = ps[i][2];
	    fwrite(val, sizeof(double), 3, outfile);
	}
	ival[0] = num_tris;
	fwrite(ival, sizeof(int), 1, outfile);

	for(int i = 0; i < num_tris; i++)
	{
	    ival[0] = tris[i][0];
	    ival[1] = tris[i][1];
	    ival[2] = tris[i][2];
	    ival[3] = comp_index[i];
	    fwrite(ival, sizeof(int), 4, outfile);
	}
	ival[0] = num_comps;
	fwrite(ival, sizeof(int), 1, outfile);

	for (int i = 0; i < num_comps; i++)
	{
	    ival[0] = comp_list[i];
	    fwrite(ival, sizeof(int), 1, outfile);
	}
    }


    fclose(outfile);

}




void DropletAnalysis::outputAllCells_binary(void)
{
    char dirname[256];
    char *out_name = pSolver->front->out_name;

    sprintf(dirname, "%s", out_name);

    char filename[256];

    sprintf(filename, "%s/bcell-nd%s", dirname, right_flush(pp_mynode(),4));

    double val[3];
    int ival[3];

    FILE *outfile = fopen(filename, "wb");

    if (hardware_is_little_endian())
    {
	ival[0] = endian_int_swap(num_cells);
	fwrite(ival, sizeof(int), 1, outfile);
	for(int i = 0; i < num_cells; i++)
	{
	    val[0] = endian_double_swap(c_centers[i][0]);
	    val[1] = endian_double_swap(c_centers[i][1]);
	    val[2] = endian_double_swap(c_centers[i][2]);
	    fwrite(val, sizeof(double), 3, outfile);
	}
	val[0] = endian_double_swap(grid_spacing[0]);
	val[1] = endian_double_swap(grid_spacing[1]);
	val[2] = endian_double_swap(grid_spacing[2]);
	fwrite(val, sizeof(double), 3, outfile);

	ival[0] = endian_int_swap(num_corners);
	fwrite(ival, sizeof(int), 1, outfile);
	for(int i = 0; i < num_corners; i++)
	{
	    val[0] = endian_double_swap(c_corners[i][0]);
	    val[1] = endian_double_swap(c_corners[i][1]);
	    val[2] = endian_double_swap(c_corners[i][2]);
	    fwrite(val, sizeof(double), 3, outfile);
	}
	ival[0] = endian_int_swap(gmax_cell[0]);
	ival[1] = endian_int_swap(gmax_cell[1]);
	ival[2] = endian_int_swap(gmax_cell[2]);
	fwrite(ival, sizeof(int), 3, outfile);
    }
    else
    {
	ival[0] = num_cells;
	fwrite(ival, sizeof(int), 1, outfile);
	for(int i = 0; i < num_cells; i++)
	{
	    val[0] = c_centers[i][0];
	    val[1] = c_centers[i][1];
	    val[2] = c_centers[i][2];
	    fwrite(val, sizeof(double), 3, outfile);
	}
	val[0] = grid_spacing[0];
	val[1] = grid_spacing[1];
	val[2] = grid_spacing[2];
	fwrite(val, sizeof(double), 3, outfile);

	ival[0] = num_corners;
	fwrite(ival, sizeof(int), 1, outfile);
	for (int i = 0; i < num_corners; i++)
	{
	    val[0] = c_corners[i][0];
	    val[1] = c_corners[i][1];
	    val[2] = c_corners[i][2];
	    fwrite(val, sizeof(double), 3, outfile);
	}
	ival[0] = gmax_cell[0];
	ival[1] = gmax_cell[1];
	ival[2] = gmax_cell[2];
	fwrite(ival, sizeof(int), 3, outfile);
    }

    fclose(outfile);
}


void DropletAnalysis::outputAllTris(void)
{
    char dirname[256];
    char *out_name = pSolver->front->out_name;

    sprintf(dirname, "%s", out_name);

    char filename[256];

    sprintf(filename, "%s/intfc-nd%s", dirname, right_flush(pp_mynode(),4));


    FILE *outfile = fopen(filename, "w");

    fprintf(outfile, "%d\n", num_points);
    for(int i = 0; i < num_points; i++)
	fprintf(outfile, "%20.16g %20.16g %20.16g\n", ps[i][0], ps[i][1], ps[i][2]);
    fprintf(outfile, "%d\n", num_tris);
    for(int i = 0;i < num_tris; i++)
	fprintf(outfile, "%d %d %d %d\n", tris[i][0], tris[i][1], tris[i][2], comp_index[i]);
    fprintf(outfile, "%d\n", num_comps);
    for(int i = 0;i < num_comps; i++)
	fprintf(outfile, "%d\n", comp_list[i]);

    fclose(outfile);

}




void DropletAnalysis::outputAllCells(void)
{
    char dirname[256];
    char *out_name = pSolver->front->out_name;

    sprintf(dirname, "%s", out_name);

    char filename[256];

    sprintf(filename, "%s/cell-nd%s", dirname, right_flush(pp_mynode(),4));


    FILE *outfile = fopen(filename, "w");

    fprintf(outfile, "%d\n", num_cells);
    for(int i = 0; i < num_cells; i++)
	fprintf(outfile, "%20.16g %20.16g %20.16g\n", c_centers[i][0], c_centers[i][1], c_centers[i][2]);
    fprintf(outfile, "%20.16g %20.16g %20.16g\n", grid_spacing[0], grid_spacing[1], grid_spacing[2]);

    fprintf(outfile, "%d\n", num_corners);
    for(int i = 0; i < num_corners; i++)
	fprintf(outfile, "%20.16g %20.16g %20.16g\n", c_corners[i][0], c_corners[i][1], c_corners[i][2]);

    fprintf(outfile, "%d %d %d\n", gmax_cell[0], gmax_cell[1], gmax_cell[2]);
    fclose(outfile);
}


void DropletAnalysis::outputHistogram(char* out_name)
{
    char filename[256];
    FILE *outfile;
    int my_id = pp_mynode();

    int i = 1;
    if (my_id == 0)
    {
	sprintf(filename, "%s/droplet-log", out_name);
	if(pSolver->front->step == 0)
	    outfile = fopen(filename, "w");
	else
	{
	    outfile = fopen(filename, "a");
	}

	(void) fprintf(outfile, "\n\n\ntime = %20.16g   step = %5d\n",
		pSolver->front->time, pSolver->front->step);

	fprintf(outfile, "\nNumber of connected physical components: %d", num_comps_global + 1);

	double total_vol_pos2 = 0.0;
	double total_vol_neg2 = 0.0;
	fprintf(outfile, "\nPositive volume:\n");

	for (map<int,double>::iterator itp = vol_comm_pos.begin(); itp != vol_comm_pos.end(); itp++)
	{
	    fprintf(outfile, "Comp Index = %d, volume = %20.16g\n", itp->first, itp->second);
	    total_vol_pos2 = total_vol_pos2 + itp->second;
	}
	fprintf(outfile, "\nNegative volume:\n");
	for (map<int,double>::iterator itn = vol_comm_neg.begin(); itn != vol_comm_neg.end(); itn++)
	{
	    fprintf(outfile, "Comp Index = %d, volume = %20.16g\n", itn->first, itn->second);
	    total_vol_neg2 = total_vol_neg2 + itn->second;
	}
	double total_vol_compare2 = total_vol_pos2 + total_vol_neg2;

	fprintf(outfile, "\nTotal volume from computation is %20.16g\n", total_vol_compare2);
   
    }
    printf("\n");
}


void DropletAnalysis::preInit(void)
{
    start_clock("setTempFrontAndRedist");
    pSolver->setTempFrontAndRedist();	
    if (debugging("trace"))
	printf("Passed setTempFrontAndRedist()\n");
    stop_clock("setTempFrontAndRedist");

    start_clock("SettingTriGlobalIndex");

    pSolver->setGlobalTriIndex();
    if (debugging("trace"))
	printf("Passed setGlobalTriIndex()\n");

    pSolver->commGlobalTriIndex();
    if (debugging("trace"))
	printf("Passed commGlobalTriIndex()\n");

    pSolver->commGlobalTriIndex2();
    if (debugging("trace"))
	printf("Passed commGlobalTriIndex2()\n");

    pSolver->checkUnsetGlobalIndex();
    if (debugging("trace"))
	printf("Passed checkUnsetGlobalIndex()\n");

    stop_clock("SettingTriGlobalIndex");
}


void DropletAnalysis::Init(void)
{
	CreatePoints();
	CreateTris();
	CreateCells();
	map_comm.clear();
	map_comm_pos.clear();
	map_comm_neg.clear();
	vol_comm_pos.clear();
	vol_comm_neg.clear();
	blank_proc = false;
	first_non_blank_cell = 0;
	num_comps_global = 0;
}

void DropletAnalysis::getCellCorners(int index, double *corners)
{
    int icoords[3];
    int id = index;

    for (int d = 0; d < 3; d++)
    {
	int G = gmax_cell[d];
	icoords[d] = id % G;
	id = (id - icoords[d])/G;
    }

    int xmin_index = icoords[0];
    int xmax_index = icoords[0] + 1;
    int ymin_index = icoords[1];
    int ymax_index = icoords[1] + 1;
    int zmin_index = icoords[2];
    int zmax_index = icoords[2] + 1;

    int Lcorner_index = (xmin_index + ymin_index*(gmax_cell[0]+1) + zmin_index*(gmax_cell[0]+1)*(gmax_cell[1]+1));
    int Ucorner_index = (xmax_index + ymax_index*(gmax_cell[0]+1) + zmax_index*(gmax_cell[0]+1)*(gmax_cell[1]+1));

    corners[0] = c_corners[Lcorner_index][0]; //xmin
    corners[1] = c_corners[Ucorner_index][0]; //xmax

    corners[2] = c_corners[Lcorner_index][1]; //ymin
    corners[3] = c_corners[Ucorner_index][1]; //ymax

    corners[4] = c_corners[Lcorner_index][2]; //zmin
    corners[5] = c_corners[Ucorner_index][2]; //zmax
}







