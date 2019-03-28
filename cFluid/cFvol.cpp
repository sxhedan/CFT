#include <cFluid.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define absmin(a,b) (((fabs(a)) < (fabs(b))) ? (a) : (b))
#define absmax(a,b) (((fabs(a)) > (fabs(b))) ? (a) : (b))

void tri_to_ctri(TRI*,CTRI*);
void copy_pt_to_cpt(POINT*,CPOINT*);
void copy_cpt_to_cpt(CPOINT*,CPOINT*);
void insert_tri_to_ctri_list(TRI*,CELL*);
void pret_tris(CELL*);
void set_polygons_in_cell(CELL*);
void init_cf_pts_and_edges(CFACE *cf);
bool set_nor(CPOINT*,CPOINT*,CPOINT*,double*);
bool tri_in_cell(CTRI*,CELL*);
bool point_in_cell(CPOINT*,CELL*);
void ctri_inter_with_cf(CTRI*,CFACE*,int*);
void set_crxps_and_edges_21(CTRI*,CFACE*,int*);
void set_crxps_and_edges_30(CTRI*,CFACE*,int*);
void set_crxps_and_edges_111(CTRI*,CFACE*,int*);
void set_crxps_and_edges_201(CTRI*,CFACE*,int*);
void add_edge(CEDGE*,CEDGE**);
void add_crxp_on_tri_edge(CPOINT*,CEDGE*);
bool edge_exists_in_edge_list(CEDGE*,CEDGE*);
bool point_on_cell_face(CPOINT*,CFACE*);
bool point_in_cell_face(CPOINT*,CFACE*);
bool point_in_tri_2d(CPOINT*,CTRI*,int);
double det2(double*,double*);
int find_crxps1(CPOINT*,CPOINT*,CFACE*,CPOINT*,CPOINT*);
int find_crxps2(CFACE*,int,CTRI*,CPOINT*,CPOINT*);
bool find_crxp_3d1(CPOINT*,CPOINT*,CFACE*,CPOINT*);
bool same_cpt(CPOINT*,CPOINT*);
bool same_cpt_double_tol(CPOINT*,CPOINT*);
bool same_undir_edge(CEDGE*,CEDGE*);
void print_edge(CEDGE*);
void print_edge_list(CEDGE*);
void print_point(CPOINT*);
void print_ctri(CTRI*);
void print_polygon(CPOLYGON*);
void construct_polygon_on_tri(CTRI*,CELL*);
void add_polyg_vertex(CPOINT*,CPOLYGON*);
bool recursively_find_next_vert(CPOINT*,CPOINT*,CPOLYGON*,CFACE*);
void push_back_vert(CPOINT*,CPOLYGON*);
void remove_last_vert(CPOLYGON*);
void add_p_as_nth_vert(CPOINT*,CPOLYGON*,int);
bool vertex_exists(CPOINT*,CPOLYGON*);
void remove_duplicated_edges(CPOLYGON*,CEDGE**);
void complete_polyg_undir_edges(CTRI*,CELL*);
void set_pic(CTRI*,CELL*,PIC*);
void add_new_edge(CPOINT*,CPOINT*,CEDGE**);
void add_new_dir_edge(CPOINT*,CPOINT*,CEDGE**);
void add_new_polyg(CPOLYGON*,CPOLYGON**);
void reverse_polyg(CPOLYGON*);
void set_directed_polyg(CPOLYGON*,double*);
void pop_edge(CEDGE**,CEDGE**);
bool find_next_endp(CPOINT,CEDGE**,CPOINT*);
bool same_nor_dir(double*,double*);
void construct_polygons_on_cfs(CELL*);
void set_polyg_edges_on_cf(CELL*,int);
int edge_dir(CEDGE*,CELL*);
bool point_on_cf_edge(CPOINT*,CFACE*,int*);
void add_crxp_on_cf_edge_sorted(CPOINT*,CFACE*,int,int);
bool pt_between_two_pts(CPOINT*,CPOINT*,CPOINT*);
void set_dir_cf(CELL*);
void set_dir_edges_by_nb_cf(CFACE*,CELL*,int);
bool is_nb_cf(int,int);
void set_polygons_on_cf(CELL*,int);
void add_next_vertices_on_ce(CPOLYGON*,CFACE*);
void add_prev_vertices_on_ce(CPOLYGON*,CFACE*);
void complete_polyg_on_cf(CPOLYGON*,CFACE*,CELL*);
void add_new_point(CPOINT*,CPOINT**);
bool the_ctri(CTRI*);
bool the_cpt(CPOINT*);
bool debug_same_cpt(CPOINT*,CPOINT*);
void set_polyhedrons_in_cell(CELL*);
void set_scps(CELL*);
void init_scp(SETOFCPOLYGS*);
void init_polygs_marks(CELL*);
void init_polygs_edges(CELL*);
void add_new_polygon(CPOLYGON*,CPOLYGON**);
void copy_polyg_to_polyg(CPOLYGON*,CPOLYGON*);
void construct_polyh(CPOLYHEDRON*,CPOLYGON*,CELL*);
void BFS_add_nb_polygs(CPOLYHEDRON*,CPOLYGON*,CELL*);
void find_nb_polygs(CPOLYGON**,CPOLYGON*,CELL*);
void add_polyg_to_scp(CPOLYGON*,SETOFCPOLYGS*);
void mark_polyg_edges(CPOLYGON*);
void BFS_add_nb_polygs(CPOLYGON*,SETOFCPOLYGS*,CPOLYGON*);
bool find_polyg_with_edge(CPOLYGON**,CEDGE*,CPOLYGON*);
void add_scp_boundary(CEDGE*,SETOFCPOLYGS*);
void insert_scp(SETOFCPOLYGS*,SETOFCPOLYGS**);
void reset_scp_boundaries(SETOFCPOLYGS*);
void complete_bdry_with_init_edge(CBOUNDARY*,CEDGE**);
void set_scp_dir(SETOFCPOLYGS*,CELL*);
void init_polyh(CPOLYHEDRON*,CELL*);
void add_scp_to_polyh(SETOFCPOLYGS*,CPOLYHEDRON*);
void add_scp_polygs_to_polyh(SETOFCPOLYGS*,CPOLYHEDRON*);
void add_scp_bdries_to_polyh(SETOFCPOLYGS*,CPOLYHEDRON*);
void add_new_bdry(CBOUNDARY*,CBOUNDARY**);
void add_scp_with_ibdries(CPOLYHEDRON*);
void inverse_bdry(CBOUNDARY*,CBOUNDARY*);
bool find_scp_with_bdry(SETOFCPOLYGS**,CBOUNDARY*,CPOLYHEDRON*);
bool find_bdry_with_edge(CBOUNDARY*,CEDGE*);
bool same_dir_edge(CEDGE*,CEDGE*);
void init_scps_bdries_marks(CELL*);
void init_scps_marks(CELL*);
void init_scpics_marks(CELL*);
void set_pf_areas(CELL*);
double area_of_polygon(CPOLYGON*);
void cross_product(double*,double*,double*);
double dot_product(double*,double*);
void set_polyh_vol(CELL*);
double vol_of_polyh(CPOLYHEDRON*);
void set_polyg_nor(CPOLYGON*);
void find_polyh_face_oncf_with_max_area(CPOLYHEDRON*,CPOLYGON**);
void find_nb_cell_with_face(CELL*,CPOLYGON*,int*);
void merge_polyh_to_cell(CPOLYHEDRON*,CELL*);
void cell_to_polyh(CELL*,CPOLYHEDRON*);
void init_new_polyg(CPOLYGON*);
void remove_short_edges_on_cf(CFACE*);
void modify_endps_of_edges(CFACE*,CPOINT*,CPOINT*);
void set_polyh_nb_cells(CPOLYHEDRON*);
void add_nbc_by_face(CPOLYHEDRON*,CPOLYGON*);
void sort_nbcs(CPOLYHEDRON*);
bool face_in_polyh(CPOLYGON*,CPOLYHEDRON*);
bool same_polyg(CPOLYGON*,CPOLYGON*);
CPOINT *prev_v(CPOLYGON*,CPOINT*);

bool DEBUGDAN = false;

void G_CARTESIAN::cft_init_cut_cells(TS_LEVEL ts)
{
	int i, j, k, l, ll, index;
	CELL *c;

	num_cells = 1;
	for (i = 0; i < dim; i++)
	    num_cells *= (top_gmax[i]+1);
	FT_VectorMemoryAlloc((POINTER*)&cells, num_cells, sizeof(CELL));

	switch (ts)
	{
	    case OLDTS:
		cells_old = cells;
		break;
	    case HALFTS:
		cells_halft = cells;
		break;
	    case NEWTS:
		cells_new = cells;
		break;
	    default:
		printf("ERROR in cft_init_cut_cells: unknown ts level.\n");
		clean_up(ERROR);
	}

	//Initialize cells in buffer zone but don't use them
	//start with (lbuf[i] ? lbuf[i] : 1)

	//FT_VectorMemoryAlloc((POINTER*)&cells_old, num_cells, sizeof(CELL));
	//FT_VectorMemoryAlloc((POINTER*)&cells_halft, num_cells, sizeof(CELL));
	//FT_VectorMemoryAlloc((POINTER*)&cells_new, num_cells, sizeof(CELL));

	//for (k = 0; k <= top_gmax[2]; k++)
	//for (j = 0; j <= top_gmax[1]; j++)
	//for (i = 0; i <= top_gmax[0]; i++)
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);

	    c = &(cells[index]);
	    c->icrds[0] = i;
	    c->icrds[1] = j;
	    c->icrds[2] = k;
	    c->ctri_polygs = NULL;
	    c->cf_polygs = NULL;
	    c->scpocs = NULL;
	    c->scpics = NULL;
	    for (l = 0; l < 3; l++)
	    {
		c->celll[l] = top_grid->L[l] + (c->icrds[l]-0.5)*top_grid->h[l];
		c->cellu[l] = c->celll[l] + top_grid->h[l];
	    }

	    //initialize cell faces
	    for (l = 0; l < 6; l++)
	    {
		for (ll = 0; ll < 3; ll++)
		{
		    c->faces[l].l[ll] = c->celll[ll];
		    c->faces[l].u[ll] = c->cellu[ll];
		}
	    }
	    c->faces[0].dir = 2;
	    c->faces[0].side = -1;
	    c->faces[0].u[2] = c->celll[2];
	    init_cf_pts_and_edges(&(c->faces[0]));
	    c->faces[1].dir = 1;
	    c->faces[1].side = -1;
	    c->faces[1].u[1] = c->celll[1];
	    init_cf_pts_and_edges(&(c->faces[1]));
	    c->faces[2].dir = 0;
	    c->faces[2].side = 1;
	    c->faces[2].l[0] = c->cellu[0];
	    init_cf_pts_and_edges(&(c->faces[2]));
	    c->faces[3].dir = 1;
	    c->faces[3].side = 1;
	    c->faces[3].l[1] = c->cellu[1];
	    init_cf_pts_and_edges(&(c->faces[3]));
	    c->faces[4].dir = 0;
	    c->faces[4].side = -1;
	    c->faces[4].u[0] = c->celll[0];
	    init_cf_pts_and_edges(&(c->faces[4]));
	    c->faces[5].dir = 2;
	    c->faces[5].side = 1;
	    c->faces[5].l[2] = c->cellu[2];
	    init_cf_pts_and_edges(&(c->faces[5]));

	    /*
	    c = &(cells_halft[index]);
	    c->icrds[0] = i;
	    c->icrds[1] = j;
	    c->icrds[2] = k;
	    c->ctri_polygs = NULL;
	    c->cf_polygs = NULL;
	    c->scpocs = NULL;
	    c->scpics = NULL;
	    for (l = 0; l < 3; l++)
	    {
		c->celll[l] = top_grid->L[l] + (c->icrds[l]-0.5)*top_grid->h[l];
		c->cellu[l] = c->celll[l] + top_grid->h[l];
	    }

	    //initialize cell faces
	    for (l = 0; l < 6; l++)
	    {
		for (ll = 0; ll < 3; ll++)
		{
		    c->faces[l].l[ll] = c->celll[ll];
		    c->faces[l].u[ll] = c->cellu[ll];
		}
	    }
	    c->faces[0].dir = 2;
	    c->faces[0].side = -1;
	    c->faces[0].u[2] = c->celll[2];
	    init_cf_pts_and_edges(&(c->faces[0]));
	    c->faces[1].dir = 1;
	    c->faces[1].side = -1;
	    c->faces[1].u[1] = c->celll[1];
	    init_cf_pts_and_edges(&(c->faces[1]));
	    c->faces[2].dir = 0;
	    c->faces[2].side = 1;
	    c->faces[2].l[0] = c->cellu[0];
	    init_cf_pts_and_edges(&(c->faces[2]));
	    c->faces[3].dir = 1;
	    c->faces[3].side = 1;
	    c->faces[3].l[1] = c->cellu[1];
	    init_cf_pts_and_edges(&(c->faces[3]));
	    c->faces[4].dir = 0;
	    c->faces[4].side = -1;
	    c->faces[4].u[0] = c->celll[0];
	    init_cf_pts_and_edges(&(c->faces[4]));
	    c->faces[5].dir = 2;
	    c->faces[5].side = 1;
	    c->faces[5].l[2] = c->cellu[2];
	    init_cf_pts_and_edges(&(c->faces[5]));

	    c = &(cells_new[index]);
	    c->icrds[0] = i;
	    c->icrds[1] = j;
	    c->icrds[2] = k;
	    c->ctri_polygs = NULL;
	    c->cf_polygs = NULL;
	    c->scpocs = NULL;
	    c->scpics = NULL;
	    for (l = 0; l < 3; l++)
	    {
		c->celll[l] = top_grid->L[l] + (c->icrds[l]-0.5)*top_grid->h[l];
		c->cellu[l] = c->celll[l] + top_grid->h[l];
	    }

	    //initialize cell faces
	    for (l = 0; l < 6; l++)
	    {
		for (ll = 0; ll < 3; ll++)
		{
		    c->faces[l].l[ll] = c->celll[ll];
		    c->faces[l].u[ll] = c->cellu[ll];
		}
	    }
	    c->faces[0].dir = 2;
	    c->faces[0].side = -1;
	    c->faces[0].u[2] = c->celll[2];
	    init_cf_pts_and_edges(&(c->faces[0]));
	    c->faces[1].dir = 1;
	    c->faces[1].side = -1;
	    c->faces[1].u[1] = c->celll[1];
	    init_cf_pts_and_edges(&(c->faces[1]));
	    c->faces[2].dir = 0;
	    c->faces[2].side = 1;
	    c->faces[2].l[0] = c->cellu[0];
	    init_cf_pts_and_edges(&(c->faces[2]));
	    c->faces[3].dir = 1;
	    c->faces[3].side = 1;
	    c->faces[3].l[1] = c->cellu[1];
	    init_cf_pts_and_edges(&(c->faces[3]));
	    c->faces[4].dir = 0;
	    c->faces[4].side = -1;
	    c->faces[4].u[0] = c->celll[0];
	    init_cf_pts_and_edges(&(c->faces[4]));
	    c->faces[5].dir = 2;
	    c->faces[5].side = 1;
	    c->faces[5].l[2] = c->cellu[2];
	    init_cf_pts_and_edges(&(c->faces[5]));
	    */
	}

	return;
}

void G_CARTESIAN::cft_set_cut_cells(TS_LEVEL ts)
{
	switch (ts)
	{
	    case OLDTS:
		cells = cells_old;
		break;
	    case HALFTS:
		cells = cells_halft;
		break;
	    case NEWTS:
		cells = cells_new;
		break;
	    default:
		printf("ERROR in cft_set_cut_cells: unknown ts level.\n");
		clean_up(ERROR);
	}

	cft_init_tris_in_cells();

	cft_set_cell_polygons();

	cft_construct_cell_polyhedrons();

	cft_set_polyhs_comps();

	cft_set_cut_cell_vol();

	return;
}

bool debugmerge = false;

void G_CARTESIAN::cft_merge_polyhs(TS_LEVEL ts)
{
	int i, j, k, index;
	CELL *c, *nb_cell;
	CPOLYHEDRON *polyh, *nb_polyh;
	CPOLYGON *max_face;
	NBCELL *nbc;
	bool merged, unmerged;

	switch (ts)
	{
	    case OLDTS:
		cells = cells_old;
		break;
	    case HALFTS:
		cells = cells_halft;
		break;
	    case NEWTS:
		cells = cells_new;
		break;
	    default:
		printf("ERROR in cft_set_cet_cells: unknown ts level.\n");
		clean_up(ERROR);
	}

	//set comp based on volumes of polyhedrons in cells
	cft_set_comp();

	//initialize pam
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);
	    c->pams = NULL;
	    c->merged = FALSE;

	    //polyh->merged is set in initialization
	}

	//first merge
	unmerged = false;
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);

	    //debugdan	FIXME
	    /*
	    if (i == 7 && j == 4 && k ==20)
		debugmerge = true;
	    else
		debugmerge = false;
	    */
	    //debugdan	FIXME

	    if (c->cut == FALSE)
	    {
		polyh = c->polyhs;
		merge_polyh_to_cell(polyh,c);
		continue;
	    }

	    polyh = c->polyhs;
	    while (polyh)
	    {
		if (polyh->comp == c->comp)
		{
		    merge_polyh_to_cell(polyh,c);
		    c->merged = TRUE;
		}
		else
		{
		    merged = false;
		    set_polyh_nb_cells(polyh);
		    nbc = polyh->sorted_nbcs;
		    while (nbc)
		    {
			index = d_index3d(nbc->inbc[0],nbc->inbc[1],nbc->inbc[2],top_gmax);
			nb_cell = &(cells[index]);
			if (nb_cell->comp == polyh->comp)
			{
			    merge_polyh_to_cell(polyh,nb_cell);
			    nb_cell->merged = TRUE;
			    merged = true;
			    break;
			}
			nbc = nbc->next;
		    }
		    if (!merged)
			unmerged = true;
		}
		polyh = polyh->next;
	    }
	}

	//second merge
	int MAX_MERGE = 3;
	int count = 0;
	while (unmerged)
	{
	    unmerged = false;

	    for (k = imin[2]; k <= imax[2]; k++)
	    for (j = imin[1]; j <= imax[1]; j++)
	    for (i = imin[0]; i <= imax[0]; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		c = &(cells[index]);

		polyh = c->polyhs;
		while (polyh)
		{
		    if (polyh->merged == true)
		    {
			polyh = polyh->next;
			continue;
		    }
		    
		    merged = false;
		    nbc = polyh->sorted_nbcs;
		    while (nbc)
		    {
			index = d_index3d(nbc->inbc[0],nbc->inbc[1],nbc->inbc[2],top_gmax);
			nb_cell = &(cells[index]);
			nb_polyh = nb_cell->polyhs;
			while (nb_polyh)
			{
			    if (nb_polyh->comp == polyh->comp)
			    {
				if (face_in_polyh(nbc->face,nb_polyh))
				{
				    if (nb_polyh->merged == true)
				    {
					merge_polyh_to_cell(polyh,nb_polyh->pam->targetc);
					merged = true;
					break;
				    }
				}
			    }
			    nb_polyh = nb_polyh->next;
			}
			nbc = nbc->next;
		    }
		    if (!merged)
			unmerged = true;

		    polyh = polyh->next;
		}
	    }

	    count++;
	    if (count > MAX_MERGE)
		break;
	}

	/*
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);

	    polyh = c->polyhs;
	    while (polyh)
	    {
		if (polyh->merged == true)
		{
		    polyh = polyh->next;
		    continue;
		}

		
		nbc = polyh->sorted_nbcs;
		while (nbc)
		{
		    index = d_index3d(nbc->inbc[0],nbc->inbc[1],nbc->inbc[2],top_gmax);
		    nb_cell = &(cells[index]);
		    nb_polyh = nb_cell->polyhs;
		    while (nb_polyh)
		    {
			if (nb_polyh->comp == polyh->comp)
			{
			    if (face_in_polyh(nbc->face,nb_polyh))
			    {
				if (nb_polyh->merged == true)
				{
				    merge_polyh_to_cell(polyh,nb_polyh->pam->targetc);
				    break;
				}
			    }
			}
			nb_polyh = nb_polyh->next;
		    }
		    nbc = nbc->next;
		}

		polyh = polyh->next;
	    }
	}
	*/

	//check if there's any polyh unmerged	FIXME
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);

	    polyh = c->polyhs;
	    while (polyh)
	    {
		if (polyh->merged == false)
		{
		    printf("%d %d %d polyh is unmerged.\n", i, j, k);
		}
		polyh = polyh->next;
	    }
	}

	//debugdan	FIXME
	//printf("END OF cft_merge_polyhs().\n");
	//exit(0);
	//debugdan	FIXME

	return;
}

bool face_in_polyh(
	CPOLYGON	*face,
	CPOLYHEDRON	*polyh)
{
	CPOLYGON *polyg = polyh->faces;

	while (polyg)
	{
	    //test?	TODO
	    if (same_polyg(face,polyg))
		return true;
	    polyg = polyg->next;
	}

	return false;
}

bool same_polyg(
	CPOLYGON	*polyg0,
	CPOLYGON	*polyg1)
{
	CPOINT *p00, *p01, *p10, *p11;
	int dir;
	bool found;

	p00 = polyg0->vertices;
	p01 = p00->next;
	if (!p00 || !p01)
	{
	    printf("ERROR in same_polyg().\n");
	    clean_up(ERROR);
	}

	p10 = polyg1->vertices;
	found = false;
	while (p10)
	{
	    if (same_cpt(p00,p10))
	    {
		found = true;
		break;
	    }
	    else
		p10 = p10->next;
	}
	if (!found)
	    return false;

	p11 = polyg1->vertices;
	found = false;
	while (p11)
	{
	    if (same_cpt(p01,p11))
	    {
		found = true;
		break;
	    }
	    else
		p11 = p11->next;
	}
	if (!found)
	    return false;

	if ((p11 == p10->next) ||
	    (p11 == polyg1->vertices && p10->next == NULL))
	{
	    dir = 1;
	}
	if ((p10 == p11->next) ||
	    (p10 == polyg1->vertices && p11->next == NULL))
	{
	    dir = -1;
	}
	else
	{
	    return false;
	}

	while (p11 != p10)
	{
	    if (!same_cpt(p01,p11))
		return false;

	    //next p01
	    if (p01->next)
		p01 = p01->next;
	    else
		p01 = polyg0->vertices;

	    //next p11
	    if (dir == 1)
	    {
		if (p11->next)
		    p11 = p11->next;
		else
		    p11 = polyg1->vertices;
	    }
	    else if (dir == -1)
	    {
		p11 = prev_v(polyg1,p11);
	    }
	}

	return true;
}

CPOINT *prev_v(
	CPOLYGON	*polyg,
	CPOINT		*p)
{
	CPOINT *prevp = polyg->vertices;

	while (prevp->next)
	{
	    if (prevp->next == p)
		return prevp;
	    else
		prevp = prevp->next;
	}

	return prevp;
}

//sort polyh's faces on cell faces based on areas
//and set corresponding neighbor cells.
void set_polyh_nb_cells(CPOLYHEDRON *polyh)
{
	double max_area = 0;
	CPOLYGON *face;

	face = polyh->faces;
	while (face)
	{
	    if (face->oncf)
	    {
		add_nbc_by_face(polyh,face);
	    }
	    face = face->next;
	}

	//debugdan	FIXME
	/*
	int count = 0;
	NBCELL *nbc = polyh->sorted_nbcs;
	while (nbc)
	{
	    count++;
	    nbc = nbc->next;
	}
	if (count > 2)
	{
	    printf("Before sort_nbcs: ");
	    nbc = polyh->sorted_nbcs;
	    while (nbc)
	    {
		printf("%e ", nbc->face->area);
		nbc = nbc->next;
	    }
	    printf("\n");
	}
	*/
	//debugdan	FIXME

	sort_nbcs(polyh);

	//debugdan	FIXME
	/*
	if (count > 2)
	{
	    printf("After sort_nbcs: ");
	    nbc = polyh->sorted_nbcs;
	    while (nbc)
	    {
		printf("%e ", nbc->face->area);
		nbc = nbc->next;
	    }
	    printf("\n");
	}
	*/
	//debugdan	FIXME

	return;
}

void sort_nbcs(CPOLYHEDRON *polyh)
{
	NBCELL *tmp, *endsort, *cur, *prev, *max, *maxprev;
	bool first = true;

	FT_ScalarMemoryAlloc((POINTER*)&(tmp),sizeof(NBCELL));
	tmp->next = polyh->sorted_nbcs;
	endsort = tmp;

	while (endsort->next)
	{
	    cur = endsort->next;
	    max = endsort->next;
	    prev = endsort;
	    maxprev = endsort;
	    while (cur)
	    {
		if (cur->face->area > max->face->area)
		{
		    max = cur;
		    maxprev = prev;
		}
		prev = cur;
		cur = cur->next;
	    }
	    if (first)
	    {
		polyh->sorted_nbcs = max;
		first = false;
	    }
	    if (max == endsort->next)
	    {
		endsort = endsort->next;
	    }
	    else
	    {
		maxprev->next = max->next;
		max->next = endsort->next;
		endsort->next = max;
		endsort = endsort->next;
	    }
	}

	free(tmp);

	return;
}

void add_nbc_by_face(
	CPOLYHEDRON	*polyh,
	CPOLYGON	*face)
{
	NBCELL *nbc;

	FT_ScalarMemoryAlloc((POINTER*)&(nbc),sizeof(NBCELL));
	nbc->face = face;
	nbc->next = polyh->sorted_nbcs;
	polyh->sorted_nbcs = nbc;
	find_nb_cell_with_face(polyh->cell,face,nbc->inbc);

	return;
}

void G_CARTESIAN::cft_set_comp()
{
	int i, j, k, index;
	double cvol = top_h[0]*top_h[1]*top_h[2];
	CELL *c;
	CPOLYHEDRON *polyh;

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);
	    c->vol[0] = 0.0;
	    c->vol[1] = 0.0;

	    polyh = c->polyhs;
	    while (polyh)
	    {
		if (polyh->comp == GAS_COMP1)
		    c->vol[0] += polyh->vol;
		else if (polyh->comp == GAS_COMP2)
		    c->vol[1] += polyh->vol;
		polyh = polyh->next;
	    }

	    if (c->vol[0] >= 0.5*cvol)
		c->comp = GAS_COMP1;
	    else
		c->comp = GAS_COMP2;
	}

	return;
}

void cell_to_polyh(
	CELL		*c,
	CPOLYHEDRON	*polyh)
{
	polyh->faces = NULL;

	if (c->comp == GAS_COMP1)
	{
	    polyh->scp_dir == OUTW;
	}
	else if (c->comp == GAS_COMP2)
	{
	    polyh->scp_dir == INW;
	}

	polyh->iscell = TRUE;
	polyh->cell = c;

	return;
}

void merge_polyh_to_cell(
	CPOLYHEDRON	*polyh,
	CELL		*c)
{
	CPAM *pam;

	FT_ScalarMemoryAlloc((POINTER*)&(pam),sizeof(CPAM));
	pam->polyh = polyh;
	pam->targetc = c;
	pam->next = NULL;
	for (int ii = 0; ii < 3; ii++)
	    pam->mdir[ii] = c->icrds[ii] - polyh->cell->icrds[ii];

	pam->next = c->pams;
	c->pams = pam;
	polyh->merged = true;
	polyh->pam = pam;
	//c->merged = TRUE;

	return;
}

void find_nb_cell_with_face(
	CELL		*c,
	CPOLYGON	*face,
	int		*nb_icrds)
{
	int i;
	CPOINT *p;
	double tol = 1e-12;

	for (i = 0; i < 3; i++)
	    nb_icrds[i] = c->icrds[i];

	for (i = 0; i < 3; i++)
	{
	    p = face->vertices;
	    while (p)
	    {
		if (fabs(p->crds[i] - c->celll[i]) > tol)
		{
		    break;
		}
		p = p->next;
	    }

	    if (p == NULL)
	    {
		nb_icrds[i]--;
		return;
	    }

	    p = face->vertices;
	    while (p)
	    {
		if (fabs(p->crds[i] - c->cellu[i]) > tol)
		{
		    break;
		}
		p = p->next;
	    }

	    if (p == NULL)
	    {
		nb_icrds[i]++;
		return;
	    }
	}

	printf("ERROR in find_nb_cell_with_face(): can't find nb cell.\n");
	clean_up(ERROR);
}

void find_polyh_face_oncf_with_max_area(
	CPOLYHEDRON	*polyh,
	CPOLYGON	**max_face)
{
	double max_area = 0;
	CPOLYGON *face;

	face = polyh->faces;
	while (face)
	{
	    if (face->oncf)
	    {
		if (face->area > max_area)
		{
		    max_area = face->area;
		    *max_face = face;
		}
	    }
	    face = face->next;
	}

	return;
}

/*
void G_CARTESIAN::cvol()
{
	printf("Enter cvol().\n");

	if (dim != 3)
	{
	    printf("cvol for 2D is NOT implemented yet.\n");
	    clean_up(ERROR);
	}

	printf("Init cells.\n");
	cft_init_grid_cells();

	printf("Init ctris.\n");
	init_tris_in_cells();

	printf("Init polygons.\n");
	set_cell_polygons();

	printf("Construct polyhedrons.\n");
	construct_cell_polyhedrons();

	printf("Calculate volume for polyhedrons.\n");
	cut_cell_vol();

	printf("Leave cvol().\n");

	return;
}
*/
/*
void G_CARTESIAN::cft_init_grid_cells()
{
	int i, j, k, l, ll, index;
	CELL *c;
	double tol;

	//FIXME
	//tol = max(max(top_grid->h[0], top_grid->h[1]), top_grid->h[2]);
	//tol = tol*tol*tol;
	//FIXME

	num_cells = 1;
	for (i = 0; i < dim; i++)
	    num_cells *= (top_gmax[i]+1);
	//Initialize cells in buffer zone but don't use them
	//start with (lbuf[i] ? lbuf[i] : 1)

	FT_VectorMemoryAlloc((POINTER*)&cells, num_cells, sizeof(CELL));

	//for (k = 0; k <= top_gmax[2]; k++)
	//for (j = 0; j <= top_gmax[1]; j++)
	//for (i = 0; i <= top_gmax[0]; i++)
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);
	    c->icrds[0] = i;
	    c->icrds[1] = j;
	    c->icrds[2] = k;
	    c->ctri_polygs = NULL;
	    c->cf_polygs = NULL;
	    c->scpocs = NULL;
	    c->scpics = NULL;
	    for (l = 0; l < 3; l++)
	    {
		c->celll[l] = top_grid->L[l] + (c->icrds[l]-0.5)*top_grid->h[l];
		c->cellu[l] = c->celll[l] + top_grid->h[l];
	    }

	    //initialize cell faces
	    for (l = 0; l < 6; l++)
	    {
		for (ll = 0; ll < 3; ll++)
		{
		    c->faces[l].l[ll] = c->celll[ll];
		    c->faces[l].u[ll] = c->cellu[ll];
		}
	    }
	    c->faces[0].dir = 2;
	    c->faces[0].side = -1;
	    c->faces[0].u[2] = c->celll[2];
	    init_cf_pts_and_edges(&(c->faces[0]));
	    c->faces[1].dir = 1;
	    c->faces[1].side = -1;
	    c->faces[1].u[1] = c->celll[1];
	    init_cf_pts_and_edges(&(c->faces[1]));
	    c->faces[2].dir = 0;
	    c->faces[2].side = 1;
	    c->faces[2].l[0] = c->cellu[0];
	    init_cf_pts_and_edges(&(c->faces[2]));
	    c->faces[3].dir = 1;
	    c->faces[3].side = 1;
	    c->faces[3].l[1] = c->cellu[1];
	    init_cf_pts_and_edges(&(c->faces[3]));
	    c->faces[4].dir = 0;
	    c->faces[4].side = -1;
	    c->faces[4].u[0] = c->celll[0];
	    init_cf_pts_and_edges(&(c->faces[4]));
	    c->faces[5].dir = 2;
	    c->faces[5].side = 1;
	    c->faces[5].l[2] = c->cellu[2];
	    init_cf_pts_and_edges(&(c->faces[5]));
	}

	return;
}
*/

void init_cf_pts_and_edges(CFACE *cf)
{
	int i, j, count;
	double tmp1[3][2];
	int tmp2[2][4] = {{0,1,1,0}, {0,0,1,1}};

	for (i = 0; i < 3; i++)
	{
	    tmp1[i][0] = cf->l[i];
	    tmp1[i][1] = cf->u[i];
	}

	for (i = 0; i < 4; i++)
	{
	    count = 0;
	    for (j = 0; j < 3; j++)
	    {
		if (j == cf->dir)
		{
		    cf->pts[i].crds[j] = cf->l[j];
		    continue;
		}
		cf->pts[i].crds[j] = tmp1[j][tmp2[count][i]];
		count++;
	    }
	}

	for (i = 0; i < 4; i++)
	{
	    cf->edges[i].endp[0] = cf->pts[i];
	    cf->edges[i].endp[1] = cf->pts[(i+1)%4];
	}

	return;
}

void G_CARTESIAN::cft_init_tris_in_cells()
{
	INTERFACE	*intfc = front->interf;
	SURFACE		**s;
	TRI		*tri;
	POINT		**pt;
	CTRI		*ctri;
	double		*pt0, *pt1, *pt2;
	double		mincrds[3], maxcrds[3];
	int		i, j, k, index;
	int		minicrds[3], maxicrds[3];
	double		tol = 1e-12;

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (Boundary(*s))
	    {
		continue;
	    }
	    for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	    {
		pt = Point_of_tri(tri);
		pt0 = Coords(pt[0]);
		pt1 = Coords(pt[1]);
		pt2 = Coords(pt[2]);

		for (i = 0; i < dim; i++)
		{
		    mincrds[i] = min(min(pt0[i], pt1[i]), pt2[i]);
		    maxcrds[i] = max(max(pt0[i], pt1[i]), pt2[i]);
		    minicrds[i] = floor((mincrds[i]-top_grid->L[i])/top_grid->h[i]+0.5-tol);
		    maxicrds[i] = floor((maxcrds[i]-top_grid->L[i])/top_grid->h[i]+0.5+tol);
		    if (minicrds[i] < 0)	//check if this happens	TODO
			minicrds[i] = 0;
		    if (maxicrds[i] > top_gmax[i])	//check if this happens	TODO
			maxicrds[i] = top_gmax[i];
		}

		for (k = minicrds[2]; k <= maxicrds[2]; k++)
		for (j = minicrds[1]; j <= maxicrds[1]; j++)
		for (i = minicrds[0]; i <= maxicrds[0]; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    insert_tri_to_ctri_list(tri,&(cells[index]));
		}
	    }
	}

	return;
}

void G_CARTESIAN::cft_set_cell_polygons()
{
	int i, j, k, index;
	CELL *c;

	//for (k = 0; k <= top_gmax[2]; k++)
	//for (j = 0; j <= top_gmax[2]; j++)
	//for (i = 0; i <= top_gmax[2]; i++)
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);

	    pret_tris(c);
	    set_polygons_in_cell(c);
	}

	return;
}

void pret_tris(CELL *c)
{
	int i, j, k;
	int pret[3];
	CTRI *ctri1, *ctri2;
	CPOINT *p1, *p2, p3;
	double tol = 1e-8;

	for (ctri1 = c->ctris; ctri1; ctri1 = ctri1->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p1 = &(ctri1->pts[i]);
		p3 = ctri1->pts[i];
		for (j = 0; j < 3; j++)
		{
		    pret[j] = 0;
		    if (fabs(p1->crds[j] - c->celll[j]) < tol)
		    {
			p1->crds[j] = c->celll[j];
			pret[j] = -1;
		    }
		    else if (fabs(p1->crds[j] - c->cellu[j]) < tol)
		    {
			p1->crds[j] = c->cellu[j];
			pret[j] = 1;
		    }
		}

		for (ctri2 = c->ctris; ctri2; ctri2 = ctri2->next)
		{
		    if (ctri2 == ctri1)
			continue;
		    for (j = 0; j < 3; j++)
		    {
			p2 = &(ctri2->pts[j]);
			if (same_cpt(p2, &p3))
			{
			    for (k = 0; k < 3; k++)
			    {
				if (pret[k] == -1)
				    p2->crds[k] = c->celll[k];
				else if (pret[k] == 1)
				    p2->crds[k] = c->cellu[k];
			    }
			}
		    }
		}
	    }
	}

	return;
}

void tri_to_ctri(
	TRI	*tri,
	CTRI	*ctri)
{
	int i;

	for (i = 0; i < 3; i++)
	    copy_pt_to_cpt(Point_of_tri(tri)[i],&(ctri->pts[i]));
	for (i = 0; i < 3; i++)
	{
	    ctri->edges[i].ncrxps = 0;
	    ctri->edges[i].endp[0] = ctri->pts[i];
	    ctri->edges[i].endp[1] = ctri->pts[(i+1)%3];
	}
	set_nor(&(ctri->pts[0]),&(ctri->pts[1]),&(ctri->pts[2]),ctri->nor);

	FT_ScalarMemoryAlloc((POINTER*)&(ctri->polyg),sizeof(CPOLYGON));
	init_new_polyg(ctri->polyg);
	//ctri->polyg->vertices = NULL;
	//ctri->polyg->undir_edges = NULL;

	return;
}

bool set_nor(
	CPOINT	*p0,
	CPOINT	*p1,
	CPOINT	*p2,
	double	*nor)
{
	int i;
	double v1[3], v2[3], l1, l2, l;
	double tol = 1e-12;

	for (i = 0; i < 3; i++)
	{
	    v1[i] = p1->crds[i] - p0->crds[i];
	    v2[i] = p2->crds[i] - p1->crds[i];
	}
	l1 = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
	l2 = sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
	for (i = 0; i < 3; i++)
	{
	    v1[i] = v1[i]/l1;
	    v2[i] = v2[i]/l2;
	}

	nor[0] = v1[1]*v2[2] - v1[2]*v2[1];
	nor[1] = -v1[0]*v2[2] + v1[2]*v2[0];
	nor[2] = v1[0]*v2[1] - v1[1]*v2[0];

	if (fabs(nor[0]) < tol && fabs(nor[1]) < tol && fabs(nor[2]) < tol)
	    return NO;

	l = sqrt(nor[0]*nor[0]+nor[1]*nor[1]+nor[2]*nor[2]);
	for (i = 0; i < 3; i++)
	    nor[i] = nor[i]/l;

	return YES;
}

void copy_pt_to_cpt(
	POINT	*pt,
	CPOINT	*cpt)
{
	int i;

	for (i = 0; i < 3; i++)
	    cpt->crds[i] = Coords(pt)[i];
	return;
}

void copy_cpt_to_cpt(
	CPOINT	*cpt1,
	CPOINT	*cpt2)
{
	int i;

	for (i = 0; i < 3; i++)
	    cpt2->crds[i] = cpt1->crds[i];
	return;
}

void insert_tri_to_ctri_list(
	TRI	*tri,
	CELL	*c)
{
	CTRI *ctri;

	FT_ScalarMemoryAlloc((POINTER*)&ctri,sizeof(CTRI));
	tri_to_ctri(tri, ctri);

	if (c->ctris)
	{
	    ctri->next = c->ctris;
	    c->ctris = ctri;
	}
	else
	    c->ctris = ctri;

	return;
}

void set_polygons_in_cell(CELL *c)
{
	int i;
	int sides[3], ss;
	CTRI *ctri;
	CFACE *cf;

	for (ctri = c->ctris; ctri; ctri = ctri->next)
	{
	    for (i = 0; i < 6; i++)
	    {
		cf = &(c->faces[i]);
		ctri_inter_with_cf(ctri,cf,sides);
		ss = sides[0] + sides[1] + sides[2];
		switch (ss)
		{
		    case 3:
		    case 300:
		    case 210:
		    case 12:
		    case 120:
		    case 30:
			break;
		    case 21:
			set_crxps_and_edges_21(ctri,cf,sides);
			break;
		    case 111:
			set_crxps_and_edges_111(ctri,cf,sides);
			break;
		    case 201:
		    case 102:
			set_crxps_and_edges_201(ctri,cf,sides);
			break;
		    default:
			printf("Unexpected case in set_polygons_in_cell().\n");
			clean_up(ERROR);
		}
	    }
	}

	//debugdan	FIXME
	/*
	if (c->icrds[0] == 9 && c->icrds[1] == 4 && c->icrds[2] == 20)
	{
	    DEBUGDAN = true;
	}
	else
	{
	    DEBUGDAN = false;
	}
	*/
	//debugdan	FIXME

	for (ctri = c->ctris; ctri; ctri = ctri->next)
	{
	    construct_polygon_on_tri(ctri,c);
	    if (ctri->polyg)
	    {
		add_new_polyg(ctri->polyg,&c->ctri_polygs);
	    }
	}

	//debugdan	FIXME
	if (DEBUGDAN)
	{
	    printf("\nDEBUGDAN:\n");
	    CPOLYGON *debug_polyg = c->ctri_polygs;
	    while (debug_polyg)
	    {
		print_polygon(debug_polyg);
		debug_polyg = debug_polyg->next;
	    }
	}
	//debugdan	FIXME

	if (c->ctri_polygs)
	    construct_polygons_on_cfs(c);

	return;
}

bool tri_in_cell(
	CTRI	*ctri,
	CELL	*c)
{
	int i, j;
	CPOINT *p;
	double tol = 1e-12;

	for (i = 0; i < 3; i++)
	{
	    if (!point_in_cell(&(ctri->pts[i]),c))
		return NO;
	}

	return YES;
}

bool point_in_cell(
	CPOINT	*p,
	CELL	*c)
{
    int i;
    double tol = 1e-12;

    for (i = 0; i < 3; i++)
    {
	if (p->crds[i] < c->celll[i]-tol || p->crds[i] > c->cellu[i]+tol)
	    return NO;
    }
    return YES;
}

/*inside: 001; on: 010; outside: 100*/
void ctri_inter_with_cf(
	CTRI	*ctri,
	CFACE	*cf,
	int	*sides)
{
	int i, dir;
	CPOINT *p;
	double tol = 1e-12;

	dir = cf->dir;
	for (i = 0; i < 3; i++)
	{
	    p = &(ctri->pts[i]);
	    if (fabs(p->crds[dir]-cf->l[dir]) <= tol)
		sides[i] = 10;	//on cell face
	    else if (p->crds[dir] > cf->l[dir])
	    {
		if (cf->side == -1)
		    sides[i] = 1;	//inside
		else
		    sides[i] = 100;	//outside
	    }
	    else
	    {
		if (cf->side == -1)
		    sides[i] = 100;	//outside
		else
		    sides[i] = 1;	//inside
	    }
	}
}

void set_crxps_and_edges_201(
	CTRI	*ctri,
	CFACE	*cf,
	int	*sides)
{
	int i, count, ncrxps, ei[2];
	bool pocf[2], found;
	CPOINT crxp[2], crxp2, crxp3;
	CEDGE *pedge;

	count = 0;
	for (i = 0; i < 3; i++)
	{
	    if (sides[i] == sides[(i+1)%3])
		continue;
	    
	    find_crxp_3d1(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),cf,&(crxp[count]));
	    ei[count] = i;
	    pocf[count] = point_on_cell_face(&(crxp[count]),cf);
	    count++;
	}
	if (count != 2)
	{
	    printf("ERROR in set_crxps_and_edges_201().\n");
	    clean_up(ERROR);
	}

	found = false;
	if (pocf[0] && pocf[1])
	{
	    FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
	    copy_cpt_to_cpt(&(crxp[0]),&(pedge->endp[0]));
	    add_crxp_on_tri_edge(&crxp[0],&(ctri->edges[ei[0]]));
	    copy_cpt_to_cpt(&(crxp[1]),&(pedge->endp[1]));
	    add_crxp_on_tri_edge(&crxp[1],&(ctri->edges[ei[1]]));
	    found = YES;
	}
	else if (pocf[0] && !pocf[1])
	{
	    ncrxps = find_crxps1(&(crxp[0]),&(crxp[1]),cf,&crxp2,&crxp3);
	    add_crxp_on_tri_edge(&crxp[0],&(ctri->edges[ei[0]]));
	    if (ncrxps == 1)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		copy_cpt_to_cpt(&(crxp[0]),&(pedge->endp[0]));
		copy_cpt_to_cpt(&crxp2,&(pedge->endp[1]));
		found = YES;
	    }
	}
	else if (!pocf[0] && pocf[1])
	{
	    ncrxps = find_crxps1(&(crxp[0]),&(crxp[1]),cf,&crxp2,&crxp3);
	    add_crxp_on_tri_edge(&crxp[1],&(ctri->edges[ei[1]]));
	    if (ncrxps == 1)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		copy_cpt_to_cpt(&crxp2,&(pedge->endp[0]));
		copy_cpt_to_cpt(&(crxp[1]),&(pedge->endp[1]));
		add_crxp_on_tri_edge(&crxp[1],&(ctri->edges[ei[1]]));
		found = YES;
	    }
	}
	else
	{
	    ncrxps = find_crxps1(&(crxp[0]),&(crxp[1]),cf,&crxp2,&crxp3);
	    if (ncrxps == 2)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		copy_cpt_to_cpt(&crxp2,&(pedge->endp[0]));
		copy_cpt_to_cpt(&crxp3,&(pedge->endp[1]));
		found = YES;
	    }
	}

	if (found)
	{
	    add_edge(pedge,&(cf->undirected_edges));
	    add_edge(pedge,&(ctri->polyg->undir_edges));
	    free(pedge);
	}
}

void set_crxps_and_edges_111(
	CTRI	*ctri,
	CFACE	*cf,
	int	*sides)
{
	int i, ncrxps;
	bool pocf1, pocf2, found;
	CPOINT *p1, *p2;
	CPOINT crxp1, crxp2, crxp3;
	CEDGE *pedge;

	for (i = 0; i < 3; i++)
	{
	    if (sides[i] == 10)
	    {
		p1 = &(ctri->pts[i]);
		break;
	    }
	}

	found = NO;
	find_crxp_3d1(&(ctri->pts[(i+1)%3]),&(ctri->pts[(i+2)%3]),cf,&crxp1);
	pocf1 = point_on_cell_face(&(ctri->pts[i]),cf);
	pocf2 = point_on_cell_face(&crxp1,cf);
	if (pocf1 && pocf2)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
	    copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
	    copy_cpt_to_cpt(&crxp1,&(pedge->endp[1]));
	    add_crxp_on_tri_edge(&crxp1,&(ctri->edges[(i+1)%3]));
	    found = YES;
	}
	else if (pocf1 && !pocf2)
	{
	    ncrxps = find_crxps1(&(ctri->pts[i]),&crxp1,cf,&crxp2,&crxp3);
	    if (ncrxps == 1)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
		copy_cpt_to_cpt(&crxp2,&(pedge->endp[1]));
		found = YES;
	    }
	}
	else if (!pocf1 && pocf2)
	{
	    ncrxps = find_crxps1(&(ctri->pts[i]),&crxp1,cf,&crxp2,&crxp3);
	    if (ncrxps == 1)	//double check	FIXME
	    {
		FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		copy_cpt_to_cpt(&crxp2,&(pedge->endp[0]));
		copy_cpt_to_cpt(&crxp1,&(pedge->endp[1]));
		add_crxp_on_tri_edge(&crxp1,&(ctri->edges[(i+1)%3]));
		found = YES;
	    }
	}
	else
	{
	    ncrxps = find_crxps1(&(ctri->pts[i]),&crxp1,cf,&crxp2,&crxp3);
	    if (ncrxps == 2)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		copy_cpt_to_cpt(&crxp2,&(pedge->endp[0]));
		copy_cpt_to_cpt(&crxp3,&(pedge->endp[1]));
		found = YES;
	    }
	}

	if (found)
	{
	    add_edge(pedge,&(cf->undirected_edges));
	    add_edge(pedge,&(ctri->polyg->undir_edges));
	    free(pedge);
	}

	return;
}

void set_crxps_and_edges_30(
	CTRI	*ctri,
	CFACE	*cf,
	int	*sides)
{
	int i, j, ncrxps;
	bool pocf1, pocf2, pit1, pit2, found;
	CPOINT crxp1, crxp2, p1, p2;
	CEDGE *pedge;

	/*cut tri edges*/
	for (i = 0; i < 3; i++)
	{
	    found = NO;
	    pocf1 = point_on_cell_face(&(ctri->pts[i]),cf);
	    pocf2 = point_on_cell_face(&(ctri->pts[(i+1)%3]),cf);
	    if (pocf1 && pocf2)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
		copy_cpt_to_cpt(&(ctri->pts[(i+1)%3]),&(pedge->endp[1]));
		found = YES;
	    }
	    else if (pocf1 && !pocf2)
	    {
		ncrxps = find_crxps1(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),cf,&crxp1,&crxp2);
		if (ncrxps == 1)
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		    copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
		    copy_cpt_to_cpt(&crxp1,&(pedge->endp[1]));
		    add_crxp_on_tri_edge(&crxp1,&(ctri->edges[i]));
		    found = YES;
		}
	    }
	    else if (!pocf1 && pocf2)
	    {
		ncrxps = find_crxps1(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),cf,&crxp1,&crxp2);
		if (ncrxps == 1)
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		    copy_cpt_to_cpt(&crxp1,&(pedge->endp[0]));
		    copy_cpt_to_cpt(&(ctri->pts[(i+1)%3]),&(pedge->endp[1]));
		    add_crxp_on_tri_edge(&crxp1,&(ctri->edges[i]));
		    found = YES;
		}
	    }
	    else if (!pocf1 && !pocf2)
	    {
		ncrxps = find_crxps1(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),cf,&crxp1,&crxp2);
		if (ncrxps == 2)
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		    copy_cpt_to_cpt(&crxp1,&(pedge->endp[0]));
		    copy_cpt_to_cpt(&crxp2,&(pedge->endp[1]));
		    add_crxp_on_tri_edge(&crxp1,&(ctri->edges[i]));
		    add_crxp_on_tri_edge(&crxp2,&(ctri->edges[i]));
		    found = YES;
		}
	    }

	    if (found)
	    {
		add_edge(pedge,&(cf->undirected_edges));
		add_edge(pedge,&(ctri->polyg->undir_edges));
		free(pedge);
	    }
	}

	/*cut cell edges*/
	for (i = 0; i < 4; i++)
	{
	    found = NO;
	    pit1 = point_in_tri_2d(&(cf->pts[i]),ctri,cf->dir);
	    pit2 = point_in_tri_2d(&(cf->pts[(i+1)%4]),ctri,cf->dir);
	    if (pit1 && pit2)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		copy_cpt_to_cpt(&(cf->pts[i]),&(pedge->endp[0]));
		copy_cpt_to_cpt(&(cf->pts[(i+1)%4]),&(pedge->endp[1]));
		found = YES;
	    }
	    else if (pit1 && !pit2)
	    {
		ncrxps = find_crxps2(cf,i,ctri,&crxp1,&crxp2);
		if (ncrxps == 1)
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		    copy_cpt_to_cpt(&(cf->pts[i]),&(pedge->endp[0]));
		    copy_cpt_to_cpt(&crxp1,&(pedge->endp[1]));
		    found = YES;
		}
	    }
	    else if (!pit1 && pit2)
	    {
		ncrxps = find_crxps2(cf,i,ctri,&crxp1,&crxp2);
		if (ncrxps == 1)
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		    copy_cpt_to_cpt(&crxp1,&(pedge->endp[0]));
		    copy_cpt_to_cpt(&(cf->pts[(i+1)%4]),&(pedge->endp[1]));
		    found = YES;
		}
	    }
	    else if (!pit1 && !pit2)
	    {
		ncrxps = find_crxps2(cf,i,ctri,&crxp1,&crxp2);
		if (ncrxps == 2)
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		    copy_cpt_to_cpt(&crxp1,&(pedge->endp[0]));
		    copy_cpt_to_cpt(&crxp2,&(pedge->endp[1]));
		    found = YES;
		}
	    }

	    if (found)
	    {
		if (!edge_exists_in_edge_list(pedge,cf->undirected_edges))
		{
		    add_edge(pedge,&(cf->undirected_edges));
		    add_edge(pedge,&(ctri->polyg->undir_edges));
		}
		free(pedge);
	    }
	}

	return;
}

/*nothing on tri; a new line on cell face*/
void set_crxps_and_edges_21(
	CTRI	*ctri,
	CFACE	*cf,
	int	*sides)
{
	int i, ncrxps;
	CPOINT crxp1, crxp2;
	CEDGE *pedge;
	bool pocf1, pocf2;
	bool found = NO;

	for (i = 0; i < 3; i++)
	{
	    if (sides[i] == 10 && sides[(i+1)%3] == 10)
	    {
//		FT_ScalarMemoryAlloc((POINTER*)&pedge,sizeof(CEDGE));
		pocf1 = point_on_cell_face(&(ctri->pts[i]),cf);
		pocf2 = point_on_cell_face(&(ctri->pts[(i+1)%3]),cf);
		if (pocf1 && pocf2)
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
		    copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
		    copy_cpt_to_cpt(&(ctri->pts[(i+1)%3]),&(pedge->endp[1]));
		    found = YES;
		}
		else if (pocf1 && !pocf2)
		{
		    ncrxps = find_crxps1(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),cf,&crxp1,&crxp2);
		    if (ncrxps == 1)
		    {
			FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
			copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
			copy_cpt_to_cpt(&crxp1,&(pedge->endp[1]));
			add_crxp_on_tri_edge(&crxp1,&(ctri->edges[i]));
			found = YES;
		    }
		}
		else if (!pocf1 && pocf2)
		{
		    ncrxps = find_crxps1(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),cf,&crxp1,&crxp2);
		    if (ncrxps == 1)
		    {
			FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
			copy_cpt_to_cpt(&crxp1,&(pedge->endp[0]));
			copy_cpt_to_cpt(&(ctri->pts[(i+1)%3]),&(pedge->endp[1]));
			add_crxp_on_tri_edge(&crxp1,&(ctri->edges[i]));
			found = YES;
		    }
		}
		else if (!pocf1 && !pocf2)
		{
		    ncrxps = find_crxps1(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),cf,&crxp1,&crxp2);
		    if (ncrxps == 2)
		    {
			FT_ScalarMemoryAlloc((POINTER*)&(pedge),sizeof(CEDGE));
			copy_cpt_to_cpt(&crxp1,&(pedge->endp[0]));
			copy_cpt_to_cpt(&crxp2,&(pedge->endp[1]));
			add_crxp_on_tri_edge(&crxp1,&(ctri->edges[i]));
			add_crxp_on_tri_edge(&crxp2,&(ctri->edges[i]));
			found = YES;
		    }
		}

		if (found)
		{
		    add_edge(pedge,&(cf->undirected_edges));
		    add_edge(pedge,&(ctri->polyg->undir_edges));
		    free(pedge);
		}
	    }
	}

	return;
}

void add_crxp_on_tri_edge(
	CPOINT	*crxp,
	CEDGE	*edge)
{
	CPOINT *newcrxp, *tmpp;

	if (same_cpt(crxp,&(edge->endp[0])) ||
	    same_cpt(crxp,&(edge->endp[1])))
	{
	    return;
	}
	tmpp = edge->crxps;
	while (tmpp)
	{
	    if (same_cpt(crxp,tmpp))
		return;
	    tmpp = tmpp->next;
	}
	FT_ScalarMemoryAlloc((POINTER*)&(newcrxp),sizeof(CPOINT));
	copy_cpt_to_cpt(crxp,newcrxp);
	newcrxp->next = edge->crxps;
	edge->crxps = newcrxp;
	edge->ncrxps++;

	return;
}

void add_edge(
	CEDGE	*edge,
	CEDGE	**edgelist)
{
	CEDGE *new_edge;
	CEDGE *tmp_edge;

	//check duplicated edges.
	tmp_edge = *edgelist;
	while (tmp_edge)
	{
	    if (same_undir_edge(tmp_edge,edge))
		return;
	    tmp_edge = tmp_edge->next;
	}

	FT_ScalarMemoryAlloc((POINTER*)&(new_edge),sizeof(CEDGE));
	copy_cpt_to_cpt(&(edge->endp[0]),&(new_edge->endp[0]));
	copy_cpt_to_cpt(&(edge->endp[1]),&(new_edge->endp[1]));
	new_edge->mark = 0;

	new_edge->next = *edgelist;
	*edgelist = new_edge;

	return;
}

bool edge_exists_in_edge_list(
	CEDGE	*edge,
	CEDGE	*edgelist)
{
	CEDGE *tmpe;

	tmpe = edgelist;
	while (tmpe)
	{
	    if (same_undir_edge(tmpe,edge))
		return YES;
	    tmpe = tmpe->next;
	}
	return NO;
}

bool point_on_cell_face(
	CPOINT	*cpt,
	CFACE	*cf)
{
	int i;
	double tol = 1e-12;

	for (i = 0; i < 3; i++)
	{
	    if (i == cf->dir)
		continue;
	    if (cpt->crds[i] < cf->l[i]-tol || cpt->crds[i] > cf->u[i]+tol)
		return NO;
	}

	return YES;
}

bool point_in_cell_face(
	CPOINT	*cpt,
	CFACE	*cf)
{
	int i;
	double tol = 1e-12;

	for (i = 0; i < 3; i++)
	{
	    if (i == cf->dir)
		continue;
	    if (cpt->crds[i] < cf->l[i]+tol || cpt->crds[i] > cf->u[i]-tol)
		return NO;
	}

	return YES;
}

bool point_in_tri_2d(
	CPOINT	*cpt,
	CTRI	*ctri,
	int	dir)
{
	int i, count;
	double a, b;
	double v1[2], v2[2], v3[2];
	double tol = 1e-12;

	count = 0;
	for (i = 0; i < 3; i++)
	{
	    if (i == dir)
		continue;
	    v1[count] = ctri->pts[1].crds[i] - ctri->pts[0].crds[i];
	    v2[count] = ctri->pts[2].crds[i] - ctri->pts[0].crds[i];
	    v3[count] = cpt->crds[i] - ctri->pts[0].crds[i];
	    count++;
	}

	/*solve a*v1 + b*v2 = v3*/
	a = det2(v3,v2)/det2(v1,v2);
	b = det2(v1,v3)/det2(v1,v2);

	if (a > -tol && a < 1+tol &&
	    b > -tol && b < 1+tol &&
	    a+b < 1+tol)
	    return YES;
	else
	    return NO;
}

double det2(
	double	*a,
	double	*b)
{
	return a[0]*b[1] - a[1]*b[0];
}

/*Intersections of a line and a cell face in 2D*/
int find_crxps1(
	CPOINT	*cpt1,
	CPOINT	*cpt2,
	CFACE	*cf,
	CPOINT	*crxp1,
	CPOINT	*crxp2)
{
	int i, j, dir, ncrxps;
	double dist1, dist2;
	CPOINT crxp;
	bool found;
	double tol = 1e-12;

	ncrxps = 0;
	dir = cf->dir;
	for (i = 0; i < 3; i++)
	{
	    if (i == dir)
		continue;

	    //intersection with the lower edge
	    found = NO;
	    dist1 = cpt1->crds[i] - cf->l[i];
	    dist2 = cpt2->crds[i] - cf->l[i];
	    //if (dist1*dist2 < -tol)
	    if ((dist1 > tol && dist2 < -tol) ||
		(dist1 < -tol && dist2 > tol))
	    {
		//set the coordinates of the crossing point
		for (j = 0; j < 3; j++)
		{
		    if (j == dir)
			crxp.crds[j] = cf->l[j];
		    else if (j == i)
			crxp.crds[j] = cf->l[j];
		    else
		    {
			dist1 = fabs(dist1);
			dist2 = fabs(dist2);
			crxp.crds[j] = dist2/(dist1+dist2)*cpt1->crds[j] +
			    	       dist1/(dist1+dist2)*cpt2->crds[j];
			if (crxp.crds[j] > cf->l[j]-tol && crxp.crds[j] < cf->u[j]+tol)
			    found = YES;
			else
			    break;
		    }
		}
		if (found)
		{
		    if (ncrxps == 0)
		    {
			ncrxps++;
			copy_cpt_to_cpt(&crxp,crxp1);
		    }
		    else if (ncrxps == 1)
		    {
			//check if it is a duplicated crxp
			if (!same_cpt(crxp1,&crxp))
			{
			    ncrxps++;
			    copy_cpt_to_cpt(&crxp,crxp2);
			}
		    }
		    else
		    {
			if (!same_cpt(crxp1,&crxp) && !same_cpt(crxp2,&crxp))
			{
			    printf("ERROR: too many crxps in find_crxps1().\n");
			    clean_up(ERROR);
			}
		    }
		}
	    }

	    //intersection with the upper edge
	    found = NO;
	    dist1 = cpt1->crds[i] - cf->u[i];
	    dist2 = cpt2->crds[i] - cf->u[i];
	    //if (dist1*dist2 < -tol)
	    if ((dist1 > tol && dist2 < -tol) ||
		(dist1 < -tol && dist2 > tol))
	    {
		//set the coordinates of the crossing point
		for (j = 0; j < 3; j++)
		{
		    if (j == dir)
			crxp.crds[j] = cf->l[j];
		    else if (j == i)
			crxp.crds[j] = cf->u[j];	//FIXME
		    else
		    {
			dist1 = fabs(dist1);
			dist2 = fabs(dist2);
			crxp.crds[j] = dist2/(dist1+dist2)*cpt1->crds[j] +
			    	       dist1/(dist1+dist2)*cpt2->crds[j];
			if (crxp.crds[j] > cf->l[j]-tol && crxp.crds[j] < cf->u[j]+tol)
			    found = YES;
			else
			    break;
		    }
		}
		if (found)
		{
		    if (ncrxps == 0)
		    {
			ncrxps++;
			copy_cpt_to_cpt(&crxp,crxp1);
		    }
		    else if (ncrxps == 1)
		    {
			//check if it is a duplicated crxp
			if (!same_cpt(crxp1,&crxp))
			{
			    ncrxps++;
			    copy_cpt_to_cpt(&crxp,crxp2);
			}
		    }
		    else
		    {
			if (!same_cpt(crxp1,&crxp) && !same_cpt(crxp2,&crxp))
			{
			    printf("ERROR: too many crxps in find_crxps1().\n");
			    clean_up(ERROR);
			}
		    }
		}
	    }
	}

	return ncrxps;
}

/*Intersections of a line and a tri in 2D*/
/*crxps are different than vertices?	TODO*/
int find_crxps2(
	CFACE	*cf,
	int	index,
	CTRI	*ctri,
	CPOINT	*crxp1,
	CPOINT	*crxp2)
{
	int i, j, edir, ncrxps;
	double f1, f2;
	CPOINT crxp;
	CPOINT *p1, *p2, *tp1, *tp2;
	bool found;
	double tol = 1e-12;

	p1 = &(cf->pts[index]);
	p2 = &(cf->pts[(index+1)%4]);

	ncrxps = 0;
	for (i = 0; i < 3; i++)
	{
	    if (fabs(p1->crds[i] - p2->crds[i]) > tol)
	    {
		edir = i;
		break;
	    }
	}

	for (i = 0; i < 3; i++)
	{
	    found = NO;
	    tp1 = &(ctri->pts[i]);
	    tp2 = &(ctri->pts[(i+1)%3]);
	    for (j = 0; j < 3; j++)
	    {
		if (j == edir)
		    continue;
		crxp.crds[j] = p1->crds[j];
		if (j != cf->dir)
		{
		    f1 = tp1->crds[j] - p1->crds[j];
		    f2 = p1->crds[j] - tp2->crds[j];
		}
	    }
	    //if (f1*f2 < tol)
	    if (!((f1 > tol && f2 > tol) || (f1 < -tol && f2 < -tol)))
		continue;
	    crxp.crds[edir] = f2/(f1+f2)*tp1->crds[edir] + f1/(f1+f2)*tp2->crds[edir];
	    //if ((p1->crds[edir]-crxp.crds[edir])*(crxp.crds[edir]-p2->crds[edir]) > tol)
	    if ((p1->crds[edir]-crxp.crds[edir] > tol && p2->crds[edir]-crxp.crds[edir] < -tol) ||
	        (p1->crds[edir]-crxp.crds[edir] < -tol && p2->crds[edir]-crxp.crds[edir] > tol))
		found = YES;
	    if (found)
	    {
		if (ncrxps == 0)
		{
		    ncrxps++;
		    copy_cpt_to_cpt(&crxp,crxp1);
		}
		else if (ncrxps == 1)
		{
		    //check if it is a duplicated crxp
		    if (!same_cpt(crxp1,&crxp))
		    {
			ncrxps++;
			copy_cpt_to_cpt(&crxp,crxp2);
		    }
		}
		else
		{
		    if (!same_cpt(crxp1,&crxp) && !same_cpt(crxp2,&crxp))
		    {
			printf("ERROR: too many crxps in find_crxps2().\n");
			clean_up(ERROR);
		    }
		}
	    }
	}
	return ncrxps;
}

bool find_crxp_3d1(
	CPOINT	*p1,
	CPOINT	*p2,
	CFACE	*cf,
	CPOINT	*crxp)
{
	int i, dir;
	double f1, f2;
	double tol = 1e-12;

	dir = cf->dir;
	f1 = p1->crds[dir] - cf->l[dir];
	f2 = cf->l[dir] - p2->crds[dir];

	//if (f1*f2 < tol)
	if (fabs(f1) < tol || fabs(f2) < tol ||
	    (f1 > tol && f2 < -tol) ||
	    (f1 < -tol && f2 > tol))
	{
	    printf("ERROR in find_crxp_3d1.\n");
	    return NO;
	}

	for (i = 0; i < 3; i++)
	{
	    if (i == dir)
	    {
		crxp->crds[i] = cf->l[i];
		continue;
	    }
	    crxp->crds[i] = f2/(f1+f2)*p1->crds[i] + f1/(f1+f2)*p2->crds[i];
	}
	return YES;
}

bool same_cpt_double_tol(
	CPOINT	*p1,
	CPOINT	*p2)
{
	int i;
	double tol = 2e-12;

	for (i = 0; i < 3; i++)
	    if (fabs(p1->crds[i] - p2->crds[i]) > tol)
		return NO;

	return YES;
}

bool same_cpt(
	CPOINT	*p1,
	CPOINT	*p2)
{
	int i;
	double tol = 1e-12;

	for (i = 0; i < 3; i++)
	    if (fabs(p1->crds[i] - p2->crds[i]) > tol)
		return NO;

	return YES;
}

bool same_undir_edge(
	CEDGE	*edge1,
	CEDGE	*edge2)
{
	if ((same_cpt(&(edge1->endp[0]),&(edge2->endp[0])) &&
	    same_cpt(&(edge1->endp[1]),&(edge2->endp[1]))) ||
	    (same_cpt(&(edge1->endp[0]),&(edge2->endp[1])) &&
	    same_cpt(&(edge1->endp[1]),&(edge2->endp[0]))))
	    return YES;

	return NO;
}

void print_point(CPOINT *p)
{
	//printf("(%lf, %lf, %lf)", p->crds[0], p->crds[1], p->crds[2]);
	printf("{%.12g, %.12g, %.12g}", p->crds[0], p->crds[1], p->crds[2]);
}

void print_edge(CEDGE *edge)
{
	CPOINT *p1, *p2;

	p1 = &(edge->endp[0]);
	p2 = &(edge->endp[1]);

	printf("Edge: ");
	print_point(p1);
	printf(", ");
	print_point(p2);
	printf(".\n");
}

void print_ctri(CTRI *ctri)
{
	printf("Tri: ");
	print_point(&(ctri->pts[0]));
	printf(", ");
	print_point(&(ctri->pts[1]));
	printf(", ");
	print_point(&(ctri->pts[2]));
	printf(".\n");
}

/*
void print_edge_list(CEDGE *edge_list)
{
    CPOINT *p1, *p2;
    CEDGE *edge = edge_list;

    while (edge)
    {
	print_edge(edge);
	edge = edge->next;
    }
}
*/
void construct_polygon_on_tri(
	CTRI	*ctri,
	CELL	*cell)
{
	int i;
	CPOINT *p;

	//tri is a polygon
	if (tri_in_cell(ctri,cell))
	{
	    for (i = 2; i >=0; i--)
	    {
		add_polyg_vertex(&(ctri->pts[i]),ctri->polyg);
	    }
	    return;
	}

	//no polygon on tri
	if (ctri->polyg->undir_edges == NULL)
	{
	    free(ctri->polyg);
	    ctri->polyg = NULL;
	    return;
	}

	//inside and outside tris have been excluded? needs double check	TODO
	//polygon on tri
	//complete ctri->polyg->undir_edges
	complete_polyg_undir_edges(ctri,cell);

	set_directed_polyg(ctri->polyg,ctri->nor);

	return;
}

void add_polyg_vertex(
	CPOINT		*p,
	CPOLYGON	*pg)
{
	CPOINT *newp;

	FT_ScalarMemoryAlloc((POINTER*)&(newp),sizeof(CPOINT));
	copy_cpt_to_cpt(p,newp);
	newp->next = pg->vertices;
	pg->vertices = newp;

	return;
}

void print_polygon(CPOLYGON *pg)
{
	CPOINT *p;
	printf("Polygon:\n");
	p = pg->vertices;
	while(p)
	{
	    print_point(p);
	    printf("\n");
	    p = p->next;
	}
}

void complete_polyg_undir_edges(
	CTRI	*ctri,
	CELL	*cell)
{
	int i, ncrxps;
	PIC pic[3];
	CPOINT crxp0, crxp1;

	set_pic(ctri,cell,pic);

	for (i = 0; i < 3; i++)
	{
	    if ((pic[i] == INC || pic[i] == ONCF) &&
		(pic[(i+1)%3] == INC || pic[(i+1)%3] ==ONCF))
	    {
		add_new_edge(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),&(ctri->polyg->undir_edges));
	    }
	    else if (pic[i] == INC && pic[(i+1)%3] == OUTC)
	    {
		if (ctri->edges[i].ncrxps != 1)
		{
		    printf("ERROR in complete_polyg_undir_edges().\n");
		    clean_up(ERROR);
		}
		add_new_edge(&(ctri->pts[i]),ctri->edges[i].crxps,&(ctri->polyg->undir_edges));
	    }
	    else if (pic[i] == OUTC && pic[(i+1)%3] == INC)
	    {
		if (ctri->edges[i].ncrxps != 1)
		{
		    printf("ERROR in complete_polyg_undir_edges().\n");
		    clean_up(ERROR);
		}
		add_new_edge(ctri->edges[i].crxps,&(ctri->pts[(i+1)%3]),&(ctri->polyg->undir_edges));
	    }
	    else if (pic[i] == ONCF && pic[(i+1)%3] == OUTC)
	    {
		if (ctri->edges[i].ncrxps == 1)
		{
		    add_new_edge(&(ctri->pts[i]),ctri->edges[i].crxps,&(ctri->polyg->undir_edges));
		}
	    }
	    else if (pic[i] == OUTC && pic[(i+1)%3] == ONCF)
	    {
		if (ctri->edges[i].ncrxps == 1)
		{
		    add_new_edge(ctri->edges[i].crxps,&(ctri->pts[(i+1)%3]),&(ctri->polyg->undir_edges));
		}
	    }
	    else if (pic[i] == OUTC && pic[(i+1)%3] == OUTC)
	    {
		if (ctri->edges[i].ncrxps == 2)
		{
		    add_new_edge(ctri->edges[i].crxps,ctri->edges[i].crxps->next,&(ctri->polyg->undir_edges));
		}
	    }
	}

	return;
}

void set_pic(
	CTRI	*ctri,
	CELL	*cell,
	PIC	*pic)
{
	int i, j;
	CPOINT	*p;
	double tol = 1e-12;

	for (i = 0; i < 3; i++)
	{
	    p = &(ctri->pts[i]);
	    pic[i] = INC;	//init
	    for (j = 0; j < 3; j++)
	    {
		if (p->crds[j] < cell->celll[j]-tol ||
		    p->crds[j] > cell->cellu[j]+tol)
		{
		    pic[i] = OUTC;
		    break;
		}
	    }
	    if (pic[i] == OUTC)
		continue;
	    else
	    {
		for (j = 0; j < 3; j++)
		{
		    if (fabs(p->crds[j] - cell->celll[j]) < tol ||
			fabs(p->crds[j] - cell->cellu[j]) < tol)
		    {
			pic[i] = ONCF;
			break;
		    }
		}
	    }
	}

	return;
}

void add_new_edge(
	CPOINT	*p0,
	CPOINT	*p1,
	CEDGE	**edgelist)
{
	CEDGE *newedge;
/*
	if (*edgelist == NULL)
	{
	    printf("ERROR in add_new_edge().\n");
	    clean_up(ERROR);
	}

	if (same_cpt(p0,p1))
	    return;
*/
	FT_ScalarMemoryAlloc((POINTER*)&(newedge),sizeof(CEDGE));
	copy_cpt_to_cpt(p0,&(newedge->endp[0]));
	copy_cpt_to_cpt(p1,&(newedge->endp[1]));
	add_edge(newedge,edgelist);

	free(newedge);
	return;
}

void add_new_dir_edge(
	CPOINT	*p0,
	CPOINT	*p1,
	CEDGE	**edgelist)
{
	CEDGE *newedge;

	FT_ScalarMemoryAlloc((POINTER*)&(newedge),sizeof(CEDGE));
	if (p0->flag == 0)
	{
	    copy_cpt_to_cpt(p1,&(newedge->endp[0]));
	    copy_cpt_to_cpt(p0,&(newedge->endp[1]));
	    //p1->flag = 1;
	}
	else if (p0->flag == 1)
	{
	    copy_cpt_to_cpt(p0,&(newedge->endp[0]));
	    copy_cpt_to_cpt(p1,&(newedge->endp[1]));
	    //p1->flag = 0;
	}
	else
	{
	    printf("ERROR in add_new_dir_edge(): flag for the first point is not set.\n");
	    clean_up(ERROR);
	}
	add_edge(newedge,edgelist);

	free(newedge);
}

void add_new_polyg(
	CPOLYGON	*polyg,
	CPOLYGON	**polyglist)
{
	CPOINT *p;
	CPOLYGON *newpolyg;

	p = polyg->vertices;
	if (p == NULL)
	    return;
	FT_ScalarMemoryAlloc((POINTER*)&(newpolyg),sizeof(CPOLYGON));
	init_new_polyg(newpolyg);
	//newpolyg->vertices = NULL;
	while (p)
	{
	    add_polyg_vertex(p,newpolyg);
	    p = p->next;
	}
	reverse_polyg(newpolyg);
	newpolyg->oncf = polyg->oncf;

	newpolyg->next = *polyglist;
	*polyglist = newpolyg;

	return;
}

void reverse_polyg(CPOLYGON *polyg)
{
	CPOINT *p, *startp, *prevp;

	startp = polyg->vertices;
	p = polyg->vertices;
	prevp = p;
	if (!p || !(p->next))
	{
	    printf("ERROR in reverse_polyg: incorrect polygon.\n");
	    clean_up(ERROR);
	}
	while (p->next)
	{
	    prevp = p;
	    p = p->next;
	}
	polyg->vertices = p;
	p->next = prevp;
	prevp->next = NULL;
	while (startp != prevp)
	{
	    p = startp;
	    while (p->next != prevp)
		p = p->next;
	    prevp->next = p;
	    p->next = NULL;
	    prevp = p;
	}

	return;
}

void set_directed_polyg(
	CPOLYGON	*polyg,
	double		*nor)
{
	int i;
	double pnor[3];
	CPOINT p0, p1, p2, tmpp;
	CEDGE *edge0;

	pop_edge(&(polyg->undir_edges),&edge0);
	copy_cpt_to_cpt(&(edge0->endp[0]),&p0);
	copy_cpt_to_cpt(&(edge0->endp[1]),&p1);
	if (!find_next_endp(p1,&(polyg->undir_edges),&p2))
	{
	    printf("Potential error in set_directed_polyg().\n");
	    free(edge0);
	    return;
	    /*
	    if (!find_next_endp(p0,&(polyg->undir_edges),&p2))
	    {
		printf("Potential error in set_directed_polyg().\n");
		free(edge0);
		return;
	    }
	    if (!set_nor(&p2,&p0,&p1,pnor))
	    {
		printf("ERROR in set_directed_polyg(): points are on the same line.\n");
		clean_up(ERROR);
	    }
	    if (same_nor_dir(pnor,nor))
	    {
		add_polyg_vertex(&p1,polyg);
		add_polyg_vertex(&p0,polyg);
		add_polyg_vertex(&p2,polyg);
	    }
	    else
	    {
		add_polyg_vertex(&p2,polyg);
		add_polyg_vertex(&p0,polyg);
		add_polyg_vertex(&p1,polyg);
	    }
	    return;
	    */
	}

	//debugdan	FIXME
	if (DEBUGDAN)
	{
	    printf("set_directed_polyg():\n");
	    printf("p0, p1, p2:\n");
	    print_point(&p0);
	    printf("\n");
	    print_point(&p1);
	    printf("\n");
	    print_point(&p2);
	    printf("\n");
	}
	//debugdan	FIXME

	if (!set_nor(&p0,&p1,&p2,pnor))
	{
	    printf("ERROR in set_directed_polyg(): points are on the same line.\n");
	    clean_up(ERROR);
	}

	if (same_nor_dir(pnor,nor))
	{
	    add_polyg_vertex(&p2,polyg);
	    add_polyg_vertex(&p1,polyg);
	    add_polyg_vertex(&p0,polyg);
	}
	else
	{
	    add_polyg_vertex(&p0,polyg);
	    add_polyg_vertex(&p1,polyg);
	    add_polyg_vertex(&p2,polyg);
	    tmpp = p0;
	    p0 = p2;
	    p2 = tmpp;
	}

	if (!find_next_endp(p0,&(polyg->undir_edges),&p1))
	    clean_up(ERROR);
	while (!same_cpt(&p1,&p2))
	{
	    add_polyg_vertex(&p1,polyg);
	    p0 = p1;
	    if (!find_next_endp(p0,&(polyg->undir_edges),&p1))
		clean_up(ERROR);
	}

	//debugdan	FIXME
	if (DEBUGDAN)
	{
	    printf("set_directed_polyg():\n");
	    print_polygon(polyg);
	    printf("pnor: %lf %lf, %lf.\n", pnor[0], pnor[1], pnor[2]);
	    printf("nor: %lf %lf, %lf.\n", nor[0], nor[1], nor[2]);
	}
	//debugdan	FIXME

	free(edge0);
	return;
}

void pop_edge(
	CEDGE	**edgelist,
	CEDGE	**edge)
{
	if (edgelist == NULL)
	{
	    printf("Empty edgelist in pop_edge().\n");
	    clean_up(ERROR);
	}

	*edge = *edgelist;
	*edgelist = (*edge)->next;

	return;
}

bool find_next_endp(
	CPOINT	p0,
	CEDGE	**edgelist,
	CPOINT	*p1)
{
	CEDGE *edge, *prev_edge;

	if (edgelist == NULL)
	{
	    printf("WARNING: Empty edgelist in find_next_endp().\n");
	    //clean_up(ERROR);
	    return NO;
	}

	prev_edge = NULL;
	edge = *edgelist;
	while (edge)
	{
	    if (same_cpt(&p0,&(edge->endp[0])))
	    {
		copy_cpt_to_cpt(&(edge->endp[1]),p1);
		if (prev_edge == NULL)
		{
		    *edgelist = edge->next;
		    free(edge);
		}
		else
		{
		    prev_edge->next = edge->next;
		    free(edge);
		}
		return YES;
	    }
	    else if (same_cpt(&p0,&(edge->endp[1])))
	    {
		copy_cpt_to_cpt(&(edge->endp[0]),p1);
		if (prev_edge == NULL)
		{
		    *edgelist = edge->next;
		    free(edge);
		}
		else
		{
		    prev_edge->next = edge->next;
		    free(edge);
		}
		return YES;
	    }

	    prev_edge = edge;
	    edge = edge->next;
	}

	printf("WARNING in find_next_endp(): can't find next endp.\n");
	return NO;
	//clean_up(ERROR);
}

bool same_nor_dir(
	double	*nor1,
	double	*nor2)
{
	int i;
	double n1, n2;
	double tol = 1e-12;

	//n1 = absmax(absmax(nor1[0],nor1[1]),nor1[2]);
	//n2 = absmax(absmax(nor2[0],nor2[1]),nor2[2]);
	if (fabs(nor1[0]) > fabs(nor1[1]) && fabs(nor1[0]) > fabs(nor1[2]))
	    i = 0;
	else if (fabs(nor1[1]) > fabs(nor1[0]) && fabs(nor1[1]) > fabs(nor1[2]))
	    i = 1;
	else
	    i = 2;

	if (nor1[i]*nor2[i] > 0)
	    return YES;
	else
	    return NO;
}

void construct_polygons_on_cfs(CELL *cell)
{
	int i;
	CEDGE *edge;

	for (i = 0; i < 6; i++)
	    set_polyg_edges_on_cf(cell,i);
	set_dir_cf(cell);

	for (i = 0; i < 6; i++)
	    set_polygons_on_cf(cell,i);

	return;
}

void set_polygons_on_cf(
	CELL	*cell,
	int	icf)
{
	CPOINT *p;
	CEDGE *edge;
	CPOLYGON *polyg;
	CFACE *cf;

	cf = &(cell->faces[icf]);

	remove_short_edges_on_cf(cf);

	//debugdan	FIXME
	if (DEBUGDAN)
	{
	    printf("\nDEBUGDAN - set_polygons_on_cf, i = %d.\n", icf);
	    printf("edges_on_ce:\n");
	    edge = cf->edges_on_ce;
	    while (edge)
	    {
		print_edge(edge);
		edge = edge->next;
	    }
	    printf("edges_in_cf:\n");
	    edge = cf->edges_in_cf;
	    while (edge)
	    {
		print_edge(edge);
		edge = edge->next;
	    }
	}
	//debugdan	FIXME
	edge = cf->edges_on_ce;
	if (!edge)
	{
	    printf("ERROR in set_polygons_on_cf(): no directed edges.\n");
	    clean_up(ERROR);
	}
	while (edge)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&(polyg),sizeof(CPOLYGON));
	    init_new_polyg(polyg);
	    add_polyg_vertex(&(edge->endp[1]),polyg);
	    add_polyg_vertex(&(edge->endp[0]),polyg);
	    complete_polyg_on_cf(polyg,cf,cell);
	    polyg->next = cell->cf_polygs;
	    cell->cf_polygs = polyg;
	    edge = cf->edges_on_ce;
	}

	return;
}

void remove_short_edges_on_cf(CFACE *cf)
{
	int i;
	CEDGE *edge, *prev_edge;
	CPOINT midp, p0, p1;

	edge = cf->edges_on_ce;
	prev_edge = NULL;
	while (edge)
	{
	    p0 = edge->endp[0];
	    p1 = edge->endp[1];
	    if (same_cpt(&p0,&p1))
	    {
		//debugdan	FIXME
		/*
		printf("Removed edge (%e, %e, %e) (%e, %e, %e).\n",
			edge->endp[0].crds[0], edge->endp[0].crds[1], edge->endp[0].crds[1],
			edge->endp[1].crds[0], edge->endp[1].crds[1], edge->endp[1].crds[1]);
		*/
		//debugdan	FIXME
		for (i = 0; i < 3; i++)
		    midp.crds[i] = (p0.crds[i]+p1.crds[i])/2;
		if (prev_edge == NULL)
		    cf->edges_on_ce = edge->next;
		else
		    prev_edge->next = edge->next;
		free(edge);

		modify_endps_of_edges(cf,&p0,&midp);
		modify_endps_of_edges(cf,&p1,&midp);
	    }
	    prev_edge = edge;
	    edge = edge->next;
	}

	edge = cf->edges_in_cf;
	prev_edge = NULL;
	while (edge)
	{
	    p0 = edge->endp[0];
	    p1 = edge->endp[1];
	    if (same_cpt(&p0,&p1))
	    {
		//debugdan	FIXME
		/*
		printf("Removed edge (%e, %e, %e) (%e, %e, %e).\n",
			edge->endp[0].crds[0], edge->endp[0].crds[1], edge->endp[0].crds[1],
			edge->endp[1].crds[0], edge->endp[1].crds[1], edge->endp[1].crds[1]);
		*/
		//debugdan	FIXME
		for (i = 0; i < 3; i++)
		    midp.crds[i] = (p0.crds[i]+p1.crds[i])/2;
		if (prev_edge == NULL)
		    cf->edges_in_cf = edge->next;
		else
		    prev_edge->next = edge->next;
		free(edge);

		modify_endps_of_edges(cf,&p0,&midp);
		modify_endps_of_edges(cf,&p1,&midp);
	    }
	    prev_edge = edge;
	    edge = edge->next;
	}

	return;
}

void modify_endps_of_edges(
	CFACE	*cf,
	CPOINT	*p,
	CPOINT	*newp)
{
	CEDGE *edge;

	edge = cf->edges_on_ce;
	while (edge)
	{
	    if (same_cpt(&(edge->endp[0]),p))
		copy_cpt_to_cpt(newp,&(edge->endp[0]));
	    if (same_cpt(&(edge->endp[1]),p))
		copy_cpt_to_cpt(newp,&(edge->endp[1]));
	    edge = edge->next;
	}
	edge = cf->edges_in_cf;
	while (edge)
	{
	    if (same_cpt(&(edge->endp[0]),p))
		copy_cpt_to_cpt(newp,&(edge->endp[0]));
	    if (same_cpt(&(edge->endp[1]),p))
		copy_cpt_to_cpt(newp,&(edge->endp[1]));
	    edge = edge->next;
	}

	return;
}

void add_new_point(
	CPOINT	*p,
	CPOINT	**plist)
{
	CPOINT *newp;

	FT_ScalarMemoryAlloc((POINTER*)&(newp),sizeof(CPOINT));
	copy_cpt_to_cpt(p,newp);
	newp->next = *plist;
	*plist = newp;

	return;
}

void push_back_vert(
	CPOINT		*p,
	CPOLYGON	*polyg)
{
	CPOINT *endv, *newv;

	FT_ScalarMemoryAlloc((POINTER*)&(newv),sizeof(CPOINT));
	copy_cpt_to_cpt(p,newv);
	newv->next = NULL;

	endv = polyg->vertices;
	if (endv == NULL)
	{
	    polyg->vertices = newv;
	    return;
	}

	while(endv->next)
	{
	    endv = endv->next;
	}
	endv->next = newv;

	return;
}

void remove_last_vert(
	CPOLYGON	*polyg)
{
	CPOINT *vert;

	vert = polyg->vertices;
	if (vert == NULL || vert->next == NULL)
	{
	    printf("ERROR in remove_last_vert(): empty polygon.\n");
	    clean_up(ERROR);
	}

	while (vert->next->next)
	{
	    vert = vert->next;
	}
	free(vert->next);
	vert->next = NULL;

	return;
}

bool recursively_find_next_vert(
	CPOINT		*currentp,
	CPOINT		*initp,
	CPOLYGON	*polyg,
	CFACE		*cf)
{
	CPOINT *nextp, *nextplist;
	CEDGE *edge;

	//set nextplist
	nextplist = NULL;
	edge = cf->edges_on_ce;
	while (edge)
	{
	    if (same_cpt(&(edge->endp[0]),currentp))
	    {
		add_new_point(&(edge->endp[1]),&nextplist);
	    }
	    edge = edge->next;
	}
	edge = cf->edges_in_cf;
	while (edge)
	{
	    if (same_cpt(&(edge->endp[0]),currentp))
	    {
		add_new_point(&(edge->endp[1]),&nextplist);
	    }
	    edge = edge->next;
	}

	push_back_vert(currentp,polyg);
	nextp = nextplist;
	while (nextp)
	{
	    if (same_cpt(nextp,initp))
	    //if (same_cpt_double_tol(nextp,initp))
	    {
		return YES;
	    }
	    else if (vertex_exists(nextp,polyg))
	    {
		nextp = nextp->next;
		continue;
	    }
	    else if (recursively_find_next_vert(nextp,initp,polyg,cf))
	    {
		return YES;
	    }
	    nextp = nextp->next;
	}

	remove_last_vert(polyg);
	return NO;
}

void add_p_as_nth_vert(
	CPOINT		*p,
	CPOLYGON	*polyg,
	int		n)
{
	int count;
	CPOINT *vert, *newvert;

	if (n < 3)
	{
	    printf("ERROR in add_p_as_nth_vert(): n should not be smaller than 3.\n");
	    clean_up(ERROR);
	}

	count = 1;
	vert = polyg->vertices;
	while (vert)
	{
	    vert = vert->next;
	    count++;
	    if (count == n-1)
		break;
	}
	if (count != n-1)
	{
	    printf("ERROR in add_p_as_nth_vert().\n");
	    clean_up(ERROR);
	}

	FT_ScalarMemoryAlloc((POINTER*)&(newvert),sizeof(CPOINT));
	copy_cpt_to_cpt(p,newvert);
	newvert->next = vert->next;
	vert->next = newvert;

	return;
}

bool vertex_exists(
	CPOINT		*p,
	CPOLYGON	*polyg)
{
	CPOINT *vert;

	vert = polyg->vertices;
	while (vert)
	{
	    if (same_cpt(vert,p))
		return YES;
	    vert = vert->next;
	}

	return NO;
}

void remove_duplicated_edges(
	CPOLYGON	*polyg,
	CEDGE		**edgelist)
{
	CPOINT *p0, *p1;
	CEDGE	*edge, *preve;
	bool found;

	p0 = polyg->vertices;
	p1 = p0->next;

	while (p1)
	{
	    found = NO;
	    preve = NULL;
	    edge = *edgelist;
	    while (edge)
	    {
		if (same_cpt(p0,&(edge->endp[0])) && 
		    same_cpt(p1,&(edge->endp[1])))
		{
		    if (preve == NULL)
		    {
			*edgelist = edge->next;
			free(edge);
		    }
		    else
		    {
			preve->next = edge->next;
			free(edge);
		    }
		    found = YES;
		    break;
		}
		preve = edge;
		edge = edge->next;
	    }
	    /*
	    if (!found)
	    {
		printf("ERROR in remove_duplicated_edges().\n");
		clean_up(ERROR);
	    }
	    */
	    p0 = p1;
	    p1 = p0->next;
	}

	//last one
	p1 = polyg->vertices;
	found = NO;
	preve = NULL;
	edge = *edgelist;
	while (edge)
	{
	    if (same_cpt(p0,&(edge->endp[0])) && 
		same_cpt(p1,&(edge->endp[1])))
	    {
		if (preve == NULL)
		{
		    *edgelist = edge->next;
		    free(edge);
		}
		else
		{
		    preve->next = edge->next;
		    free(edge);
		}
		found = YES;
		break;
	    }
	    preve = edge;
	    edge = edge->next;
	}
	/*
	if (!found)
	{
	    printf("ERROR in remove_duplicated_edges().\n");
	    clean_up(ERROR);
	}
	*/

	return;
}

void complete_polyg_on_cf(
	CPOLYGON	*polyg,
	CFACE		*cf,
	CELL		*c)
{
	CPOINT *initp, *secondp, *nextp;
	CPOINT *nextplist;
	CEDGE *edge;
	bool found;

	initp = polyg->vertices;
	secondp = initp->next;
	nextplist = NULL;

	edge = cf->edges_on_ce;
	while (edge)
	{
	    if (same_cpt(&(edge->endp[0]),secondp))
	    {
		add_new_point(&(edge->endp[1]),&nextplist);
	    }
	    edge = edge->next;
	}
	edge = cf->edges_in_cf;
	while (edge)
	{
	    if (same_cpt(&(edge->endp[0]),secondp))
	    {
		add_new_point(&(edge->endp[1]),&nextplist);
	    }
	    edge = edge->next;
	}

	found = NO;
	nextp = nextplist;
	while (nextp)
	{
	    if (recursively_find_next_vert(nextp,initp,polyg,cf))
	    {
		found = YES;
		break;
	    }
	    nextp = nextp->next;
	}

	if (!found)
	{
	    printf("ERROR in complete_polyg_on_cf(): function failed.\n");
	    //debugdan	FIXME
	    printf("ctri_polygs:\n");
	    CPOLYGON *polyg = c->ctri_polygs;
	    double nor[3];
	    CPOINT *p0, *p1, *p2;
	    while (polyg)
	    {
		print_polygon(polyg);
		p0 = polyg->vertices;
		p1 = p0->next;
		p2 = p1->next;
		set_nor(p0,p1,p2,nor);
		printf("nor: %e %e %e.\n", nor[0], nor[1], nor[2]);
		polyg = polyg->next;
	    }
	    printf("ctri nor:\n");
	    CTRI *ctri = c->ctris;
	    while (ctri)
	    {
		printf("%e %e %e.\n", ctri->nor[0], ctri->nor[1], ctri->nor[2]);
		ctri = ctri->next;
	    }
	    nextp = nextplist;
	    while (nextp)
	    {
		if (recursively_find_next_vert(nextp,initp,polyg,cf))
		{
		    found = YES;
		    break;
		}
		nextp = nextp->next;
	    }
	    //debugdan	FIXME
	    clean_up(ERROR);
	}

	//remove used edges in edges_on_ce
	remove_duplicated_edges(polyg,&(cf->edges_on_ce));

	return;
}
/*
void complete_polyg_on_cf(
	CPOLYGON	*polyg,
	CFACE		*cf)
{
	int niter;
	CPOINT *startv, *endv, *newv;
	CEDGE *edge, *preve;
	bool found;

	startv = polyg->vertices;
	endv = polyg->vertices;
	if (!endv)
	{
	    printf("ERROR in complete_polyg_on_cf(): no vertex in poly.\n");
	    clean_up(ERROR);
	}
	while (endv->next)
	    endv = endv->next;

	niter = 0;
	while (!same_cpt(startv,endv))
	{
	    niter++;

	    //find prev vertex on edges_on_ce
	    edge = cf->edges_on_ce;
	    preve = NULL;
	    found = NO;
	    while (edge)
	    {
		if (same_cpt(&(edge->endp[1]),startv))
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(newv),sizeof(CPOINT));
		    copy_cpt_to_cpt(&(edge->endp[0]),newv);
		    newv->next = startv;
		    polyg->vertices = newv;
		    startv = newv;
		    //remove edge
		    if (preve == NULL)
		    {
			cf->edges_on_ce = edge->next;
			free(edge);
		    }
		    else
		    {
			preve->next = edge->next;
			free(edge);
		    }
		    found = YES;
		    break;
		}
		else
		{
		    preve = edge;
		    edge = edge->next;
		}
	    }
	    if (found)
		continue;

	    //if not found in edges_on_ce, search edges_in_cf
	    edge = cf->edges_in_cf;
	    while (edge)
	    {
		if (same_cpt(&(edge->endp[1]),startv))
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(newv),sizeof(CPOINT));
		    copy_cpt_to_cpt(&(edge->endp[0]),newv);
		    newv->next = startv;
		    polyg->vertices = newv;
		    startv = newv;
		    break;
		}
		else
		    edge = edge->next;
	    }

	    if (niter > 5000)
	    {
		printf("ERROR in complete_polyg_on_cf(): too many iteration steps.\n");
		printf("Polygon vertices already set:\n");
		newv = polyg->vertices;
		while (newv)
		{
		    print_point(newv);
		    printf("\n");
		    newv = newv->next;
		}
		printf("Edges in cell face:\n");
		edge = cf->edges_in_cf;
		while (edge)
		{
		    print_edge(edge);
		    edge = edge->next;
		}
		clean_up(ERROR);
	    }
	}

	return;
}
*/
/*
void complete_polyg_on_cf(
	CPOLYGON	*polyg,
	CFACE		*cf)
{
	int niter;
	CPOINT *startv, *endv, *newv;
	CEDGE *edge;

	startv = polyg->vertices;
	endv = polyg->vertices;
	if (!endv)
	{
	    printf("ERROR in complete_polyg_on_cf(): no vertex in poly.\n");
	    clean_up(ERROR);
	}
	while (endv->next)
	    endv = endv->next;

	niter = 0;
	while (!same_cpt(startv,endv))
	{
	    niter++;
	    edge = cf->edges_in_cf;
	    if (!edge)
	    {
		printf("ERROR in complete_polyg_on_cf(): no directed edges in cell face.\n");
		clean_up(ERROR);
	    }
	    while (edge)
	    {
		if (same_cpt(&(edge->endp[1]),startv))
		{
		    newv = (CPOINT*)malloc(sizeof(CPOINT));
		    copy_cpt_to_cpt(&(edge->endp[0]),newv);
		    newv->next = startv;
		    polyg->vertices = newv;
		    startv = newv;
		    break;
		}
		else
		    edge = edge->next;
	    }
	    if (niter > 5000)
	    {
		printf("ERROR in complete_polyg_on_cf(): too many iteration steps.\n");
		printf("Polygon vertices already set:\n");
		newv = polyg->vertices;
		while (newv)
		{
		    print_point(newv);
		    printf("\n");
		    newv = newv->next;
		}
		printf("Edges in cell face:\n");
		edge = cf->edges_in_cf;
		while (edge)
		{
		    print_edge(edge);
		    edge = edge->next;
		}
		clean_up(ERROR);
	    }
	}

	return;
}
*/
void add_next_vertices_on_ce(
	CPOLYGON	*polyg,
	CFACE		*cf)
{
	CPOINT *startv, *endv, *newv;
	CEDGE *edge, *preve;
	bool found;

	startv = polyg->vertices;
	endv = polyg->vertices;
	if (!endv)
	{
	    printf("ERROR in add_next_vertices_on_ce(): no vertices found in polyg.\n");
	    clean_up(ERROR);
	}
	while (endv->next)
	    endv = endv->next;

	while (cf->edges_on_ce)
	{
	    edge = cf->edges_on_ce;
	    preve = NULL;
	    if (!edge)
	    {
		printf("ERROR in add_next_vertices_on_ce(): no directed edges.\n");
		clean_up(ERROR);
	    }
	    found = NO;
	    while (edge)
	    {
		if (same_cpt(&(edge->endp[0]),endv))
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(newv),sizeof(CPOINT));
		    newv->next = NULL;
		    copy_cpt_to_cpt(&(edge->endp[1]),newv);
		    endv->next = newv;
		    endv = newv;
		    //delete edge
		    if (!preve)
		    {
			cf->edges_on_ce = edge->next;
			free(edge);
			found = YES;
			break;
		    }
		    else
		    {
			preve->next = edge->next;
			free(edge);
			found = YES;
			break;
		    }
		}
		else
		{
		    preve = edge;
		    edge = edge->next;
		}
	    }
	    if (!found || same_cpt(startv,endv))
		break;
	}
}

void add_prev_vertices_on_ce(
	CPOLYGON	*polyg,
	CFACE		*cf)
{
	CPOINT *startv, *endv, *newv;
	CEDGE *edge, *preve;
	bool found;

	startv = polyg->vertices;
	endv = polyg->vertices;
	if (!endv)
	{
	    printf("ERROR in add_prev_vertices_on_ce(): no vertices found in polyg.\n");
	    clean_up(ERROR);
	}
	while (endv->next)
	    endv = endv->next;

	while (cf->edges_on_ce)
	{
	    edge = cf->edges_on_ce;
	    preve = NULL;
	    if (!edge)
	    {
		printf("ERROR in add_prev_vertices_on_ce(): no directed edges.\n");
		clean_up(ERROR);
	    }
	    found = NO;
	    while (edge)
	    {
		if (same_cpt(&(edge->endp[1]),startv))
		{
		    FT_ScalarMemoryAlloc((POINTER*)&(newv),sizeof(CPOINT));
		    newv->next = NULL;
		    copy_cpt_to_cpt(&(edge->endp[0]),newv);
		    newv->next = startv;
		    polyg->vertices = newv;
		    startv = newv;
		    //delete edge
		    if (!preve)
		    {
			cf->edges_on_ce = edge->next;
			free(edge);
			found = YES;
			break;
		    }
		    else
		    {
			preve->next = edge->next;
			free(edge);
			found = YES;
			break;
		    }
		}
		else
		{
		    preve = edge;
		    edge = edge->next;
		}
	    }
	    if (!found || same_cpt(startv,endv))
		break;
	}
}

void set_polyg_edges_on_cf(
	CELL	*cell,
	int	index)
{
	int i, j, edir;
	CPOINT *p, *startp, *prevp;
	CPOINT midp;
	CEDGE *edge;
	CFACE *cf = &(cell->faces[index]);
	bool eicf;	//edge in cell face

	//debugdan	FIXME
	if (DEBUGDAN)
	{
	    printf("\nDEBUGDAN - set_polyg_edges_on_cf, i = %d.\n", index);
	    printf("undirected_edges:\n");
	}
	//debugdan	FIXME

	//set directed edges from undirected edges
	edge = cf->undirected_edges;
	while (edge)
	{
	    //debugdan	FIXME
	    if (DEBUGDAN)
		print_edge(edge);
	    //debugdan	FIXME
	    edir = edge_dir(edge,cell);
	    for (i = 0; i < 3; i++)
	    {
		midp.crds[i] = (edge->endp[0].crds[i] + edge->endp[1].crds[i])/2.0;
	    }
	    if (point_in_cell_face(&(edge->endp[0]),cf) ||
		point_in_cell_face(&(edge->endp[1]),cf) ||
		point_in_cell_face(&midp,cf))
	    {
		if (edir == -1)
		    add_new_edge(&(edge->endp[0]),&(edge->endp[1]),&(cf->edges_in_cf));
		else
		    add_new_edge(&(edge->endp[1]),&(edge->endp[0]),&(cf->edges_in_cf));
	    }
	    else
	    {
		if (edir == -1)
		    add_new_edge(&(edge->endp[0]),&(edge->endp[1]),&(cf->edges_on_ce));
		else
		    add_new_edge(&(edge->endp[1]),&(edge->endp[0]),&(cf->edges_on_ce));
	    }
	    edge = edge->next;
	}

	//complete directed edges on cf edges

	//set crxp list on cf edges
	//edge = cf->undirected_edges;
	edge = cf->edges_in_cf;
	while (edge)
	{
	    p = &(edge->endp[0]);
	    if (point_on_cf_edge(p,cf,&i))
		add_crxp_on_cf_edge_sorted(p,cf,i,0);

	    p = &(edge->endp[1]);
	    if (point_on_cf_edge(p,cf,&i))
		add_crxp_on_cf_edge_sorted(p,cf,i,1);

	    edge = edge->next;
	}
	edge = cf->edges_on_ce;
	while (edge)
	{
	    p = &(edge->endp[0]);
	    if (point_on_cf_edge(p,cf,&i))
		add_crxp_on_cf_edge_sorted(p,cf,i,0);

	    p = &(edge->endp[1]);
	    if (point_on_cf_edge(p,cf,&i))
		add_crxp_on_cf_edge_sorted(p,cf,i,1);

	    edge = edge->next;
	}

	//set directed edges on edges of cell face from a starting point
	//set directed edges on the first edge which has at least one crxp on it
	for (i = 0; i < 4; i++)
	{
	    if (cf->edges[i].crxps)
	    {
		startp = cf->edges[i].crxps;
		break;
	    }
	}
	if (i == 4)
	{
	    return;
	}
	if (!same_cpt(&(cf->pts[i]),startp))
	{
	    add_new_dir_edge(startp,&(cf->pts[i]),&(cf->edges_on_ce));
	    cf->pts[i].flag = startp->flag;
	}
	prevp = startp;
	p = startp->next;
	while (p)
	{
	    add_new_dir_edge(prevp,p,&(cf->edges_on_ce));
	    p->flag = (prevp->flag+1)%2;
	    prevp = p;
	    p = p->next;
	}
	add_new_dir_edge(prevp,&(cf->pts[(i+1)%4]),&(cf->edges_on_ce));
	cf->pts[(i+1)%4].flag = prevp->flag;

	//set directed edges for rest edges
	for (j = i+1; j < i+4; j++)
	{
	    prevp = &(cf->pts[j%4]);
	    p = cf->edges[j%4].crxps;
	    if (!p)
	    {
		add_new_dir_edge(prevp,&(cf->pts[(j+1)%4]),&(cf->edges_on_ce));
		cf->pts[(j+1)%4].flag = prevp->flag;
		continue;
	    }
	    if (same_cpt(prevp,p))
	    {
		prevp = p;
		p = p->next;
	    }
	    while (p)
	    {
		add_new_dir_edge(prevp,p,&(cf->edges_on_ce));
		p->flag = (prevp->flag+1)%2;
		prevp = p;
		p = p->next;
	    }
	    add_new_dir_edge(prevp,&(cf->pts[(j+1)%4]),&(cf->edges_on_ce));
	    cf->pts[(j+1)%4].flag = prevp->flag;
	}

	//if there're crxing proints in cf but no crxing points on cf edges
	//TODO???

	return;
}

void set_dir_cf(CELL *cell)
{
	int starti, i, j;
	CEDGE *edge;
	CFACE *cf;

	for (starti = 0; starti < 6; starti++)
	{
	    cf = &(cell->faces[starti]);
	    if (cf->edges_in_cf || cf->edges_on_ce)
		break;
	}
	//start from the neighbor of a cell face which has directed edges
	for (i = starti+3; i < starti+9; i++)
	{
	    j = i%6;
	    cf = &(cell->faces[j]);
	    if (cf->edges_in_cf == NULL && cf->edges_on_ce == NULL)
		set_dir_edges_by_nb_cf(cf,cell,j);
	}
}

void set_dir_edges_by_nb_cf(
	CFACE	*cf,
	CELL	*cell,
	int	icf)
{
	int i, j, dir;
	CPOINT *p0, *p1;
	CEDGE *edge;
	CFACE *ncf;
	bool found = NO;

	for (i = 0; i < 6; i++)
	{
	    if (i == icf)
		continue;
	    ncf = &(cell->faces[i]);
	    if (is_nb_cf(icf,i) && ncf->edges_on_ce)
	    {
		for (j = 0; j < 4; j++)
		{
		    p0 = &(cf->pts[j]);
		    p1 = &(cf->pts[(j+1)%4]);
		    if (point_on_cell_face(p0,ncf) &&
			point_on_cell_face(p1,ncf))
		    {
			edge = ncf->edges_on_ce;
			while (edge)
			{
			    if (same_cpt(p0,&(edge->endp[0])) &&
				same_cpt(p1,&(edge->endp[1])))
			    {
				add_new_edge(p1,p0,&(cf->edges_on_ce));
				dir = -1;
				found = YES;
			    }
			    else if (same_cpt(p0,&(edge->endp[1])) &&
				     same_cpt(p1,&(edge->endp[0])))
			    {
				add_new_edge(p0,p1,&(cf->edges_on_ce));
				dir = 1;
				found = YES;
			    }
			    if (found)
				break;
			    edge = edge->next;
			}
			if (found)
			    break;
		    }
		}
		if (found)
		    break;
	    }
	    if (found)
		break;
	}

	if (!found)
	{
	    printf("ERROR in set_dir_edges_by_nb_cf(): can't find directed neighboring cf.\n");
	    clean_up(ERROR);
	}

	for (i = j+1; i < j+4; i++)
	{
	    if (dir == -1)
		add_new_edge(&(cf->pts[(i+1)%4]),&(cf->pts[i%4]),&(cf->edges_on_ce));
	    else
		add_new_edge(&(cf->pts[i%4]),&(cf->pts[(i+1)%4]),&(cf->edges_on_ce));
	}

	return;
}

bool is_nb_cf(
	int	i,
	int	j)
{
	if ((i == 0 && j == 5) ||
	    (i == 1 && j == 3) ||
	    (i == 2 && j == 4) ||
	    (i == 3 && j == 1) ||
	    (i == 4 && j == 2) ||
	    (i == 5 && j == 0))
	    return NO;
	return YES;
}

int edge_dir(
	CEDGE	*edge,
	CELL	*cell)
{
	CPOLYGON *polyg;
	CPOINT *p, *startp;

	polyg = cell->ctri_polygs;
	while (polyg)
	{
	    startp = polyg->vertices;
	    p = polyg->vertices;
	    if (!p)
	    {
		printf("ERROR in edge_dir(): empty polyg.\n");
		clean_up(ERROR);
	    }
	    while (p)
	    {
		if (same_cpt(&(edge->endp[0]),p))
		{
		    if (p->next)
		    {
			if (same_cpt(&(edge->endp[1]),p->next))
			    return 1;
		    }
		    else
		    {
			if (same_cpt(&(edge->endp[1]),startp))
			    return 1;
		    }
		}
		else if (same_cpt(&(edge->endp[1]),p))
		{
		    if (p->next)
		    {
			if (same_cpt(&(edge->endp[0]),p->next))
			    return -1;
		    }
		    else
		    {
			if (same_cpt(&(edge->endp[0]),startp))
			    return -1;
		    }
		}
		p = p->next;
	    }
	    polyg = polyg->next;
	}
}

bool point_on_cf_edge(
	CPOINT	*p,
	CFACE	*cf,
	int	*index)
{
	int i, j, k, dir;
	double tol = 1e-12;
	bool found = NO;

	dir = cf->dir;

	if (fabs(p->crds[dir] - cf->l[dir]) > tol)
	{
	    printf("ERROR in point_on_cf_edge(): point is not on cell face.\n");
	    clean_up(ERROR);
	}

	if (point_in_cell_face(p,cf))
	    return NO;

	for (i = 0; i < 4; i++)
	{
	    if (pt_between_two_pts(p,&(cf->pts[i]),&(cf->pts[(i+1)%4])))
	    {
		*index = i;
		return YES;
	    }
	}

	return NO;
}

bool pt_between_two_pts(
	CPOINT	*p0,
	CPOINT	*p1,
	CPOINT	*p2)
{
	int i;
	double v1[3], v2[3];
	double tol = 1e-12;

	for (i = 0; i < 3; i++)
	{
	    v1[i] = p1->crds[i] - p0->crds[i];
	    v2[i] = p2->crds[i] - p0->crds[i];
	}

	/*p0 can be the same as p1, but can't be the same as p2.*/
	if (fabs(v2[0]) < tol && fabs(v2[1]) < tol && fabs(v2[2]) < tol)
	    return NO;

	if (fabs(v1[1]*v2[2]-v1[2]*v2[1]) < tol &&
	    fabs(v1[2]*v2[0]-v1[0]*v2[2]) < tol &&
	    fabs(v1[0]*v2[1]-v1[1]*v2[0]) < tol)
	{
	    return YES;
	}

	return NO;
}

void add_crxp_on_cf_edge_sorted(
	CPOINT	*p,
	CFACE	*cf,
	int	index,
	int	flag)
{
	int dir;
	double f1, f2;
	CPOINT *crxp, *prevp, *newp;
	CEDGE *edge;
	double tol = 1e-12;

	FT_ScalarMemoryAlloc((POINTER*)&(newp),sizeof(CPOINT));
	newp->next = NULL;
	copy_cpt_to_cpt(p,newp);
	newp->flag = flag;

	edge = &(cf->edges[index]);
	crxp = edge->crxps;

	if (edge->crxps == NULL)
	{
	    newp->next = edge->crxps;
	    edge->crxps = newp;
	    return;
	}

	for (dir = 0; dir < 3; dir++)
	{
	    if (fabs(newp->crds[dir] - cf->pts[index].crds[dir]) > tol)
		break;
	}

	/*p is the same as cf->pts[index]*/
	if (dir == 3)
	{
	    newp->next = edge->crxps;
	    edge->crxps = newp;
	    return;
	}

	f1 = newp->crds[dir] - cf->pts[index].crds[dir];
	prevp = NULL;
	while (crxp)
	{
	    f2 = crxp->crds[dir] - cf->pts[index].crds[dir];
	    if (fabs(f1) < fabs(f2))
	    {
		newp->next = crxp;
		if (!prevp)
		    edge->crxps = newp;
		else
		    prevp->next = newp;
		return;
	    }
	    prevp = crxp;
	    crxp = crxp->next;
	}
	newp->next = crxp;
	prevp->next = newp;

	return;
}

bool the_ctri(CTRI *ctri)
{
	CPOINT p0, p1, p2;

	p0.crds[0] = -0.0251126955278;
	p0.crds[1] = -0.0251126954227;
	p0.crds[2] = 2.19961122047;

	p1.crds[0] = 0.0251126950244;
	p1.crds[1] = -0.0251126955093;
	p1.crds[2] = 2.19961122039;
/*
	p2.crds[0] = 0.525;
	p2.crds[1] = 0.375;
	p2.crds[2] = 1.83052048782;

	if (debug_same_cpt(&p0,&(ctri->pts[0])) &&
	    debug_same_cpt(&p1,&(ctri->pts[1])) &&
	    debug_same_cpt(&p2,&(ctri->pts[2])))
*/
	if (debug_same_cpt(&p0,&(ctri->pts[0])) &&
	    debug_same_cpt(&p1,&(ctri->pts[1])))
	    return YES;

	return NO;
}

bool the_cpt(CPOINT *p)
{
	CPOINT thep;

	thep.crds[0] = 0.9;
	thep.crds[1] = 0.026051612160407842;
	thep.crds[2] = 2.2;

	if (debug_same_cpt(&thep,p))
	    return YES;
	return NO;
}

bool debug_same_cpt(
	CPOINT	*p1,
	CPOINT	*p2)
{
	int i;
	double tol = 1e-8;

	for (i = 0; i < 3; i++)
	    if (fabs(p1->crds[i] - p2->crds[i]) > tol)
		return NO;

	return YES;
}

void G_CARTESIAN::cft_construct_cell_polyhedrons()
{
	int i, j, k, index;
	CELL *c;
	CPOLYHEDRON *polyh;

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);

	    if (c->ctri_polygs == NULL)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&polyh,sizeof(CPOLYHEDRON));
		init_polyh(polyh,c);
		cell_to_polyh(c,polyh);
		c->polyhs = polyh;
		continue;
	    }
	    /*
	    if (c->ctri_polygs == NULL)
		continue;
	    */

	    //set sets of connected polygons
	    set_scps(c);

	    //set polyhedrons based on scps
	    set_polyhedrons_in_cell(c);
	}

	return;
}

void set_polyhedrons_in_cell(CELL *c)
{
	CPOLYHEDRON *polyh;
	SETOFCPOLYGS *scp;
	CBOUNDARY *bdry;

	if (c->scpocs == NULL)
	    return;

	//init_scps_bdries_marks(c);
	init_scps_marks(c);

	scp = c->scpocs;
	while (scp)
	{
	    if (scp->inpolyh)
	    {
		scp = scp->next;
		continue;
	    }

	    init_scpics_marks(c);

	    FT_ScalarMemoryAlloc((POINTER*)&polyh,sizeof(CPOLYHEDRON));
	    init_polyh(polyh,c);
	    polyh->scp_dir = scp->dir;
	    add_scp_to_polyh(scp,polyh);

	    if (polyh->boundaries)
	    {
		//find scp with inverse polyh->boundaries (should be unique) and add scp
		add_scp_with_ibdries(polyh);
	    }

	    polyh->next = c->polyhs;
	    c->polyhs = polyh;

	    scp = scp->next;
	}

	return;
}

void init_scpics_marks(CELL *c)
{
	SETOFCPOLYGS *scp;

	scp = c->scpics;
	while (scp)
	{
	    scp->inpolyh = NO;
	    scp = scp->next;
	}

	return;
}

void init_scps_marks(CELL *c)
{
	SETOFCPOLYGS *scp;

	scp = c->scpocs;
	while (scp)
	{
	    scp->oncf = YES;
	    scp = scp->next;
	}

	scp = c->scpics;
	while (scp)
	{
	    scp->oncf = NO;
	    scp = scp->next;
	}

	return;
}

void init_scps_bdries_marks(CELL *c)
{
	SETOFCPOLYGS *scp;
	CBOUNDARY *bdry;

	scp = c->scpocs;
	while (scp)
	{
	    bdry = scp->boundaries;
	    while (bdry)
	    {
		bdry->mark = 0;
		bdry = bdry->next;
	    }
	    scp = scp->next;
	}

	scp = c->scpics;
	while (scp)
	{
	    bdry = scp->boundaries;
	    while (bdry)
	    {
		bdry->mark = 0;
		bdry = bdry->next;
	    }
	    scp = scp->next;
	}

	return;
}

void add_scp_with_ibdries(CPOLYHEDRON *polyh)
{
	CBOUNDARY *bdry, *ibdry;
	SETOFCPOLYGS *scp;

	bdry = polyh->boundaries;
	while (bdry)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&ibdry,sizeof(CBOUNDARY));
	    inverse_bdry(bdry,ibdry);
	    if (find_scp_with_bdry(&scp,ibdry,polyh))
	    {
		//bdry->mark++;
		add_scp_to_polyh(scp,polyh);
	    }

	    polyh->boundaries = bdry->next;
	    free(bdry);
	    bdry = polyh->boundaries;
	}

	return;
}

bool find_scp_with_bdry(
	SETOFCPOLYGS	**scp,
	CBOUNDARY	*bdry,
	CPOLYHEDRON	*polyh)
{
	CEDGE *edge;
	SETOFCPOLYGS *scp0;

	edge = bdry->edges;

	scp0 = polyh->cell->scpocs;
	while (scp0)
	{
	    if (scp0->inpolyh || (scp0->dir != polyh->scp_dir))
	    {
		scp0 = scp0->next;
		continue;
	    }

	    if (find_bdry_with_edge(scp0->boundaries,edge))
	    {
		*scp = scp0;
		return YES;
	    }

	    scp0 = scp0->next;
	}

	scp0 = polyh->cell->scpics;
	while (scp0)
	{
	    if (scp0->inpolyh)
	    {
		scp0 = scp0->next;
		continue;
	    }

	    if (find_bdry_with_edge(scp0->boundaries,edge))
	    {
		*scp = scp0;
		return YES;
	    }

	    scp0 = scp0->next;
	}

	return NO;
}

bool find_bdry_with_edge(
	CBOUNDARY	*bdries,
	CEDGE		*edge)
{
	CEDGE *edge0;
	CBOUNDARY *bdry;

	bdry = bdries;
	while (bdry)
	{
	    edge0 = bdry->edges;
	    while (edge0)
	    {
		if (same_dir_edge(edge0,edge))
		{
		    bdry->mark++;
		    return YES;
		}
		edge0 = edge0->next;
	    }
	    bdry = bdry->next;
	}

	return NO;
}

bool same_dir_edge(
	CEDGE	*edge1,
	CEDGE	*edge2)
{
	if ((same_cpt(&(edge1->endp[0]),&(edge2->endp[0])) &&
	    same_cpt(&(edge1->endp[1]),&(edge2->endp[1]))))
	    return YES;

	return NO;
}

void inverse_bdry(
	CBOUNDARY	*bdry,
	CBOUNDARY	*ibdry)
{
	CEDGE *edge;

	edge = bdry->edges;
	while (edge)
	{
	    add_new_edge(&(edge->endp[1]),&(edge->endp[0]),&(ibdry->edges));
	    edge = edge->next;
	}

	//ibdry->mark = bdry->mark;

	return;
}

void add_scp_to_polyh(
	SETOFCPOLYGS	*scp,
	CPOLYHEDRON	*polyh)
{
	//add scp->polygs to polyh
	add_scp_polygs_to_polyh(scp,polyh);
	//add scp->boundaries to polyh
	add_scp_bdries_to_polyh(scp,polyh);
	scp->inpolyh = YES;

	return;
}

void add_scp_polygs_to_polyh(
	SETOFCPOLYGS	*scp,
	CPOLYHEDRON	*polyh)
{
	CPOLYGON *polyg;

	polyg = scp->polygs;
	while (polyg)
	{
	    add_new_polyg(polyg,&(polyh->faces));
	    polyg = polyg->next;
	}

	return;
}

void add_scp_bdries_to_polyh(
	SETOFCPOLYGS	*scp,
	CPOLYHEDRON	*polyh)
{
	CBOUNDARY *bdry;

	bdry = scp->boundaries;
	while (bdry)
	{
	    if (scp->oncf && (bdry->mark > 0))
	    {
		bdry = bdry->next;
	    }

	    add_new_bdry(bdry,&(polyh->boundaries));
	    bdry = bdry->next;
	}

	return;
}

void add_new_bdry(
	CBOUNDARY	*bdry,
	CBOUNDARY	**bdrylist)
{
	CEDGE *edge;
	CBOUNDARY *newbdry;

	FT_ScalarMemoryAlloc((POINTER*)&newbdry,sizeof(CBOUNDARY));
	newbdry->next = NULL;

	edge = bdry->edges;
	while (edge)
	{
	    add_edge(edge,&(newbdry->edges));
	    edge = edge->next;
	}

	newbdry->mark = bdry->mark;

	newbdry->next = *bdrylist;
	*bdrylist = newbdry;

	return;
}

void init_polyh(
	CPOLYHEDRON	*polyh,
	CELL		*c)
{
	polyh->faces = NULL;
	polyh->boundaries = NULL;
	polyh->cell = c;
	polyh->sorted_nbcs = NULL;
	polyh->pam = NULL;
	polyh->next = NULL;
	polyh->iscell = false;
	polyh->merged = false;

	return;
}

void set_scps(CELL *c)
{
	CEDGE *edge;
	CPOLYGON *polyg;
	SETOFCPOLYGS *scp;

	init_polygs_marks(c);
	init_polygs_edges(c);

	//sets of connected polygons on cell faces
	polyg = c->cf_polygs;
	while (polyg)
	{
	    if (polyg->inscp)
	    {
		polyg = polyg->next;
		continue;
	    }

	    FT_ScalarMemoryAlloc((POINTER*)&scp,sizeof(SETOFCPOLYGS));
	    init_scp(scp);
	    insert_scp(scp,&(c->scpocs));

	    add_polyg_to_scp(polyg,scp);
	    polyg->inscp = YES;

	    BFS_add_nb_polygs(polyg,scp,c->cf_polygs);

	    reset_scp_boundaries(scp);

	    //set direction for scp on cell face
	    set_scp_dir(scp,c);

	    polyg = polyg->next;
	}

	//sets of connected polygons in cell
	polyg = c->ctri_polygs;
	while (polyg)
	{
	    if (polyg->inscp)
	    {
		polyg = polyg->next;
		continue;
	    }

	    FT_ScalarMemoryAlloc((POINTER*)&scp,sizeof(SETOFCPOLYGS));
	    init_scp(scp);
	    insert_scp(scp,&(c->scpics));

	    add_polyg_to_scp(polyg,scp);
	    polyg->inscp = YES;

	    BFS_add_nb_polygs(polyg,scp,c->ctri_polygs);

	    reset_scp_boundaries(scp);

	    polyg = polyg->next;
	}

	return;
}

void set_scp_dir(
	SETOFCPOLYGS	*scp,
	CELL		*c)
{
	int i;
	double v1[3], v2[3], nor[3];
	CPOINT *p0, *p1, *p2;
	double tol = 1e-12;

	p0 = scp->polygs->vertices;
	p1 = p0->next;
	p2 = p1->next;

	while (true)
	{
	    if (same_cpt(p0,p1) || same_cpt(p1,p2))
	    {
		printf("ERROR in set_scp_dir(): duplicated points in a polygon.\n");
		clean_up(ERROR);
	    }

	    if (!set_nor(p0,p1,p2,nor))
	    {
		p0 = p1;
		p1 = p2;
		p2 = p2->next;
		if (p2 == NULL)
		{
		    printf("ERROR in set_scp_dir().\n");
		    clean_up(ERROR);
		}
		continue;
	    }
	    break;
	}

	for (i = 0; i < 3; i++)
	{
	    if (fabs(p0->crds[i] - p1->crds[i]) < tol &&
		fabs(p1->crds[i] - p2->crds[i]) < tol)
		break;
	}

	if (fabs(p0->crds[i] - c->celll[i]) < tol)
	{
	    if (nor[i] > 0)
		scp->dir = INW;
	    else
		scp->dir = OUTW;
	}
	else if (fabs(p0->crds[i] - c->cellu[i]) < tol)
	{
	    if (nor[i] > 0)
		scp->dir = OUTW;
	    else
		scp->dir = INW;
	}
	else
	{
	    printf("ERROR in set_scp_dir(): polygon is not on cell face.\n");
	    clean_up(ERROR);
	}

	return;
}

void reset_scp_boundaries(SETOFCPOLYGS *scp)
{
	CBOUNDARY *bdry, *newbdry;
	int count = 0;

	bdry = scp->boundaries;
	scp->boundaries = NULL;

	if (bdry->edges == NULL)
	{
	    printf("ERROR in reset_scp_boundaries(): empty boundary.\n");
	    clean_up(ERROR);
	}
	while (bdry->edges)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&newbdry,sizeof(CBOUNDARY));
	    complete_bdry_with_init_edge(newbdry,&(bdry->edges));

	    newbdry->next = scp->boundaries;
	    scp->boundaries = newbdry;

	    count++;
	    if (count > 1000)
	    {
		printf("ERROR in reset_scp_boundaries().\n");
		clean_up(ERROR);
	    }
	}

	free(bdry);
	return;
}

void complete_bdry_with_init_edge(
	CBOUNDARY	*bdry,
	CEDGE		**init_edge)
{
	CPOINT p0, p1;
	CEDGE *edge, *preve;
	int count = 0;

	edge = *init_edge;
	copy_cpt_to_cpt(&(edge->endp[0]),&p0);
	copy_cpt_to_cpt(&(edge->endp[1]),&p1);
	add_edge(edge,&(bdry->edges));
	//remove first edge in bdry->edges
	*init_edge = edge->next;
	free(edge);

	while (true)
	{
	    edge = *init_edge;
	    preve = NULL;
	    while (edge)
	    {
		if (same_cpt(&(edge->endp[1]),&p0))
		{
		    add_edge(edge,&(bdry->edges));
		    copy_cpt_to_cpt(&(edge->endp[0]),&p0);
		    if (preve == NULL)
		    {
			*init_edge = edge->next;
		    }
		    else
		    {
			preve->next = edge->next;
		    }
		    break;
		}
		else
		{
		    preve = edge;
		    edge = edge->next;
		}
	    }

	    if (same_cpt(&(edge->endp[0]),&p1))
	    {
		free(edge);
		break;
	    }
	    else
	    {
		free(edge);
	    }

	    count++;
	    if (count > 10000)
	    {
		printf("ERROR in complete_bdry_with_init_edge().\n");
		clean_up(ERROR);
	    }
	}

	return;
}

void insert_scp(
	SETOFCPOLYGS	*scp,
	SETOFCPOLYGS	**scplist)
{
	scp->next = *scplist;
	*scplist = scp;

	return;
}

void init_scp(SETOFCPOLYGS *scp)
{
	scp->polygs = NULL;
	scp->boundaries = NULL;
	scp->dir = UNSET;
	scp->inpolyh = NO;

	return;
}

void BFS_add_nb_polygs(
	CPOLYGON	*polyg,
	SETOFCPOLYGS	*scp,
	CPOLYGON	*polyglist)
{
	CEDGE *edge;
	CEDGE ie;
	CPOLYGON *nb_polyg;

	edge = polyg->edges;
	while (edge)
	{
	    ie.endp[0] = edge->endp[1];
	    ie.endp[1] = edge->endp[0];
	    if (!find_polyg_with_edge(&nb_polyg,&ie,polyglist))
	    {
		//edge is boundary
		add_scp_boundary(edge,scp);
		edge = edge->next;
		continue;
	    }
	    else
	    {
		if (!nb_polyg->inscp)
		{
		    add_polyg_to_scp(nb_polyg,scp);
		    nb_polyg->inscp = YES;
		    BFS_add_nb_polygs(nb_polyg,scp,polyglist);
		}
		edge = edge->next;
	    }
	}

	return;
}

void add_scp_boundary(
	CEDGE		*edge,
	SETOFCPOLYGS	*scp)
{
	if (scp->boundaries == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&(scp->boundaries),sizeof(CBOUNDARY));
	    scp->boundaries->edges = NULL;
	    scp->boundaries->next = NULL;
	}
	add_edge(edge,&(scp->boundaries->edges));

	return;
}

bool find_polyg_with_edge(
	CPOLYGON	**polyg,
	CEDGE		*edge,
	CPOLYGON	*polyglist)
{
	CEDGE *cpedge;
	CPOLYGON *cpolyg;

	cpolyg = polyglist;
	while (cpolyg)
	{
	    cpedge = cpolyg->edges;
	    if (cpedge == NULL)
	    {
		printf("ERROR in find_polyg_with_edge(): empty polygon.\n");
		clean_up(ERROR);
	    }
	    while (cpedge)
	    {
		if (same_cpt(&(cpedge->endp[0]),&(edge->endp[0])) &&
		    same_cpt(&(cpedge->endp[1]),&(edge->endp[1])))
		{
		    //copy_polyg_to_polyg(cpolyg,polyg);
		    *polyg = cpolyg;
		    return YES;
		}
		cpedge = cpedge->next;
	    }

	    cpolyg = cpolyg->next;
	}

	*polyg = NULL;
	return NO;
}

void mark_polyg_edges(CPOLYGON *polyg)
{
	CEDGE *edge;

	edge = polyg->edges;
	while (edge)
	{
	    edge->mark++;
	    edge = edge->next;
	}

	return;
}

void add_polyg_to_scp(
	CPOLYGON	*polyg,
	SETOFCPOLYGS	*scp)
{
	add_new_polygon(polyg,&(scp->polygs));

	return;
}

void init_polygs_marks(CELL *c)
{
	CPOLYGON *polyg;

	polyg = c->ctri_polygs;
	while (polyg)
	{
	    polyg->mark = 0;
	    polyg->inscp = FALSE;
	    polyg->oncf = FALSE;
	    polyg = polyg->next;
	}

	polyg = c->cf_polygs;
	while (polyg)
	{
	    polyg->mark = 0;
	    polyg->inscp = FALSE;
	    polyg->oncf = TRUE;
	    polyg = polyg->next;
	}

	return;
}

void init_polygs_edges(CELL *c)
{
	CPOLYGON *polyg;
	CEDGE *edge;
	CPOINT *p0, *p1;

	polyg = c->ctri_polygs;
	while (polyg)
	{
	    p0 = polyg->vertices;
	    p1 = p0->next;
	    while (p1)
	    {
		add_new_edge(p0,p1,&(polyg->edges));
		p0 = p1;
		p1 = p0->next;
	    }
	    p1 = polyg->vertices;
	    add_new_edge(p0,p1,&(polyg->edges));
	    polyg = polyg->next;
	}

	polyg = c->cf_polygs;
	while (polyg)
	{
	    p0 = polyg->vertices;
	    p1 = p0->next;
	    while (p1)
	    {
		add_new_edge(p0,p1,&(polyg->edges));
		p0 = p1;
		p1 = p0->next;
	    }
	    p1 = polyg->vertices;
	    add_new_edge(p0,p1,&(polyg->edges));
	    polyg = polyg->next;
	}

	return;
}

void add_new_polygon(
	CPOLYGON	*polyg,
	CPOLYGON	**polyglist)
{
	CPOLYGON *newpolyg;

	FT_ScalarMemoryAlloc((POINTER*)&(newpolyg),sizeof(CPOLYGON));
	init_new_polyg(newpolyg);
	copy_polyg_to_polyg(polyg,newpolyg);

	newpolyg->next = *polyglist;
	*polyglist = newpolyg;

	return;
}

void copy_polyg_to_polyg(
	CPOLYGON	*polyg,
	CPOLYGON	*newpolyg)
{
	CPOINT *p, *newp, *prevp;
	CEDGE *edge;

	p = polyg->vertices;
	prevp = NULL;
	while (p)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&newp,sizeof(CPOINT));
	    newp->next = NULL;
	    copy_cpt_to_cpt(p,newp);
	    if (prevp == NULL)
	    {
		newpolyg->vertices = newp;
	    }
	    else
	    {
		prevp->next = newp;
	    }
	    prevp = newp;
	    p = p->next;
	}

	edge = polyg->edges;
	while (edge)
	{
	    add_edge(edge,&(newpolyg->edges));
	    edge = edge->next;
	}

	newpolyg->oncf = polyg->oncf;

	return;
}
/*
void construct_polyh(
	CPOLYHEDRON	*polyh,
	CPOLYGON	*polyg,
	CELL		*c)
{
	BFS_add_nb_polygs(polyh,polyg,c);
	return;
}

void BFS_add_nb_polygs(
	CPOLYHEDRON	*polyh,
	CPOLYGON	*polyg,
	CELL		*c)
{
	CPOLYGON *nb_polyg_list;

	find_nb_polygs(&nb_polyg_list,polyg,c);
	//if (nb_polyg_list == NULL)

	return;
}

void find_nb_polygs(
	CPOLYGON	**nb_polyg_list,
	CPOLYGON	*polyg,
	CELL		*c)
{
	CPOINT *p1, *p2;
	CPOLYGON *nb_epolygs;

	p1 = polyg->vertices;
	if (p1 == NULL)
	{
	    printf("ERROR in find_nb_polygs(): empty polygon.\n");
	    clean_up(ERROR);
	}
	p2 = p1->next;

	while (p2)
	{
	    find_polyg_with_edge(&nb_epolygs,p2,p1,c);
	    if (nb_epolygs != NULL)
	    {
		//add nb_polyg
		;
	    }
	    p1 = p2;
	    p2 = p1->next;
	}

	return;
}

void find_polyg_with_edge(
	CPOLYGON	**polygs,
	CPOINT		*p1,
	CPOINT		*p2,
	CELL		*c)
{
	CPOINT *v1, *v2;
	CPOLYGON *cpolyg;

	cpolyg = c->ctri_polygs;
	while (cpolyg)
	{
	    v1 = cpolyg->vertices;
	    v2 = v1->next;
	    if (v1 == NULL)
	    {
		printf("ERROR in find_polyg_with_edge(): empty polygon.\n");
		clean_up(ERROR);
	    }
	    while (v2)
	    {
		if (same_cpt(v1,p1) && same_cpt(v2,p2))
		{
		    //could be more than one
		    //return;
		}
		v1 = v2;
		v2 = v1->next;
	    }
	    //last one
	    v2 = cpolyg->vertices;

	    cpolyg = cpolyg->next;
	}

	//cpolyg = c->cf_polygs;

	return;
}
*/

void G_CARTESIAN::cft_set_cut_cell_vol()
{
	int i, j, k, index;
	double cvol;
	CELL *c;
	CPOLYHEDRON *polyh;

	cvol = top_h[0]*top_h[1]*top_h[2];

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);

	    /*
	    if (c->polyhs == NULL)
		continue;
	    */
	    if (c->polyhs->iscell == TRUE)
	    {
		c->polyhs->vol = cvol;
		continue;
	    }

	    //set areas of polyhedrons' faces
	    set_pf_areas(c);

	    //set volume of polyhedrons
	    set_polyh_vol(c);
	}

	return;
}

void set_pf_areas(CELL *c)
{
	CPOLYHEDRON *polyh;
	CPOLYGON *polyg;

	polyh = c->polyhs;
	while (polyh)
	{
	    polyg = polyh->faces;
	    while (polyg)
	    {
		polyg->area = area_of_polygon(polyg);
		polyg = polyg->next;
	    }
	    polyh = polyh->next;
	}

	return;
}

double area_of_polygon(CPOLYGON *polyg)
{
	int i;
	double area = 0;
	double vcp[3], sumofvcp[3], nor[3];
	CPOINT *v0, *v1, *v2;

	v0 = polyg->vertices;
	v1 = v0->next;
	v2 = v1->next;
	if ((v1 == NULL) || (v2 == NULL))
	{
	    printf("ERROR in area_of_polygon(): incorrect polygon.\n");
	    clean_up(ERROR);
	}
	set_nor(v0,v1,v2,nor);
	for (i = 0; i < 3; i++)
	    sumofvcp[i] = 0.0;
	while (v1)
	{
	    cross_product(v0->crds,v1->crds,vcp);
	    for (i = 0; i < 3; i++)
		sumofvcp[i] += vcp[i];
	    v0 = v1;
	    v1 = v0->next;
	}
	v1 = polyg->vertices;
	cross_product(v0->crds,v1->crds,vcp);
	for (i = 0; i < 3; i++)
	    sumofvcp[i] += vcp[i];

	area = fabs(dot_product(sumofvcp,nor)/2.0);

	return area;
}

void cross_product(
	double	*v0,
	double	*v1,
	double	*p)
{
	p[0] = v0[1]*v1[2] - v0[2]*v1[1];
	p[1] = -v0[0]*v1[2] + v0[2]*v1[0];
	p[2] = v0[0]*v1[1] - v0[1]*v1[0];

	return;
}

double dot_product(
	double	*v0,
	double	*v1)
{
	return (v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]);
}

void set_polyh_vol(CELL *c)
{
	CPOLYHEDRON *polyh;

	polyh = c->polyhs;
	while (polyh)
	{
	    polyh->vol = vol_of_polyh(polyh);

	    polyh = polyh->next;
	}

	return;
}

double vol_of_polyh(CPOLYHEDRON *polyh)
{
	double vol = 0.0;
	CPOLYGON *polyg;

	polyg = polyh->faces;
	while (polyg)
	{
	    set_polyg_nor(polyg);

	    vol += (polyg->area * dot_product(polyg->vertices->crds,polyg->nor))/3.0;

	    polyg = polyg->next;
	}

	return fabs(vol);
}

void set_polyg_nor(CPOLYGON *polyg)
{
	CPOINT *p0, *p1, *p2;

	p0 = polyg->vertices;
	if (p0 == NULL)
	{
	    printf("ERROR in set_polyg_nor().\n");
	    clean_up(ERROR);
	}
	p1 = p0->next;
	if (p1 == NULL)
	{
	    printf("ERROR in set_polyg_nor().\n");
	    clean_up(ERROR);
	}
	p2 = p1->next;
	if (p2 == NULL)
	{
	    printf("ERROR in set_polyg_nor().\n");
	    clean_up(ERROR);
	}

	set_nor(p0,p1,p2,polyg->nor);

	return;
}

//set cell->cut, polyh->comp
void G_CARTESIAN::cft_set_polyhs_comps()
{
	int i, j, k, index;
	COMPONENT comp;
	CELL *c;
	CPOLYHEDRON *polyh;

	setDomain();

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);
	    comp = top_comp[index];
	    polyh = c->polyhs;

	    if (polyh == NULL)
	    {
		printf("ERROR in cft_set_polyhs_comps(): polyh doesn't exist.\n");
		clean_up(ERROR);
	    }

	    if (polyh->iscell)
	    {
		polyh->comp = comp;
		c->cut = FALSE;
		c->comp = comp;
		continue;
	    }
	    
	    c->cut = TRUE;
	    while (polyh)
	    {
		if (polyh->scp_dir == INW)
		{
		    polyh->comp = GAS_COMP2;
		}
		else if (polyh->scp_dir = OUTW)
		{
		    polyh->comp = GAS_COMP1;
		}
		else
		{
		    printf("ERROR in cft_set_polyhs_comps(): unset scp_dir.\n");
		    clean_up(ERROR);
		}
		polyh = polyh->next;
	    }
	}

	return;
}

void init_new_polyg(CPOLYGON *polyg)
{
	polyg->vertices = NULL;
	polyg->edges = NULL;
	polyg->undir_edges = NULL;
	polyg->next = NULL;

	return;
}
