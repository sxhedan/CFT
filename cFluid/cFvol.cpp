#include <cFluid.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))
/*
LOCAL void set_nor_of_ctri(CTRI*);
LOCAL bool tri_in_cell(CTRI*,CELL*);
LOCAL void ctri_to_polygon(CTRI*,POLYGON*);
LOCAL bool crxp_to_polygon(CTRI*,POLYGON*,CELL*);
LOCAL void uedges_to_edges(CTRI*,POLYGON*);
LOCAL bool pt_in_cell(CPOINT*,CELL*);
LOCAL bool pt_on_cell_face(CPOINT*,CFACE*);
LOCAL bool edge_on_cell_face(CEDGE*,CFACE*);
LOCAL void set_line_seg_on_cell_face(CEDGE*,CFACE*);
LOCAL void set_ls_on_edges(CFACE*);
LOCAL LINESEG *make_line_seg(CPOINT*,CPOINT*);
LOCAL void set_face_polygons_with_ls(CFACE*);
LOCAL bool point_on_edge(CPOINT*,CEDGE*);
LOCAL void insert_crxp_into_crxplist(CRXPOINT*,CRXPOINT*);
LOCAL void find_tri_edge_crx_with_cell(CEDGE*,CELL*,CTRI*,int*);
LOCAL bool find_cell_edge_crx_with_tri(CEDGE*,CTRI*);
LOCAL bool duplicated_crxp_on_edge(CPOINT*,CEDGE*);
LOCAL int side_of_cell_face(CPOINT*,CFACE*);
LOCAL bool set_crxps_and_edges_201(CTRI*,CFACE*,int*);
LOCAL bool set_crxps_and_edges_111(CTRI*,CFACE*,int*);
LOCAL bool set_crxps_and_edges_120(CTRI*,CFACE*,int*);
LOCAL bool set_crxps_and_edges_030(CTRI*,CFACE*,int*);
LOCAL void find_crxp_on_cf_plane(CPOINT*,CPOINT*,CFACE*,CPOINT*);
LOCAL bool crxp_in_cf(CPOINT*,CFACE*);
LOCAL void find_crxp_in_tri(CTRI*,CFACE*,int*);
LOCAL void find_crxp_on_tri_edge(CEDGE*,CFACE*,int*);
LOCAL bool crxp_on_edges_of_cell_face(CPOINT*,CFACE*);
LOCAL bool duplicated_crxps(CRXPOINT*,CRXPOINT*);
LOCAL void debug_print_cell(CELL*);
LOCAL void debug_print_ctri(CTRI*);
LOCAL void debug_print_cpoint(CPOINT*);
*/

LOCAL void copy_tri_to_ctri(TRI*,CTRI*);
LOCAL void insert_tri_to_ctri_list(TRI*,CELL*);
LOCAL void set_polygons_in_cell(CELL*);

//volume calculation
void G_CARTESIAN::cvol()
{
	printf("Enter cvol().\n");

	if (dim != 3)
	{
	    printf("cvol for 2D is NOT implemented yet.\n");
	    clean_up(ERROR);
	}

	printf("Init cells.\n");
	/*Initialize cell structures.*/
	init_grid_cells();	//We might not need this.	FIXME

	printf("Init ctris.\n");
        /*Initialize triangles related to cells.*/
	init_tris_in_cells();

	printf("Init polygons.\n");
	set_cell_polygons();

	printf("Calculate volume for polyhedrons.\n");
//	cut_cell_vol();

	//free_ctri_list(cells[i]->ctri_list);
	free(cells);

	printf("Leave cvol().\n");

	return;
}

void G_CARTESIAN::init_grid_cells()
{
	int i;

	num_cells = 1;
	for (i = 0; i < dim; i++)
	    num_cells *= (top_gmax[i]+1);
	//Initialize cells in buffer zone but don't use them?	FIXME?
	//start with (lbuf[i] ? lbuf[i] : 1)

	FT_VectorMemoryAlloc((POINTER*)&cells, num_cells, sizeof(CELL));

	return;
}

void G_CARTESIAN::init_tris_in_cells()
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
		printf("Boundary.\n");	//debugdan	FIXME
		continue;
	    }
	    for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	    {
		pt = Point_of_tri(tri);
		pt0 = Coords(pt[0]);
		pt1 = Coords(pt[1]);
		pt2 = Coords(pt[2]);

		//debugdan	FIXME
		printf("pt0: (%lf, %lf, %lf).\n", pt0[0], pt0[1], pt0[2]);
		printf("pt1: (%lf, %lf, %lf).\n", pt1[0], pt1[1], pt1[2]);
		printf("pt2: (%lf, %lf, %lf).\n", pt2[0], pt2[1], pt2[2]);
		//debugdan	FIXME

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

		//debugdan	FIXME
		printf("mincrds: (%lf, %lf, %lf)\n", mincrds[0], mincrds[1], mincrds[2]);
		printf("maxcrds: (%lf, %lf, %lf)\n", maxcrds[0], maxcrds[1], maxcrds[2]);
		printf("minicrds: (%d, %d, %d)\n", minicrds[0], minicrds[1], minicrds[2]);
		printf("maxicrds: (%d, %d, %d)\n", maxicrds[0], maxicrds[1], maxicrds[2]);
		//debugdan	FIXME

		for (k = minicrds[2]; k <= maxicrds[2]; k++)
		for (j = minicrds[1]; j <= maxicrds[1]; j++)
		for (i = minicrds[0]; i <= maxicrds[0]; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    insert_tri_to_ctri_list(tri,&(cells[index]));
		}

		break;
	    }
	    printf("Surf.\n");	//debugdan	FIXME
	}

	return;
}

void G_CARTESIAN::set_cell_polygons()
{
	int i, j, k, index;
	CELL *c;

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[2]; j++)
	for (i = 0; i <= top_gmax[2]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);

	    if (c->num_of_ctris != 0)
	    {
		set_polygons_in_cell(c);;
	    }
	}

	return;
}

//remove this?	Dan	FIXME
LOCAL void copy_tri_to_ctri(TRI *tri, CTRI *ctri)
{
	int i;

	return;
}

LOCAL void insert_tri_to_ctri_list(TRI *tri, CELL *c)
{
	CTRI *ctri;

	FT_ScalarMemoryAlloc((POINTER*)&ctri,sizeof(CTRI));
	//copy_tri_to_ctri(tri, ctri);

	if (c->ctri_list_head)
	{
	    ctri->next = c->ctri_list_head;
	    c->ctri_list_head = ctri;
	}
	else
	    c->ctri_list_head = ctri;

	return;
}

LOCAL void set_polygons_in_cell(CELL *c)
{
	CTRI *ctri;

	c->num_of_ctris = 0;
	for (ctri = c->ctri_list_head; c; ctri=ctri->next)
	{
	    c->num_of_ctris++;
	}

	return;
}
