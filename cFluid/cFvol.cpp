#include <cFluid.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

void tri_to_ctri(TRI*,CTRI*);
void copy_pt_to_cpt(POINT*,CPOINT*);
void copy_cpt_to_cpt(CPOINT*,CPOINT*);
void insert_tri_to_ctri_list(TRI*,CELL*);
void set_polygons_in_cell(CELL*);
void init_cf_pts_and_edges(CFACE *cf);
void set_nor(CPOINT*,CPOINT*,CPOINT*,double*);
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
bool same_edge(CEDGE*,CEDGE*);
void print_edge(CEDGE*);
void print_edge_list(CEDGE*);
void print_point(CPOINT*);
void print_polygon(CPOLYGON*);
void construct_polygon_on_tri(CTRI*,CELL*);
void add_polyg_vertex(CPOINT*,CPOLYGON*);
void complete_polyg_undir_edges(CTRI*,CELL*);
void set_pic(CTRI*,CELL*,PIC*);
void add_new_edge(CPOINT*,CPOINT*,CEDGE**);
void add_new_dir_edge(CPOINT*,CPOINT*,CEDGE**);
void add_new_polyg(CPOLYGON*,CPOLYGON**);
void reverse_polyg(CPOLYGON*);
void set_directed_polyg(CPOLYGON*,double*);
void pop_edge(CEDGE**,CEDGE**);
void find_next_endp(CPOINT,CEDGE**,CPOINT*);
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
void complete_polyg_on_cf(CPOLYGON*,CFACE*);

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
	init_grid_cells();

	printf("Init ctris.\n");
        /*Initialize triangles related to cells.*/
	init_tris_in_cells();

	printf("Init polygons.\n");
	set_cell_polygons();

	printf("Calculate volume for polyhedrons.\n");
	//cut_cell_vol();

	////free_ctri_list(cells[i]->ctri_list);
	////free(cells);

	printf("Leave cvol().\n");

	return;
}

void G_CARTESIAN::init_grid_cells()
{
	int i, j, k, l, ll, index;
	CELL *c;

	num_cells = 1;
	for (i = 0; i < dim; i++)
	    num_cells *= (top_gmax[i]+1);
	//Initialize cells in buffer zone but don't use them?	FIXME?
	//start with (lbuf[i] ? lbuf[i] : 1)

	FT_VectorMemoryAlloc((POINTER*)&cells, num_cells, sizeof(CELL));

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    c = &(cells[index]);
	    c->icrds[0] = i;
	    c->icrds[1] = j;
	    c->icrds[2] = k;
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

	    set_polygons_in_cell(c);;
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

	ctri->polyg = (CPOLYGON*)malloc(sizeof(CPOLYGON));

	return;
}

void set_nor(
	CPOINT	*p0,
	CPOINT	*p1,
	CPOINT	*p2,
	double	*nor)
{
	int i;
	double v1[3], v2[3];

	for (i = 0; i < 3; i++)
	{
	    v1[i] = p1->crds[i] - p0->crds[i];
	    v2[i] = p2->crds[i] - p0->crds[i];
	}

	nor[0] = v1[1]*v2[2] - v1[2]*v2[1];
	nor[1] = -v1[0]*v2[2] + v1[2]*v2[0];
	nor[2] = v1[0]*v2[1] - v1[1]*v2[0];
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
	CPOLYGON *cpg;

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
	for (ctri = c->ctris; ctri; ctri = ctri->next)
	{
	    construct_polygon_on_tri(ctri,c);
	    if (ctri->polyg->vertices)
		add_new_polyg(ctri->polyg,&c->ctri_polygs);
	}

	CPOLYGON *ctri_polygs = c->ctri_polygs;
	while (ctri_polygs)
	{
	    //print_polygon(ctri_polygs);
	    ctri_polygs = ctri_polygs->next;
	}

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
	    pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		pedge = (CEDGE*)malloc(sizeof(CEDGE));
		copy_cpt_to_cpt(&(crxp[0]),&(pedge->endp[0]));
		//add_crxp_on_tri_edge(&crxp[0],&(ctri->edges[ei[0]]));
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
		pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		pedge = (CEDGE*)malloc(sizeof(CEDGE));
		copy_cpt_to_cpt(&crxp2,&(pedge->endp[0]));
		copy_cpt_to_cpt(&crxp3,&(pedge->endp[1]));
		found = YES;
	    }
	}

	if (found)
	{
	    add_edge(pedge,&(cf->undirected_edges));
	    add_edge(pedge,&(ctri->polyg->undir_edges));
	    //free(pedge);
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
	    pedge = (CEDGE*)malloc(sizeof(CEDGE));
	    copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
	    copy_cpt_to_cpt(&crxp1,&(pedge->endp[1]));
	    found = YES;
	}
	else if (pocf1 && !pocf2)
	{
	    ncrxps = find_crxps1(&(ctri->pts[i]),&crxp1,cf,&crxp2,&crxp3);
	    if (ncrxps == 1)
	    {
		pedge = (CEDGE*)malloc(sizeof(CEDGE));
		copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
		copy_cpt_to_cpt(&crxp2,&(pedge->endp[1]));
		add_crxp_on_tri_edge(&crxp2,&(ctri->edges[(i+1)%3]));
		found = YES;
	    }
	}
	else if (!pocf1 && pocf2)
	{
	    ncrxps = find_crxps1(&(ctri->pts[i]),&crxp1,cf,&crxp2,&crxp3);
	    if (ncrxps == 1)	//double check	FIXME
	    {
		pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		pedge = (CEDGE*)malloc(sizeof(CEDGE));
		copy_cpt_to_cpt(&crxp2,&(pedge->endp[0]));
		copy_cpt_to_cpt(&crxp3,&(pedge->endp[1]));
		found = YES;
	    }
	}

	if (found)
	{
	    add_edge(pedge,&(cf->undirected_edges));
	    add_edge(pedge,&(ctri->polyg->undir_edges));
	    //free(pedge);
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
		pedge = (CEDGE*)malloc(sizeof(CEDGE));
		copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
		copy_cpt_to_cpt(&(ctri->pts[(i+1)%3]),&(pedge->endp[1]));
		found = YES;
	    }
	    else if (pocf1 && !pocf2)
	    {
		ncrxps = find_crxps1(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),cf,&crxp1,&crxp2);
		if (ncrxps == 1)
		{
		    pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		    pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		    pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		//free(pedge);
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
		pedge = (CEDGE*)malloc(sizeof(CEDGE));
		copy_cpt_to_cpt(&(cf->pts[i]),&(pedge->endp[0]));
		copy_cpt_to_cpt(&(cf->pts[(i+1)%4]),&(pedge->endp[1]));
		found = YES;
	    }
	    else if (pit1 && !pit2)
	    {
		ncrxps = find_crxps2(cf,i,ctri,&crxp1,&crxp2);
		if (ncrxps == 1)
		{
		    pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		    pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		    pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		//free(pedge);
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
		    pedge = (CEDGE*)malloc(sizeof(CEDGE));
		    copy_cpt_to_cpt(&(ctri->pts[i]),&(pedge->endp[0]));
		    copy_cpt_to_cpt(&(ctri->pts[(i+1)%3]),&(pedge->endp[1]));
		    found = YES;
		}
		else if (pocf1 && !pocf2)
		{
		    ncrxps = find_crxps1(&(ctri->pts[i]),&(ctri->pts[(i+1)%3]),cf,&crxp1,&crxp2);
		    if (ncrxps == 1)
		    {
			pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
			pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
			pedge = (CEDGE*)malloc(sizeof(CEDGE));
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
		    //free(pedge);
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
	newcrxp = (CPOINT*)malloc(sizeof(CPOINT));
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
	    if (same_edge(tmp_edge,edge))
		return;
	    tmp_edge = tmp_edge->next;
	}

	new_edge = (CEDGE*)malloc(sizeof(CEDGE));
	copy_cpt_to_cpt(&(edge->endp[0]),&(new_edge->endp[0]));
	copy_cpt_to_cpt(&(edge->endp[1]),&(new_edge->endp[1]));

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
	    if (same_edge(tmpe,edge))
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
	    if (dist1*dist2 < -tol)
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
	    if (dist1*dist2 < -tol)
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
	    if (f1*f2 < tol)
		continue;
	    crxp.crds[edir] = f2/(f1+f2)*tp1->crds[edir] + f1/(f1+f2)*tp2->crds[edir];
	    if ((p1->crds[edir]-crxp.crds[edir])*(crxp.crds[edir]-p2->crds[edir]) > tol)
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

	if (f1*f2 < tol)
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

bool same_edge(
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
/*
void print_point(CPOINT *p)
{
    cout << "("
	 << p->crds[0] << ", "
	 << p->crds[1] << ", "
	 << p->crds[2] << ")";
}

void print_edge(CEDGE *edge)
{
    CPOINT *p1, *p2;

    p1 = &(edge->endp[0]);
    p2 = &(edge->endp[1]);

    cout << "Edge: ";
    print_point(p1);
    cout << ", ";
    print_point(p2);
    cout << "\n";
}

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
		//check direction	TODO
		add_polyg_vertex(&(ctri->pts[i]),ctri->polyg);
	    }
	    return;
	}

	//no polygon on tri
	if (ctri->polyg->undir_edges == NULL)
	    return;

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

	newp = (CPOINT*)malloc(sizeof(CPOINT));
	copy_cpt_to_cpt(p,newp);
	newp->next = pg->vertices;
	pg->vertices = newp;

	return;
}
/*
void print_polygon(CPOLYGON *pg)
{
    CPOINT *p;
    cout << "Polygon: \n";
    p = pg->vertices;
    while(p)
    {
	print_point(p);
	cout << " ";
	p = p->next;
    }
    cout << endl;
}
*/
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
*/
	newedge = (CEDGE*)malloc(sizeof(CEDGE));
	copy_cpt_to_cpt(p0,&(newedge->endp[0]));
	copy_cpt_to_cpt(p1,&(newedge->endp[1]));
	add_edge(newedge,edgelist);

	//free(newedge);
}

void add_new_dir_edge(
	CPOINT	*p0,
	CPOINT	*p1,
	CEDGE	**edgelist)
{
	CEDGE *newedge;

	newedge = (CEDGE*)malloc(sizeof(CEDGE));
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

	//free(newedge);
}

void add_new_polyg(
	CPOLYGON	*polyg,
	CPOLYGON	**polyglist)
{
	CPOINT *p;
	CPOLYGON *newpolyg;

	p = polyg->vertices;
	if (p)
	    newpolyg = (CPOLYGON*)malloc(sizeof(CPOLYGON));
	while (p)
	{
	    add_polyg_vertex(p,newpolyg);
	    p = p->next;
	}
	reverse_polyg(newpolyg);

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

	//edge0 = (CEDGE*)malloc(sizeof(CEDGE));
	pop_edge(&(polyg->undir_edges),&edge0);
	copy_cpt_to_cpt(&(edge0->endp[0]),&p0);
	copy_cpt_to_cpt(&(edge0->endp[1]),&p1);
	find_next_endp(p1,&(polyg->undir_edges),&p2);
	set_nor(&p0,&p1,&p2,pnor);
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
	    p0 = p2;	//FIXME
	    p2 = tmpp;
	}

	find_next_endp(p0,&(polyg->undir_edges),&p1);
	while (!same_cpt(&p1,&p2))
	{
	    add_polyg_vertex(&p1,polyg);
	    p0 = p1;
	    find_next_endp(p0,&(polyg->undir_edges),&p1);
	}

	//free(edge0);
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

void find_next_endp(
	CPOINT	p0,
	CEDGE	**edgelist,
	CPOINT	*p1)
{
	CEDGE *edge, *prev_edge;

	if (edgelist == NULL)
	{
	    printf("Empty edgelist in find_next_endp().\n");
	    clean_up(ERROR);
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
		    //free(edge);
		}
		else
		{
		    prev_edge->next = edge->next;
		    //free(edge);
		}
		return;
	    }
	    else if (same_cpt(&p0,&(edge->endp[1])))
	    {
		copy_cpt_to_cpt(&(edge->endp[0]),p1);
		if (prev_edge == NULL)
		{
		    *edgelist = edge->next;
		    //free(edge);
		}
		else
		{
		    prev_edge->next = edge->next;
		    //free(edge);
		}
		return;
	    }

	    prev_edge = edge;
	    edge = edge->next;
	}

	printf("ERROR: can't find next endp.\n");
	clean_up(ERROR);
}

bool same_nor_dir(
	double	*nor1,
	double	*nor2)
{
	int i;
	double tol = 1e-12;

	for (i = 0; i < 3; i++)
	    if (nor1[i]*nor2[i] < -tol)
		return NO;

	return YES;
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
	pop_edge(&(cf->edges_on_ce),&edge);
	if (!edge)
	{
	    printf("ERROR in set_polygons_on_cf(): no directed edges.\n");
	    clean_up(ERROR);
	}
	while (edge)
	{
	    polyg = (CPOLYGON*)malloc(sizeof(CPOLYGON));
	    add_polyg_vertex(&(edge->endp[1]),polyg);
	    add_polyg_vertex(&(edge->endp[0]),polyg);
	    add_next_vertices_on_ce(polyg,cf);
	    add_prev_vertices_on_ce(polyg,cf);
	    complete_polyg_on_cf(polyg,cf);
	    //cf->edges_in_cf is not //freed after this	FIXME
	    //remove duplicated point	TODO
	    p = polyg->vertices;
	    polyg->vertices = p->next;
	    polyg->next = cell->cf_polygs;
	    cell->cf_polygs = polyg;
	    //free(p);
	    //free(edge);
	    if (!cf->edges_on_ce)
		break;
	    pop_edge(&(cf->edges_on_ce),&edge);
	}

	return;
}

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
	    if (niter > 1000)
	    {
		printf("ERROR in complete_polyg_on_cf(): too many iteration steps.\n");
		clean_up(ERROR);
	    }
	}

	return;
}

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
		    newv = (CPOINT*)malloc(sizeof(CPOINT));
		    copy_cpt_to_cpt(&(edge->endp[1]),newv);
		    endv->next = newv;
		    endv = newv;
		    //delete edge
		    if (!preve)
		    {
			cf->edges_on_ce = edge->next;
			//free(edge);
			found = YES;
			break;
		    }
		    else
		    {
			preve->next = edge->next;
			//free(edge);
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
		    newv = (CPOINT*)malloc(sizeof(CPOINT));
		    copy_cpt_to_cpt(&(edge->endp[0]),newv);
		    newv->next = startv;
		    polyg->vertices = newv;
		    startv = newv;
		    //delete edge
		    if (!preve)
		    {
			cf->edges_on_ce = edge->next;
			//free(edge);
			found = YES;
			break;
		    }
		    else
		    {
			preve->next = edge->next;
			//free(edge);
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

	//set directed edges from undirected edges
	edge = cf->undirected_edges;
	while (edge)
	{
	    edir = edge_dir(edge,cell);
	    for (i = 0; i < 3; i++)
	    {
		midp.crds[i] = (edge->endp[0].crds[i] + edge->endp[1].crds[i])/2.0;
	    }
	    if (point_in_cell_face(&midp,cf))
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

	newp = (CPOINT*)malloc(sizeof(CPOINT));
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
	if (f1 < -tol)
	{
	    printf("ERROR in add_crxp_on_cf_edge_sorted(): negative f1.\n");
	    clean_up(ERROR);
	}
	prevp = NULL;
	while (crxp)
	{
	    f2 = crxp->crds[dir] - cf->pts[index].crds[dir];
	    if (f1 < f2)
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
