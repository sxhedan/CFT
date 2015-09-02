/*
 * compute_volume_usingtets.c
 *
 * Code generation for function 'compute_volume_usingtets'
 *
 * C source code generated on: Fri Jul 20 16:35:57 2012
 *
 */

/* Include files */
#include "compute_volume_usingtets.h"

/* Custom Source Code */
#ifdef BUILD_MEX

/* Define macros to support building function into MATLAB executable. */
#include "mex.h"
#define malloc                         mxMalloc
#define calloc                         mxCalloc
#define realloc                        mxRealloc
#define free                           mxFree
#endif

#define M2C_ASSIGN(dst,src)            dst=src

/* Type Definitions */
#ifndef struct_emxArray__common
#define struct_emxArray__common

typedef struct emxArray__common
{
  void *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
} emxArray__common;

#endif

#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T

typedef struct emxArray_boolean_T
{
  boolean_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
} emxArray_boolean_T;

#endif

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void assign_id_basedon_xieta(real_T xi, real_T eta, real_T xp, real_T yp,
  real_T cx, real_T cxmin, real_T cxmax, real_T cy, real_T cymin, real_T cymax,
  real_T faceid, const real_T fcormap[4], const real_T fedgmap[4], real_T p[3],
  real_T id[6]);
static void b_eml_null_assignment(emxArray_real_T *x, int32_T idx);
static boolean_T b_eml_sort_le(const emxArray_int32_T *v, const int32_T col[2],
  int32_T irow1, int32_T irow2);
static void b_emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T
  numDimensions);
static void b_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);
static void b_gridcell_data(const real_T cell_bnds[6], real_T corners[24],
  real_T nrm[18], real_T face_equ[6], int32_T face_equID[6], real_T facebnd[30]);
static void b_mod(const int32_T x[3], real_T y, int32_T r[3]);
static void b_unique(const emxArray_int32_T *a, emxArray_int32_T *b);
static void bitshift(const uint32_T a[3], real_T k, uint32_T c[3]);
static void c_eml_null_assignment(emxArray_real_T *x, const emxArray_int32_T
  *idx);
static void c_emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T
  numDimensions);
static void c_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);
static void c_unique(const int32_T a[6], emxArray_int32_T *b);
extern void call_tetgen(const emxArray_real_T *Fps, const emxArray_int32_T
  *Elems, emxArray_real_T *Fps_new, emxArray_int32_T *Ftris);
extern void call_triangle(const emxArray_real_T *Fps, const emxArray_int32_T
  *Elems, emxArray_real_T *Fps_new, emxArray_int32_T *Ftris);
static void chk_pnt_inside_MB(const real_T cell_bnds[6], const real_T tri[9],
  int32_T flag_v[3]);
static void chk_tri_orientation(const int32_T Tri1[3], const int32_T Tri2[3],
  real_T *SameVert, real_T *SameOrient);
static void compute_intersections(const real_T cell_bnds[6], const
  emxArray_int32_T *opphes, const int32_T trivertID[3], int32_T triID, const
  real_T tri[9], emxArray_real_T *ps_local, emxArray_int32_T *ps_localID,
  int32_T *ps_idx);
static void d_unique(const emxArray_real_T *a, emxArray_real_T *b);
static void determine_opposite_halfedge_tri(int32_T nv, const emxArray_int32_T
  *tris, emxArray_int32_T *opphes);
static void determine_opposite_halfface(int32_T nv, const emxArray_int32_T
  *elems, emxArray_int32_T *opphfs);
static void determine_opposite_halfface_tet(int32_T nv, const emxArray_int32_T
  *elems, emxArray_int32_T *opphfs);
static int32_T div_s32_floor(int32_T numerator, int32_T denominator);
static void e_unique(const emxArray_int32_T *a, emxArray_int32_T *b);
static void eml_null_assignment(emxArray_int32_T *x, int32_T idx);
static boolean_T eml_sort_le(const emxArray_int32_T *v, int32_T col, int32_T
  irow1, int32_T irow2);
static real_T eml_xdot(int32_T n, const real_T x[3], int32_T ix0, int32_T incx,
  const real_T y[3], int32_T iy0, int32_T incy);
static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize);
static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int32_T
  numDimensions);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);
static void extract_cellface_edges(real_T FaceId, const int32_T Faces[4], const
  int32_T FaceEdgeMap[4], const emxArray_int32_T *cor_map, int32_T ps_index,
  const emxArray_real_T *ps_local, const emxArray_int32_T *ps_localID, const
  emxArray_int32_T *tris_labelled, const emxArray_int32_T *opphes,
  emxArray_int32_T *Face_ps, emxArray_int32_T *Face_elems);
static int32_T extract_trisedges_on_cellface(const emxArray_int32_T *Face_ps,
  int32_T npnts, const emxArray_int32_T *tris, const emxArray_int32_T *opphes,
  emxArray_int32_T *Face_elems);
static void get_celledgelems_eight(const int32_T Faces[4], const int32_T
  FaceEdgeMap[4], int32_T b_index, const emxArray_int32_T *Face_ps,
  emxArray_int32_T *Face_elems, int32_T *fe_idx, int32_T ps_index, const
  emxArray_real_T *ps_local, const emxArray_int32_T *ps_localID);
static void get_celledgelems_noteight(const int32_T Faces[4], const int32_T
  FaceEdgeMap[4], const emxArray_int32_T *cor_map, int32_T b_index, const
  emxArray_int32_T *Face_ps, emxArray_int32_T *Face_elems, int32_T *fe_idx,
  int32_T ps_index, const emxArray_real_T *ps_local, const emxArray_int32_T
  *ps_localID);
static void get_pntlist(const real_T cell_bnds[6], const emxArray_int32_T
  *opphes, int32_T triID, const int32_T trivertID[3], const real_T tri[9], const
  int32_T flag_v[3], emxArray_real_T *ps_local, emxArray_int32_T *ps_localID,
  int32_T *ps_idx);
static void get_trisincell(const real_T cell_bnds[6], const emxArray_real_T *ps,
  int32_T num_tris, const emxArray_int32_T *tris, const emxArray_int32_T
  *comp_index, emxArray_int32_T *tris_cell, emxArray_int32_T *comp_index_cell,
  int32_T *flag_cell, int32_T *tris_idx);
static void gridcell_data(const real_T cell_bnds[6], real_T corners[24], real_T
  nrm[18], real_T face_equ[6], int32_T face_equID[6]);
static void intersection_triedge_cellface_2(const real_T cell_bnds[6], int32_T
  triID, real_T edgID, const int32_T edgvertID[2], const real_T edgps[6],
  int32_T *ps_idx, emxArray_real_T *ps_local, emxArray_int32_T *ps_localID);
static void label_cellface(const emxArray_int32_T *ps_localID, const
  emxArray_int32_T *tris_local, const emxArray_int32_T *tris_bndID, const
  emxArray_int32_T *tets, const emxArray_int32_T *TETcomp, const
  emxArray_int32_T *TETsign, int32_T num_comps, const emxArray_int32_T
  *comp_list, const emxArray_int32_T *pos_map, const emxArray_int32_T *neg_map,
  int32_T face_type[6], boolean_T face_sign[6]);
static void label_tets(int32_T nv, const emxArray_int32_T *tets, const
  emxArray_int32_T *tris_labelled, const emxArray_int32_T *tris_comp, int32_T
  num_comps, const emxArray_int32_T *comp_list, emxArray_int32_T *pos_map,
  emxArray_int32_T *neg_map, emxArray_int32_T *TetComp, emxArray_int32_T
  *TetSign);
static real_T nnz(const real_T s[4]);
static void nonzeros(const emxArray_int32_T *s, emxArray_int32_T *v);
static int32_T re_triangulate(const emxArray_real_T *ps, const emxArray_int32_T *
  tris_cell, const emxArray_int32_T *opphes, int32_T tris_idx, const
  emxArray_int32_T *comp_index_cell, const emxArray_real_T *ps_local, const
  emxArray_int32_T *ps_localID, emxArray_int32_T *tris_labelled,
  emxArray_int32_T *tris_comp);
static void searchstar(int32_T key, int32_T idx, int32_T num_comps, const
  emxArray_int32_T *comp_list, const emxArray_int32_T *map, emxArray_real_T *L,
  real_T *p);
static void tetrahedralize_cell(const real_T cell_bnds[6], const emxArray_real_T
  *ps, int32_T num_points, const emxArray_int32_T *tris_cell, int32_T tris_idx,
  const emxArray_int32_T *comp_index_cell, emxArray_real_T *ps_vol,
  emxArray_int32_T *tets, emxArray_int32_T *tris_labelled, emxArray_int32_T
  *tris_comp, emxArray_int32_T *tris_local, emxArray_int32_T *tris_bndID,
  emxArray_int32_T *ps_localID, real_T *flag_tri, int32_T *nv);
static void triangulate_cellface(const real_T Face_nrm[18], const real_T
  Face_equ[6], const int32_T Face_equID[6], const emxArray_int32_T *cor_map,
  int32_T ps_index, emxArray_real_T *ps_local, const emxArray_int32_T
  *ps_localID, int32_T tris_index, emxArray_int32_T *tris_local,
  emxArray_int32_T *tris_bndID, const emxArray_int32_T *tris_labelled);
static void unique(const emxArray_int32_T *a, emxArray_int32_T *b);
static void update_map_component(const emxArray_int32_T *comp, int32_T num_comps,
  const emxArray_int32_T *comp_list, emxArray_int32_T *map);

/* Function Definitions */

/*
 * function [p,id] = assign_id_basedon_xieta(xi,eta,xp,yp,cx,cxmin,cxmax,cy,cymin,cymax,faceid,fcormap,fedgmap,p,id)
 */
static void assign_id_basedon_xieta(real_T xi, real_T eta, real_T xp, real_T yp,
  real_T cx, real_T cxmin, real_T cxmax, real_T cy, real_T cymin, real_T cymax,
  real_T faceid, const real_T fcormap[4], const real_T fedgmap[4], real_T p[3],
  real_T id[6])
{
  /* 'assign_id_basedon_xieta:2' if (xi ==0 && eta == 0) */
  if ((xi == 0.0) && (eta == 0.0)) {
    /* 'assign_id_basedon_xieta:3' p(cx) = cxmin; */
    p[(int32_T)cx - 1] = cxmin;

    /* 'assign_id_basedon_xieta:4' p(cy) = cymin; */
    p[(int32_T)cy - 1] = cymin;

    /* 'assign_id_basedon_xieta:5' id(6) = fcormap(1); */
    id[5] = fcormap[0];
  } else if ((xi == 1.0) && (eta == 0.0)) {
    /* 'assign_id_basedon_xieta:6' elseif (xi ==1 && eta == 0) */
    /* 'assign_id_basedon_xieta:7' p(cx) = cxmax; */
    p[(int32_T)cx - 1] = cxmax;

    /* 'assign_id_basedon_xieta:8' p(cy) = cymin; */
    p[(int32_T)cy - 1] = cymin;

    /* 'assign_id_basedon_xieta:9' id(6) = fcormap(2); */
    id[5] = fcormap[1];
  } else if ((xi == 1.0) && (eta == 1.0)) {
    /* 'assign_id_basedon_xieta:10' elseif (xi ==1 && eta == 1) */
    /* 'assign_id_basedon_xieta:11' p(cx) = cxmax; */
    p[(int32_T)cx - 1] = cxmax;

    /* 'assign_id_basedon_xieta:12' p(cy) = cymax; */
    p[(int32_T)cy - 1] = cymax;

    /* 'assign_id_basedon_xieta:13' id(6) = fcormap(3); */
    id[5] = fcormap[2];
  } else if ((xi == 0.0) && (eta == 1.0)) {
    /* 'assign_id_basedon_xieta:14' elseif (xi ==0 && eta == 1) */
    /* 'assign_id_basedon_xieta:15' p(cx) = cxmin; */
    p[(int32_T)cx - 1] = cxmin;

    /* 'assign_id_basedon_xieta:16' p(cy) = cymax; */
    p[(int32_T)cy - 1] = cymax;

    /* 'assign_id_basedon_xieta:17' id(6) =fcormap(4); */
    id[5] = fcormap[3];
  } else if ((xi > 0.0) && (xi < 1.0) && (eta == 0.0)) {
    /* 'assign_id_basedon_xieta:18' elseif (xi >0 && xi < 1 && eta == 0) */
    /* 'assign_id_basedon_xieta:19' p(cx) = xp; */
    p[(int32_T)cx - 1] = xp;

    /* 'assign_id_basedon_xieta:20' p(cy) = cymin; */
    p[(int32_T)cy - 1] = cymin;

    /* 'assign_id_basedon_xieta:21' id(5) = fedgmap(1); */
    id[4] = fedgmap[0];
  } else if ((xi == 1.0) && (eta > 0.0) && (eta < 1.0)) {
    /* 'assign_id_basedon_xieta:22' elseif (xi ==1 && eta > 0 && eta <1) */
    /* 'assign_id_basedon_xieta:23' p(cx) = cxmax; */
    p[(int32_T)cx - 1] = cxmax;

    /* 'assign_id_basedon_xieta:24' p(cy) = yp; */
    p[(int32_T)cy - 1] = yp;

    /* 'assign_id_basedon_xieta:25' id(5) = fedgmap(2); */
    id[4] = fedgmap[1];
  } else if ((xi > 0.0) && (xi < 1.0) && (eta == 1.0)) {
    /* 'assign_id_basedon_xieta:26' elseif (xi >0 && xi <1 && eta == 1) */
    /* 'assign_id_basedon_xieta:27' p(cx) = xp; */
    p[(int32_T)cx - 1] = xp;

    /* 'assign_id_basedon_xieta:28' p(cy) = cymax; */
    p[(int32_T)cy - 1] = cymax;

    /* 'assign_id_basedon_xieta:29' id(5) = fedgmap(3); */
    id[4] = fedgmap[2];
  } else if ((xi == 0.0) && (eta > 0.0) && (eta < 1.0)) {
    /* 'assign_id_basedon_xieta:30' elseif (xi ==0 && eta >0 && eta<1) */
    /* 'assign_id_basedon_xieta:31' p(cx) = cxmin; */
    p[(int32_T)cx - 1] = cxmin;

    /* 'assign_id_basedon_xieta:32' p(cy) = yp; */
    p[(int32_T)cy - 1] = yp;

    /* 'assign_id_basedon_xieta:33' id(5) = fedgmap(4); */
    id[4] = fedgmap[3];
  } else {
    if ((xi > 0.0) && (xi < 1.0) && (eta > 0.0) && (eta < 1.0)) {
      /* 'assign_id_basedon_xieta:34' elseif (xi >0 && xi <1 && eta >0 && eta<1) */
      /* 'assign_id_basedon_xieta:35' p(cx)=xp; */
      p[(int32_T)cx - 1] = xp;

      /* 'assign_id_basedon_xieta:36' p(cy)=yp; */
      p[(int32_T)cy - 1] = yp;

      /* 'assign_id_basedon_xieta:37' id(4) = faceid; */
      id[3] = faceid;
    }
  }
}

/*
 *
 */
static void b_eml_null_assignment(emxArray_real_T *x, int32_T idx)
{
  int32_T nrows;
  int32_T j;
  int32_T loop_ub;
  emxArray_real_T *b_x;
  int32_T i15;
  nrows = x->size[0] - 1;
  for (j = 0; j < 3; j++) {
    for (loop_ub = idx; loop_ub <= nrows; loop_ub++) {
      x->data[(loop_ub + x->size[0] * j) - 1] = x->data[loop_ub + x->size[0] * j];
    }
  }

  if (1 > nrows) {
    nrows = 0;
  }

  emxInit_real_T(&b_x, 2);
  j = b_x->size[0] * b_x->size[1];
  b_x->size[0] = nrows;
  b_x->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)b_x, j, (int32_T)sizeof(real_T));
  for (j = 0; j < 3; j++) {
    loop_ub = nrows - 1;
    for (i15 = 0; i15 <= loop_ub; i15++) {
      b_x->data[i15 + b_x->size[0] * j] = x->data[i15 + x->size[0] * j];
    }
  }

  j = x->size[0] * x->size[1];
  x->size[0] = b_x->size[0];
  x->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)x, j, (int32_T)sizeof(real_T));
  for (j = 0; j < 3; j++) {
    loop_ub = b_x->size[0] - 1;
    for (i15 = 0; i15 <= loop_ub; i15++) {
      x->data[i15 + x->size[0] * j] = b_x->data[i15 + b_x->size[0] * j];
    }
  }

  emxFree_real_T(&b_x);
}

/*
 *
 */
static boolean_T b_eml_sort_le(const emxArray_int32_T *v, const int32_T col[2],
  int32_T irow1, int32_T irow2)
{
  boolean_T p;
  int32_T k;
  boolean_T exitg1;
  int32_T abscolk;
  boolean_T b2;
  p = TRUE;
  k = 0;
  exitg1 = 0U;
  while ((exitg1 == 0U) && (k < 2)) {
    if (col[k] < 0) {
      abscolk = -col[k];
    } else {
      abscolk = col[k];
    }

    abscolk = (abscolk - 1) * v->size[0] - 1;
    if (v->data[abscolk + irow1] == v->data[abscolk + irow2]) {
      b2 = TRUE;
    } else {
      b2 = FALSE;
    }

    if (!b2) {
      if (col[k] < 0) {
        if (v->data[abscolk + irow1] >= v->data[abscolk + irow2]) {
          p = TRUE;
        } else {
          p = FALSE;
        }
      } else if (v->data[abscolk + irow1] <= v->data[abscolk + irow2]) {
        p = TRUE;
      } else {
        p = FALSE;
      }

      exitg1 = 1U;
    } else {
      k++;
    }
  }

  return p;
}

static void b_emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T
  numDimensions)
{
  emxArray_int32_T *emxArray;
  int32_T loop_ub;
  int32_T i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int32_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    emxArray->size[i] = 0;
  }
}

static void b_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions)
{
  emxArray_real_T *emxArray;
  int32_T loop_ub;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * function [faces,corners,nrm,face_equ,face_equID,facebnd,fbmap,fcormap,fedgmap] = gridcell_data(cell_bnds)
 *  Bounds of the cell
 */
static void b_gridcell_data(const real_T cell_bnds[6], real_T corners[24],
  real_T nrm[18], real_T face_equ[6], int32_T face_equID[6], real_T facebnd[30])
{
  int32_T i;
  static const int8_T iv9[3] = { 0, 0, -1 };

  static const int8_T iv10[3] = { 1, 0, 0 };

  static const int8_T iv11[3] = { 0, 1, 0 };

  static const int8_T iv12[3] = { -1, 0, 0 };

  static const int8_T iv13[3] = { 0, -1, 0 };

  static const int8_T iv14[3] = { 0, 0, 1 };

  /* 'gridcell_data:4' xmin = cell_bnds(1); */
  /* 'gridcell_data:5' xmax = cell_bnds(2); */
  /* 'gridcell_data:6' ymin = cell_bnds(3); */
  /* 'gridcell_data:7' ymax = cell_bnds(4); */
  /* 'gridcell_data:8' zmin = cell_bnds(5); */
  /* 'gridcell_data:9' zmax = cell_bnds(6); */
  /*  Face-Corner Map */
  /* 'gridcell_data:12' faces = int32([1 4 3 2;1 2 6 5;2 3 7 6;3 4 8 7;1 5 8 4;5 6 7 8]); */
  /*  Corners of the cell */
  /* 'gridcell_data:15' corners = zeros(8,3); */
  for (i = 0; i < 24; i++) {
    corners[i] = 0.0;
  }

  /* 'gridcell_data:16' corners(1,:) = [xmax ymin zmin]; */
  corners[0] = cell_bnds[1];
  corners[8] = cell_bnds[2];
  corners[16] = cell_bnds[4];

  /* 'gridcell_data:17' corners(2,:) = [xmax ymax zmin]; */
  corners[1] = cell_bnds[1];
  corners[9] = cell_bnds[3];
  corners[17] = cell_bnds[4];

  /* 'gridcell_data:18' corners(3,:) = [xmin ymax zmin]; */
  corners[2] = cell_bnds[0];
  corners[10] = cell_bnds[3];
  corners[18] = cell_bnds[4];

  /* 'gridcell_data:19' corners(4,:) = [xmin ymin zmin]; */
  corners[3] = cell_bnds[0];
  corners[11] = cell_bnds[2];
  corners[19] = cell_bnds[4];

  /* 'gridcell_data:20' corners(5,:) = [xmax ymin zmax]; */
  corners[4] = cell_bnds[1];
  corners[12] = cell_bnds[2];
  corners[20] = cell_bnds[5];

  /* 'gridcell_data:21' corners(6,:) = [xmax ymax zmax]; */
  corners[5] = cell_bnds[1];
  corners[13] = cell_bnds[3];
  corners[21] = cell_bnds[5];

  /* 'gridcell_data:22' corners(7,:) = [xmin ymax zmax]; */
  corners[6] = cell_bnds[0];
  corners[14] = cell_bnds[3];
  corners[22] = cell_bnds[5];

  /* 'gridcell_data:23' corners(8,:) = [xmin ymin zmax]; */
  corners[7] = cell_bnds[0];
  corners[15] = cell_bnds[2];
  corners[23] = cell_bnds[5];

  /*  Normals to the cell faces */
  /* 'gridcell_data:26' nrm = zeros(6,3); */
  for (i = 0; i < 18; i++) {
    nrm[i] = 0.0;
  }

  /* 'gridcell_data:27' nrm(1,:) = [0 0 -1]; */
  for (i = 0; i < 3; i++) {
    nrm[6 * i] = (real_T)iv9[i];

    /* 'gridcell_data:28' nrm(2,:) = [1 0 0]; */
    nrm[1 + 6 * i] = (real_T)iv10[i];

    /* 'gridcell_data:29' nrm(3,:) = [0 1 0]; */
    nrm[2 + 6 * i] = (real_T)iv11[i];

    /* 'gridcell_data:30' nrm(4,:) = [-1 0 0]; */
    nrm[3 + 6 * i] = (real_T)iv12[i];

    /* 'gridcell_data:31' nrm(5,:) = [0 -1 0]; */
    nrm[4 + 6 * i] = (real_T)iv13[i];

    /* 'gridcell_data:32' nrm(6,:) = [0 0 1]; */
    nrm[5 + 6 * i] = (real_T)iv14[i];
  }

  /*  Equation of the cell faces */
  /* 'gridcell_data:35' face_equ = zeros(6,1); */
  /* 'gridcell_data:36' face_equ(1) = zmin; */
  for (i = 0; i < 6; i++) {
    face_equ[i] = 0.0;

    /*  ID of the face to distinguish if it is the x/y/z-plane */
    /* 'gridcell_data:44' face_equID = zeros(6,1,'int32'); */
    face_equID[i] = 0;
  }

  face_equ[0] = cell_bnds[4];

  /* 'gridcell_data:37' face_equ(2) = xmax; */
  face_equ[1] = cell_bnds[1];

  /* 'gridcell_data:38' face_equ(3) = ymax; */
  face_equ[2] = cell_bnds[3];

  /* 'gridcell_data:39' face_equ(4) = xmin; */
  face_equ[3] = cell_bnds[0];

  /* 'gridcell_data:40' face_equ(5) = ymin; */
  face_equ[4] = cell_bnds[2];

  /* 'gridcell_data:41' face_equ(6) = zmax; */
  face_equ[5] = cell_bnds[5];

  /* 'gridcell_data:45' face_equID(1) = 3; */
  face_equID[0] = 3;

  /* 'gridcell_data:46' face_equID(2) = 1; */
  face_equID[1] = 1;

  /* 'gridcell_data:47' face_equID(3) = 2; */
  face_equID[2] = 2;

  /* 'gridcell_data:48' face_equID(4) = 1; */
  face_equID[3] = 1;

  /* 'gridcell_data:49' face_equID(5) = 2; */
  face_equID[4] = 2;

  /* 'gridcell_data:50' face_equID(6) = 3; */
  face_equID[5] = 3;

  /*  Face bounds and map with coordinate for each cell face. */
  /*  For a typical face, */
  /*  fbmap(1),fbmap(2) : The coordinates that vary for the face, their bounds are */
  /*  given by facebnd(1:2),facebnd(3:4), respectively. fbmap(3) is the */
  /*  coordinate that remains constant and this constant value is facebnd(5). */
  /* 'gridcell_data:58' facebnd = [xmin xmax ymin ymax zmin; */
  /* 'gridcell_data:59'     ymin ymax zmin zmax xmax; */
  /* 'gridcell_data:60'     xmin xmax zmin zmax ymax; */
  /* 'gridcell_data:61'     ymin ymax zmin zmax xmin; */
  /* 'gridcell_data:62'     xmin xmax zmin zmax ymin; */
  /* 'gridcell_data:63'     xmin xmax ymin ymax zmax]; */
  facebnd[0] = cell_bnds[0];
  facebnd[6] = cell_bnds[1];
  facebnd[12] = cell_bnds[2];
  facebnd[18] = cell_bnds[3];
  facebnd[24] = cell_bnds[4];
  facebnd[1] = cell_bnds[2];
  facebnd[7] = cell_bnds[3];
  facebnd[13] = cell_bnds[4];
  facebnd[19] = cell_bnds[5];
  facebnd[25] = cell_bnds[1];
  facebnd[2] = cell_bnds[0];
  facebnd[8] = cell_bnds[1];
  facebnd[14] = cell_bnds[4];
  facebnd[20] = cell_bnds[5];
  facebnd[26] = cell_bnds[3];
  facebnd[3] = cell_bnds[2];
  facebnd[9] = cell_bnds[3];
  facebnd[15] = cell_bnds[4];
  facebnd[21] = cell_bnds[5];
  facebnd[27] = cell_bnds[0];
  facebnd[4] = cell_bnds[0];
  facebnd[10] = cell_bnds[1];
  facebnd[16] = cell_bnds[4];
  facebnd[22] = cell_bnds[5];
  facebnd[28] = cell_bnds[2];
  facebnd[5] = cell_bnds[0];
  facebnd[11] = cell_bnds[1];
  facebnd[17] = cell_bnds[2];
  facebnd[23] = cell_bnds[3];
  facebnd[29] = cell_bnds[5];

  /* 'gridcell_data:65' fbmap =[1 2 3; 2 3 1; 1 3 2; 2 3 1;1 3 2;1 2 3]; */
  /*  Mapping of the barycentric coordinates with the corners and cell edges. */
  /*  Corners: fcormap(1:4) corresponds to (xi,eta) = (0,0),(1,0),(1,1),(0,1) */
  /*  Edges : fedgmap(1:4) corresponds to(xi,eta) = ((0,1),0), (1,(0,1)),((0,1),1), (0,(0,1)) */
  /* 'gridcell_data:70' fcormap =[4 1 2 3; 1 2 6 5; 3 2 6 7; 4 3 7 8; 4 1 5 8; 8 5 6 7]; */
  /* 'gridcell_data:71' fedgmap = [4 1 2 3; 1 6 9 5; 2 6 10 7; 3 7 11 8; 4 5 12 8; 12 9 10 11]; */
}

/*
 *
 */
static void b_mod(const int32_T x[3], real_T y, int32_T r[3])
{
  int32_T k;
  real_T d1;
  int32_T yk;
  for (k = 0; k < 3; k++) {
    d1 = y;
    d1 = d1 < 0.0 ? ceil(d1 - 0.5) : floor(d1 + 0.5);
    yk = (int32_T)d1;
    if (yk == 0) {
      yk = x[k];
    } else {
      yk = x[k] - div_s32_floor(x[k], yk) * yk;
    }

    r[k] = yk;
  }
}

/*
 *
 */
static void b_unique(const emxArray_int32_T *a, emxArray_int32_T *b)
{
  int32_T p;
  int32_T k0;
  int32_T col[2];
  int32_T k;
  emxArray_int32_T *idx;
  int32_T n;
  int32_T np1;
  emxArray_int32_T *idx0;
  int32_T i;
  int32_T m;
  int32_T j;
  int32_T nb;
  int32_T qEnd;
  int32_T kEnd;
  int32_T exitg1;
  boolean_T b_p;
  boolean_T exitg2;
  emxArray_int32_T *b_b;
  if (a->size[0] == 0) {
    p = b->size[0] * b->size[1];
    b->size[0] = a->size[0];
    b->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)b, p, (int32_T)sizeof(int32_T));
    k0 = a->size[0] * a->size[1] - 1;
    for (p = 0; p <= k0; p++) {
      b->data[p] = a->data[p];
    }
  } else {
    for (k = 0; k < 2; k++) {
      col[k] = k + 1;
    }

    b_emxInit_int32_T(&idx, 1);
    n = a->size[0];
    p = idx->size[0];
    idx->size[0] = n;
    emxEnsureCapacity((emxArray__common *)idx, p, (int32_T)sizeof(int32_T));
    np1 = n + 1;
    for (k = 1; k <= n; k++) {
      idx->data[k - 1] = k;
    }

    for (k = 1; k <= n - 1; k += 2) {
      if (b_eml_sort_le(a, col, k, k + 1)) {
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    b_emxInit_int32_T(&idx0, 1);
    p = idx0->size[0];
    idx0->size[0] = n;
    emxEnsureCapacity((emxArray__common *)idx0, p, (int32_T)sizeof(int32_T));
    k0 = n - 1;
    for (p = 0; p <= k0; p++) {
      idx0->data[p] = 1;
    }

    i = 2;
    while (i < n) {
      m = i << 1;
      j = 1;
      for (k0 = 1 + i; k0 < np1; k0 = qEnd + i) {
        p = j;
        nb = k0;
        qEnd = j + m;
        if (qEnd > np1) {
          qEnd = np1;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          if (b_eml_sort_le(a, col, idx->data[p - 1], idx->data[nb - 1])) {
            idx0->data[k] = idx->data[p - 1];
            p++;
            if (p == k0) {
              while (nb < qEnd) {
                k++;
                idx0->data[k] = idx->data[nb - 1];
                nb++;
              }
            }
          } else {
            idx0->data[k] = idx->data[nb - 1];
            nb++;
            if (nb == qEnd) {
              while (p < k0) {
                k++;
                idx0->data[k] = idx->data[p - 1];
                p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k + 1 <= kEnd; k++) {
          idx->data[(j + k) - 1] = idx0->data[k];
        }

        j = qEnd;
      }

      i = m;
    }

    p = b->size[0] * b->size[1];
    b->size[0] = a->size[0];
    b->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)b, p, (int32_T)sizeof(int32_T));
    k0 = a->size[0] * a->size[1] - 1;
    for (p = 0; p <= k0; p++) {
      b->data[p] = a->data[p];
    }

    m = a->size[0];
    col[0] = m;
    col[1] = 1;
    p = idx0->size[0];
    idx0->size[0] = col[0];
    emxEnsureCapacity((emxArray__common *)idx0, p, (int32_T)sizeof(int32_T));
    for (j = 0; j < 2; j++) {
      for (i = 0; i + 1 <= m; i++) {
        idx0->data[i] = b->data[(idx->data[i] + b->size[0] * j) - 1];
      }

      for (i = 0; i + 1 <= m; i++) {
        b->data[i + b->size[0] * j] = idx0->data[i];
      }
    }

    emxFree_int32_T(&idx0);
    emxFree_int32_T(&idx);
    nb = 0;
    m = a->size[0];
    k = 0;
    while (k + 1 <= m) {
      k0 = k;
      do {
        exitg1 = 0U;
        k++;
        if (k + 1 > m) {
          exitg1 = 1U;
        } else {
          b_p = FALSE;
          j = 0;
          exitg2 = 0U;
          while ((exitg2 == 0U) && (j < 2)) {
            if (!(b->data[k0 + b->size[0] * j] == b->data[k + b->size[0] * j]))
            {
              b_p = TRUE;
              exitg2 = 1U;
            } else {
              j++;
            }
          }

          if (b_p) {
            exitg1 = 1U;
          }
        }
      } while (exitg1 == 0U);

      nb++;
      for (j = 0; j < 2; j++) {
        b->data[(nb + b->size[0] * j) - 1] = b->data[k0 + b->size[0] * j];
      }
    }

    if (1 > nb) {
      nb = 0;
    }

    emxInit_int32_T(&b_b, 2);
    p = b_b->size[0] * b_b->size[1];
    b_b->size[0] = nb;
    b_b->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)b_b, p, (int32_T)sizeof(int32_T));
    for (p = 0; p < 2; p++) {
      k0 = nb - 1;
      for (m = 0; m <= k0; m++) {
        b_b->data[m + b_b->size[0] * p] = b->data[m + b->size[0] * p];
      }
    }

    p = b->size[0] * b->size[1];
    b->size[0] = b_b->size[0];
    b->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)b, p, (int32_T)sizeof(int32_T));
    for (p = 0; p < 2; p++) {
      k0 = b_b->size[0] - 1;
      for (m = 0; m <= k0; m++) {
        b->data[m + b->size[0] * p] = b_b->data[m + b_b->size[0] * p];
      }
    }

    emxFree_int32_T(&b_b);
  }
}

/*
 *
 */
static void bitshift(const uint32_T a[3], real_T k, uint32_T c[3])
{
  int32_T m;
  real_T d0;
  uint8_T absk1;
  for (m = 0; m < 3; m++) {
    c[m] = 0U;
    if (k < 0.0) {
      d0 = -k;
      d0 = floor(d0 + 0.5);
      absk1 = (uint8_T)d0;
      if (absk1 < 32) {
        c[m] = a[m] >> (uint32_T)absk1;
      }
    } else {
      d0 = k;
      d0 = floor(d0 + 0.5);
      absk1 = (uint8_T)d0;
      if (absk1 < 32) {
        c[m] = a[m] << (uint32_T)absk1;
      }
    }
  }
}

/*
 *
 */
static void c_eml_null_assignment(emxArray_real_T *x, const emxArray_int32_T
  *idx)
{
  int32_T nrows;
  int32_T j;
  int32_T i;
  emxArray_boolean_T *b;
  int32_T i16;
  int32_T nb;
  int32_T k;
  emxArray_real_T *b_x;
  if (idx->size[1] == 1) {
    nrows = x->size[0] - 1;
    for (j = 0; j < 3; j++) {
      for (i = idx->data[0]; i <= nrows; i++) {
        x->data[(i + x->size[0] * j) - 1] = x->data[i + x->size[0] * j];
      }
    }
  } else {
    emxInit_boolean_T(&b, 2);
    i16 = b->size[0] * b->size[1];
    b->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)b, i16, (int32_T)sizeof(boolean_T));
    nb = x->size[0];
    i16 = b->size[0] * b->size[1];
    b->size[1] = nb;
    emxEnsureCapacity((emxArray__common *)b, i16, (int32_T)sizeof(boolean_T));
    nb = x->size[0] - 1;
    for (i16 = 0; i16 <= nb; i16++) {
      b->data[i16] = FALSE;
    }

    for (k = 1; k <= idx->size[1]; k++) {
      b->data[idx->data[k - 1] - 1] = TRUE;
    }

    nb = 0;
    for (k = 1; k <= b->size[1]; k++) {
      i = b->data[k - 1];
      nb += i;
    }

    nrows = x->size[0] - nb;
    nb = b->size[1];
    i = 0;
    i16 = x->size[0];
    for (k = 1; k <= i16; k++) {
      if ((k > nb) || (!b->data[k - 1])) {
        for (j = 0; j < 3; j++) {
          x->data[i + x->size[0] * j] = x->data[(k + x->size[0] * j) - 1];
        }

        i++;
      }
    }

    emxFree_boolean_T(&b);
  }

  if (1 > nrows) {
    nrows = 0;
  }

  emxInit_real_T(&b_x, 2);
  i16 = b_x->size[0] * b_x->size[1];
  b_x->size[0] = nrows;
  b_x->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)b_x, i16, (int32_T)sizeof(real_T));
  for (i16 = 0; i16 < 3; i16++) {
    nb = nrows - 1;
    for (i = 0; i <= nb; i++) {
      b_x->data[i + b_x->size[0] * i16] = x->data[i + x->size[0] * i16];
    }
  }

  i16 = x->size[0] * x->size[1];
  x->size[0] = b_x->size[0];
  x->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)x, i16, (int32_T)sizeof(real_T));
  for (i16 = 0; i16 < 3; i16++) {
    nb = b_x->size[0] - 1;
    for (i = 0; i <= nb; i++) {
      x->data[i + x->size[0] * i16] = b_x->data[i + b_x->size[0] * i16];
    }
  }

  emxFree_real_T(&b_x);
}

static void c_emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T
  numDimensions)
{
  emxArray_int32_T *emxArray;
  int32_T loop_ub;
  int32_T i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int32_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    emxArray->size[i] = 0;
  }
}

static void c_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions)
{
  emxArray_real_T *emxArray;
  int32_T loop_ub;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 *
 */
static void c_unique(const int32_T a[6], emxArray_int32_T *b)
{
  int8_T idx[6];
  int32_T k;
  boolean_T p;
  int8_T idx0[6];
  int32_T i;
  int32_T nb;
  int32_T j;
  int32_T pEnd;
  int32_T b_p;
  int32_T q;
  int32_T qEnd;
  int32_T kEnd;
  emxArray_int32_T *r5;
  emxArray_int32_T *b_b;
  int32_T iv16[2];
  emxArray_int32_T r6;
  for (k = 0; k < 6; k++) {
    idx[k] = (int8_T)(k + 1);
  }

  for (k = 0; k < 6; k += 2) {
    if (a[k] <= a[k + 1]) {
      p = TRUE;
    } else {
      p = FALSE;
    }

    if (p) {
    } else {
      idx[k] = (int8_T)(k + 2);
      idx[k + 1] = (int8_T)(k + 1);
    }
  }

  for (i = 0; i < 6; i++) {
    idx0[i] = 1;
  }

  i = 2;
  while (i < 6) {
    nb = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < 7; pEnd = qEnd + i) {
      b_p = j;
      q = pEnd;
      qEnd = j + nb;
      if (qEnd > 7) {
        qEnd = 7;
      }

      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        if (a[idx[b_p - 1] - 1] <= a[idx[q - 1] - 1]) {
          p = TRUE;
        } else {
          p = FALSE;
        }

        if (p) {
          idx0[k] = idx[b_p - 1];
          b_p++;
          if (b_p == pEnd) {
            while (q < qEnd) {
              k++;
              idx0[k] = idx[q - 1];
              q++;
            }
          }
        } else {
          idx0[k] = idx[q - 1];
          q++;
          if (q == qEnd) {
            while (b_p < pEnd) {
              k++;
              idx0[k] = idx[b_p - 1];
              b_p++;
            }
          }
        }

        k++;
      }

      for (k = 1; k <= kEnd; k++) {
        idx[(j + k) - 2] = idx0[k - 1];
      }

      j = qEnd;
    }

    i = nb;
  }

  for (pEnd = 0; pEnd < 2; pEnd++) {
    j = b->size[0] * b->size[1];
    b->size[pEnd] = 1 + 5 * pEnd;
    emxEnsureCapacity((emxArray__common *)b, j, (int32_T)sizeof(int32_T));
  }

  for (k = 0; k < 6; k++) {
    b->data[k] = a[idx[k] - 1];
  }

  nb = 0;
  k = 1;
  while (k <= 6) {
    i = b->data[k - 1];
    do {
      k++;
    } while (!((k > 6) || (!(b->data[k - 1] == i))));

    nb++;
    b->data[nb - 1] = i;
  }

  if (1 > nb) {
    nb = 0;
  }

  b_emxInit_int32_T(&r5, 1);
  pEnd = r5->size[0];
  r5->size[0] = nb;
  emxEnsureCapacity((emxArray__common *)r5, pEnd, (int32_T)sizeof(int32_T));
  i = nb - 1;
  for (pEnd = 0; pEnd <= i; pEnd++) {
    r5->data[pEnd] = 1 + pEnd;
  }

  emxInit_int32_T(&b_b, 2);
  iv16[0] = 1;
  iv16[1] = r5->size[0];
  pEnd = b_b->size[0] * b_b->size[1];
  b_b->size[0] = iv16[0];
  b_b->size[1] = iv16[1];
  emxEnsureCapacity((emxArray__common *)b_b, pEnd, (int32_T)sizeof(int32_T));
  i = iv16[1] - 1;
  for (pEnd = 0; pEnd <= i; pEnd++) {
    nb = iv16[0] - 1;
    for (j = 0; j <= nb; j++) {
      r6 = *r5;
      r6.size = (int32_T *)&iv16;
      r6.numDimensions = 1;
      b_b->data[j + b_b->size[0] * pEnd] = b->data[r6.data[j + r6.size[0] * pEnd]
        - 1];
    }
  }

  emxFree_int32_T(&r5);
  pEnd = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = b_b->size[1];
  emxEnsureCapacity((emxArray__common *)b, pEnd, (int32_T)sizeof(int32_T));
  i = b_b->size[1] - 1;
  for (pEnd = 0; pEnd <= i; pEnd++) {
    b->data[b->size[0] * pEnd] = b_b->data[b_b->size[0] * pEnd];
  }

  emxFree_int32_T(&b_b);
}

/*
 * function [ps_vol,tets] = call_tetgen(ps_local,tris_local)
 */

/*
 * function [Fps_new,Ftris]=call_triangle(Fps,Elems)
 */

/*
 * function [flag_v] = chk_pnt_inside_MB(cell_bnds,tri,flag_v )
 *  hx = grid_spacing(1);
 *  hy = grid_spacing(2);
 *  hz = grid_spacing(3);
 *  x_plane_pos = c_centers(1)+hx/2;
 *  x_plane_neg = c_centers(1)-hx/2;
 *  y_plane_pos = c_centers(2)+hy/2;
 *  y_plane_neg = c_centers(2)-hy/2;
 *  z_plane_pos = c_centers(3)+hz/2;
 *  z_plane_neg = c_centers(3)-hz/2;
 */
static void chk_pnt_inside_MB(const real_T cell_bnds[6], const real_T tri[9],
  int32_T flag_v[3])
{
  /* 'chk_pnt_inside_MB:12' xmin = cell_bnds(1); */
  /* 'chk_pnt_inside_MB:13' xmax = cell_bnds(2); */
  /* 'chk_pnt_inside_MB:14' ymin = cell_bnds(3); */
  /* 'chk_pnt_inside_MB:15' ymax = cell_bnds(4); */
  /* 'chk_pnt_inside_MB:16' zmin = cell_bnds(5); */
  /* 'chk_pnt_inside_MB:17' zmax = cell_bnds(6); */
  /* 'chk_pnt_inside_MB:19' v_1 = tri(1,:); */
  /* 'chk_pnt_inside_MB:19' v_2 = tri(2,:); */
  /* 'chk_pnt_inside_MB:19' v_3 = tri(3,:); */
  /* 'chk_pnt_inside_MB:21' if ((v_1(1)<= xmax) && (v_1(1)>= xmin)&& ... */
  /* 'chk_pnt_inside_MB:22'         (v_1(2)<= ymax) && (v_1(2)>= ymin)&& ... */
  /* 'chk_pnt_inside_MB:23'         (v_1(3)<= zmax) && (v_1(3)>= zmin)) */
  if ((tri[0] <= cell_bnds[1]) && (tri[0] >= cell_bnds[0]) && (tri[3] <=
       cell_bnds[3]) && (tri[3] >= cell_bnds[2]) && (tri[6] <= cell_bnds[5]) &&
      (tri[6] >= cell_bnds[4])) {
    /* 'chk_pnt_inside_MB:24' flag_v(1) = 7; */
    flag_v[0] = 7;
  }

  /* 'chk_pnt_inside_MB:27' if ((v_2(1)<= xmax) && (v_2(1)>= xmin)&& ... */
  /* 'chk_pnt_inside_MB:28'         (v_2(2)<= ymax) && (v_2(2)>= ymin)&& ... */
  /* 'chk_pnt_inside_MB:29'         (v_2(3)<= zmax) && (v_2(3)>= zmin)) */
  if ((tri[1] <= cell_bnds[1]) && (tri[1] >= cell_bnds[0]) && (tri[4] <=
       cell_bnds[3]) && (tri[4] >= cell_bnds[2]) && (tri[7] <= cell_bnds[5]) &&
      (tri[7] >= cell_bnds[4])) {
    /* 'chk_pnt_inside_MB:30' flag_v(2) = 7; */
    flag_v[1] = 7;
  }

  /* 'chk_pnt_inside_MB:33' if ((v_3(1)<= xmax) && (v_3(1)>= xmin)&& ... */
  /* 'chk_pnt_inside_MB:34'         (v_3(2)<= ymax) && (v_3(2)>= ymin)&& ... */
  /* 'chk_pnt_inside_MB:35'         (v_3(3)<= zmax) && (v_3(3)>= zmin)) */
  if ((tri[2] <= cell_bnds[1]) && (tri[2] >= cell_bnds[0]) && (tri[5] <=
       cell_bnds[3]) && (tri[5] >= cell_bnds[2]) && (tri[8] <= cell_bnds[5]) &&
      (tri[8] >= cell_bnds[4])) {
    /* 'chk_pnt_inside_MB:36' flag_v(3) = 7; */
    flag_v[2] = 7;
  }
}

/*
 * function [SameVert, SameOrient]=chk_tri_orientation(Tri1,Tri2)
 *  Given two triangles, this function checks 1. if they have the same set of
 *  vertices, and 2. if they have the same orientation.
 */
static void chk_tri_orientation(const int32_T Tri1[3], const int32_T Tri2[3],
  real_T *SameVert, real_T *SameOrient)
{
  int32_T b_Tri1[6];
  int32_T i9;
  emxArray_int32_T *T;

  /* 'chk_tri_orientation:4' T = [Tri1 Tri2]; */
  /* 'chk_tri_orientation:5' T = unique(T); */
  for (i9 = 0; i9 < 3; i9++) {
    b_Tri1[i9] = Tri1[i9];
  }

  for (i9 = 0; i9 < 3; i9++) {
    b_Tri1[i9 + 3] = Tri2[i9];
  }

  emxInit_int32_T(&T, 2);
  c_unique(b_Tri1, T);

  /* 'chk_tri_orientation:6' if (size(T,2) == 3) */
  if (T->size[1] == 3) {
    /* 'chk_tri_orientation:7' SameVert = 1; */
    *SameVert = 1.0;
  } else {
    /* 'chk_tri_orientation:8' else */
    /* 'chk_tri_orientation:9' SameVert = 0; */
    *SameVert = 0.0;
  }

  emxFree_int32_T(&T);

  /* 'chk_tri_orientation:12' comb1 = [Tri1(1) Tri1(2) Tri1(3)]; */
  /* 'chk_tri_orientation:13' comb2 = [Tri1(2) Tri1(3) Tri1(1)]; */
  /* 'chk_tri_orientation:14' comb3 = [Tri1(3) Tri1(1) Tri1(2)]; */
  /* 'chk_tri_orientation:17' if (((Tri2(1) == comb1(1)) && (Tri2(2) == comb1(2)) && (Tri2(3) == comb1(3))) || ... */
  /* 'chk_tri_orientation:18'     ((Tri2(1) == comb2(1)) && (Tri2(2) == comb2(2)) && (Tri2(3) == comb2(3))) || ... */
  /* 'chk_tri_orientation:19'     ((Tri2(1) == comb3(1)) && (Tri2(2) == comb3(2)) && (Tri2(3) == comb3(3)))) */
  if (((Tri2[0] == Tri1[0]) && (Tri2[1] == Tri1[1]) && (Tri2[2] == Tri1[2])) ||
      ((Tri2[0] == Tri1[1]) && (Tri2[1] == Tri1[2]) && (Tri2[2] == Tri1[0])) ||
      ((Tri2[0] == Tri1[2]) && (Tri2[1] == Tri1[0]) && (Tri2[2] == Tri1[1]))) {
    /* 'chk_tri_orientation:20' SameOrient = 1; */
    *SameOrient = 1.0;
  } else {
    /* 'chk_tri_orientation:21' else */
    /* 'chk_tri_orientation:22' SameOrient = 0; */
    *SameOrient = 0.0;
  }
}

/*
 * function [ps_local,ps_localID,ps_idx]=compute_intersections(cell_bnds,opphes,trivertID,triID,tri,ps_local,ps_localID,ps_idx)
 *  Face details
 * [~,corner] = cell_detail(c_centers,grid_spacing);
 */
static void compute_intersections(const real_T cell_bnds[6], const
  emxArray_int32_T *opphes, const int32_T trivertID[3], int32_T triID, const
  real_T tri[9], emxArray_real_T *ps_local, emxArray_int32_T *ps_localID,
  int32_T *ps_idx)
{
  uint32_T b_opphes[3];
  int32_T i12;
  uint32_T uv6[3];
  int32_T ngbtriID[3];
  int32_T c_opphes[3];
  int32_T ngbtriEdgeID[3];
  int32_T kk;
  int32_T count;
  int32_T i;
  boolean_T exitg4;
  int32_T b_ps_idx;
  int32_T b_trivertID[2];
  static const int8_T elems[6] = { 1, 2, 3, 2, 3, 1 };

  real_T b_tri[6];
  real_T corner[24];
  real_T E1[3];
  real_T E2[3];
  static const int8_T cell_elems[24] = { 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 2,
    3, 4, 1, 5, 6, 7, 8, 6, 7, 8, 5 };

  real_T D[3];
  real_T T[3];
  real_T b_D[3];
  real_T deno;
  real_T t;
  real_T u;
  boolean_T exitg3;
  boolean_T exitg2;
  boolean_T exitg1;

  /* 'compute_intersections:4' elems = [1 2; 2 3; 3 1]; */
  /* 'compute_intersections:5' ngbtriID = int32( bitshift( uint32(opphes(triID,:)),-2)); */
  for (i12 = 0; i12 < 3; i12++) {
    b_opphes[i12] = (uint32_T)opphes->data[(triID + opphes->size[0] * i12) - 1];
  }

  bitshift(b_opphes, -2.0, uv6);
  for (i12 = 0; i12 < 3; i12++) {
    ngbtriID[i12] = (int32_T)uv6[i12];
  }

  /* 'compute_intersections:6' ngbtriEdgeID = mod(opphes(triID,:),4)+1; */
  for (i12 = 0; i12 < 3; i12++) {
    c_opphes[i12] = opphes->data[(triID + opphes->size[0] * i12) - 1];
  }

  b_mod(c_opphes, 4.0, ngbtriEdgeID);
  for (i12 = 0; i12 < 3; i12++) {
    ngbtriEdgeID[i12]++;
  }

  /* Intersection of Triangle Edges with Cell Faces */
  /* 'compute_intersections:9' for kk=1:3 */
  for (kk = 0; kk < 3; kk++) {
    /* 'compute_intersections:10' count = 0; */
    count = 0;

    /* 'compute_intersections:11' if (ngbtriID(kk) ~=0) */
    if (ngbtriID[kk] != 0) {
      /* 'compute_intersections:12' for ii=1:(ps_idx-1) */
      i12 = *ps_idx - 1;
      i = 0;
      exitg4 = 0U;
      while ((exitg4 == 0U) && (i + 1 <= i12)) {
        /* 'compute_intersections:13' if ((ps_localID(ii,2)==ngbtriID(kk)) && (ps_localID(ii,3)==ngbtriEdgeID(kk))) */
        if ((ps_localID->data[i + ps_localID->size[0]] == ngbtriID[kk]) &&
            (ps_localID->data[i + (ps_localID->size[0] << 1)] == ngbtriEdgeID[kk]))
        {
          /* 'compute_intersections:14' count = 1; */
          count = 1;
          exitg4 = 1U;
        } else {
          i++;
        }
      }
    }

    /* 'compute_intersections:19' if ((count == 0)||(ngbtriID(kk) == 0)) */
    if ((count == 0) || (ngbtriID[kk] == 0)) {
      /* 'compute_intersections:20' [ps_idx,ps_local,ps_localID] =intersection_triedge_cellface_2(cell_bnds,triID,kk,trivertID(elems(kk,:)),tri(elems(kk,:),:),ps_idx,ps_local,ps_localID); */
      b_ps_idx = *ps_idx;
      for (i12 = 0; i12 < 2; i12++) {
        b_trivertID[i12] = trivertID[elems[kk + 3 * i12] - 1];
      }

      for (i12 = 0; i12 < 3; i12++) {
        for (i = 0; i < 2; i++) {
          b_tri[i + (i12 << 1)] = tri[(elems[kk + 3 * i] + 3 * i12) - 1];
        }
      }

      intersection_triedge_cellface_2(cell_bnds, triID, 1.0 + (real_T)kk,
        b_trivertID, b_tri, &b_ps_idx, ps_local, ps_localID);
      *ps_idx = b_ps_idx;
    }
  }

  /*  Intersection of cell edges with triangle */
  /* 'compute_intersections:25' [ps_idx,ps_local,ps_localID]=intersection_celledg_triface_2(cell_bnds,tri,triID,ps_idx,ps_local,ps_localID); */
  b_ps_idx = *ps_idx - 1;

  /* 'intersection_celledg_triface_2:2' [~,corner] = gridcell_data(cell_bnds); */
  /*  Bounds of the cell */
  /* 'gridcell_data:4' xmin = cell_bnds(1); */
  /* 'gridcell_data:5' xmax = cell_bnds(2); */
  /* 'gridcell_data:6' ymin = cell_bnds(3); */
  /* 'gridcell_data:7' ymax = cell_bnds(4); */
  /* 'gridcell_data:8' zmin = cell_bnds(5); */
  /* 'gridcell_data:9' zmax = cell_bnds(6); */
  /*  Face-Corner Map */
  /* 'gridcell_data:12' faces = int32([1 4 3 2;1 2 6 5;2 3 7 6;3 4 8 7;1 5 8 4;5 6 7 8]); */
  /*  Corners of the cell */
  /* 'gridcell_data:15' corners = zeros(8,3); */
  for (i12 = 0; i12 < 24; i12++) {
    corner[i12] = 0.0;
  }

  /* 'gridcell_data:16' corners(1,:) = [xmax ymin zmin]; */
  corner[0] = cell_bnds[1];
  corner[8] = cell_bnds[2];
  corner[16] = cell_bnds[4];

  /* 'gridcell_data:17' corners(2,:) = [xmax ymax zmin]; */
  corner[1] = cell_bnds[1];
  corner[9] = cell_bnds[3];
  corner[17] = cell_bnds[4];

  /* 'gridcell_data:18' corners(3,:) = [xmin ymax zmin]; */
  corner[2] = cell_bnds[0];
  corner[10] = cell_bnds[3];
  corner[18] = cell_bnds[4];

  /* 'gridcell_data:19' corners(4,:) = [xmin ymin zmin]; */
  corner[3] = cell_bnds[0];
  corner[11] = cell_bnds[2];
  corner[19] = cell_bnds[4];

  /* 'gridcell_data:20' corners(5,:) = [xmax ymin zmax]; */
  corner[4] = cell_bnds[1];
  corner[12] = cell_bnds[2];
  corner[20] = cell_bnds[5];

  /* 'gridcell_data:21' corners(6,:) = [xmax ymax zmax]; */
  corner[5] = cell_bnds[1];
  corner[13] = cell_bnds[3];
  corner[21] = cell_bnds[5];

  /* 'gridcell_data:22' corners(7,:) = [xmin ymax zmax]; */
  corner[6] = cell_bnds[0];
  corner[14] = cell_bnds[3];
  corner[22] = cell_bnds[5];

  /* 'gridcell_data:23' corners(8,:) = [xmin ymin zmax]; */
  corner[7] = cell_bnds[0];
  corner[15] = cell_bnds[2];
  corner[23] = cell_bnds[5];

  /*  Normals to the cell faces */
  /* 'gridcell_data:26' nrm = zeros(6,3); */
  /* 'gridcell_data:27' nrm(1,:) = [0 0 -1]; */
  /* 'gridcell_data:28' nrm(2,:) = [1 0 0]; */
  /* 'gridcell_data:29' nrm(3,:) = [0 1 0]; */
  /* 'gridcell_data:30' nrm(4,:) = [-1 0 0]; */
  /* 'gridcell_data:31' nrm(5,:) = [0 -1 0]; */
  /* 'gridcell_data:32' nrm(6,:) = [0 0 1]; */
  /*  Equation of the cell faces */
  /* 'gridcell_data:35' face_equ = zeros(6,1); */
  /* 'gridcell_data:36' face_equ(1) = zmin; */
  /* 'gridcell_data:37' face_equ(2) = xmax; */
  /* 'gridcell_data:38' face_equ(3) = ymax; */
  /* 'gridcell_data:39' face_equ(4) = xmin; */
  /* 'gridcell_data:40' face_equ(5) = ymin; */
  /* 'gridcell_data:41' face_equ(6) = zmax; */
  /*  ID of the face to distinguish if it is the x/y/z-plane */
  /* 'gridcell_data:44' face_equID = zeros(6,1,'int32'); */
  /* 'gridcell_data:45' face_equID(1) = 3; */
  /* 'gridcell_data:46' face_equID(2) = 1; */
  /* 'gridcell_data:47' face_equID(3) = 2; */
  /* 'gridcell_data:48' face_equID(4) = 1; */
  /* 'gridcell_data:49' face_equID(5) = 2; */
  /* 'gridcell_data:50' face_equID(6) = 3; */
  /*  Face bounds and map with coordinate for each cell face. */
  /*  For a typical face, */
  /*  fbmap(1),fbmap(2) : The coordinates that vary for the face, their bounds are */
  /*  given by facebnd(1:2),facebnd(3:4), respectively. fbmap(3) is the */
  /*  coordinate that remains constant and this constant value is facebnd(5). */
  /* 'gridcell_data:58' facebnd = [xmin xmax ymin ymax zmin; */
  /* 'gridcell_data:59'     ymin ymax zmin zmax xmax; */
  /* 'gridcell_data:60'     xmin xmax zmin zmax ymax; */
  /* 'gridcell_data:61'     ymin ymax zmin zmax xmin; */
  /* 'gridcell_data:62'     xmin xmax zmin zmax ymin; */
  /* 'gridcell_data:63'     xmin xmax ymin ymax zmax]; */
  /* 'gridcell_data:65' fbmap =[1 2 3; 2 3 1; 1 3 2; 2 3 1;1 3 2;1 2 3]; */
  /*  Mapping of the barycentric coordinates with the corners and cell edges. */
  /*  Corners: fcormap(1:4) corresponds to (xi,eta) = (0,0),(1,0),(1,1),(0,1) */
  /*  Edges : fedgmap(1:4) corresponds to(xi,eta) = ((0,1),0), (1,(0,1)),((0,1),1), (0,(0,1)) */
  /* 'gridcell_data:70' fcormap =[4 1 2 3; 1 2 6 5; 3 2 6 7; 4 3 7 8; 4 1 5 8; 8 5 6 7]; */
  /* 'gridcell_data:71' fedgmap = [4 1 2 3; 1 6 9 5; 2 6 10 7; 3 7 11 8; 4 5 12 8; 12 9 10 11]; */
  /* 'intersection_celledg_triface_2:3' cell_elems =[1 2; 2 3; 3 4; 4 1; 1 5; 2 6; 3 7; 4 8; 5 6; 6 7; 7 8; 8 5]; */
  /* 'intersection_celledg_triface_2:4' E1 = tri(2,:)-tri(1,:); */
  for (i12 = 0; i12 < 3; i12++) {
    E1[i12] = tri[1 + 3 * i12] - tri[3 * i12];

    /* 'intersection_celledg_triface_2:5' E2 = tri(3,:)-tri(1,:); */
    E2[i12] = tri[2 + 3 * i12] - tri[3 * i12];
  }

  /* 'intersection_celledg_triface_2:6' tol = 1e-6; */
  /* 'intersection_celledg_triface_2:7' for kk=1:12 */
  for (kk = 0; kk < 12; kk++) {
    /* 'intersection_celledg_triface_2:8' D = corner(cell_elems(kk,2),:) - corner(cell_elems(kk,1),:); */
    /* 'intersection_celledg_triface_2:9' T = corner(cell_elems(kk,1),:) - tri(1,:); */
    for (i12 = 0; i12 < 3; i12++) {
      D[i12] = corner[(cell_elems[12 + kk] + (i12 << 3)) - 1] - corner
        [(cell_elems[kk] + (i12 << 3)) - 1];
      T[i12] = corner[(cell_elems[kk] + (i12 << 3)) - 1] - tri[3 * i12];
    }

    /* 'intersection_celledg_triface_2:10' deno = dot(E1,cross_col(D,E2)); */
    /* CROSS_COL Efficient routine for computing cross product of two  */
    /* 3-dimensional column vectors. */
    /*  CROSS_COL(A,B) Efficiently computes the cross product between */
    /*  3-dimensional column vector A, and 3-dimensional column vector B. */
    /* 'cross_col:7' c = [a(2)*b(3)-a(3)*b(2); a(3)*b(1)-a(1)*b(3); a(1)*b(2)-a(2)*b(1)]; */
    b_D[0] = D[1] * E2[2] - D[2] * E2[1];
    b_D[1] = D[2] * E2[0] - D[0] * E2[2];
    b_D[2] = D[0] * E2[1] - D[1] * E2[0];
    deno = eml_xdot(3, E1, 1, 1, b_D, 1, 1);

    /* 'intersection_celledg_triface_2:12' t = dot(E2,cross_col(T,E1))/deno; */
    /* CROSS_COL Efficient routine for computing cross product of two  */
    /* 3-dimensional column vectors. */
    /*  CROSS_COL(A,B) Efficiently computes the cross product between */
    /*  3-dimensional column vector A, and 3-dimensional column vector B. */
    /* 'cross_col:7' c = [a(2)*b(3)-a(3)*b(2); a(3)*b(1)-a(1)*b(3); a(1)*b(2)-a(2)*b(1)]; */
    b_D[0] = T[1] * E1[2] - T[2] * E1[1];
    b_D[1] = T[2] * E1[0] - T[0] * E1[2];
    b_D[2] = T[0] * E1[1] - T[1] * E1[0];
    t = eml_xdot(3, E2, 1, 1, b_D, 1, 1) / deno;

    /* 'intersection_celledg_triface_2:13' u = dot(T,cross_col(D,E2))/deno; */
    /* CROSS_COL Efficient routine for computing cross product of two  */
    /* 3-dimensional column vectors. */
    /*  CROSS_COL(A,B) Efficiently computes the cross product between */
    /*  3-dimensional column vector A, and 3-dimensional column vector B. */
    /* 'cross_col:7' c = [a(2)*b(3)-a(3)*b(2); a(3)*b(1)-a(1)*b(3); a(1)*b(2)-a(2)*b(1)]; */
    b_D[0] = D[1] * E2[2] - D[2] * E2[1];
    b_D[1] = D[2] * E2[0] - D[0] * E2[2];
    b_D[2] = D[0] * E2[1] - D[1] * E2[0];
    u = eml_xdot(3, T, 1, 1, b_D, 1, 1) / deno;

    /* 'intersection_celledg_triface_2:14' v = dot(D,cross_col(T,E1))/deno; */
    /* CROSS_COL Efficient routine for computing cross product of two  */
    /* 3-dimensional column vectors. */
    /*  CROSS_COL(A,B) Efficiently computes the cross product between */
    /*  3-dimensional column vector A, and 3-dimensional column vector B. */
    /* 'cross_col:7' c = [a(2)*b(3)-a(3)*b(2); a(3)*b(1)-a(1)*b(3); a(1)*b(2)-a(2)*b(1)]; */
    b_D[0] = T[1] * E1[2] - T[2] * E1[1];
    b_D[1] = T[2] * E1[0] - T[0] * E1[2];
    b_D[2] = T[0] * E1[1] - T[1] * E1[0];
    deno = eml_xdot(3, D, 1, 1, b_D, 1, 1) / deno;

    /* 'intersection_celledg_triface_2:16' if(abs(u)<=tol) */
    if (fabs(u) <= 1.0E-6) {
      /* 'intersection_celledg_triface_2:17' u = 0; */
      u = 0.0;
    } else {
      if (fabs(u - 1.0) <= 1.0E-6) {
        /* 'intersection_celledg_triface_2:18' elseif (abs(u-1)<=tol) */
        /* 'intersection_celledg_triface_2:19' u = 1; */
        u = 1.0;
      }
    }

    /* 'intersection_celledg_triface_2:21' if(abs(v)<=tol) */
    if (fabs(deno) <= 1.0E-6) {
      /* 'intersection_celledg_triface_2:22' v = 0; */
      deno = 0.0;
    } else {
      if (fabs(deno - 1.0) <= 1.0E-6) {
        /* 'intersection_celledg_triface_2:23' elseif (abs(v-1)<=tol) */
        /* 'intersection_celledg_triface_2:24' v = 1; */
        deno = 1.0;
      }
    }

    /* 'intersection_celledg_triface_2:26' if (abs(u+v-1)<=tol) */
    if (fabs((u + deno) - 1.0) <= 1.0E-6) {
      /* 'intersection_celledg_triface_2:27' u = 1-v; */
      u = 1.0 - deno;
    }

    /* 'intersection_celledg_triface_2:29' if(abs(t)<=tol) */
    if (fabs(t) <= 1.0E-6) {
      /* 'intersection_celledg_triface_2:30' t = 0; */
      t = 0.0;
    } else {
      if (fabs(t - 1.0) <= 1.0E-6) {
        /* 'intersection_celledg_triface_2:31' elseif (abs(t-1)<=tol) */
        /* 'intersection_celledg_triface_2:32' t = 1; */
        t = 1.0;
      }
    }

    /* 'intersection_celledg_triface_2:36' if (t >=0) && (t <=1) &&(u>0) && (v>0) && (u+v<1) */
    if ((t >= 0.0) && (t <= 1.0) && (u > 0.0) && (deno > 0.0) && (u + deno < 1.0))
    {
      /* 'intersection_celledg_triface_2:37' if (t ==0) */
      if (t == 0.0) {
        /* 'intersection_celledg_triface_2:38' p = corner(cell_elems(kk,1),:); */
        /* 'intersection_celledg_triface_2:39' count = 0; */
        count = 0;

        /* 'intersection_celledg_triface_2:40' for i=1:ps_idx-1 */
        i = 1;
        exitg3 = 0U;
        while ((exitg3 == 0U) && (i <= b_ps_idx)) {
          /* 'intersection_celledg_triface_2:41' if (ps_localID(i,6)==cell_elems(kk,1)) */
          if (ps_localID->data[(i + ps_localID->size[0] * 5) - 1] ==
              cell_elems[kk]) {
            /* 'intersection_celledg_triface_2:42' count = 1; */
            count = 1;
            exitg3 = 1U;
          } else {
            i++;
          }
        }

        /* 'intersection_celledg_triface_2:45' if (count ==0) */
        if (count == 0) {
          /* 'intersection_celledg_triface_2:46' ps_local(ps_idx,:) = p; */
          for (i12 = 0; i12 < 3; i12++) {
            ps_local->data[b_ps_idx + ps_local->size[0] * i12] = corner
              [(cell_elems[kk] + (i12 << 3)) - 1];
          }

          /* 'intersection_celledg_triface_2:47' ps_localID(ps_idx,:) = [0,triID,0,0,0,cell_elems(kk,1)]; */
          ps_localID->data[b_ps_idx] = 0;
          ps_localID->data[b_ps_idx + ps_localID->size[0]] = triID;
          ps_localID->data[b_ps_idx + (ps_localID->size[0] << 1)] = 0;
          ps_localID->data[b_ps_idx + ps_localID->size[0] * 3] = 0;
          ps_localID->data[b_ps_idx + (ps_localID->size[0] << 2)] = 0;
          ps_localID->data[b_ps_idx + ps_localID->size[0] * 5] = (int32_T)
            cell_elems[kk];

          /* 'intersection_celledg_triface_2:48' ps_idx = ps_idx + 1; */
          b_ps_idx++;
        }
      } else if (t == 1.0) {
        /* 'intersection_celledg_triface_2:50' elseif (t==1) */
        /* 'intersection_celledg_triface_2:51' p = corner(cell_elems(kk,2),:); */
        /* 'intersection_celledg_triface_2:52' count = 0; */
        count = 0;

        /* 'intersection_celledg_triface_2:53' for i=1:ps_idx-1 */
        i = 1;
        exitg2 = 0U;
        while ((exitg2 == 0U) && (i <= b_ps_idx)) {
          /* 'intersection_celledg_triface_2:54' if (ps_localID(i,6)==cell_elems(kk,2)) */
          if (ps_localID->data[(i + ps_localID->size[0] * 5) - 1] == cell_elems
              [12 + kk]) {
            /* 'intersection_celledg_triface_2:55' count = 1; */
            count = 1;
            exitg2 = 1U;
          } else {
            i++;
          }
        }

        /* 'intersection_celledg_triface_2:58' if (count ==0) */
        if (count == 0) {
          /* 'intersection_celledg_triface_2:59' ps_local(ps_idx,:) = p; */
          for (i12 = 0; i12 < 3; i12++) {
            ps_local->data[b_ps_idx + ps_local->size[0] * i12] = corner
              [(cell_elems[12 + kk] + (i12 << 3)) - 1];
          }

          /* 'intersection_celledg_triface_2:60' ps_localID(ps_idx,:) = [0,triID,0,0,0,cell_elems(kk,2)]; */
          ps_localID->data[b_ps_idx] = 0;
          ps_localID->data[b_ps_idx + ps_localID->size[0]] = triID;
          ps_localID->data[b_ps_idx + (ps_localID->size[0] << 1)] = 0;
          ps_localID->data[b_ps_idx + ps_localID->size[0] * 3] = 0;
          ps_localID->data[b_ps_idx + (ps_localID->size[0] << 2)] = 0;
          ps_localID->data[b_ps_idx + ps_localID->size[0] * 5] = (int32_T)
            cell_elems[12 + kk];

          /* 'intersection_celledg_triface_2:61' ps_idx = ps_idx + 1; */
          b_ps_idx++;
        }
      } else {
        /* 'intersection_celledg_triface_2:63' else */
        /* 'intersection_celledg_triface_2:64' count = 0; */
        count = 0;

        /* 'intersection_celledg_triface_2:65' for i=1:ps_idx-1 */
        i = 0;
        exitg1 = 0U;
        while ((exitg1 == 0U) && (i + 1 <= b_ps_idx)) {
          /* 'intersection_celledg_triface_2:66' if (ps_localID(i,2) == triID && ps_localID(i,5)==kk) */
          if ((ps_localID->data[i + ps_localID->size[0]] == triID) &&
              (ps_localID->data[i + (ps_localID->size[0] << 2)] == 1 + kk)) {
            /* 'intersection_celledg_triface_2:67' count = 1; */
            count = 1;
            exitg1 = 1U;
          } else {
            i++;
          }
        }

        /* 'intersection_celledg_triface_2:70' if (count ==0) */
        if (count == 0) {
          /* 'intersection_celledg_triface_2:71' p = corner(cell_elems(kk,1),:)+ t*(corner(cell_elems(kk,2),:)-corner(cell_elems(kk,1),:)); */
          for (i12 = 0; i12 < 3; i12++) {
            D[i12] = corner[(cell_elems[kk] + (i12 << 3)) - 1] + t * (corner
              [(cell_elems[12 + kk] + (i12 << 3)) - 1] - corner[(cell_elems[kk]
              + (i12 << 3)) - 1]);
          }

          /* 'intersection_celledg_triface_2:72' ps_local(ps_idx,:) = p; */
          for (i12 = 0; i12 < 3; i12++) {
            ps_local->data[b_ps_idx + ps_local->size[0] * i12] = D[i12];
          }

          /* 'intersection_celledg_triface_2:73' ps_localID(ps_idx,:) = [0,triID,0,0,kk,0]; */
          ps_localID->data[b_ps_idx] = 0;
          ps_localID->data[b_ps_idx + ps_localID->size[0]] = triID;
          ps_localID->data[b_ps_idx + (ps_localID->size[0] << 1)] = 0;
          ps_localID->data[b_ps_idx + ps_localID->size[0] * 3] = 0;
          ps_localID->data[b_ps_idx + (ps_localID->size[0] << 2)] = 1 + kk;
          ps_localID->data[b_ps_idx + ps_localID->size[0] * 5] = 0;

          /* 'intersection_celledg_triface_2:74' ps_idx = ps_idx + 1; */
          b_ps_idx++;
        }
      }
    }
  }

  *ps_idx = b_ps_idx + 1;
}

/*
 *
 */
static void d_unique(const emxArray_real_T *a, emxArray_real_T *b)
{
  int32_T na;
  int32_T n;
  uint32_T uv4[2];
  int32_T pEnd;
  emxArray_int32_T *idx;
  int32_T np1;
  int32_T k;
  boolean_T p;
  emxArray_int32_T *idx0;
  int32_T j;
  int32_T i;
  int32_T i2;
  int32_T b_p;
  int32_T q;
  int32_T qEnd;
  int32_T kEnd;
  real_T x;
  int32_T exitg1;
  real_T absxk;
  int32_T exponent;
  emxArray_int32_T *r7;
  emxArray_real_T *b_b;
  int32_T iv17[2];
  emxArray_int32_T r8;
  na = a->size[1];
  n = a->size[1];
  for (pEnd = 0; pEnd < 2; pEnd++) {
    uv4[pEnd] = (uint32_T)a->size[pEnd];
  }

  emxInit_int32_T(&idx, 2);
  pEnd = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  idx->size[1] = (int32_T)uv4[1];
  emxEnsureCapacity((emxArray__common *)idx, pEnd, (int32_T)sizeof(int32_T));
  np1 = n + 1;
  if (a->size[1] == 0) {
    for (k = 1; k <= n; k++) {
      idx->data[k - 1] = k;
    }
  } else {
    for (k = 1; k <= n; k++) {
      idx->data[k - 1] = k;
    }

    for (k = 1; k <= n - 1; k += 2) {
      if (a->data[k - 1] <= a->data[k]) {
        p = TRUE;
      } else {
        p = FALSE;
      }

      if (p) {
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    b_emxInit_int32_T(&idx0, 1);
    pEnd = idx0->size[0];
    idx0->size[0] = n;
    emxEnsureCapacity((emxArray__common *)idx0, pEnd, (int32_T)sizeof(int32_T));
    j = n - 1;
    for (pEnd = 0; pEnd <= j; pEnd++) {
      idx0->data[pEnd] = 1;
    }

    i = 2;
    while (i < n) {
      i2 = i << 1;
      j = 1;
      for (pEnd = 1 + i; pEnd < np1; pEnd = qEnd + i) {
        b_p = j;
        q = pEnd;
        qEnd = j + i2;
        if (qEnd > np1) {
          qEnd = np1;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          if (a->data[idx->data[b_p - 1] - 1] <= a->data[idx->data[q - 1] - 1])
          {
            p = TRUE;
          } else {
            p = FALSE;
          }

          if (p) {
            idx0->data[k] = idx->data[b_p - 1];
            b_p++;
            if (b_p == pEnd) {
              while (q < qEnd) {
                k++;
                idx0->data[k] = idx->data[q - 1];
                q++;
              }
            }
          } else {
            idx0->data[k] = idx->data[q - 1];
            q++;
            if (q == qEnd) {
              while (b_p < pEnd) {
                k++;
                idx0->data[k] = idx->data[b_p - 1];
                b_p++;
              }
            }
          }

          k++;
        }

        for (k = -1; k + 2 <= kEnd; k++) {
          idx->data[j + k] = idx0->data[k + 1];
        }

        j = qEnd;
      }

      i = i2;
    }

    emxFree_int32_T(&idx0);
  }

  for (pEnd = 0; pEnd < 2; pEnd++) {
    uv4[pEnd] = (uint32_T)a->size[pEnd];
  }

  pEnd = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = (int32_T)uv4[1];
  emxEnsureCapacity((emxArray__common *)b, pEnd, (int32_T)sizeof(real_T));
  for (k = 0; k + 1 <= na; k++) {
    b->data[k] = a->data[idx->data[k] - 1];
  }

  emxFree_int32_T(&idx);
  i2 = 0;
  k = 1;
  while (k <= na) {
    x = b->data[k - 1];
    do {
      exitg1 = 0U;
      k++;
      if (k > na) {
        exitg1 = 1U;
      } else {
        absxk = fabs(x / 2.0);
        if (absxk <= 2.2250738585072014E-308) {
          absxk = 4.94065645841247E-324;
        } else {
          frexp(absxk, &exponent);
          absxk = ldexp(1.0, exponent - 53);
        }

        if (fabs(x - b->data[k - 1]) < absxk) {
          p = TRUE;
        } else {
          p = FALSE;
        }

        if (!p) {
          exitg1 = 1U;
        }
      }
    } while (exitg1 == 0U);

    i2++;
    b->data[i2 - 1] = x;
  }

  if (1 > i2) {
    i2 = 0;
  }

  b_emxInit_int32_T(&r7, 1);
  pEnd = r7->size[0];
  r7->size[0] = i2;
  emxEnsureCapacity((emxArray__common *)r7, pEnd, (int32_T)sizeof(int32_T));
  j = i2 - 1;
  for (pEnd = 0; pEnd <= j; pEnd++) {
    r7->data[pEnd] = 1 + pEnd;
  }

  emxInit_real_T(&b_b, 2);
  iv17[0] = 1;
  iv17[1] = r7->size[0];
  pEnd = b_b->size[0] * b_b->size[1];
  b_b->size[0] = iv17[0];
  b_b->size[1] = iv17[1];
  emxEnsureCapacity((emxArray__common *)b_b, pEnd, (int32_T)sizeof(real_T));
  j = iv17[1] - 1;
  for (pEnd = 0; pEnd <= j; pEnd++) {
    i2 = iv17[0] - 1;
    for (i = 0; i <= i2; i++) {
      r8 = *r7;
      r8.size = (int32_T *)&iv17;
      r8.numDimensions = 1;
      b_b->data[i + b_b->size[0] * pEnd] = b->data[r8.data[i + r8.size[0] * pEnd]
        - 1];
    }
  }

  emxFree_int32_T(&r7);
  pEnd = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = b_b->size[1];
  emxEnsureCapacity((emxArray__common *)b, pEnd, (int32_T)sizeof(real_T));
  j = b_b->size[1] - 1;
  for (pEnd = 0; pEnd <= j; pEnd++) {
    b->data[b->size[0] * pEnd] = b_b->data[b_b->size[0] * pEnd];
  }

  emxFree_real_T(&b_b);
}

/*
 * function opphes = determine_opposite_halfedge_tri(nv, tris, opphes)
 * DETERMINE_OPPOSITE_HALFEDGE_TRI Determine opposite half-edges for triangle
 * mesh.
 *  DETERMINE_OPPOSITE_HALFEDGE_TRI(NV,TRIS,OPPHES) Determines
 *  opposite half-edges for triangle mesh. The following explains the input
 *  and output arguments.
 *
 *  OPPHES = DETERMINE_OPPOSITE_HALFEDGE_TRI(NV,TRIS)
 *  OPPHES = DETERMINE_OPPOSITE_HALFEDGE_TRI(NV,TRIS,OPPHES)
 *  Computes mapping from each half-edge to its opposite half-edge for
 *  triangle mesh.
 *
 *  Convention: Each half-edge is indicated by <face_id,local_edge_id>.
 *  We assign 2 bits to local_edge_id (starts from 0).
 *
 *  See also DETERMINE_OPPOSITE_HALFEDGE
 */
static void determine_opposite_halfedge_tri(int32_T nv, const emxArray_int32_T
  *tris, emxArray_int32_T *opphes)
{
  emxArray_int32_T *is_index;
  int32_T ntris;
  int32_T i3;
  int32_T ii;
  boolean_T exitg4;
  int32_T b_is_index[3];
  emxArray_int32_T *v2nv;
  emxArray_int32_T *v2he;
  int32_T ne;
  static const int8_T iv7[3] = { 1, 2, 0 };

  boolean_T exitg3;
  int32_T exitg2;
  boolean_T guard1 = FALSE;
  int32_T found;
  static const int8_T iv8[3] = { 2, 3, 1 };

  int32_T b_index;
  int32_T exitg1;
  boolean_T guard2 = FALSE;
  uint32_T a;
  b_emxInit_int32_T(&is_index, 1);

  /* 'determine_opposite_halfedge_tri:18' coder.inline('never'); */
  /* 'determine_opposite_halfedge_tri:20' nepE = int32(3); */
  /*  Number of edges per element */
  /* 'determine_opposite_halfedge_tri:21' next = int32([2,3,1]); */
  /* 'determine_opposite_halfedge_tri:22' inds = int32(1:3); */
  /* 'determine_opposite_halfedge_tri:24' ntris = int32(size(tris,1)); */
  ntris = tris->size[0];

  /* % First, build is_index to store starting position for each vertex. */
  /* 'determine_opposite_halfedge_tri:26' is_index = zeros(nv+1,1,'int32'); */
  i3 = is_index->size[0];
  is_index->size[0] = nv + 1;
  emxEnsureCapacity((emxArray__common *)is_index, i3, (int32_T)sizeof(int32_T));
  for (i3 = 0; i3 <= nv; i3++) {
    is_index->data[i3] = 0;
  }

  /* 'determine_opposite_halfedge_tri:27' for ii=1:ntris */
  ii = 0;
  exitg4 = 0U;
  while ((exitg4 == 0U) && (ii + 1 <= ntris)) {
    /* 'determine_opposite_halfedge_tri:28' if tris(ii,1)==0 */
    if (tris->data[ii] == 0) {
      /* 'determine_opposite_halfedge_tri:28' ntris=ii-1; */
      ntris = ii;
      exitg4 = 1U;
    } else {
      /* 'determine_opposite_halfedge_tri:29' is_index(tris(ii,inds)+1) = is_index(tris(ii,inds)+1) + 1; */
      for (i3 = 0; i3 < 3; i3++) {
        b_is_index[i3] = is_index->data[tris->data[ii + tris->size[0] * i3]] + 1;
      }

      for (i3 = 0; i3 < 3; i3++) {
        is_index->data[tris->data[ii + tris->size[0] * i3]] = b_is_index[i3];
      }

      ii++;
    }
  }

  /* 'determine_opposite_halfedge_tri:31' is_index(1) = 1; */
  is_index->data[0] = 1;

  /* 'determine_opposite_halfedge_tri:32' for ii=1:nv */
  for (ii = 1; ii <= nv; ii++) {
    /* 'determine_opposite_halfedge_tri:33' is_index(ii+1) = is_index(ii) + is_index(ii+1); */
    is_index->data[ii] += is_index->data[ii - 1];
  }

  b_emxInit_int32_T(&v2nv, 1);
  b_emxInit_int32_T(&v2he, 1);

  /* 'determine_opposite_halfedge_tri:36' ne = ntris*nepE; */
  ne = ntris * 3;

  /* 'determine_opposite_halfedge_tri:37' v2nv = coder.nullcopy(zeros( ne,1, 'int32')); */
  i3 = v2nv->size[0];
  v2nv->size[0] = ne;
  emxEnsureCapacity((emxArray__common *)v2nv, i3, (int32_T)sizeof(int32_T));

  /*  Vertex to next vertex in each halfedge. */
  /* 'determine_opposite_halfedge_tri:38' v2he = coder.nullcopy(zeros( ne,1, 'int32')); */
  i3 = v2he->size[0];
  v2he->size[0] = ne;
  emxEnsureCapacity((emxArray__common *)v2he, i3, (int32_T)sizeof(int32_T));

  /*  Vertex to half-edge. */
  /* 'determine_opposite_halfedge_tri:39' for ii=1:ntris */
  for (ii = 0; ii + 1 <= ntris; ii++) {
    /* 'determine_opposite_halfedge_tri:40' v2nv(is_index( tris(ii,inds))) = tris(ii,next); */
    for (i3 = 0; i3 < 3; i3++) {
      v2nv->data[is_index->data[tris->data[ii + tris->size[0] * i3] - 1] - 1] =
        tris->data[ii + tris->size[0] * iv7[i3]];
    }

    /* 'determine_opposite_halfedge_tri:41' v2he(is_index( tris(ii,inds))) = 4*ii-1+inds; */
    ne = (ii + 1) << 2;
    for (i3 = 0; i3 < 3; i3++) {
      v2he->data[is_index->data[tris->data[ii + tris->size[0] * i3] - 1] - 1] =
        i3 + ne;
    }

    /* 'determine_opposite_halfedge_tri:42' is_index(tris(ii,inds)) = is_index(tris(ii,inds)) + 1; */
    for (i3 = 0; i3 < 3; i3++) {
      b_is_index[i3] = is_index->data[tris->data[ii + tris->size[0] * i3] - 1] +
        1;
    }

    for (i3 = 0; i3 < 3; i3++) {
      is_index->data[tris->data[ii + tris->size[0] * i3] - 1] = b_is_index[i3];
    }
  }

  /* 'determine_opposite_halfedge_tri:44' for ii=nv-1:-1:1 */
  for (ii = nv - 1; ii > 0; ii--) {
    /* 'determine_opposite_halfedge_tri:44' is_index(ii+1) = is_index(ii); */
    is_index->data[ii] = is_index->data[ii - 1];
  }

  /* 'determine_opposite_halfedge_tri:45' is_index(1)=1; */
  is_index->data[0] = 1;

  /* % Set opphes */
  /* 'determine_opposite_halfedge_tri:47' if nargin<3 || isempty(opphes) */
  /* 'determine_opposite_halfedge_tri:48' opphes = zeros(size(tris,1), nepE, 'int32'); */
  ne = tris->size[0];
  i3 = opphes->size[0] * opphes->size[1];
  opphes->size[0] = ne;
  opphes->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)opphes, i3, (int32_T)sizeof(int32_T));
  ne = tris->size[0] * 3 - 1;
  for (i3 = 0; i3 <= ne; i3++) {
    opphes->data[i3] = 0;
  }

  /* 'determine_opposite_halfedge_tri:54' for ii=1:ntris */
  ii = 0;
  exitg3 = 0U;
  while ((exitg3 == 0U) && (ii + 1 <= ntris)) {
    /* 'determine_opposite_halfedge_tri:55' for jj=int32(1):3 */
    ne = 0;
    do {
      exitg2 = 0U;
      if (ne + 1 < 4) {
        /* 'determine_opposite_halfedge_tri:56' if opphes(ii,jj) */
        guard1 = FALSE;
        if (opphes->data[ii + opphes->size[0] * ne] != 0) {
          guard1 = TRUE;
        } else {
          /* 'determine_opposite_halfedge_tri:57' v = tris(ii,jj); */
          /* 'determine_opposite_halfedge_tri:57' vn = tris(ii,next(jj)); */
          /*  LOCATE: Locate index col in v2nv(first:last) */
          /* 'determine_opposite_halfedge_tri:60' found = int32(0); */
          found = 0;

          /* 'determine_opposite_halfedge_tri:61' for index = is_index(vn):is_index(vn+1)-1 */
          i3 = is_index->data[tris->data[ii + tris->size[0] * (iv8[ne] - 1)]] -
            1;
          for (b_index = is_index->data[tris->data[ii + tris->size[0] * (iv8[ne]
                - 1)] - 1] - 1; b_index + 1 <= i3; b_index++) {
            /* 'determine_opposite_halfedge_tri:62' if v2nv(index)==v */
            if (v2nv->data[b_index] == tris->data[ii + tris->size[0] * ne]) {
              /* 'determine_opposite_halfedge_tri:63' opp = v2he(index); */
              /* 'determine_opposite_halfedge_tri:64' opphes(ii,jj) = opp; */
              opphes->data[ii + opphes->size[0] * ne] = v2he->data[b_index];

              /* opphes(heid2fid(opp),heid2leid(opp)) = ii*4+jj-1; */
              /* 'determine_opposite_halfedge_tri:66' opphes(bitshift(uint32(opp),-2),mod(opp,4)+1) = ii*4+jj-1; */
              opphes->data[((int32_T)((uint32_T)v2he->data[b_index] >> 2U) +
                            opphes->size[0] * (v2he->data[b_index] -
                ((v2he->data[b_index] >> 2) << 2))) - 1] = ((ii + 1) << 2) + ne;

              /* 'determine_opposite_halfedge_tri:68' found = found + 1; */
              found++;
            }
          }

          /*  Check for consistency */
          /* 'determine_opposite_halfedge_tri:73' if found>1 */
          if ((found > 1) || (found != 0)) {
            /* 'determine_opposite_halfedge_tri:74' error( 'Input mesh is not an oriented manifold.'); */
            guard1 = TRUE;
          } else {
            /* 'determine_opposite_halfedge_tri:75' elseif ~found */
            /* 'determine_opposite_halfedge_tri:76' for index = is_index(v):is_index(v+1)-1 */
            i3 = is_index->data[tris->data[ii + tris->size[0] * ne]] - 1;
            b_index = is_index->data[tris->data[ii + tris->size[0] * ne] - 1];
            do {
              exitg1 = 0U;
              if (b_index <= i3) {
                /* 'determine_opposite_halfedge_tri:77' if v2nv(index)==vn && int32(bitshift( uint32(v2he(index)),-2))~=ii */
                guard2 = FALSE;
                if (v2nv->data[b_index - 1] == tris->data[ii + tris->size[0] *
                    (iv8[ne] - 1)]) {
                  a = (uint32_T)v2he->data[b_index - 1];
                  if ((int32_T)(a >> 2U) != ii + 1) {
                    /* 'determine_opposite_halfedge_tri:78' if nargin==3 */
                    /* 'determine_opposite_halfedge_tri:80' else */
                    /* 'determine_opposite_halfedge_tri:81' opphes = zeros(0,3, 'int32'); */
                    i3 = opphes->size[0] * opphes->size[1];
                    opphes->size[0] = 0;
                    opphes->size[1] = 3;
                    emxEnsureCapacity((emxArray__common *)opphes, i3, (int32_T)
                                      sizeof(int32_T));
                    exitg1 = 2U;
                  } else {
                    guard2 = TRUE;
                  }
                } else {
                  guard2 = TRUE;
                }

                if (guard2 == TRUE) {
                  b_index++;
                }
              } else {
                exitg1 = 1U;
              }
            } while (exitg1 == 0U);

            if (exitg1 == 1U) {
              guard1 = TRUE;
            } else {
              exitg2 = 2U;
            }
          }
        }

        if (guard1 == TRUE) {
          ne++;
        }
      } else {
        ii++;
        exitg2 = 1U;
      }
    } while (exitg2 == 0U);

    if (exitg2 == 1U) {
    } else {
      exitg3 = 1U;
    }
  }

  emxFree_int32_T(&v2he);
  emxFree_int32_T(&v2nv);
  emxFree_int32_T(&is_index);
}

/*
 * function opphfs = determine_opposite_halfface( nv, elems, opphfs)
 * DETERMINE_OPPOSITE_HALFFACE determines the opposite half-face of
 *  each half-face of an oriented, manifold volume mesh with or
 *  without boundary.
 *
 *  OPPHFS = DETERMINE_OPPOSITE_HALFFACE(NV,ELEMS)
 *  OPPHFS = DETERMINE_OPPOSITE_HALFFACE(NV,ELEMS,OPPHFS)
 *  Computes mapping from each half-face to its opposite half-face.
 *
 *  See also DETERMINE_NEXTPAGE_VOL, DETERMINE_INCIDENT_HALFFACES
 */
static void determine_opposite_halfface(int32_T nv, const emxArray_int32_T
  *elems, emxArray_int32_T *opphfs)
{
  uint32_T uv3[2];
  int32_T i8;

  /* 'determine_opposite_halfface:12' if nargin<3 */
  /* 'determine_opposite_halfface:12' opphfs = coder.nullcopy(zeros(size(elems),'int32')); */
  for (i8 = 0; i8 < 2; i8++) {
    uv3[i8] = (uint32_T)elems->size[i8];
  }

  i8 = opphfs->size[0] * opphfs->size[1];
  opphfs->size[0] = (int32_T)uv3[0];
  opphfs->size[1] = 4;
  emxEnsureCapacity((emxArray__common *)opphfs, i8, (int32_T)sizeof(int32_T));

  /* 'determine_opposite_halfface:14' switch size(elems,2) */
  /* 'determine_opposite_halfface:17' case {4,10} % tet */
  /*  tet */
  /* 'determine_opposite_halfface:18' opphfs = determine_opposite_halfface_tet(nv, elems, opphfs); */
  determine_opposite_halfface_tet(nv, elems, opphfs);
}

/*
 * function opphfs = determine_opposite_halfface_tet( nv, elems, opphfs)
 * DETERMINE_OPPOSITE_HALFFACE_TET Determine the opposite half-face.
 *  DETERMINE_OPPOSITE_HALFFACE_TET(NV,ELEMS,OPPHFS) Determines the
 *  opposite half-face.
 *
 *  OPPHFS = DETERMINE_OPPOSITE_HALFFACE_TET(NV,ELEMS)
 *  OPPHFS = DETERMINE_OPPOSITE_HALFFACE_TET(NV,ELEMS,OPPHFS)
 *  computes mapping from each half-face to its opposite half-face.
 *
 *  We assign three bits to local_face_id.
 */
static void determine_opposite_halfface_tet(int32_T nv, const emxArray_int32_T
  *elems, emxArray_int32_T *opphfs)
{
  emxArray_int32_T *is_index;
  int32_T i18;
  int32_T nelems;
  int32_T ii;
  boolean_T exitg2;
  int32_T jj;
  static const int8_T hf_tet[12] = { 1, 1, 2, 3, 3, 2, 3, 1, 2, 4, 4, 4 };

  int32_T mtmp;
  int32_T ix;
  emxArray_int32_T *v2hf;
  emxArray_int32_T *v2oe_v1;
  emxArray_int32_T *v2oe_v2;
  int32_T itmp;
  static const int8_T iv19[3] = { 2, 3, 1 };

  static const int8_T iv20[3] = { 3, 1, 2 };

  int32_T b_index;
  boolean_T exitg1;
  b_emxInit_int32_T(&is_index, 1);

  /*  Note: See http://www.grc.nasa.gov/WWW/cgns/CGNS_docs_current/sids/conv.html for numbering */
  /*        convention of faces. */
  /*  Table for vertices of each face. */
  /* 'determine_opposite_halfface_tet:16' hf_tet    = int32([1 3 2; 1 2 4; 2 3 4; 3 1 4]); */
  /* 'determine_opposite_halfface_tet:18' next = int32([2,3,1]); */
  /* 'determine_opposite_halfface_tet:19' prev = int32([3 1 2]); */
  /* % First, build is_index to store starting position for each vertex. */
  /* 'determine_opposite_halfface_tet:22' is_index = zeros(nv+1,1,'int32'); */
  i18 = is_index->size[0];
  is_index->size[0] = nv + 1;
  emxEnsureCapacity((emxArray__common *)is_index, i18, (int32_T)sizeof(int32_T));
  for (i18 = 0; i18 <= nv; i18++) {
    is_index->data[i18] = 0;
  }

  /* 'determine_opposite_halfface_tet:23' nelems = int32(size(elems,1)); */
  nelems = elems->size[0];

  /* 'determine_opposite_halfface_tet:24' for ii=1:nelems */
  ii = 0;
  exitg2 = 0U;
  while ((exitg2 == 0U) && (ii + 1 <= nelems)) {
    /* 'determine_opposite_halfface_tet:25' if elems(ii,1)==0 */
    if (elems->data[ii] == 0) {
      /* 'determine_opposite_halfface_tet:25' nelems=ii-1; */
      nelems = ii;
      exitg2 = 1U;
    } else {
      /* 'determine_opposite_halfface_tet:27' for jj=1:4 */
      for (jj = 0; jj < 4; jj++) {
        /* 'determine_opposite_halfface_tet:28' vs = elems(ii,hf_tet(jj,:)); */
        /* 'determine_opposite_halfface_tet:29' v = max( vs, [], 2); */
        mtmp = elems->data[ii + elems->size[0] * (hf_tet[jj] - 1)];
        for (ix = 0; ix < 2; ix++) {
          if (elems->data[ii + elems->size[0] * (hf_tet[jj + ((ix + 1) << 2)] -
               1)] > mtmp) {
            mtmp = elems->data[ii + elems->size[0] * (hf_tet[jj + ((ix + 1) << 2)]
              - 1)];
          }
        }

        /* 'determine_opposite_halfface_tet:30' is_index(v+1) = is_index(v+1)+1; */
        is_index->data[mtmp]++;
      }

      ii++;
    }
  }

  /* 'determine_opposite_halfface_tet:33' is_index(1) = 1; */
  is_index->data[0] = 1;

  /* 'determine_opposite_halfface_tet:34' for ii=1:nv */
  for (ii = 1; ii <= nv; ii++) {
    /* 'determine_opposite_halfface_tet:34' is_index(ii+1) = is_index(ii) + is_index(ii+1); */
    is_index->data[ii] += is_index->data[ii - 1];
  }

  b_emxInit_int32_T(&v2hf, 1);
  b_emxInit_int32_T(&v2oe_v1, 1);
  b_emxInit_int32_T(&v2oe_v2, 1);

  /*  v2hf stores mapping from each vertex to half-face ID. */
  /*  v2oe stores mapping from each vertex to the encoding of the opposite */
  /*      edge of each half-face.. */
  /* 'determine_opposite_halfface_tet:39' v2hf = coder.nullcopy(zeros(is_index(nv+1),1,'int32')); */
  i18 = v2hf->size[0];
  v2hf->size[0] = is_index->data[nv];
  emxEnsureCapacity((emxArray__common *)v2hf, i18, (int32_T)sizeof(int32_T));

  /* 'determine_opposite_halfface_tet:40' v2oe_v1 = coder.nullcopy(zeros(is_index(nv+1),1, 'int32')); */
  i18 = v2oe_v1->size[0];
  v2oe_v1->size[0] = is_index->data[nv];
  emxEnsureCapacity((emxArray__common *)v2oe_v1, i18, (int32_T)sizeof(int32_T));

  /* 'determine_opposite_halfface_tet:41' v2oe_v2 = coder.nullcopy(zeros(is_index(nv+1),1, 'int32')); */
  i18 = v2oe_v2->size[0];
  v2oe_v2->size[0] = is_index->data[nv];
  emxEnsureCapacity((emxArray__common *)v2oe_v2, i18, (int32_T)sizeof(int32_T));

  /* 'determine_opposite_halfface_tet:43' for ii=1:nelems */
  for (ii = 0; ii + 1 <= nelems; ii++) {
    /* 'determine_opposite_halfface_tet:44' for jj=int32(1):4 */
    for (jj = 0; jj < 4; jj++) {
      /* 'determine_opposite_halfface_tet:45' vs = elems(ii,hf_tet(jj,:)); */
      /* 'determine_opposite_halfface_tet:46' [v,kk] = max( vs, [], 2); */
      mtmp = elems->data[ii + elems->size[0] * (hf_tet[jj] - 1)] - 1;
      itmp = 0;
      for (ix = 0; ix < 2; ix++) {
        if (elems->data[ii + elems->size[0] * (hf_tet[jj + ((ix + 1) << 2)] - 1)]
            > mtmp + 1) {
          mtmp = elems->data[ii + elems->size[0] * (hf_tet[jj + ((ix + 1) << 2)]
            - 1)] - 1;
          itmp = ix + 1;
        }
      }

      /* 'determine_opposite_halfface_tet:48' v2oe_v1(is_index(v)) = vs( next(kk)); */
      v2oe_v1->data[is_index->data[mtmp] - 1] = elems->data[ii + elems->size[0] *
        (hf_tet[jj + ((iv19[itmp] - 1) << 2)] - 1)];

      /* 'determine_opposite_halfface_tet:49' v2oe_v2(is_index(v)) = vs( prev(kk)); */
      v2oe_v2->data[is_index->data[mtmp] - 1] = elems->data[ii + elems->size[0] *
        (hf_tet[jj + ((iv20[itmp] - 1) << 2)] - 1)];

      /* 'determine_opposite_halfface_tet:50' v2hf(is_index(v)) = ii*8 + jj-1; */
      v2hf->data[is_index->data[mtmp] - 1] = ((ii + 1) << 3) + jj;

      /* 'determine_opposite_halfface_tet:51' is_index(v) = is_index(v)+1; */
      is_index->data[mtmp]++;
    }
  }

  /* 'determine_opposite_halfface_tet:54' for ii=nv-1:-1:1 */
  for (ii = nv - 1; ii > 0; ii--) {
    /* 'determine_opposite_halfface_tet:54' is_index(ii+1) = is_index(ii); */
    is_index->data[ii] = is_index->data[ii - 1];
  }

  /* 'determine_opposite_halfface_tet:55' is_index(1)=1; */
  is_index->data[0] = 1;

  /*  Fill in opphfs for each half-face. */
  /* 'determine_opposite_halfface_tet:58' if nargin<3 || isempty(opphfs) */
  if (opphfs->size[0] == 0) {
    /* 'determine_opposite_halfface_tet:59' opphfs = zeros(size(elems,1),4,'int32'); */
    mtmp = elems->size[0];
    i18 = opphfs->size[0] * opphfs->size[1];
    opphfs->size[0] = mtmp;
    opphfs->size[1] = 4;
    emxEnsureCapacity((emxArray__common *)opphfs, i18, (int32_T)sizeof(int32_T));
    for (i18 = 0; i18 < 4; i18++) {
      ix = mtmp - 1;
      for (b_index = 0; b_index <= ix; b_index++) {
        opphfs->data[b_index + opphfs->size[0] * i18] = 0;
      }
    }
  } else {
    /* 'determine_opposite_halfface_tet:60' else */
    /* 'determine_opposite_halfface_tet:61' assert( size(opphfs,1)>=nelems && size(opphfs,2)>=4); */
    /* 'determine_opposite_halfface_tet:62' opphfs(:) = 0; */
    i18 = opphfs->size[0] * opphfs->size[1];
    opphfs->size[1] = 4;
    emxEnsureCapacity((emxArray__common *)opphfs, i18, (int32_T)sizeof(int32_T));
    for (i18 = 0; i18 < 4; i18++) {
      ix = opphfs->size[0] - 1;
      for (b_index = 0; b_index <= ix; b_index++) {
        opphfs->data[b_index + opphfs->size[0] * i18] = 0;
      }
    }
  }

  /* 'determine_opposite_halfface_tet:65' for ii=1:nelems */
  for (ii = 0; ii + 1 <= nelems; ii++) {
    /* 'determine_opposite_halfface_tet:66' for jj=int32(1):4 */
    for (jj = 0; jj < 4; jj++) {
      /*  local face ID */
      /* 'determine_opposite_halfface_tet:67' if opphfs(ii,jj) */
      if (opphfs->data[ii + opphfs->size[0] * jj] != 0) {
      } else {
        /* 'determine_opposite_halfface_tet:68' vs = elems(ii, hf_tet(jj,:)); */
        /*  list of vertices of face */
        /* 'determine_opposite_halfface_tet:69' [v,imax] = max( vs, [], 2); */
        mtmp = elems->data[ii + elems->size[0] * (hf_tet[jj] - 1)];
        itmp = 0;
        for (ix = 0; ix < 2; ix++) {
          if (elems->data[ii + elems->size[0] * (hf_tet[jj + ((ix + 1) << 2)] -
               1)] > mtmp) {
            mtmp = elems->data[ii + elems->size[0] * (hf_tet[jj + ((ix + 1) << 2)]
              - 1)];
            itmp = ix + 1;
          }
        }

        /* 'determine_opposite_halfface_tet:71' found = false; */
        /* 'determine_opposite_halfface_tet:72' v1 = vs(prev(imax)); */
        /* 'determine_opposite_halfface_tet:72' v2 = vs(next(imax)); */
        /*  Search for opposite half-face. */
        /* 'determine_opposite_halfface_tet:74' for index = is_index( v):is_index( v+1)-1 */
        i18 = is_index->data[mtmp] - 1;
        b_index = is_index->data[mtmp - 1] - 1;
        exitg1 = 0U;
        while ((exitg1 == 0U) && (b_index + 1 <= i18)) {
          /* 'determine_opposite_halfface_tet:75' if v2oe_v1(index) == v1 && v2oe_v2(index) == v2 */
          if ((v2oe_v1->data[b_index] == elems->data[ii + elems->size[0] *
               (hf_tet[jj + ((iv20[itmp] - 1) << 2)] - 1)]) && (v2oe_v2->
               data[b_index] == elems->data[ii + elems->size[0] * (hf_tet[jj +
                ((iv19[itmp] - 1) << 2)] - 1)])) {
            /* 'determine_opposite_halfface_tet:76' opp = v2hf(index); */
            /* 'determine_opposite_halfface_tet:77' opphfs(ii,jj) = opp; */
            opphfs->data[ii + opphfs->size[0] * jj] = v2hf->data[b_index];

            /*  opphfs(hfid2cid(opp),hfid2lfid(opp)) = ii*8+jj-1; */
            /* 'determine_opposite_halfface_tet:80' lfid0=mod(opp,8); */
            /* 'determine_opposite_halfface_tet:80' opphfs(bitshift(uint32(opp),-3),lfid0+1) = ii*8+jj-1; */
            opphfs->data[((int32_T)((uint32_T)v2hf->data[b_index] >> 3U) +
                          opphfs->size[0] * (v2hf->data[b_index] - ((v2hf->
              data[b_index] >> 3) << 3))) - 1] = ((ii + 1) << 3) + jj;

            /* 'determine_opposite_halfface_tet:82' found = true; */
            exitg1 = 1U;
          } else {
            b_index++;
          }
        }

        /* 'determine_opposite_halfface_tet:87' if ~found */
      }
    }
  }

  emxFree_int32_T(&v2oe_v2);
  emxFree_int32_T(&v2oe_v1);
  emxFree_int32_T(&v2hf);
  emxFree_int32_T(&is_index);
}

static int32_T div_s32_floor(int32_T numerator, int32_T denominator)
{
  int32_T quotient;
  uint32_T absNumerator;
  uint32_T absDenominator;
  int32_T quotientNeedsNegation;
  uint32_T tempAbsQuotient;
  if (denominator == 0) {
    quotient = numerator >= 0 ? MAX_int32_T : MIN_int32_T;
  } else {
    absNumerator = (uint32_T)(numerator >= 0 ? numerator : -numerator);
    absDenominator = (uint32_T)(denominator >= 0 ? denominator : -denominator);
    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    tempAbsQuotient = absNumerator / absDenominator;
    if ((uint32_T)quotientNeedsNegation) {
      absNumerator %= absDenominator;
      if (absNumerator > (uint32_T)0) {
        tempAbsQuotient++;
      }
    }

    quotient = (uint32_T)quotientNeedsNegation ? -(int32_T)tempAbsQuotient :
      (int32_T)tempAbsQuotient;
  }

  return quotient;
}

/*
 *
 */
static void e_unique(const emxArray_int32_T *a, emxArray_int32_T *b)
{
  int32_T na;
  int32_T n;
  uint32_T uv5[2];
  int32_T pEnd;
  emxArray_int32_T *idx;
  int32_T np1;
  int32_T k;
  boolean_T p;
  emxArray_int32_T *idx0;
  int32_T j;
  int32_T i;
  int32_T i2;
  int32_T b_p;
  int32_T q;
  int32_T qEnd;
  int32_T kEnd;
  emxArray_int32_T *r9;
  emxArray_int32_T *b_b;
  int32_T iv18[2];
  emxArray_int32_T r10;
  na = a->size[1];
  n = a->size[1];
  for (pEnd = 0; pEnd < 2; pEnd++) {
    uv5[pEnd] = (uint32_T)a->size[pEnd];
  }

  emxInit_int32_T(&idx, 2);
  pEnd = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  idx->size[1] = (int32_T)uv5[1];
  emxEnsureCapacity((emxArray__common *)idx, pEnd, (int32_T)sizeof(int32_T));
  np1 = n + 1;
  if (a->size[1] == 0) {
    for (k = 1; k <= n; k++) {
      idx->data[k - 1] = k;
    }
  } else {
    for (k = 1; k <= n; k++) {
      idx->data[k - 1] = k;
    }

    for (k = 1; k <= n - 1; k += 2) {
      if (a->data[k - 1] <= a->data[k]) {
        p = TRUE;
      } else {
        p = FALSE;
      }

      if (p) {
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    b_emxInit_int32_T(&idx0, 1);
    pEnd = idx0->size[0];
    idx0->size[0] = n;
    emxEnsureCapacity((emxArray__common *)idx0, pEnd, (int32_T)sizeof(int32_T));
    j = n - 1;
    for (pEnd = 0; pEnd <= j; pEnd++) {
      idx0->data[pEnd] = 1;
    }

    i = 2;
    while (i < n) {
      i2 = i << 1;
      j = 1;
      for (pEnd = 1 + i; pEnd < np1; pEnd = qEnd + i) {
        b_p = j;
        q = pEnd;
        qEnd = j + i2;
        if (qEnd > np1) {
          qEnd = np1;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          if (a->data[idx->data[b_p - 1] - 1] <= a->data[idx->data[q - 1] - 1])
          {
            p = TRUE;
          } else {
            p = FALSE;
          }

          if (p) {
            idx0->data[k] = idx->data[b_p - 1];
            b_p++;
            if (b_p == pEnd) {
              while (q < qEnd) {
                k++;
                idx0->data[k] = idx->data[q - 1];
                q++;
              }
            }
          } else {
            idx0->data[k] = idx->data[q - 1];
            q++;
            if (q == qEnd) {
              while (b_p < pEnd) {
                k++;
                idx0->data[k] = idx->data[b_p - 1];
                b_p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k + 1 <= kEnd; k++) {
          idx->data[(j + k) - 1] = idx0->data[k];
        }

        j = qEnd;
      }

      i = i2;
    }

    emxFree_int32_T(&idx0);
  }

  for (pEnd = 0; pEnd < 2; pEnd++) {
    uv5[pEnd] = (uint32_T)a->size[pEnd];
  }

  pEnd = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = (int32_T)uv5[1];
  emxEnsureCapacity((emxArray__common *)b, pEnd, (int32_T)sizeof(int32_T));
  for (k = 0; k + 1 <= na; k++) {
    b->data[k] = a->data[idx->data[k] - 1];
  }

  emxFree_int32_T(&idx);
  i2 = 0;
  k = 1;
  while (k <= na) {
    i = b->data[k - 1];
    do {
      k++;
    } while (!((k > na) || (!(b->data[k - 1] == i))));

    i2++;
    b->data[i2 - 1] = i;
  }

  if (1 > i2) {
    i2 = 0;
  }

  b_emxInit_int32_T(&r9, 1);
  pEnd = r9->size[0];
  r9->size[0] = i2;
  emxEnsureCapacity((emxArray__common *)r9, pEnd, (int32_T)sizeof(int32_T));
  j = i2 - 1;
  for (pEnd = 0; pEnd <= j; pEnd++) {
    r9->data[pEnd] = 1 + pEnd;
  }

  emxInit_int32_T(&b_b, 2);
  iv18[0] = 1;
  iv18[1] = r9->size[0];
  pEnd = b_b->size[0] * b_b->size[1];
  b_b->size[0] = iv18[0];
  b_b->size[1] = iv18[1];
  emxEnsureCapacity((emxArray__common *)b_b, pEnd, (int32_T)sizeof(int32_T));
  j = iv18[1] - 1;
  for (pEnd = 0; pEnd <= j; pEnd++) {
    i = iv18[0] - 1;
    for (i2 = 0; i2 <= i; i2++) {
      r10 = *r9;
      r10.size = (int32_T *)&iv18;
      r10.numDimensions = 1;
      b_b->data[i2 + b_b->size[0] * pEnd] = b->data[r10.data[i2 + r10.size[0] *
        pEnd] - 1];
    }
  }

  emxFree_int32_T(&r9);
  pEnd = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = b_b->size[1];
  emxEnsureCapacity((emxArray__common *)b, pEnd, (int32_T)sizeof(int32_T));
  j = b_b->size[1] - 1;
  for (pEnd = 0; pEnd <= j; pEnd++) {
    b->data[b->size[0] * pEnd] = b_b->data[b_b->size[0] * pEnd];
  }

  emxFree_int32_T(&b_b);
}

/*
 *
 */
static void eml_null_assignment(emxArray_int32_T *x, int32_T idx)
{
  int32_T nrows;
  int32_T j;
  int32_T loop_ub;
  emxArray_int32_T *b_x;
  int32_T i14;
  nrows = x->size[0] - 1;
  for (j = 0; j < 6; j++) {
    for (loop_ub = idx; loop_ub <= nrows; loop_ub++) {
      x->data[(loop_ub + x->size[0] * j) - 1] = x->data[loop_ub + x->size[0] * j];
    }
  }

  if (1 > nrows) {
    nrows = 0;
  }

  emxInit_int32_T(&b_x, 2);
  j = b_x->size[0] * b_x->size[1];
  b_x->size[0] = nrows;
  b_x->size[1] = 6;
  emxEnsureCapacity((emxArray__common *)b_x, j, (int32_T)sizeof(int32_T));
  for (j = 0; j < 6; j++) {
    loop_ub = nrows - 1;
    for (i14 = 0; i14 <= loop_ub; i14++) {
      b_x->data[i14 + b_x->size[0] * j] = x->data[i14 + x->size[0] * j];
    }
  }

  j = x->size[0] * x->size[1];
  x->size[0] = b_x->size[0];
  x->size[1] = 6;
  emxEnsureCapacity((emxArray__common *)x, j, (int32_T)sizeof(int32_T));
  for (j = 0; j < 6; j++) {
    loop_ub = b_x->size[0] - 1;
    for (i14 = 0; i14 <= loop_ub; i14++) {
      x->data[i14 + x->size[0] * j] = b_x->data[i14 + b_x->size[0] * j];
    }
  }

  emxFree_int32_T(&b_x);
}

/*
 *
 */
static boolean_T eml_sort_le(const emxArray_int32_T *v, int32_T col, int32_T
  irow1, int32_T irow2)
{
  boolean_T p;
  int32_T abscolk;
  boolean_T b1;
  p = TRUE;
  if (col < 0) {
    abscolk = -col;
  } else {
    abscolk = col;
  }

  abscolk = (abscolk - 1) * v->size[0] - 1;
  if (v->data[abscolk + irow1] == v->data[abscolk + irow2]) {
    b1 = TRUE;
  } else {
    b1 = FALSE;
  }

  if (!b1) {
    if (col < 0) {
      if (v->data[abscolk + irow1] >= v->data[abscolk + irow2]) {
        p = TRUE;
      } else {
        p = FALSE;
      }
    } else if (v->data[abscolk + irow1] <= v->data[abscolk + irow2]) {
      p = TRUE;
    } else {
      p = FALSE;
    }
  }

  return p;
}

/*
 *
 */
static real_T eml_xdot(int32_T n, const real_T x[3], int32_T ix0, int32_T incx,
  const real_T y[3], int32_T iy0, int32_T incy)
{
  real_T d;
  int32_T ix;
  int32_T iy;
  int32_T k;
  d = 0.0;
  if (n < 1) {
  } else {
    ix = ix0;
    iy = iy0;
    for (k = 1; k <= n; k++) {
      d += x[ix - 1] * y[iy - 1];
      if (incx < 0) {
        ix += incx;
      } else {
        ix += incx;
      }

      if (incy < 0) {
        iy += incy;
      } else {
        iy += incy;
      }
    }
  }

  return d;
}

static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize)
{
  int32_T newNumel;
  int32_T loop_ub;
  int32_T i;
  void *newData;
  newNumel = 1;
  loop_ub = emxArray->numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    if (emxArray->allocatedSize) {
      loop_ub = emxArray->allocatedSize;
      if (loop_ub < 16) {
        loop_ub = 16;
       }

      while (loop_ub < newNumel) {
        loop_ub <<= 1;
      }
    } else {
      loop_ub = newNumel;
    }

    newData = calloc((uint32_T)loop_ub, (uint32_T)elementSize);
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, (uint32_T)(elementSize * oldNumel));
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }

    emxArray->data = newData;
    emxArray->allocatedSize = loop_ub;
    emxArray->canFreeData = TRUE;
  }
}

static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_boolean_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_boolean_T *)NULL;
  }
}

static void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int32_T
  numDimensions)
{
  emxArray_boolean_T *emxArray;
  int32_T loop_ub;
  int32_T i;
  *pEmxArray = (emxArray_boolean_T *)malloc(sizeof(emxArray_boolean_T));
  emxArray = *pEmxArray;
  emxArray->data = (boolean_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    emxArray->size[i] = 0;
  }
}

static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions)
{
  emxArray_int32_T *emxArray;
  int32_T loop_ub;
  int32_T i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int32_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    emxArray->size[i] = 0;
  }
}

static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions)
{
  emxArray_real_T *emxArray;
  int32_T loop_ub;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * function [Face_ps,Face_elems] = extract_cellface_edges(FaceId,Faces,FaceEdgeMap,cor_map,ps_index,ps_local,ps_localID,tris_labelled,opphes)
 *  Face_ps : nx1 matrix containing the id's of all points (including
 *  the corners) in ps_local that lie on the given face.
 *  Face_elems : mx2 matrix containing the ordered list of edges(includes cell edge discretization as well).
 */
static void extract_cellface_edges(real_T FaceId, const int32_T Faces[4], const
  int32_T FaceEdgeMap[4], const emxArray_int32_T *cor_map, int32_T ps_index,
  const emxArray_real_T *ps_local, const emxArray_int32_T *ps_localID, const
  emxArray_int32_T *tris_labelled, const emxArray_int32_T *opphes,
  emxArray_int32_T *Face_ps, emxArray_int32_T *Face_elems)
{
  int32_T i5;
  int32_T loop_ub;
  int32_T f_idx;
  int32_T ps_local_idx_0;
  int32_T fe_idx;
  int32_T ii;
  int32_T b_index;
  int32_T common_corner[4];
  int32_T count;
  int32_T jj;
  boolean_T exitg3;
  emxArray_int32_T *r2;
  emxArray_int32_T *r3;
  emxArray_int32_T *r4;
  emxArray_int32_T *b_common_corner;
  int32_T exitg1;
  int32_T iv15[2];
  emxArray_int32_T b_Face_ps;
  boolean_T exitg2;
  real_T d2;
  emxArray_int32_T *c_Face_ps;
  int32_T d_Face_ps[2];
  emxArray_int32_T *b_Face_elems;
  emxArray_int32_T *c_Face_elems;

  /* 'extract_cellface_edges:6' Face_ps = zeros(size(ps_local,1),1,'int32'); */
  i5 = Face_ps->size[0];
  Face_ps->size[0] = ps_local->size[0];
  emxEnsureCapacity((emxArray__common *)Face_ps, i5, (int32_T)sizeof(int32_T));
  loop_ub = ps_local->size[0] - 1;
  for (i5 = 0; i5 <= loop_ub; i5++) {
    Face_ps->data[i5] = 0;
  }

  /* 'extract_cellface_edges:7' f_idx = int32(1); */
  f_idx = 0;

  /* 'extract_cellface_edges:8' Face_elems = zeros(size(ps_local,1),2,'int32'); */
  ps_local_idx_0 = ps_local->size[0];
  i5 = Face_elems->size[0] * Face_elems->size[1];
  Face_elems->size[0] = ps_local_idx_0;
  Face_elems->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)Face_elems, i5, (int32_T)sizeof(int32_T));
  loop_ub = (ps_local->size[0] << 1) - 1;
  for (i5 = 0; i5 <= loop_ub; i5++) {
    Face_elems->data[i5] = 0;
  }

  /* 'extract_cellface_edges:8' fe_idx = int32(1); */
  fe_idx = 1;

  /*  COLLECT ID'S OF VERTICES LYING ON THE FACE */
  /*  Points From Embedded Surface */
  /* 'extract_cellface_edges:11' for ii=1:ps_index */
  for (ii = 0; ii + 1 <= ps_index; ii++) {
    /* 'extract_cellface_edges:12' if (ps_localID(ii,4)== FaceId) */
    if ((real_T)ps_localID->data[ii + ps_localID->size[0] * 3] == FaceId) {
      /* 'extract_cellface_edges:13' Face_ps(f_idx) = ii; */
      Face_ps->data[f_idx] = ii + 1;

      /* 'extract_cellface_edges:14' f_idx = f_idx + 1; */
      f_idx++;
    } else if ((ps_localID->data[ii + (ps_localID->size[0] << 2)] ==
                FaceEdgeMap[0]) || (ps_localID->data[ii + (ps_localID->size[0] <<
      2)] == FaceEdgeMap[1]) || (ps_localID->data[ii + (ps_localID->size[0] << 2)]
                == FaceEdgeMap[2]) || (ps_localID->data[ii + (ps_localID->size[0]
      << 2)] == FaceEdgeMap[3])) {
      /* 'extract_cellface_edges:15' elseif ((ps_localID(ii,5) == FaceEdgeMap(1))||(ps_localID(ii,5) == FaceEdgeMap(2))|| ... */
      /* 'extract_cellface_edges:16'             (ps_localID(ii,5) == FaceEdgeMap(3))||(ps_localID(ii,5) == FaceEdgeMap(4))) */
      /* 'extract_cellface_edges:17' Face_ps(f_idx) = ii; */
      Face_ps->data[f_idx] = ii + 1;

      /* 'extract_cellface_edges:18' f_idx = f_idx + 1; */
      f_idx++;
    } else {
      if ((ps_localID->data[ii + ps_localID->size[0] * 5] == Faces[0]) ||
          (ps_localID->data[ii + ps_localID->size[0] * 5] == Faces[1]) ||
          (ps_localID->data[ii + ps_localID->size[0] * 5] == Faces[2]) ||
          (ps_localID->data[ii + ps_localID->size[0] * 5] == Faces[3])) {
        /* 'extract_cellface_edges:19' elseif ((ps_localID(ii,6) == Faces(1))||(ps_localID(ii,6) == Faces(2))|| ... */
        /* 'extract_cellface_edges:20'             (ps_localID(ii,6) == Faces(3))||(ps_localID(ii,6) == Faces(4))) */
        /* 'extract_cellface_edges:21' Face_ps(f_idx) = ii; */
        Face_ps->data[f_idx] = ii + 1;

        /* 'extract_cellface_edges:22' f_idx = f_idx + 1; */
        f_idx++;
      }
    }
  }

  /* 'extract_cellface_edges:25' index = f_idx - 1; */
  b_index = f_idx;

  /* 'extract_cellface_edges:26' if index~=0 */
  if (f_idx != 0) {
    /* 'extract_cellface_edges:27' [Face_elems,fe_idx] = extract_trisedges_on_cellface(Face_ps,index,tris_labelled,opphes,Face_elems,fe_idx); */
    ps_local_idx_0 = ps_local->size[0];
    i5 = Face_elems->size[0] * Face_elems->size[1];
    Face_elems->size[0] = ps_local_idx_0;
    Face_elems->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)Face_elems, i5, (int32_T)sizeof
                      (int32_T));
    loop_ub = (ps_local->size[0] << 1) - 1;
    for (i5 = 0; i5 <= loop_ub; i5++) {
      Face_elems->data[i5] = 0;
    }

    fe_idx = extract_trisedges_on_cellface(Face_ps, f_idx, tris_labelled, opphes,
      Face_elems);
  }

  /*  ADD THE CORNERS AND THE ELEMS ON THE CELL EDGES */
  /* 'extract_cellface_edges:30' if (size(cor_map,1)==8) */
  if (cor_map->size[0] == 8) {
    /* 'extract_cellface_edges:31' for ii=1:4 */
    for (ii = 0; ii < 4; ii++) {
      /* 'extract_cellface_edges:32' Face_ps(f_idx) = Faces(ii)+ps_index; */
      Face_ps->data[f_idx] = Faces[ii] + ps_index;

      /* 'extract_cellface_edges:33' f_idx = f_idx + 1; */
      f_idx++;
    }
  } else {
    /* 'extract_cellface_edges:35' else */
    /*  common_corner = intersect(cor_map',Faces); */
    /* 'extract_cellface_edges:37' common_corner = getcommoncorners(cor_map,Faces); */
    /* 'extract_cellface_edges:172' common_corner = coder.nullcopy(zeros(1,4,'int32')); */
    /* 'extract_cellface_edges:172' count = 1; */
    count = 0;

    /* 'extract_cellface_edges:173' for i=1:4 */
    for (ps_local_idx_0 = 0; ps_local_idx_0 < 4; ps_local_idx_0++) {
      /* 'extract_cellface_edges:174' for j=1:size(cor_map,1) */
      jj = 0;
      exitg3 = 0U;
      while ((exitg3 == 0U) && (jj <= cor_map->size[0] - 1)) {
        /* 'extract_cellface_edges:175' if (cor_map(j)==Faces(i)) */
        if (cor_map->data[jj] == Faces[ps_local_idx_0]) {
          /* 'extract_cellface_edges:176' common_corner(count) = cor_map(j); */
          common_corner[count] = cor_map->data[jj];

          /* 'extract_cellface_edges:177' count = count + 1; */
          count++;
          exitg3 = 1U;
        } else {
          jj++;
        }
      }
    }

    /* 'extract_cellface_edges:182' common_corner = common_corner(1:count-1); */
    if (1 > count) {
      count = 0;
    }

    /* 'extract_cellface_edges:38' for ii=1:size(common_corner,2) */
    ii = 0;
    emxInit_int32_T(&r2, 2);
    b_emxInit_int32_T(&r3, 1);
    b_emxInit_int32_T(&r4, 1);
    emxInit_int32_T(&b_common_corner, 2);
    do {
      exitg1 = 0U;
      i5 = r3->size[0];
      r3->size[0] = count;
      emxEnsureCapacity((emxArray__common *)r3, i5, (int32_T)sizeof(int32_T));
      loop_ub = count - 1;
      for (i5 = 0; i5 <= loop_ub; i5++) {
        r3->data[i5] = 1 + i5;
      }

      iv15[0] = 1;
      iv15[1] = r3->size[0];
      i5 = b_common_corner->size[0] * b_common_corner->size[1];
      b_common_corner->size[0] = iv15[0];
      b_common_corner->size[1] = iv15[1];
      emxEnsureCapacity((emxArray__common *)b_common_corner, i5, (int32_T)sizeof
                        (int32_T));
      loop_ub = iv15[1] - 1;
      for (i5 = 0; i5 <= loop_ub; i5++) {
        ps_local_idx_0 = iv15[0] - 1;
        for (jj = 0; jj <= ps_local_idx_0; jj++) {
          b_Face_ps = *r3;
          b_Face_ps.size = (int32_T *)&iv15;
          b_Face_ps.numDimensions = 1;
          b_common_corner->data[jj + b_common_corner->size[0] * i5] =
            common_corner[b_Face_ps.data[jj + b_Face_ps.size[0] * i5] - 1];
        }
      }

      if (ii <= b_common_corner->size[1] - 1) {
        /* 'extract_cellface_edges:39' for jj=1:size(cor_map,1) */
        jj = 0;
        exitg2 = 0U;
        while ((exitg2 == 0U) && (jj <= cor_map->size[0] - 1)) {
          /* 'extract_cellface_edges:40' if (cor_map(jj)==common_corner(ii)) */
          i5 = r4->size[0];
          r4->size[0] = count;
          emxEnsureCapacity((emxArray__common *)r4, i5, (int32_T)sizeof(int32_T));
          loop_ub = count - 1;
          for (i5 = 0; i5 <= loop_ub; i5++) {
            r4->data[i5] = 1 + i5;
          }

          i5 = r2->size[0] * r2->size[1];
          r2->size[0] = 1;
          emxEnsureCapacity((emxArray__common *)r2, i5, (int32_T)sizeof(int32_T));
          ps_local_idx_0 = r4->size[0];
          i5 = r2->size[0] * r2->size[1];
          r2->size[1] = ps_local_idx_0;
          emxEnsureCapacity((emxArray__common *)r2, i5, (int32_T)sizeof(int32_T));
          loop_ub = r4->size[0] - 1;
          for (i5 = 0; i5 <= loop_ub; i5++) {
            r2->data[i5] = r4->data[i5];
          }

          if (cor_map->data[jj] == common_corner[r2->data[r2->size[0] * ii] - 1])
          {
            /* 'extract_cellface_edges:41' Face_ps(f_idx) = jj+ps_index; */
            d2 = (1.0 + (real_T)jj) + (real_T)ps_index;
            d2 = d2 < 0.0 ? ceil(d2 - 0.5) : floor(d2 + 0.5);
            Face_ps->data[f_idx] = (int32_T)d2;

            /* 'extract_cellface_edges:42' f_idx = f_idx + 1; */
            f_idx++;
            exitg2 = 1U;
          } else {
            jj++;
          }
        }

        ii++;
      } else {
        exitg1 = 1U;
      }
    } while (exitg1 == 0U);

    emxFree_int32_T(&b_common_corner);
    emxFree_int32_T(&r4);
    emxFree_int32_T(&r3);
    emxFree_int32_T(&r2);
  }

  /* 'extract_cellface_edges:49' if  (size(cor_map,1)==8) */
  if (cor_map->size[0] == 8) {
    /* 'extract_cellface_edges:50' [Face_elems,fe_idx] = get_celledgelems_eight(Faces,FaceEdgeMap,index,Face_ps,Face_elems,fe_idx,ps_index,ps_local,ps_localID); */
    get_celledgelems_eight(Faces, FaceEdgeMap, b_index, Face_ps, Face_elems,
      &fe_idx, ps_index, ps_local, ps_localID);
  } else {
    /* 'extract_cellface_edges:51' else */
    /* 'extract_cellface_edges:52' [Face_elems,fe_idx] = get_celledgelems_noteight(Faces,FaceEdgeMap,cor_map,index,Face_ps,Face_elems,fe_idx,ps_index,ps_local,ps_localID); */
    get_celledgelems_noteight(Faces, FaceEdgeMap, cor_map, b_index, Face_ps,
      Face_elems, &fe_idx, ps_index, ps_local, ps_localID);
  }

  /* 'extract_cellface_edges:55' Face_ps = Face_ps(1:f_idx-1,:); */
  if (1 > f_idx) {
    f_idx = 0;
  }

  b_emxInit_int32_T(&c_Face_ps, 1);
  ps_local_idx_0 = Face_ps->size[0];
  d_Face_ps[0] = ps_local_idx_0;
  d_Face_ps[1] = 1;
  i5 = c_Face_ps->size[0];
  c_Face_ps->size[0] = f_idx;
  emxEnsureCapacity((emxArray__common *)c_Face_ps, i5, (int32_T)sizeof(int32_T));
  loop_ub = f_idx - 1;
  for (i5 = 0; i5 <= loop_ub; i5++) {
    b_Face_ps = *Face_ps;
    b_Face_ps.size = (int32_T *)&d_Face_ps;
    b_Face_ps.numDimensions = 1;
    c_Face_ps->data[i5] = b_Face_ps.data[i5];
  }

  i5 = Face_ps->size[0];
  Face_ps->size[0] = c_Face_ps->size[0];
  emxEnsureCapacity((emxArray__common *)Face_ps, i5, (int32_T)sizeof(int32_T));
  loop_ub = c_Face_ps->size[0] - 1;
  for (i5 = 0; i5 <= loop_ub; i5++) {
    Face_ps->data[i5] = c_Face_ps->data[i5];
  }

  emxFree_int32_T(&c_Face_ps);

  /* 'extract_cellface_edges:56' Face_elems = Face_elems(1:fe_idx-1,:); */
  i5 = fe_idx - 1;
  if (1 > i5) {
    i5 = 0;
  }

  emxInit_int32_T(&b_Face_elems, 2);
  jj = b_Face_elems->size[0] * b_Face_elems->size[1];
  b_Face_elems->size[0] = i5;
  b_Face_elems->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)b_Face_elems, jj, (int32_T)sizeof
                    (int32_T));
  for (jj = 0; jj < 2; jj++) {
    loop_ub = i5 - 1;
    for (ps_local_idx_0 = 0; ps_local_idx_0 <= loop_ub; ps_local_idx_0++) {
      b_Face_elems->data[ps_local_idx_0 + b_Face_elems->size[0] * jj] =
        Face_elems->data[ps_local_idx_0 + Face_elems->size[0] * jj];
    }
  }

  i5 = Face_elems->size[0] * Face_elems->size[1];
  Face_elems->size[0] = b_Face_elems->size[0];
  Face_elems->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)Face_elems, i5, (int32_T)sizeof(int32_T));
  for (i5 = 0; i5 < 2; i5++) {
    loop_ub = b_Face_elems->size[0] - 1;
    for (jj = 0; jj <= loop_ub; jj++) {
      Face_elems->data[jj + Face_elems->size[0] * i5] = b_Face_elems->data[jj +
        b_Face_elems->size[0] * i5];
    }
  }

  emxFree_int32_T(&b_Face_elems);
  emxInit_int32_T(&c_Face_elems, 2);

  /* 'extract_cellface_edges:57' Face_elems = unique(Face_elems,'rows'); */
  i5 = c_Face_elems->size[0] * c_Face_elems->size[1];
  c_Face_elems->size[0] = Face_elems->size[0];
  c_Face_elems->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)c_Face_elems, i5, (int32_T)sizeof
                    (int32_T));
  loop_ub = Face_elems->size[0] * Face_elems->size[1] - 1;
  for (i5 = 0; i5 <= loop_ub; i5++) {
    c_Face_elems->data[i5] = Face_elems->data[i5];
  }

  b_unique(c_Face_elems, Face_elems);
  emxFree_int32_T(&c_Face_elems);
}

/*
 * function [Face_elems,fe_idx] = extract_trisedges_on_cellface(Face_ps,npnts,tris,opphes,Face_elems,fe_idx)
 */
static int32_T extract_trisedges_on_cellface(const emxArray_int32_T *Face_ps,
  int32_T npnts, const emxArray_int32_T *tris, const emxArray_int32_T *opphes,
  emxArray_int32_T *Face_elems)
{
  int32_T fe_idx;
  emxArray_int32_T *visited;
  int32_T count;
  int32_T ngbtriID;
  emxArray_int32_T *stack;
  int32_T stack_idx;
  int32_T counter;
  int32_T jj;
  int32_T kk;
  static const int8_T emap[6] = { 1, 2, 3, 2, 3, 1 };

  int32_T next;
  boolean_T exitg2;
  uint32_T b_opphes[3];
  uint32_T uv7[3];
  int32_T b_ngbtriID[3];
  boolean_T exitg1;
  uint32_T a;
  b_emxInit_int32_T(&visited, 1);
  fe_idx = 1;

  /* 'extract_trisedges_on_cellface:3' ntris = size(tris,1); */
  /* 'extract_trisedges_on_cellface:4' emap = int32([1 2; 2 3; 3 1]); */
  /* 'extract_trisedges_on_cellface:5' visited = zeros(ntris,1,'int32'); */
  count = visited->size[0];
  visited->size[0] = tris->size[0];
  emxEnsureCapacity((emxArray__common *)visited, count, (int32_T)sizeof(int32_T));
  ngbtriID = tris->size[0] - 1;
  for (count = 0; count <= ngbtriID; count++) {
    visited->data[count] = 0;
  }

  b_emxInit_int32_T(&stack, 1);

  /* 'extract_trisedges_on_cellface:6' stack = zeros(ntris,1,'int32'); */
  count = stack->size[0];
  stack->size[0] = tris->size[0];
  emxEnsureCapacity((emxArray__common *)stack, count, (int32_T)sizeof(int32_T));
  ngbtriID = tris->size[0] - 1;
  for (count = 0; count <= ngbtriID; count++) {
    stack->data[count] = 0;
  }

  /* 'extract_trisedges_on_cellface:6' stack_idx = int32(1); */
  stack_idx = 1;

  /* 'extract_trisedges_on_cellface:7' counter = int32(1); */
  counter = 1;

  /*  SET THE FIRST ELEMENT OF THE STACK */
  /* 'extract_trisedges_on_cellface:10' stack(counter)=1; */
  stack->data[0] = 1;

  /* 'extract_trisedges_on_cellface:11' for jj=1:3 */
  for (jj = 0; jj < 3; jj++) {
    /* 'extract_trisedges_on_cellface:12' c1 = int32(0); */
    ngbtriID = 0;

    /* 'extract_trisedges_on_cellface:12' c2 = int32(0); */
    count = 0;

    /* 'extract_trisedges_on_cellface:13' for kk=1:npnts */
    for (kk = 0; kk + 1 <= npnts; kk++) {
      /* 'extract_trisedges_on_cellface:14' if (tris(1,emap(jj,1))==Face_ps(kk)) */
      if (tris->data[tris->size[0] * (emap[jj] - 1)] == Face_ps->data[kk]) {
        /* 'extract_trisedges_on_cellface:15' c1 = int32(1); */
        ngbtriID = 1;
      }

      /* 'extract_trisedges_on_cellface:17' if (tris(1,emap(jj,2))==Face_ps(kk)) */
      if (tris->data[tris->size[0] * (emap[3 + jj] - 1)] == Face_ps->data[kk]) {
        /* 'extract_trisedges_on_cellface:18' c2 = int32(1); */
        count = 1;
      }
    }

    /* 'extract_trisedges_on_cellface:21' if (c1==1)&&(c2==1) */
    if ((ngbtriID == 1) && (count == 1)) {
      /* 'extract_trisedges_on_cellface:22' Face_elems(fe_idx,:)=tris(1,emap(jj,:)); */
      for (count = 0; count < 2; count++) {
        Face_elems->data[(fe_idx + Face_elems->size[0] * count) - 1] =
          tris->data[tris->size[0] * (emap[jj + 3 * count] - 1)];
      }

      /* 'extract_trisedges_on_cellface:23' fe_idx = fe_idx + 1; */
      fe_idx++;
    }
  }

  /* 'extract_trisedges_on_cellface:26' visited(1) = int32(1); */
  visited->data[0] = 1;

  /*  START THE TRAVERSAL */
  /* 'extract_trisedges_on_cellface:29' while(1) */
  do {
    /* 'extract_trisedges_on_cellface:30' next = int32(0); */
    next = -1;

    /* 'extract_trisedges_on_cellface:31' if (stack_idx ==0) && (counter ~= ntris) */
    if ((stack_idx == 0) && (counter != tris->size[0])) {
      /* 'extract_trisedges_on_cellface:32' for ii=1:ntris */
      ngbtriID = 0;
      exitg2 = 0U;
      while ((exitg2 == 0U) && (ngbtriID <= tris->size[0] - 1)) {
        /* 'extract_trisedges_on_cellface:33' count = int32(0); */
        count = 0;

        /* 'extract_trisedges_on_cellface:34' for jj=1:counter */
        for (jj = 1; jj <= counter; jj++) {
          /* 'extract_trisedges_on_cellface:35' if (ii~=stack(jj)) */
          if (1 + ngbtriID != stack->data[jj - 1]) {
            /* 'extract_trisedges_on_cellface:36' count = count + 1; */
            count++;
          }
        }

        /* 'extract_trisedges_on_cellface:39' if (count ==counter) */
        if (count == counter) {
          /* 'extract_trisedges_on_cellface:40' next = int32(ii); */
          next = ngbtriID;
          exitg2 = 1U;
        } else {
          ngbtriID++;
        }
      }
    } else {
      /* 'extract_trisedges_on_cellface:44' else */
      /* 'extract_trisedges_on_cellface:45' top = stack(stack_idx); */
      /* 'extract_trisedges_on_cellface:46' next = getnexttri(top,opphes,visited,next); */
      /* 'extract_trisedges_on_cellface:85' ngbtriID = int32( bitshift( uint32(opphes(triID,:)),-2)); */
      ngbtriID = stack->data[stack_idx - 1];
      for (count = 0; count < 3; count++) {
        b_opphes[count] = (uint32_T)opphes->data[(ngbtriID + opphes->size[0] *
          count) - 1];
      }

      bitshift(b_opphes, -2.0, uv7);
      for (count = 0; count < 3; count++) {
        b_ngbtriID[count] = (int32_T)uv7[count];
      }

      /* 'extract_trisedges_on_cellface:86' count = int32(0); */
      count = 0;

      /* 'extract_trisedges_on_cellface:87' for kk=1:3 */
      kk = 0;
      exitg1 = 0U;
      while ((exitg1 == 0U) && (kk < 3)) {
        /* 'extract_trisedges_on_cellface:88' if (ngbtriID(kk) ~=0) && (visited(ngbtriID(kk)) ==0) */
        if ((b_ngbtriID[kk] != 0) && (visited->data[b_ngbtriID[kk] - 1] == 0)) {
          /* 'extract_trisedges_on_cellface:89' next = ngbtriID(kk); */
          next = b_ngbtriID[kk] - 1;
          exitg1 = 1U;
        } else {
          if ((b_ngbtriID[kk] == 0) || (visited->data[b_ngbtriID[kk] - 1] == 1))
          {
            /* 'extract_trisedges_on_cellface:91' elseif (ngbtriID(kk) ==0)||(visited(ngbtriID(kk)) ==1) */
            /* 'extract_trisedges_on_cellface:92' count = count+1; */
            count++;
          }

          kk++;
        }
      }

      /* 'extract_trisedges_on_cellface:95' if (count==3) */
      if (count == 3) {
        /* 'extract_trisedges_on_cellface:96' next = int32(-1); */
        next = -2;
      }
    }

    /* 'extract_trisedges_on_cellface:49' if (next == -1) */
    if (next + 1 == -1) {
      /* 'extract_trisedges_on_cellface:50' stack_idx = stack_idx - 1; */
      stack_idx--;
    } else {
      /* 'extract_trisedges_on_cellface:51' else */
      /* 'extract_trisedges_on_cellface:52' counter = counter + 1; */
      counter++;

      /* 'extract_trisedges_on_cellface:53' stack_idx = counter; */
      stack_idx = counter;

      /* 'extract_trisedges_on_cellface:54' stack(counter) = next; */
      stack->data[counter - 1] = next + 1;

      /*  display(stack(1:counter)) */
      /* 'extract_trisedges_on_cellface:57' for jj=1:3 */
      for (jj = 0; jj < 3; jj++) {
        /* 'extract_trisedges_on_cellface:58' ngbtriID = int32( bitshift( uint32(opphes(next,jj)),-2)); */
        a = (uint32_T)opphes->data[next + opphes->size[0] * jj];
        ngbtriID = (int32_T)(a >> 2U);

        /* 'extract_trisedges_on_cellface:59' if (ngbtriID == 0) || (visited(ngbtriID)==0) */
        if ((ngbtriID == 0) || (visited->data[ngbtriID - 1] == 0)) {
          /* 'extract_trisedges_on_cellface:60' c1 = int32(0); */
          ngbtriID = 0;

          /* 'extract_trisedges_on_cellface:60' c2 = int32(0); */
          count = 0;

          /* 'extract_trisedges_on_cellface:61' for kk=1:npnts */
          for (kk = 0; kk + 1 <= npnts; kk++) {
            /* 'extract_trisedges_on_cellface:62' if (tris(next,emap(jj,1))==Face_ps(kk)) */
            if (tris->data[next + tris->size[0] * (emap[jj] - 1)] ==
                Face_ps->data[kk]) {
              /* 'extract_trisedges_on_cellface:63' c1 = int32(1); */
              ngbtriID = 1;
            }

            /* 'extract_trisedges_on_cellface:65' if (tris(next,emap(jj,2))==Face_ps(kk)) */
            if (tris->data[next + tris->size[0] * (emap[3 + jj] - 1)] ==
                Face_ps->data[kk]) {
              /* 'extract_trisedges_on_cellface:66' c2 = int32(1); */
              count = 1;
            }
          }

          /* 'extract_trisedges_on_cellface:69' if (c1==1)&&(c2==1) */
          if ((ngbtriID == 1) && (count == 1)) {
            /* 'extract_trisedges_on_cellface:70' Face_elems(fe_idx,:)=tris(next,emap(jj,:)); */
            for (count = 0; count < 2; count++) {
              Face_elems->data[(fe_idx + Face_elems->size[0] * count) - 1] =
                tris->data[next + tris->size[0] * (emap[jj + 3 * count] - 1)];
            }

            /* 'extract_trisedges_on_cellface:71' fe_idx = fe_idx + 1; */
            fe_idx++;
          }
        }
      }

      /* 'extract_trisedges_on_cellface:75' visited(next) = 1; */
      visited->data[next] = 1;
    }

    /* 'extract_trisedges_on_cellface:78' if (counter == ntris) */
  } while (!(counter == tris->size[0]));

  emxFree_int32_T(&stack);
  emxFree_int32_T(&visited);
  return fe_idx;
}

/*
 * function [Face_elems,fe_idx] = get_celledgelems_eight(Faces,FaceEdgeMap,index,Face_ps,Face_elems,fe_idx,ps_index,ps_local,ps_localID)
 */
static void get_celledgelems_eight(const int32_T Faces[4], const int32_T
  FaceEdgeMap[4], int32_T b_index, const emxArray_int32_T *Face_ps,
  emxArray_int32_T *Face_elems, int32_T *fe_idx, int32_T ps_index, const
  emxArray_real_T *ps_local, const emxArray_int32_T *ps_localID)
{
  emxArray_real_T *dist;
  emxArray_int32_T *idx;
  emxArray_real_T *vwork;
  emxArray_int32_T *iidx;
  emxArray_int32_T *idx0;
  int32_T ii;
  int32_T i;
  int32_T P[10];
  int32_T count;
  int32_T j;
  static const int8_T map[8] = { 1, 2, 3, 4, 2, 3, 4, 1 };

  real_T a;
  real_T b_a;
  real_T c_a;
  real_T y;
  int32_T d_a;
  int8_T unnamed_idx_0;
  int32_T i1;
  int32_T k;
  int32_T n;
  int32_T np1;
  boolean_T p;
  int32_T i2;
  int32_T pEnd;
  int32_T b_p;
  int32_T q;
  int32_T qEnd;
  int32_T kEnd;

  /* 'extract_cellface_edges:62' map = [1 2; 2 3; 3 4; 4 1]; */
  /* 'extract_cellface_edges:63' for ii=1:4 */
  b_emxInit_real_T(&dist, 1);
  b_emxInit_int32_T(&idx, 1);
  b_emxInit_real_T(&vwork, 1);
  b_emxInit_int32_T(&iidx, 1);
  b_emxInit_int32_T(&idx0, 1);
  for (ii = 0; ii < 4; ii++) {
    /* 'extract_cellface_edges:64' P = zeros(10,1); */
    for (i = 0; i < 10; i++) {
      P[i] = 0;
    }

    /* 'extract_cellface_edges:64' count = 1; */
    count = -1;

    /* 'extract_cellface_edges:65' for jj=1:index */
    for (i = 0; i + 1 <= b_index; i++) {
      /* 'extract_cellface_edges:66' if (ps_localID(Face_ps(jj),5) == FaceEdgeMap(ii)) */
      if (ps_localID->data[(Face_ps->data[i] + (ps_localID->size[0] << 2)) - 1] ==
          FaceEdgeMap[ii]) {
        /* 'extract_cellface_edges:67' P(count) = Face_ps(jj); */
        P[count + 1] = Face_ps->data[i];

        /* 'extract_cellface_edges:68' count = count + 1; */
        count++;
      }
    }

    /* 'extract_cellface_edges:71' if (count == 1) */
    if (count + 2 == 1) {
      /* 'extract_cellface_edges:72' Face_elems(fe_idx,:) = Faces(map(ii,:))+ps_index; */
      i = *fe_idx;
      for (j = 0; j < 2; j++) {
        Face_elems->data[(i + Face_elems->size[0] * j) - 1] = Faces[map[ii + (j <<
          2)] - 1] + ps_index;
      }

      /* 'extract_cellface_edges:73' fe_idx = fe_idx + 1; */
      (*fe_idx)++;
    } else {
      /* 'extract_cellface_edges:74' else */
      /* 'extract_cellface_edges:75' count = count-1; */
      /* 'extract_cellface_edges:76' if (count ==1) */
      if (count + 1 == 1) {
        /* 'extract_cellface_edges:77' Face_elems(fe_idx,1) =Faces(ii)+ps_index; */
        Face_elems->data[*fe_idx - 1] = Faces[ii] + ps_index;

        /* 'extract_cellface_edges:78' Face_elems(fe_idx,2) = P(1); */
        Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = P[0];

        /* 'extract_cellface_edges:79' fe_idx = fe_idx + 1; */
        (*fe_idx)++;

        /* 'extract_cellface_edges:80' Face_elems(fe_idx,1) = P(1); */
        Face_elems->data[*fe_idx - 1] = P[0];

        /* 'extract_cellface_edges:81' Face_elems(fe_idx,2) = Faces(map(ii,2))+ps_index; */
        Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = Faces[map[4 + ii]
          - 1] + ps_index;

        /* 'extract_cellface_edges:82' fe_idx = fe_idx + 1; */
        (*fe_idx)++;
      } else {
        /* 'extract_cellface_edges:83' else */
        /* 'extract_cellface_edges:84' dist = zeros(count,1); */
        j = dist->size[0];
        dist->size[0] = count + 1;
        emxEnsureCapacity((emxArray__common *)dist, j, (int32_T)sizeof(real_T));
        for (j = 0; j <= count; j++) {
          dist->data[j] = 0.0;
        }

        /* 'extract_cellface_edges:85' C = ps_local(Faces(ii)+ps_index,:); */
        /* 'extract_cellface_edges:86' for jj=1:count */
        for (i = 0; i <= count; i++) {
          /* 'extract_cellface_edges:87' dist(jj) = sqrt((ps_local(P(jj),1)-C(1))^2+(ps_local(P(jj),2)-C(2))^2+(ps_local(P(jj),3)-C(3))^2); */
          a = ps_local->data[P[i] - 1] - ps_local->data[(Faces[ii] + ps_index) -
            1];
          b_a = ps_local->data[(P[i] + ps_local->size[0]) - 1] - ps_local->data
            [((Faces[ii] + ps_index) + ps_local->size[0]) - 1];
          c_a = ps_local->data[(P[i] + (ps_local->size[0] << 1)) - 1] -
            ps_local->data[((Faces[ii] + ps_index) + (ps_local->size[0] << 1)) -
            1];
          y = pow(a, 2.0);
          b_a = pow(b_a, 2.0);
          a = pow(c_a, 2.0);
          dist->data[i] = sqrt((y + b_a) + a);
        }

        /* 'extract_cellface_edges:89' [~,id] = sort(dist,'ascend'); */
        d_a = dist->size[0] - 1;
        j = vwork->size[0];
        vwork->size[0] = (int32_T)(int8_T)(d_a + 1);
        emxEnsureCapacity((emxArray__common *)vwork, j, (int32_T)sizeof(real_T));
        unnamed_idx_0 = (int8_T)dist->size[0];
        j = idx->size[0];
        idx->size[0] = (int32_T)unnamed_idx_0;
        emxEnsureCapacity((emxArray__common *)idx, j, (int32_T)sizeof(int32_T));
        i1 = 0;
        j = 1;
        while (j <= 1) {
          i1++;
          i = i1;
          for (k = 0; k <= d_a; k++) {
            vwork->data[k] = dist->data[i - 1];
            i++;
          }

          n = vwork->size[0];
          unnamed_idx_0 = (int8_T)vwork->size[0];
          j = iidx->size[0];
          iidx->size[0] = (int32_T)unnamed_idx_0;
          emxEnsureCapacity((emxArray__common *)iidx, j, (int32_T)sizeof(int32_T));
          np1 = n + 1;
          for (k = 1; k <= n; k++) {
            iidx->data[k - 1] = k;
          }

          for (k = 1; k <= n - 1; k += 2) {
            if (vwork->data[k - 1] <= vwork->data[k]) {
              p = TRUE;
            } else {
              p = FALSE;
            }

            if (p) {
            } else {
              iidx->data[k - 1] = k + 1;
              iidx->data[k] = k;
            }
          }

          j = idx0->size[0];
          idx0->size[0] = n;
          emxEnsureCapacity((emxArray__common *)idx0, j, (int32_T)sizeof(int32_T));
          i = n - 1;
          for (j = 0; j <= i; j++) {
            idx0->data[j] = 1;
          }

          i = 2;
          while (i < n) {
            i2 = i << 1;
            j = 1;
            for (pEnd = 1 + i; pEnd < np1; pEnd = qEnd + i) {
              b_p = j;
              q = pEnd;
              qEnd = j + i2;
              if (qEnd > np1) {
                qEnd = np1;
              }

              k = 0;
              kEnd = qEnd - j;
              while (k + 1 <= kEnd) {
                if (vwork->data[iidx->data[b_p - 1] - 1] <= vwork->data
                    [iidx->data[q - 1] - 1]) {
                  p = TRUE;
                } else {
                  p = FALSE;
                }

                if (p) {
                  idx0->data[k] = iidx->data[b_p - 1];
                  b_p++;
                  if (b_p == pEnd) {
                    while (q < qEnd) {
                      k++;
                      idx0->data[k] = iidx->data[q - 1];
                      q++;
                    }
                  }
                } else {
                  idx0->data[k] = iidx->data[q - 1];
                  q++;
                  if (q == qEnd) {
                    while (b_p < pEnd) {
                      k++;
                      idx0->data[k] = iidx->data[b_p - 1];
                      b_p++;
                    }
                  }
                }

                k++;
              }

              for (k = 0; k + 1 <= kEnd; k++) {
                iidx->data[(j + k) - 1] = idx0->data[k];
              }

              j = qEnd;
            }

            i = i2;
          }

          i = i1;
          for (k = 0; k <= d_a; k++) {
            idx->data[i - 1] = iidx->data[k];
            i++;
          }

          j = 2;
        }

        j = dist->size[0];
        dist->size[0] = idx->size[0];
        emxEnsureCapacity((emxArray__common *)dist, j, (int32_T)sizeof(real_T));
        i = idx->size[0] - 1;
        for (j = 0; j <= i; j++) {
          dist->data[j] = (real_T)idx->data[j];
        }

        /* 'extract_cellface_edges:90' Face_elems(fe_idx,1) =Faces(ii)+ps_index; */
        Face_elems->data[*fe_idx - 1] = Faces[ii] + ps_index;

        /* 'extract_cellface_edges:91' Face_elems(fe_idx,2) = P(id(1)); */
        Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = P[(int32_T)
          dist->data[0] - 1];

        /* 'extract_cellface_edges:92' fe_idx = fe_idx + 1; */
        (*fe_idx)++;

        /* 'extract_cellface_edges:93' for jj=1:count-1 */
        for (i = 0; i <= count - 1; i++) {
          /* 'extract_cellface_edges:94' Face_elems(fe_idx,1) = P(id(jj)); */
          Face_elems->data[*fe_idx - 1] = P[(int32_T)dist->data[i] - 1];

          /* 'extract_cellface_edges:95' Face_elems(fe_idx,2) = P(id(jj+1)); */
          Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = P[(int32_T)
            dist->data[i + 1] - 1];

          /* 'extract_cellface_edges:96' fe_idx = fe_idx + 1; */
          (*fe_idx)++;
        }

        /* 'extract_cellface_edges:98' Face_elems(fe_idx,1) = P(id(count)); */
        Face_elems->data[*fe_idx - 1] = P[(int32_T)dist->data[count] - 1];

        /* 'extract_cellface_edges:99' Face_elems(fe_idx,2) = Faces(map(ii,2))+ps_index; */
        Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = Faces[map[4 + ii]
          - 1] + ps_index;

        /* 'extract_cellface_edges:100' fe_idx = fe_idx + 1; */
        (*fe_idx)++;
      }
    }
  }

  emxFree_int32_T(&idx0);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
  emxFree_int32_T(&idx);
  emxFree_real_T(&dist);
}

/*
 * function [Face_elems,fe_idx] = get_celledgelems_noteight(Faces,FaceEdgeMap,cor_map,index,Face_ps,Face_elems,fe_idx,ps_index,ps_local,ps_localID)
 */
static void get_celledgelems_noteight(const int32_T Faces[4], const int32_T
  FaceEdgeMap[4], const emxArray_int32_T *cor_map, int32_T b_index, const
  emxArray_int32_T *Face_ps, emxArray_int32_T *Face_elems, int32_T *fe_idx,
  int32_T ps_index, const emxArray_real_T *ps_local, const emxArray_int32_T
  *ps_localID)
{
  int32_T cor_idx[4];
  int32_T i2;
  int32_T ii;
  int32_T ix;
  int32_T i;
  boolean_T exitg2;
  real_T a;
  boolean_T exitg1;
  emxArray_real_T *dist;
  emxArray_int32_T *idx;
  emxArray_real_T *vwork;
  emxArray_int32_T *iidx;
  emxArray_int32_T *idx0;
  int32_T P[10];
  int32_T count;
  static const int8_T map[8] = { 1, 2, 3, 4, 2, 3, 4, 1 };

  real_T b_a;
  real_T c_a;
  real_T y;
  int32_T d_a;
  int8_T unnamed_idx_0;
  int32_T i1;
  int32_T k;
  int32_T n;
  int32_T np1;
  boolean_T p;
  int32_T pEnd;
  int32_T b_p;
  int32_T q;
  int32_T qEnd;
  int32_T kEnd;

  /* 'extract_cellface_edges:107' map = [1 2; 2 3; 3 4; 4 1]; */
  /* 'extract_cellface_edges:108' cor_idx =int32([0 0 0 0]); */
  for (i2 = 0; i2 < 4; i2++) {
    cor_idx[i2] = 0;
  }

  /* 'extract_cellface_edges:109' for ii=1:4 */
  for (ii = 0; ii < 4; ii++) {
    /* 'extract_cellface_edges:110' flag = 0; */
    ix = 0;

    /* 'extract_cellface_edges:111' for jj=1:size(cor_map,1) */
    i = 0;
    exitg2 = 0U;
    while ((exitg2 == 0U) && (i <= cor_map->size[0] - 1)) {
      /* 'extract_cellface_edges:112' if (Faces(ii)==cor_map(jj)) */
      if (Faces[ii] == cor_map->data[i]) {
        /* 'extract_cellface_edges:113' cor_idx(ii) = jj + ps_index; */
        a = (1.0 + (real_T)i) + (real_T)ps_index;
        a = a < 0.0 ? ceil(a - 0.5) : floor(a + 0.5);
        cor_idx[ii] = (int32_T)a;

        /* 'extract_cellface_edges:114' flag = 1; */
        ix = 1;
        exitg2 = 1U;
      } else {
        i++;
      }
    }

    /* 'extract_cellface_edges:118' if (flag==0) */
    if (ix == 0) {
      /* 'extract_cellface_edges:119' for jj=1:ps_index */
      i = 1;
      exitg1 = 0U;
      while ((exitg1 == 0U) && (i <= ps_index)) {
        /* 'extract_cellface_edges:120' if (ps_localID(jj,6)==Faces(ii)) */
        if (ps_localID->data[(i + ps_localID->size[0] * 5) - 1] == Faces[ii]) {
          /* 'extract_cellface_edges:121' cor_idx(ii) = jj; */
          cor_idx[ii] = i;
          exitg1 = 1U;
        } else {
          i++;
        }
      }
    }
  }

  /* 'extract_cellface_edges:128' for ii=1:4 */
  b_emxInit_real_T(&dist, 1);
  b_emxInit_int32_T(&idx, 1);
  b_emxInit_real_T(&vwork, 1);
  b_emxInit_int32_T(&iidx, 1);
  b_emxInit_int32_T(&idx0, 1);
  for (ii = 0; ii < 4; ii++) {
    /* 'extract_cellface_edges:129' P = zeros(10,1); */
    for (i = 0; i < 10; i++) {
      P[i] = 0;
    }

    /* 'extract_cellface_edges:129' count = 1; */
    count = -1;

    /* 'extract_cellface_edges:130' for jj=1:index */
    for (i = 0; i + 1 <= b_index; i++) {
      /* 'extract_cellface_edges:131' if ((ps_localID(Face_ps(jj),5) == FaceEdgeMap(ii))) */
      if (ps_localID->data[(Face_ps->data[i] + (ps_localID->size[0] << 2)) - 1] ==
          FaceEdgeMap[ii]) {
        /* 'extract_cellface_edges:132' P(count) = Face_ps(jj); */
        P[count + 1] = Face_ps->data[i];

        /* 'extract_cellface_edges:133' count = count + 1; */
        count++;
      }
    }

    /* 'extract_cellface_edges:136' if (count == 1) */
    if (count + 2 == 1) {
      /* 'extract_cellface_edges:137' Face_elems(fe_idx,:) = cor_idx(map(ii,:)); */
      ix = *fe_idx;
      for (i2 = 0; i2 < 2; i2++) {
        Face_elems->data[(ix + Face_elems->size[0] * i2) - 1] = cor_idx[map[ii +
          (i2 << 2)] - 1];
      }

      /* 'extract_cellface_edges:138' fe_idx = fe_idx + 1; */
      (*fe_idx)++;
    } else {
      /* 'extract_cellface_edges:139' else */
      /* 'extract_cellface_edges:140' count = count-1; */
      /* 'extract_cellface_edges:141' if (count ==1) */
      if (count + 1 == 1) {
        /* 'extract_cellface_edges:142' Face_elems(fe_idx,1) =cor_idx(ii); */
        Face_elems->data[*fe_idx - 1] = cor_idx[ii];

        /* 'extract_cellface_edges:143' Face_elems(fe_idx,2) = P(1); */
        Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = P[0];

        /* 'extract_cellface_edges:144' fe_idx = fe_idx + 1; */
        (*fe_idx)++;

        /* 'extract_cellface_edges:145' Face_elems(fe_idx,1) = P(1); */
        Face_elems->data[*fe_idx - 1] = P[0];

        /* 'extract_cellface_edges:146' Face_elems(fe_idx,2) = cor_idx(map(ii,2)); */
        Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = cor_idx[map[4 +
          ii] - 1];

        /* 'extract_cellface_edges:147' fe_idx = fe_idx + 1; */
        (*fe_idx)++;
      } else {
        /* 'extract_cellface_edges:148' else */
        /* 'extract_cellface_edges:149' dist = zeros(count,1); */
        i2 = dist->size[0];
        dist->size[0] = count + 1;
        emxEnsureCapacity((emxArray__common *)dist, i2, (int32_T)sizeof(real_T));
        for (i2 = 0; i2 <= count; i2++) {
          dist->data[i2] = 0.0;
        }

        /* 'extract_cellface_edges:150' C = ps_local(cor_idx(ii),:); */
        /* 'extract_cellface_edges:151' for jj=1:count */
        for (i = 0; i <= count; i++) {
          /* 'extract_cellface_edges:152' dist(jj) = sqrt((ps_local(P(jj),1)-C(1))^2+(ps_local(P(jj),2)-C(2))^2+(ps_local(P(jj),3)-C(3))^2); */
          a = ps_local->data[P[i] - 1] - ps_local->data[cor_idx[ii] - 1];
          b_a = ps_local->data[(P[i] + ps_local->size[0]) - 1] - ps_local->data
            [(cor_idx[ii] + ps_local->size[0]) - 1];
          c_a = ps_local->data[(P[i] + (ps_local->size[0] << 1)) - 1] -
            ps_local->data[(cor_idx[ii] + (ps_local->size[0] << 1)) - 1];
          y = pow(a, 2.0);
          b_a = pow(b_a, 2.0);
          a = pow(c_a, 2.0);
          dist->data[i] = sqrt((y + b_a) + a);
        }

        /* 'extract_cellface_edges:154' [~,id] = sort(dist,'ascend'); */
        d_a = dist->size[0] - 1;
        i2 = vwork->size[0];
        vwork->size[0] = (int32_T)(int8_T)(d_a + 1);
        emxEnsureCapacity((emxArray__common *)vwork, i2, (int32_T)sizeof(real_T));
        unnamed_idx_0 = (int8_T)dist->size[0];
        i2 = idx->size[0];
        idx->size[0] = (int32_T)unnamed_idx_0;
        emxEnsureCapacity((emxArray__common *)idx, i2, (int32_T)sizeof(int32_T));
        i1 = 0;
        ix = 1;
        while (ix <= 1) {
          i1++;
          ix = i1;
          for (k = 0; k <= d_a; k++) {
            vwork->data[k] = dist->data[ix - 1];
            ix++;
          }

          n = vwork->size[0];
          unnamed_idx_0 = (int8_T)vwork->size[0];
          i2 = iidx->size[0];
          iidx->size[0] = (int32_T)unnamed_idx_0;
          emxEnsureCapacity((emxArray__common *)iidx, i2, (int32_T)sizeof
                            (int32_T));
          np1 = n + 1;
          for (k = 1; k <= n; k++) {
            iidx->data[k - 1] = k;
          }

          for (k = 1; k <= n - 1; k += 2) {
            if (vwork->data[k - 1] <= vwork->data[k]) {
              p = TRUE;
            } else {
              p = FALSE;
            }

            if (p) {
            } else {
              iidx->data[k - 1] = k + 1;
              iidx->data[k] = k;
            }
          }

          i2 = idx0->size[0];
          idx0->size[0] = n;
          emxEnsureCapacity((emxArray__common *)idx0, i2, (int32_T)sizeof
                            (int32_T));
          i = n - 1;
          for (i2 = 0; i2 <= i; i2++) {
            idx0->data[i2] = 1;
          }

          i = 2;
          while (i < n) {
            i2 = i << 1;
            ix = 1;
            for (pEnd = 1 + i; pEnd < np1; pEnd = qEnd + i) {
              b_p = ix;
              q = pEnd;
              qEnd = ix + i2;
              if (qEnd > np1) {
                qEnd = np1;
              }

              k = 0;
              kEnd = qEnd - ix;
              while (k + 1 <= kEnd) {
                if (vwork->data[iidx->data[b_p - 1] - 1] <= vwork->data
                    [iidx->data[q - 1] - 1]) {
                  p = TRUE;
                } else {
                  p = FALSE;
                }

                if (p) {
                  idx0->data[k] = iidx->data[b_p - 1];
                  b_p++;
                  if (b_p == pEnd) {
                    while (q < qEnd) {
                      k++;
                      idx0->data[k] = iidx->data[q - 1];
                      q++;
                    }
                  }
                } else {
                  idx0->data[k] = iidx->data[q - 1];
                  q++;
                  if (q == qEnd) {
                    while (b_p < pEnd) {
                      k++;
                      idx0->data[k] = iidx->data[b_p - 1];
                      b_p++;
                    }
                  }
                }

                k++;
              }

              for (k = 0; k + 1 <= kEnd; k++) {
                iidx->data[(ix + k) - 1] = idx0->data[k];
              }

              ix = qEnd;
            }

            i = i2;
          }

          ix = i1;
          for (k = 0; k <= d_a; k++) {
            idx->data[ix - 1] = iidx->data[k];
            ix++;
          }

          ix = 2;
        }

        i2 = dist->size[0];
        dist->size[0] = idx->size[0];
        emxEnsureCapacity((emxArray__common *)dist, i2, (int32_T)sizeof(real_T));
        i = idx->size[0] - 1;
        for (i2 = 0; i2 <= i; i2++) {
          dist->data[i2] = (real_T)idx->data[i2];
        }

        /* 'extract_cellface_edges:155' Face_elems(fe_idx,1) =cor_idx(ii); */
        Face_elems->data[*fe_idx - 1] = cor_idx[ii];

        /* 'extract_cellface_edges:156' Face_elems(fe_idx,2) = P(id(1)); */
        Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = P[(int32_T)
          dist->data[0] - 1];

        /* 'extract_cellface_edges:157' fe_idx = fe_idx + 1; */
        (*fe_idx)++;

        /* 'extract_cellface_edges:158' for jj=1:count-1 */
        for (i = 0; i <= count - 1; i++) {
          /* 'extract_cellface_edges:159' Face_elems(fe_idx,1) = P(id(jj)); */
          Face_elems->data[*fe_idx - 1] = P[(int32_T)dist->data[i] - 1];

          /* 'extract_cellface_edges:160' Face_elems(fe_idx,2) = P(id(jj+1)); */
          Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = P[(int32_T)
            dist->data[i + 1] - 1];

          /* 'extract_cellface_edges:161' fe_idx = fe_idx + 1; */
          (*fe_idx)++;
        }

        /* 'extract_cellface_edges:163' Face_elems(fe_idx,1) = P(id(count)); */
        Face_elems->data[*fe_idx - 1] = P[(int32_T)dist->data[count] - 1];

        /* 'extract_cellface_edges:164' Face_elems(fe_idx,2) = cor_idx(map(ii,2)); */
        Face_elems->data[(*fe_idx + Face_elems->size[0]) - 1] = cor_idx[map[4 +
          ii] - 1];

        /* 'extract_cellface_edges:165' fe_idx = fe_idx + 1; */
        (*fe_idx)++;
      }
    }
  }

  emxFree_int32_T(&idx0);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
  emxFree_int32_T(&idx);
  emxFree_real_T(&dist);
}

/*
 * function [ps_local,ps_localID,ps_idx]= get_pntlist(cell_bnds,opphes,triID,trivertID,tri,flag_v,ps_local,ps_localID,ps_idx)
 */
static void get_pntlist(const real_T cell_bnds[6], const emxArray_int32_T
  *opphes, int32_T triID, const int32_T trivertID[3], const real_T tri[9], const
  int32_T flag_v[3], emxArray_real_T *ps_local, emxArray_int32_T *ps_localID,
  int32_T *ps_idx)
{
  int32_T count;
  int32_T ii;
  int32_T i11;
  int32_T jj;
  boolean_T exitg2;
  real_T id[6];
  real_T v[3];
  real_T facebnd[30];
  int32_T unusedU4[6];
  real_T unusedU3[6];
  real_T unusedU2[18];
  real_T unusedU1[24];
  boolean_T exitg1;
  boolean_T guard1 = FALSE;
  static const int8_T fbmap[18] = { 1, 2, 1, 2, 1, 1, 2, 3, 3, 3, 3, 2, 3, 1, 2,
    1, 2, 3 };

  real_T xi;
  real_T eta;
  real_T fcormap[4];
  real_T fedgmap[4];
  static const int8_T b_fcormap[24] = { 4, 1, 3, 4, 4, 8, 1, 2, 2, 3, 1, 5, 2, 6,
    6, 7, 5, 6, 3, 5, 7, 8, 8, 7 };

  static const int8_T b_fedgmap[24] = { 4, 1, 2, 3, 4, 12, 1, 6, 6, 7, 5, 9, 2,
    9, 10, 11, 12, 10, 3, 5, 7, 8, 8, 11 };

  /* 'get_pntlist:3' [ps_local,ps_localID,ps_idx]=compute_intersections(cell_bnds,opphes,trivertID,triID,tri,ps_local,ps_localID,ps_idx); */
  count = *ps_idx;
  compute_intersections(cell_bnds, opphes, trivertID, triID, tri, ps_local,
                        ps_localID, &count);
  *ps_idx = count;

  /* 'get_pntlist:4' for ii=1:3 */
  for (ii = 0; ii < 3; ii++) {
    /* 'get_pntlist:5' if (flag_v(ii)==7) */
    if (flag_v[ii] == 7) {
      /* 'get_pntlist:6' v = tri(ii,:); */
      /* 'get_pntlist:7' count = int32(0); */
      count = 0;

      /* 'get_pntlist:8' for jj=1:(ps_idx-1) */
      i11 = *ps_idx - 1;
      jj = 1;
      exitg2 = 0U;
      while ((exitg2 == 0U) && (jj <= i11)) {
        /* 'get_pntlist:9' if (trivertID(ii) == ps_localID(jj,1)) */
        if (trivertID[ii] == ps_localID->data[jj - 1]) {
          /* 'get_pntlist:10' count = int32(1); */
          count = 1;
          exitg2 = 1U;
        } else {
          jj++;
        }
      }

      /* 'get_pntlist:14' if (count ==0) */
      if (count == 0) {
        /* 'get_pntlist:15' id = [0 0 0 0 0 0]; */
        for (i11 = 0; i11 < 6; i11++) {
          id[i11] = 0.0;
        }

        /* 'get_pntlist:16' id(1)=trivertID(ii); */
        id[0] = (real_T)trivertID[ii];

        /* 'get_pntlist:17' [v,id]=assign_id(cell_bnds,v,id); */
        for (i11 = 0; i11 < 3; i11++) {
          v[i11] = tri[ii + 3 * i11];
        }

        /* 'get_pntlist:26' tol = 1e-6; */
        /* 'get_pntlist:27' [~,~,~,~,~,facebnd,fbmap,fcormap,fedgmap] = gridcell_data(cell_bnds); */
        b_gridcell_data(cell_bnds, unusedU1, unusedU2, unusedU3, unusedU4,
                        facebnd);

        /* 'get_pntlist:29' for jj=1:6 */
        jj = 0;
        exitg1 = 0U;
        while ((exitg1 == 0U) && (jj < 6)) {
          /* 'get_pntlist:30' cx = fbmap(jj,1); */
          /* 'get_pntlist:30' cxmin=facebnd(jj,1); */
          /* 'get_pntlist:30' cxmax = facebnd(jj,2); */
          /* 'get_pntlist:31' cy = fbmap(jj,2); */
          /* 'get_pntlist:31' cymin=facebnd(jj,3); */
          /* 'get_pntlist:31' cymax = facebnd(jj,4); */
          /* 'get_pntlist:32' cz = fbmap(jj,3); */
          /* 'get_pntlist:32' czval = facebnd(jj,5); */
          /* 'get_pntlist:33' eps = v(cz)-czval; */
          /* 'get_pntlist:34' if (eps ==0) */
          guard1 = FALSE;
          if (tri[ii + 3 * (fbmap[12 + jj] - 1)] - facebnd[24 + jj] == 0.0) {
            /* 'get_pntlist:35' xp = v(cx); */
            /* 'get_pntlist:35' yp = v(cy); */
            /* 'get_pntlist:36' xi = (xp-cxmin)/(cxmax-cxmin); */
            xi = (tri[ii + 3 * (fbmap[jj] - 1)] - facebnd[jj]) / (facebnd[6 + jj]
              - facebnd[jj]);

            /* 'get_pntlist:37' eta = (yp-cymin)/(cymax-cymin); */
            eta = (tri[ii + 3 * (fbmap[6 + jj] - 1)] - facebnd[12 + jj]) /
              (facebnd[18 + jj] - facebnd[12 + jj]);

            /* 'get_pntlist:39' if (abs(xi)<=tol) */
            if (fabs(xi) <= 1.0E-6) {
              /* 'get_pntlist:40' xi =0 ; */
              xi = 0.0;
            } else {
              if (fabs(xi - 1.0) <= 1.0E-6) {
                /* 'get_pntlist:41' elseif (abs(xi-1)<=tol) */
                /* 'get_pntlist:42' xi = 1; */
                xi = 1.0;
              }
            }

            /* 'get_pntlist:44' if (abs(eta)<=tol) */
            if (fabs(eta) <= 1.0E-6) {
              /* 'get_pntlist:45' eta =0 ; */
              eta = 0.0;
            } else {
              if (fabs(eta - 1.0) <= 1.0E-6) {
                /* 'get_pntlist:46' elseif (abs(eta-1)<=tol) */
                /* 'get_pntlist:47' eta = 1; */
                eta = 1.0;
              }
            }

            /* 'get_pntlist:50' if ((xi>=0) && (xi<=1) && (eta>=0) && (eta<=1)) */
            if ((xi >= 0.0) && (xi <= 1.0) && (eta >= 0.0) && (eta <= 1.0)) {
              /* 'get_pntlist:51' [v,id] = assign_id_basedon_xieta(xi,eta,xp,yp,cx,cxmin,cxmax,cy,cymin,cymax,jj,fcormap(jj,:),fedgmap(jj,:),v,id); */
              for (i11 = 0; i11 < 3; i11++) {
                v[i11] = tri[ii + 3 * i11];
              }

              for (i11 = 0; i11 < 4; i11++) {
                fcormap[i11] = (real_T)b_fcormap[jj + 6 * i11];
                fedgmap[i11] = (real_T)b_fedgmap[jj + 6 * i11];
              }

              assign_id_basedon_xieta(xi, eta, tri[ii + 3 * (fbmap[jj] - 1)],
                tri[ii + 3 * (fbmap[6 + jj] - 1)], (real_T)fbmap[jj], facebnd[jj],
                facebnd[6 + jj], (real_T)fbmap[6 + jj], facebnd[12 + jj],
                facebnd[18 + jj], 1.0 + (real_T)jj, fcormap, fedgmap, v, id);
              exitg1 = 1U;
            } else {
              guard1 = TRUE;
            }
          } else {
            guard1 = TRUE;
          }

          if (guard1 == TRUE) {
            jj++;
          }
        }

        /* 'get_pntlist:18' ps_local(ps_idx,:) = v; */
        count = *ps_idx;
        for (i11 = 0; i11 < 3; i11++) {
          ps_local->data[(count + ps_local->size[0] * i11) - 1] = v[i11];
        }

        /* 'get_pntlist:19' ps_localID(ps_idx,:)=id; */
        count = *ps_idx;
        for (i11 = 0; i11 < 6; i11++) {
          xi = id[i11];
          xi = xi < 0.0 ? ceil(xi - 0.5) : floor(xi + 0.5);
          ps_localID->data[(count + ps_localID->size[0] * i11) - 1] = (int32_T)
            xi;
        }

        /* 'get_pntlist:20' ps_idx = ps_idx + 1; */
        (*ps_idx)++;
      }
    }
  }
}

/*
 * function [flag_cell,tris_cell,tris_idx,comp_index_cell] = get_trisincell(cell_bnds,ps,num_tris,tris,comp_index)
 */
static void get_trisincell(const real_T cell_bnds[6], const emxArray_real_T *ps,
  int32_T num_tris, const emxArray_int32_T *tris, const emxArray_int32_T
  *comp_index, emxArray_int32_T *tris_cell, emxArray_int32_T *comp_index_cell,
  int32_T *flag_cell, int32_T *tris_idx)
{
  int32_T i0;
  int32_T flag;
  real_T tri_bnds[6];
  int32_T jj;
  real_T mtmp;
  boolean_T b0;
  int32_T b_comp_index_cell[2];
  int32_T b_comp_index[2];
  emxArray_int32_T c_comp_index_cell;
  emxArray_int32_T c_comp_index;
  emxArray_int32_T *b_tris_cell;
  int32_T i1;
  emxArray_int32_T *d_comp_index_cell;

  /* 'get_trisincell:3' coder.inline('never') */
  /* 'get_trisincell:4' tris_cell = zeros(num_tris,3,'int32'); */
  i0 = tris_cell->size[0] * tris_cell->size[1];
  tris_cell->size[0] = num_tris;
  tris_cell->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)tris_cell, i0, (int32_T)sizeof(int32_T));
  flag = num_tris * 3 - 1;
  for (i0 = 0; i0 <= flag; i0++) {
    tris_cell->data[i0] = 0;
  }

  /* 'get_trisincell:5' comp_index_cell = zeros(num_tris,1,'int32'); */
  i0 = comp_index_cell->size[0];
  comp_index_cell->size[0] = num_tris;
  emxEnsureCapacity((emxArray__common *)comp_index_cell, i0, (int32_T)sizeof
                    (int32_T));
  flag = num_tris - 1;
  for (i0 = 0; i0 <= flag; i0++) {
    comp_index_cell->data[i0] = 0;
  }

  /* 'get_trisincell:6' l = int32(1); */
  *tris_idx = 0;

  /* 'get_trisincell:7' flag_cell = int32(0); */
  *flag_cell = 0;

  /* 'get_trisincell:8' tri_bnds = zeros(1,6); */
  for (i0 = 0; i0 < 6; i0++) {
    tri_bnds[i0] = 0.0;
  }

  /* 'get_trisincell:10' for jj = 1:num_tris */
  for (jj = 0; jj + 1 <= num_tris; jj++) {
    /* 'get_trisincell:11' tri_bnds(1) = min(ps(tris(jj,:),1)); */
    mtmp = ps->data[tris->data[jj] - 1];
    for (flag = 0; flag < 2; flag++) {
      if (ps->data[tris->data[jj + tris->size[0] * (flag + 1)] - 1] < mtmp) {
        mtmp = ps->data[tris->data[jj + tris->size[0] * (flag + 1)] - 1];
      }
    }

    tri_bnds[0] = mtmp;

    /* 'get_trisincell:12' tri_bnds(2) = max(ps(tris(jj,:),1)); */
    mtmp = ps->data[tris->data[jj] - 1];
    for (flag = 0; flag < 2; flag++) {
      if (ps->data[tris->data[jj + tris->size[0] * (flag + 1)] - 1] > mtmp) {
        mtmp = ps->data[tris->data[jj + tris->size[0] * (flag + 1)] - 1];
      }
    }

    tri_bnds[1] = mtmp;

    /* 'get_trisincell:13' tri_bnds(3) = min(ps(tris(jj,:),2)); */
    mtmp = ps->data[(tris->data[jj] + ps->size[0]) - 1];
    for (flag = 0; flag < 2; flag++) {
      if (ps->data[(tris->data[jj + tris->size[0] * (flag + 1)] + ps->size[0]) -
          1] < mtmp) {
        mtmp = ps->data[(tris->data[jj + tris->size[0] * (flag + 1)] + ps->size
                         [0]) - 1];
      }
    }

    tri_bnds[2] = mtmp;

    /* 'get_trisincell:14' tri_bnds(4) = max(ps(tris(jj,:),2)); */
    mtmp = ps->data[(tris->data[jj] + ps->size[0]) - 1];
    for (flag = 0; flag < 2; flag++) {
      if (ps->data[(tris->data[jj + tris->size[0] * (flag + 1)] + ps->size[0]) -
          1] > mtmp) {
        mtmp = ps->data[(tris->data[jj + tris->size[0] * (flag + 1)] + ps->size
                         [0]) - 1];
      }
    }

    tri_bnds[3] = mtmp;

    /* 'get_trisincell:15' tri_bnds(5) = min(ps(tris(jj,:),3)); */
    mtmp = ps->data[(tris->data[jj] + (ps->size[0] << 1)) - 1];
    for (flag = 0; flag < 2; flag++) {
      if (ps->data[(tris->data[jj + tris->size[0] * (flag + 1)] + (ps->size[0] <<
            1)) - 1] < mtmp) {
        mtmp = ps->data[(tris->data[jj + tris->size[0] * (flag + 1)] + (ps->
          size[0] << 1)) - 1];
      }
    }

    tri_bnds[4] = mtmp;

    /* 'get_trisincell:16' tri_bnds(6) = max(ps(tris(jj,:),3)); */
    mtmp = ps->data[(tris->data[jj] + (ps->size[0] << 1)) - 1];
    for (flag = 0; flag < 2; flag++) {
      if (ps->data[(tris->data[jj + tris->size[0] * (flag + 1)] + (ps->size[0] <<
            1)) - 1] > mtmp) {
        mtmp = ps->data[(tris->data[jj + tris->size[0] * (flag + 1)] + (ps->
          size[0] << 1)) - 1];
      }
    }

    tri_bnds[5] = mtmp;

    /* 'get_trisincell:17' [flag]=bndbox_overlap(tri_bnds,cell_bnds); */
    /* 'bndbox_overlap:2' flag = int32(0); */
    flag = 0;

    /* 'bndbox_overlap:3' txmin = tri_bnds(1); */
    /* 'bndbox_overlap:4' txmax = tri_bnds(2); */
    /* 'bndbox_overlap:5' tymin = tri_bnds(3); */
    /* 'bndbox_overlap:6' tymax = tri_bnds(4); */
    /* 'bndbox_overlap:7' tzmin = tri_bnds(5); */
    /* 'bndbox_overlap:8' tzmax = tri_bnds(6); */
    /* 'bndbox_overlap:10' cxmin = cell_bnds(1); */
    /* 'bndbox_overlap:11' cxmax = cell_bnds(2); */
    /* 'bndbox_overlap:12' cymin = cell_bnds(3); */
    /* 'bndbox_overlap:13' cymax = cell_bnds(4); */
    /* 'bndbox_overlap:14' czmin = cell_bnds(5); */
    /* 'bndbox_overlap:15' czmax = cell_bnds(6); */
    /* 'bndbox_overlap:18' if ~((txmin >cxmax)||(txmax <cxmin)||(tymin >cymax)||(tymax <cymin)||(tzmin >czmax)||(tzmax <czmin)) */
    if ((tri_bnds[0] > cell_bnds[1]) || (tri_bnds[1] < cell_bnds[0]) ||
        (tri_bnds[2] > cell_bnds[3]) || (tri_bnds[3] < cell_bnds[2]) ||
        (tri_bnds[4] > cell_bnds[5]) || (tri_bnds[5] < cell_bnds[4])) {
      b0 = TRUE;
    } else {
      b0 = FALSE;
    }

    if (!b0) {
      /* 'bndbox_overlap:19' flag = int32(1); */
      flag = 1;
    }

    /* 'get_trisincell:18' if (flag ==1) */
    if (flag == 1) {
      /* 'get_trisincell:19' tris_cell(l,:) = tris(jj,:); */
      for (i0 = 0; i0 < 3; i0++) {
        tris_cell->data[*tris_idx + tris_cell->size[0] * i0] = tris->data[jj +
          tris->size[0] * i0];
      }

      /* 'get_trisincell:20' comp_index_cell(l,:) = comp_index(jj,:); */
      b_comp_index_cell[0] = comp_index_cell->size[0];
      b_comp_index_cell[1] = 1;
      b_comp_index[0] = comp_index->size[0];
      b_comp_index[1] = 1;
      c_comp_index_cell = *comp_index_cell;
      c_comp_index_cell.size = (int32_T *)&b_comp_index_cell;
      c_comp_index_cell.numDimensions = 1;
      c_comp_index = *comp_index;
      c_comp_index.size = (int32_T *)&b_comp_index;
      c_comp_index.numDimensions = 1;
      c_comp_index_cell.data[*tris_idx] = c_comp_index.data[jj];

      /* fprintf('Tris_cell %d is Tris Global %d\n',l,jj); */
      /* 'get_trisincell:22' l = l+1; */
      (*tris_idx)++;

      /* 'get_trisincell:23' flag_cell = int32(1); */
      *flag_cell = 1;
    }
  }

  /* 'get_trisincell:27' tris_idx = l-1; */
  /* 'get_trisincell:28' tris_cell = tris_cell(1:tris_idx,:); */
  if (1 > *tris_idx) {
    i0 = 0;
  } else {
    i0 = *tris_idx;
  }

  emxInit_int32_T(&b_tris_cell, 2);
  jj = b_tris_cell->size[0] * b_tris_cell->size[1];
  b_tris_cell->size[0] = i0;
  b_tris_cell->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)b_tris_cell, jj, (int32_T)sizeof(int32_T));
  for (jj = 0; jj < 3; jj++) {
    flag = i0 - 1;
    for (i1 = 0; i1 <= flag; i1++) {
      b_tris_cell->data[i1 + b_tris_cell->size[0] * jj] = tris_cell->data[i1 +
        tris_cell->size[0] * jj];
    }
  }

  i0 = tris_cell->size[0] * tris_cell->size[1];
  tris_cell->size[0] = b_tris_cell->size[0];
  tris_cell->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)tris_cell, i0, (int32_T)sizeof(int32_T));
  for (i0 = 0; i0 < 3; i0++) {
    flag = b_tris_cell->size[0] - 1;
    for (jj = 0; jj <= flag; jj++) {
      tris_cell->data[jj + tris_cell->size[0] * i0] = b_tris_cell->data[jj +
        b_tris_cell->size[0] * i0];
    }
  }

  emxFree_int32_T(&b_tris_cell);

  /* 'get_trisincell:29' comp_index_cell = comp_index_cell(1:tris_idx); */
  if (1 > *tris_idx) {
    i0 = 0;
  } else {
    i0 = *tris_idx;
  }

  b_emxInit_int32_T(&d_comp_index_cell, 1);
  jj = d_comp_index_cell->size[0];
  d_comp_index_cell->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)d_comp_index_cell, jj, (int32_T)sizeof
                    (int32_T));
  flag = i0 - 1;
  for (i0 = 0; i0 <= flag; i0++) {
    d_comp_index_cell->data[i0] = comp_index_cell->data[i0];
  }

  i0 = comp_index_cell->size[0];
  comp_index_cell->size[0] = d_comp_index_cell->size[0];
  emxEnsureCapacity((emxArray__common *)comp_index_cell, i0, (int32_T)sizeof
                    (int32_T));
  flag = d_comp_index_cell->size[0] - 1;
  for (i0 = 0; i0 <= flag; i0++) {
    comp_index_cell->data[i0] = d_comp_index_cell->data[i0];
  }

  emxFree_int32_T(&d_comp_index_cell);
}

/*
 * function [faces,corners,nrm,face_equ,face_equID,facebnd,fbmap,fcormap,fedgmap] = gridcell_data(cell_bnds)
 *  Bounds of the cell
 */
static void gridcell_data(const real_T cell_bnds[6], real_T corners[24], real_T
  nrm[18], real_T face_equ[6], int32_T face_equID[6])
{
  int32_T i;
  static const int8_T iv0[3] = { 0, 0, -1 };

  static const int8_T iv1[3] = { 1, 0, 0 };

  static const int8_T iv2[3] = { 0, 1, 0 };

  static const int8_T iv3[3] = { -1, 0, 0 };

  static const int8_T iv4[3] = { 0, -1, 0 };

  static const int8_T iv5[3] = { 0, 0, 1 };

  /* 'gridcell_data:4' xmin = cell_bnds(1); */
  /* 'gridcell_data:5' xmax = cell_bnds(2); */
  /* 'gridcell_data:6' ymin = cell_bnds(3); */
  /* 'gridcell_data:7' ymax = cell_bnds(4); */
  /* 'gridcell_data:8' zmin = cell_bnds(5); */
  /* 'gridcell_data:9' zmax = cell_bnds(6); */
  /*  Face-Corner Map */
  /* 'gridcell_data:12' faces = int32([1 4 3 2;1 2 6 5;2 3 7 6;3 4 8 7;1 5 8 4;5 6 7 8]); */
  /*  Corners of the cell */
  /* 'gridcell_data:15' corners = zeros(8,3); */
  for (i = 0; i < 24; i++) {
    corners[i] = 0.0;
  }

  /* 'gridcell_data:16' corners(1,:) = [xmax ymin zmin]; */
  corners[0] = cell_bnds[1];
  corners[8] = cell_bnds[2];
  corners[16] = cell_bnds[4];

  /* 'gridcell_data:17' corners(2,:) = [xmax ymax zmin]; */
  corners[1] = cell_bnds[1];
  corners[9] = cell_bnds[3];
  corners[17] = cell_bnds[4];

  /* 'gridcell_data:18' corners(3,:) = [xmin ymax zmin]; */
  corners[2] = cell_bnds[0];
  corners[10] = cell_bnds[3];
  corners[18] = cell_bnds[4];

  /* 'gridcell_data:19' corners(4,:) = [xmin ymin zmin]; */
  corners[3] = cell_bnds[0];
  corners[11] = cell_bnds[2];
  corners[19] = cell_bnds[4];

  /* 'gridcell_data:20' corners(5,:) = [xmax ymin zmax]; */
  corners[4] = cell_bnds[1];
  corners[12] = cell_bnds[2];
  corners[20] = cell_bnds[5];

  /* 'gridcell_data:21' corners(6,:) = [xmax ymax zmax]; */
  corners[5] = cell_bnds[1];
  corners[13] = cell_bnds[3];
  corners[21] = cell_bnds[5];

  /* 'gridcell_data:22' corners(7,:) = [xmin ymax zmax]; */
  corners[6] = cell_bnds[0];
  corners[14] = cell_bnds[3];
  corners[22] = cell_bnds[5];

  /* 'gridcell_data:23' corners(8,:) = [xmin ymin zmax]; */
  corners[7] = cell_bnds[0];
  corners[15] = cell_bnds[2];
  corners[23] = cell_bnds[5];

  /*  Normals to the cell faces */
  /* 'gridcell_data:26' nrm = zeros(6,3); */
  for (i = 0; i < 18; i++) {
    nrm[i] = 0.0;
  }

  /* 'gridcell_data:27' nrm(1,:) = [0 0 -1]; */
  for (i = 0; i < 3; i++) {
    nrm[6 * i] = (real_T)iv0[i];

    /* 'gridcell_data:28' nrm(2,:) = [1 0 0]; */
    nrm[1 + 6 * i] = (real_T)iv1[i];

    /* 'gridcell_data:29' nrm(3,:) = [0 1 0]; */
    nrm[2 + 6 * i] = (real_T)iv2[i];

    /* 'gridcell_data:30' nrm(4,:) = [-1 0 0]; */
    nrm[3 + 6 * i] = (real_T)iv3[i];

    /* 'gridcell_data:31' nrm(5,:) = [0 -1 0]; */
    nrm[4 + 6 * i] = (real_T)iv4[i];

    /* 'gridcell_data:32' nrm(6,:) = [0 0 1]; */
    nrm[5 + 6 * i] = (real_T)iv5[i];
  }

  /*  Equation of the cell faces */
  /* 'gridcell_data:35' face_equ = zeros(6,1); */
  /* 'gridcell_data:36' face_equ(1) = zmin; */
  for (i = 0; i < 6; i++) {
    face_equ[i] = 0.0;

    /*  ID of the face to distinguish if it is the x/y/z-plane */
    /* 'gridcell_data:44' face_equID = zeros(6,1,'int32'); */
    face_equID[i] = 0;
  }

  face_equ[0] = cell_bnds[4];

  /* 'gridcell_data:37' face_equ(2) = xmax; */
  face_equ[1] = cell_bnds[1];

  /* 'gridcell_data:38' face_equ(3) = ymax; */
  face_equ[2] = cell_bnds[3];

  /* 'gridcell_data:39' face_equ(4) = xmin; */
  face_equ[3] = cell_bnds[0];

  /* 'gridcell_data:40' face_equ(5) = ymin; */
  face_equ[4] = cell_bnds[2];

  /* 'gridcell_data:41' face_equ(6) = zmax; */
  face_equ[5] = cell_bnds[5];

  /* 'gridcell_data:45' face_equID(1) = 3; */
  face_equID[0] = 3;

  /* 'gridcell_data:46' face_equID(2) = 1; */
  face_equID[1] = 1;

  /* 'gridcell_data:47' face_equID(3) = 2; */
  face_equID[2] = 2;

  /* 'gridcell_data:48' face_equID(4) = 1; */
  face_equID[3] = 1;

  /* 'gridcell_data:49' face_equID(5) = 2; */
  face_equID[4] = 2;

  /* 'gridcell_data:50' face_equID(6) = 3; */
  face_equID[5] = 3;

  /*  Face bounds and map with coordinate for each cell face. */
  /*  For a typical face, */
  /*  fbmap(1),fbmap(2) : The coordinates that vary for the face, their bounds are */
  /*  given by facebnd(1:2),facebnd(3:4), respectively. fbmap(3) is the */
  /*  coordinate that remains constant and this constant value is facebnd(5). */
  /* 'gridcell_data:58' facebnd = [xmin xmax ymin ymax zmin; */
  /* 'gridcell_data:59'     ymin ymax zmin zmax xmax; */
  /* 'gridcell_data:60'     xmin xmax zmin zmax ymax; */
  /* 'gridcell_data:61'     ymin ymax zmin zmax xmin; */
  /* 'gridcell_data:62'     xmin xmax zmin zmax ymin; */
  /* 'gridcell_data:63'     xmin xmax ymin ymax zmax]; */
  /* 'gridcell_data:65' fbmap =[1 2 3; 2 3 1; 1 3 2; 2 3 1;1 3 2;1 2 3]; */
  /*  Mapping of the barycentric coordinates with the corners and cell edges. */
  /*  Corners: fcormap(1:4) corresponds to (xi,eta) = (0,0),(1,0),(1,1),(0,1) */
  /*  Edges : fedgmap(1:4) corresponds to(xi,eta) = ((0,1),0), (1,(0,1)),((0,1),1), (0,(0,1)) */
  /* 'gridcell_data:70' fcormap =[4 1 2 3; 1 2 6 5; 3 2 6 7; 4 3 7 8; 4 1 5 8; 8 5 6 7]; */
  /* 'gridcell_data:71' fedgmap = [4 1 2 3; 1 6 9 5; 2 6 10 7; 3 7 11 8; 4 5 12 8; 12 9 10 11]; */
}

/*
 * function [ps_idx,ps_local,ps_localID] =intersection_triedge_cellface_2(cell_bnds,triID,edgID,edgvertID,edgps,ps_idx,ps_local,ps_localID)
 *  This function computes the intersection of a triangle edge with edgID with a cell
 *  face. Uses exact arithmetic to find intersections.
 */
static void intersection_triedge_cellface_2(const real_T cell_bnds[6], int32_T
  triID, real_T edgID, const int32_T edgvertID[2], const real_T edgps[6],
  int32_T *ps_idx, emxArray_real_T *ps_local, emxArray_int32_T *ps_localID)
{
  real_T facebnd[30];
  int32_T unusedU4[6];
  real_T unusedU3[6];
  real_T unusedU2[18];
  real_T unusedU1[24];
  real_T intpoints[15];
  int32_T i13;
  int32_T idx;
  int32_T intpoints_id[30];
  int32_T jj;
  real_T counter;
  real_T yp;
  static const int8_T fbmap[18] = { 1, 2, 1, 2, 1, 1, 2, 3, 3, 3, 3, 2, 3, 1, 2,
    1, 2, 3 };

  real_T alpha;
  real_T xi;
  real_T eta;
  int8_T b_fbmap[3];
  static const int8_T fcormap[24] = { 4, 1, 3, 4, 4, 8, 1, 2, 2, 3, 1, 5, 2, 6,
    6, 7, 5, 6, 3, 5, 7, 8, 8, 7 };

  real_T facecorID[4];
  static const int8_T fedgmap[24] = { 4, 1, 2, 3, 4, 12, 1, 6, 6, 7, 5, 9, 2, 9,
    10, 11, 12, 10, 3, 5, 7, 8, 8, 11 };

  real_T faceedgID[4];
  real_T id[6];
  real_T p[3];
  int32_T count;
  int32_T b_ps_idx;
  boolean_T exitg2;
  emxArray_real_T *b_intpoints;
  emxArray_int32_T *b_intpoints_id;
  emxArray_int32_T *samefacepnts;
  emxArray_int32_T *b_samefacepnts;
  emxArray_int32_T *c_samefacepnts;
  int32_T ii;
  boolean_T guard1 = FALSE;
  boolean_T exitg1;

  /* 'intersection_triedge_cellface_2:6' [~,~,~,~,~,facebnd,fbmap,fcormap,fedgmap] = gridcell_data(cell_bnds); */
  b_gridcell_data(cell_bnds, unusedU1, unusedU2, unusedU3, unusedU4, facebnd);

  /* 'intersection_triedge_cellface_2:8' tol = 1e-6; */
  /* 'intersection_triedge_cellface_2:9' v1 = edgps(1,:); */
  /* 'intersection_triedge_cellface_2:10' v2 = edgps(2,:); */
  /* 'intersection_triedge_cellface_2:11' intpoints = zeros(5,3); */
  for (i13 = 0; i13 < 15; i13++) {
    intpoints[i13] = 0.0;
  }

  /* 'intersection_triedge_cellface_2:11' idx = int32(1); */
  idx = 0;

  /* 'intersection_triedge_cellface_2:12' intpoints_id = zeros(5,6,'int32'); */
  memset((void *)&intpoints_id[0], 0, 30U * sizeof(int32_T));

  /* 'intersection_triedge_cellface_2:13' for jj=1:6 */
  for (jj = 0; jj < 6; jj++) {
    /* 'intersection_triedge_cellface_2:14' cx = fbmap(jj,1); */
    /* 'intersection_triedge_cellface_2:14' cxmin=facebnd(jj,1); */
    /* 'intersection_triedge_cellface_2:14' cxmax = facebnd(jj,2); */
    /* 'intersection_triedge_cellface_2:15' cy = fbmap(jj,2); */
    /* 'intersection_triedge_cellface_2:15' cymin=facebnd(jj,3); */
    /* 'intersection_triedge_cellface_2:15' cymax = facebnd(jj,4); */
    /* 'intersection_triedge_cellface_2:16' cz = fbmap(jj,3); */
    /* 'intersection_triedge_cellface_2:16' czval = facebnd(jj,5); */
    /* 'intersection_triedge_cellface_2:17' xp = 0; */
    counter = 0.0;

    /* 'intersection_triedge_cellface_2:17' yp = 0; */
    yp = 0.0;

    /* 'intersection_triedge_cellface_2:18' eps = (v2(cz)-v1(cz)); */
    /* 'intersection_triedge_cellface_2:20' if (eps ~= 0) */
    if (edgps[1 + ((fbmap[12 + jj] - 1) << 1)] - edgps[(fbmap[12 + jj] - 1) << 1]
        != 0.0) {
      /* 'intersection_triedge_cellface_2:21' alpha = (czval-v1(cz))/(v2(cz)-v1(cz)); */
      alpha = (facebnd[24 + jj] - edgps[(fbmap[12 + jj] - 1) << 1]) / (edgps[1 +
        ((fbmap[12 + jj] - 1) << 1)] - edgps[(fbmap[12 + jj] - 1) << 1]);

      /* 'intersection_triedge_cellface_2:22' if (abs(alpha)<=tol) */
      if (fabs(alpha) <= 1.0E-6) {
        /* 'intersection_triedge_cellface_2:23' alpha =0 ; */
        alpha = 0.0;
      } else {
        if (fabs(alpha - 1.0) <= 1.0E-6) {
          /* 'intersection_triedge_cellface_2:24' elseif (abs(alpha-1)<=tol) */
          /* 'intersection_triedge_cellface_2:25' alpha = 1; */
          alpha = 1.0;
        }
      }

      /* 'intersection_triedge_cellface_2:27' if (alpha>=0 && alpha<=1) */
      if ((alpha >= 0.0) && (alpha <= 1.0)) {
        /* 'intersection_triedge_cellface_2:28' if (alpha ==0) */
        if (alpha == 0.0) {
          /* 'intersection_triedge_cellface_2:29' xp = v1(cx); */
          counter = edgps[(fbmap[jj] - 1) << 1];

          /* 'intersection_triedge_cellface_2:30' yp = v1(cy); */
          yp = edgps[(fbmap[6 + jj] - 1) << 1];
        } else if (alpha == 1.0) {
          /* 'intersection_triedge_cellface_2:31' elseif (alpha == 1) */
          /* 'intersection_triedge_cellface_2:32' xp = v2(cx); */
          counter = edgps[1 + ((fbmap[jj] - 1) << 1)];

          /* 'intersection_triedge_cellface_2:33' yp = v2(cy); */
          yp = edgps[1 + ((fbmap[6 + jj] - 1) << 1)];
        } else {
          if ((alpha > 0.0) && (alpha < 1.0)) {
            /* 'intersection_triedge_cellface_2:34' elseif ((alpha>0) && (alpha<1)) */
            /* 'intersection_triedge_cellface_2:35' xp = v1(cx)+alpha*(v2(cx)-v1(cx)); */
            counter = edgps[(fbmap[jj] - 1) << 1] + alpha * (edgps[1 +
              ((fbmap[jj] - 1) << 1)] - edgps[(fbmap[jj] - 1) << 1]);

            /* 'intersection_triedge_cellface_2:36' yp = v1(cy)+alpha*(v2(cy)-v1(cy)); */
            yp = edgps[(fbmap[6 + jj] - 1) << 1] + alpha * (edgps[1 + ((fbmap[6
              + jj] - 1) << 1)] - edgps[(fbmap[6 + jj] - 1) << 1]);
          }
        }

        /* 'intersection_triedge_cellface_2:39' xi = (xp-cxmin)/(cxmax-cxmin); */
        xi = (counter - facebnd[jj]) / (facebnd[6 + jj] - facebnd[jj]);

        /* 'intersection_triedge_cellface_2:40' eta = (yp-cymin)/(cymax-cymin); */
        eta = (yp - facebnd[12 + jj]) / (facebnd[18 + jj] - facebnd[12 + jj]);

        /* 'intersection_triedge_cellface_2:42' if (abs(xi)<=tol) */
        if (fabs(xi) <= 1.0E-6) {
          /* 'intersection_triedge_cellface_2:43' xi =0 ; */
          xi = 0.0;
        } else {
          if (fabs(xi - 1.0) <= 1.0E-6) {
            /* 'intersection_triedge_cellface_2:44' elseif (abs(xi-1)<=tol) */
            /* 'intersection_triedge_cellface_2:45' xi = 1; */
            xi = 1.0;
          }
        }

        /* 'intersection_triedge_cellface_2:47' if (abs(eta)<=tol) */
        if (fabs(eta) <= 1.0E-6) {
          /* 'intersection_triedge_cellface_2:48' eta =0 ; */
          eta = 0.0;
        } else {
          if (fabs(eta - 1.0) <= 1.0E-6) {
            /* 'intersection_triedge_cellface_2:49' elseif (abs(eta-1)<=tol) */
            /* 'intersection_triedge_cellface_2:50' eta = 1; */
            eta = 1.0;
          }
        }

        /* 'intersection_triedge_cellface_2:53' if ((xi>=0) && (xi<=1) && (eta>=0) && (eta<=1)) */
        if ((xi >= 0.0) && (xi <= 1.0) && (eta >= 0.0) && (eta <= 1.0)) {
          /* 'intersection_triedge_cellface_2:54' [p,id]=assign_pntid(triID,edgID,edgvertID,alpha,xi,eta,xp,yp,jj,fbmap(jj,:),... */
          /* 'intersection_triedge_cellface_2:55'                     facebnd(jj,:),fcormap(jj,:),fedgmap(jj,:)); */
          for (i13 = 0; i13 < 3; i13++) {
            b_fbmap[i13] = fbmap[jj + 6 * i13];
          }

          for (i13 = 0; i13 < 4; i13++) {
            facecorID[i13] = (real_T)fcormap[jj + 6 * i13];
            faceedgID[i13] = (real_T)fedgmap[jj + 6 * i13];
          }

          /* % */
          /* 'intersection_triedge_cellface_2:134' cx = fbmap(1); */
          /* 'intersection_triedge_cellface_2:134' cxmin=facebnd(1); */
          /* 'intersection_triedge_cellface_2:134' cxmax = facebnd(2); */
          /* 'intersection_triedge_cellface_2:135' cy = fbmap(2); */
          /* 'intersection_triedge_cellface_2:135' cymin=facebnd(3); */
          /* 'intersection_triedge_cellface_2:135' cymax = facebnd(4); */
          /* 'intersection_triedge_cellface_2:136' cz = fbmap(3); */
          /* 'intersection_triedge_cellface_2:136' czval = facebnd(5); */
          /* 'intersection_triedge_cellface_2:137' id = [0 0 0 0 0 0]; */
          for (i13 = 0; i13 < 6; i13++) {
            id[i13] = 0.0;
          }

          /* 'intersection_triedge_cellface_2:138' p = [0 0 0]; */
          for (i13 = 0; i13 < 3; i13++) {
            p[i13] = 0.0;
          }

          /* 'intersection_triedge_cellface_2:139' p(cz) = czval; */
          p[b_fbmap[2] - 1] = facebnd[24 + jj];

          /* 'intersection_triedge_cellface_2:140' if (alpha ==0) */
          if (alpha == 0.0) {
            /* 'intersection_triedge_cellface_2:141' id(1) = edgvertID(1); */
            id[0] = (real_T)edgvertID[0];

            /* 'intersection_triedge_cellface_2:142' [p,id] = assign_id_basedon_xieta(xi,eta,xp,yp,cx,cxmin,cxmax,cy,cymin,cymax,faceid,fcormap,fedgmap,p,id); */
            assign_id_basedon_xieta(xi, eta, counter, yp, (real_T)b_fbmap[0],
              facebnd[jj], facebnd[6 + jj], (real_T)b_fbmap[1], facebnd[12 + jj],
              facebnd[18 + jj], 1.0 + (real_T)jj, facecorID, faceedgID, p, id);
          } else if (alpha == 1.0) {
            /* 'intersection_triedge_cellface_2:144' elseif (alpha == 1) */
            /* 'intersection_triedge_cellface_2:145' id(1) = edgvertID(2); */
            id[0] = (real_T)edgvertID[1];

            /* 'intersection_triedge_cellface_2:146' [p,id] = assign_id_basedon_xieta(xi,eta,xp,yp,cx,cxmin,cxmax,cy,cymin,cymax,faceid,fcormap,fedgmap,p,id); */
            assign_id_basedon_xieta(xi, eta, counter, yp, (real_T)b_fbmap[0],
              facebnd[jj], facebnd[6 + jj], (real_T)b_fbmap[1], facebnd[12 + jj],
              facebnd[18 + jj], 1.0 + (real_T)jj, facecorID, faceedgID, p, id);
          } else {
            if ((alpha > 0.0) && (alpha < 1.0)) {
              /* 'intersection_triedge_cellface_2:148' elseif (alpha >0 &&  alpha < 1) */
              /* 'intersection_triedge_cellface_2:149' id(2) = triID; */
              id[1] = (real_T)triID;

              /* 'intersection_triedge_cellface_2:150' id(3) = edgID; */
              id[2] = edgID;

              /* 'intersection_triedge_cellface_2:151' [p,id] = assign_id_basedon_xieta(xi,eta,xp,yp,cx,cxmin,cxmax,cy,cymin,cymax,faceid,fcormap,fedgmap,p,id); */
              assign_id_basedon_xieta(xi, eta, counter, yp, (real_T)b_fbmap[0],
                facebnd[jj], facebnd[6 + jj], (real_T)b_fbmap[1], facebnd[12 +
                jj], facebnd[18 + jj], 1.0 + (real_T)jj, facecorID, faceedgID, p,
                id);
            }
          }

          /* 'intersection_triedge_cellface_2:56' count = int32(0); */
          count = 0;

          /* 'intersection_triedge_cellface_2:57' for i=1:size(intpoints,1) */
          b_ps_idx = 0;
          exitg2 = 0U;
          while ((exitg2 == 0U) && (b_ps_idx < 5)) {
            /* 'intersection_triedge_cellface_2:58' if  id(1) == intpoints_id(i,1) && id(2) == intpoints_id(i,2) && id(3) == intpoints_id(i,3) && id(4) == intpoints_id(i,4) && ... */
            /* 'intersection_triedge_cellface_2:59'                             id(5) == intpoints_id(i,5) && id(6) == intpoints_id(i,6) */
            if ((id[0] == (real_T)intpoints_id[b_ps_idx]) && (id[1] == (real_T)
                 intpoints_id[5 + b_ps_idx]) && (id[2] == (real_T)intpoints_id
                 [10 + b_ps_idx]) && (id[3] == (real_T)intpoints_id[15 +
                 b_ps_idx]) && (id[4] == (real_T)intpoints_id[20 + b_ps_idx]) &&
                (id[5] == (real_T)intpoints_id[25 + b_ps_idx])) {
              /* 'intersection_triedge_cellface_2:60' count = int32(1); */
              count = 1;
              exitg2 = 1U;
            } else {
              b_ps_idx++;
            }
          }

          /* 'intersection_triedge_cellface_2:63' if (count == 0) */
          if (count == 0) {
            /* 'intersection_triedge_cellface_2:64' intpoints(idx,:) = p; */
            for (i13 = 0; i13 < 3; i13++) {
              intpoints[idx + 5 * i13] = p[i13];
            }

            /* 'intersection_triedge_cellface_2:65' intpoints_id(idx,:) = id; */
            for (i13 = 0; i13 < 6; i13++) {
              counter = id[i13];
              counter = counter < 0.0 ? ceil(counter - 0.5) : floor(counter +
                0.5);
              intpoints_id[idx + 5 * i13] = (int32_T)counter;
            }

            /* 'intersection_triedge_cellface_2:66' idx = idx + 1; */
            idx++;
          }
        }
      }
    }
  }

  /*  Checking for duplicates */
  /* 'intersection_triedge_cellface_2:75' intpoints = intpoints(1:idx-1,:); */
  if (1 > idx) {
    i13 = 0;
  } else {
    i13 = idx;
  }

  emxInit_real_T(&b_intpoints, 2);
  b_ps_idx = b_intpoints->size[0] * b_intpoints->size[1];
  b_intpoints->size[0] = i13;
  b_intpoints->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)b_intpoints, b_ps_idx, (int32_T)sizeof
                    (real_T));
  for (b_ps_idx = 0; b_ps_idx < 3; b_ps_idx++) {
    count = i13 - 1;
    for (jj = 0; jj <= count; jj++) {
      b_intpoints->data[jj + b_intpoints->size[0] * b_ps_idx] = intpoints[jj + 5
        * b_ps_idx];
    }
  }

  /* 'intersection_triedge_cellface_2:76' intpoints_id = intpoints_id(1:idx-1,:); */
  if (1 > idx) {
    i13 = 0;
  } else {
    i13 = idx;
  }

  emxInit_int32_T(&b_intpoints_id, 2);
  b_ps_idx = b_intpoints_id->size[0] * b_intpoints_id->size[1];
  b_intpoints_id->size[0] = i13;
  b_intpoints_id->size[1] = 6;
  emxEnsureCapacity((emxArray__common *)b_intpoints_id, b_ps_idx, (int32_T)
                    sizeof(int32_T));
  for (b_ps_idx = 0; b_ps_idx < 6; b_ps_idx++) {
    count = i13 - 1;
    for (jj = 0; jj <= count; jj++) {
      b_intpoints_id->data[jj + b_intpoints_id->size[0] * b_ps_idx] =
        intpoints_id[jj + 5 * b_ps_idx];
    }
  }

  /* 'intersection_triedge_cellface_2:77' index = idx-1; */
  /* 'intersection_triedge_cellface_2:78' if (index ~= 0) */
  if (idx != 0) {
    b_emxInit_int32_T(&samefacepnts, 1);

    /* 'intersection_triedge_cellface_2:79' samefacepnts = zeros(index,1,'int32'); */
    i13 = samefacepnts->size[0];
    samefacepnts->size[0] = idx;
    emxEnsureCapacity((emxArray__common *)samefacepnts, i13, (int32_T)sizeof
                      (int32_T));
    count = idx - 1;
    for (i13 = 0; i13 <= count; i13++) {
      samefacepnts->data[i13] = 0;
    }

    /* 'intersection_triedge_cellface_2:80' for ii=1:6 */
    b_emxInit_int32_T(&b_samefacepnts, 1);
    b_emxInit_int32_T(&c_samefacepnts, 1);
    for (ii = 0; ii < 6; ii++) {
      /* 'intersection_triedge_cellface_2:81' facecorID = fcormap(ii,:); */
      for (i13 = 0; i13 < 4; i13++) {
        facecorID[i13] = (real_T)fcormap[ii + 6 * i13];

        /* 'intersection_triedge_cellface_2:82' faceedgID = fedgmap(ii,:); */
        faceedgID[i13] = (real_T)fedgmap[ii + 6 * i13];
      }

      /* 'intersection_triedge_cellface_2:83' counter = 0; */
      counter = 0.0;

      /* 'intersection_triedge_cellface_2:84' for jj=1:index */
      for (jj = 0; jj + 1 <= idx; jj++) {
        /* 'intersection_triedge_cellface_2:85' if (intpoints_id(jj,4) == ii || intpoints_id(jj,5)==faceedgID(1) || intpoints_id(jj,5)==faceedgID(2) || intpoints_id(jj,5)==faceedgID(3) || intpoints_id(jj,5)==faceedgID(4)|| ... */
        /* 'intersection_triedge_cellface_2:86'                     intpoints_id(jj,6) == facecorID(1) || intpoints_id(jj,6) == facecorID(2)||intpoints_id(jj,6) == facecorID(3)||intpoints_id(jj,6) == facecorID(4)) */
        if ((b_intpoints_id->data[jj + b_intpoints_id->size[0] * 3] == 1 + ii) ||
            (b_intpoints_id->data[jj + (b_intpoints_id->size[0] << 2)] ==
             (int32_T)faceedgID[0]) || (b_intpoints_id->data[jj +
             (b_intpoints_id->size[0] << 2)] == (int32_T)faceedgID[1]) ||
            (b_intpoints_id->data[jj + (b_intpoints_id->size[0] << 2)] ==
             (int32_T)faceedgID[2]) || (b_intpoints_id->data[jj +
             (b_intpoints_id->size[0] << 2)] == (int32_T)faceedgID[3]) ||
            (b_intpoints_id->data[jj + b_intpoints_id->size[0] * 5] == (int32_T)
             facecorID[0]) || (b_intpoints_id->data[jj + b_intpoints_id->size[0]
             * 5] == (int32_T)facecorID[1]) || (b_intpoints_id->data[jj +
             b_intpoints_id->size[0] * 5] == (int32_T)facecorID[2]) ||
            (b_intpoints_id->data[jj + b_intpoints_id->size[0] * 5] == (int32_T)
             facecorID[3])) {
          /* 'intersection_triedge_cellface_2:87' samefacepnts(jj) = int32(jj); */
          samefacepnts->data[jj] = jj + 1;

          /* 'intersection_triedge_cellface_2:88' counter = counter + 1; */
          counter++;
        }
      }

      /* 'intersection_triedge_cellface_2:91' samefacepnts = unique(nonzeros(samefacepnts),'rows'); */
      i13 = c_samefacepnts->size[0];
      c_samefacepnts->size[0] = samefacepnts->size[0];
      emxEnsureCapacity((emxArray__common *)c_samefacepnts, i13, (int32_T)sizeof
                        (int32_T));
      count = samefacepnts->size[0] - 1;
      for (i13 = 0; i13 <= count; i13++) {
        c_samefacepnts->data[i13] = samefacepnts->data[i13];
      }

      nonzeros(c_samefacepnts, samefacepnts);
      i13 = b_samefacepnts->size[0];
      b_samefacepnts->size[0] = samefacepnts->size[0];
      emxEnsureCapacity((emxArray__common *)b_samefacepnts, i13, (int32_T)sizeof
                        (int32_T));
      count = samefacepnts->size[0] - 1;
      for (i13 = 0; i13 <= count; i13++) {
        b_samefacepnts->data[i13] = samefacepnts->data[i13];
      }

      unique(b_samefacepnts, samefacepnts);

      /* 'intersection_triedge_cellface_2:92' if (counter == 2) */
      if (counter == 2.0) {
        /* 'intersection_triedge_cellface_2:93' v1 = intpoints_id(samefacepnts(1),:); */
        /* 'intersection_triedge_cellface_2:94' v2 = intpoints_id(samefacepnts(2),:); */
        /* 'intersection_triedge_cellface_2:95' if (v1(4) ~= 0 && v2(4) == 0) */
        b_ps_idx = samefacepnts->data[0];
        guard1 = FALSE;
        if (b_intpoints_id->data[(b_ps_idx + b_intpoints_id->size[0] * 3) - 1]
            != 0) {
          b_ps_idx = samefacepnts->data[1];
          if (b_intpoints_id->data[(b_ps_idx + b_intpoints_id->size[0] * 3) - 1]
              == 0) {
            /* 'intersection_triedge_cellface_2:96' intpoints(samefacepnts(1),:) = []; */
            b_eml_null_assignment(b_intpoints, samefacepnts->data[0]);

            /* 'intersection_triedge_cellface_2:97' intpoints_id(samefacepnts(1),:) = []; */
            eml_null_assignment(b_intpoints_id, samefacepnts->data[0]);

            /* 'intersection_triedge_cellface_2:98' index = index - 1; */
            idx--;
          } else {
            guard1 = TRUE;
          }
        } else {
          guard1 = TRUE;
        }

        if (guard1 == TRUE) {
          b_ps_idx = samefacepnts->data[0];
          if (b_intpoints_id->data[(b_ps_idx + b_intpoints_id->size[0] * 3) - 1]
              == 0) {
            b_ps_idx = samefacepnts->data[1];
            if (b_intpoints_id->data[(b_ps_idx + b_intpoints_id->size[0] * 3) -
                1] != 0) {
              /* 'intersection_triedge_cellface_2:99' elseif (v1(4) == 0 && v2(4) ~= 0) */
              /* 'intersection_triedge_cellface_2:100' intpoints(samefacepnts(2),:) = []; */
              b_eml_null_assignment(b_intpoints, samefacepnts->data[1]);

              /* 'intersection_triedge_cellface_2:101' intpoints_id(samefacepnts(2),:) = []; */
              eml_null_assignment(b_intpoints_id, samefacepnts->data[1]);

              /* 'intersection_triedge_cellface_2:102' index = index - 1; */
              idx--;
            }
          }
        }
      }
    }

    emxFree_int32_T(&c_samefacepnts);
    emxFree_int32_T(&b_samefacepnts);
    emxFree_int32_T(&samefacepnts);

    /*  Putting points into ps_local */
    /* 'intersection_triedge_cellface_2:108' for kk=1:index */
    for (jj = 0; jj + 1 <= idx; jj++) {
      /* 'intersection_triedge_cellface_2:109' p = intpoints(kk,:); */
      /* 'intersection_triedge_cellface_2:110' id = intpoints_id(kk,:); */
      /* 'intersection_triedge_cellface_2:111' count = int32(0); */
      count = 0;

      /* 'intersection_triedge_cellface_2:112' for ii=1:ps_idx-1 */
      i13 = *ps_idx - 1;
      ii = 0;
      exitg1 = 0U;
      while ((exitg1 == 0U) && (ii + 1 <= i13)) {
        /* 'intersection_triedge_cellface_2:113' if (id(1)~=0)&&(id(1) == ps_localID(ii,1)) */
        if ((b_intpoints_id->data[jj] != 0) && (b_intpoints_id->data[jj] ==
             ps_localID->data[ii])) {
          /* 'intersection_triedge_cellface_2:114' count = int32(1); */
          count = 1;
          exitg1 = 1U;
        } else if ((b_intpoints_id->data[jj + b_intpoints_id->size[0]] != 0) &&
                   (b_intpoints_id->data[jj + (b_intpoints_id->size[0] << 1)] !=
                    0) && (b_intpoints_id->data[jj + (b_intpoints_id->size[0] <<
          2)] != 0) && (b_intpoints_id->data[jj + b_intpoints_id->size[0]] ==
                        ps_localID->data[ii + ps_localID->size[0]]) &&
                   (b_intpoints_id->data[jj + (b_intpoints_id->size[0] << 1)] ==
                    ps_localID->data[ii + (ps_localID->size[0] << 1)]) &&
                   (b_intpoints_id->data[jj + (b_intpoints_id->size[0] << 2)] ==
                    ps_localID->data[ii + (ps_localID->size[0] << 2)])) {
          /* 'intersection_triedge_cellface_2:116' elseif (id(2)~=0) && (id(3)~=0)&& (id(5)~=0) && (id(2)==ps_localID(ii,2) && id(3)==ps_localID(ii,3) && id(5) == ps_localID(ii,5)) */
          /* 'intersection_triedge_cellface_2:117' count = int32(1); */
          count = 1;
          exitg1 = 1U;
        } else if ((b_intpoints_id->data[jj + b_intpoints_id->size[0]] != 0) &&
                   (b_intpoints_id->data[jj + (b_intpoints_id->size[0] << 1)] !=
                    0) && (b_intpoints_id->data[jj + b_intpoints_id->size[0] * 5]
                           != 0) && (b_intpoints_id->data[jj +
                    b_intpoints_id->size[0]] == ps_localID->data[ii +
                    ps_localID->size[0]]) && (b_intpoints_id->data[jj +
                    (b_intpoints_id->size[0] << 1)] == ps_localID->data[ii +
                    (ps_localID->size[0] << 1)]) && (b_intpoints_id->data[jj +
                    b_intpoints_id->size[0] * 5] == ps_localID->data[ii +
                    ps_localID->size[0] * 5])) {
          /* 'intersection_triedge_cellface_2:119' elseif (id(2)~=0) && (id(3)~=0)&& (id(6)~=0) && (id(2)==ps_localID(ii,2) && id(3)==ps_localID(ii,3) && id(6) == ps_localID(ii,6)) */
          /* 'intersection_triedge_cellface_2:120' count = int32(1); */
          count = 1;
          exitg1 = 1U;
        } else {
          ii++;
        }
      }

      /* 'intersection_triedge_cellface_2:124' if (count ==0) */
      if (count == 0) {
        /* 'intersection_triedge_cellface_2:125' ps_local(ps_idx,:) = p; */
        b_ps_idx = *ps_idx;
        for (i13 = 0; i13 < 3; i13++) {
          ps_local->data[(b_ps_idx + ps_local->size[0] * i13) - 1] =
            b_intpoints->data[jj + b_intpoints->size[0] * i13];
        }

        /* 'intersection_triedge_cellface_2:126' ps_localID(ps_idx,:)=id; */
        b_ps_idx = *ps_idx;
        for (i13 = 0; i13 < 6; i13++) {
          ps_localID->data[(b_ps_idx + ps_localID->size[0] * i13) - 1] =
            b_intpoints_id->data[jj + b_intpoints_id->size[0] * i13];
        }

        /* 'intersection_triedge_cellface_2:127' ps_idx = ps_idx + 1; */
        (*ps_idx)++;
      }
    }
  }

  emxFree_int32_T(&b_intpoints_id);
  emxFree_real_T(&b_intpoints);
}

/*
 * function [face_type,face_sign] = label_cellface(ps_localID,tris_local,tris_bndID,tets,TETcomp,TETsign,num_comps,comp_list,pos_map,neg_map,face_type,face_sign)
 * This function outputs the labels of the cell faces
 * face_type : 1x6 array
 *  face_type = -2 by default,
 *            = -1 if it is a crossed face,
 *            = othervalue, otherwise
 *  face_sign : bool
 *  face_sign = true/1 : default or positive sign
 *            = false/0 : negative sign
 */
static void label_cellface(const emxArray_int32_T *ps_localID, const
  emxArray_int32_T *tris_local, const emxArray_int32_T *tris_bndID, const
  emxArray_int32_T *tets, const emxArray_int32_T *TETcomp, const
  emxArray_int32_T *TETsign, int32_T num_comps, const emxArray_int32_T
  *comp_list, const emxArray_int32_T *pos_map, const emxArray_int32_T *neg_map,
  int32_T face_type[6], boolean_T face_sign[6])
{
  int32_T tetnum;
  emxArray_int32_T *tlist;
  emxArray_int32_T *tcomp;
  emxArray_int32_T *tsign;
  emxArray_int32_T *SIGN;
  emxArray_int32_T *b_tcomp;
  emxArray_int32_T *c_tcomp;
  emxArray_int32_T *b_tsign;
  int32_T ii;
  int32_T idx;
  int32_T jj;
  static const int8_T FaceEdgeMap[24] = { 4, 1, 2, 3, 5, 9, 3, 6, 7, 8, 12, 10,
    2, 9, 10, 11, 8, 11, 1, 5, 6, 7, 4, 12 };

  static const int8_T Faces[24] = { 1, 1, 2, 3, 1, 5, 4, 2, 3, 4, 5, 6, 3, 6, 7,
    8, 8, 7, 2, 5, 6, 7, 4, 8 };

  int32_T i19;
  int32_T p;
  int32_T mtmp;
  boolean_T exitg1;
  int32_T wall[12];
  int32_T ix;
  boolean_T exitg2;
  int32_T b_wall[3];
  int32_T b_tris_local[3];
  int32_T b_tlist;
  real_T SO;
  real_T SV;

  /* 'label_cellface:10' coder.inline('never') */
  /* 'label_cellface:11' max_no_tris = size(tris_local,1); */
  /* 'label_cellface:12' tetnum = int32(size(tets,1)); */
  tetnum = tets->size[0];

  /* 'label_cellface:14' Faces = int32([1 4 3 2;1 2 6 5;2 3 7 6;3 4 8 7;1 5 8 4;5 6 7 8]); */
  /* 'label_cellface:15' FaceEdgeMap = int32([4,3,2,1;1,6,9,5;2,7,10,6;3,8,11,7;5,12,8,4;9,10,11,12]); */
  /* 'label_cellface:17' for ii=1:6 */
  b_emxInit_int32_T(&tlist, 1);
  b_emxInit_int32_T(&tcomp, 1);
  b_emxInit_int32_T(&tsign, 1);
  emxInit_int32_T(&SIGN, 2);
  emxInit_int32_T(&b_tcomp, 2);
  b_emxInit_int32_T(&c_tcomp, 1);
  emxInit_int32_T(&b_tsign, 2);
  for (ii = 0; ii < 6; ii++) {
    /* 'label_cellface:18' f_idx = int32(1); */
    idx = 0;

    /* 'label_cellface:19' for jj=1:size(ps_localID,1) */
    for (jj = 0; jj <= ps_localID->size[0] - 1; jj++) {
      /* 'label_cellface:20' if (ps_localID(jj,4)== ii) */
      if (ps_localID->data[jj + ps_localID->size[0] * 3] == 1 + ii) {
        /* 'label_cellface:21' f_idx = f_idx + 1; */
        idx++;
      } else if ((ps_localID->data[jj + (ps_localID->size[0] << 2)] ==
                  FaceEdgeMap[ii]) || (ps_localID->data[jj + (ps_localID->size[0]
        << 2)] == FaceEdgeMap[6 + ii]) || (ps_localID->data[jj +
                  (ps_localID->size[0] << 2)] == FaceEdgeMap[12 + ii]) ||
                 (ps_localID->data[jj + (ps_localID->size[0] << 2)] ==
                  FaceEdgeMap[18 + ii])) {
        /* 'label_cellface:22' elseif ((ps_localID(jj,5) == FaceEdgeMap(ii,1))||(ps_localID(jj,5) == FaceEdgeMap(ii,2))|| ... */
        /* 'label_cellface:23'                 (ps_localID(jj,5) == FaceEdgeMap(ii,3))||(ps_localID(jj,5) == FaceEdgeMap(ii,4))) */
        /* 'label_cellface:24' f_idx = f_idx + 1; */
        idx++;
      } else {
        if ((ps_localID->data[jj + ps_localID->size[0] * 5] == Faces[ii]) ||
            (ps_localID->data[jj + ps_localID->size[0] * 5] == Faces[6 + ii]) ||
            (ps_localID->data[jj + ps_localID->size[0] * 5] == Faces[12 + ii]) ||
            (ps_localID->data[jj + ps_localID->size[0] * 5] == Faces[18 + ii]))
        {
          /* 'label_cellface:25' elseif ((ps_localID(jj,6) == Faces(ii,1))||(ps_localID(jj,6) == Faces(ii,2))|| ... */
          /* 'label_cellface:26'                 (ps_localID(jj,6) == Faces(ii,3))||(ps_localID(jj,6) == Faces(ii,4))) */
          /* 'label_cellface:27' f_idx = f_idx + 1; */
          idx++;
        }
      }
    }

    /* 'label_cellface:30' index = f_idx - 1; */
    /* 'label_cellface:31' if (index >= 3) */
    if (idx >= 3) {
      /* 'label_cellface:32' face_type(ii) = -1; */
      face_type[ii] = -1;
    } else {
      /* 'label_cellface:33' else */
      /* 'label_cellface:34' tlist =zeros(max_no_tris,1,'int32'); */
      i19 = tlist->size[0];
      tlist->size[0] = tris_local->size[0];
      emxEnsureCapacity((emxArray__common *)tlist, i19, (int32_T)sizeof(int32_T));
      idx = tris_local->size[0] - 1;
      for (i19 = 0; i19 <= idx; i19++) {
        tlist->data[i19] = 0;
      }

      /* 'label_cellface:34' p = int32(1); */
      p = 0;

      /* 'label_cellface:35' for jj=1:size(tris_bndID,1) */
      for (jj = 0; jj <= tris_bndID->size[0] - 1; jj++) {
        /* 'label_cellface:36' if (tris_bndID(jj) == ii) */
        if (tris_bndID->data[jj] == 1 + ii) {
          /* 'label_cellface:37' tlist(p) = jj; */
          tlist->data[p] = jj + 1;

          /* 'label_cellface:38' p = p+1; */
          p++;
        }
      }

      /* 'label_cellface:42' tcomp = zeros(p-1,1,'int32'); */
      i19 = tcomp->size[0];
      tcomp->size[0] = p;
      emxEnsureCapacity((emxArray__common *)tcomp, i19, (int32_T)sizeof(int32_T));
      idx = p - 1;
      for (i19 = 0; i19 <= idx; i19++) {
        tcomp->data[i19] = 0;
      }

      /* 'label_cellface:43' tsign = zeros(p-1,1,'int32'); */
      i19 = tsign->size[0];
      tsign->size[0] = p;
      emxEnsureCapacity((emxArray__common *)tsign, i19, (int32_T)sizeof(int32_T));
      idx = p - 1;
      for (i19 = 0; i19 <= idx; i19++) {
        tsign->data[i19] = 0;
      }

      /* 'label_cellface:45' for jj=1:(p-1) */
      for (jj = 0; jj + 1 <= p; jj++) {
        /* 'label_cellface:46' if (tlist(jj) ==0) */
        if (tlist->data[jj] == 0) {
        } else {
          /* 'label_cellface:47' else */
          /* 'label_cellface:48' flag = 0; */
          idx = 0;

          /* 'label_cellface:49' for kk=1:tetnum */
          mtmp = 0;
          exitg1 = 0U;
          while ((exitg1 == 0U) && (mtmp + 1 <= tetnum)) {
            /* 'label_cellface:50' wall = int32(zeros(4,3)); */
            for (i19 = 0; i19 < 12; i19++) {
              wall[i19] = 0;
            }

            /* 'label_cellface:51' wall(1,:) = [tets(kk,1) tets(kk,3) tets(kk,2)]; */
            wall[0] = tets->data[mtmp];
            wall[4] = tets->data[mtmp + (tets->size[0] << 1)];
            wall[8] = tets->data[mtmp + tets->size[0]];

            /* 'label_cellface:52' wall(2,:) = [tets(kk,1) tets(kk,2) tets(kk,4)]; */
            wall[1] = tets->data[mtmp];
            wall[5] = tets->data[mtmp + tets->size[0]];
            wall[9] = tets->data[mtmp + tets->size[0] * 3];

            /* 'label_cellface:53' wall(3,:) = [tets(kk,2) tets(kk,3) tets(kk,4)]; */
            wall[2] = tets->data[mtmp + tets->size[0]];
            wall[6] = tets->data[mtmp + (tets->size[0] << 1)];
            wall[10] = tets->data[mtmp + tets->size[0] * 3];

            /* 'label_cellface:54' wall(4,:) = [tets(kk,3) tets(kk,1) tets(kk,4)]; */
            wall[3] = tets->data[mtmp + (tets->size[0] << 1)];
            wall[7] = tets->data[mtmp];
            wall[11] = tets->data[mtmp + tets->size[0] * 3];

            /* 'label_cellface:55' for i=1:4 */
            ix = 0;
            exitg2 = 0U;
            while ((exitg2 == 0U) && (ix < 4)) {
              /* 'label_cellface:56' [SV,SO] = chk_tri_orientation(wall(i,:),tris_local(tlist(jj),:)); */
              for (i19 = 0; i19 < 3; i19++) {
                b_wall[i19] = wall[ix + (i19 << 2)];
              }

              b_tlist = tlist->data[jj];
              for (i19 = 0; i19 < 3; i19++) {
                b_tris_local[i19] = tris_local->data[(b_tlist + tris_local->
                  size[0] * i19) - 1];
              }

              chk_tri_orientation(b_wall, b_tris_local, &SV, &SO);

              /* 'label_cellface:57' if ((SV == 1)&&(SO == 1)) */
              if ((SV == 1.0) && (SO == 1.0)) {
                /* 'label_cellface:58' tcomp(jj) = TETcomp(kk,5); */
                tcomp->data[jj] = TETcomp->data[mtmp + (TETcomp->size[0] << 2)];

                /* 'label_cellface:59' tsign(jj) = TETsign(kk); */
                tsign->data[jj] = TETsign->data[mtmp];

                /* 'label_cellface:60' flag = 1; */
                idx = 1;
                exitg2 = 1U;
              } else {
                ix++;
              }
            }

            /* 'label_cellface:64' if (flag == 1) */
            if (idx == 1) {
              exitg1 = 1U;
            } else {
              mtmp++;
            }
          }
        }
      }

      /* 'label_cellface:70' SIGN = unique(tsign'); */
      i19 = b_tsign->size[0] * b_tsign->size[1];
      b_tsign->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)b_tsign, i19, (int32_T)sizeof
                        (int32_T));
      idx = tsign->size[0];
      i19 = b_tsign->size[0] * b_tsign->size[1];
      b_tsign->size[1] = idx;
      emxEnsureCapacity((emxArray__common *)b_tsign, i19, (int32_T)sizeof
                        (int32_T));
      idx = tsign->size[0] - 1;
      for (i19 = 0; i19 <= idx; i19++) {
        b_tsign->data[i19] = tsign->data[i19];
      }

      e_unique(b_tsign, SIGN);

      /* 'label_cellface:71' if ((size(SIGN,1)==1) && (size(SIGN,2)==1)) */
      if (SIGN->size[1] == 1) {
        /* 'label_cellface:72' if (SIGN(1,1) == 1) */
        if (SIGN->data[0] == 1) {
          /* 'label_cellface:73' face_sign(ii)= 1; */
          face_sign[ii] = TRUE;
        } else {
          if (SIGN->data[0] == -1) {
            /* 'label_cellface:74' elseif (SIGN(1,1) == -1) */
            /* 'label_cellface:75' face_sign(ii)= 0; */
            face_sign[ii] = FALSE;
          }
        }

        /* 'label_cellface:77' tcomp = nonzeros(tcomp); */
        i19 = c_tcomp->size[0];
        c_tcomp->size[0] = tcomp->size[0];
        emxEnsureCapacity((emxArray__common *)c_tcomp, i19, (int32_T)sizeof
                          (int32_T));
        idx = tcomp->size[0] - 1;
        for (i19 = 0; i19 <= idx; i19++) {
          c_tcomp->data[i19] = tcomp->data[i19];
        }

        nonzeros(c_tcomp, tcomp);

        /* 'label_cellface:78' comp = unique(tcomp'); */
        i19 = b_tcomp->size[0] * b_tcomp->size[1];
        b_tcomp->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)b_tcomp, i19, (int32_T)sizeof
                          (int32_T));
        idx = tcomp->size[0];
        i19 = b_tcomp->size[0] * b_tcomp->size[1];
        b_tcomp->size[1] = idx;
        emxEnsureCapacity((emxArray__common *)b_tcomp, i19, (int32_T)sizeof
                          (int32_T));
        idx = tcomp->size[0] - 1;
        for (i19 = 0; i19 <= idx; i19++) {
          b_tcomp->data[i19] = tcomp->data[i19];
        }

        e_unique(b_tcomp, SIGN);

        /* 'label_cellface:79' face_type(ii) = min(comp); */
        idx = SIGN->size[1];
        mtmp = SIGN->data[0];
        if (idx > 1) {
          for (ix = 1; ix + 1 <= idx; ix++) {
            if (SIGN->data[ix] < mtmp) {
              mtmp = SIGN->data[ix];
            }
          }
        }

        face_type[ii] = mtmp;
      } else {
        /* 'label_cellface:80' else */
        /* 'label_cellface:81' face_type(ii) = -1; */
        face_type[ii] = -1;
      }
    }
  }

  emxFree_int32_T(&b_tsign);
  emxFree_int32_T(&c_tcomp);
  emxFree_int32_T(&b_tcomp);
  emxFree_int32_T(&SIGN);
  emxFree_int32_T(&tsign);
  emxFree_int32_T(&tcomp);
  emxFree_int32_T(&tlist);

  /*  Finally use the comp_list, pos_map, neg_map to get the comps of the */
  /*  faces. */
  /* 'label_cellface:88' idx = int32(0); */
  idx = -1;

  /* 'label_cellface:89' for ii=1:6 */
  for (ii = 0; ii < 6; ii++) {
    /* 'label_cellface:90' if (face_type(ii) ~= -1) */
    if (face_type[ii] != -1) {
      /* 'label_cellface:91' for jj=1:num_comps */
      for (jj = 0; jj + 1 <= num_comps; jj++) {
        /* 'label_cellface:92' if (comp_list(jj) == face_type(ii)) */
        if (comp_list->data[jj] == face_type[ii]) {
          /* 'label_cellface:93' idx = jj; */
          idx = jj;
        }
      }

      /* 'label_cellface:96' if (face_sign(ii) == 1) */
      if (face_sign[ii] == 1) {
        /* 'label_cellface:97' face_type(ii) = pos_map(idx); */
        face_type[ii] = pos_map->data[idx];
      } else {
        /* 'label_cellface:98' else */
        /* 'label_cellface:99' face_type(ii) = neg_map(idx); */
        face_type[ii] = neg_map->data[idx];
      }
    }
  }
}

/*
 * function [TetComp,TetSign,pos_map,neg_map] = label_tets(nv,tets,tris_labelled,tris_comp,num_comps,comp_list,pos_map,neg_map)
 * This function assigns labels to all the tets.
 *  The output should be labelled tets, and lists of local components
 *  All the tets are labelled using a depth first search.
 *  Tetsign = -1, for negative components,
 *          = 1, for positive components,
 *          = 0, otherwise : To distinguish that this tet has not been
 *          assigned a sign yet.
 * coder.extrinsic('fprintf')
 */
static void label_tets(int32_T nv, const emxArray_int32_T *tets, const
  emxArray_int32_T *tris_labelled, const emxArray_int32_T *tris_comp, int32_T
  num_comps, const emxArray_int32_T *comp_list, emxArray_int32_T *pos_map,
  emxArray_int32_T *neg_map, emxArray_int32_T *TetComp, emxArray_int32_T
  *TetSign)
{
  int32_T tnum;
  int32_T ntris;
  int32_T k;
  int32_T j;
  emxArray_int32_T *tetflag;
  emxArray_int32_T *stack;
  emxArray_int32_T *opphfs;
  int32_T stack_idx;
  int32_T counter;
  int32_T ii;
  boolean_T exitg4;
  int32_T F[12];
  real_T SV[4];
  real_T SO[4];
  int32_T b_F[3];
  int32_T b_tris_labelled[3];
  real_T b_SO;
  real_T b_SV;
  boolean_T exitg3;
  boolean_T guard2 = FALSE;
  emxArray_real_T *SIGN;
  emxArray_int32_T *T;
  emxArray_int32_T *comp;
  emxArray_real_T *L;
  emxArray_int32_T *b_T;
  emxArray_real_T *b_L;
  int32_T ngbtets[4];
  uint32_T a;
  int32_T next;
  boolean_T exitg2;
  boolean_T guard1 = FALSE;
  boolean_T exitg1;

  /* 'label_tets:11' tnum = int32(size(tets,1)); */
  tnum = tets->size[0];

  /* 'label_tets:12' ntris = int32(size(tris_labelled,1)); */
  ntris = tris_labelled->size[0];

  /* 'label_tets:13' TetComp = zeros(tnum,5,'int32'); */
  k = TetComp->size[0] * TetComp->size[1];
  TetComp->size[0] = tnum;
  TetComp->size[1] = 5;
  emxEnsureCapacity((emxArray__common *)TetComp, k, (int32_T)sizeof(int32_T));
  j = tnum * 5 - 1;
  for (k = 0; k <= j; k++) {
    TetComp->data[k] = 0;
  }

  /* 'label_tets:14' TetSign = zeros(tnum,1,'int32'); */
  k = TetSign->size[0];
  TetSign->size[0] = tnum;
  emxEnsureCapacity((emxArray__common *)TetSign, k, (int32_T)sizeof(int32_T));
  j = tnum - 1;
  for (k = 0; k <= j; k++) {
    TetSign->data[k] = 0;
  }

  b_emxInit_int32_T(&tetflag, 1);

  /* 'label_tets:15' tetflag = zeros(tnum,1,'int32'); */
  k = tetflag->size[0];
  tetflag->size[0] = tnum;
  emxEnsureCapacity((emxArray__common *)tetflag, k, (int32_T)sizeof(int32_T));
  j = tnum - 1;
  for (k = 0; k <= j; k++) {
    tetflag->data[k] = 0;
  }

  b_emxInit_int32_T(&stack, 1);

  /* 'label_tets:16' stack = zeros(tnum,1,'int32'); */
  k = stack->size[0];
  stack->size[0] = tnum;
  emxEnsureCapacity((emxArray__common *)stack, k, (int32_T)sizeof(int32_T));
  j = tnum - 1;
  for (k = 0; k <= j; k++) {
    stack->data[k] = 0;
  }

  emxInit_int32_T(&opphfs, 2);

  /* 'label_tets:16' stack_idx = int32(1); */
  stack_idx = 1;

  /* 'label_tets:17' counter = int32(1); */
  counter = 1;

  /* 'label_tets:19' opphfs = determine_opposite_halfface(nv,tets); */
  determine_opposite_halfface(nv, tets, opphfs);

  /* % Find the first element of the stack */
  /* 'label_tets:21' for ii=1:tnum */
  ii = 0;
  exitg4 = 0U;
  while ((exitg4 == 0U) && (ii + 1 <= tnum)) {
    /* 'label_tets:22' F = int32(zeros(4,3)); */
    for (k = 0; k < 12; k++) {
      F[k] = 0;
    }

    /* 'label_tets:23' F(1,:) = [tets(ii,1) tets(ii,3) tets(ii,2)]; */
    F[0] = tets->data[ii];
    F[4] = tets->data[ii + (tets->size[0] << 1)];
    F[8] = tets->data[ii + tets->size[0]];

    /* 'label_tets:24' F(2,:) = [tets(ii,1) tets(ii,2) tets(ii,4)]; */
    F[1] = tets->data[ii];
    F[5] = tets->data[ii + tets->size[0]];
    F[9] = tets->data[ii + tets->size[0] * 3];

    /* 'label_tets:25' F(3,:) = [tets(ii,2) tets(ii,3) tets(ii,4)]; */
    F[2] = tets->data[ii + tets->size[0]];
    F[6] = tets->data[ii + (tets->size[0] << 1)];
    F[10] = tets->data[ii + tets->size[0] * 3];

    /* 'label_tets:26' F(4,:) = [tets(ii,3) tets(ii,1) tets(ii,4)]; */
    F[3] = tets->data[ii + (tets->size[0] << 1)];
    F[7] = tets->data[ii];
    F[11] = tets->data[ii + tets->size[0] * 3];

    /* 'label_tets:28' SV = zeros(4,1); */
    /* 'label_tets:28' SO = zeros(4,1); */
    /* 'label_tets:29' for kk =1:4 */
    for (j = 0; j < 4; j++) {
      /* 'label_tets:30' [SV(kk),SO(kk)] = chk_tri_orientation(F(kk,:),tris_labelled(1,:)); */
      for (k = 0; k < 3; k++) {
        b_F[k] = F[j + (k << 2)];
      }

      for (k = 0; k < 3; k++) {
        b_tris_labelled[k] = tris_labelled->data[tris_labelled->size[0] * k];
      }

      chk_tri_orientation(b_F, b_tris_labelled, &b_SV, &b_SO);
      SV[j] = b_SV;
      SO[j] = b_SO;
    }

    /* 'label_tets:32' if ((nnz(SV) == 1) && (nnz(SO) == 1)) */
    if ((nnz(SV) == 1.0) && (nnz(SO) == 1.0)) {
      /* 'label_tets:33' TetComp(ii,:) = tris_comp(1); */
      for (k = 0; k < 5; k++) {
        TetComp->data[ii + TetComp->size[0] * k] = tris_comp->data[0];
      }

      /* 'label_tets:34' TetSign(ii) = -1; */
      TetSign->data[ii] = -1;

      /* 'label_tets:35' tetflag(ii) = 1; */
      tetflag->data[ii] = 1;

      /* 'label_tets:36' stack(counter) = ii; */
      stack->data[0] = ii + 1;
      exitg4 = 1U;
    } else {
      ii++;
    }
  }

  /* % If stack is still empty */
  /* 'label_tets:42' if (stack(1) == 0) */
  if (stack->data[0] == 0) {
    /* 'label_tets:43' for ii=1:tnum */
    ii = 0;
    exitg3 = 0U;
    while ((exitg3 == 0U) && (ii + 1 <= tnum)) {
      /* 'label_tets:44' F = int32(zeros(4,3)); */
      for (k = 0; k < 12; k++) {
        F[k] = 0;
      }

      /* 'label_tets:45' F(1,:) = [tets(ii,1) tets(ii,3) tets(ii,2)]; */
      F[0] = tets->data[ii];
      F[4] = tets->data[ii + (tets->size[0] << 1)];
      F[8] = tets->data[ii + tets->size[0]];

      /* 'label_tets:46' F(2,:) = [tets(ii,1) tets(ii,2) tets(ii,4)]; */
      F[1] = tets->data[ii];
      F[5] = tets->data[ii + tets->size[0]];
      F[9] = tets->data[ii + tets->size[0] * 3];

      /* 'label_tets:47' F(3,:) = [tets(ii,2) tets(ii,3) tets(ii,4)]; */
      F[2] = tets->data[ii + tets->size[0]];
      F[6] = tets->data[ii + (tets->size[0] << 1)];
      F[10] = tets->data[ii + tets->size[0] * 3];

      /* 'label_tets:48' F(4,:) = [tets(ii,3) tets(ii,1) tets(ii,4)]; */
      F[3] = tets->data[ii + (tets->size[0] << 1)];
      F[7] = tets->data[ii];
      F[11] = tets->data[ii + tets->size[0] * 3];

      /* 'label_tets:50' SV = zeros(4,1); */
      /* 'label_tets:50' SO = zeros(4,1); */
      /* 'label_tets:51' for kk =1:4 */
      for (j = 0; j < 4; j++) {
        /* 'label_tets:52' [SV(kk),SO(kk)] = chk_tri_orientation(F(kk,:),tris_labelled(1,:)); */
        for (k = 0; k < 3; k++) {
          b_F[k] = F[j + (k << 2)];
        }

        for (k = 0; k < 3; k++) {
          b_tris_labelled[k] = tris_labelled->data[tris_labelled->size[0] * k];
        }

        chk_tri_orientation(b_F, b_tris_labelled, &b_SV, &b_SO);
        SV[j] = b_SV;
        SO[j] = b_SO;
      }

      /* 'label_tets:54' if ((nnz(SV) == 1) && (nnz(SO) == 0)) */
      guard2 = FALSE;
      if (nnz(SV) == 1.0) {
        j = 0;
        for (k = 0; k < 4; k++) {
          if (SO[k] != 0.0) {
            j++;
          }
        }

        if (j == 0) {
          /* 'label_tets:55' TetComp(ii,:) = tris_comp(1); */
          for (k = 0; k < 5; k++) {
            TetComp->data[ii + TetComp->size[0] * k] = tris_comp->data[0];
          }

          /* 'label_tets:56' TetSign(ii) = 1; */
          TetSign->data[ii] = 1;

          /* 'label_tets:57' tetflag(ii) = 1; */
          tetflag->data[ii] = 1;

          /* 'label_tets:58' stack(counter) = ii; */
          stack->data[0] = ii + 1;
          exitg3 = 1U;
        } else {
          guard2 = TRUE;
        }
      } else {
        guard2 = TRUE;
      }

      if (guard2 == TRUE) {
        ii++;
      }
    }
  }

  /* % Begin Depth First Search */
  /* 'label_tets:67' while (1) */
  emxInit_real_T(&SIGN, 2);
  b_emxInit_int32_T(&T, 1);
  emxInit_int32_T(&comp, 2);
  b_emxInit_real_T(&L, 1);
  emxInit_int32_T(&b_T, 2);
  emxInit_real_T(&b_L, 2);
  do {
    /* 'label_tets:68' top = stack(stack_idx); */
    /* 'label_tets:69' [ngbtets] = get_ngb_tets(top,opphfs); */
    /* 'get_ngb_tets:2' ngbtets = zeros(1,4,'int32'); */
    /* 'get_ngb_tets:3' ngbtetfs = zeros(1,4,'int32'); */
    /* 'get_ngb_tets:4' for ii=1:4 */
    for (ii = 0; ii < 4; ii++) {
      /* 'get_ngb_tets:5' ngbtets(ii) = int32(bitshift( uint32(opphfs(keyID,ii)),-3)); */
      a = (uint32_T)opphfs->data[(stack->data[stack_idx - 1] + opphfs->size[0] *
        ii) - 1];
      ngbtets[ii] = (int32_T)(a >> 3U);

      /* 'get_ngb_tets:6' if (ngbtets(ii) == 0) */
    }

    /* 'label_tets:71' next = get_next_tet(ngbtets,tetflag); */
    /*  ngbtets -- Id of the four neighbor tets */
    /* 'get_next_tet:4' count = int32(0); */
    j = 0;

    /* 'get_next_tet:4' next = int32(1); */
    next = 0;

    /* 'get_next_tet:5' for ii=1:size(ngbtets,2) */
    ii = 0;
    exitg2 = 0U;
    while ((exitg2 == 0U) && (ii < 4)) {
      /* 'get_next_tet:6' if (ngbtets(ii) == 0) */
      guard1 = FALSE;
      if (ngbtets[ii] == 0) {
        /* 'get_next_tet:7' count = count + 1; */
        j++;
        guard1 = TRUE;
      } else if (tetflag->data[ngbtets[ii] - 1] == 1) {
        /* 'get_next_tet:8' elseif (tetflag(ngbtets(ii))==1) */
        /* 'get_next_tet:9' count = count + 1; */
        j++;
        guard1 = TRUE;
      } else {
        /* 'get_next_tet:10' else */
        /* 'get_next_tet:11' next = ngbtets(ii); */
        next = ngbtets[ii] - 1;
        exitg2 = 1U;
      }

      if (guard1 == TRUE) {
        ii++;
      }
    }

    /* 'get_next_tet:14' if (count == size(ngbtets,2)) */
    if (j == 4) {
      /* 'get_next_tet:15' next = int32(-1); */
      next = -2;
    }

    /*  next = int32(-1); */
    /*  for ii=1:size(ngbtets,2) */
    /*      if (ngbtets(ii) == 0)         */
    /*      elseif (tetflag(ngbtets(ii))==1)         */
    /*      else */
    /*          next = ngbtets(ii);break */
    /*      end */
    /*  end  */
    /*  end */
    /* 'label_tets:72' if (next == -1) */
    if (next + 1 == -1) {
      /* 'label_tets:73' stack_idx = stack_idx - 1; */
      stack_idx--;
    } else {
      /* 'label_tets:74' else */
      /* 'label_tets:75' counter = counter + 1; */
      counter++;

      /* 'label_tets:76' stack_idx = counter; */
      stack_idx = counter;

      /* 'label_tets:77' stack(counter) = next; */
      stack->data[counter - 1] = next + 1;

      /*  Color the walls of this tet */
      /* 'label_tets:80' wall = int32(zeros(4,3)); */
      for (k = 0; k < 12; k++) {
        F[k] = 0;
      }

      /* 'label_tets:81' wall(1,:) = [tets(next,1) tets(next,3) tets(next,2)]; */
      F[0] = tets->data[next];
      F[4] = tets->data[next + (tets->size[0] << 1)];
      F[8] = tets->data[next + tets->size[0]];

      /* 'label_tets:82' wall(2,:) = [tets(next,1) tets(next,2) tets(next,4)]; */
      F[1] = tets->data[next];
      F[5] = tets->data[next + tets->size[0]];
      F[9] = tets->data[next + tets->size[0] * 3];

      /* 'label_tets:83' wall(3,:) = [tets(next,2) tets(next,3) tets(next,4)]; */
      F[2] = tets->data[next + tets->size[0]];
      F[6] = tets->data[next + (tets->size[0] << 1)];
      F[10] = tets->data[next + tets->size[0] * 3];

      /* 'label_tets:84' wall(4,:) = [tets(next,3) tets(next,1) tets(next,4)]; */
      F[3] = tets->data[next + (tets->size[0] << 1)];
      F[7] = tets->data[next];
      F[11] = tets->data[next + tets->size[0] * 3];

      /* 'label_tets:85' L = zeros(1,4); */
      for (k = 0; k < 4; k++) {
        SV[k] = 0.0;
      }

      /* 'label_tets:87' for kk =1:4 */
      for (j = 0; j < 4; j++) {
        /* 'label_tets:88' for ii=1:ntris */
        ii = 0;
        exitg1 = 0U;
        while ((exitg1 == 0U) && (ii + 1 <= ntris)) {
          /* 'label_tets:89' [SV,SO] = chk_tri_orientation(wall(kk,:),tris_labelled(ii,:)); */
          for (k = 0; k < 3; k++) {
            b_F[k] = F[j + (k << 2)];
          }

          for (k = 0; k < 3; k++) {
            b_tris_labelled[k] = tris_labelled->data[ii + tris_labelled->size[0]
              * k];
          }

          chk_tri_orientation(b_F, b_tris_labelled, &b_SV, &b_SO);

          /* 'label_tets:90' if ((SV == 1)&&(SO == 1)) */
          if ((b_SV == 1.0) && (b_SO == 1.0)) {
            /* 'label_tets:91' TetComp(next,kk) = tris_comp(ii); */
            TetComp->data[next + TetComp->size[0] * j] = tris_comp->data[ii];

            /* 'label_tets:92' L(kk) = -1; */
            SV[j] = -1.0;
            exitg1 = 1U;
          } else if ((b_SV == 1.0) && (b_SO == 0.0)) {
            /* 'label_tets:94' elseif ((SV == 1)&&(SO ==0)) */
            /* 'label_tets:95' TetComp(next,kk) = tris_comp(ii); */
            TetComp->data[next + TetComp->size[0] * j] = tris_comp->data[ii];

            /* 'label_tets:96' L(kk) = 1; */
            SV[j] = 1.0;
            exitg1 = 1U;
          } else {
            ii++;
          }
        }
      }

      /* 'label_tets:102' for ii=1:4 */
      for (ii = 0; ii < 4; ii++) {
        /* 'label_tets:103' if (TetComp(next,ii)==0) */
        if (TetComp->data[next + TetComp->size[0] * ii] == 0) {
          /* 'label_tets:104' [Opp,Oppfs] = get_opposite_tet(next,ii,opphfs); */
          /*  Input : Index of a TET, a face of the given TET, and list of all TETS */
          /*  Output: Index of the TET opposite to the given face, and Index of the */
          /*  face which is the opposite face. */
          /* 'get_opposite_tet:5' Opptetfs = int32(0); */
          j = -1;

          /* 'get_opposite_tet:6' Opptet = int32(bitshift( uint32(opphfs(tetid,faceid)),-3)); */
          a = (uint32_T)opphfs->data[next + opphfs->size[0] * ii];
          k = (int32_T)(a >> 3U);

          /* 'get_opposite_tet:7' if (Opptet == 0) */
          if (k == 0) {
          } else {
            /* 'get_opposite_tet:8' else */
            /* 'get_opposite_tet:9' Opptetfs = mod(opphfs(tetid,faceid),8)+1; */
            j = opphfs->data[next + opphfs->size[0] * ii] - ((opphfs->data[next
              + opphfs->size[0] * ii] >> 3) << 3);
          }

          /* 'label_tets:105' if (Opp == 0) */
          if (k == 0) {
          } else {
            /* 'label_tets:106' else */
            /* 'label_tets:107' TetComp(next,ii) = TetComp(Opp,Oppfs); */
            TetComp->data[next + TetComp->size[0] * ii] = TetComp->data[(k +
              TetComp->size[0] * j) - 1];

            /* 'label_tets:108' L(ii)=TetSign(Opp); */
            SV[ii] = (real_T)TetSign->data[k - 1];
          }
        }
      }

      /*  Assign sign */
      /* 'label_tets:114' L = nonzeros(L); */
      j = 0;
      for (k = 0; k < 4; k++) {
        if (SV[k] != 0.0) {
          j++;
        }
      }

      k = L->size[0];
      L->size[0] = j;
      emxEnsureCapacity((emxArray__common *)L, k, (int32_T)sizeof(real_T));
      j = -1;
      for (k = 0; k < 4; k++) {
        if (SV[k] != 0.0) {
          j++;
          L->data[j] = SV[k];
        }
      }

      /* 'label_tets:115' SIGN = unique(L'); */
      k = b_L->size[0] * b_L->size[1];
      b_L->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)b_L, k, (int32_T)sizeof(real_T));
      j = L->size[0];
      k = b_L->size[0] * b_L->size[1];
      b_L->size[1] = j;
      emxEnsureCapacity((emxArray__common *)b_L, k, (int32_T)sizeof(real_T));
      j = L->size[0] - 1;
      for (k = 0; k <= j; k++) {
        b_L->data[k] = L->data[k];
      }

      d_unique(b_L, SIGN);

      /* 'label_tets:116' if (size(SIGN,1)==1) */
      /* 'label_tets:117' TetSign(next) = int32(SIGN(1,1)); */
      b_SV = SIGN->data[0];
      b_SV = b_SV < 0.0 ? ceil(b_SV - 0.5) : floor(b_SV + 0.5);
      TetSign->data[next] = (int32_T)b_SV;

      /*  Assign label */
      /* 'label_tets:123' T = nonzeros(TetComp(next,1:4)); */
      j = 0;
      for (k = 0; k < 4; k++) {
        if (TetComp->data[next + TetComp->size[0] * k] != 0) {
          j++;
        }
      }

      k = T->size[0];
      T->size[0] = j;
      emxEnsureCapacity((emxArray__common *)T, k, (int32_T)sizeof(int32_T));
      j = -1;
      for (k = 0; k < 4; k++) {
        if (TetComp->data[next + TetComp->size[0] * k] != 0) {
          j++;
          T->data[j] = TetComp->data[next + TetComp->size[0] * k];
        }
      }

      /* 'label_tets:124' comp = unique(T'); */
      k = b_T->size[0] * b_T->size[1];
      b_T->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)b_T, k, (int32_T)sizeof(int32_T));
      j = T->size[0];
      k = b_T->size[0] * b_T->size[1];
      b_T->size[1] = j;
      emxEnsureCapacity((emxArray__common *)b_T, k, (int32_T)sizeof(int32_T));
      j = T->size[0] - 1;
      for (k = 0; k <= j; k++) {
        b_T->data[k] = T->data[k];
      }

      e_unique(b_T, comp);

      /* 'label_tets:125' mincomp = min(comp); */
      j = comp->size[1];
      ii = comp->data[0];
      if (j > 1) {
        for (k = 1; k + 1 <= j; k++) {
          if (comp->data[k] < ii) {
            ii = comp->data[k];
          }
        }
      }

      /* 'label_tets:126' TetComp(next,:) = mincomp; */
      for (k = 0; k < 5; k++) {
        TetComp->data[next + TetComp->size[0] * k] = ii;
      }

      /* 'label_tets:127' tetflag(next) = 1; */
      tetflag->data[next] = 1;

      /*  Get mapping */
      /* 'label_tets:130' if (SIGN(1,1) == 1) */
      if (SIGN->data[0] == 1.0) {
        /* 'label_tets:131' [pos_map] = update_map_component(comp,num_comps,comp_list,pos_map); */
        update_map_component(comp, num_comps, comp_list, pos_map);
      } else {
        if (SIGN->data[0] == -1.0) {
          /* 'label_tets:132' elseif (SIGN(1,1) == -1) */
          /* 'label_tets:133' [neg_map] = update_map_component(comp,num_comps,comp_list,neg_map); */
          update_map_component(comp, num_comps, comp_list, neg_map);
        }
      }
    }

    /* 'label_tets:136' if (counter == tnum) */
  } while (!(counter == tnum));

  emxFree_real_T(&b_L);
  emxFree_int32_T(&b_T);
  emxFree_real_T(&L);
  emxFree_int32_T(&comp);
  emxFree_int32_T(&T);
  emxFree_real_T(&SIGN);
  emxFree_int32_T(&opphfs);
  emxFree_int32_T(&stack);
  emxFree_int32_T(&tetflag);
}

/*
 *
 */
static real_T nnz(const real_T s[4])
{
  int32_T j;
  int32_T k;
  j = 0;
  for (k = 0; k < 4; k++) {
    if (s[k] != 0.0) {
      j++;
    }
  }

  return (real_T)j;
}

/*
 *
 */
static void nonzeros(const emxArray_int32_T *s, emxArray_int32_T *v)
{
  int32_T n;
  int32_T j;
  int32_T k;
  int32_T i;
  n = s->size[0];
  j = 0;
  for (k = 0; k <= s->size[0] - 1; k++) {
    if (s->data[k] != 0) {
      j++;
    }
  }

  i = v->size[0];
  v->size[0] = j;
  emxEnsureCapacity((emxArray__common *)v, i, (int32_T)sizeof(int32_T));
  i = -1;
  for (k = 0; k + 1 <= n; k++) {
    if (s->data[k] != 0) {
      i++;
      v->data[i] = s->data[k];
    }
  }
}

/*
 * function [tris_index ,tris_labelled,tris_comp]= re_triangulate(ps,tris_cell,opphes,tris_idx,comp_index_cell,ps_local,ps_localID)
 *  This function retriangulates a given list of triangles (tris_cell) using
 *  the local list of points "ps_local"
 */
static int32_T re_triangulate(const emxArray_real_T *ps, const emxArray_int32_T *
  tris_cell, const emxArray_int32_T *opphes, int32_T tris_idx, const
  emxArray_int32_T *comp_index_cell, const emxArray_real_T *ps_local, const
  emxArray_int32_T *ps_localID, emxArray_int32_T *tris_labelled,
  emxArray_int32_T *tris_comp)
{
  int32_T tris_index;
  int32_T idx;
  int32_T i4;
  int32_T ii;
  emxArray_real_T *pntlist;
  emxArray_real_T *indls;
  emxArray_real_T *cyl;
  int32_T b_tris_cell;
  real_T p[3];
  real_T E2[3];
  real_T nrm[3];
  real_T count;
  uint32_T b_opphes[3];
  uint32_T uv0[3];
  int32_T ngbtriID[3];
  int32_T c_opphes[3];
  int32_T ngbtriEdgeID[3];
  real_T b_pntlist[10];
  int32_T i;
  int32_T pidx;
  int32_T indls_idx;
  int32_T cyl_idx;
  int32_T j;
  real_T b_p[3];
  int32_T c_pntlist[2];
  emxArray_real_T d_pntlist;
  int32_T e_pntlist[2];
  int32_T b_tris_comp[2];
  emxArray_int32_T c_tris_comp;

  /* 're_triangulate:4' factor = int32(5); */
  /* 're_triangulate:5' tris_labelled = zeros(tris_idx*factor,3,'int32'); */
  idx = tris_idx * 5;
  i4 = tris_labelled->size[0] * tris_labelled->size[1];
  tris_labelled->size[0] = idx;
  tris_labelled->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)tris_labelled, i4, (int32_T)sizeof
                    (int32_T));
  idx = idx * 3 - 1;
  for (i4 = 0; i4 <= idx; i4++) {
    tris_labelled->data[i4] = 0;
  }

  /* 're_triangulate:5' trilocal_idx = int32(1); */
  tris_index = 0;

  /* 're_triangulate:6' tris_comp = zeros(tris_idx*factor,1,'int32'); */
  idx = tris_idx * 5;
  i4 = tris_comp->size[0];
  tris_comp->size[0] = idx;
  emxEnsureCapacity((emxArray__common *)tris_comp, i4, (int32_T)sizeof(int32_T));
  idx--;
  for (i4 = 0; i4 <= idx; i4++) {
    tris_comp->data[i4] = 0;
  }

  /* 're_triangulate:7' for ii=1:tris_idx */
  ii = 0;
  b_emxInit_real_T(&pntlist, 1);
  emxInit_real_T(&indls, 2);
  b_emxInit_real_T(&cyl, 1);
  while (ii + 1 <= tris_idx) {
    /*  Triangle Normal */
    /* 're_triangulate:9' nrm = cross_col(ps(tris_cell(ii,2),:)-ps(tris_cell(ii,1),:),ps(tris_cell(ii,3),:)-ps(tris_cell(ii,1),:)); */
    idx = tris_cell->data[ii + tris_cell->size[0]];
    b_tris_cell = tris_cell->data[ii];
    for (i4 = 0; i4 < 3; i4++) {
      p[i4] = ps->data[(idx + ps->size[0] * i4) - 1] - ps->data[(b_tris_cell +
        ps->size[0] * i4) - 1];
    }

    idx = tris_cell->data[ii + (tris_cell->size[0] << 1)];
    b_tris_cell = tris_cell->data[ii];
    for (i4 = 0; i4 < 3; i4++) {
      E2[i4] = ps->data[(idx + ps->size[0] * i4) - 1] - ps->data[(b_tris_cell +
        ps->size[0] * i4) - 1];
    }

    /* CROSS_COL Efficient routine for computing cross product of two  */
    /* 3-dimensional column vectors. */
    /*  CROSS_COL(A,B) Efficiently computes the cross product between */
    /*  3-dimensional column vector A, and 3-dimensional column vector B. */
    /* 'cross_col:7' c = [a(2)*b(3)-a(3)*b(2); a(3)*b(1)-a(1)*b(3); a(1)*b(2)-a(2)*b(1)]; */
    nrm[0] = p[1] * E2[2] - p[2] * E2[1];
    nrm[1] = p[2] * E2[0] - p[0] * E2[2];
    nrm[2] = p[0] * E2[1] - p[1] * E2[0];

    /* 're_triangulate:10' nrm = nrm/sqrt(nrm'*nrm); */
    count = sqrt(eml_xdot(3, nrm, 1, 1, nrm, 1, 1));
    for (i4 = 0; i4 < 3; i4++) {
      nrm[i4] /= count;
    }

    /*  Get triangle neighbors */
    /* 're_triangulate:13' ngbtriID = int32( bitshift( uint32(opphes(ii,:)),-2)); */
    for (i4 = 0; i4 < 3; i4++) {
      b_opphes[i4] = (uint32_T)opphes->data[ii + opphes->size[0] * i4];
    }

    bitshift(b_opphes, -2.0, uv0);
    for (i4 = 0; i4 < 3; i4++) {
      ngbtriID[i4] = (int32_T)uv0[i4];
    }

    /* 're_triangulate:14' ngbtriEdgeID = mod(opphes(ii,:),4)+1; */
    for (i4 = 0; i4 < 3; i4++) {
      c_opphes[i4] = opphes->data[ii + opphes->size[0] * i4];
    }

    b_mod(c_opphes, 4.0, ngbtriEdgeID);
    for (i4 = 0; i4 < 3; i4++) {
      ngbtriEdgeID[i4]++;
    }

    /*  Collect indices of all points (intersections/vertices) for this triangle */
    /* 're_triangulate:17' pntlist = zeros(10,1); */
    for (i = 0; i < 10; i++) {
      b_pntlist[i] = 0.0;
    }

    /* 're_triangulate:17' pidx = int32(1); */
    pidx = 0;

    /* 're_triangulate:18' for jj=1:size(ps_local,1) */
    for (b_tris_cell = 0; b_tris_cell <= ps_local->size[0] - 1; b_tris_cell++) {
      /* 're_triangulate:19' for kk=1:3 */
      for (idx = 0; idx < 3; idx++) {
        /* 're_triangulate:20' if (ps_localID(jj,1) == tris_cell(ii,kk)) */
        if (ps_localID->data[b_tris_cell] == tris_cell->data[ii +
            tris_cell->size[0] * idx]) {
          /* 're_triangulate:21' pntlist(pidx) = jj; */
          b_pntlist[pidx] = 1.0 + (real_T)b_tris_cell;

          /* 're_triangulate:22' pidx = pidx + 1; */
          pidx++;
        }
      }

      /* 're_triangulate:25' if (ps_localID(jj,2) == ii) */
      if (ps_localID->data[b_tris_cell + ps_localID->size[0]] == ii + 1) {
        /* 're_triangulate:26' pntlist(pidx) = jj; */
        b_pntlist[pidx] = 1.0 + (real_T)b_tris_cell;

        /* 're_triangulate:27' pidx = pidx + 1; */
        pidx++;
      }

      /* 're_triangulate:29' for i=1:3 */
      for (i = 0; i < 3; i++) {
        /* 're_triangulate:30' if (ps_localID(jj,2) ==ngbtriID(i)) */
        if ((ps_localID->data[b_tris_cell + ps_localID->size[0]] == ngbtriID[i])
            && (ps_localID->data[b_tris_cell + (ps_localID->size[0] << 1)] ==
                ngbtriEdgeID[i])) {
          /* 're_triangulate:31' if (ps_localID(jj,3)==ngbtriEdgeID(i)) */
          /* 're_triangulate:32' pntlist(pidx) = jj; */
          b_pntlist[pidx] = 1.0 + (real_T)b_tris_cell;

          /* 're_triangulate:33' pidx = pidx + 1; */
          pidx++;
        }
      }
    }

    /* 're_triangulate:39' np = pidx-1; */
    /* 're_triangulate:40' if (np > 2) */
    if (pidx > 2) {
      /* 're_triangulate:41' pntlist = orient_pnts(np,pntlist,ps_local,nrm); */
      /* This function orients a given set of points along a given direction. It is */
      /* assumed that the set of points all lie on a plane and 'dir' is the normal */
      /* of this plane. */
      /* 'orient_pnts:5' ps_new = zeros(np,1); */
      i4 = pntlist->size[0];
      pntlist->size[0] = pidx;
      emxEnsureCapacity((emxArray__common *)pntlist, i4, (int32_T)sizeof(real_T));
      idx = pidx - 1;
      for (i4 = 0; i4 <= idx; i4++) {
        pntlist->data[i4] = 0.0;
      }

      /* 'orient_pnts:6' indls = zeros(np,2); */
      i4 = indls->size[0] * indls->size[1];
      indls->size[0] = pidx;
      indls->size[1] = 2;
      emxEnsureCapacity((emxArray__common *)indls, i4, (int32_T)sizeof(real_T));
      idx = (pidx << 1) - 1;
      for (i4 = 0; i4 <= idx; i4++) {
        indls->data[i4] = 0.0;
      }

      /* 'orient_pnts:6' indls_idx = 1; */
      indls_idx = 0;

      /* 'orient_pnts:7' cyl = zeros(np,1); */
      i4 = cyl->size[0];
      cyl->size[0] = pidx;
      emxEnsureCapacity((emxArray__common *)cyl, i4, (int32_T)sizeof(real_T));
      idx = pidx - 1;
      for (i4 = 0; i4 <= idx; i4++) {
        cyl->data[i4] = 0.0;
      }

      /* 'orient_pnts:7' cyl_idx = 1; */
      cyl_idx = 0;

      /* 'orient_pnts:8' for i=1:np */
      for (i = 1; i <= pidx; i++) {
        /* 'orient_pnts:9' for j=1:np */
        for (j = 1; j <= pidx; j++) {
          /* 'orient_pnts:10' if (j ~= i) */
          if (j != i) {
            /* 'orient_pnts:11' E1 = ps(pntlist(j),:)- ps(pntlist(i),:); */
            for (i4 = 0; i4 < 3; i4++) {
              p[i4] = ps_local->data[((int32_T)b_pntlist[j - 1] + ps_local->
                size[0] * i4) - 1] - ps_local->data[((int32_T)b_pntlist[i - 1] +
                ps_local->size[0] * i4) - 1];
            }

            /* 'orient_pnts:12' count = 0; */
            count = 0.0;

            /* 'orient_pnts:13' for k=1:np */
            for (b_tris_cell = 1; b_tris_cell <= pidx; b_tris_cell++) {
              /* 'orient_pnts:14' if ((k ~= j)&&(k ~= i)) */
              if ((b_tris_cell != j) && (b_tris_cell != i)) {
                /* 'orient_pnts:15' E2 = ps(pntlist(k),:)-ps(pntlist(i),:); */
                for (i4 = 0; i4 < 3; i4++) {
                  E2[i4] = ps_local->data[((int32_T)b_pntlist[b_tris_cell - 1] +
                    ps_local->size[0] * i4) - 1] - ps_local->data[((int32_T)
                    b_pntlist[i - 1] + ps_local->size[0] * i4) - 1];
                }

                /* 'orient_pnts:16' pdir = cross(E1,E2); */
                /* 'orient_pnts:17' sign = dot(pdir,dir); */
                b_p[0] = p[1] * E2[2] - p[2] * E2[1];
                b_p[1] = p[2] * E2[0] - p[0] * E2[2];
                b_p[2] = p[0] * E2[1] - p[1] * E2[0];

                /* 'orient_pnts:18' if (sign >= 0) */
                if (eml_xdot(3, b_p, 1, 1, nrm, 1, 1) >= 0.0) {
                  /* 'orient_pnts:19' count = count + 1; */
                  count++;
                }
              }
            }

            /* 'orient_pnts:23' if (count == np-2) */
            if (count == (real_T)pidx - 2.0) {
              /* 'orient_pnts:24' indls(indls_idx,:) = [i j]; */
              indls->data[indls_idx] = (real_T)(int8_T)i;
              indls->data[indls_idx + indls->size[0]] = (real_T)(int8_T)j;

              /* 'orient_pnts:25' indls_idx = indls_idx + 1; */
              indls_idx++;
            }
          }
        }
      }

      /*  Extract a cycle from the indls */
      /* 'orient_pnts:31' i = int32(1); */
      i = 0;

      /* 'orient_pnts:31' idx = int32(0); */
      idx = -1;

      /* 'orient_pnts:32' while (cyl_idx <= np) */
      while (cyl_idx + 1 <= pidx) {
        /* 'orient_pnts:33' cyl(cyl_idx) = indls(i,1); */
        cyl->data[cyl_idx] = indls->data[i];

        /* 'orient_pnts:34' cyl_idx = cyl_idx + 1; */
        cyl_idx++;

        /* 'orient_pnts:35' for j=1:np */
        for (j = 1; j <= pidx; j++) {
          /* 'orient_pnts:36' if (j~=i) */
          if ((j != i + 1) && (indls->data[j - 1] == indls->data[i + indls->
                               size[0]])) {
            /* 'orient_pnts:37' if (indls(j,1) == indls(i,2)) */
            /* 'orient_pnts:38' idx = j; */
            idx = j - 1;
          }
        }

        /* 'orient_pnts:42' i = idx; */
        i = idx;
      }

      /*  Now put the points in order in the new list */
      /* 'orient_pnts:45' for i=1:np */
      for (i = 0; i + 1 <= pidx; i++) {
        /* 'orient_pnts:46' ps_new(i,:)=pntlist(cyl(i),:); */
        c_pntlist[0] = pntlist->size[0];
        c_pntlist[1] = 1;
        d_pntlist = *pntlist;
        d_pntlist.size = (int32_T *)&c_pntlist;
        d_pntlist.numDimensions = 1;
        d_pntlist.data[i] = b_pntlist[(int32_T)cyl->data[i] - 1];
      }

      /* 're_triangulate:42' v = pntlist(1,:); */
      /* 're_triangulate:43' for i=2:(np-1) */
      for (i = 2; i <= pidx - 1; i++) {
        /* 're_triangulate:44' p = [v pntlist(i) pntlist(i+1)]; */
        e_pntlist[0] = pntlist->size[0];
        e_pntlist[1] = 1;
        d_pntlist = *pntlist;
        d_pntlist.size = (int32_T *)&e_pntlist;
        d_pntlist.numDimensions = 1;
        p[0] = d_pntlist.data[0];
        p[1] = pntlist->data[i - 1];
        p[2] = pntlist->data[i];

        /* 're_triangulate:45' tris_labelled(trilocal_idx,:) = p; */
        for (i4 = 0; i4 < 3; i4++) {
          count = p[i4];
          count = count < 0.0 ? ceil(count - 0.5) : floor(count + 0.5);
          tris_labelled->data[tris_index + tris_labelled->size[0] * i4] =
            (int32_T)count;
        }

        /* 're_triangulate:46' tris_comp(trilocal_idx,:) = comp_index_cell(ii); */
        b_tris_comp[0] = tris_comp->size[0];
        b_tris_comp[1] = 1;
        c_tris_comp = *tris_comp;
        c_tris_comp.size = (int32_T *)&b_tris_comp;
        c_tris_comp.numDimensions = 1;
        c_tris_comp.data[tris_index] = comp_index_cell->data[ii];

        /* 're_triangulate:47' trilocal_idx = trilocal_idx + 1; */
        tris_index++;
      }
    }

    ii++;
  }

  emxFree_real_T(&cyl);
  emxFree_real_T(&indls);
  emxFree_real_T(&pntlist);

  /* 're_triangulate:51' tris_index = trilocal_idx - 1; */
  /*  function visualize(ps_local,tid,tris_local,tidx_before,tidx_after) */
  /*   */
  /*  trisurf(tris_local(tidx_before:tidx_after,:),ps_local(:,1),ps_local(:,2),ps_local(:,3)) */
  /*  % mid = (ps(tris(ii,1),:) +ps(tris(ii,2),:)+ps(tris_local(ii,3),:))/3;  */
  /*   %text(mid(1),mid(2),mid(3),['\leftarrow',num2str(tid)]); */
  /*  for i=tidx_before:tidx_after */
  /*      text(ps_local(tris_local(i,1),1),ps_local(tris_local(i,1),2),ps_local(tris_local(i,1),3),['\leftarrow',num2str(tris_local(i,1))]); */
  /*      text(ps_local(tris_local(i,2),1),ps_local(tris_local(i,2),2),ps_local(tris_local(i,2),3),['\leftarrow',num2str(tris_local(i,2))]); */
  /*      text(ps_local(tris_local(i,3),1),ps_local(tris_local(i,3),2),ps_local(tris_local(i,3),3),['\leftarrow',num2str(tris_local(i,3))]); */
  /*  end */
  /*  end */
  /*  function vis_tris(ps,tris,ps_local,pntlist) */
  /*  figure */
  /*  trisurf(tris,ps(:,1),ps(:,2),ps(:,3)); */
  /*  for i=1:size(pntlist,1) */
  /*      text(ps_local(pntlist(i),1),ps_local(pntlist(i),2),ps_local(pntlist(i),3),['\leftarrow',num2str(pntlist(i))]); */
  /*  end */
  /*  end */
  return tris_index;
}

/*
 * function [L,p] = searchstar(key,idx,num_comps,comp_list,map,L,p)
 */
static void searchstar(int32_T key, int32_T idx, int32_T num_comps, const
  emxArray_int32_T *comp_list, const emxArray_int32_T *map, emxArray_real_T *L,
  real_T *p)
{
  int32_T jj;
  real_T count;
  int32_T i;

  /* 'update_map_component:73' for jj=1:num_comps */
  for (jj = 0; jj + 1 <= num_comps; jj++) {
    /* 'update_map_component:74' if (jj ~= idx) */
    if ((jj + 1 != idx) && (map->data[jj] == key)) {
      /* 'update_map_component:75' if (map(jj) == key) */
      /* 'update_map_component:76' count = 0; */
      count = 0.0;

      /* 'update_map_component:77' for i=1:num_comps */
      i = 1;
      while ((i <= num_comps) && (!(L->data[i - 1] == (real_T)comp_list->data[jj])))
      {
        /* 'update_map_component:78' if (L(i) == comp_list(jj)) */
        /* 'update_map_component:80' else */
        /* 'update_map_component:81' count = count+1; */
        count++;
        i++;
      }

      /* 'update_map_component:84' if (count == num_comps) */
      if (count == (real_T)num_comps) {
        /* 'update_map_component:85' L(p) = comp_list(jj); */
        L->data[(int32_T)*p - 1] = (real_T)comp_list->data[jj];

        /* 'update_map_component:86' p = p+1; */
        (*p)++;
      }
    }
  }
}

/*
 * function[flag_tri,nv,ps_vol,tets,tris_labelled,tris_comp,tris_local,tris_bndID,ps_localID] = tetrahedralize_cell( cell_bnds,...
 *     ps,num_points,tris_cell,tris_idx,comp_index_cell,flag_tri)
 *  This function performs the tetrahedralization of a mesh block with the
 *  embedded surface mesh as a constrain.
 *  Outputs :
 *  tets : kx4 matrix,
 *  ps_vol : positions, nv : no of points,
 *  tris_labelled : Surface triangles
 *  tris_comp: their corresponding components,
 *  tris_local: tris_labelled + all triangulated cell faces
 *  tris_bndID : identification of cell face triangles
 */
static void tetrahedralize_cell(const real_T cell_bnds[6], const emxArray_real_T
  *ps, int32_T num_points, const emxArray_int32_T *tris_cell, int32_T tris_idx,
  const emxArray_int32_T *comp_index_cell, emxArray_real_T *ps_vol,
  emxArray_int32_T *tets, emxArray_int32_T *tris_labelled, emxArray_int32_T
  *tris_comp, emxArray_int32_T *tris_local, emxArray_int32_T *tris_bndID,
  emxArray_int32_T *ps_localID, real_T *flag_tri, int32_T *nv)
{
  emxArray_real_T *ps_local;
  int32_T Face_equID[6];
  real_T Face_equ[6];
  real_T Face_nrm[18];
  real_T corners[24];
  int32_T nb;
  int32_T i2;
  emxArray_int32_T *opphes;
  int32_T ps_idx;
  real_T tri[9];
  int32_T nrows;
  int32_T iv6[3];
  int32_T b_tris_cell[3];
  emxArray_int32_T *r0;
  int32_T i;
  emxArray_int32_T *idx;
  emxArray_int32_T *r1;
  int32_T j;
  emxArray_boolean_T *b;
  int32_T k;
  emxArray_int32_T *b_ps_localID;
  emxArray_int32_T *b_opphes;
  emxArray_int32_T *cor_map;
  int32_T tris_index;
  int8_T b_cor_map[8];
  boolean_T exitg1;
  emxArray_real_T *b_ps_local;
  emxArray_int32_T *c_opphes;
  emxInit_real_T(&ps_local, 2);
  *flag_tri = 0.0;

  /* [Faces,corners,Face_nrm,~,Face_equ,Face_equID] = cell_detail(c_centers,grid_spacing); */
  /* 'tetrahedralize_cell:14' [Faces,corners,Face_nrm,Face_equ,Face_equID] =  gridcell_data(cell_bnds); */
  gridcell_data(cell_bnds, corners, Face_nrm, Face_equ, Face_equID);

  /* 'tetrahedralize_cell:15' factor = int32(6); */
  /* 'tetrahedralize_cell:16' ps_local = zeros(tris_idx*factor+8,3); */
  nb = tris_idx * 6 + 8;
  i2 = ps_local->size[0] * ps_local->size[1];
  ps_local->size[0] = nb;
  ps_local->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)ps_local, i2, (int32_T)sizeof(real_T));
  nb = nb * 3 - 1;
  for (i2 = 0; i2 <= nb; i2++) {
    ps_local->data[i2] = 0.0;
  }

  /* 'tetrahedralize_cell:17' ps_localID = zeros(tris_idx*factor+8,6,'int32'); */
  nb = tris_idx * 6 + 8;
  i2 = ps_localID->size[0] * ps_localID->size[1];
  ps_localID->size[0] = nb;
  ps_localID->size[1] = 6;
  emxEnsureCapacity((emxArray__common *)ps_localID, i2, (int32_T)sizeof(int32_T));
  nb = nb * 6 - 1;
  for (i2 = 0; i2 <= nb; i2++) {
    ps_localID->data[i2] = 0;
  }

  emxInit_int32_T(&opphes, 2);

  /* 'tetrahedralize_cell:18' ps_idx = int32(1); */
  ps_idx = 1;

  /* 'tetrahedralize_cell:19' opphes = determine_opposite_halfedge(num_points, tris_cell); */
  /* DETERMINE_OPPOSITE_HALFEDGE determines the opposite half-edge of  */
  /*  each halfedge for an oriented, manifold surface mesh with or */
  /*  without boundary. It works for both triangle and quadrilateral  */
  /*  meshes that are either linear and quadratic. */
  /*  */
  /*  OPPHES = DETERMINE_OPPOSITE_HALFEDGE(NV,ELEMS) */
  /*  OPPHES = DETERMINE_OPPOSITE_HALFEDGE(NV,ELEMS,OPPHES) */
  /*  computes mapping from each half-edge to its opposite half-edge. This  */
  /*  function supports triangular, quadrilateral, and mixed meshes. */
  /*  */
  /*  Convention: Each half-edge is indicated by <face_id,local_edge_id>. */
  /*     We assign 2 bits to local_edge_id. */
  /*  */
  /*  See also DETERMINE_NEXTPAGE_SURF, DETERMINE_INCIDENT_HALFEDGES */
  /* 'determine_opposite_halfedge:18' if nargin<3 */
  /* 'determine_opposite_halfedge:19' switch size(elems,2) */
  /* 'determine_opposite_halfedge:20' case {3,6} % tri */
  /*  tri */
  /* 'determine_opposite_halfedge:21' opphes = determine_opposite_halfedge_tri(nv, elems); */
  determine_opposite_halfedge_tri(num_points, tris_cell, opphes);

  /* 'tetrahedralize_cell:20' tri = zeros(3,3); */
  for (i2 = 0; i2 < 9; i2++) {
    tri[i2] = 0.0;
  }

  /* LOCAL POINT LIST */
  /* 'tetrahedralize_cell:23' for ii=1:tris_idx */
  for (nrows = 0; nrows + 1 <= tris_idx; nrows++) {
    /* 'tetrahedralize_cell:24' tri_ID = int32(ii); */
    /* 'tetrahedralize_cell:25' tri(1,:) = ps(tris_cell(ii,1),:); */
    nb = tris_cell->data[nrows];
    for (i2 = 0; i2 < 3; i2++) {
      tri[3 * i2] = ps->data[(nb + ps->size[0] * i2) - 1];
    }

    /* 'tetrahedralize_cell:26' tri(2,:) = ps(tris_cell(ii,2),:); */
    nb = tris_cell->data[nrows + tris_cell->size[0]];
    for (i2 = 0; i2 < 3; i2++) {
      tri[1 + 3 * i2] = ps->data[(nb + ps->size[0] * i2) - 1];
    }

    /* 'tetrahedralize_cell:27' tri(3,:) = ps(tris_cell(ii,3),:); */
    nb = tris_cell->data[nrows + (tris_cell->size[0] << 1)];
    for (i2 = 0; i2 < 3; i2++) {
      tri[2 + 3 * i2] = ps->data[(nb + ps->size[0] * i2) - 1];
    }

    /* 'tetrahedralize_cell:28' flag_v = zeros(1,3,'int32'); */
    /* 'tetrahedralize_cell:29' [flag_v] = chk_pnt_inside_MB(cell_bnds,tri,flag_v); */
    /* 'tetrahedralize_cell:30' [ps_local,ps_localID,ps_idx]=  get_pntlist(cell_bnds,opphes,tri_ID,tris_cell(ii,:),tri,flag_v,ps_local,ps_localID,ps_idx); */
    for (i2 = 0; i2 < 3; i2++) {
      iv6[i2] = 0;
    }

    chk_pnt_inside_MB(cell_bnds, tri, iv6);
    for (i2 = 0; i2 < 3; i2++) {
      b_tris_cell[i2] = tris_cell->data[nrows + tris_cell->size[0] * i2];
    }

    get_pntlist(cell_bnds, opphes, nrows + 1, b_tris_cell, tri, iv6, ps_local,
                ps_localID, &ps_idx);
  }

  b_emxInit_int32_T(&r0, 1);

  /* 'tetrahedralize_cell:33' ps_local(ps_idx:end,:) = []; */
  i2 = ps_local->size[0];
  i = r0->size[0];
  r0->size[0] = (i2 - ps_idx) + 1;
  emxEnsureCapacity((emxArray__common *)r0, i, (int32_T)sizeof(int32_T));
  nb = i2 - ps_idx;
  for (i2 = 0; i2 <= nb; i2++) {
    r0->data[i2] = ps_idx + i2;
  }

  emxInit_int32_T(&idx, 2);
  i2 = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)idx, i2, (int32_T)sizeof(int32_T));
  nb = r0->size[0];
  i2 = idx->size[0] * idx->size[1];
  idx->size[1] = nb;
  emxEnsureCapacity((emxArray__common *)idx, i2, (int32_T)sizeof(int32_T));
  nb = r0->size[0] - 1;
  for (i2 = 0; i2 <= nb; i2++) {
    idx->data[i2] = r0->data[i2];
  }

  emxFree_int32_T(&r0);
  b_emxInit_int32_T(&r1, 1);
  c_eml_null_assignment(ps_local, idx);

  /* 'tetrahedralize_cell:34' ps_localID(ps_idx:end,:) = []; */
  i2 = ps_localID->size[0];
  i = r1->size[0];
  r1->size[0] = (i2 - ps_idx) + 1;
  emxEnsureCapacity((emxArray__common *)r1, i, (int32_T)sizeof(int32_T));
  nb = i2 - ps_idx;
  for (i2 = 0; i2 <= nb; i2++) {
    r1->data[i2] = ps_idx + i2;
  }

  i2 = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)idx, i2, (int32_T)sizeof(int32_T));
  nb = r1->size[0];
  i2 = idx->size[0] * idx->size[1];
  idx->size[1] = nb;
  emxEnsureCapacity((emxArray__common *)idx, i2, (int32_T)sizeof(int32_T));
  nb = r1->size[0] - 1;
  for (i2 = 0; i2 <= nb; i2++) {
    idx->data[i2] = r1->data[i2];
  }

  emxFree_int32_T(&r1);
  if (idx->size[1] == 1) {
    nrows = ps_localID->size[0] - 1;
    for (j = 0; j < 6; j++) {
      for (i = idx->data[0]; i <= nrows; i++) {
        ps_localID->data[(i + ps_localID->size[0] * j) - 1] = ps_localID->data[i
          + ps_localID->size[0] * j];
      }
    }
  } else {
    emxInit_boolean_T(&b, 2);
    i2 = b->size[0] * b->size[1];
    b->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)b, i2, (int32_T)sizeof(boolean_T));
    nb = ps_localID->size[0];
    i2 = b->size[0] * b->size[1];
    b->size[1] = nb;
    emxEnsureCapacity((emxArray__common *)b, i2, (int32_T)sizeof(boolean_T));
    nb = ps_localID->size[0] - 1;
    for (i2 = 0; i2 <= nb; i2++) {
      b->data[i2] = FALSE;
    }

    for (k = 1; k <= idx->size[1]; k++) {
      b->data[idx->data[k - 1] - 1] = TRUE;
    }

    nb = 0;
    for (k = 1; k <= b->size[1]; k++) {
      nrows = b->data[k - 1];
      nb += nrows;
    }

    nrows = ps_localID->size[0] - nb;
    nb = b->size[1];
    i = 0;
    i2 = ps_localID->size[0];
    for (k = 1; k <= i2; k++) {
      if ((k > nb) || (!b->data[k - 1])) {
        for (j = 0; j < 6; j++) {
          ps_localID->data[i + ps_localID->size[0] * j] = ps_localID->data[(k +
            ps_localID->size[0] * j) - 1];
        }

        i++;
      }
    }

    emxFree_boolean_T(&b);
  }

  emxFree_int32_T(&idx);
  if (1 > nrows) {
    nrows = 0;
  }

  emxInit_int32_T(&b_ps_localID, 2);
  i2 = b_ps_localID->size[0] * b_ps_localID->size[1];
  b_ps_localID->size[0] = nrows;
  b_ps_localID->size[1] = 6;
  emxEnsureCapacity((emxArray__common *)b_ps_localID, i2, (int32_T)sizeof
                    (int32_T));
  for (i2 = 0; i2 < 6; i2++) {
    nb = nrows - 1;
    for (i = 0; i <= nb; i++) {
      b_ps_localID->data[i + b_ps_localID->size[0] * i2] = ps_localID->data[i +
        ps_localID->size[0] * i2];
    }
  }

  i2 = ps_localID->size[0] * ps_localID->size[1];
  ps_localID->size[0] = b_ps_localID->size[0];
  ps_localID->size[1] = 6;
  emxEnsureCapacity((emxArray__common *)ps_localID, i2, (int32_T)sizeof(int32_T));
  for (i2 = 0; i2 < 6; i2++) {
    nb = b_ps_localID->size[0] - 1;
    for (i = 0; i <= nb; i++) {
      ps_localID->data[i + ps_localID->size[0] * i2] = b_ps_localID->data[i +
        b_ps_localID->size[0] * i2];
    }
  }

  emxFree_int32_T(&b_ps_localID);
  emxInit_int32_T(&b_opphes, 2);

  /* %% DEBUG VISUALIZATION */
  /*  figure */
  /*   hold on */
  /*  elems = [1 2; 2 3; 3 4; 4 1; 1 5; 2 6; 3 7; 4 8; 5 6; 6 7; 7 8; 5 8]; */
  /*  for ii=1:12 */
  /*      x = [corners(elems(ii,1),1);corners(elems(ii,2),1)]; */
  /*      y = [corners(elems(ii,1),2);corners(elems(ii,2),2)]; */
  /*      z = [corners(elems(ii,1),3);corners(elems(ii,2),3)]; */
  /*      plot3(x,y,z); */
  /*  end */
  /*  RETRIANGULATION */
  /* 'tetrahedralize_cell:48' [tris_index,tris_labelled,tris_comp]= re_triangulate(ps,tris_cell,opphes,tris_idx,comp_index_cell,ps_local,ps_localID); */
  i2 = b_opphes->size[0] * b_opphes->size[1];
  b_opphes->size[0] = opphes->size[0];
  b_opphes->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)b_opphes, i2, (int32_T)sizeof(int32_T));
  nb = opphes->size[0] * opphes->size[1] - 1;
  for (i2 = 0; i2 <= nb; i2++) {
    b_opphes->data[i2] = opphes->data[i2];
  }

  b_emxInit_int32_T(&cor_map, 1);
  tris_index = re_triangulate(ps, tris_cell, b_opphes, tris_idx, comp_index_cell,
    ps_local, ps_localID, opphes, cor_map);

  /* 'tetrahedralize_cell:49' tris_labelled = tris_labelled(1:tris_index,:); */
  emxFree_int32_T(&b_opphes);
  if (1 > tris_index) {
    i2 = -1;
  } else {
    i2 = tris_index - 1;
  }

  i = tris_labelled->size[0] * tris_labelled->size[1];
  tris_labelled->size[0] = i2 + 1;
  tris_labelled->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)tris_labelled, i, (int32_T)sizeof
                    (int32_T));
  for (i = 0; i < 3; i++) {
    for (nrows = 0; nrows <= i2; nrows++) {
      tris_labelled->data[nrows + tris_labelled->size[0] * i] = opphes->
        data[nrows + opphes->size[0] * i];
    }
  }

  /* 'tetrahedralize_cell:50' tris_comp = tris_comp(1:tris_index); */
  if (1 > tris_index) {
    i = 0;
  } else {
    i = tris_index;
  }

  nrows = tris_comp->size[0];
  tris_comp->size[0] = i;
  emxEnsureCapacity((emxArray__common *)tris_comp, nrows, (int32_T)sizeof
                    (int32_T));
  nb = i - 1;
  for (i = 0; i <= nb; i++) {
    tris_comp->data[i] = cor_map->data[i];
  }

  /* 'tetrahedralize_cell:52' if (tris_index ~=0) */
  if (tris_index != 0) {
    /* 'tetrahedralize_cell:53' flag_tri = 1; */
    *flag_tri = 1.0;

    /*  TRIANGULATE THE CELL FACES */
    /* 'tetrahedralize_cell:55' MAX_TCF = 100*6; */
    /*  Max no of triangles on all 6 cell faces */
    /* 'tetrahedralize_cell:56' tris_local = zeros(tris_index+MAX_TCF,3,'int32'); */
    i = tris_local->size[0] * tris_local->size[1];
    tris_local->size[0] = tris_index + 600;
    tris_local->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)tris_local, i, (int32_T)sizeof(int32_T));
    nb = (tris_index + 600) * 3 - 1;
    for (i = 0; i <= nb; i++) {
      tris_local->data[i] = 0;
    }

    /* 'tetrahedralize_cell:57' tris_local(1:tris_index,:) = tris_labelled; */
    for (i = 0; i < 3; i++) {
      for (nrows = 0; nrows <= i2; nrows++) {
        tris_local->data[nrows + tris_local->size[0] * i] = opphes->data[nrows +
          opphes->size[0] * i];
      }
    }

    /* 'tetrahedralize_cell:58' tris_bndID = zeros(tris_index+MAX_TCF,1,'int32'); */
    /* 'tetrahedralize_cell:59' ps_index = ps_idx-1; */
    ps_idx--;

    /* 'tetrahedralize_cell:61' cor_map =int32([1; 2; 3; 4; 5; 6; 7; 8]); */
    for (i = 0; i < 8; i++) {
      b_cor_map[i] = (int8_T)(1 + i);
    }

    /* 'tetrahedralize_cell:62' for jj=1:8 */
    for (nb = 0; nb < 8; nb++) {
      /* 'tetrahedralize_cell:63' for ii=1:ps_index */
      nrows = 1;
      exitg1 = 0U;
      while ((exitg1 == 0U) && (nrows <= ps_idx)) {
        /* 'tetrahedralize_cell:64' if (ps_localID(ii,6)==cor_map(jj)) */
        if (ps_localID->data[(nrows + ps_localID->size[0] * 5) - 1] ==
            b_cor_map[nb]) {
          /* 'tetrahedralize_cell:65' cor_map(jj)=0; */
          b_cor_map[nb] = 0;
          exitg1 = 1U;
        } else {
          nrows++;
        }
      }
    }

    /* 'tetrahedralize_cell:69' cor_map = nonzeros(cor_map); */
    j = 0;
    for (k = 0; k < 8; k++) {
      if (b_cor_map[k] != 0) {
        j++;
      }
    }

    i = cor_map->size[0];
    cor_map->size[0] = j;
    emxEnsureCapacity((emxArray__common *)cor_map, i, (int32_T)sizeof(int32_T));
    i = -1;
    for (k = 0; k < 8; k++) {
      if (b_cor_map[k] != 0) {
        i++;
        cor_map->data[i] = (int32_T)b_cor_map[k];
      }
    }

    /* 'tetrahedralize_cell:70' if (size(cor_map,1)==8) */
    if (cor_map->size[0] == 8) {
      emxInit_real_T(&b_ps_local, 2);

      /* 'tetrahedralize_cell:71' ps_local = [ps_local; corners]; */
      i = b_ps_local->size[0] * b_ps_local->size[1];
      b_ps_local->size[0] = ps_local->size[0] + 8;
      b_ps_local->size[1] = 3;
      emxEnsureCapacity((emxArray__common *)b_ps_local, i, (int32_T)sizeof
                        (real_T));
      for (i = 0; i < 3; i++) {
        nb = ps_local->size[0] - 1;
        for (nrows = 0; nrows <= nb; nrows++) {
          b_ps_local->data[nrows + b_ps_local->size[0] * i] = ps_local->
            data[nrows + ps_local->size[0] * i];
        }
      }

      for (i = 0; i < 3; i++) {
        for (nrows = 0; nrows < 8; nrows++) {
          b_ps_local->data[(nrows + ps_local->size[0]) + b_ps_local->size[0] * i]
            = corners[nrows + (i << 3)];
        }
      }

      i = ps_local->size[0] * ps_local->size[1];
      ps_local->size[0] = b_ps_local->size[0];
      ps_local->size[1] = 3;
      emxEnsureCapacity((emxArray__common *)ps_local, i, (int32_T)sizeof(real_T));
      for (i = 0; i < 3; i++) {
        nb = b_ps_local->size[0] - 1;
        for (nrows = 0; nrows <= nb; nrows++) {
          ps_local->data[nrows + ps_local->size[0] * i] = b_ps_local->data[nrows
            + b_ps_local->size[0] * i];
        }
      }

      emxFree_real_T(&b_ps_local);
    } else {
      emxInit_real_T(&b_ps_local, 2);

      /* 'tetrahedralize_cell:72' else */
      /* 'tetrahedralize_cell:73' ps_local = [ps_local; corners(cor_map,:)]; */
      i = b_ps_local->size[0] * b_ps_local->size[1];
      b_ps_local->size[0] = ps_local->size[0] + cor_map->size[0];
      b_ps_local->size[1] = 3;
      emxEnsureCapacity((emxArray__common *)b_ps_local, i, (int32_T)sizeof
                        (real_T));
      for (i = 0; i < 3; i++) {
        nb = ps_local->size[0] - 1;
        for (nrows = 0; nrows <= nb; nrows++) {
          b_ps_local->data[nrows + b_ps_local->size[0] * i] = ps_local->
            data[nrows + ps_local->size[0] * i];
        }
      }

      for (i = 0; i < 3; i++) {
        nb = cor_map->size[0] - 1;
        for (nrows = 0; nrows <= nb; nrows++) {
          b_ps_local->data[(nrows + ps_local->size[0]) + b_ps_local->size[0] * i]
            = corners[(cor_map->data[nrows] + (i << 3)) - 1];
        }
      }

      i = ps_local->size[0] * ps_local->size[1];
      ps_local->size[0] = b_ps_local->size[0];
      ps_local->size[1] = 3;
      emxEnsureCapacity((emxArray__common *)ps_local, i, (int32_T)sizeof(real_T));
      for (i = 0; i < 3; i++) {
        nb = b_ps_local->size[0] - 1;
        for (nrows = 0; nrows <= nb; nrows++) {
          ps_local->data[nrows + ps_local->size[0] * i] = b_ps_local->data[nrows
            + b_ps_local->size[0] * i];
        }
      }

      emxFree_real_T(&b_ps_local);
    }

    /* figure */
    /* 'tetrahedralize_cell:77' [ps_local,tris_local,tris_bndID]=triangulate_cellface(Faces,Face_nrm,Face_equ,Face_equID,... */
    /* 'tetrahedralize_cell:78'         cor_map,ps_index,ps_local,ps_localID,tris_index,tris_local,tris_bndID,tris_labelled); */
    i = tris_bndID->size[0];
    tris_bndID->size[0] = tris_index + 600;
    emxEnsureCapacity((emxArray__common *)tris_bndID, i, (int32_T)sizeof(int32_T));
    nb = tris_index + 599;
    for (i = 0; i <= nb; i++) {
      tris_bndID->data[i] = 0;
    }

    emxInit_int32_T(&c_opphes, 2);
    i = c_opphes->size[0] * c_opphes->size[1];
    c_opphes->size[0] = i2 + 1;
    c_opphes->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)c_opphes, i, (int32_T)sizeof(int32_T));
    for (i = 0; i < 3; i++) {
      for (nrows = 0; nrows <= i2; nrows++) {
        c_opphes->data[nrows + c_opphes->size[0] * i] = opphes->data[nrows +
          opphes->size[0] * i];
      }
    }

    triangulate_cellface(Face_nrm, Face_equ, Face_equID, cor_map, ps_idx,
                         ps_local, ps_localID, tris_index, tris_local,
                         tris_bndID, c_opphes);

    /*    fprintf('Passed Triangulate Cell Face\n'); */
    /* vis_triincell(ps_local,tris_local,1,corners) */
    /*  TETRAHEDRALIZATION */
    /* 'tetrahedralize_cell:84' [ps_vol,tets] = call_tetgen(ps_local, tris_local); */
    call_tetgen(ps_local, tris_local, ps_vol, tets);

    /* fprintf('Passed Tetgen\n'); */
    /* 'tetrahedralize_cell:85' nv = int32(size(ps_vol,1)); */
    *nv = ps_vol->size[0];
    emxFree_int32_T(&c_opphes);
  } else {
    /* 'tetrahedralize_cell:86' else */
    /* 'tetrahedralize_cell:87' nv = int32(0); */
    *nv = 0;

    /* 'tetrahedralize_cell:88' ps_vol = [0 0 0]; */
    i2 = ps_vol->size[0] * ps_vol->size[1];
    ps_vol->size[0] = 1;
    ps_vol->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)ps_vol, i2, (int32_T)sizeof(real_T));
    for (i2 = 0; i2 < 3; i2++) {
      ps_vol->data[i2] = 0.0;
    }

    /* 'tetrahedralize_cell:89' tets = int32([0 0 0 0]); */
    i2 = tets->size[0] * tets->size[1];
    tets->size[0] = 1;
    tets->size[1] = 4;
    emxEnsureCapacity((emxArray__common *)tets, i2, (int32_T)sizeof(int32_T));
    for (i2 = 0; i2 < 4; i2++) {
      tets->data[i2] = 0;
    }

    /* 'tetrahedralize_cell:90' tris_local = int32([0 0 0]); */
    i2 = tris_local->size[0] * tris_local->size[1];
    tris_local->size[0] = 1;
    tris_local->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)tris_local, i2, (int32_T)sizeof
                      (int32_T));
    for (i2 = 0; i2 < 3; i2++) {
      tris_local->data[i2] = 0;
    }

    /* 'tetrahedralize_cell:91' tris_bndID = int32(0); */
    i2 = tris_bndID->size[0];
    tris_bndID->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)tris_bndID, i2, (int32_T)sizeof
                      (int32_T));
    tris_bndID->data[0] = 0;
  }

  emxFree_int32_T(&cor_map);
  emxFree_int32_T(&opphes);
  emxFree_real_T(&ps_local);
}

/*
 * function [ps_local,tris_local,tris_bndID] = triangulate_cellface(Faces,Face_nrm,...
 *     Face_equ,Face_equID,cor_map,ps_index,ps_local,ps_localID,tris_index,tris_local,tris_bndID,tris_labelled)
 *  This function triangulates the cell faces using 'triangle' and adds the
 *  extra vertices and triangles generated to the corresponding local list.
 */
static void triangulate_cellface(const real_T Face_nrm[18], const real_T
  Face_equ[6], const int32_T Face_equID[6], const emxArray_int32_T *cor_map,
  int32_T ps_index, emxArray_real_T *ps_local, const emxArray_int32_T
  *ps_localID, int32_T tris_index, emxArray_int32_T *tris_local,
  emxArray_int32_T *tris_bndID, const emxArray_int32_T *tris_labelled)
{
  emxArray_int32_T *opphes;
  real_T chkdir;
  int32_T ps_idx;
  int32_T tidx;
  emxArray_int32_T *Face_ps;
  emxArray_real_T *Fps;
  emxArray_int32_T *Elems;
  emxArray_int32_T *b_Face_ps;
  emxArray_int32_T *Face_elems;
  emxArray_real_T *Fps_new;
  emxArray_int32_T *Ftris;
  emxArray_int32_T *T;
  emxArray_real_T *b_ps_local;
  emxArray_int32_T *b_tris_local;
  emxArray_int32_T *b_tris_bndID;
  int32_T ii;
  int32_T Faces[4];
  int32_T FaceEdgeMap[4];
  int32_T i17;
  static const int8_T b_Faces[24] = { 1, 1, 2, 3, 1, 5, 4, 2, 3, 4, 5, 6, 3, 6,
    7, 8, 8, 7, 2, 5, 6, 7, 4, 8 };

  static const int8_T b_FaceEdgeMap[24] = { 4, 1, 2, 3, 5, 9, 3, 6, 7, 8, 12, 10,
    2, 9, 10, 11, 8, 11, 1, 5, 6, 7, 4, 12 };

  int32_T loop_ub;
  real_T fps_idx;
  int32_T SV;
  int32_T i;
  real_T numextra;
  static const int8_T Coord_twodim[6] = { 2, 1, 1, 3, 3, 2 };

  int32_T chk;
  int32_T tri[3];
  boolean_T exitg1;
  int32_T b_tri[6];
  real_T a[3];
  real_T b[3];
  real_T nrm[3];
  emxArray_int32_T *c_tris_local;
  emxArray_int32_T *c_tris_bndID;
  emxInit_int32_T(&opphes, 2);

  /* 'triangulate_cellface:6' coder.inline('never') */
  /* 'triangulate_cellface:7' Coord_twodim =int32([2 3; 1 3; 1 2]); */
  /* 'triangulate_cellface:8' FaceEdgeMap = int32([4,3,2,1;1,6,9,5;2,7,10,6;3,8,11,7;5,12,8,4;9,10,11,12]); */
  /* [b2v, bdedgs] = extract_border_curv(ps_index, tris_local); */
  /* 'triangulate_cellface:10' opphes = determine_opposite_halfedge(ps_index, tris_labelled); */
  /* DETERMINE_OPPOSITE_HALFEDGE determines the opposite half-edge of  */
  /*  each halfedge for an oriented, manifold surface mesh with or */
  /*  without boundary. It works for both triangle and quadrilateral  */
  /*  meshes that are either linear and quadratic. */
  /*  */
  /*  OPPHES = DETERMINE_OPPOSITE_HALFEDGE(NV,ELEMS) */
  /*  OPPHES = DETERMINE_OPPOSITE_HALFEDGE(NV,ELEMS,OPPHES) */
  /*  computes mapping from each half-edge to its opposite half-edge. This  */
  /*  function supports triangular, quadrilateral, and mixed meshes. */
  /*  */
  /*  Convention: Each half-edge is indicated by <face_id,local_edge_id>. */
  /*     We assign 2 bits to local_edge_id. */
  /*  */
  /*  See also DETERMINE_NEXTPAGE_SURF, DETERMINE_INCIDENT_HALFEDGES */
  /* 'determine_opposite_halfedge:18' if nargin<3 */
  /* 'determine_opposite_halfedge:19' switch size(elems,2) */
  /* 'determine_opposite_halfedge:20' case {3,6} % tri */
  /*  tri */
  /* 'determine_opposite_halfedge:21' opphes = determine_opposite_halfedge_tri(nv, elems); */
  determine_opposite_halfedge_tri(ps_index, tris_labelled, opphes);

  /* 'triangulate_cellface:11' ps_idx = int32(size(ps_local,1)+1); */
  chkdir = (real_T)ps_local->size[0] + 1.0;
  chkdir = chkdir < 0.0 ? ceil(chkdir - 0.5) : floor(chkdir + 0.5);
  ps_idx = (int32_T)chkdir - 1;

  /* 'triangulate_cellface:12' tidx = tris_index + 1; */
  tidx = tris_index;

  /* 'triangulate_cellface:13' for ii=1:6 */
  b_emxInit_int32_T(&Face_ps, 1);
  emxInit_real_T(&Fps, 2);
  emxInit_int32_T(&Elems, 2);
  b_emxInit_int32_T(&b_Face_ps, 1);
  emxInit_int32_T(&Face_elems, 2);
  emxInit_real_T(&Fps_new, 2);
  emxInit_int32_T(&Ftris, 2);
  emxInit_int32_T(&T, 2);
  emxInit_real_T(&b_ps_local, 2);
  emxInit_int32_T(&b_tris_local, 2);
  b_emxInit_int32_T(&b_tris_bndID, 1);
  for (ii = 0; ii < 6; ii++) {
    /*  EXTRACTION OF PNTS/EDGES ON CELL BOUNDARIES */
    /* 'triangulate_cellface:15' [Face_ps,Face_elems] = extract_cellface_edges(ii,Faces(ii,:),FaceEdgeMap(ii,:),cor_map,ps_index,ps_local,ps_localID,tris_labelled,opphes); */
    for (i17 = 0; i17 < 4; i17++) {
      Faces[i17] = (int32_T)b_Faces[ii + 6 * i17];
      FaceEdgeMap[i17] = (int32_T)b_FaceEdgeMap[ii + 6 * i17];
    }

    extract_cellface_edges(1.0 + (real_T)ii, Faces, FaceEdgeMap, cor_map,
      ps_index, ps_local, ps_localID, tris_labelled, opphes, b_Face_ps,
      Face_elems);
    i17 = Face_ps->size[0];
    Face_ps->size[0] = b_Face_ps->size[0];
    emxEnsureCapacity((emxArray__common *)Face_ps, i17, (int32_T)sizeof(int32_T));
    loop_ub = b_Face_ps->size[0] - 1;
    for (i17 = 0; i17 <= loop_ub; i17++) {
      Face_ps->data[i17] = b_Face_ps->data[i17];
    }

    /*     visualize_FaceEdges(ps_local,Face_elems) */
    /*  TRANSFORMATION TO LOCAL ID */
    /* 'triangulate_cellface:19' numvert = size(Face_ps,1); */
    /* 'triangulate_cellface:19' fps_idx = numvert+1; */
    fps_idx = (real_T)b_Face_ps->size[0] + 1.0;

    /* 'triangulate_cellface:20' Fps = zeros(numvert,2); */
    SV = b_Face_ps->size[0];
    i17 = Fps->size[0] * Fps->size[1];
    Fps->size[0] = SV;
    Fps->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)Fps, i17, (int32_T)sizeof(real_T));
    loop_ub = (b_Face_ps->size[0] << 1) - 1;
    for (i17 = 0; i17 <= loop_ub; i17++) {
      Fps->data[i17] = 0.0;
    }

    /* 'triangulate_cellface:21' Elems = zeros(size(Face_elems,1),2, 'int32'); */
    SV = Face_elems->size[0];
    i17 = Elems->size[0] * Elems->size[1];
    Elems->size[0] = SV;
    Elems->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)Elems, i17, (int32_T)sizeof(int32_T));
    loop_ub = (Face_elems->size[0] << 1) - 1;
    for (i17 = 0; i17 <= loop_ub; i17++) {
      Elems->data[i17] = 0;
    }

    /* 'triangulate_cellface:22' for i=1:size(Face_elems,1) */
    for (i = 0; i <= Face_elems->size[0] - 1; i++) {
      /* 'triangulate_cellface:23' v1 = Face_elems(i,1); */
      /* 'triangulate_cellface:23' v1id = 0; */
      chkdir = 0.0;

      /* 'triangulate_cellface:24' v2 = Face_elems(i,2); */
      /* 'triangulate_cellface:24' v2id = 0; */
      numextra = 0.0;

      /* 'triangulate_cellface:25' for j=1:numvert */
      for (SV = 0; SV <= b_Face_ps->size[0] - 1; SV++) {
        /* 'triangulate_cellface:26' if (v1 == Face_ps(j)) */
        if (Face_elems->data[(int32_T)(1.0 + (real_T)i) - 1] == b_Face_ps->data
            [(int32_T)(1.0 + (real_T)SV) - 1]) {
          /* 'triangulate_cellface:27' v1id = j; */
          chkdir = 1.0 + (real_T)SV;
        }

        /* 'triangulate_cellface:29' if (v2 == Face_ps(j)) */
        if (Face_elems->data[((int32_T)(1.0 + (real_T)i) + Face_elems->size[0])
            - 1] == b_Face_ps->data[(int32_T)(1.0 + (real_T)SV) - 1]) {
          /* 'triangulate_cellface:30' v2id = j; */
          numextra = 1.0 + (real_T)SV;
        }
      }

      /* 'triangulate_cellface:33' Elems(i,1) = v1id; */
      chkdir = chkdir < 0.0 ? ceil(chkdir - 0.5) : floor(chkdir + 0.5);
      Elems->data[(int32_T)(1.0 + (real_T)i) - 1] = (int32_T)chkdir;

      /* 'triangulate_cellface:34' Elems(i,2) = v2id; */
      numextra = numextra < 0.0 ? ceil(numextra - 0.5) : floor(numextra + 0.5);
      Elems->data[((int32_T)(1.0 + (real_T)i) + Elems->size[0]) - 1] = (int32_T)
        numextra;
    }

    /*  TRANSFORMATION TO 2D, 'TRIANGLE', BACK TO 3D */
    /* 'triangulate_cellface:39' kk = Face_equID(ii); */
    /* 'triangulate_cellface:40' for i=1:numvert */
    for (i = 0; i <= b_Face_ps->size[0] - 1; i++) {
      /* 'triangulate_cellface:41' Fps(i,1)= ps_local(Face_ps(i),Coord_twodim(kk,1)); */
      Fps->data[(int32_T)(1.0 + (real_T)i) - 1] = ps_local->data
        [(b_Face_ps->data[(int32_T)(1.0 + (real_T)i) - 1] + ps_local->size[0] *
          (Coord_twodim[Face_equID[ii] - 1] - 1)) - 1];

      /* 'triangulate_cellface:42' Fps(i,2)= ps_local(Face_ps(i),Coord_twodim(kk,2)); */
      Fps->data[((int32_T)(1.0 + (real_T)i) + Fps->size[0]) - 1] =
        ps_local->data[(b_Face_ps->data[(int32_T)(1.0 + (real_T)i) - 1] +
                        ps_local->size[0] * (Coord_twodim[Face_equID[ii] + 2] -
        1)) - 1];
    }

    /* 'triangulate_cellface:45' [Fps_new,Ftris]=call_triangle(Fps,Elems); */
    call_triangle(Fps, Elems, Fps_new, Ftris);

    /* 'triangulate_cellface:46' if (size(Fps_new,1) > size(Fps,1)) */
    if (Fps_new->size[0] > Fps->size[0]) {
      /* 'triangulate_cellface:47' numextra = size(Fps_new,1)-numvert; */
      numextra = (real_T)Fps_new->size[0] - (real_T)b_Face_ps->size[0];

      /* 'triangulate_cellface:48' if size(ps_local,1)<ps_idx+numextra-1 */
      chkdir = (real_T)(ps_idx + 1) + numextra;
      chkdir = chkdir < 0.0 ? ceil(chkdir - 0.5) : floor(chkdir + 0.5);
      if (ps_local->size[0] < (int32_T)chkdir - 1) {
        /* 'triangulate_cellface:49' ps_local = [ps_local; zeros(numextra,3)]; */
        i17 = b_ps_local->size[0] * b_ps_local->size[1];
        b_ps_local->size[0] = ps_local->size[0] + (int32_T)numextra;
        b_ps_local->size[1] = 3;
        emxEnsureCapacity((emxArray__common *)b_ps_local, i17, (int32_T)sizeof
                          (real_T));
        for (i17 = 0; i17 < 3; i17++) {
          loop_ub = ps_local->size[0] - 1;
          for (chk = 0; chk <= loop_ub; chk++) {
            b_ps_local->data[chk + b_ps_local->size[0] * i17] = ps_local->
              data[chk + ps_local->size[0] * i17];
          }
        }

        for (i17 = 0; i17 < 3; i17++) {
          loop_ub = (int32_T)numextra - 1;
          for (chk = 0; chk <= loop_ub; chk++) {
            b_ps_local->data[(chk + ps_local->size[0]) + b_ps_local->size[0] *
              i17] = 0.0;
          }
        }

        i17 = ps_local->size[0] * ps_local->size[1];
        ps_local->size[0] = b_ps_local->size[0];
        ps_local->size[1] = 3;
        emxEnsureCapacity((emxArray__common *)ps_local, i17, (int32_T)sizeof
                          (real_T));
        for (i17 = 0; i17 < 3; i17++) {
          loop_ub = b_ps_local->size[0] - 1;
          for (chk = 0; chk <= loop_ub; chk++) {
            ps_local->data[chk + ps_local->size[0] * i17] = b_ps_local->data[chk
              + b_ps_local->size[0] * i17];
          }
        }

        /* 'triangulate_cellface:50' Face_ps = [Face_ps; zeros(numextra,1)]; */
        i17 = Face_ps->size[0];
        Face_ps->size[0] = b_Face_ps->size[0] + (int32_T)numextra;
        emxEnsureCapacity((emxArray__common *)Face_ps, i17, (int32_T)sizeof
                          (int32_T));
        loop_ub = b_Face_ps->size[0] - 1;
        for (i17 = 0; i17 <= loop_ub; i17++) {
          Face_ps->data[i17] = b_Face_ps->data[i17];
        }

        loop_ub = (int32_T)numextra - 1;
        for (i17 = 0; i17 <= loop_ub; i17++) {
          Face_ps->data[i17 + b_Face_ps->size[0]] = 0;
        }
      }

      /* 'triangulate_cellface:52' for i=1:numextra */
      for (i = 0; i <= (int32_T)numextra - 1; i++) {
        /* 'triangulate_cellface:53' ps_local(ps_idx,kk) = Face_equ(ii); */
        ps_local->data[ps_idx + ps_local->size[0] * (Face_equID[ii] - 1)] =
          Face_equ[ii];

        /* 'triangulate_cellface:54' ps_local(ps_idx,Coord_twodim(kk,1)) = Fps_new(i+numvert,1); */
        ps_local->data[ps_idx + ps_local->size[0] * (Coord_twodim[Face_equID[ii]
          - 1] - 1)] = Fps_new->data[(int32_T)((1.0 + (real_T)i) + (real_T)
          b_Face_ps->size[0]) - 1];

        /* 'triangulate_cellface:55' ps_local(ps_idx,Coord_twodim(kk,2)) = Fps_new(i+numvert,2); */
        ps_local->data[ps_idx + ps_local->size[0] * (Coord_twodim[Face_equID[ii]
          + 2] - 1)] = Fps_new->data[((int32_T)((1.0 + (real_T)i) + (real_T)
          b_Face_ps->size[0]) + Fps_new->size[0]) - 1];

        /* 'triangulate_cellface:56' ps_idx = ps_idx+1; */
        ps_idx++;

        /* 'triangulate_cellface:57' Face_ps(fps_idx) = i+numvert; */
        chkdir = (1.0 + (real_T)i) + (real_T)b_Face_ps->size[0];
        chkdir = chkdir < 0.0 ? ceil(chkdir - 0.5) : floor(chkdir + 0.5);
        Face_ps->data[(int32_T)fps_idx - 1] = (int32_T)chkdir;

        /* 'triangulate_cellface:58' fps_idx = fps_idx + 1; */
        fps_idx++;
      }
    }

    /* 'triangulate_cellface:62' if size(tris_local,1) < tidx+size(Ftris,1)-1 */
    chkdir = (real_T)(tidx + 1) + (real_T)Ftris->size[0];
    chkdir = chkdir < 0.0 ? ceil(chkdir - 0.5) : floor(chkdir + 0.5);
    if (tris_local->size[0] < (int32_T)chkdir - 1) {
      /* 'triangulate_cellface:63' tris_local = [tris_local; zeros(size(Ftris,1),3,'int32')]; */
      SV = Ftris->size[0];
      i17 = b_tris_local->size[0] * b_tris_local->size[1];
      b_tris_local->size[0] = tris_local->size[0] + SV;
      b_tris_local->size[1] = 3;
      emxEnsureCapacity((emxArray__common *)b_tris_local, i17, (int32_T)sizeof
                        (int32_T));
      for (i17 = 0; i17 < 3; i17++) {
        loop_ub = tris_local->size[0] - 1;
        for (chk = 0; chk <= loop_ub; chk++) {
          b_tris_local->data[chk + b_tris_local->size[0] * i17] =
            tris_local->data[chk + tris_local->size[0] * i17];
        }
      }

      for (i17 = 0; i17 < 3; i17++) {
        loop_ub = SV - 1;
        for (chk = 0; chk <= loop_ub; chk++) {
          b_tris_local->data[(chk + tris_local->size[0]) + b_tris_local->size[0]
            * i17] = 0;
        }
      }

      i17 = tris_local->size[0] * tris_local->size[1];
      tris_local->size[0] = b_tris_local->size[0];
      tris_local->size[1] = 3;
      emxEnsureCapacity((emxArray__common *)tris_local, i17, (int32_T)sizeof
                        (int32_T));
      for (i17 = 0; i17 < 3; i17++) {
        loop_ub = b_tris_local->size[0] - 1;
        for (chk = 0; chk <= loop_ub; chk++) {
          tris_local->data[chk + tris_local->size[0] * i17] = b_tris_local->
            data[chk + b_tris_local->size[0] * i17];
        }
      }

      /* 'triangulate_cellface:64' tris_bndID = [tris_bndID; zeros(size(Ftris,1),1,'int32')]; */
      SV = Ftris->size[0];
      i17 = b_tris_bndID->size[0];
      b_tris_bndID->size[0] = tris_bndID->size[0] + SV;
      emxEnsureCapacity((emxArray__common *)b_tris_bndID, i17, (int32_T)sizeof
                        (int32_T));
      loop_ub = tris_bndID->size[0] - 1;
      for (i17 = 0; i17 <= loop_ub; i17++) {
        b_tris_bndID->data[i17] = tris_bndID->data[i17];
      }

      loop_ub = SV - 1;
      for (i17 = 0; i17 <= loop_ub; i17++) {
        b_tris_bndID->data[i17 + tris_bndID->size[0]] = 0;
      }

      i17 = tris_bndID->size[0];
      tris_bndID->size[0] = b_tris_bndID->size[0];
      emxEnsureCapacity((emxArray__common *)tris_bndID, i17, (int32_T)sizeof
                        (int32_T));
      loop_ub = b_tris_bndID->size[0] - 1;
      for (i17 = 0; i17 <= loop_ub; i17++) {
        tris_bndID->data[i17] = b_tris_bndID->data[i17];
      }
    }

    /* 'triangulate_cellface:67' tri = zeros(1,3); */
    for (i17 = 0; i17 < 3; i17++) {
      tri[i17] = 0;
    }

    /* 'triangulate_cellface:68' for i=1:size(Ftris,1) */
    for (i = 0; i <= Ftris->size[0] - 1; i++) {
      /* 'triangulate_cellface:69' tri(1)=Face_ps(Ftris(i,1)); */
      tri[0] = Face_ps->data[Ftris->data[(int32_T)(1.0 + (real_T)i) - 1] - 1];

      /* 'triangulate_cellface:70' tri(2)=Face_ps(Ftris(i,2)); */
      tri[1] = Face_ps->data[Ftris->data[((int32_T)(1.0 + (real_T)i) +
        Ftris->size[0]) - 1] - 1];

      /* 'triangulate_cellface:71' tri(3)=Face_ps(Ftris(i,3)); */
      tri[2] = Face_ps->data[Ftris->data[((int32_T)(1.0 + (real_T)i) +
        (Ftris->size[0] << 1)) - 1] - 1];

      /* 'triangulate_cellface:73' chk = chktri_in_trislocal(tri,tris_local,tidx); */
      /* 'triangulate_cellface:95' chk = int32(0); */
      chk = 0;

      /* 'triangulate_cellface:96' for ii=1:tidx-1 */
      loop_ub = 1;
      exitg1 = 0U;
      while ((exitg1 == 0U) && (loop_ub <= tidx)) {
        /* 'triangulate_cellface:97' [SV] = chk_tri_orientation(tri,tris_local(ii,:)); */
        /*  Given two triangles, this function checks 1. if they have the same set of */
        /*  vertices, and 2. if they have the same orientation. */
        /* 'chk_tri_orientation:4' T = [Tri1 Tri2]; */
        /* 'chk_tri_orientation:5' T = unique(T); */
        for (i17 = 0; i17 < 3; i17++) {
          b_tri[i17] = tri[i17];
        }

        for (i17 = 0; i17 < 3; i17++) {
          b_tri[i17 + 3] = tris_local->data[(loop_ub + tris_local->size[0] * i17)
            - 1];
        }

        c_unique(b_tri, T);

        /* 'chk_tri_orientation:6' if (size(T,2) == 3) */
        if (T->size[1] == 3) {
          /* 'chk_tri_orientation:7' SameVert = 1; */
          SV = 1;
        } else {
          /* 'chk_tri_orientation:8' else */
          /* 'chk_tri_orientation:9' SameVert = 0; */
          SV = 0;
        }

        /* 'chk_tri_orientation:12' comb1 = [Tri1(1) Tri1(2) Tri1(3)]; */
        /* 'chk_tri_orientation:13' comb2 = [Tri1(2) Tri1(3) Tri1(1)]; */
        /* 'chk_tri_orientation:14' comb3 = [Tri1(3) Tri1(1) Tri1(2)]; */
        /* 'chk_tri_orientation:17' if (((Tri2(1) == comb1(1)) && (Tri2(2) == comb1(2)) && (Tri2(3) == comb1(3))) || ... */
        /* 'chk_tri_orientation:18'     ((Tri2(1) == comb2(1)) && (Tri2(2) == comb2(2)) && (Tri2(3) == comb2(3))) || ... */
        /* 'chk_tri_orientation:19'     ((Tri2(1) == comb3(1)) && (Tri2(2) == comb3(2)) && (Tri2(3) == comb3(3)))) */
        /* 'triangulate_cellface:98' if (SV == 1) */
        if (SV == 1) {
          /* 'triangulate_cellface:99' chk = int32(1); */
          chk = 1;
          exitg1 = 1U;
        } else {
          loop_ub++;
        }
      }

      /*  function visualize_FaceEdges(ps,edges) */
      /*  for i=1:size(edges,1) */
      /*      x = [ps(edges(i,1),1);ps(edges(i,2),1)]; */
      /*      y = [ps(edges(i,1),2);ps(edges(i,2),2)]; */
      /*      z = [ps(edges(i,1),3);ps(edges(i,2),3)]; */
      /*      plot3(x,y,z); */
      /*       text(x(1),y(1),z(1),['\leftarrow',num2str(edges(i,1))]); */
      /*      hold on */
      /*      text(x(2),y(2),z(2),['\rightarrow',num2str(edges(i,2))]); */
      /*      scatter3(ps(edges(i,1),1),ps(edges(i,1),2),ps(edges(i,1),3),'filled'); */
      /*      scatter3(ps(edges(i,2),1),ps(edges(i,2),2),ps(edges(i,2),3),'filled'); */
      /*  end */
      /*  end */
      /* 'triangulate_cellface:74' if (chk == 0) */
      if (chk == 0) {
        /* 'triangulate_cellface:75' nrm = cross_col(ps_local(tri(2),:)-ps_local(tri(1),:),ps_local(tri(3),:)-ps_local(tri(1),:)); */
        for (i17 = 0; i17 < 3; i17++) {
          a[i17] = ps_local->data[(tri[1] + ps_local->size[0] * i17) - 1] -
            ps_local->data[(tri[0] + ps_local->size[0] * i17) - 1];
        }

        for (i17 = 0; i17 < 3; i17++) {
          b[i17] = ps_local->data[(tri[2] + ps_local->size[0] * i17) - 1] -
            ps_local->data[(tri[0] + ps_local->size[0] * i17) - 1];
        }

        /* CROSS_COL Efficient routine for computing cross product of two  */
        /* 3-dimensional column vectors. */
        /*  CROSS_COL(A,B) Efficiently computes the cross product between */
        /*  3-dimensional column vector A, and 3-dimensional column vector B. */
        /* 'cross_col:7' c = [a(2)*b(3)-a(3)*b(2); a(3)*b(1)-a(1)*b(3); a(1)*b(2)-a(2)*b(1)]; */
        nrm[0] = a[1] * b[2] - a[2] * b[1];
        nrm[1] = a[2] * b[0] - a[0] * b[2];
        nrm[2] = a[0] * b[1] - a[1] * b[0];

        /* 'triangulate_cellface:76' nrm = nrm/sqrt(nrm'*nrm); */
        chkdir = sqrt(eml_xdot(3, nrm, 1, 1, nrm, 1, 1));
        for (i17 = 0; i17 < 3; i17++) {
          nrm[i17] /= chkdir;
        }

        /* 'triangulate_cellface:77' chkdir = dot(nrm,Face_nrm(ii,:)); */
        chkdir = 0.0;
        SV = 0;
        chk = 0;
        for (loop_ub = 0; loop_ub < 3; loop_ub++) {
          chkdir += nrm[SV] * Face_nrm[ii + 6 * chk];
          SV++;
          chk++;
        }

        /* 'triangulate_cellface:78' if (chkdir ~= 1) */
        if (chkdir != 1.0) {
          /* 'triangulate_cellface:79' a = tri(1); */
          SV = tri[0];

          /* 'triangulate_cellface:80' tri(1) = tri(2); */
          tri[0] = tri[1];

          /* 'triangulate_cellface:81' tri(2) = a; */
          tri[1] = SV;
        }

        /* 'triangulate_cellface:84' tris_local(tidx,1:3)=tri; */
        for (i17 = 0; i17 < 3; i17++) {
          tris_local->data[tidx + tris_local->size[0] * i17] = tri[i17];
        }

        /* 'triangulate_cellface:85' tris_bndID(tidx)=ii; */
        tris_bndID->data[tidx] = 1 + ii;

        /* 'triangulate_cellface:86' tidx = tidx + 1; */
        tidx++;
      }
    }
  }

  emxFree_int32_T(&b_tris_bndID);
  emxFree_int32_T(&b_tris_local);
  emxFree_real_T(&b_ps_local);
  emxFree_int32_T(&T);
  emxFree_int32_T(&Ftris);
  emxFree_real_T(&Fps_new);
  emxFree_int32_T(&Face_elems);
  emxFree_int32_T(&b_Face_ps);
  emxFree_int32_T(&Elems);
  emxFree_real_T(&Fps);
  emxFree_int32_T(&Face_ps);
  emxFree_int32_T(&opphes);

  /* 'triangulate_cellface:91' tris_local = tris_local(1:tidx-1,:); */
  if (1 > tidx) {
    i17 = 0;
  } else {
    i17 = tidx;
  }

  emxInit_int32_T(&c_tris_local, 2);
  chk = c_tris_local->size[0] * c_tris_local->size[1];
  c_tris_local->size[0] = i17;
  c_tris_local->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)c_tris_local, chk, (int32_T)sizeof
                    (int32_T));
  for (chk = 0; chk < 3; chk++) {
    loop_ub = i17 - 1;
    for (SV = 0; SV <= loop_ub; SV++) {
      c_tris_local->data[SV + c_tris_local->size[0] * chk] = tris_local->data[SV
        + tris_local->size[0] * chk];
    }
  }

  i17 = tris_local->size[0] * tris_local->size[1];
  tris_local->size[0] = c_tris_local->size[0];
  tris_local->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)tris_local, i17, (int32_T)sizeof(int32_T));
  for (i17 = 0; i17 < 3; i17++) {
    loop_ub = c_tris_local->size[0] - 1;
    for (chk = 0; chk <= loop_ub; chk++) {
      tris_local->data[chk + tris_local->size[0] * i17] = c_tris_local->data[chk
        + c_tris_local->size[0] * i17];
    }
  }

  emxFree_int32_T(&c_tris_local);

  /* 'triangulate_cellface:92' tris_bndID = tris_bndID(1:tidx-1); */
  if (1 > tidx) {
    tidx = 0;
  }

  b_emxInit_int32_T(&c_tris_bndID, 1);
  i17 = c_tris_bndID->size[0];
  c_tris_bndID->size[0] = tidx;
  emxEnsureCapacity((emxArray__common *)c_tris_bndID, i17, (int32_T)sizeof
                    (int32_T));
  loop_ub = tidx - 1;
  for (i17 = 0; i17 <= loop_ub; i17++) {
    c_tris_bndID->data[i17] = tris_bndID->data[i17];
  }

  i17 = tris_bndID->size[0];
  tris_bndID->size[0] = c_tris_bndID->size[0];
  emxEnsureCapacity((emxArray__common *)tris_bndID, i17, (int32_T)sizeof(int32_T));
  loop_ub = c_tris_bndID->size[0] - 1;
  for (i17 = 0; i17 <= loop_ub; i17++) {
    tris_bndID->data[i17] = c_tris_bndID->data[i17];
  }

  emxFree_int32_T(&c_tris_bndID);
}

/*
 *
 */
static void unique(const emxArray_int32_T *a, emxArray_int32_T *b)
{
  int32_T k0;
  int32_T m;
  emxArray_int32_T *idx;
  int32_T n;
  int32_T np1;
  int32_T k;
  emxArray_int32_T *idx0;
  int32_T i;
  int32_T nb;
  int32_T p;
  int32_T q;
  int32_T qEnd;
  int32_T kEnd;
  int32_T b_a[2];
  emxArray_int32_T b_b;
  int32_T c_b[2];
  int32_T exitg1;
  boolean_T b_p;
  int32_T d_b[2];
  int32_T e_b[2];
  emxArray_int32_T f_b;
  int32_T g_b[2];
  int32_T h_b[2];
  emxArray_int32_T *i_b;
  int32_T j_b[2];
  if (a->size[0] == 0) {
    k0 = b->size[0];
    b->size[0] = a->size[0];
    emxEnsureCapacity((emxArray__common *)b, k0, (int32_T)sizeof(int32_T));
    m = a->size[0] - 1;
    for (k0 = 0; k0 <= m; k0++) {
      b->data[k0] = a->data[k0];
    }
  } else {
    b_emxInit_int32_T(&idx, 1);
    n = a->size[0];
    k0 = idx->size[0];
    idx->size[0] = n;
    emxEnsureCapacity((emxArray__common *)idx, k0, (int32_T)sizeof(int32_T));
    np1 = n + 1;
    for (k = 1; k <= n; k++) {
      idx->data[k - 1] = k;
    }

    for (k = 1; k <= n - 1; k += 2) {
      if (eml_sort_le(a, 1, k, k + 1)) {
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    b_emxInit_int32_T(&idx0, 1);
    k0 = idx0->size[0];
    idx0->size[0] = n;
    emxEnsureCapacity((emxArray__common *)idx0, k0, (int32_T)sizeof(int32_T));
    m = n - 1;
    for (k0 = 0; k0 <= m; k0++) {
      idx0->data[k0] = 1;
    }

    i = 2;
    while (i < n) {
      m = i << 1;
      k0 = 1;
      for (nb = 1 + i; nb < np1; nb = qEnd + i) {
        p = k0;
        q = nb;
        qEnd = k0 + m;
        if (qEnd > np1) {
          qEnd = np1;
        }

        k = 0;
        kEnd = qEnd - k0;
        while (k + 1 <= kEnd) {
          if (eml_sort_le(a, 1, idx->data[p - 1], idx->data[q - 1])) {
            idx0->data[k] = idx->data[p - 1];
            p++;
            if (p == nb) {
              while (q < qEnd) {
                k++;
                idx0->data[k] = idx->data[q - 1];
                q++;
              }
            }
          } else {
            idx0->data[k] = idx->data[q - 1];
            q++;
            if (q == qEnd) {
              while (p < nb) {
                k++;
                idx0->data[k] = idx->data[p - 1];
                p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k + 1 <= kEnd; k++) {
          idx->data[(k0 + k) - 1] = idx0->data[k];
        }

        k0 = qEnd;
      }

      i = m;
    }

    k0 = b->size[0];
    b->size[0] = a->size[0];
    emxEnsureCapacity((emxArray__common *)b, k0, (int32_T)sizeof(int32_T));
    m = a->size[0] - 1;
    for (k0 = 0; k0 <= m; k0++) {
      b->data[k0] = a->data[k0];
    }

    m = a->size[0];
    k0 = idx0->size[0];
    idx0->size[0] = m;
    emxEnsureCapacity((emxArray__common *)idx0, k0, (int32_T)sizeof(int32_T));
    for (i = 0; i + 1 <= m; i++) {
      b_a[0] = a->size[0];
      b_a[1] = 1;
      b_b = *a;
      b_b.size = (int32_T *)&b_a;
      b_b.numDimensions = 1;
      idx0->data[i] = b_b.data[idx->data[i] - 1];
    }

    emxFree_int32_T(&idx);
    for (i = 0; i + 1 <= m; i++) {
      c_b[0] = b->size[0];
      c_b[1] = 1;
      b_b = *b;
      b_b.size = (int32_T *)&c_b;
      b_b.numDimensions = 1;
      b_b.data[i] = idx0->data[i];
    }

    emxFree_int32_T(&idx0);
    nb = 0;
    m = a->size[0];
    k = 0;
    while (k + 1 <= m) {
      k0 = k;
      do {
        exitg1 = 0U;
        k++;
        if (k + 1 > m) {
          exitg1 = 1U;
        } else {
          b_p = FALSE;
          d_b[0] = b->size[0];
          d_b[1] = 1;
          e_b[0] = b->size[0];
          e_b[1] = 1;
          b_b = *b;
          b_b.size = (int32_T *)&d_b;
          b_b.numDimensions = 1;
          f_b = *b;
          f_b.size = (int32_T *)&e_b;
          f_b.numDimensions = 1;
          if (!(b_b.data[k0] == f_b.data[k])) {
            b_p = TRUE;
          }

          if (b_p) {
            exitg1 = 1U;
          }
        }
      } while (exitg1 == 0U);

      nb++;
      g_b[0] = b->size[0];
      g_b[1] = 1;
      h_b[0] = b->size[0];
      h_b[1] = 1;
      b_b = *b;
      b_b.size = (int32_T *)&g_b;
      b_b.numDimensions = 1;
      f_b = *b;
      f_b.size = (int32_T *)&h_b;
      f_b.numDimensions = 1;
      b_b.data[nb - 1] = f_b.data[k0];
    }

    if (1 > nb) {
      nb = 0;
    }

    b_emxInit_int32_T(&i_b, 1);
    m = b->size[0];
    j_b[0] = m;
    j_b[1] = 1;
    k0 = i_b->size[0];
    i_b->size[0] = nb;
    emxEnsureCapacity((emxArray__common *)i_b, k0, (int32_T)sizeof(int32_T));
    m = nb - 1;
    for (k0 = 0; k0 <= m; k0++) {
      b_b = *b;
      b_b.size = (int32_T *)&j_b;
      b_b.numDimensions = 1;
      i_b->data[k0] = b_b.data[k0];
    }

    k0 = b->size[0];
    b->size[0] = i_b->size[0];
    emxEnsureCapacity((emxArray__common *)b, k0, (int32_T)sizeof(int32_T));
    m = i_b->size[0] - 1;
    for (k0 = 0; k0 <= m; k0++) {
      b->data[k0] = i_b->data[k0];
    }

    emxFree_int32_T(&i_b);
  }
}

/*
 * function [map] = update_map_component(comp,num_comps,comp_list,map)
 * Input : comp = 1x4 array, the list of components of the 4 ngb tets
 *  mincomp = minimum value of comp
 *  sign = sign of the component
 *  comp_list,pos_map,neg_map are global lists.
 */
static void update_map_component(const emxArray_int32_T *comp, int32_T num_comps,
  const emxArray_int32_T *comp_list, emxArray_int32_T *map)
{
  int32_T idx;
  int32_T idx_val;
  int32_T b_index;
  emxArray_real_T *L;
  int32_T k;
  int32_T j;
  real_T p;
  real_T count;
  int32_T i;
  emxArray_real_T *LL;

  /* 'update_map_component:6' idx = int32(0); */
  idx = -1;

  /* 'update_map_component:6' idx_val = int32(0); */
  idx_val = 0;

  /* 'update_map_component:6' index = int32(0); */
  b_index = -1;

  /* 'update_map_component:7' if (size(comp,2)==1) */
  if (comp->size[1] == 1) {
  } else {
    emxInit_real_T(&L, 2);

    /* 'update_map_component:8' else */
    /* 'update_map_component:9' L = zeros(1,num_comps); */
    k = L->size[0] * L->size[1];
    L->size[0] = 1;
    L->size[1] = num_comps;
    emxEnsureCapacity((emxArray__common *)L, k, (int32_T)sizeof(real_T));
    j = num_comps - 1;
    for (k = 0; k <= j; k++) {
      L->data[k] = 0.0;
    }

    /* 'update_map_component:9' p = 1; */
    p = 1.0;

    /* 'update_map_component:10' for ii=1:size(comp,2) */
    for (j = 0; j <= comp->size[1] - 1; j++) {
      /*  Check if this component is already in L. If not put it into */
      /*  L. */
      /* 'update_map_component:13' key = comp(ii); */
      /* 'update_map_component:13' count = 0; */
      count = 0.0;

      /* 'update_map_component:14' for i=1:num_comps */
      i = 1;
      while ((i <= num_comps) && (!(L->data[i - 1] == (real_T)comp->data[j]))) {
        /* 'update_map_component:15' if (L(i) == key) */
        /* 'update_map_component:17' else */
        /* 'update_map_component:18' count = count+1; */
        count++;
        i++;
      }

      /* 'update_map_component:21' if (count == num_comps) */
      if (count == (real_T)num_comps) {
        /* 'update_map_component:22' L(p) = key; */
        L->data[(int32_T)p - 1] = (real_T)comp->data[j];

        /* 'update_map_component:23' p = p+1; */
        p++;
      }

      /*  Check if this component is mapped to itself: If yes, leave it */
      /*  as it is. Otherwise it belongs to a star. */
      /* 'update_map_component:27' for jj=1:num_comps */
      for (k = 0; k + 1 <= num_comps; k++) {
        /* 'update_map_component:28' if (comp_list(jj)==key) */
        if (comp_list->data[k] == comp->data[j]) {
          /* 'update_map_component:29' idx = jj; */
          idx = k;
        }
      }

      /* 'update_map_component:32' if (map(idx)==key) */
      if (map->data[idx] == comp->data[j]) {
        /*  Find all comps that are mapped to key */
        /* 'update_map_component:34' [L,p] = searchstar(key,idx,num_comps,comp_list,map,L,p); */
        searchstar(comp->data[j], idx + 1, num_comps, comp_list, map, L, &p);
      } else {
        /* 'update_map_component:35' else */
        /*  Find the comps that are also mapped to 'val' */
        /* 'update_map_component:37' val = map(idx); */
        /* 'update_map_component:37' count = 0; */
        count = 0.0;

        /* 'update_map_component:38' for i=1:num_comps */
        i = 1;
        while ((i <= num_comps) && (!(L->data[i - 1] == (real_T)map->data[idx])))
        {
          /* 'update_map_component:39' if (L(i) == val) */
          /* 'update_map_component:41' else */
          /* 'update_map_component:42' count = count+1; */
          count++;
          i++;
        }

        /* 'update_map_component:45' if (count == num_comps) */
        if (count == (real_T)num_comps) {
          /* 'update_map_component:46' L(p) = val; */
          L->data[(int32_T)p - 1] = (real_T)map->data[idx];

          /* 'update_map_component:47' p = p+1; */
          p++;
        }

        /* 'update_map_component:49' for jj=1:num_comps */
        for (k = 1; k <= num_comps; k++) {
          /* 'update_map_component:50' if (comp_list(jj)==val) */
          if (comp_list->data[k - 1] == map->data[idx]) {
            /* 'update_map_component:51' idx_val = jj; */
            idx_val = k;
          }
        }

        /* 'update_map_component:54' [L,p] = searchstar(val,idx_val,num_comps,comp_list,map,L,p); */
        searchstar(map->data[idx], idx_val, num_comps, comp_list, map, L, &p);
      }
    }

    /*  Finally map all nonzero elements of L to the minimum. */
    /* 'update_map_component:58' LL = nonzeros(L); */
    idx = L->size[1];
    j = 0;
    for (k = 0; k <= L->size[1] - 1; k++) {
      if (L->data[(int32_T)(1.0 + (real_T)k) - 1] != 0.0) {
        j++;
      }
    }

    b_emxInit_real_T(&LL, 1);
    k = LL->size[0];
    LL->size[0] = j;
    emxEnsureCapacity((emxArray__common *)LL, k, (int32_T)sizeof(real_T));
    i = -1;
    for (k = 0; k + 1 <= idx; k++) {
      if (L->data[k] != 0.0) {
        i++;
        LL->data[i] = L->data[k];
      }
    }

    emxFree_real_T(&L);

    /* 'update_map_component:59' minL = min(LL); */
    idx = LL->size[0];
    p = LL->data[0];
    if (idx > 1) {
      for (j = 1; j + 1 <= idx; j++) {
        if (LL->data[j] < p) {
          p = LL->data[j];
        }
      }
    }

    /* 'update_map_component:60' for jj=1:size(LL,1) */
    for (k = 0; k <= LL->size[0] - 1; k++) {
      /* 'update_map_component:61' if (LL(jj) ~= minL) */
      if (LL->data[k] != p) {
        /* 'update_map_component:62' for kk=1:num_comps */
        for (j = 0; j + 1 <= num_comps; j++) {
          /* 'update_map_component:63' if (comp_list(kk)==LL(jj)) */
          if ((real_T)comp_list->data[j] == LL->data[k]) {
            /* 'update_map_component:64' index = kk; */
            b_index = j;
          }
        }

        /* 'update_map_component:67' map(index) = minL; */
        count = p;
        count = count < 0.0 ? ceil(count - 0.5) : floor(count + 0.5);
        map->data[b_index] = (int32_T)count;
      }
    }

    emxFree_real_T(&LL);
  }
}

/*
 * function [face_type,face_sign,pos_map,neg_map,pos_vol,neg_vol] = compute_volume_usingtets(cell_bnds,face_type,face_sign,...
 *     ps,num_points,tris,num_tris,comp_index,num_comps,comp_list,pos_map,neg_map,pos_vol,neg_vol)
 * % INPUT DESCRIPTION :
 *  CELL DATA : c_centers = cell center, grid_spacing = [hx,hy,hz],
 *  face_type = -2 by default,
 *            = -1 if it is a crossed face,
 *            = othervalue, otherwise
 *  face_sign = true/1 : default or positive sign
 *            = false/0 : negative sign
 *  TRIANGLE DATA : ps = position of vertices, tris = triangle connectivity,
 *  comp_index = component of each triangle in tris
 *  COMPONENT DATA : comp_list = list of components in one processor,
 *  pos_map = mapping of the positive components, neg_map = mapping of the
 *  negative components.
 *  comp_list, pos_map, neg_map are kx1 arrays
 *  VOLUME DATA : pos_vol, neg_vol : default value is 0 and their dimension
 *  is kx1 arrays.
 * % INPUT TYPE
 * assert((size(c_centers,1)==1)&& (size(c_centers,2)==3) && isa(c_centers,'double'));
 * assert((size(grid_spacing,1)==1)&& (size(grid_spacing,2)==3) && isa(grid_spacing,'double'));
 * assert(isscalar(flag_cell)&& isa(flag_cell,'int32'));
 */
void compute_volume_usingtets(const real_T cell_bnds[6], int32_T face_type[6],
  boolean_T face_sign[6], const emxArray_real_T *ps, int32_T num_points, const
  emxArray_int32_T *tris, int32_T num_tris, const emxArray_int32_T *comp_index,
  int32_T num_comps, const emxArray_int32_T *comp_list, emxArray_int32_T
  *pos_map, emxArray_int32_T *neg_map, emxArray_real_T *pos_vol, emxArray_real_T
  *neg_vol)
{
  emxArray_int32_T *tris_cell;
  emxArray_int32_T *comp_index_cell;
  int32_T tris_idx;
  int32_T ix;
  emxArray_int32_T *b_tris_cell;
  int32_T i10;
  emxArray_int32_T *b_comp_index_cell;
  emxArray_real_T *ps_vol;
  emxArray_int32_T *tets;
  emxArray_int32_T *tris_local;
  emxArray_int32_T *tris_bndID;
  emxArray_int32_T *ps_localID;
  int32_T iy;
  real_T volume;
  emxArray_int32_T *c_comp_index_cell;
  emxArray_int32_T *TETcomp;
  int32_T idx;
  int32_T idx_val;
  int32_T ii;
  int32_T b_tets;
  real_T AD[3];
  real_T BD[3];
  real_T CD[3];
  real_T c[3];
  emxInit_int32_T(&tris_cell, 2);
  b_emxInit_int32_T(&comp_index_cell, 1);

  /* 'compute_volume_usingtets:22' assert((size(cell_bnds,1)==1)&& (size(cell_bnds,2)==6) && isa(cell_bnds,'double')); */
  /* 'compute_volume_usingtets:23' assert((size(face_type,1)==1)&& (size(face_type,2)==6) && isa(face_type,'int32')); */
  /* 'compute_volume_usingtets:24' assert((isa(face_sign,'logical'))&&(size(face_sign,1)==1)&& (size(face_sign,2)==6)); */
  /* 'compute_volume_usingtets:26' assert(isscalar(num_points)&& isscalar(num_tris)); */
  /* 'compute_volume_usingtets:27' assert(isa(num_points,'int32')&& isa(num_tris,'int32')); */
  /* 'compute_volume_usingtets:28' assert((size(ps,2)==3) && (size(ps,1)>=1) && isa(ps,'double')); */
  /* 'compute_volume_usingtets:29' assert((size(tris,2)==3) && (size(tris,1)>=1) && isa(tris,'int32')); */
  /* 'compute_volume_usingtets:30' assert((size(comp_index,2)==1) && (size(comp_index,1)>=1) && isa(comp_index,'int32')); */
  /* 'compute_volume_usingtets:32' assert(isscalar(num_comps)&& isa(num_comps,'int32')); */
  /* 'compute_volume_usingtets:33' assert((size(comp_list,2)==1)&&(size(comp_list,1)>=1) && isa(comp_list,'int32')); */
  /* 'compute_volume_usingtets:34' assert((size(pos_map,2)==1)&&(size(pos_map,1)>=1) && isa(pos_map,'int32')); */
  /* 'compute_volume_usingtets:35' assert((size(neg_map,2)==1)&&(size(neg_map,1)>=1) && isa(neg_map,'int32')); */
  /* 'compute_volume_usingtets:36' assert((size(pos_vol,2)==1)&&(size(pos_vol,1)>=1) && isa(pos_vol,'double')); */
  /* 'compute_volume_usingtets:37' assert((size(neg_vol,2)==1)&&(size(neg_vol,1)>=1) && isa(neg_vol,'double')); */
  /* % GET A LIST OF TRIANGLES THAT INTERSECT THE CELL OR LIE TOTALLY INSIDE IT */
  /* TODO: Make get_trilist_inMB into a separate preprocessing step. */
  /* [tris_idx,tris_cell,comp_index_cell,flag_cell] = get_trilist_inMB(c_centers,grid_spacing,ps,tris,num_tris,comp_index); */
  /* 'compute_volume_usingtets:42' [flag_cell,tris_cell,tris_idx,comp_index_cell] = get_trisincell(cell_bnds,ps,num_tris,tris,comp_index) ; */
  get_trisincell(cell_bnds, ps, num_tris, tris, comp_index, tris_cell,
                 comp_index_cell, &ix, &tris_idx);

  /* 'compute_volume_usingtets:43' if (flag_cell == 1) */
  if (ix == 1) {
    emxInit_int32_T(&b_tris_cell, 2);

    /*     %% TETRAHEDRALIZATION */
    /* 'compute_volume_usingtets:45' flag_tri = 0; */
    /* 'compute_volume_usingtets:46' [flag_tri,nv,ps_vol,tets,tris_labelled,tris_comp,tris_local,tris_bndID,ps_localID] = tetrahedralize_cell(cell_bnds,... */
    /* 'compute_volume_usingtets:47'         ps,num_points,tris_cell,tris_idx,comp_index_cell,flag_tri); */
    i10 = b_tris_cell->size[0] * b_tris_cell->size[1];
    b_tris_cell->size[0] = tris_cell->size[0];
    b_tris_cell->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)b_tris_cell, i10, (int32_T)sizeof
                      (int32_T));
    ix = tris_cell->size[0] * tris_cell->size[1] - 1;
    for (i10 = 0; i10 <= ix; i10++) {
      b_tris_cell->data[i10] = tris_cell->data[i10];
    }

    b_emxInit_int32_T(&b_comp_index_cell, 1);
    i10 = b_comp_index_cell->size[0];
    b_comp_index_cell->size[0] = comp_index_cell->size[0];
    emxEnsureCapacity((emxArray__common *)b_comp_index_cell, i10, (int32_T)
                      sizeof(int32_T));
    ix = comp_index_cell->size[0] - 1;
    for (i10 = 0; i10 <= ix; i10++) {
      b_comp_index_cell->data[i10] = comp_index_cell->data[i10];
    }

    emxInit_real_T(&ps_vol, 2);
    emxInit_int32_T(&tets, 2);
    emxInit_int32_T(&tris_local, 2);
    b_emxInit_int32_T(&tris_bndID, 1);
    emxInit_int32_T(&ps_localID, 2);
    tetrahedralize_cell(cell_bnds, ps, num_points, b_tris_cell, tris_idx,
                        b_comp_index_cell, ps_vol, tets, tris_cell,
                        comp_index_cell, tris_local, tris_bndID, ps_localID,
                        &volume, &iy);

    /* 'compute_volume_usingtets:49' if (flag_tri ==1) */
    emxFree_int32_T(&b_comp_index_cell);
    emxFree_int32_T(&b_tris_cell);
    if (volume == 1.0) {
      b_emxInit_int32_T(&c_comp_index_cell, 1);

      /*         %% LABEL TETS AND GET CONNECTED COMPONENT MAPS */
      /* 'compute_volume_usingtets:52' [TETcomp,TETsign,pos_map,neg_map] = label_tets(nv,tets,tris_labelled,tris_comp,num_comps,comp_list,pos_map,neg_map); */
      i10 = c_comp_index_cell->size[0];
      c_comp_index_cell->size[0] = comp_index_cell->size[0];
      emxEnsureCapacity((emxArray__common *)c_comp_index_cell, i10, (int32_T)
                        sizeof(int32_T));
      ix = comp_index_cell->size[0] - 1;
      for (i10 = 0; i10 <= ix; i10++) {
        c_comp_index_cell->data[i10] = comp_index_cell->data[i10];
      }

      emxInit_int32_T(&TETcomp, 2);
      label_tets(iy, tets, tris_cell, c_comp_index_cell, num_comps, comp_list,
                 pos_map, neg_map, TETcomp, comp_index_cell);

      /*         %% GET FACE LABELS */
      /* 'compute_volume_usingtets:55' [face_type,face_sign] = label_cellface(ps_localID,tris_local,tris_bndID,tets,TETcomp,TETsign,num_comps,comp_list,pos_map,neg_map,face_type,face_sign); */
      label_cellface(ps_localID, tris_local, tris_bndID, tets, TETcomp,
                     comp_index_cell, num_comps, comp_list, pos_map, neg_map,
                     face_type, face_sign);

      /*         %% COMPUTE VOLUME OF CONNECTED COMPONENTS */
      /* 'compute_volume_usingtets:58' tnum = size(tets,1); */
      /* 'compute_volume_usingtets:59' idx = int32(0); */
      idx = -1;

      /* 'compute_volume_usingtets:59' idx_val = int32(0); */
      idx_val = -1;

      /* 'compute_volume_usingtets:59' volume = 0; */
      /* 'compute_volume_usingtets:60' for ii=1:tnum */
      ii = 0;
      emxFree_int32_T(&c_comp_index_cell);
      while (ii <= tets->size[0] - 1) {
        /* 'compute_volume_usingtets:61' a = ps_vol(tets(ii,1),:); */
        /* 'compute_volume_usingtets:62' b = ps_vol(tets(ii,2),:); */
        /* 'compute_volume_usingtets:63' c = ps_vol(tets(ii,3),:); */
        /* 'compute_volume_usingtets:64' d = ps_vol(tets(ii,4),:); */
        /* 'compute_volume_usingtets:65' AD = a - d; */
        iy = tets->data[(int32_T)(1.0 + (real_T)ii) - 1];
        b_tets = tets->data[((int32_T)(1.0 + (real_T)ii) + tets->size[0] * 3) -
          1];
        for (i10 = 0; i10 < 3; i10++) {
          AD[i10] = ps_vol->data[(iy + ps_vol->size[0] * i10) - 1] -
            ps_vol->data[(b_tets + ps_vol->size[0] * i10) - 1];
        }

        /* 'compute_volume_usingtets:66' BD = b - d; */
        iy = tets->data[((int32_T)(1.0 + (real_T)ii) + tets->size[0]) - 1];
        b_tets = tets->data[((int32_T)(1.0 + (real_T)ii) + tets->size[0] * 3) -
          1];
        for (i10 = 0; i10 < 3; i10++) {
          BD[i10] = ps_vol->data[(iy + ps_vol->size[0] * i10) - 1] -
            ps_vol->data[(b_tets + ps_vol->size[0] * i10) - 1];
        }

        /* 'compute_volume_usingtets:67' CD = c - d; */
        iy = tets->data[((int32_T)(1.0 + (real_T)ii) + (tets->size[0] << 1)) - 1];
        b_tets = tets->data[((int32_T)(1.0 + (real_T)ii) + tets->size[0] * 3) -
          1];
        for (i10 = 0; i10 < 3; i10++) {
          CD[i10] = ps_vol->data[(iy + ps_vol->size[0] * i10) - 1] -
            ps_vol->data[(b_tets + ps_vol->size[0] * i10) - 1];
        }

        /*  r_centroid = (a(1)+b(1)+c(1)+d(1))/4; % In frontier: Order of */
        /*  coordinates is theta, z, r, hence */
        /* 'compute_volume_usingtets:70' r_centroid = (a(3)+b(3)+c(3)+d(3))/4; */
        /* 'compute_volume_usingtets:71' volume =  r_centroid*abs(dot(AD,cross(BD,CD))/6); */
        c[0] = BD[1] * CD[2] - BD[2] * CD[1];
        c[1] = BD[2] * CD[0] - BD[0] * CD[2];
        c[2] = BD[0] * CD[1] - BD[1] * CD[0];
        volume = 0.0;
        ix = 0;
        iy = 0;
        for (tris_idx = 0; tris_idx < 3; tris_idx++) {
          volume += AD[ix] * c[iy];
          ix++;
          iy++;
        }

        iy = tets->data[(int32_T)(1.0 + (real_T)ii) - 1];
        b_tets = tets->data[((int32_T)(1.0 + (real_T)ii) + tets->size[0]) - 1];
        tris_idx = tets->data[((int32_T)(1.0 + (real_T)ii) + (tets->size[0] << 1))
          - 1];
        ix = tets->data[((int32_T)(1.0 + (real_T)ii) + tets->size[0] * 3) - 1];
        volume = (((ps_vol->data[(iy + (ps_vol->size[0] << 1)) - 1] +
                    ps_vol->data[(b_tets + (ps_vol->size[0] << 1)) - 1]) +
                   ps_vol->data[(tris_idx + (ps_vol->size[0] << 1)) - 1]) +
                  ps_vol->data[(ix + (ps_vol->size[0] << 1)) - 1]) / 4.0 * fabs
          (volume / 6.0);

        /*  Cartesian Volume */
        /* volume = abs(dot(AD,cross(BD,CD))/6); % Cylindrical Volume */
        /* 'compute_volume_usingtets:73' if (TETsign(ii) == 1) */
        if (comp_index_cell->data[(int32_T)(1.0 + (real_T)ii) - 1] == 1) {
          /* 'compute_volume_usingtets:74' comp = TETcomp(ii,5); */
          /* 'compute_volume_usingtets:75' for jj=1:num_comps */
          for (ix = 0; ix + 1 <= num_comps; ix++) {
            /* 'compute_volume_usingtets:76' if (comp_list(jj)==comp) */
            if (comp_list->data[ix] == TETcomp->data[((int32_T)(1.0 + (real_T)ii)
                 + (TETcomp->size[0] << 2)) - 1]) {
              /* 'compute_volume_usingtets:77' idx = jj; */
              idx = ix;
            }
          }

          /* 'compute_volume_usingtets:80' val = pos_map(idx); */
          /* 'compute_volume_usingtets:81' for jj=1:num_comps */
          for (ix = 0; ix + 1 <= num_comps; ix++) {
            /* 'compute_volume_usingtets:82' if (comp_list(jj)==val) */
            if (comp_list->data[ix] == pos_map->data[idx]) {
              /* 'compute_volume_usingtets:83' idx_val = jj; */
              idx_val = ix;
            }
          }

          /* 'compute_volume_usingtets:86' if ((pos_vol(idx) ~= 0) && (val ~= comp)) */
          if ((pos_vol->data[idx] != 0.0) && (pos_map->data[idx] !=
               TETcomp->data[((int32_T)(1.0 + (real_T)ii) + (TETcomp->size[0] <<
                 2)) - 1])) {
            /* 'compute_volume_usingtets:87' pos_vol(idx_val) = pos_vol(idx_val) + pos_vol(idx); */
            pos_vol->data[idx_val] += pos_vol->data[idx];

            /* 'compute_volume_usingtets:88' pos_vol(idx) = 0; */
            pos_vol->data[idx] = 0.0;
          }

          /* 'compute_volume_usingtets:90' pos_vol(idx_val) = pos_vol(idx_val) + volume; */
          pos_vol->data[idx_val] += volume;
        } else {
          if (comp_index_cell->data[(int32_T)(1.0 + (real_T)ii) - 1] == -1) {
            /* 'compute_volume_usingtets:91' elseif (TETsign(ii) == -1) */
            /* 'compute_volume_usingtets:92' comp = TETcomp(ii,5); */
            /* 'compute_volume_usingtets:93' for jj=1:num_comps */
            for (ix = 0; ix + 1 <= num_comps; ix++) {
              /* 'compute_volume_usingtets:94' if (comp_list(jj)==comp) */
              if (comp_list->data[ix] == TETcomp->data[((int32_T)(1.0 + (real_T)
                    ii) + (TETcomp->size[0] << 2)) - 1]) {
                /* 'compute_volume_usingtets:95' idx = jj; */
                idx = ix;
              }
            }

            /* 'compute_volume_usingtets:98' val = neg_map(idx); */
            /* 'compute_volume_usingtets:99' for jj=1:num_comps */
            for (ix = 0; ix + 1 <= num_comps; ix++) {
              /* 'compute_volume_usingtets:100' if (comp_list(jj)==val) */
              if (comp_list->data[ix] == neg_map->data[idx]) {
                /* 'compute_volume_usingtets:101' idx_val = jj; */
                idx_val = ix;
              }
            }

            /* 'compute_volume_usingtets:104' if ((neg_vol(idx) ~= 0) && (val ~= comp)) */
            if ((neg_vol->data[idx] != 0.0) && (neg_map->data[idx] !=
                 TETcomp->data[((int32_T)(1.0 + (real_T)ii) + (TETcomp->size[0] <<
                   2)) - 1])) {
              /* 'compute_volume_usingtets:105' neg_vol(idx_val) = neg_vol(idx_val) + neg_vol(idx); */
              neg_vol->data[idx_val] += neg_vol->data[idx];

              /* 'compute_volume_usingtets:106' neg_vol(idx) = 0; */
              neg_vol->data[idx] = 0.0;
            }

            /* 'compute_volume_usingtets:108' neg_vol(idx_val) = neg_vol(idx_val) + volume; */
            neg_vol->data[idx_val] += volume;
          }
        }

        ii++;
      }

      emxFree_int32_T(&TETcomp);
    } else {
      /* 'compute_volume_usingtets:111' else */
    }

    emxFree_int32_T(&ps_localID);
    emxFree_int32_T(&tris_bndID);
    emxFree_int32_T(&tris_local);
    emxFree_int32_T(&tets);
    emxFree_real_T(&ps_vol);
  } else {
    /* 'compute_volume_usingtets:114' else */
  }

  emxFree_int32_T(&comp_index_cell);
  emxFree_int32_T(&tris_cell);
}

void compute_volume_usingtets_initialize(void)
{
}

void compute_volume_usingtets_terminate(void)
{
  /* (no terminate code required) */
}

emxArray_int32_T *emxCreateND_int32_T(int32_T numDimensions, int32_T *size)
{
  emxArray_int32_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  c_emxInit_int32_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (int32_T *)calloc((uint32_T)numEl, sizeof(int32_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  c_emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_int32_T *emxCreateWrapperND_int32_T(int32_T *data, int32_T
  numDimensions, int32_T *size)
{
  emxArray_int32_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  c_emxInit_int32_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T numDimensions,
  int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  c_emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_int32_T *emxCreateWrapper_int32_T(int32_T *data, int32_T rows, int32_T
  cols)
{
  emxArray_int32_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  c_emxInit_int32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_real_T *emxCreateWrapper_real_T(real_T *data, int32_T rows, int32_T
  cols)
{
  emxArray_real_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  c_emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_int32_T *emxCreate_int32_T(int32_T rows, int32_T cols)
{
  emxArray_int32_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  c_emxInit_int32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (int32_T *)calloc((uint32_T)numEl, sizeof(int32_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols)
{
  emxArray_real_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  c_emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

void emxDestroyArray_int32_T(emxArray_int32_T *emxArray)
{
  emxFree_int32_T(&emxArray);
}

void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
}

/* End of code generation (compute_volume_usingtets.c) */
