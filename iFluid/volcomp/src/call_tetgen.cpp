#define TETLIBRARY
#define SELF_CHECK
#include "tetgen.cpp"
//#include "tetgen.cxx"
#include "rtwtypes.h"

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

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T
typedef struct emxArray_int32_T
{
    int32_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray_int32_T;
#endif

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
typedef struct emxArray_real_T
{
    real_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray_real_T;
#endif

static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize)
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

extern "C" void call_tetgen(const emxArray_real_T *ps_local, const emxArray_int32_T
        *tris_local, emxArray_real_T *ps_vol, emxArray_int32_T *tets)
{
    /* INPUT DECLARATION */
    int32_T nv;
    int32_T ne;
    int32_T i;
    int32_T j;
    
    nv = ps_local->size[0];
    ne = tris_local->size[0];
    
    tetgenio in, out;
    in.firstnumber = 1;
    
    /* VERTICES */
    in.numberofpoints = nv;
    in.pointlist = new REAL[in.numberofpoints * 3];
    
    for (i = 0; i<nv; i++){
        for (j = 0; j<3; j++){
            in.pointlist[3*i+j] = ps_local->data[nv*j+i];
        }
    }
    
    /* TRIANGLE */
    in.numberoffacets = ne;
    in.facetlist = new tetgenio::facet[ne];
    
    for (i=0; i<ne; i++){
        tetgenio::facet *f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[1];
        f->numberofholes = 0;
        f->holelist = NULL;
        tetgenio::polygon *p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[3];
        p->vertexlist[0] = tris_local->data[i];
        p->vertexlist[1] = tris_local->data[ne + i];
        p->vertexlist[2] = tris_local->data[ne*2 + i];
    }
    
    /* CALL TETGEN */
    const char *switches = "pYQ";
    tetrahedralize((char*)switches, &in, &out);
    
    /* OUTPUT */
    /* Obtain nodal coordinates*/
    int32_T out_nv ;
    int32_T out_ne ;
    out_nv = out.numberofpoints;
    out_ne = out.numberoftetrahedra;
    
    ps_vol->size[0] = out_nv;
    ps_vol->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)ps_vol, 0, (real_T)sizeof(real_T));
    
    for (i=0; i<out.numberofpoints; i++) {
        ps_vol->data[i] = out.pointlist[3*i];
        ps_vol->data[out.numberofpoints + i] = out.pointlist[3*i + 1];
        ps_vol->data[2*out.numberofpoints + i] = out.pointlist[3*i + 2];
    }
    
    tets->size[0] = out_ne;
    tets->size[1] = 4;
    emxEnsureCapacity((emxArray__common *)tets, 0, (int32_T)sizeof(int32_T));
       
    for (i=0; i<out.numberoftetrahedra; i++){
        tets->data[i] = out.tetrahedronlist[4*i];
        tets->data[out.numberoftetrahedra + i] = out.tetrahedronlist[4*i + 1];
        tets->data[2*out.numberoftetrahedra +i] = out.tetrahedronlist[4*i + 2];
        tets->data[3*out.numberoftetrahedra +i] = out.tetrahedronlist[4*i + 3];
    }
}
