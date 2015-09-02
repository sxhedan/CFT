
/* Custom Source Code */
#ifdef BUILD_MEX

/* Define macros to support building function into MATLAB executable. */
#include "mex.h"
#define malloc                         mxMalloc
#define calloc                         mxCalloc
#define realloc                        mxRealloc
#define free                           mxFree
#define printf                         mexPrintf
#endif

#define TRILIBRARY
#define NO_TIMER
#define REAL double
#include "triangle.c"

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

void call_triangle(const emxArray_real_T *Fps, const emxArray_int32_T *Elems,
        emxArray_real_T *Fps_new, emxArray_int32_T *Ftris)
{
    int32_T m1 = Fps->size[0];
    int32_T m2 = Elems->size[0];
    int32_T i;
    int32_T j;
    
    /* Declare and Initialize two 1D arrays : nodes and edges which should be passed as the input to "triangle" */
    real_T  *nodes = (real_T*)calloc( 2*m1, sizeof(real_T));
    int32_T *edges = (int32_T*)calloc( 2*m2, sizeof(int32_T));
    
    int32_T out_m1, out_m2;
    
    /* Define the I/O structure for "triangle"*/
    struct triangulateio in;
    struct triangulateio out;
    
    /* Change 2D matrices Fps and Elems to 1D arrays nodes and edges */
    for (i = 0; i<m1; i++){
        for (j = 0; j<2; j++){
            nodes[2*i+j] = Fps->data[m1*j+i];
        }
    }
    
    for (i = 0; i<m2; i++){
        for (j = 0; j<2; j++){
            edges[2*i+j] = Elems->data[m2*j+i];
        }
    }
    
    memset( &in, 0, sizeof(struct triangulateio));
    memset( &out, 0, sizeof(struct triangulateio));
    
    in.pointlist = nodes;
    in.numberofpoints = m1;
    in.segmentlist = edges;
    in.numberofsegments = m2;
    
    triangulate_l("pYQ", &in, &out, NULL);
    
    free(nodes);
    free(edges);
    
    /* Collect the output */   
    out_m1 = out.numberofpoints;
    out_m2 = out.numberoftriangles;
    
    Fps_new->size[0] = out_m1;
    Fps_new->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)Fps_new, 0, (int32_T)sizeof(real_T));
        
    for (i=0; i<out.numberofpoints; i++){
        Fps_new->data[i] = out.pointlist[2*i];
        Fps_new->data[out.numberofpoints + i] = out.pointlist[2*i + 1];
    }
    
    Ftris->size[0] = out_m2;
    Ftris->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)Ftris, 0, (int32_T)sizeof(int32_T));  
    
    for (i=0; i<out.numberoftriangles; i++){
        Ftris->data[i] = out.trianglelist[3*i];
        Ftris->data[out.numberoftriangles + i] = out.trianglelist[3*i + 1];
        Ftris->data[2*out.numberoftriangles +i] = out.trianglelist[3*i + 2];
    }
}
