
#include "compute_volume_usingtets.h"

void volumecompfnc(const double cell_bndsIN[6], 
        int face_typeIN[6], unsigned char face_signIN[6], 
        const double *psIN, int num_pointsIN,
        const int *trisIN, int num_trisIN, 
        const int *comp_indexIN, int num_compsIN, 
        int *comp_listIN, int *pos_mapIN, int *neg_mapIN, double *pos_volIN, double *neg_volIN)
{
    int i;
    
    /* Process the given input into a target function compatible form */
    emxArray_real_T *ps;
    emxArray_int32_T *tris;
    emxArray_int32_T *comp_index;
    emxArray_int32_T *comp_list;
    emxArray_int32_T *pos_map;
    emxArray_int32_T *neg_map;
    emxArray_real_T *pos_vol;
    emxArray_real_T *neg_vol;
    
    /*    const real_T *c_centers = c_centersIN;
    const real_T *grid_spacing = grid_spacingIN;*/
    const real_T *cell_bnds = cell_bndsIN;
    int32_T *face_type = face_typeIN ;
    boolean_T *face_sign= face_signIN;
        
    int32_T num_points = num_pointsIN;
    int32_T num_tris = num_trisIN;
    int32_T num_comps = num_compsIN;
    
    /* Allocate memory for emxArrays */ 
    ps = emxCreate_real_T(num_pointsIN, 3);
    tris = emxCreate_int32_T(num_trisIN, 3);

    comp_index = emxCreateWrapperND_int32_T((int*)comp_indexIN, 1, &num_trisIN);
    comp_list  = emxCreateWrapperND_int32_T(comp_listIN, 1, &num_compsIN);
    pos_map = emxCreateWrapperND_int32_T(pos_mapIN, 1, &num_compsIN);
    neg_map = emxCreateWrapperND_int32_T(neg_mapIN, 1, &num_compsIN);
    pos_vol = emxCreateWrapperND_real_T(pos_volIN, 1, &num_compsIN);
    neg_vol = emxCreateWrapperND_real_T(neg_volIN, 1, &num_compsIN);
    
    /* Copy data into points and triangles */
    for (i=0;i<num_points;i++) {
      ps->data[i]=psIN[i*3];
      ps->data[num_points+i]=psIN[i*3+1];
      ps->data[2*num_points+i]=psIN[i*3+2];
    }

    for (i=0;i<num_tris;i++) {
      tris->data[i]=trisIN[3*i]+1;
      tris->data[i+num_tris]=trisIN[3*i+1]+1;
      tris->data[i+2*num_tris]=trisIN[3*i+2]+1;      
    }

    /* Invoke the target function */
    compute_volume_usingtets_initialize();
    compute_volume_usingtets(cell_bnds, face_type, face_sign, ps, num_points, tris, num_tris, comp_index, num_comps, comp_list, pos_map, neg_map, pos_vol, neg_vol);
    compute_volume_usingtets_terminate();
    
    /* Free memory */
    emxDestroyArray_real_T(ps);
    emxDestroyArray_int32_T(tris);
    emxDestroyArray_int32_T(comp_index);
    emxDestroyArray_int32_T(comp_list);
    emxDestroyArray_int32_T(pos_map);
    emxDestroyArray_int32_T(neg_map);
    emxDestroyArray_real_T(pos_vol);
    emxDestroyArray_real_T(neg_vol);
}
