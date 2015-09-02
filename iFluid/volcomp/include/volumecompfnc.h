
#ifndef __VOLUMECOMPFNC_H__
#define __VOLUMECOMPFNC_H__

/* Only define EXTERN_C if it hasn't been defined already. This allows
 * individual modules to have more control over managing their exports.
 */
#ifndef EXTERN_C

#ifdef __cplusplus
  #define EXTERN_C extern "C"
#else
  #define EXTERN_C extern
#endif

#endif

EXTERN_C void volumecompfnc(const double *cell_bnds,
        int *face_type, unsigned char *face_sign, double *ps,
        int num_points, int *tris,
        int num_tris, const int *comp_index, int num_comps,
        int *comp_list, int *pos_map, int *neg_map,double *pos_vol, double *neg_vol);


#endif
