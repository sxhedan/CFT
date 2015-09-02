/*
 * compute_volume_usingtets.h
 *
 * Code generation for function 'compute_volume_usingtets'
 *
 * C source code generated on: Sat May 12 14:38:40 2012
 *
 */

#ifndef __COMPUTE_VOLUME_USINGTETS_H__
#define __COMPUTE_VOLUME_USINGTETS_H__
/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rtwtypes.h"
#include "compute_volume_usingtets_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void compute_volume_usingtets(const real_T cell_bnds[6], int32_T face_type[6], boolean_T face_sign[6], const emxArray_real_T *ps, int32_T num_points, const emxArray_int32_T *tris, int32_T num_tris, const emxArray_int32_T *comp_index, int32_T num_comps, const emxArray_int32_T *comp_list, emxArray_int32_T *pos_map, emxArray_int32_T *neg_map, emxArray_real_T *pos_vol, emxArray_real_T *neg_vol);
extern void compute_volume_usingtets_initialize(void);
extern void compute_volume_usingtets_terminate(void);
extern emxArray_int32_T *emxCreateND_int32_T(int32_T numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size);
extern emxArray_int32_T *emxCreateWrapperND_int32_T(int32_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_int32_T *emxCreateWrapper_int32_T(int32_T *data, int32_T rows, int32_T cols);
extern emxArray_real_T *emxCreateWrapper_real_T(real_T *data, int32_T rows, int32_T cols);
extern emxArray_int32_T *emxCreate_int32_T(int32_T rows, int32_T cols);
extern emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_int32_T(emxArray_int32_T *emxArray);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
#endif
/* End of code generation (compute_volume_usingtets.h) */
