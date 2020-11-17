//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: arcs_cal_emxAPI.h
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 16-Sep-2016 12:12:10
//
#ifndef __ARCS_CAL_EMXAPI_H__
#define __ARCS_CAL_EMXAPI_H__

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "arcs_cal_types.h"

// Function Declarations
extern emxArray_real32_T *emxCreateND_real32_T(int numDimensions, int *size);
extern emxArray_real32_T *emxCreateWrapperND_real32_T(float *data, int
  numDimensions, int *size);
extern emxArray_real32_T *emxCreateWrapper_real32_T(float *data, int rows, int
  cols);
extern emxArray_real32_T *emxCreate_real32_T(int rows, int cols);
extern void emxDestroyArray_real32_T(emxArray_real32_T *emxArray);
extern void emxInitArray_real32_T(emxArray_real32_T **pEmxArray, int
  numDimensions);

#endif

//
// File trailer for arcs_cal_emxAPI.h
//
// [EOF]
//
