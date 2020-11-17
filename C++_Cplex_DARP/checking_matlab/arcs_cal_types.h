//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: arcs_cal_types.h
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 16-Sep-2016 12:12:10
//
#ifndef __ARCS_CAL_TYPES_H__
#define __ARCS_CAL_TYPES_H__

// Include Files
#include "rtwtypes.h"

// Type Definitions
#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
  void *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray__common

#ifndef struct_emxArray_real32_T
#define struct_emxArray_real32_T

struct emxArray_real32_T
{
  float *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_real32_T
#endif

//
// File trailer for arcs_cal_types.h
//
// [EOF]
//
