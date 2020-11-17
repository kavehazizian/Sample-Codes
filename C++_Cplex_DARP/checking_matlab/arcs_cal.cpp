//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: arcs_cal.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 16-Sep-2016 12:12:10
//

// Include Files
#include "rt_nonfinite.h"
#include "arcs_cal.h"
#include "arcs_cal_emxutil.h"

// Function Definitions

//
// Arguments    : short n
//                emxArray_real32_T *arcs_set
// Return Type  : void
//
void arcs_cal(short n, emxArray_real32_T *arcs_set)
{
  short i0;
  int i1;
  short c;
  int j;
  int counter_arc;
  short ip;
  short b_c;
  short jp;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  float counterp;
  int i;
  float b_arcs_set[2];
  emxArray_real32_T *x;
  emxArray_real32_T *b_x;
  int b_i;
  i0 = n;
  if (i0 > 8191) {
    i0 = MAX_int16_T;
  } else if (i0 <= -8192) {
    i0 = MIN_int16_T;
  } else {
    i0 = (short)(i0 << 2);
  }

  i1 = i0 * n;
  if (i1 > 32767) {
    i1 = 32767;
  } else {
    if (i1 < -32768) {
      i1 = -32768;
    }
  }

  i1 -= n;
  if (i1 > 32767) {
    i1 = 32767;
  } else {
    if (i1 < -32768) {
      i1 = -32768;
    }
  }

  c = (short)i1;
  i1 = arcs_set->size[0] * arcs_set->size[1];
  arcs_set->size[0] = c;
  arcs_set->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)arcs_set, i1, (int)sizeof(float));
  j = c << 1;
  for (i1 = 0; i1 < j; i1++) {
    arcs_set->data[i1] = 0.0F;
  }

  counter_arc = -1;
  i0 = n;
  if (i0 > 16383) {
    i0 = MAX_int16_T;
  } else if (i0 <= -16384) {
    i0 = MIN_int16_T;
  } else {
    i0 <<= 1;
  }

  i1 = i0 + 2;
  if (i1 > 32767) {
    i1 = 32767;
  }

  c = (short)i1;
  for (ip = 1; ip <= c; ip++) {
    i0 = n;
    if (i0 > 16383) {
      i0 = MAX_int16_T;
    } else if (i0 <= -16384) {
      i0 = MIN_int16_T;
    } else {
      i0 <<= 1;
    }

    i1 = i0 + 2;
    if (i1 > 32767) {
      i1 = 32767;
    }

    b_c = (short)i1;
    for (jp = 1; jp <= b_c; jp++) {
      guard1 = false;
      guard2 = false;
      guard3 = false;
      if (ip == jp) {
        //          elseif jp==2*n+2
        //              countp=counter_arc+(n-(ip-n-1))*(2*n-2)+(ip-n-1-1);
        //                  arcs_set(countp,1:2)=[ip-1 jp-1];
      } else {
        i1 = jp + n;
        if (i1 > 32767) {
          i1 = 32767;
        } else {
          if (i1 < -32768) {
            i1 = -32768;
          }
        }

        if (ip == i1) {
          //          elseif jp==2*n+2
          //              countp=counter_arc+(n-(ip-n-1))*(2*n-2)+(ip-n-1-1);
          //                  arcs_set(countp,1:2)=[ip-1 jp-1];
        } else if (ip == 1) {
          i1 = n + 1;
          if (i1 > 32767) {
            i1 = 32767;
          }

          if (jp > i1) {
            //          elseif jp==2*n+2
            //              countp=counter_arc+(n-(ip-n-1))*(2*n-2)+(ip-n-1-1);
            //                  arcs_set(countp,1:2)=[ip-1 jp-1];
          } else {
            guard3 = true;
          }
        } else {
          guard3 = true;
        }
      }

      if (guard3) {
        i0 = n;
        if (i0 > 16383) {
          i0 = MAX_int16_T;
        } else if (i0 <= -16384) {
          i0 = MIN_int16_T;
        } else {
          i0 <<= 1;
        }

        i1 = i0 + 2;
        if (i1 > 32767) {
          i1 = 32767;
        }

        if (jp == i1) {
          i1 = n + 1;
          if (i1 > 32767) {
            i1 = 32767;
          }

          if (ip <= i1) {
            //          elseif jp==2*n+2
            //              countp=counter_arc+(n-(ip-n-1))*(2*n-2)+(ip-n-1-1);
            //                  arcs_set(countp,1:2)=[ip-1 jp-1];
          } else {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }
      }

      if (guard2) {
        if ((jp == 1) && (ip > 1)) {
          //          elseif jp==2*n+2
          //              countp=counter_arc+(n-(ip-n-1))*(2*n-2)+(ip-n-1-1);
          //                  arcs_set(countp,1:2)=[ip-1 jp-1];
        } else {
          i1 = n + 1;
          if (i1 > 32767) {
            i1 = 32767;
          }

          if (ip < i1) {
            i0 = n;
            if (i0 > 16383) {
              i0 = MAX_int16_T;
            } else if (i0 <= -16384) {
              i0 = MIN_int16_T;
            } else {
              i0 <<= 1;
            }

            i1 = i0 + 2;
            if (i1 > 32767) {
              i1 = 32767;
            }

            if (jp == i1) {
              //          elseif jp==2*n+2
              //              countp=counter_arc+(n-(ip-n-1))*(2*n-2)+(ip-n-1-1); 
              //                  arcs_set(countp,1:2)=[ip-1 jp-1];
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }
        }
      }

      if (guard1) {
        i0 = n;
        if (i0 > 16383) {
          i0 = MAX_int16_T;
        } else if (i0 <= -16384) {
          i0 = MIN_int16_T;
        } else {
          i0 <<= 1;
        }

        i1 = i0 + 2;
        if (i1 > 32767) {
          i1 = 32767;
        }

        if (ip == i1) {
          //          elseif jp==2*n+2
          //              countp=counter_arc+(n-(ip-n-1))*(2*n-2)+(ip-n-1-1);
          //                  arcs_set(countp,1:2)=[ip-1 jp-1];
        } else {
          counter_arc++;
          i1 = ip - 1;
          if (i1 < -32768) {
            i1 = -32768;
          }

          arcs_set->data[counter_arc] = (float)i1;
          i1 = jp - 1;
          if (i1 < -32768) {
            i1 = -32768;
          }

          arcs_set->data[counter_arc + arcs_set->size[0]] = (float)i1;
        }
      }
    }
  }

  counterp = 0.0F;
  for (i = 0; i <= counter_arc; i++) {
    i1 = n + 1;
    if (i1 > 32767) {
      i1 = 32767;
    }

    if ((int)arcs_set->data[i] >= i1) {
      i0 = n;
      if (i0 > 16383) {
        i0 = MAX_int16_T;
      } else if (i0 <= -16384) {
        i0 = MIN_int16_T;
      } else {
        i0 <<= 1;
      }

      i1 = i0 + 1;
      if (i1 > 32767) {
        i1 = 32767;
      }

      if ((int)arcs_set->data[i + arcs_set->size[0]] == i1) {
        counterp++;
        j = (int)((float)(counter_arc + 1) + counterp);
        for (i1 = 0; i1 < 2; i1++) {
          b_arcs_set[i1] = arcs_set->data[i + arcs_set->size[0] * i1];
        }

        for (i1 = 0; i1 < 2; i1++) {
          arcs_set->data[(j + arcs_set->size[0] * i1) - 1] = b_arcs_set[i1];
        }

        for (i1 = 0; i1 < 2; i1++) {
          arcs_set->data[i + arcs_set->size[0] * i1] = 0.0F;
        }
      }
    }
  }

  i = 0;
  emxInit_real32_T(&x, 2);
  emxInit_real32_T(&b_x, 2);
  while (i <= counter_arc) {
    if ((arcs_set->data[i] == 0.0F) && (arcs_set->data[i + arcs_set->size[0]] ==
         0.0F)) {
      i1 = x->size[0] * x->size[1];
      x->size[0] = arcs_set->size[0];
      x->size[1] = 2;
      emxEnsureCapacity((emxArray__common *)x, i1, (int)sizeof(float));
      j = arcs_set->size[0] * arcs_set->size[1];
      for (i1 = 0; i1 < j; i1++) {
        x->data[i1] = arcs_set->data[i1];
      }

      for (j = 0; j < 2; j++) {
        for (b_i = i + 1; b_i < arcs_set->size[0]; b_i++) {
          x->data[(b_i + x->size[0] * j) - 1] = x->data[b_i + x->size[0] * j];
        }
      }

      if (1 > arcs_set->size[0] - 1) {
        j = 0;
      } else {
        j = arcs_set->size[0] - 1;
      }

      i1 = b_x->size[0] * b_x->size[1];
      b_x->size[0] = j;
      b_x->size[1] = 2;
      emxEnsureCapacity((emxArray__common *)b_x, i1, (int)sizeof(float));
      for (i1 = 0; i1 < 2; i1++) {
        for (b_i = 0; b_i < j; b_i++) {
          b_x->data[b_i + b_x->size[0] * i1] = x->data[b_i + x->size[0] * i1];
        }
      }

      i1 = x->size[0] * x->size[1];
      x->size[0] = b_x->size[0];
      x->size[1] = 2;
      emxEnsureCapacity((emxArray__common *)x, i1, (int)sizeof(float));
      for (i1 = 0; i1 < 2; i1++) {
        j = b_x->size[0];
        for (b_i = 0; b_i < j; b_i++) {
          x->data[b_i + x->size[0] * i1] = b_x->data[b_i + b_x->size[0] * i1];
        }
      }

      i1 = arcs_set->size[0] * arcs_set->size[1];
      arcs_set->size[0] = x->size[0];
      arcs_set->size[1] = 2;
      emxEnsureCapacity((emxArray__common *)arcs_set, i1, (int)sizeof(float));
      j = x->size[0] * x->size[1];
      for (i1 = 0; i1 < j; i1++) {
        arcs_set->data[i1] = x->data[i1];
      }
    }

    i++;
  }

  emxFree_real32_T(&b_x);
  emxFree_real32_T(&x);
}

//
// File trailer for arcs_cal.cpp
//
// [EOF]
//
