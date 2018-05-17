/*! \file GaIn.h
    \brief function declarations for linking with libGaIn
*/
#ifndef LIB_GAIN_H
#define LIB_GAIN_H

#include "GaIn.config.h"

#define C_Overlap_C     FC_FUNC_(c_overlap_c)
#define C_Laplacian_C   FC_FUNC_(c_laplacian_c)
#define CC_Coulomb_Ion  FC_FUNC_(cc_coulomb_ion)
#define CC_Coulomb_CC   FC_FUNC_(cc_coulomb_cc)
#define CC_Coulomb_C    FC_FUNC_(cc_coulomb_c)
#define C_Coulomb_C     FC_FUNC_(c_coulomb_c)
#define Y_Overlap_Y     FC_FUNC_(y_overlap_y)
#define Y_Laplacian_Y   FC_FUNC_(y_laplacian_y)
#define YY_Coulomb_Ion  FC_FUNC_(yy_coulomb_ion)
#define YY_Coulomb_YY   FC_FUNC_(yy_coulomb_yy)
#define YY_Coulomb_Y    FC_FUNC_(yy_coulomb_y)
#define Y_Coulomb_Y     FC_FUNC_(y_coulomb_y)

#ifdef __cplusplus
extern "C"
  {
#endif
      double C_Overlap_C   (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2);
      double CC_Overlap_C  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2,
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double C_Laplacian_C (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2);
      double C_Coulomb_C   (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2);
      double CC_Coulomb_C  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double CC_Coulomb_CC (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1,
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3, 
                            const double *a4, const double *r4, const int *nx4, const int *ny4, const int *nz4);
      double CC_Coulomb_Ion(const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                            const double *rion);
      double Y_Value       (const double *a1, const int *l1, const int *m1, const double *r1);
      double Y_Overlap_Y   (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2);
      double YY_Overlap_Y  (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2,
                            const double *a3, const double *r3, const int *l3, const int *m3);
      double Y_Laplacian_Y (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2);
      double Y_Coulomb_Y   (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2);
      double YY_Coulomb_Y  (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                            const double *a3, const double *r3, const int *l3, const int *m3);
      double YY_Coulomb_YY (const double *a1, const double *r1, const int *l1, const int *m1,
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                            const double *a3, const double *r3, const int *l3, const int *m3, 
                            const double *a4, const double *r4, const int *l4, const int *m4);
      double YY_Coulomb_Ion(const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                            const double *rion);
#ifdef __cplusplus
  }
#endif

#endif
