#ifndef RMATOP_H
#define RMATOP_H

/*
  1. RSpectra 0.17-0 contains both the old Spectra v0.y.z and the new v1.y.z.
  2. For future code that wants to link to RSpectra:
      (a) If the user indeed wants to use Spectra v1.y.z, first define
          a macro called PREFER_SPECTRA_1YZ before including <SymEigs.h>
          (or <GenEigs.h> etc.).
          - If the RSpectra package that is linked to indeed contains v1.y.z,
            it will define a macro called SPECTRA_1YZ_INCLUDED, and the
            user can query this macro to verify the existence of v1.y.z.
          - If the RSpectra package that is linked to is an old one
            (this is the case when RSpectra is in the process of updaing on CRAN),
            then the macro SUPPORT_SPECTRA_1YZ is undefined, and the user can
            decide what to do next based on this status. Typically this means
            that the old API needs to be used.
      (b) If the user does not want to use Spectra v1.y.z, then just use the
          old API, since RSpectra will be back compatible for some time.
          But note that this is discouraged, and eventually RSpectra will
          move to the new interface.
*/

#ifdef PREFER_SPECTRA_1YZ
    /* Indicates that this version of RSpectra contains Spectra 1.y.z */
    #ifndef SPECTRA_1YZ_INCLUDED
    #define SPECTRA_1YZ_INCLUDED
    #endif
#endif

#ifdef __cplusplus

#include "RMatOp/MatProd.h"
#include "RMatOp/MatProd_matrix.h"
#include "RMatOp/MatProd_sym_matrix.h"
#include "RMatOp/MatProd_dgeMatrix.h"
#include "RMatOp/MatProd_sym_dgeMatrix.h"
#include "RMatOp/MatProd_dsyMatrix.h"
#include "RMatOp/MatProd_sparseMatrix.h"
#include "RMatOp/MatProd_sym_sparseMatrix.h"
#include "RMatOp/MatProd_function.h"

#include "RMatOp/RealShift.h"
#include "RMatOp/RealShift_matrix.h"
#include "RMatOp/RealShift_sym_matrix.h"
#include "RMatOp/RealShift_dgeMatrix.h"
#include "RMatOp/RealShift_sym_dgeMatrix.h"
#include "RMatOp/RealShift_dsyMatrix.h"
#include "RMatOp/RealShift_sparseMatrix.h"
#include "RMatOp/RealShift_sym_sparseMatrix.h"

#include "RMatOp/ComplexShift.h"
#include "RMatOp/ComplexShift_matrix.h"
#include "RMatOp/ComplexShift_dgeMatrix.h"
#include "RMatOp/ComplexShift_sparseMatrix.h"

#include "RMatOp/SVDOp.h"

#endif /* __cplusplus */

#endif /* RMATOP_H */
