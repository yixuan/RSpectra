#ifndef SYMEIGS_H
#define SYMEIGS_H

#ifdef __cplusplus

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
    #include <next/Spectra/SymEigsSolver.h>
    #include <next/Spectra/SymEigsShiftSolver.h>

    #include <next/Spectra/MatOp/DenseSymMatProd.h>
    #include <next/Spectra/MatOp/SparseSymMatProd.h>
    #include <next/Spectra/MatOp/DenseSymShiftSolve.h>
    #include <next/Spectra/MatOp/SparseSymShiftSolve.h>

    /* Indicates that this version of RSpectra contains Spectra 1.y.z */
    #ifndef SPECTRA_1YZ_INCLUDED
    #define SPECTRA_1YZ_INCLUDED
    #endif
#else
    #include <Spectra/SymEigsSolver.h>
    #include <Spectra/SymEigsShiftSolver.h>
#endif

#endif /* __cplusplus */

#endif /* SYMEIGS_H */
