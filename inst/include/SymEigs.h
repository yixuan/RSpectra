#ifndef SYMEIGS_H
#define SYMEIGS_H

#ifdef __cplusplus

#ifdef USE_SPECTRA_1YZ
    #include <next/Spectra/SymEigsSolver.h>
    #include <next/Spectra/SymEigsShiftSolver.h>

    #include <next/Spectra/MatOp/DenseSymMatProd.h>
    #include <next/Spectra/MatOp/SparseSymMatProd.h>
    #include <next/Spectra/MatOp/DenseSymShiftSolve.h>
    #include <next/Spectra/MatOp/SparseSymShiftSolve.h>
#else
    #include <Spectra/SymEigsSolver.h>
    #include <Spectra/SymEigsShiftSolver.h>
#endif

#endif

#endif /* SYMEIGS_H */
