#ifndef GENEIGS_H
#define GENEIGS_H

#ifdef __cplusplus

#ifdef USE_SPECTRA_1YZ
    #include <next/Spectra/GenEigsSolver.h>
    #include <next/Spectra/GenEigsRealShiftSolver.h>
    #include <next/Spectra/GenEigsComplexShiftSolver.h>

    #include <next/Spectra/MatOp/DenseGenMatProd.h>
    #include <next/Spectra/MatOp/SparseGenMatProd.h>
    #include <next/Spectra/MatOp/DenseGenRealShiftSolve.h>
    #include <next/Spectra/MatOp/SparseGenRealShiftSolve.h>
    #include <next/Spectra/MatOp/DenseGenComplexShiftSolve.h>
    #include <next/Spectra/MatOp/SparseGenComplexShiftSolve.h>
#else
    #include <Spectra/GenEigsSolver.h>
    #include <Spectra/GenEigsRealShiftSolver.h>
    #include <Spectra/GenEigsComplexShiftSolver.h>
#endif

#endif

#endif /* GENEIGS_H */
