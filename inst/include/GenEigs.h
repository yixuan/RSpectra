#ifndef GENEIGS_H
#define GENEIGS_H

#ifdef __cplusplus

#ifdef USE_SPECTRA_1YZ
    #include <next/Spectra/GenEigsSolver.h>
    #include <next/Spectra/GenEigsRealShiftSolver.h>
    #include <next/Spectra/GenEigsComplexShiftSolver.h>
#else
    #include <Spectra/GenEigsSolver.h>
    #include <Spectra/GenEigsRealShiftSolver.h>
    #include <Spectra/GenEigsComplexShiftSolver.h>
#endif

#endif

#endif /* GENEIGS_H */
