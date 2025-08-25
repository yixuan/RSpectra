#ifndef SPECTRAC_H
#define SPECTRAC_H

#ifdef __cplusplus
extern "C" {
#endif

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

/* Options related to Spectra */
typedef struct {
    int    rule;    /* 0-LM, 1-LR, 2-LI, 3-LA, 4-SM, 5-SR, 6-SI, 7-SA, 8-BE  */
    int    ncv;     /* number of Ritz values in iteration */
    double tol;     /* precision parameter */
    int    maxitr;  /* maximum number of iterations */
    int    retvec;  /* 0 - do not return eigenvectors, 1 - return eigenvectors */
} spectra_opts;

/* Function to represent matrix operation */
typedef void (*mat_op)(const double* x_in, double* y_out, int n, void* data);

/*
 * op:    Function pointer of matrix operation (in)
 * n:     Dimension of the matrix              (in)
 * k:     Number of eigenvalues requested      (in)
 * opts:  Additional options                   (in)
 * nconv: Number of converged eigenvalues      (out)
 * niter: Number of iterations used            (out)
 * nops:  Number of matrix operations used     (out)
 * evals: The first nconv elements will be overwritten by converged eigenvalues.
 *        On entry this array must contain at least k elements. (out)
 * evecs: If opts->retvec is nonzero, the first n*nconv elements will be
 *        overwritten by converged eigenvectors column by column.
 *        If opts->retvec is zero, this array will be intact.   (out)
 * info:  Nonzero if error occurs.             (out)
 */
void eigs_sym_c(
    mat_op op, int n, int k,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals, double* evecs, int* info
);

typedef void (*eigs_sym_c_funtype)(
    mat_op op, int n, int k,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals, double* evecs, int* info
);


/*
 * sigma: The shift. (in)
 */
void eigs_sym_shift_c(
    mat_op op, int n, int k, double sigma,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals, double* evecs, int* info
);

typedef void (*eigs_sym_shift_c_funtype)(
    mat_op op, int n, int k, double sigma,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals, double* evecs, int* info
);

/*
 * evals_r: Real part of the eigenvalues.       (out)
 * evals_i: Imaginary part of the eigenvalues.  (out)
 * evecs_r: Real part of the eigenvectors.      (out)
 * evecs_i: Imaginary part of the eigenvectors. (out)
 */
void eigs_gen_c(
    mat_op op, int n, int k,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals_r, double* evals_i, double* evecs_r, double* evecs_i, int* info
);

typedef void (*eigs_gen_c_funtype)(
    mat_op op, int n, int k,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals_r, double* evals_i, double* evecs_r, double* evecs_i, int* info
);

/*
 * sigmar: Real part of the shift.      (in)
 * sigmai: Imaginary part of the shift. (in)
 */
void eigs_gen_shift_c(
    mat_op op, int n, int k, double sigmar, double sigmai,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals_r, double* evals_i, double* evecs_r, double* evecs_i, int* info
);

typedef void (*eigs_gen_shift_c_funtype)(
    mat_op op, int n, int k, double sigmar, double sigmai,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals_r, double* evals_i, double* evecs_r, double* evecs_i, int* info
);


#ifdef __cplusplus
}
#endif

#endif /* SPECTRAC_H */
