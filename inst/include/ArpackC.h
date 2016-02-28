#ifndef ARPACKC_H
#define ARPACKC_H

#ifdef __cplusplus
extern "C" {
#endif



/* Options related to ARPACK */
typedef struct {
    int rule;     /* 0-LM, 1-LR, 2-LI, 3-LA, 4-SM, 5-SR, 6-SI, 7-SA, 8-BE  */
    int ncv;      /* number of Ritz values in iteration */
    double tol;   /* precision parameter */
    int maxitr;   /* maximum number of iterations */
    int retvec;   /* 0 - do not return eigenvectors, 1 - return eigenvectors */
} arpack_opts;

/* Function to represent matrix operation */
typedef void (*mat_op)(double *x_in, double *y_out, int n, void *data);

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
    const arpack_opts *opts, void *data,
    int *nconv, int *niter, int *nops,
    double *evals, double *evecs, int *info
);

typedef void (*eigs_sym_c_funtype)(
        mat_op op, int n, int k,
        const arpack_opts *opts, void *data,
        int *nconv, int *niter, int *nops,
        double *evals, double *evecs, int *info
);


/*
 * sigma: The shift. (in)
 */
void eigs_sym_shift_c(
    mat_op op, int n, int k, double sigma,
    const arpack_opts *opts, void *data,
    int *nconv, int *niter, int *nops,
    double *evals, double *evecs, int *info
);

typedef void (*eigs_sym_shift_c_funtype)(
        mat_op op, int n, int k, double sigma,
        const arpack_opts *opts, void *data,
        int *nconv, int *niter, int *nops,
        double *evals, double *evecs, int *info
);

/*
 * evals_r: Real part of the eigenvalues.       (out)
 * evals_i: Imaginary part of the eigenvalues.  (out)
 * evecs_r: Real part of the eigenvectors.      (out)
 * evecs_i: Imaginary part of the eigenvectors. (out)
 */
void eigs_gen_c(
    mat_op op, int n, int k,
    const arpack_opts *opts, void *data,
    int *nconv, int *niter, int *nops,
    double *evals_r, double *evals_i, double *evecs_r, double *evecs_i, int *info
);

typedef void (*eigs_gen_c_funtype)(
        mat_op op, int n, int k,
        const arpack_opts *opts, void *data,
        int *nconv, int *niter, int *nops,
        double *evals_r, double *evals_i, double *evecs_r, double *evecs_i, int *info
);

/*
 * sigmar: Real part of the shift.      (in)
 * sigmai: Imaginary part of the shift. (in)
 */
void eigs_gen_shift_c(
    mat_op op, int n, int k, double sigmar, double sigmai,
    const arpack_opts *opts, void *data,
    int *nconv, int *niter, int *nops,
    double *evals_r, double *evals_i, double *evecs_r, double *evecs_i, int *info
);

typedef void (*eigs_gen_shift_c_funtype)(
        mat_op op, int n, int k, double sigmar, double sigmai,
        const arpack_opts *opts, void *data,
        int *nconv, int *niter, int *nops,
        double *evals_r, double *evals_i, double *evecs_r, double *evecs_i, int *info
);


#ifdef __cplusplus
}
#endif

#endif /* ARPACKC_H */
