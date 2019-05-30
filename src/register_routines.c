#include <SpectraC.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Exported functions */
SEXP eigs_sym(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
);

SEXP eigs_shift_sym(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
);

SEXP eigs_gen(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
);

SEXP eigs_real_shift_gen(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
);

SEXP eigs_complex_shift_gen(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
);

SEXP svds_sym(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP nu_scalar_r, SEXP nv_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
);

SEXP svds_gen(
    SEXP A_mat_r, SEXP m_scalar_r, SEXP n_scalar_r,
    SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
);

SEXP is_sym_dgCMatrix(SEXP mat, SEXP tol);

SEXP is_sym_dgRMatrix(SEXP mat, SEXP tol);

static const R_CallMethodDef CallEntries[] = {
    {"eigs_sym",               (DL_FUNC) &eigs_sym,               5},
    {"eigs_shift_sym",         (DL_FUNC) &eigs_shift_sym,         5},
    {"eigs_gen",               (DL_FUNC) &eigs_gen,               5},
    {"eigs_real_shift_gen",    (DL_FUNC) &eigs_real_shift_gen,    5},
    {"eigs_complex_shift_gen", (DL_FUNC) &eigs_complex_shift_gen, 5},
    {"svds_sym",               (DL_FUNC) &svds_sym,               7},
    {"svds_gen",               (DL_FUNC) &svds_gen,               8},
    {"is_sym_dgCMatrix",       (DL_FUNC) &is_sym_dgCMatrix,       2},
    {"is_sym_dgRMatrix",       (DL_FUNC) &is_sym_dgRMatrix,       2},
    {NULL, NULL, 0}
};

void R_init_RSpectra(DllInfo* info)
{
    /* Register C interface */
    R_RegisterCCallable("RSpectra", "eigs_sym_c",       (DL_FUNC) eigs_sym_c);
    R_RegisterCCallable("RSpectra", "eigs_sym_shift_c", (DL_FUNC) eigs_sym_shift_c);
    R_RegisterCCallable("RSpectra", "eigs_gen_c",       (DL_FUNC) eigs_gen_c);
    R_RegisterCCallable("RSpectra", "eigs_gen_shift_c", (DL_FUNC) eigs_gen_shift_c);

    /* Register R .Call functions */
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
