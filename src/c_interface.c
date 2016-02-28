#include <SpectraC.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_rARPACK(DllInfo *info)
{
    R_RegisterCCallable("RSpectra", "eigs_sym_c",       (DL_FUNC) eigs_sym_c);
    R_RegisterCCallable("RSpectra", "eigs_sym_shift_c", (DL_FUNC) eigs_sym_shift_c);
    R_RegisterCCallable("RSpectra", "eigs_gen_c",       (DL_FUNC) eigs_gen_c);
    R_RegisterCCallable("RSpectra", "eigs_gen_shift_c", (DL_FUNC) eigs_gen_shift_c);
}

