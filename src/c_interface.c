#include <ArpackC.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_rARPACK(DllInfo *info)
{
    R_RegisterCCallable("rARPACK", "eigs_sym_c",       (DL_FUNC) eigs_sym_c);
    R_RegisterCCallable("rARPACK", "eigs_sym_shift_c", (DL_FUNC) eigs_sym_shift_c);
    R_RegisterCCallable("rARPACK", "eigs_gen_c",       (DL_FUNC) eigs_gen_c);
    R_RegisterCCallable("rARPACK", "eigs_gen_shift_c", (DL_FUNC) eigs_gen_shift_c);
}

