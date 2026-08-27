#include "kdml.h"

#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

static const R_CallMethodDef CallEntries[] = {
    {"kdml_mcmc_call", (DL_FUNC) &kdml_mcmc_call, 1},
    {"kdml_score_call", (DL_FUNC) &kdml_score_call, 1},
    {"kdml_distance_call", (DL_FUNC) &kdml_distance_call, 1},
    {"kdml_cuda_info_call", (DL_FUNC) &kdml_cuda_info_call, 0},
    {NULL, NULL, 0}
};

void attribute_visible R_init_kdml(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
