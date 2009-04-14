#include <R.h>
#include <R_ext/Rdynload.h>
#include "flowCore.h"

static const R_CallMethodDef CallEntries[] = {
   // {"biexponential_transform", (DL_FUNC)&biexponential_transform, 9},
     {"logicle_transform", (DL_FUNC)&logicle_transform, 9},
     {"biexponential_transform", (DL_FUNC)&biexponential_transform, 9},
    {"inPolygon", (DL_FUNC)&inPolygon, 2},
    {"inPolytope", (DL_FUNC)&inPolytope, 3},
    {NULL, NULL, 0}
};

void R_init_flowCore(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    //R_useDynamicSymbols(dll, FALSE);
}