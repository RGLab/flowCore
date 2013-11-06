#include <R.h>
#include <R_ext/Rdynload.h>
#include "flowCore.h"

static const R_CallMethodDef CallEntries[] = {
   // {"biexponential_transform", (DL_FUNC)&biexponential_transform, 9},
     {"logicle_transform", (DL_FUNC)&logicle_transform, 5},
     {"invLogicle_transform",(DL_FUNC)&invLogicle_transform, 5},
     {"hyperlog_transform",(DL_FUNC)&hyperlog_transform, 5},
     {"invHyperlog_transform",(DL_FUNC)&invHyperlog_transform, 5},
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
