#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void dtweediep1(void *, void *, void *, void *, void *);
extern void predictExtrapDown(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void predictExtrapUp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void predictInterp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void truncatedRat(void *, void *, void *, void *, void *, void *);
extern void truncatedWalk(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"dtweediep1",        (DL_FUNC) &dtweediep1,         5},
    {"predictExtrapDown", (DL_FUNC) &predictExtrapDown, 10},
    {"predictExtrapUp",   (DL_FUNC) &predictExtrapUp,   10},
    {"predictInterp",     (DL_FUNC) &predictInterp,     11},
    {"truncatedRat",      (DL_FUNC) &truncatedRat,       6},
    {"truncatedWalk",     (DL_FUNC) &truncatedWalk,      5},
    {NULL, NULL, 0}
};

void R_init_WaverideR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}