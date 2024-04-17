#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void ar2cor(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"ar2cor", (DL_FUNC) &ar2cor, 5},
  {NULL, NULL, 0}
};

void R_init_INLAspacetime(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
