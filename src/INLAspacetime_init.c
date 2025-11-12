
/* INLAspacetime_init.c
 *
 * Copyright (C) 2024 Elias T Krainski
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Elias Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */

#include "INLAspacetime.h"
#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[] = {
	{"ar2covk", (DL_FUNC) &ar2covk, 5},
	{"CSphere_gamma_alpha", (DL_FUNC) &CSphere_gamma_alpha, 3},
	{NULL, NULL, 0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
    __attribute__((visibility("default")))
#endif
    R_init_INLAspacetime(DllInfo *dll)
{
	R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
