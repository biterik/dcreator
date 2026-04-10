/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* atoms.c - dc_atoms_t lifecycle */

#include "internal.h"
#include <stdlib.h>
#include <string.h>

void dc_atoms_init(dc_atoms_t *a)
{
    memset(a, 0, sizeof(*a));
    a->periodic[0] = a->periodic[1] = a->periodic[2] = 1;
    a->style = DC_ATOM_STYLE_ATOMIC;
    a->source_format = DC_FORMAT_UNKNOWN;
}

void dc_atoms_free(dc_atoms_t *a)
{
    if (!a) return;
    free(a->id);
    free(a->type);
    free(a->x);
    free(a->y);
    free(a->z);
    free(a->mass);
    free(a->q);
    free(a->molecule);
    free(a->vx);
    free(a->vy);
    free(a->vz);
    if (a->extra) {
        for (size_t i = 0; i < a->n; ++i) free(a->extra[i]);
        free(a->extra);
    }
    free(a->header_blob);
    free(a->trailer_blob);
    dc_atoms_init(a);
}

static void *xrealloc(void *p, size_t sz)
{
    return realloc(p, sz);
}

dc_status_t dc_atoms_reserve(dc_atoms_t *a, size_t n)
{
    if (n <= a->cap) return DC_OK;

    int64_t *id  = (int64_t *)xrealloc(a->id,   n * sizeof(int64_t));
    int32_t *ty  = (int32_t *)xrealloc(a->type, n * sizeof(int32_t));
    double  *x   = (double  *)xrealloc(a->x,    n * sizeof(double));
    double  *y   = (double  *)xrealloc(a->y,    n * sizeof(double));
    double  *z   = (double  *)xrealloc(a->z,    n * sizeof(double));

    if (!id || !ty || !x || !y || !z) {
        /* On failure the originals are still valid; drop new ones only. */
        free(id != a->id ? id : NULL);
        free(ty != a->type ? ty : NULL);
        free(x  != a->x ? x : NULL);
        free(y  != a->y ? y : NULL);
        free(z  != a->z ? z : NULL);
        return DC_ERR_MEMORY;
    }
    a->id = id; a->type = ty; a->x = x; a->y = y; a->z = z;

    /* Optional columns: grow only if already allocated. */
    if (a->mass)     a->mass     = (double *)xrealloc(a->mass,     n * sizeof(double));
    if (a->q)        a->q        = (double *)xrealloc(a->q,        n * sizeof(double));
    if (a->molecule) a->molecule = (int64_t *)xrealloc(a->molecule, n * sizeof(int64_t));
    if (a->vx)       a->vx       = (double *)xrealloc(a->vx,       n * sizeof(double));
    if (a->vy)       a->vy       = (double *)xrealloc(a->vy,       n * sizeof(double));
    if (a->vz)       a->vz       = (double *)xrealloc(a->vz,       n * sizeof(double));
    if (a->extra)    a->extra    = (char  **)xrealloc(a->extra,    n * sizeof(char *));

    a->cap = n;
    return DC_OK;
}
