/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* internal.h - shared helpers not part of the public API */
#ifndef DCREATOR_INTERNAL_H
#define DCREATOR_INTERNAL_H

#include "dcreator.h"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

/* Safe snprintf into errbuf; no-op if errbuf is NULL or zero. */
__attribute__((format(printf, 3, 4)))
static inline void dc_errf(char *errbuf, size_t sz, const char *fmt, ...)
{
    if (!errbuf || sz == 0) return;
    va_list ap;
    va_start(ap, fmt);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
    vsnprintf(errbuf, sz, fmt, ap);
#pragma GCC diagnostic pop
    va_end(ap);
}

/* ---- vector helpers ---- */
static inline double dc_vabs(dc_vec3_t a)
{
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}
static inline double dc_dot(dc_vec3_t a, dc_vec3_t b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static inline dc_vec3_t dc_cross(dc_vec3_t a, dc_vec3_t b)
{
    dc_vec3_t r = { a.y*b.z - a.z*b.y,
                    a.z*b.x - a.x*b.z,
                    a.x*b.y - a.y*b.x };
    return r;
}
static inline dc_vec3_t dc_scale(dc_vec3_t a, double s)
{
    dc_vec3_t r = { a.x*s, a.y*s, a.z*s };
    return r;
}
static inline int dc_is_zero(dc_vec3_t a)
{
    return a.x == 0.0 && a.y == 0.0 && a.z == 0.0;
}
static inline dc_vec3_t dc_normalize(dc_vec3_t a)
{
    double n = dc_vabs(a);
    if (n > 0.0) { a.x /= n; a.y /= n; a.z /= n; }
    return a;
}

/* matrix[row][col] * vec */
static inline dc_vec3_t dc_mv(const double m[3][3], dc_vec3_t v)
{
    dc_vec3_t r = {
        m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
        m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
        m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z,
    };
    return r;
}

static inline double dc_det(const double m[3][3])
{
    return m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
         - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
         + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
}

/* ---- reading helpers ---- */

/* Read whole file into malloc'd buffer, NUL-terminated. Caller frees.
 * Sets *out_len to the number of bytes (not counting the NUL). */
int dc_slurp(const char *filename, char **out_buf, size_t *out_len,
             char *errbuf, size_t errbuf_sz);

/* Trim leading/trailing whitespace in place, returning the new start. */
char *dc_strstrip(char *s);

/* Case-insensitive compare of two strings up to n chars; returns 0 on match. */
int dc_strncasecmp_(const char *a, const char *b, size_t n);

#endif /* DCREATOR_INTERNAL_H */
