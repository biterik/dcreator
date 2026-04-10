/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* lammps_dump.c - reader for LAMMPS dump (text) files.
 *
 * Only the first snapshot in the file is consumed.
 *
 * Supported column names on "ITEM: ATOMS ..." line:
 *   id, type, x, y, z            (unscaled)
 *   xu, yu, zu                   (unwrapped)
 *   xs, ys, zs                   (scaled; converted to unscaled on read)
 *   ix, iy, iz                   (image flags; preserved as trailing text)
 *   anything else                preserved as trailing text
 *
 * Writing LAMMPS dump is not implemented in v1 — output is lammps_data or imd.
 */

#include "../core/internal.h"
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

enum col_kind {
    COL_IGNORE = 0,
    COL_ID,
    COL_TYPE,
    COL_X,  COL_Y,  COL_Z,
    COL_XS, COL_YS, COL_ZS,
    COL_XU, COL_YU, COL_ZU,
};

typedef struct {
    int kind;
    char name[32];
} col_spec_t;

static int classify_col(const char *name)
{
    if (strcmp(name, "id") == 0)      return COL_ID;
    if (strcmp(name, "type") == 0)    return COL_TYPE;
    if (strcmp(name, "x") == 0)       return COL_X;
    if (strcmp(name, "y") == 0)       return COL_Y;
    if (strcmp(name, "z") == 0)       return COL_Z;
    if (strcmp(name, "xu") == 0)      return COL_XU;
    if (strcmp(name, "yu") == 0)      return COL_YU;
    if (strcmp(name, "zu") == 0)      return COL_ZU;
    if (strcmp(name, "xs") == 0)      return COL_XS;
    if (strcmp(name, "ys") == 0)      return COL_YS;
    if (strcmp(name, "zs") == 0)      return COL_ZS;
    return COL_IGNORE;
}

static int fgets_line(FILE *f, char *buf, size_t sz)
{
    if (!fgets(buf, (int)sz, f)) return 0;
    size_t n = strlen(buf);
    while (n > 0 && (buf[n-1] == '\n' || buf[n-1] == '\r')) buf[--n] = '\0';
    return 1;
}

dc_status_t dc_read_lammps_dump(const char *filename, dc_atoms_t *a,
                                char *errbuf, size_t errbuf_sz)
{
    FILE *f = fopen(filename, "r");
    if (!f) {
        dc_errf(errbuf, errbuf_sz, "cannot open '%s': %s",
                filename, strerror(errno));
        return DC_ERR_IO;
    }

    dc_atoms_init(a);
    a->source_format = DC_FORMAT_LAMMPS_DUMP;

    char line[8192];

    /* ITEM: TIMESTEP */
    if (!fgets_line(f, line, sizeof(line)) ||
        strncmp(line, "ITEM: TIMESTEP", 14) != 0) {
        dc_errf(errbuf, errbuf_sz, "LAMMPS dump: missing ITEM: TIMESTEP");
        fclose(f);
        return DC_ERR_FORMAT;
    }
    if (!fgets_line(f, line, sizeof(line))) { fclose(f); return DC_ERR_FORMAT; }
    /* timestep ignored */

    /* ITEM: NUMBER OF ATOMS */
    if (!fgets_line(f, line, sizeof(line)) ||
        strncmp(line, "ITEM: NUMBER OF ATOMS", 21) != 0) {
        dc_errf(errbuf, errbuf_sz, "LAMMPS dump: missing ITEM: NUMBER OF ATOMS");
        fclose(f);
        return DC_ERR_FORMAT;
    }
    size_t n_atoms;
    if (!fgets_line(f, line, sizeof(line)) ||
        sscanf(line, "%zu", &n_atoms) != 1) {
        dc_errf(errbuf, errbuf_sz, "LAMMPS dump: malformed atom count");
        fclose(f);
        return DC_ERR_FORMAT;
    }

    /* ITEM: BOX BOUNDS ... */
    if (!fgets_line(f, line, sizeof(line)) ||
        strncmp(line, "ITEM: BOX BOUNDS", 16) != 0) {
        dc_errf(errbuf, errbuf_sz, "LAMMPS dump: missing ITEM: BOX BOUNDS");
        fclose(f);
        return DC_ERR_FORMAT;
    }
    int triclinic = (strstr(line, "xy") != NULL);

    double b[3][3] = {{0}};
    int cols_box = triclinic ? 3 : 2;
    for (int i = 0; i < 3; ++i) {
        if (!fgets_line(f, line, sizeof(line))) {
            dc_errf(errbuf, errbuf_sz, "LAMMPS dump: short BOX BOUNDS");
            fclose(f);
            return DC_ERR_FORMAT;
        }
        int got = sscanf(line, "%lf %lf %lf", &b[i][0], &b[i][1], &b[i][2]);
        if (got < cols_box) {
            dc_errf(errbuf, errbuf_sz, "LAMMPS dump: bad BOX BOUNDS row %d", i);
            fclose(f);
            return DC_ERR_FORMAT;
        }
    }
    if (triclinic) {
        /* LAMMPS dump encodes BOX BOUNDS for triclinic as:
         *   xlo_bound xhi_bound xy
         *   ylo_bound yhi_bound xz
         *   zlo_bound zhi_bound yz
         * where {xlo,xhi} = {xlo_bound - MIN(0,xy,xz,xy+xz), xhi_bound - MAX(0,xy,xz,xy+xz)}
         * For v1 we store xlo_bound/xhi_bound directly. */
        a->box.xlo = b[0][0]; a->box.xhi = b[0][1];
        a->box.ylo = b[1][0]; a->box.yhi = b[1][1];
        a->box.zlo = b[2][0]; a->box.zhi = b[2][1];
        a->box.xy  = b[0][2];
        a->box.xz  = b[1][2];
        a->box.yz  = b[2][2];
        a->box.triclinic = 1;
    } else {
        a->box.xlo = b[0][0]; a->box.xhi = b[0][1];
        a->box.ylo = b[1][0]; a->box.yhi = b[1][1];
        a->box.zlo = b[2][0]; a->box.zhi = b[2][1];
    }

    /* ITEM: ATOMS <column names...> */
    if (!fgets_line(f, line, sizeof(line)) ||
        strncmp(line, "ITEM: ATOMS", 11) != 0) {
        dc_errf(errbuf, errbuf_sz, "LAMMPS dump: missing ITEM: ATOMS");
        fclose(f);
        return DC_ERR_FORMAT;
    }
    /* Parse column names after "ITEM: ATOMS ". */
    col_spec_t cols[64];
    int ncols = 0;
    int id_col = -1, type_col = -1;
    int x_col = -1, y_col = -1, z_col = -1;
    int scale_mode = 0;   /* 0 unscaled, 1 scaled */

    char *save = NULL;
    char *rest = line + 11;
    char *tok = strtok_r(rest, " \t", &save);
    while (tok && ncols < (int)(sizeof(cols)/sizeof(cols[0]))) {
        int kind = classify_col(tok);
        cols[ncols].kind = kind;
        strncpy(cols[ncols].name, tok, sizeof(cols[ncols].name) - 1);
        cols[ncols].name[sizeof(cols[ncols].name) - 1] = '\0';
        if (kind == COL_ID)                         id_col = ncols;
        else if (kind == COL_TYPE)                  type_col = ncols;
        else if (kind == COL_X || kind == COL_XU)   x_col = ncols;
        else if (kind == COL_Y || kind == COL_YU)   y_col = ncols;
        else if (kind == COL_Z || kind == COL_ZU)   z_col = ncols;
        else if (kind == COL_XS) { x_col = ncols; scale_mode = 1; }
        else if (kind == COL_YS) { y_col = ncols; scale_mode = 1; }
        else if (kind == COL_ZS) { z_col = ncols; scale_mode = 1; }
        ncols++;
        tok = strtok_r(NULL, " \t", &save);
    }
    if (id_col < 0 || type_col < 0 || x_col < 0 || y_col < 0 || z_col < 0) {
        dc_errf(errbuf, errbuf_sz,
                "LAMMPS dump: ATOMS header must contain id, type, x/y/z");
        fclose(f);
        return DC_ERR_FORMAT;
    }

    a->n = n_atoms;
    a->style = DC_ATOM_STYLE_ATOMIC;
    dc_status_t rc = dc_atoms_reserve(a, n_atoms);
    if (rc != DC_OK) { fclose(f); return rc; }
    a->extra = (char **)calloc(n_atoms, sizeof(char *));
    if (!a->extra) { fclose(f); return DC_ERR_MEMORY; }

    /* Box extents for scaled-coordinate conversion. */
    double lx = a->box.xhi - a->box.xlo;
    double ly = a->box.yhi - a->box.ylo;
    double lz = a->box.zhi - a->box.zlo;

    /* Tokenize each atom line, plucking the values by column index. */
    char *values[64];
    for (size_t i = 0; i < n_atoms; ++i) {
        if (!fgets_line(f, line, sizeof(line))) {
            dc_errf(errbuf, errbuf_sz,
                    "LAMMPS dump: truncated atoms section at line %zu", i + 1);
            fclose(f);
            return DC_ERR_FORMAT;
        }
        /* Split into up to ncols tokens; remaining tail after the last
         * core column is kept as trailing text. */
        int nv = 0;
        char *sp = line;
        while (*sp && nv < ncols) {
            while (*sp == ' ' || *sp == '\t') sp++;
            if (!*sp) break;
            values[nv++] = sp;
            while (*sp && *sp != ' ' && *sp != '\t') sp++;
            if (*sp) { *sp = '\0'; sp++; }
        }
        if (nv < ncols) {
            /* Not all columns present — try to proceed if core columns ok. */
        }

        int64_t id_v = 0;
        int32_t type_v = 0;
        double x_v = 0, y_v = 0, z_v = 0;

        if (id_col < nv)   id_v   = (int64_t)strtoll(values[id_col], NULL, 10);
        if (type_col < nv) type_v = (int32_t)strtol(values[type_col], NULL, 10);
        if (x_col < nv)    x_v    = strtod(values[x_col], NULL);
        if (y_col < nv)    y_v    = strtod(values[y_col], NULL);
        if (z_col < nv)    z_v    = strtod(values[z_col], NULL);

        if (scale_mode) {
            x_v = a->box.xlo + x_v * lx;
            y_v = a->box.ylo + y_v * ly;
            z_v = a->box.zlo + z_v * lz;
        }

        a->id[i]   = id_v;
        a->type[i] = type_v;
        a->x[i] = x_v;
        a->y[i] = y_v;
        a->z[i] = z_v;

        /* Build trailing text: all non-core columns joined by spaces.
         * This preserves extra fields on round-trip to LAMMPS data. */
        char tail[2048];
        size_t tl = 0;
        tail[0] = '\0';
        for (int k = 0; k < nv; ++k) {
            if (k == id_col || k == type_col ||
                k == x_col  || k == y_col  || k == z_col) continue;
            size_t vl = strlen(values[k]);
            if (tl + vl + 2 >= sizeof(tail)) break;
            if (tl > 0) tail[tl++] = ' ';
            memcpy(tail + tl, values[k], vl);
            tl += vl;
            tail[tl] = '\0';
        }
        if (tl > 0) {
            a->extra[i] = strdup(tail);
            if (!a->extra[i]) { fclose(f); return DC_ERR_MEMORY; }
        } else {
            a->extra[i] = NULL;
        }
    }

    fclose(f);
    return DC_OK;
}
