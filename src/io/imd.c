/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* imd.c - reader and writer for the IMD ASCII checkpoint format.
 *
 * Header:
 *   #F <fmt> <n_id> <n_type> <n_mass> <n_pos> <n_vel> <n_data>
 *   #C <column names...>
 *   #X <ax ay az>     cell vector a
 *   #Y <bx by bz>     cell vector b
 *   #Z <cx cy cz>     cell vector c
 *   ## <free text>    (ignored, preserved on write)
 *   #E                end of header
 *   <atom lines>
 *
 * fmt is 'A' for ASCII. Binary formats are NOT supported in v1.
 *
 * Each atom line: id type [mass] x y z [vx vy vz] [extra columns...]
 * The presence of id/type/mass/vel/extras is driven by the #F fields.
 */

#include "../core/internal.h"
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IMD_LINE_BUF  4096

typedef struct {
    int fmt_char;
    int n_id, n_type, n_mass, n_pos, n_vel, n_data;
    dc_vec3_t cell_a, cell_b, cell_c;
    int has_cell_a, has_cell_b, has_cell_c;
    /* free-form comment lines, captured verbatim for round-trip */
    char **comments;
    size_t n_comments;
    size_t cap_comments;
} imd_header_t;

static void imd_header_init(imd_header_t *h) { memset(h, 0, sizeof(*h)); }
static void imd_header_free(imd_header_t *h)
{
    for (size_t i = 0; i < h->n_comments; ++i) free(h->comments[i]);
    free(h->comments);
    memset(h, 0, sizeof(*h));
}

static int imd_push_comment(imd_header_t *h, const char *line)
{
    if (h->n_comments == h->cap_comments) {
        size_t newcap = h->cap_comments ? h->cap_comments * 2 : 8;
        char **p = (char **)realloc(h->comments, newcap * sizeof(char *));
        if (!p) return 0;
        h->comments = p;
        h->cap_comments = newcap;
    }
    h->comments[h->n_comments] = strdup(line);
    if (!h->comments[h->n_comments]) return 0;
    h->n_comments++;
    return 1;
}

static int parse_F(const char *s, imd_header_t *h, char *errbuf, size_t sz)
{
    /* skip "#F" */
    while (*s && !isspace((unsigned char)*s)) s++;
    while (*s && isspace((unsigned char)*s)) s++;
    if (!*s) { dc_errf(errbuf, sz, "IMD: empty #F line"); return 0; }

    h->fmt_char = *s++;
    if (h->fmt_char != 'A' && h->fmt_char != 'a') {
        dc_errf(errbuf, sz, "IMD: only ASCII ('A') format supported, got '%c'",
                h->fmt_char);
        return 0;
    }
    int n_id=0, n_type=0, n_mass=0, n_pos=0, n_vel=0, n_data=0;
    if (sscanf(s, "%d %d %d %d %d %d",
               &n_id, &n_type, &n_mass, &n_pos, &n_vel, &n_data) < 5) {
        dc_errf(errbuf, sz, "IMD: malformed #F line");
        return 0;
    }
    if (n_pos != 3) {
        dc_errf(errbuf, sz, "IMD: expected n_pos=3, got %d", n_pos);
        return 0;
    }
    if (n_vel != 0 && n_vel != 3) {
        dc_errf(errbuf, sz, "IMD: expected n_vel 0 or 3, got %d", n_vel);
        return 0;
    }
    h->n_id   = n_id;
    h->n_type = n_type;
    h->n_mass = n_mass;
    h->n_pos  = n_pos;
    h->n_vel  = n_vel;
    h->n_data = n_data;
    return 1;
}

static int parse_vec3_after_tag(const char *s, dc_vec3_t *v)
{
    /* skip tag */
    while (*s && !isspace((unsigned char)*s)) s++;
    return sscanf(s, "%lf %lf %lf", &v->x, &v->y, &v->z) == 3;
}

static void imd_cell_to_box(const imd_header_t *h, dc_box_t *box)
{
    memset(box, 0, sizeof(*box));
    /* Orthogonal: origin at 0, extents = diagonal components. */
    box->xlo = 0.0; box->xhi = h->cell_a.x;
    box->ylo = 0.0; box->yhi = h->cell_b.y;
    box->zlo = 0.0; box->zhi = h->cell_c.z;
    /* Triclinic tilt factors (LAMMPS convention lower-triangular). */
    box->xy = h->cell_b.x;
    box->xz = h->cell_c.x;
    box->yz = h->cell_c.y;
    if (box->xy != 0.0 || box->xz != 0.0 || box->yz != 0.0) {
        box->triclinic = 1;
    }
    /* Assume off-diagonal parts of cell_a, cell_b.z, cell_c.y handled;
     * cell_a.y, cell_a.z, cell_b.z should be zero for a standard IMD box. */
}

static dc_status_t parse_imd_header(FILE *f, imd_header_t *h, int *have_F,
                                    char *errbuf, size_t errbuf_sz)
{
    char line[IMD_LINE_BUF];
    *have_F = 0;
    while (fgets(line, sizeof(line), f)) {
        /* must start with '#' */
        if (line[0] != '#') {
            dc_errf(errbuf, errbuf_sz,
                    "IMD: header ended before '#E' marker");
            return DC_ERR_FORMAT;
        }
        if (line[1] == 'E' || line[1] == 'e') return DC_OK;

        if (line[1] == 'F' || line[1] == 'f') {
            if (!parse_F(line, h, errbuf, errbuf_sz)) return DC_ERR_FORMAT;
            *have_F = 1;
        } else if (line[1] == 'X' || line[1] == 'x') {
            if (!parse_vec3_after_tag(line, &h->cell_a)) {
                dc_errf(errbuf, errbuf_sz, "IMD: malformed #X line");
                return DC_ERR_FORMAT;
            }
            h->has_cell_a = 1;
        } else if (line[1] == 'Y' || line[1] == 'y') {
            if (!parse_vec3_after_tag(line, &h->cell_b)) {
                dc_errf(errbuf, errbuf_sz, "IMD: malformed #Y line");
                return DC_ERR_FORMAT;
            }
            h->has_cell_b = 1;
        } else if (line[1] == 'Z' || line[1] == 'z') {
            if (!parse_vec3_after_tag(line, &h->cell_c)) {
                dc_errf(errbuf, errbuf_sz, "IMD: malformed #Z line");
                return DC_ERR_FORMAT;
            }
            h->has_cell_c = 1;
        } else if (line[1] == '#') {
            /* user comment line; store verbatim without trailing newline */
            size_t n = strlen(line);
            if (n > 0 && line[n-1] == '\n') line[n-1] = '\0';
            if (!imd_push_comment(h, line)) {
                dc_errf(errbuf, errbuf_sz, "IMD: out of memory");
                return DC_ERR_MEMORY;
            }
        } else {
            /* #C and other header lines: ignore silently */
        }
    }
    dc_errf(errbuf, errbuf_sz, "IMD: reached EOF before '#E' marker");
    return DC_ERR_FORMAT;
}

dc_status_t dc_read_imd(const char *filename, dc_atoms_t *a,
                        char *errbuf, size_t errbuf_sz)
{
    FILE *f = fopen(filename, "r");
    if (!f) {
        dc_errf(errbuf, errbuf_sz, "cannot open IMD file '%s': %s",
                filename, strerror(errno));
        return DC_ERR_IO;
    }
    imd_header_t h;
    imd_header_init(&h);
    int have_F = 0;
    dc_status_t st = parse_imd_header(f, &h, &have_F, errbuf, errbuf_sz);
    if (st != DC_OK) { imd_header_free(&h); fclose(f); return st; }
    if (!have_F) {
        dc_errf(errbuf, errbuf_sz, "IMD: missing #F header line");
        imd_header_free(&h); fclose(f); return DC_ERR_FORMAT;
    }

    /* Pre-count atom lines so we can allocate exactly once.  For huge
     * files we pay one read pass; alternative is amortized-doubling
     * realloc which is also fine. Keep it simple: count. */
    long atoms_start = ftell(f);
    size_t n_atoms = 0;
    char line[IMD_LINE_BUF];
    while (fgets(line, sizeof(line), f)) {
        const char *s = line;
        while (*s && isspace((unsigned char)*s)) s++;
        if (*s == '\0' || *s == '#') continue;
        n_atoms++;
    }
    fseek(f, atoms_start, SEEK_SET);

    dc_atoms_init(a);
    a->source_format = DC_FORMAT_IMD;
    a->n = n_atoms;
    st = dc_atoms_reserve(a, n_atoms);
    if (st != DC_OK) {
        imd_header_free(&h); fclose(f); return st;
    }
    if (h.n_mass == 1) {
        a->mass = (double *)calloc(n_atoms, sizeof(double));
        if (!a->mass) { imd_header_free(&h); fclose(f); return DC_ERR_MEMORY; }
    }
    if (h.n_vel == 3) {
        a->vx = (double *)calloc(n_atoms, sizeof(double));
        a->vy = (double *)calloc(n_atoms, sizeof(double));
        a->vz = (double *)calloc(n_atoms, sizeof(double));
        if (!a->vx || !a->vy || !a->vz) {
            imd_header_free(&h); fclose(f); return DC_ERR_MEMORY;
        }
    }

    /* Read atoms. */
    size_t idx = 0;
    int lineno = 0;
    while (fgets(line, sizeof(line), f)) {
        lineno++;
        const char *s = line;
        while (*s && isspace((unsigned char)*s)) s++;
        if (*s == '\0' || *s == '#') continue;

        char *end;
        long id = 0;
        long type = 0;
        double mass = 0.0;
        double x, y, z;
        double vx = 0.0, vy = 0.0, vz = 0.0;

        char *p = (char *)s;
        if (h.n_id) {
            id = strtol(p, &end, 10);
            if (end == p) goto parse_err;
            p = end;
        }
        if (h.n_type) {
            type = strtol(p, &end, 10);
            if (end == p) goto parse_err;
            p = end;
        }
        if (h.n_mass) {
            mass = strtod(p, &end);
            if (end == p) goto parse_err;
            p = end;
        }
        x = strtod(p, &end); if (end == p) goto parse_err; p = end;
        y = strtod(p, &end); if (end == p) goto parse_err; p = end;
        z = strtod(p, &end); if (end == p) goto parse_err; p = end;
        if (h.n_vel == 3) {
            vx = strtod(p, &end); if (end == p) goto parse_err; p = end;
            vy = strtod(p, &end); if (end == p) goto parse_err; p = end;
            vz = strtod(p, &end); if (end == p) goto parse_err; p = end;
        }

        a->id[idx]   = (int64_t)id;
        a->type[idx] = (int32_t)type;
        a->x[idx] = x;
        a->y[idx] = y;
        a->z[idx] = z;
        if (a->mass) a->mass[idx] = mass;
        if (a->vx)   a->vx[idx] = vx;
        if (a->vy)   a->vy[idx] = vy;
        if (a->vz)   a->vz[idx] = vz;
        idx++;
        continue;

parse_err:
        dc_errf(errbuf, errbuf_sz,
                "IMD: parse error in atom line %d of '%s'", lineno, filename);
        fclose(f);
        dc_atoms_free(a);
        imd_header_free(&h);
        return DC_ERR_FORMAT;
    }
    fclose(f);

    /* Build the box record. */
    imd_cell_to_box(&h, &a->box);

    /* Serialize the header blob (comments + #F + #C + #X/#Y/#Z) so the
     * writer can round-trip it later. We stash it as a printable string. */
    {
        size_t cap = 1024;
        char *blob = (char *)malloc(cap);
        if (!blob) { imd_header_free(&h); return DC_ERR_MEMORY; }
        size_t len = 0;
        #define APPEND(fmt_, ...) do { \
            int needed = snprintf(NULL, 0, fmt_, ##__VA_ARGS__) + 1; \
            while (len + (size_t)needed >= cap) { \
                cap *= 2; \
                char *nb = (char *)realloc(blob, cap); \
                if (!nb) { free(blob); imd_header_free(&h); return DC_ERR_MEMORY; } \
                blob = nb; \
            } \
            len += (size_t)snprintf(blob + len, cap - len, fmt_, ##__VA_ARGS__); \
        } while (0)

        APPEND("#F %c %d %d %d %d %d %d\n",
               h.fmt_char, h.n_id, h.n_type, h.n_mass, h.n_pos, h.n_vel, h.n_data);
        /* Column names (reconstructed) */
        APPEND("#C");
        if (h.n_id)   APPEND(" number");
        if (h.n_type) APPEND(" type");
        if (h.n_mass) APPEND(" mass");
        APPEND(" x y z");
        if (h.n_vel)  APPEND(" vx vy vz");
        APPEND("\n");
        APPEND("#X %.16g %.16g %.16g\n", h.cell_a.x, h.cell_a.y, h.cell_a.z);
        APPEND("#Y %.16g %.16g %.16g\n", h.cell_b.x, h.cell_b.y, h.cell_b.z);
        APPEND("#Z %.16g %.16g %.16g\n", h.cell_c.x, h.cell_c.y, h.cell_c.z);
        for (size_t i = 0; i < h.n_comments; ++i) {
            APPEND("%s\n", h.comments[i]);
        }
        APPEND("#E\n");
        #undef APPEND

        a->header_blob = blob;
        a->header_blob_len = len;
    }

    imd_header_free(&h);
    return DC_OK;
}

dc_status_t dc_write_imd(const char *filename, const dc_atoms_t *a,
                         char *errbuf, size_t errbuf_sz)
{
    FILE *f = fopen(filename, "w");
    if (!f) {
        dc_errf(errbuf, errbuf_sz, "cannot open '%s' for writing: %s",
                filename, strerror(errno));
        return DC_ERR_IO;
    }

    /* Reuse captured header only when the source format matches; otherwise
     * the blob is for a different format and must not be emitted. */
    if (a->header_blob && a->header_blob_len > 0 &&
        a->source_format == DC_FORMAT_IMD) {
        fwrite(a->header_blob, 1, a->header_blob_len, f);
    } else {
        /* Synthesize a header from box info and column presence. */
        int n_mass = (a->mass != NULL) ? 1 : 0;
        int n_vel  = (a->vx   != NULL) ? 3 : 0;
        fprintf(f, "#F A 1 1 %d 3 %d 0\n", n_mass, n_vel);
        fprintf(f, "#C number type");
        if (n_mass) fprintf(f, " mass");
        fprintf(f, " x y z");
        if (n_vel) fprintf(f, " vx vy vz");
        fprintf(f, "\n");
        /* Convert LAMMPS-style box back to IMD cell vectors. */
        double ax = a->box.xhi - a->box.xlo;
        double by = a->box.yhi - a->box.ylo;
        double cz = a->box.zhi - a->box.zlo;
        fprintf(f, "#X %.16g 0 0\n", ax);
        fprintf(f, "#Y %.16g %.16g 0\n", a->box.xy, by);
        fprintf(f, "#Z %.16g %.16g %.16g\n", a->box.xz, a->box.yz, cz);
        fprintf(f, "## written by dcreator 2.0\n");
        fprintf(f, "#E\n");
    }

    for (size_t i = 0; i < a->n; ++i) {
        fprintf(f, "%lld %d",
                (long long)a->id[i], (int)a->type[i]);
        if (a->mass) fprintf(f, " %.16g", a->mass[i]);
        fprintf(f, " %.16f %.16f %.16f", a->x[i], a->y[i], a->z[i]);
        if (a->vx)   fprintf(f, " %.16g %.16g %.16g", a->vx[i], a->vy[i], a->vz[i]);
        fputc('\n', f);
    }
    if (fclose(f) != 0) {
        dc_errf(errbuf, errbuf_sz, "error closing '%s': %s",
                filename, strerror(errno));
        return DC_ERR_IO;
    }
    return DC_OK;
}
