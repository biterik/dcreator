/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* lammps_data.c - reader and writer for LAMMPS data files.
 *
 * The strategy is "transcribe with surgery": we pass every header
 * section and every non-Atoms section through verbatim, only touching
 * the Atoms section. That keeps Masses, Velocities, Bonds, Pair Coeffs
 * and any other LAMMPS-specific content intact.
 *
 * Within the Atoms section we parse the style-dependent columns and
 * store any trailing text (image flags, arbitrary extras) per-atom
 * for untouched round-trip.
 *
 * Supported atom_style on read: atomic, charge, molecular, full, bond.
 * Image flags and extra trailing columns are preserved verbatim.
 */

#include "../core/internal.h"
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

/* Return pointer to first non-space char, NULL if line is empty/commented. */
static const char *skip_ws(const char *s)
{
    while (*s && isspace((unsigned char)*s)) s++;
    return *s ? s : NULL;
}

/* Does `line` (after whitespace) start with `keyword` followed by end-of-line
 * or whitespace+anything? Case-sensitive — LAMMPS section names are fixed. */
static int line_is_section(const char *line, const char *keyword)
{
    line = skip_ws(line);
    if (!line) return 0;
    size_t klen = strlen(keyword);
    if (strncmp(line, keyword, klen) != 0) return 0;
    char c = line[klen];
    return c == '\0' || c == '\n' || c == ' ' || c == '\t' || c == '#';
}

/* LAMMPS section names that we recognize. Any line matching one of these
 * (at start after whitespace) is a section header. */
static const char *const kSections[] = {
    "Atoms", "Masses", "Velocities", "Bonds", "Angles", "Dihedrals",
    "Impropers", "Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs",
    "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs", "Ellipsoids",
    "Lines", "Triangles", "Bodies", "BondBond Coeffs", "BondAngle Coeffs",
    "MiddleBondTorsion Coeffs", "EndBondTorsion Coeffs",
    "AngleTorsion Coeffs", "AngleAngleTorsion Coeffs", "BondBond13 Coeffs",
    "AngleAngle Coeffs",
    NULL
};

static int line_is_any_section(const char *line)
{
    for (int i = 0; kSections[i]; ++i)
        if (line_is_section(line, kSections[i])) return 1;
    return 0;
}

/* Detect atom_style from an "Atoms # <style>" header line.
 * Returns the detected style or -1 if none found. */
static int detect_style_comment(const char *line)
{
    const char *h = strchr(line, '#');
    if (!h) return -1;
    h++;
    while (*h && isspace((unsigned char)*h)) h++;
    if (strncasecmp(h, "atomic",    6) == 0) return DC_ATOM_STYLE_ATOMIC;
    if (strncasecmp(h, "charge",    6) == 0) return DC_ATOM_STYLE_CHARGE;
    if (strncasecmp(h, "molecular", 9) == 0) return DC_ATOM_STYLE_MOLECULAR;
    if (strncasecmp(h, "full",      4) == 0) return DC_ATOM_STYLE_FULL;
    if (strncasecmp(h, "bond",      4) == 0) return DC_ATOM_STYLE_BOND;
    return -1;
}

static void parse_header_counts(const char *buf, size_t len, size_t *n_atoms,
                                dc_box_t *box, int *have_n)
{
    /* Walk line by line over the header portion. */
    const char *p = buf;
    const char *end = buf + len;
    char line[1024];
    *have_n = 0;
    memset(box, 0, sizeof(*box));

    while (p < end) {
        size_t m = 0;
        while (p + m < end && p[m] != '\n' && m + 1 < sizeof(line)) m++;
        memcpy(line, p, m);
        line[m] = '\0';
        p += m;
        if (p < end && *p == '\n') p++;

        const char *s = skip_ws(line);
        if (!s || *s == '#') continue;

        /* "N atoms" */
        {
            long v;
            char rest[64];
            if (sscanf(s, "%ld %63s", &v, rest) == 2 &&
                strcmp(rest, "atoms") == 0) {
                *n_atoms = (size_t)v;
                *have_n = 1;
                continue;
            }
        }
        /* box bounds */
        {
            double a, b;
            char t1[16], t2[16];
            if (sscanf(s, "%lf %lf %15s %15s", &a, &b, t1, t2) == 4) {
                if (strcmp(t1, "xlo") == 0 && strcmp(t2, "xhi") == 0) {
                    box->xlo = a; box->xhi = b; continue;
                }
                if (strcmp(t1, "ylo") == 0 && strcmp(t2, "yhi") == 0) {
                    box->ylo = a; box->yhi = b; continue;
                }
                if (strcmp(t1, "zlo") == 0 && strcmp(t2, "zhi") == 0) {
                    box->zlo = a; box->zhi = b; continue;
                }
            }
        }
        /* tilt factors "xy xz yz" */
        {
            double xy, xz, yz;
            char t1[8], t2[8], t3[8];
            if (sscanf(s, "%lf %lf %lf %7s %7s %7s",
                       &xy, &xz, &yz, t1, t2, t3) == 6 &&
                strcmp(t1, "xy") == 0 && strcmp(t2, "xz") == 0 &&
                strcmp(t3, "yz") == 0) {
                box->xy = xy; box->xz = xz; box->yz = yz;
                box->triclinic = 1;
                continue;
            }
        }
    }
}

/* Parse one atom line according to style. On success returns 1 and
 * advances *pp past the core columns; the remaining trailing text
 * (image flags etc.) is written (copied) to *out_trailing if non-NULL. */
static int parse_atom_line(const char *line, dc_atom_style_t style,
                           int64_t *id, int32_t *type, int64_t *mol,
                           double *q, double *x, double *y, double *z,
                           const char **out_trailing)
{
    char *p = (char *)line;
    char *end;

    *id = (int64_t)strtoll(p, &end, 10);
    if (end == p) return 0;
    p = end;

    if (style == DC_ATOM_STYLE_MOLECULAR || style == DC_ATOM_STYLE_FULL ||
        style == DC_ATOM_STYLE_BOND) {
        *mol = (int64_t)strtoll(p, &end, 10);
        if (end == p) return 0;
        p = end;
    } else {
        *mol = 0;
    }

    *type = (int32_t)strtol(p, &end, 10);
    if (end == p) return 0;
    p = end;

    if (style == DC_ATOM_STYLE_CHARGE || style == DC_ATOM_STYLE_FULL) {
        *q = strtod(p, &end);
        if (end == p) return 0;
        p = end;
    } else {
        *q = 0.0;
    }

    *x = strtod(p, &end); if (end == p) return 0; p = end;
    *y = strtod(p, &end); if (end == p) return 0; p = end;
    *z = strtod(p, &end); if (end == p) return 0; p = end;

    if (out_trailing) {
        while (*p == ' ' || *p == '\t') p++;
        *out_trailing = (*p && *p != '\n') ? p : NULL;
    }
    return 1;
}

dc_status_t dc_read_lammps_data(const char *filename, dc_atom_style_t style,
                                dc_atoms_t *a,
                                char *errbuf, size_t errbuf_sz)
{
    char *buf = NULL;
    size_t len = 0;
    int rc = dc_slurp(filename, &buf, &len, errbuf, errbuf_sz);
    if (rc != DC_OK) return rc;

    dc_atoms_init(a);
    a->source_format = DC_FORMAT_LAMMPS_DATA;
    a->style = style;

    /* First pass: walk lines, find the "Atoms" section. Everything
     * before it becomes header_blob. Also parse header counts/box. */
    const char *p = buf;
    const char *end = buf + len;
    const char *atoms_hdr_begin = NULL; /* start of the line "Atoms" */
    const char *atoms_data_begin = NULL;/* start of the first atom line */
    int style_from_comment = -1;

    while (p < end) {
        const char *lp = p;
        while (p < end && *p != '\n') p++;
        size_t llen = (size_t)(p - lp);
        if (p < end) p++;  /* consume newline */

        /* skip leading whitespace without copying */
        const char *s = lp;
        while (s < lp + llen && isspace((unsigned char)*s)) s++;
        if (s == lp + llen) continue;
        if (*s == '#') continue;

        /* Section header? */
        if ((size_t)(lp + llen - s) >= 5 && strncmp(s, "Atoms", 5) == 0 &&
            (s[5] == '\0' || s[5] == '\n' || s[5] == ' ' || s[5] == '\t' || s[5] == '#')) {
            atoms_hdr_begin = lp;
            /* detect comment */
            char linebuf[256];
            size_t ll = llen < sizeof(linebuf) - 1 ? llen : sizeof(linebuf) - 1;
            memcpy(linebuf, lp, ll); linebuf[ll] = '\0';
            style_from_comment = detect_style_comment(linebuf);

            /* Skip blank lines after "Atoms" header. */
            while (p < end) {
                const char *lp2 = p;
                while (p < end && *p != '\n') p++;
                size_t ll2 = (size_t)(p - lp2);
                if (p < end) p++;
                const char *s2 = lp2;
                while (s2 < lp2 + ll2 && isspace((unsigned char)*s2)) s2++;
                if (s2 == lp2 + ll2) continue; /* blank */
                atoms_data_begin = lp2;
                /* rewind p to start of this line so the count loop starts here */
                p = lp2;
                break;
            }
            break;
        }
    }

    if (!atoms_hdr_begin || !atoms_data_begin) {
        dc_errf(errbuf, errbuf_sz, "LAMMPS data: no 'Atoms' section found");
        free(buf);
        return DC_ERR_FORMAT;
    }

    /* Capture header blob verbatim. */
    size_t hdr_len = (size_t)(atoms_hdr_begin - buf);
    a->header_blob = (char *)malloc(hdr_len + 1);
    if (!a->header_blob) { free(buf); return DC_ERR_MEMORY; }
    memcpy(a->header_blob, buf, hdr_len);
    a->header_blob[hdr_len] = '\0';
    a->header_blob_len = hdr_len;

    /* Parse counts and box from the header blob. */
    int have_n = 0;
    size_t n_atoms = 0;
    parse_header_counts(a->header_blob, a->header_blob_len, &n_atoms, &a->box, &have_n);
    if (!have_n) {
        dc_errf(errbuf, errbuf_sz, "LAMMPS data: could not parse atom count");
        free(buf);
        dc_atoms_free(a);
        return DC_ERR_FORMAT;
    }

    /* If the Atoms line carried a style comment, honor it. */
    if (style_from_comment >= 0) {
        a->style = (dc_atom_style_t)style_from_comment;
    }

    /* Allocate atom SoA. Optional columns allocated on demand. */
    a->n = n_atoms;
    rc = dc_atoms_reserve(a, n_atoms);
    if (rc != DC_OK) { free(buf); return rc; }
    if (a->style == DC_ATOM_STYLE_CHARGE || a->style == DC_ATOM_STYLE_FULL) {
        a->q = (double *)calloc(n_atoms, sizeof(double));
        if (!a->q) { free(buf); return DC_ERR_MEMORY; }
    }
    if (a->style == DC_ATOM_STYLE_MOLECULAR || a->style == DC_ATOM_STYLE_FULL ||
        a->style == DC_ATOM_STYLE_BOND) {
        a->molecule = (int64_t *)calloc(n_atoms, sizeof(int64_t));
        if (!a->molecule) { free(buf); return DC_ERR_MEMORY; }
    }
    a->extra = (char **)calloc(n_atoms, sizeof(char *));
    if (!a->extra) { free(buf); return DC_ERR_MEMORY; }

    /* Parse atom lines. After N lines we capture the trailer blob. */
    size_t idx = 0;
    char linebuf[4096];
    while (idx < n_atoms && p < end) {
        const char *lp = p;
        while (p < end && *p != '\n') p++;
        size_t llen = (size_t)(p - lp);
        const char *nl = p;
        if (p < end) p++;

        const char *s = lp;
        while (s < lp + llen && isspace((unsigned char)*s)) s++;
        if (s == lp + llen) continue;
        if (*s == '#') continue;

        size_t ll = llen < sizeof(linebuf) - 1 ? llen : sizeof(linebuf) - 1;
        memcpy(linebuf, lp, ll); linebuf[ll] = '\0';

        int64_t id = 0, mol = 0;
        int32_t type = 0;
        double q = 0, x, y, z;
        const char *trailing = NULL;
        if (!parse_atom_line(linebuf, a->style, &id, &type, &mol, &q,
                             &x, &y, &z, &trailing)) {
            dc_errf(errbuf, errbuf_sz,
                    "LAMMPS data: parse error in atom line %zu", idx + 1);
            free(buf);
            dc_atoms_free(a);
            return DC_ERR_FORMAT;
        }
        a->id[idx]   = id;
        a->type[idx] = type;
        a->x[idx] = x;
        a->y[idx] = y;
        a->z[idx] = z;
        if (a->q)        a->q[idx]        = q;
        if (a->molecule) a->molecule[idx] = mol;
        if (trailing && *trailing) {
            size_t t_len = strlen(trailing);
            /* strip trailing CR/LF */
            while (t_len > 0 && (trailing[t_len-1] == '\r' ||
                                 trailing[t_len-1] == '\n')) t_len--;
            char *t = (char *)malloc(t_len + 1);
            if (!t) { free(buf); dc_atoms_free(a); return DC_ERR_MEMORY; }
            memcpy(t, trailing, t_len);
            t[t_len] = '\0';
            a->extra[idx] = t;
        } else {
            a->extra[idx] = NULL;
        }
        idx++;

        (void)nl;
    }
    if (idx != n_atoms) {
        dc_errf(errbuf, errbuf_sz,
                "LAMMPS data: expected %zu atoms, got %zu", n_atoms, idx);
        free(buf);
        dc_atoms_free(a);
        return DC_ERR_FORMAT;
    }

    /* Remaining content -> trailer blob. */
    size_t tail_len = (size_t)(end - p);
    a->trailer_blob = (char *)malloc(tail_len + 1);
    if (!a->trailer_blob) { free(buf); return DC_ERR_MEMORY; }
    memcpy(a->trailer_blob, p, tail_len);
    a->trailer_blob[tail_len] = '\0';
    a->trailer_blob_len = tail_len;

    free(buf);
    return DC_OK;
}

/* -------- writer -------- */

static void write_atom_line(FILE *f, const dc_atoms_t *a, size_t i)
{
    switch (a->style) {
    case DC_ATOM_STYLE_ATOMIC:
        fprintf(f, "%lld %d %.16g %.16g %.16g",
                (long long)a->id[i], (int)a->type[i],
                a->x[i], a->y[i], a->z[i]);
        break;
    case DC_ATOM_STYLE_CHARGE:
        fprintf(f, "%lld %d %.16g %.16g %.16g %.16g",
                (long long)a->id[i], (int)a->type[i],
                a->q ? a->q[i] : 0.0,
                a->x[i], a->y[i], a->z[i]);
        break;
    case DC_ATOM_STYLE_MOLECULAR:
    case DC_ATOM_STYLE_BOND:
        fprintf(f, "%lld %lld %d %.16g %.16g %.16g",
                (long long)a->id[i],
                (long long)(a->molecule ? a->molecule[i] : 0),
                (int)a->type[i],
                a->x[i], a->y[i], a->z[i]);
        break;
    case DC_ATOM_STYLE_FULL:
        fprintf(f, "%lld %lld %d %.16g %.16g %.16g %.16g",
                (long long)a->id[i],
                (long long)(a->molecule ? a->molecule[i] : 0),
                (int)a->type[i],
                a->q ? a->q[i] : 0.0,
                a->x[i], a->y[i], a->z[i]);
        break;
    }
    if (a->extra && a->extra[i]) {
        fputc(' ', f);
        fputs(a->extra[i], f);
    }
    fputc('\n', f);
}

dc_status_t dc_write_lammps_data(const char *filename, const dc_atoms_t *a,
                                 char *errbuf, size_t errbuf_sz)
{
    FILE *f = fopen(filename, "w");
    if (!f) {
        dc_errf(errbuf, errbuf_sz, "cannot open '%s' for writing: %s",
                filename, strerror(errno));
        return DC_ERR_IO;
    }

    /* Only reuse the captured header blob when the source was also a
     * LAMMPS data file; otherwise the blob is in a foreign format. */
    if (a->header_blob && a->header_blob_len > 0 &&
        a->source_format == DC_FORMAT_LAMMPS_DATA) {
        fwrite(a->header_blob, 1, a->header_blob_len, f);
    } else {
        /* Synthesize a minimal header (IMD -> LAMMPS path). */
        fprintf(f, "LAMMPS data file written by dcreator 2.0\n\n");
        fprintf(f, "%zu atoms\n", a->n);
        /* LAMMPS requires atom types >= 1, but IMD files often use
         * 0-based types. Compute the minimum and shift everything up
         * so the smallest stored type becomes 1. */
        int32_t mint =  2147483647;
        int32_t maxt = -2147483647;
        for (size_t i = 0; i < a->n; ++i) {
            if (a->type[i] < mint) mint = a->type[i];
            if (a->type[i] > maxt) maxt = a->type[i];
        }
        int shift = (mint <= 0) ? (1 - mint) : 0;
        int ntypes = (int)(maxt + shift);
        if (ntypes < 1) ntypes = 1;
        fprintf(f, "%d atom types\n\n", ntypes);
        fprintf(f, "%.16g %.16g xlo xhi\n", a->box.xlo, a->box.xhi);
        fprintf(f, "%.16g %.16g ylo yhi\n", a->box.ylo, a->box.yhi);
        fprintf(f, "%.16g %.16g zlo zhi\n", a->box.zlo, a->box.zhi);
        if (a->box.triclinic) {
            fprintf(f, "%.16g %.16g %.16g xy xz yz\n",
                    a->box.xy, a->box.xz, a->box.yz);
        }
        fprintf(f, "\n");
        /* Masses: synthesize 1.0 if we have no per-atom mass information. */
        if (a->mass) {
            double *by_type = (double *)calloc((size_t)ntypes + 1, sizeof(double));
            if (by_type) {
                for (size_t i = 0; i < a->n; ++i) {
                    int t = (int)a->type[i] + shift;
                    if (t >= 1 && t <= ntypes && by_type[t] == 0.0)
                        by_type[t] = a->mass[i];
                }
                fprintf(f, "Masses\n\n");
                for (int t = 1; t <= ntypes; ++t) {
                    fprintf(f, "%d %.16g\n", t,
                            by_type[t] != 0.0 ? by_type[t] : 1.0);
                }
                fprintf(f, "\n");
                free(by_type);
            }
        }
        fprintf(f, "Atoms # %s\n\n", dc_atom_style_name(a->style));

        /* Write atoms with optional type shift. We briefly mutate the
         * (mutable) type array rather than carrying shift into the
         * writer. This function owns a mutable alias to `a` for this
         * reason; see the cast at the top of the else branch. */
        dc_atoms_t *aw = (dc_atoms_t *)a;
        if (shift) {
            for (size_t i = 0; i < aw->n; ++i) aw->type[i] += shift;
        }
        for (size_t i = 0; i < a->n; ++i) write_atom_line(f, a, i);
        if (shift) {
            for (size_t i = 0; i < aw->n; ++i) aw->type[i] -= shift;
        }
        if (a->vx) {
            fprintf(f, "\nVelocities\n\n");
            for (size_t i = 0; i < a->n; ++i) {
                fprintf(f, "%lld %.16g %.16g %.16g\n",
                        (long long)a->id[i], a->vx[i], a->vy[i], a->vz[i]);
            }
        }
        if (fclose(f) != 0) {
            dc_errf(errbuf, errbuf_sz, "error closing '%s': %s",
                    filename, strerror(errno));
            return DC_ERR_IO;
        }
        return DC_OK;
    }

    /* Header blob path: emit Atoms section then the trailer blob. */
    fprintf(f, "Atoms # %s\n\n", dc_atom_style_name(a->style));
    for (size_t i = 0; i < a->n; ++i) write_atom_line(f, a, i);

    if (a->trailer_blob && a->trailer_blob_len > 0) {
        /* Ensure a blank line between Atoms and trailer if trailer
         * doesn't already start with a blank/newline. */
        if (a->trailer_blob[0] != '\n') fputc('\n', f);
        fwrite(a->trailer_blob, 1, a->trailer_blob_len, f);
    }

    if (fclose(f) != 0) {
        dc_errf(errbuf, errbuf_sz, "error closing '%s': %s",
                filename, strerror(errno));
        return DC_ERR_IO;
    }
    return DC_OK;
}

/* line_is_any_section is reserved for future use (multi-section surgery). */
__attribute__((unused))
static int dc__keep_line_is_any_section(const char *line) {
    return line_is_any_section(line);
}
