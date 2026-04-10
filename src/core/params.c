/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* params.c - parameter file parsing and validation.
 *
 * Accepts the legacy dcreator .dissparam format plus the new keys
 * introduced by dcreator 2.0.  Unknown keys produce a warning to
 * stderr and are ignored (matches legacy behavior).
 *
 * Recognized keys:
 *   configfile       <path>
 *   outfile          <path>
 *   new_x/new_y/new_z <dx dy dz>
 *   burgersv         <bx by bz>          (magnitude encoded)
 *   linev            <lx ly lz>
 *   glideplane       <nx ny nz>
 *   point            <ox oy oz>
 *   stripe           <width>
 *
 *   displacement_mode  upper_half | lower_half | symmetric    (default: upper_half)
 *   rigid_side         right | left                           (default: right)
 *   input_format       auto | lammps_data | lammps_dump | imd (default: auto)
 *   output_format      lammps_data | imd                      (default: lammps_data)
 *   lammps_atom_style  atomic | charge | molecular | full | bond  (default: atomic)
 *   num_threads        <N>                                    (default: 0 = OMP env)
 */

#include "internal.h"
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

void dc_params_init(dc_params_t *p)
{
    memset(p, 0, sizeof(*p));
    p->input_format  = DC_FORMAT_AUTO;
    p->output_format = DC_FORMAT_LAMMPS_DATA;
    p->lammps_style  = DC_ATOM_STYLE_ATOMIC;
    p->mode          = DC_MODE_UPPER_HALF;
    p->rigid_side    = DC_RIGID_RIGHT;
    p->num_threads   = 0;
}

const char *dc_format_name(dc_format_t f)
{
    switch (f) {
    case DC_FORMAT_AUTO:        return "auto";
    case DC_FORMAT_LAMMPS_DATA: return "lammps_data";
    case DC_FORMAT_LAMMPS_DUMP: return "lammps_dump";
    case DC_FORMAT_IMD:         return "imd";
    case DC_FORMAT_UNKNOWN:     return "unknown";
    }
    return "?";
}

const char *dc_mode_name(dc_mode_t m)
{
    switch (m) {
    case DC_MODE_UPPER_HALF: return "upper_half";
    case DC_MODE_LOWER_HALF: return "lower_half";
    case DC_MODE_SYMMETRIC:  return "symmetric";
    }
    return "?";
}

const char *dc_rigid_side_name(dc_rigid_side_t r)
{
    return r == DC_RIGID_RIGHT ? "right" : "left";
}

const char *dc_atom_style_name(dc_atom_style_t s)
{
    switch (s) {
    case DC_ATOM_STYLE_ATOMIC:    return "atomic";
    case DC_ATOM_STYLE_CHARGE:    return "charge";
    case DC_ATOM_STYLE_MOLECULAR: return "molecular";
    case DC_ATOM_STYLE_FULL:      return "full";
    case DC_ATOM_STYLE_BOND:      return "bond";
    }
    return "?";
}

/* --- token helpers --- */

static int next_real(char **saveptr, double *out)
{
    char *tok = strtok_r(NULL, " \t\n\r", saveptr);
    if (!tok) return 0;
    char *end = NULL;
    double v = strtod(tok, &end);
    if (end == tok) return 0;
    *out = v;
    return 1;
}

static int next_int(char **saveptr, long *out)
{
    char *tok = strtok_r(NULL, " \t\n\r", saveptr);
    if (!tok) return 0;
    char *end = NULL;
    long v = strtol(tok, &end, 10);
    if (end == tok) return 0;
    *out = v;
    return 1;
}

static int next_str(char **saveptr, char *out, size_t out_sz)
{
    char *tok = strtok_r(NULL, " \t\n\r", saveptr);
    if (!tok) return 0;
    strncpy(out, tok, out_sz - 1);
    out[out_sz - 1] = '\0';
    return 1;
}

static int read_vec3(char **saveptr, dc_vec3_t *v)
{
    return next_real(saveptr, &v->x) &&
           next_real(saveptr, &v->y) &&
           next_real(saveptr, &v->z);
}

/* Strip '#' comments in place (anything after first '#' on a line). */
static void strip_comment(char *line)
{
    for (char *c = line; *c; ++c) {
        if (*c == '#') { *c = '\0'; break; }
    }
}

dc_status_t dc_params_load(dc_params_t *p, const char *filename,
                           char *errbuf, size_t errbuf_sz)
{
    FILE *f = fopen(filename, "r");
    if (!f) {
        dc_errf(errbuf, errbuf_sz, "cannot open parameter file '%s': %s",
                filename, strerror(errno));
        return DC_ERR_IO;
    }

    char buf[2048];
    int lineno = 0;
    dc_status_t status = DC_OK;

    while (fgets(buf, sizeof(buf), f)) {
        lineno++;
        strip_comment(buf);
        char *saveptr = NULL;
        char *key = strtok_r(buf, " \t\n\r", &saveptr);
        if (!key) continue;

        int ok = 1;
        const char *why = "";

        if (strcasecmp(key, "configfile") == 0) {
            if (!next_str(&saveptr, p->input_file, sizeof(p->input_file))) {
                ok = 0; why = "expected filename";
            }
        } else if (strcasecmp(key, "outfile") == 0) {
            if (!next_str(&saveptr, p->output_file, sizeof(p->output_file))) {
                ok = 0; why = "expected filename";
            }
        } else if (strcasecmp(key, "new_x") == 0) {
            if (!read_vec3(&saveptr, &p->new_x)) { ok = 0; why = "expected 3 reals"; }
        } else if (strcasecmp(key, "new_y") == 0) {
            if (!read_vec3(&saveptr, &p->new_y)) { ok = 0; why = "expected 3 reals"; }
        } else if (strcasecmp(key, "new_z") == 0) {
            if (!read_vec3(&saveptr, &p->new_z)) { ok = 0; why = "expected 3 reals"; }
        } else if (strcasecmp(key, "burgersv") == 0 ||
                   strcasecmp(key, "burgers")  == 0) {
            if (!read_vec3(&saveptr, &p->burgers)) { ok = 0; why = "expected 3 reals"; }
        } else if (strcasecmp(key, "linev") == 0 ||
                   strcasecmp(key, "line")  == 0) {
            if (!read_vec3(&saveptr, &p->line)) { ok = 0; why = "expected 3 reals"; }
        } else if (strcasecmp(key, "glideplane") == 0) {
            if (!read_vec3(&saveptr, &p->glide_plane_normal)) {
                ok = 0; why = "expected 3 reals";
            }
        } else if (strcasecmp(key, "point") == 0 ||
                   strcasecmp(key, "origin") == 0) {
            if (!read_vec3(&saveptr, &p->origin)) { ok = 0; why = "expected 3 reals"; }
        } else if (strcasecmp(key, "stripe") == 0) {
            if (!next_real(&saveptr, &p->stripe)) { ok = 0; why = "expected real"; }
        }
        /* --- new keys --- */
        else if (strcasecmp(key, "displacement_mode") == 0 ||
                 strcasecmp(key, "mode") == 0) {
            char val[64];
            if (!next_str(&saveptr, val, sizeof(val))) {
                ok = 0; why = "expected mode name";
            } else if (strcasecmp(val, "upper_half") == 0) p->mode = DC_MODE_UPPER_HALF;
              else if (strcasecmp(val, "lower_half") == 0) p->mode = DC_MODE_LOWER_HALF;
              else if (strcasecmp(val, "symmetric")  == 0 ||
                       strcasecmp(val, "legacy")     == 0) p->mode = DC_MODE_SYMMETRIC;
              else { ok = 0; why = "unknown displacement_mode"; }
        } else if (strcasecmp(key, "rigid_side") == 0) {
            char val[64];
            if (!next_str(&saveptr, val, sizeof(val))) {
                ok = 0; why = "expected side";
            } else if (strcasecmp(val, "right") == 0) p->rigid_side = DC_RIGID_RIGHT;
              else if (strcasecmp(val, "left")  == 0) p->rigid_side = DC_RIGID_LEFT;
              else { ok = 0; why = "expected 'right' or 'left'"; }
        } else if (strcasecmp(key, "input_format") == 0) {
            char val[64];
            if (!next_str(&saveptr, val, sizeof(val))) {
                ok = 0; why = "expected format name";
            } else if (strcasecmp(val, "auto")        == 0) p->input_format = DC_FORMAT_AUTO;
              else if (strcasecmp(val, "lammps_data") == 0) p->input_format = DC_FORMAT_LAMMPS_DATA;
              else if (strcasecmp(val, "lammps_dump") == 0) p->input_format = DC_FORMAT_LAMMPS_DUMP;
              else if (strcasecmp(val, "imd")         == 0) p->input_format = DC_FORMAT_IMD;
              else { ok = 0; why = "unknown input_format"; }
        } else if (strcasecmp(key, "output_format") == 0) {
            char val[64];
            if (!next_str(&saveptr, val, sizeof(val))) {
                ok = 0; why = "expected format name";
            } else if (strcasecmp(val, "lammps_data") == 0) p->output_format = DC_FORMAT_LAMMPS_DATA;
              else if (strcasecmp(val, "imd")         == 0) p->output_format = DC_FORMAT_IMD;
              else { ok = 0; why = "output_format must be lammps_data or imd"; }
        } else if (strcasecmp(key, "lammps_atom_style") == 0) {
            char val[64];
            if (!next_str(&saveptr, val, sizeof(val))) {
                ok = 0; why = "expected atom style";
            } else if (strcasecmp(val, "atomic")    == 0) p->lammps_style = DC_ATOM_STYLE_ATOMIC;
              else if (strcasecmp(val, "charge")    == 0) p->lammps_style = DC_ATOM_STYLE_CHARGE;
              else if (strcasecmp(val, "molecular") == 0) p->lammps_style = DC_ATOM_STYLE_MOLECULAR;
              else if (strcasecmp(val, "full")      == 0) p->lammps_style = DC_ATOM_STYLE_FULL;
              else if (strcasecmp(val, "bond")      == 0) p->lammps_style = DC_ATOM_STYLE_BOND;
              else { ok = 0; why = "unknown lammps_atom_style"; }
        } else if (strcasecmp(key, "num_threads") == 0) {
            long v;
            if (!next_int(&saveptr, &v)) { ok = 0; why = "expected integer"; }
            else p->num_threads = (int)v;
        } else {
            fprintf(stderr, "dcreator: warning: unknown key '%s' at line %d ignored\n",
                    key, lineno);
            continue;
        }

        if (!ok) {
            dc_errf(errbuf, errbuf_sz, "parameter file '%s' line %d: %s (%s)",
                    filename, lineno, why, key);
            status = DC_ERR_PARAM;
            break;
        }
    }
    fclose(f);
    return status;
}

dc_status_t dc_params_validate(const dc_params_t *p,
                               char *errbuf, size_t errbuf_sz)
{
    if (p->input_file[0] == '\0') {
        dc_errf(errbuf, errbuf_sz, "configfile/input file not specified");
        return DC_ERR_PARAM;
    }
    if (p->output_file[0] == '\0') {
        dc_errf(errbuf, errbuf_sz, "outfile/output file not specified");
        return DC_ERR_PARAM;
    }
    int nn = (!dc_is_zero(p->new_x)) + (!dc_is_zero(p->new_y)) + (!dc_is_zero(p->new_z));
    if (nn < 2) {
        dc_errf(errbuf, errbuf_sz,
                "need at least two of new_x/new_y/new_z to be non-zero");
        return DC_ERR_PARAM;
    }
    if (dc_is_zero(p->burgers)) {
        dc_errf(errbuf, errbuf_sz, "burgersv is zero");
        return DC_ERR_PARAM;
    }
    if (dc_is_zero(p->line)) {
        dc_errf(errbuf, errbuf_sz, "linev is zero");
        return DC_ERR_PARAM;
    }
    if (dc_is_zero(p->glide_plane_normal)) {
        dc_errf(errbuf, errbuf_sz, "glideplane is zero");
        return DC_ERR_PARAM;
    }
    if (p->stripe < 0.0) {
        dc_errf(errbuf, errbuf_sz, "stripe must be >= 0");
        return DC_ERR_PARAM;
    }
    if (p->output_format != DC_FORMAT_LAMMPS_DATA &&
        p->output_format != DC_FORMAT_IMD) {
        dc_errf(errbuf, errbuf_sz,
                "output_format must be lammps_data or imd");
        return DC_ERR_PARAM;
    }
    return DC_OK;
}
