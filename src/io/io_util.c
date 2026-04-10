/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* io_util.c - small I/O helpers shared by readers */

#include "../core/internal.h"
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

const char *dc_strerror(dc_status_t s)
{
    switch (s) {
    case DC_OK:              return "ok";
    case DC_ERR_IO:          return "I/O error";
    case DC_ERR_FORMAT:      return "format error";
    case DC_ERR_PARAM:       return "parameter error";
    case DC_ERR_MEMORY:      return "out of memory";
    case DC_ERR_GEOMETRY:    return "geometry error";
    case DC_ERR_UNSUPPORTED: return "unsupported";
    }
    return "unknown";
}

int dc_slurp(const char *filename, char **out_buf, size_t *out_len,
             char *errbuf, size_t errbuf_sz)
{
    *out_buf = NULL;
    *out_len = 0;

    FILE *f = fopen(filename, "rb");
    if (!f) {
        dc_errf(errbuf, errbuf_sz,
                "cannot open '%s': %s", filename, strerror(errno));
        return DC_ERR_IO;
    }
    if (fseek(f, 0, SEEK_END) != 0) {
        dc_errf(errbuf, errbuf_sz, "fseek failed on '%s'", filename);
        fclose(f);
        return DC_ERR_IO;
    }
    long sz = ftell(f);
    if (sz < 0) {
        dc_errf(errbuf, errbuf_sz, "ftell failed on '%s'", filename);
        fclose(f);
        return DC_ERR_IO;
    }
    rewind(f);

    char *buf = (char *)malloc((size_t)sz + 1);
    if (!buf) {
        dc_errf(errbuf, errbuf_sz, "out of memory slurping '%s'", filename);
        fclose(f);
        return DC_ERR_MEMORY;
    }
    size_t nr = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    if (nr != (size_t)sz) {
        dc_errf(errbuf, errbuf_sz,
                "short read on '%s' (%zu of %ld)", filename, nr, sz);
        free(buf);
        return DC_ERR_IO;
    }
    buf[sz] = '\0';

    *out_buf = buf;
    *out_len = (size_t)sz;
    return DC_OK;
}

char *dc_strstrip(char *s)
{
    if (!s) return s;
    while (*s && isspace((unsigned char)*s)) s++;
    if (*s == '\0') return s;
    char *end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) *end-- = '\0';
    return s;
}

int dc_strncasecmp_(const char *a, const char *b, size_t n)
{
    return strncasecmp(a, b, n);
}

dc_status_t dc_read(const char *filename, dc_format_t fmt,
                    dc_atom_style_t style, dc_atoms_t *a,
                    char *errbuf, size_t errbuf_sz)
{
    if (fmt == DC_FORMAT_AUTO) {
        fmt = dc_detect_format(filename);
        if (fmt == DC_FORMAT_UNKNOWN) {
            dc_errf(errbuf, errbuf_sz,
                    "cannot auto-detect format of '%s'", filename);
            return DC_ERR_FORMAT;
        }
    }
    switch (fmt) {
    case DC_FORMAT_LAMMPS_DATA:
        return dc_read_lammps_data(filename, style, a, errbuf, errbuf_sz);
    case DC_FORMAT_LAMMPS_DUMP:
        return dc_read_lammps_dump(filename, a, errbuf, errbuf_sz);
    case DC_FORMAT_IMD:
        return dc_read_imd(filename, a, errbuf, errbuf_sz);
    default:
        dc_errf(errbuf, errbuf_sz, "unsupported input format");
        return DC_ERR_UNSUPPORTED;
    }
}

dc_status_t dc_write(const char *filename, dc_format_t fmt,
                     const dc_atoms_t *a,
                     char *errbuf, size_t errbuf_sz)
{
    switch (fmt) {
    case DC_FORMAT_LAMMPS_DATA:
        return dc_write_lammps_data(filename, a, errbuf, errbuf_sz);
    case DC_FORMAT_IMD:
        return dc_write_imd(filename, a, errbuf, errbuf_sz);
    default:
        dc_errf(errbuf, errbuf_sz,
                "unsupported output format (%s)", dc_format_name(fmt));
        return DC_ERR_UNSUPPORTED;
    }
}
