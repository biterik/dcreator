/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* format_detect.c - sniff the file format from the first non-blank line. */

#include "../core/internal.h"
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

dc_format_t dc_detect_format(const char *filename)
{
    FILE *f = fopen(filename, "rb");
    if (!f) return DC_FORMAT_UNKNOWN;

    char buf[1024];
    dc_format_t result = DC_FORMAT_UNKNOWN;

    /* Scan up to the first ~50 lines looking for a distinguishing
     * marker. LAMMPS data files start with a free-form comment line,
     * so we can't rely on the first line alone. */
    for (int i = 0; i < 50; ++i) {
        if (!fgets(buf, sizeof(buf), f)) break;

        /* skip pure whitespace */
        const char *s = buf;
        while (*s && isspace((unsigned char)*s)) s++;
        if (*s == '\0') continue;

        /* IMD: first non-blank starts with "#F " */
        if (i == 0 && s[0] == '#' && (s[1] == 'F' || s[1] == 'f') &&
            (s[2] == ' ' || s[2] == '\t')) {
            result = DC_FORMAT_IMD;
            break;
        }

        /* LAMMPS dump */
        if (strncasecmp(s, "ITEM:", 5) == 0) {
            result = DC_FORMAT_LAMMPS_DUMP;
            break;
        }

        /* LAMMPS data: the header has lines like "N atoms", "N bonds",
         * "xlo xhi", "Atoms". */
        if (strstr(s, " atoms") || strstr(s, " bonds") ||
            strstr(s, "xlo xhi") || strstr(s, "ylo yhi") ||
            strstr(s, "zlo zhi") || strncmp(s, "Atoms", 5) == 0) {
            result = DC_FORMAT_LAMMPS_DATA;
            break;
        }
    }
    fclose(f);
    return result;
}
