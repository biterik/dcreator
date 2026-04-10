/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* dcreator CLI main entry point. */

#include "dcreator.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef DCREATOR_HAVE_OPENMP
#  include <omp.h>
#endif

static dc_status_t dispatch_read(const dc_params_t *p, dc_atoms_t *a,
                                 char *err, size_t errsz)
{
    dc_format_t fmt = p->input_format;
    if (fmt == DC_FORMAT_AUTO) {
        fmt = dc_detect_format(p->input_file);
        if (fmt == DC_FORMAT_UNKNOWN) {
            snprintf(err, errsz, "cannot auto-detect format of '%s'",
                     p->input_file);
            return DC_ERR_FORMAT;
        }
    }
    switch (fmt) {
    case DC_FORMAT_LAMMPS_DATA:
        return dc_read_lammps_data(p->input_file, p->lammps_style, a, err, errsz);
    case DC_FORMAT_LAMMPS_DUMP:
        return dc_read_lammps_dump(p->input_file, a, err, errsz);
    case DC_FORMAT_IMD:
        return dc_read_imd(p->input_file, a, err, errsz);
    default:
        snprintf(err, errsz, "unsupported input format");
        return DC_ERR_UNSUPPORTED;
    }
}

static dc_status_t dispatch_write(const dc_params_t *p, const dc_atoms_t *a,
                                  char *err, size_t errsz)
{
    switch (p->output_format) {
    case DC_FORMAT_LAMMPS_DATA:
        return dc_write_lammps_data(p->output_file, a, err, errsz);
    case DC_FORMAT_IMD:
        return dc_write_imd(p->output_file, a, err, errsz);
    default:
        snprintf(err, errsz, "unsupported output format");
        return DC_ERR_UNSUPPORTED;
    }
}

static void print_usage(const char *argv0)
{
    fprintf(stderr,
        "dcreator 2.0 - parallel Volterra dislocation inserter\n"
        "(C) 2000-2026 Erik Bitzek, GPL-3.0-or-later, NO WARRANTY.\n"
        "Cite: Vaid et al., Acta Mater. 236 (2022) 118095\n"
        "      https://doi.org/10.1016/j.actamat.2022.118095\n\n"
        "Usage: %s <parameter_file>\n\n"
        "Parameter keys (superset of legacy dcreator):\n"
        "  configfile, outfile\n"
        "  new_x, new_y, new_z           sample orientation (any two required)\n"
        "  burgersv, linev, glideplane   dislocation vectors\n"
        "  point                         point on dislocation line\n"
        "  stripe                        width of linear ramp region\n"
        "  displacement_mode             upper_half | lower_half | symmetric\n"
        "                                (default: upper_half)\n"
        "  rigid_side                    right | left  (default: right)\n"
        "  input_format                  auto | lammps_data | lammps_dump | imd\n"
        "  output_format                 lammps_data | imd\n"
        "  lammps_atom_style             atomic | charge | molecular | full | bond\n"
        "  num_threads                   OpenMP thread count (0 = default)\n"
        , argv0);
}

int main(int argc, char **argv)
{
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_usage(argv[0]);
        return (argc < 2) ? 2 : 0;
    }

    dc_params_t params;
    dc_params_init(&params);

    char err[1024] = {0};
    dc_status_t rc;

    rc = dc_params_load(&params, argv[1], err, sizeof(err));
    if (rc != DC_OK) {
        fprintf(stderr, "dcreator: %s\n", err[0] ? err : dc_strerror(rc));
        return 2;
    }
    rc = dc_params_validate(&params, err, sizeof(err));
    if (rc != DC_OK) {
        fprintf(stderr, "dcreator: %s\n", err[0] ? err : dc_strerror(rc));
        return 2;
    }

    printf("dcreator 2.0 - parallel Volterra dislocation inserter\n");
    printf("  (C) 2000-2026 Erik Bitzek, GPL-3.0-or-later, NO WARRANTY.\n");
    printf("  If you publish results obtained with dcreator please cite:\n");
    printf("    Vaid, Wei, Bitzek, Nasiri, Zaiser, Acta Mater. 236 (2022) 118095.\n");
    printf("    https://doi.org/10.1016/j.actamat.2022.118095\n\n");
    printf("  input  : %s\n", params.input_file);
    printf("  output : %s (%s)\n", params.output_file,
           dc_format_name(params.output_format));
    printf("  mode   : %s\n", dc_mode_name(params.mode));
    printf("  rigid  : %s\n", dc_rigid_side_name(params.rigid_side));
    printf("  stripe : %.6g\n", params.stripe);
#ifdef DCREATOR_HAVE_OPENMP
    int nthr = params.num_threads > 0 ? params.num_threads : omp_get_max_threads();
    printf("  threads: %d (OpenMP)\n", nthr);
#else
    printf("  threads: 1 (serial)\n");
#endif

    /* Build geometry. */
    dc_geometry_t geo;
    memset(&geo, 0, sizeof(geo));
    rc = dc_geometry_build(&params, &geo, err, sizeof(err));
    if (rc != DC_OK) {
        fprintf(stderr, "dcreator: %s\n", err[0] ? err : dc_strerror(rc));
        return 3;
    }
    printf("  |b|    : %.6g\n", geo.bmag);
    printf("  char   : %s (angle(b,l)=%.3f deg)\n",
           geo.is_screw ? "screw" : geo.is_edge ? "edge" : "mixed",
           geo.angle_bl_deg);

    /* Read atoms. */
    dc_atoms_t atoms;
    dc_atoms_init(&atoms);
    rc = dispatch_read(&params, &atoms, err, sizeof(err));
    if (rc != DC_OK) {
        fprintf(stderr, "dcreator: %s\n", err[0] ? err : dc_strerror(rc));
        dc_atoms_free(&atoms);
        return 4;
    }
    printf("  atoms  : %zu (source format: %s)\n",
           atoms.n, dc_format_name(atoms.source_format));

    /* Apply displacement. */
    dc_apply_displacement(&atoms, &params, &geo);

    /* Write. */
    rc = dispatch_write(&params, &atoms, err, sizeof(err));
    if (rc != DC_OK) {
        fprintf(stderr, "dcreator: %s\n", err[0] ? err : dc_strerror(rc));
        dc_atoms_free(&atoms);
        return 5;
    }

    dc_atoms_free(&atoms);
    printf("done.\n");
    return 0;
}
