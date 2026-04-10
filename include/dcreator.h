/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * See LICENSE for full license text and COPYING.HEADER for citation info.
 */
/* dcreator.h - public API for libdcreator
 *
 * Insert a Volterra dislocation into an atomistic configuration by
 * applying a piecewise displacement field in the dislocation frame:
 *
 *     untouched region -> linear ramp (width `stripe`) -> rigid shift
 *
 * See README.md for the physical model.
 */

#ifndef DCREATOR_H
#define DCREATOR_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/*  Status codes                                                       */
/* ------------------------------------------------------------------ */
typedef enum {
    DC_OK             = 0,
    DC_ERR_IO         = 1,
    DC_ERR_FORMAT     = 2,
    DC_ERR_PARAM      = 3,
    DC_ERR_MEMORY     = 4,
    DC_ERR_GEOMETRY   = 5,
    DC_ERR_UNSUPPORTED= 6,
} dc_status_t;

const char *dc_strerror(dc_status_t s);

/* ------------------------------------------------------------------ */
/*  Enums                                                              */
/* ------------------------------------------------------------------ */
typedef enum {
    DC_FORMAT_AUTO = 0,
    DC_FORMAT_LAMMPS_DATA,
    DC_FORMAT_LAMMPS_DUMP,
    DC_FORMAT_IMD,
    DC_FORMAT_UNKNOWN,
} dc_format_t;

typedef enum {
    DC_MODE_UPPER_HALF = 0,  /* default: upper half (y_d > 0) by +b, lower untouched */
    DC_MODE_LOWER_HALF,      /* lower half by -b, upper untouched */
    DC_MODE_SYMMETRIC,       /* legacy: upper by -b/2, lower by +b/2 */
} dc_mode_t;

typedef enum {
    DC_RIGID_RIGHT = 0,      /* default: rigid shift at pos_d.x > +stripe/2 */
    DC_RIGID_LEFT,           /* rigid shift at pos_d.x < -stripe/2 */
} dc_rigid_side_t;

/* LAMMPS atom_style (we only need to know the column layout). */
typedef enum {
    DC_ATOM_STYLE_ATOMIC = 0,  /* id type x y z [extras] */
    DC_ATOM_STYLE_CHARGE,      /* id type q x y z [extras] */
    DC_ATOM_STYLE_MOLECULAR,   /* id mol type x y z [extras] */
    DC_ATOM_STYLE_FULL,        /* id mol type q x y z [extras] */
    DC_ATOM_STYLE_BOND,        /* id mol type x y z [extras] */
} dc_atom_style_t;

/* ------------------------------------------------------------------ */
/*  Geometry primitives                                                */
/* ------------------------------------------------------------------ */
typedef struct {
    double x, y, z;
} dc_vec3_t;

typedef struct {
    double xlo, xhi;
    double ylo, yhi;
    double zlo, zhi;
    double xy, xz, yz;       /* tilt factors, 0 if orthogonal */
    int    triclinic;        /* 1 if any tilt factor is non-zero */
} dc_box_t;

/* ------------------------------------------------------------------ */
/*  Parameters                                                         */
/* ------------------------------------------------------------------ */
typedef struct {
    /* IO */
    char        input_file[4096];
    char        output_file[4096];
    dc_format_t input_format;    /* DC_FORMAT_AUTO detects from content */
    dc_format_t output_format;   /* DC_FORMAT_LAMMPS_DATA or DC_FORMAT_IMD */
    dc_atom_style_t lammps_style;/* relevant for LAMMPS data I/O */

    /* sample orientation: at least two of new_x/new_y/new_z non-zero */
    dc_vec3_t   new_x;
    dc_vec3_t   new_y;
    dc_vec3_t   new_z;

    /* dislocation specification, in reference/crystal frame */
    dc_vec3_t   burgers;              /* full Burgers vector (magnitude encoded) */
    dc_vec3_t   line;                 /* line direction (unnormalized ok) */
    dc_vec3_t   glide_plane_normal;   /* glide plane normal (unnormalized ok) */
    dc_vec3_t   origin;               /* point on dislocation line, sample frame */

    double      stripe;               /* width of the linear ramp region */

    dc_mode_t       mode;
    dc_rigid_side_t rigid_side;

    /* OpenMP thread count: 0 = use OMP_NUM_THREADS / default */
    int         num_threads;
} dc_params_t;

void        dc_params_init(dc_params_t *p);
dc_status_t dc_params_load(dc_params_t *p, const char *filename,
                           char *errbuf, size_t errbuf_sz);
dc_status_t dc_params_validate(const dc_params_t *p,
                               char *errbuf, size_t errbuf_sz);

const char *dc_format_name(dc_format_t f);
const char *dc_mode_name(dc_mode_t m);
const char *dc_rigid_side_name(dc_rigid_side_t r);
const char *dc_atom_style_name(dc_atom_style_t s);

/* ------------------------------------------------------------------ */
/*  Atom storage (structure of arrays)                                 */
/* ------------------------------------------------------------------ */
typedef struct {
    size_t   n;                  /* number of atoms */
    size_t   cap;                /* allocated capacity */

    /* Core columns (always present after a read) */
    int64_t *id;                 /* atom id */
    int32_t *type;               /* atom type */
    double  *x, *y, *z;          /* positions (sample frame) */

    /* Optional columns (NULL if not present in source) */
    double  *mass;               /* per-atom mass (IMD) */
    double  *q;                  /* charge (LAMMPS charge/full) */
    int64_t *molecule;           /* molecule id (LAMMPS molecular/full/bond) */
    double  *vx, *vy, *vz;       /* velocities */

    /* Per-atom trailing text (everything after the core columns) so
     * arbitrary LAMMPS dump / data columns round-trip untouched.
     * NULL if no trailing text was present. */
    char   **extra;

    dc_box_t box;
    int      periodic[3];        /* 1 = periodic, 0 = shrink/fixed (pass-through) */

    /* Format-specific header blob for round-trip fidelity (LAMMPS data
     * header sections before "Atoms", IMD #-prefixed lines, etc.).
     * Owned by the atoms struct. */
    char    *header_blob;
    size_t   header_blob_len;

    /* Parsed LAMMPS data sub-section blob: Masses, Velocities, Bonds
     * (everything that isn't the Atoms section) preserved verbatim. */
    char    *trailer_blob;
    size_t   trailer_blob_len;

    dc_atom_style_t style;       /* layout used in the source file */
    dc_format_t     source_format;
} dc_atoms_t;

void        dc_atoms_init(dc_atoms_t *a);
void        dc_atoms_free(dc_atoms_t *a);
dc_status_t dc_atoms_reserve(dc_atoms_t *a, size_t n);

/* ------------------------------------------------------------------ */
/*  Format auto-detection                                              */
/* ------------------------------------------------------------------ */
dc_format_t dc_detect_format(const char *filename);

/* ------------------------------------------------------------------ */
/*  Readers / writers                                                  */
/* ------------------------------------------------------------------ */
dc_status_t dc_read(const char *filename, dc_format_t fmt,
                    dc_atom_style_t style, dc_atoms_t *a,
                    char *errbuf, size_t errbuf_sz);

dc_status_t dc_write(const char *filename, dc_format_t fmt,
                     const dc_atoms_t *a,
                     char *errbuf, size_t errbuf_sz);

dc_status_t dc_read_lammps_data(const char *filename, dc_atom_style_t style,
                                dc_atoms_t *a,
                                char *errbuf, size_t errbuf_sz);
dc_status_t dc_read_lammps_dump(const char *filename, dc_atoms_t *a,
                                char *errbuf, size_t errbuf_sz);
dc_status_t dc_read_imd(const char *filename, dc_atoms_t *a,
                        char *errbuf, size_t errbuf_sz);

dc_status_t dc_write_lammps_data(const char *filename, const dc_atoms_t *a,
                                 char *errbuf, size_t errbuf_sz);
dc_status_t dc_write_imd(const char *filename, const dc_atoms_t *a,
                         char *errbuf, size_t errbuf_sz);

/* ------------------------------------------------------------------ */
/*  Geometry (matrix setup)                                            */
/* ------------------------------------------------------------------ */
typedef struct {
    double    R_ref2act[3][3];   /* reference/crystal -> sample frame */
    double    R_act2dis[3][3];   /* sample -> dislocation frame */
    double    R_dis2act[3][3];   /* dislocation -> sample (transpose) */
    dc_vec3_t b_sample;          /* Burgers vector in sample frame (full magnitude) */
    dc_vec3_t origin_sample;     /* point on dislocation line in sample frame */
    double    bmag;              /* |b_sample| */
    int       is_screw;
    int       is_edge;
    double    angle_bl_deg;      /* angle(b,l) in degrees */
} dc_geometry_t;

dc_status_t dc_geometry_build(const dc_params_t *p, dc_geometry_t *g,
                              char *errbuf, size_t errbuf_sz);

/* ------------------------------------------------------------------ */
/*  Displacement kernel                                                */
/*                                                                     */
/*  Mutates a->x, a->y, a->z in place according to `p` and `g`.        */
/*  OpenMP-parallel; thread count honored from p->num_threads.          */
/* ------------------------------------------------------------------ */
void dc_apply_displacement(dc_atoms_t *a,
                           const dc_params_t *p,
                           const dc_geometry_t *g);

#ifdef __cplusplus
}
#endif

#endif /* DCREATOR_H */
