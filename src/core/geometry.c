/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* geometry.c - build the reference->sample and sample->dislocation matrices */

#include "internal.h"
#include <math.h>

/* Fill the 3x3 R_ref2act matrix from two non-zero basis vectors among
 * new_x, new_y, new_z. The third is reconstructed by cross product so
 * the resulting basis is right-handed and orthonormal. */
static dc_status_t build_ref2act(const dc_params_t *p, double R[3][3],
                                 char *errbuf, size_t errbuf_sz)
{
    dc_vec3_t nx = p->new_x;
    dc_vec3_t ny = p->new_y;
    dc_vec3_t nz = p->new_z;

    int have_x = !dc_is_zero(nx);
    int have_y = !dc_is_zero(ny);
    int have_z = !dc_is_zero(nz);
    int nhave = have_x + have_y + have_z;

    if (nhave < 2) {
        dc_errf(errbuf, errbuf_sz,
                "need at least two of new_x/new_y/new_z to be non-zero");
        return DC_ERR_PARAM;
    }

    if (have_x && have_z) {
        nx = dc_normalize(nx);
        nz = dc_normalize(nz);
        ny = dc_cross(nz, nx);
    } else if (have_x && have_y) {
        nx = dc_normalize(nx);
        ny = dc_normalize(ny);
        nz = dc_cross(nx, ny);
    } else /* have_y && have_z */ {
        ny = dc_normalize(ny);
        nz = dc_normalize(nz);
        nx = dc_cross(ny, nz);
    }

    /* Even if the third axis was explicitly provided, recomputing it
     * guarantees a right-handed orthonormal basis when the user passed
     * consistent but slightly non-unit vectors. */
    ny = dc_normalize(ny);
    nx = dc_normalize(nx);
    nz = dc_normalize(nz);

    /* Check orthogonality of the provided pair to catch bogus inputs. */
    double cross_check = 0.0;
    if (have_x && have_z) cross_check = fabs(dc_dot(p->new_x, p->new_z)
                                             / (dc_vabs(p->new_x) * dc_vabs(p->new_z)));
    else if (have_x && have_y) cross_check = fabs(dc_dot(p->new_x, p->new_y)
                                             / (dc_vabs(p->new_x) * dc_vabs(p->new_y)));
    else cross_check = fabs(dc_dot(p->new_y, p->new_z)
                             / (dc_vabs(p->new_y) * dc_vabs(p->new_z)));
    if (cross_check > 1.0e-6) {
        dc_errf(errbuf, errbuf_sz,
                "supplied sample axes are not orthogonal (cos=%.3e)", cross_check);
        return DC_ERR_GEOMETRY;
    }

    R[0][0] = nx.x; R[0][1] = nx.y; R[0][2] = nx.z;
    R[1][0] = ny.x; R[1][1] = ny.y; R[1][2] = ny.z;
    R[2][0] = nz.x; R[2][1] = nz.y; R[2][2] = nz.z;

    if (fabs(dc_det(R) - 1.0) > 1.0e-6) {
        dc_errf(errbuf, errbuf_sz,
                "R_ref2act determinant = %.6e, expected 1.0", dc_det(R));
        return DC_ERR_GEOMETRY;
    }
    return DC_OK;
}

dc_status_t dc_geometry_build(const dc_params_t *p, dc_geometry_t *g,
                              char *errbuf, size_t errbuf_sz)
{
    dc_status_t s = build_ref2act(p, g->R_ref2act, errbuf, errbuf_sz);
    if (s != DC_OK) return s;

    /* Burgers vector: keep full magnitude. Line and normal: unit. */
    dc_vec3_t b_ref = p->burgers;
    double bmag = dc_vabs(b_ref);
    if (bmag == 0.0) {
        dc_errf(errbuf, errbuf_sz, "burgers vector is zero");
        return DC_ERR_PARAM;
    }
    dc_vec3_t l_ref = dc_normalize(p->line);
    dc_vec3_t n_ref = dc_normalize(p->glide_plane_normal);
    if (dc_is_zero(l_ref) || dc_is_zero(n_ref)) {
        dc_errf(errbuf, errbuf_sz, "line vector or glide plane normal is zero");
        return DC_ERR_PARAM;
    }

    /* Transform to sample frame. */
    dc_vec3_t b_s = dc_mv(g->R_ref2act, b_ref);
    dc_vec3_t l_s = dc_mv(g->R_ref2act, l_ref);
    dc_vec3_t n_s = dc_mv(g->R_ref2act, n_ref);

    /* Dislocation frame: z = line, y = glide plane normal, x = y x z. */
    dc_vec3_t z_d = dc_normalize(l_s);
    dc_vec3_t y_d = dc_normalize(n_s);
    dc_vec3_t x_d = dc_cross(y_d, z_d);
    x_d = dc_normalize(x_d);

    /* Sanity: y and z should be perpendicular. */
    if (fabs(dc_dot(y_d, z_d)) > 1.0e-6) {
        dc_errf(errbuf, errbuf_sz,
                "line vector is not perpendicular to glide plane normal "
                "(cos=%.3e)", dc_dot(y_d, z_d));
        return DC_ERR_GEOMETRY;
    }

    g->R_act2dis[0][0] = x_d.x; g->R_act2dis[0][1] = x_d.y; g->R_act2dis[0][2] = x_d.z;
    g->R_act2dis[1][0] = y_d.x; g->R_act2dis[1][1] = y_d.y; g->R_act2dis[1][2] = y_d.z;
    g->R_act2dis[2][0] = z_d.x; g->R_act2dis[2][1] = z_d.y; g->R_act2dis[2][2] = z_d.z;

    /* Transpose = inverse for orthonormal matrix. */
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            g->R_dis2act[i][j] = g->R_act2dis[j][i];

    if (fabs(dc_det(g->R_act2dis) - 1.0) > 1.0e-6) {
        dc_errf(errbuf, errbuf_sz,
                "R_act2dis determinant = %.6e, expected 1.0",
                dc_det(g->R_act2dis));
        return DC_ERR_GEOMETRY;
    }

    g->b_sample     = b_s;  /* full magnitude preserved */
    g->bmag         = bmag;
    g->origin_sample= p->origin;

    /* Character classification via unit b and unit l. */
    dc_vec3_t b_unit = dc_scale(b_ref, 1.0 / bmag);
    double bl = dc_dot(b_unit, l_ref);
    g->is_edge  = (fabs(bl) < 1.0e-6);
    g->is_screw = (fabs(fabs(bl) - 1.0) < 1.0e-6);
    double c = bl;
    if (c >  1.0) c =  1.0;
    if (c < -1.0) c = -1.0;
    g->angle_bl_deg = acos(c) * 180.0 / 3.14159265358979323846;

    return DC_OK;
}
