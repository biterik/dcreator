/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 *
 * The piecewise displacement construction implemented here was inspired
 * by D. Rodney, Acta Mater. 52 (3) (2004) 607-614.  When using dcreator
 * in published work, please cite:
 *   A. Vaid, D. Wei, E. Bitzek, S. Nasiri, M. Zaiser,
 *   Acta Materialia 236 (2022) 118095
 *   https://doi.org/10.1016/j.actamat.2022.118095
 */
/* kernel.c - the displacement kernel, OpenMP parallel.
 *
 * For each atom:
 *
 *   p' = R_act2dis * (p - O)               (sample -> dislocation frame)
 *   t  = sgn * p'.x                        (sgn = +1 if rigid_side=right, -1 if left)
 *
 *            t < -stripe/2  :  u_d = 0
 *   |t| <= stripe/2         :  u_d = (t + stripe/2) / stripe * U(p'.y)
 *            t >  stripe/2  :  u_d =                            U(p'.y)
 *
 * where U(y) depends on the displacement mode:
 *
 *   upper_half  (default): U(y>0) = +b_d,        U(y<=0) = 0
 *   lower_half           : U(y>0) = 0,           U(y<=0) = -b_d
 *   symmetric (legacy)   : U(y>0) = -b_d/2,      U(y<=0) = +b_d/2
 *
 * b_d is the Burgers vector expressed in the dislocation frame.
 * Finally p_new = R_dis2act * (p' + u_d) + O.
 */

#include "internal.h"

#ifdef DCREATOR_HAVE_OPENMP
#  include <omp.h>
#endif

/* Evaluate U(y) * factor and accumulate into the dislocation-frame
 * displacement components. */
static inline void apply_U(double y, double factor, dc_mode_t mode,
                           const dc_vec3_t b_d,
                           double *ux, double *uy, double *uz)
{
    double s = 0.0;
    switch (mode) {
    case DC_MODE_UPPER_HALF:
        s = (y > 0.0) ? factor : 0.0;
        break;
    case DC_MODE_LOWER_HALF:
        s = (y > 0.0) ? 0.0 : -factor;
        break;
    case DC_MODE_SYMMETRIC:
        s = (y > 0.0) ? -0.5 * factor : 0.5 * factor;
        break;
    }
    *ux += s * b_d.x;
    *uy += s * b_d.y;
    *uz += s * b_d.z;
}

void dc_apply_displacement(dc_atoms_t *a,
                           const dc_params_t *p,
                           const dc_geometry_t *g)
{
    const double m00 = g->R_act2dis[0][0];
    const double m01 = g->R_act2dis[0][1];
    const double m02 = g->R_act2dis[0][2];
    const double m10 = g->R_act2dis[1][0];
    const double m11 = g->R_act2dis[1][1];
    const double m12 = g->R_act2dis[1][2];
    const double m20 = g->R_act2dis[2][0];
    const double m21 = g->R_act2dis[2][1];
    const double m22 = g->R_act2dis[2][2];

    const double t00 = g->R_dis2act[0][0];
    const double t01 = g->R_dis2act[0][1];
    const double t02 = g->R_dis2act[0][2];
    const double t10 = g->R_dis2act[1][0];
    const double t11 = g->R_dis2act[1][1];
    const double t12 = g->R_dis2act[1][2];
    const double t20 = g->R_dis2act[2][0];
    const double t21 = g->R_dis2act[2][1];
    const double t22 = g->R_dis2act[2][2];

    const double Ox = g->origin_sample.x;
    const double Oy = g->origin_sample.y;
    const double Oz = g->origin_sample.z;

    /* Burgers vector in the dislocation frame (full magnitude). */
    const dc_vec3_t b_d = dc_mv(g->R_act2dis, g->b_sample);

    const double stripe      = p->stripe;
    const double half_stripe = 0.5 * stripe;
    const double inv_stripe  = (stripe > 0.0) ? (1.0 / stripe) : 0.0;
    const double sgn         = (p->rigid_side == DC_RIGID_RIGHT) ? +1.0 : -1.0;
    const dc_mode_t mode     = p->mode;

#ifdef DCREATOR_HAVE_OPENMP
    if (p->num_threads > 0) omp_set_num_threads(p->num_threads);
#endif

    const size_t N = a->n;
    double * restrict X = a->x;
    double * restrict Y = a->y;
    double * restrict Z = a->z;

#ifdef DCREATOR_HAVE_OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < N; ++i) {
        /* shift to origin */
        double dx = X[i] - Ox;
        double dy = Y[i] - Oy;
        double dz = Z[i] - Oz;

        /* rotate to dislocation frame */
        double px = m00*dx + m01*dy + m02*dz;
        double py = m10*dx + m11*dy + m12*dz;
        double pz = m20*dx + m21*dy + m22*dz;

        /* region selector */
        double t = sgn * px;

        double ux = 0.0, uy = 0.0, uz = 0.0;
        if (t > half_stripe) {
            apply_U(py, 1.0, mode, b_d, &ux, &uy, &uz);
        } else if (t >= -half_stripe) {
            double frac = (t + half_stripe) * inv_stripe;   /* 0 .. 1 */
            apply_U(py, frac, mode, b_d, &ux, &uy, &uz);
        }
        /* else: untouched region */

        px += ux;
        py += uy;
        pz += uz;

        /* rotate back and undo origin shift */
        double nx = t00*px + t01*py + t02*pz + Ox;
        double ny = t10*px + t11*py + t12*pz + Oy;
        double nz = t20*px + t21*py + t22*pz + Oz;

        X[i] = nx;
        Y[i] = ny;
        Z[i] = nz;
    }
}
