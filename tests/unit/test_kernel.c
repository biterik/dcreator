/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* test_kernel.c - verifies the displacement kernel produces the
 * expected rigid-shift result for atoms far from the core in a simple
 * axis-aligned configuration.
 *
 * Geometry:
 *   - reference axes == sample axes (identity rotation)
 *   - line vector    = (0, 0, 1)  (along z)
 *   - glide plane n  = (0, 1, 0)  (plane y = 0)
 *   - burgers vector = (b, 0, 0)  with b = 2.5
 *   - origin         = (0, 0, 0)
 *
 * Dislocation frame: z_d = l = (0,0,1), y_d = n = (0,1,0),
 *                    x_d = y_d x z_d = (1,0,0).
 * So the disloc frame coincides with the sample frame.
 *
 * For `upper_half` mode with `rigid_side=right` (defaults):
 *   atom at (+10, +5, 0)  -> +b  in x   (upper half, rigid region)
 *   atom at (+10, -5, 0)  -> unchanged  (lower half untouched)
 *   atom at (-10, +5, 0)  -> unchanged  (region I)
 *   atom at ( 0,  +5, 0)  -> +b/2 in x  (middle of linear ramp; t=0, frac=0.5)
 */

#include "dcreator.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(cond, msg) do {                                          \
    if (!(cond)) {                                                     \
        fprintf(stderr, "FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); \
        return 1;                                                      \
    }                                                                  \
} while (0)

static int close_eq(double a, double b) { return fabs(a - b) < 1e-9; }

static int run_case(dc_mode_t mode, dc_rigid_side_t side,
                    const double expected_dx[4], const double expected_dy[4])
{
    const double b = 2.5;
    const double stripe = 10.0;

    dc_params_t p;
    dc_params_init(&p);
    p.new_x = (dc_vec3_t){ 1, 0, 0 };
    p.new_z = (dc_vec3_t){ 0, 0, 1 };
    p.burgers = (dc_vec3_t){ b, 0, 0 };
    p.line    = (dc_vec3_t){ 0, 0, 1 };
    p.glide_plane_normal = (dc_vec3_t){ 0, 1, 0 };
    p.origin = (dc_vec3_t){ 0, 0, 0 };
    p.stripe = stripe;
    p.mode = mode;
    p.rigid_side = side;

    dc_geometry_t g;
    char err[512] = {0};
    dc_status_t rc = dc_geometry_build(&p, &g, err, sizeof(err));
    CHECK(rc == DC_OK, err);

    dc_atoms_t a;
    dc_atoms_init(&a);
    a.n = 4;
    dc_atoms_reserve(&a, 4);
    /* (x, y, z) for 4 atoms */
    double x0[4] = { +10.0, +10.0, -10.0,   0.0 };
    double y0[4] = {  +5.0,  -5.0,  +5.0,  +5.0 };
    double z0[4] = {   0.0,   0.0,   0.0,   0.0 };
    for (int i = 0; i < 4; ++i) {
        a.id[i] = i + 1;
        a.type[i] = 1;
        a.x[i] = x0[i];
        a.y[i] = y0[i];
        a.z[i] = z0[i];
    }

    dc_apply_displacement(&a, &p, &g);

    for (int i = 0; i < 4; ++i) {
        double dx = a.x[i] - x0[i];
        double dy = a.y[i] - y0[i];
        double dz = a.z[i] - z0[i];
        if (!close_eq(dx, expected_dx[i]) ||
            !close_eq(dy, expected_dy[i]) ||
            !close_eq(dz, 0.0)) {
            fprintf(stderr,
                "FAIL: atom %d expected dx=%.3g dy=%.3g got dx=%.6g dy=%.6g dz=%.6g\n",
                i, expected_dx[i], expected_dy[i], dx, dy, dz);
            dc_atoms_free(&a);
            return 1;
        }
    }
    dc_atoms_free(&a);
    return 0;
}

int main(void)
{
    const double b = 2.5;
    /* upper_half + right: upper(+y) shifted by +b in rigid region */
    double exp_dx_upper[4] = { +b, 0.0, 0.0, +0.5*b };
    double exp_dy[4]       = { 0.0, 0.0, 0.0,  0.0 };
    if (run_case(DC_MODE_UPPER_HALF, DC_RIGID_RIGHT, exp_dx_upper, exp_dy)) return 1;

    /* lower_half + right: lower(-y) shifted by -b */
    double exp_dx_lower[4] = { 0.0, -b, 0.0, 0.0 };
    if (run_case(DC_MODE_LOWER_HALF, DC_RIGID_RIGHT, exp_dx_lower, exp_dy)) return 1;

    /* symmetric + right: upper -b/2, lower +b/2 (legacy) */
    double exp_dx_sym[4] = { -0.5*b, +0.5*b, 0.0, -0.25*b };
    if (run_case(DC_MODE_SYMMETRIC, DC_RIGID_RIGHT, exp_dx_sym, exp_dy)) return 1;

    /* upper_half + left: the rigid region is now on the LEFT, so the
     * atom at x=-10 gets the rigid shift and x=+10 is untouched. */
    double exp_dx_left[4] = { 0.0, 0.0, +b, +0.5*b };
    if (run_case(DC_MODE_UPPER_HALF, DC_RIGID_LEFT, exp_dx_left, exp_dy)) return 1;

    printf("OK\n");
    return 0;
}
