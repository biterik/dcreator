/*
 * dcreator 2.0 - Parallel Volterra dislocation inserter
 * Copyright (C) 2000-2026 Erik Bitzek
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See LICENSE and COPYING.HEADER for full license text and citation info.
 */
/* test_geometry.c - sanity checks for the matrix / frame setup. */

#include "dcreator.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(cond, msg) do {                                         \
    if (!(cond)) {                                                    \
        fprintf(stderr, "FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__);\
        return 1;                                                     \
    }                                                                 \
} while (0)

static int close_eq(double a, double b) { return fabs(a - b) < 1e-10; }

/* Edge-dislocation-like setup matching the DCREATOR_Example parameter file. */
static int test_edge_frame(void)
{
    dc_params_t p;
    dc_params_init(&p);
    p.new_x = (dc_vec3_t){ 1.0, -1.0, 0.0 };
    p.new_z = (dc_vec3_t){ 1.0,  1.0, 1.0 };
    p.burgers = (dc_vec3_t){ 1.76, -1.76, 0.0 };
    p.line    = (dc_vec3_t){ 1.0,  1.0, -2.0 };
    p.glide_plane_normal = (dc_vec3_t){ -1.0, -1.0, -1.0 };
    p.origin = (dc_vec3_t){ 0, 0, 0 };
    p.stripe = 50.0;

    dc_geometry_t g;
    char err[512] = {0};
    dc_status_t rc = dc_geometry_build(&p, &g, err, sizeof(err));
    CHECK(rc == DC_OK, err);

    /* R_ref2act determinant == 1 */
    double det_ref = 0;
    for (int i = 0; i < 3; ++i) {
        double c = 0;
        c += g.R_ref2act[0][i] * (g.R_ref2act[1][(i+1)%3]*g.R_ref2act[2][(i+2)%3]
                                - g.R_ref2act[1][(i+2)%3]*g.R_ref2act[2][(i+1)%3]);
        det_ref += c;
    }
    CHECK(close_eq(det_ref, 1.0), "det(R_ref2act) != 1");

    /* R_act2dis * R_dis2act = I */
    double I[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                I[i][j] += g.R_act2dis[i][k] * g.R_dis2act[k][j];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double want = (i == j) ? 1.0 : 0.0;
            CHECK(close_eq(I[i][j], want), "R_act2dis * R_dis2act != I");
        }

    /* |b_sample| == |b| */
    double b_orig = sqrt(p.burgers.x*p.burgers.x +
                         p.burgers.y*p.burgers.y +
                         p.burgers.z*p.burgers.z);
    CHECK(close_eq(g.bmag, b_orig), "burgers magnitude mismatch");

    /* For an edge dislocation b.l == 0 */
    double bl = p.burgers.x*p.line.x + p.burgers.y*p.line.y + p.burgers.z*p.line.z;
    CHECK(close_eq(bl, 0.0), "test setup should be an edge dislocation");
    CHECK(g.is_edge == 1, "edge classification");

    return 0;
}

/* Screw: b parallel to l. */
static int test_screw_frame(void)
{
    dc_params_t p;
    dc_params_init(&p);
    p.new_x = (dc_vec3_t){ 1.0, 0.0, 1.0 };
    p.new_z = (dc_vec3_t){ 1.0,-1.0,-1.0 };
    p.burgers = (dc_vec3_t){ -1.76, -3.52, 1.76 };
    p.line    = (dc_vec3_t){ -1.0,  -2.0, 1.0 };
    p.glide_plane_normal = (dc_vec3_t){ 1.0, -1.0, -1.0 };
    p.origin = (dc_vec3_t){ 0, 0, 0 };
    p.stripe = 24.0;

    dc_geometry_t g;
    char err[512] = {0};
    dc_status_t rc = dc_geometry_build(&p, &g, err, sizeof(err));
    CHECK(rc == DC_OK, err);
    CHECK(g.is_screw == 1, "screw classification");
    return 0;
}

int main(void)
{
    if (test_edge_frame())  return 1;
    if (test_screw_frame()) return 1;
    printf("OK\n");
    return 0;
}
