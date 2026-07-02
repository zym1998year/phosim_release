///
/// @package phosim
/// @file raytrace.h
/// @brief raytrace header file
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include "vector_type.h"

void vectorCopy(Vector vectorIn, Vector *vectorOut);
void setup_tangent (double ra, double dec, Vector *tpx, Vector *tpy, Vector *tpz);
void tangent (double ra, double dec, double *x, double *y, Vector *tpx, Vector *tpy, Vector *tpz);
double vectorDot( Vector *vector1, Vector *vector2);
void vectorAdd(Vector vectorA, Vector vectorB, Vector *vectorOut);
void vectorInit (Vector *vector);
void vectorSubtract( Vector vectorA, Vector vectorB, Vector *vectorOut);
void normalize (Vector *vector);
double modulus (Vector *vector);
void propagate(Vector *position, Vector angle, double distance);
void reflect(Vector *vectorIn, Vector normal, double *twicevidotvn);
void refract (Vector *vectorIn, Vector normal, double n_1, double n_2);
double rotateInverseX(Vector *vector, double angle);
double rotateInverseY(Vector *vector, double angle);
void shift_mu(Vector *vector, double mu, double phi);
void setupGrids(double *zernikeR, double *zernikePhi, double *zernikeX, double *zernikeY, double *zernikeA, double *zernikeB, long numelements, int coordinate);
void zernikes(double *zernike_r, double *zernike_phi, double *zernike_r_grid,
              double *zernike_phi_grid, double *zernike_normal_r, double *zernike_normal_phi,
              long numelements, long nzernikes, int coordinate);
double chebyshevT(int n, double x);
void chebyshevT2D(int n, double x, double y, double *t);
void zernikeValue(int n, double x, double y, double *t);
void chebyshevs(double *r_grid, double *phi_grid, double *chebyshev, double *chebyshev_r, double *chebyshev_phi, long nPoint, long nTerm, int coordinate);
void polys (double *r_grid, double *phi_grid, double *poly, double *poly_r, double *poly_phi, long nPoint, long nTerm, int coordinate);
