///
/// @package phosim
/// @file perturbation.h
/// @brief header file for perturbation class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#ifndef Perturbation_H
#define Perturbation_H

class Perturbation {

 public:
    std::vector<double> eulerPhi;
    std::vector<double> eulerPsi;
    std::vector<double> eulerTheta;
    std::vector<double> decenterX;
    std::vector<double> decenterY;
    std::vector<double> defocus;
    double *zernike_r_grid;
    double *zernike_phi_grid;
    double *zernike_x_grid;
    double *zernike_y_grid;
    double *zernike_a_grid;
    double *zernike_b_grid;
    double *zernike_coeff;
    double *zernike_summed;
    int *zernike_summed_flag;
    double *zernike_summed_nr_p;
    double *zernike_summed_np_r;
    double *scatteringprofile;
    double *scatteringsigma;
    double *scatteringb;
    double *scatteringc;
    double *rotationmatrix;
    double *inverserotationmatrix;
    double *jitterrot;
    double *jitterele;
    double *jitterazi;
    double *jittertime;
    double *windshake;
    std::vector<int> zernikeflag;
    std::vector<double> rmax;

};

#endif
