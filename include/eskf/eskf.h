//
// Created by Haumin Zhu on 5/25/18.
//

#ifndef FILTEROPTICALFLOW_ESKF_H
#define FILTEROPTICALFLOW_ESKF_H

#include "include/estimate_velocity/estimate_velocity.h"
#include <eigen3/Eigen/Dense>
#include <tuple>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <set>
#include <tuple>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace Eigen;

class ESKF{
public:

    ESKF(){
        VectorXf q(12);
        VectorXf cov(15);
        VectorXf r(3);
        VectorXf nom_state(16);
        VectorXf err_state(15);
        VectorXf true_state(16);

        q << 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02;
        cov << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.025, 0.025, 0.005, 0.001, 0.001, 0.001;
        r << 0.1, 0.1, 0.1;
        nom_state << 0, 0, 0, 0, 0, 0, 0.9073, -0.0369, -0.0165, -0.4185, 0, 0, 0, 0, 0, 0;
        err_state << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        true_state << 0, 0, 0, 0, 0, 0, 0.9073, -0.0369, -0.0165, -0.4185, 0, 0, 0, 0, 0, 0;

        _nom_state = nom_state;
        _err_state = err_state;
        _true_state = true_state;
        _Q = q.asDiagonal();
        _cov = cov.asDiagonal();
        _R = r.asDiagonal();
        estimateVelocity_ = EstimateVelocity();
        _ransac_iteration = 100;
        _ransac_threshold = 0.3;
    }

    ~ESKF(){

    }

    void filter(float dt, MatrixXf imu, MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position, bool is_ready);

private:

    VectorXf _nom_state;
    VectorXf _err_state;
    VectorXf _true_state;
    MatrixXf _Q;
    MatrixXf _cov;
    MatrixXf _R;

    EstimateVelocity estimateVelocity_;

    int _ransac_iteration;
    float _ransac_threshold;

    void prediction(float dt, MatrixXf imu);

    void measurement(MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position, MatrixXf imu);

    Vector4f quatMultiply(Vector4f q, Vector4f r);

    Matrix3f skew(Vector3f vec);

};

#endif //FILTEROPTICALFLOW_ESKF_H
