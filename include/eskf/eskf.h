//
// Created by Haumin Zhu on 5/25/18.
//

#ifndef FILTEROPTICALFLOW_ESKF_H
#define FILTEROPTICALFLOW_ESKF_H

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <set>
#include <tuple>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include "include/estimate_velocity/estimate_velocity.h"

using namespace Eigen;

class ESKF{
public:

    ESKF(){

    }

    ~ESKF(){

    }

    void filter(float dt, VectorXf imu, MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position);

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

    void prediction(float dt, VectorXf imu);

    void measurement(MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position);

    Vector4f quatMultiply(Vector4f q, Vector4f r);

    Matrix3f skew(Vector3f vec);

};

#endif //FILTEROPTICALFLOW_ESKF_H
