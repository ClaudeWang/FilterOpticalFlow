//
// Created by claude on 4/21/18.
//

#ifndef FILTEROPTICALFLOW_ESTIMATEVELOCITY_H
#define FILTEROPTICALFLOW_ESTIMATEVELOCITY_H

#include <eigen3/Eigen/Dense>
#include <tuple>
#include <iostream>

using namespace Eigen;

class EstimateVelocity {

public:

    EstimateVelocity() {

    }
    ~EstimateVelocity() {
    }

    MatrixXf velocityFromFlow(MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position);

};


#endif //FILTEROPTICALFLOW_ESTIMATEVELOCITY_H
