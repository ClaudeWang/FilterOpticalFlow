//
// Created by claude on 4/21/18.
//

#ifndef FILTEROPTICALFLOW_ESTIMATEVELOCITY_H
#define FILTEROPTICALFLOW_ESTIMATEVELOCITY_H

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

class EstimateVelocity {
public:

    EstimateVelocity() {
    }

    ~EstimateVelocity() {
    }

    std::tuple<MatrixXf, MatrixXf, MatrixXf> velocityFromFlow(MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position);
    MatrixXf velocityFromRANSAC(int iterations, float threshold, MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position);
    std::set<int> generateRandomNum(int num, int range);
};


#endif //FILTEROPTICALFLOW_ESTIMATEVELOCITY_H
