#include <iostream>
//#include <EstimateVelocity.h>
#include "EstimateVelocity.h"

int main() {
    EstimateVelocity a;
    MatrixXf pos(4, 2);
    MatrixXf flow(4, 2);
    MatrixXf depth(4, 1);
    MatrixXf K(3, 3);

    pos << 1, 2, 3, 4, 5, 6, 7, 8;
    flow << 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4;
    depth << 1, 1, 1, 1;
    K << 1, 0, 0, 0, 1, 0, 0, 0, 1;

    MatrixXf result;
    result = a.velocityFromFlow(flow, K, depth, pos);
    return 0;
}