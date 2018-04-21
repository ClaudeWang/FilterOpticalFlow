//
// Created by claude on 4/21/18.
//

#include "EstimateVelocity.h"

MatrixXf EstimateVelocity::velocityFromFlow(MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position) {

    // convert coordinates into ideal coordinates.
    int num_pts = position.rows();

    MatrixXf h_pts(3, num_pts);
    MatrixXf new_pos(num_pts, 2);
    h_pts.block(0, 0, 2, num_pts) = position.transpose().block(0, 0, 2, num_pts);
    h_pts.block(2, 0, 1, num_pts) =  MatrixXf::Ones(1, num_pts).block(0, 0, 1, num_pts);
    h_pts = K.inverse() * h_pts;
    h_pts = h_pts.array().rowwise() / h_pts.row(2).array();

    // extract the points in ideal coordinates.
    new_pos = h_pts.transpose().block(0, 0, num_pts, 2);

    // convert the flow into ideal coordinates.
    flow.col(0) = flow.col(0) / K(0, 0);
    flow.col(1) = flow.col(1) / K(1, 1);

    // Construct A matrix and b matrix to solve the linear/angular velocity.
    MatrixXf ones;
    MatrixXf zeros;
    ones = MatrixXf::Ones(num_pts, 1);
    zeros = MatrixXf::Zero(num_pts, 1);
    MatrixXf A_upper(num_pts, 6);
    MatrixXf A_lower(num_pts, 6);

    // upper matrix
    A_upper.col(0) = -1 * (ones.array().colwise() / depth.col(0).array());
    A_upper.col(1) = zeros.col(0);
    A_upper.col(2) = new_pos.col(0).array().colwise() / depth.col(0).array();
    A_upper.col(3) = new_pos.col(0).array().colwise() * new_pos.col(1).array();
    A_upper.col(4) = -1 * (1 + (new_pos.col(0).array().colwise() * new_pos.col(0).array()));
    A_upper.col(5) = new_pos.col(1);

    // lower matrix
    A_lower.col(0) = zeros.col(0);
    A_lower.col(1) = -1 * (ones.array().colwise() / depth.col(0).array());
    A_lower.col(2) = new_pos.col(1).array().colwise() / depth.col(0).array();
    A_lower.col(3) = 1 + (new_pos.col(1).array().colwise() * new_pos.col(1).array());
    A_lower.col(4) = -1 * (new_pos.col(0).array().colwise() * new_pos.col(1).array());
    A_lower.col(5) = -1 * new_pos.col(0);

    MatrixXf A(2 * num_pts, 6);
    A.block(0, 0, num_pts, 6) = A_upper.block(0, 0, num_pts, 6);
    A.block(num_pts, 0, num_pts, 6) = A_lower.block(0, 0, num_pts, 6);

    // b matrix
    MatrixXf b(2 * num_pts, 1);
    b.block(0, 0, num_pts, 1) = flow.col(0).block(0, 0, num_pts, 1);
    b.block(num_pts, 0, num_pts, 1) = flow.col(1).block(0, 0, num_pts, 1);

    // solve
    int start = clock();
    MatrixXf X(6, 1);
    ColPivHouseholderQR<MatrixXf> dec(A);
    X.col(0) = dec.solve(b);

    std::cout << (clock() - start) / double(CLOCKS_PER_SEC) << std::endl;

    return X;
}
