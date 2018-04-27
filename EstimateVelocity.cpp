//
// Created by Haumin Zhu on 4/23/18.
//

#include "EstimateVelocity.h"

std::tuple<MatrixXf, MatrixXf, MatrixXf> EstimateVelocity::velocityFromFlow(MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position) {
    /**
     * INPUT:
     * - flow: Nx2 matrix for optical flow
     * - K: 3x3 calibration matrix for camera
     * - depth: Nx1 matrix for depth
     * - position: Nx2 matrix for position
     * OUTPUT:
     * - X: 6x1 velocity matrix
     * - A: 2N x 6 matrix
     * - b: 2N x 1 matrix
     */
    // convert coordinates into ideal coordinates.
    int num_pts = position.rows();

    MatrixXf h_pts(3, num_pts);
    MatrixXf new_pos(num_pts, 2);
    h_pts.block(0, 0, 2, num_pts) = position.transpose().block(0, 0, 2, num_pts);
    h_pts.block(2, 0, 1, num_pts) = MatrixXf::Ones(1, num_pts).block(0, 0, 1, num_pts);
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

//    std::cout << (clock() - start) / double(CLOCKS_PER_SEC) << std::endl;
//    std::cout << X << std::endl;

    return std::make_tuple(X, A, b);
}

std::vector<int> EstimateVelocity::RANSAC(int iterations, float threshold, MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position) {
    /**
     * INPUT:
     * - iterations: number of iterations the RANSAC will take
     * - threshold: threshold for rejecting outliers
     * - flow: Nx2 matrix for optical flow
     * - K: 3x3 calibration matrix for camera
     * - depth: Nx1 matrix for depth
     * - position: Nx2 matrix for position
     */
    int num_pts = position.rows();
    int num_rand_pts = 3;
    int max_inlier = -1;
    std::vector<int> result;

    MatrixXf X, A, b; // X is not used
    std::tie(X, A, b) = velocityFromFlow(flow, K, depth, position);

    // RANSAC
    for(int i = 0; i < iterations; i ++){
        // randomly select points
        std::set<int> rand_idx = generateRandomNum(num_rand_pts, num_pts);
        std::set<int>::iterator iterator;

        MatrixXf select_postition, select_flow, select_depth;
        select_postition = MatrixXf::Zero(num_rand_pts, 2);
        select_flow = MatrixXf::Zero(num_rand_pts, 2);
        select_depth = MatrixXf::Zero(num_rand_pts, 1);

        int c = 0;
        for(iterator = rand_idx.begin(); iterator != rand_idx.end() && c < num_rand_pts; iterator++){
            int idx = *iterator;
            select_postition.row(c) = position.row(idx);
            select_flow.row(c) = flow.row(idx);
            select_depth.row(c) = depth.row(idx);
            c ++;
        }

        MatrixXf vel, A_rand, b_rand;
        std::tie(vel, A_rand, b_rand) = velocityFromFlow(select_flow, K, select_depth, select_postition);

        MatrixXf diff = A * vel - b;
        MatrixXf x_diff, y_diff;
        x_diff = MatrixXf::Zero(num_pts, 1);
        y_diff = MatrixXf::Zero(num_pts, 1);
        x_diff.block(0, 0, num_pts, 1) = diff.block(0, 0, num_pts, 1);
        y_diff.block(0, 0, num_pts, 1) = diff.block(num_pts, 0, num_pts, 1);

        int count = 0;
        std::vector<int> result_tmp;
        for(int j = 0; j < num_pts; j ++){
            double candidate = sqrt(pow(x_diff(j, 0), 2.0) + pow(y_diff(j, 0), 2.0));
            if(candidate < threshold){
                count ++;
                result_tmp.push_back(j);
            }
        }
        if(count > max_inlier){
            max_inlier = count;
            result = result_tmp;
            if(max_inlier > (num_pts * 0.95)){
                break;
            };
        }
        std::cout << "Iteration: "<< i << " inliers number: " << count << std::endl;
    }

    return result;

}

std::set<int> EstimateVelocity::generateRandomNum(int num_pts, int range){
    /**
     * num: number of random integer
     * range: range of random number (upper bound)
     */
    std::set<int> output;
    srand((unsigned) time(NULL));
    while(output.size() < num_pts){
        int number = rand() % range;
        output.insert(number);
    }
    return output;
}
