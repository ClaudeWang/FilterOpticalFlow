//
// Created by Haumin Zhu on 4/23/18.
//

#include "estimate_velocity/estimate_velocity.h"

namespace filter_optical_flow{

    EstimateVelocity::EstimateVelocity(ros::NodeHandle nh, ros::NodeHandle nh_priv):
        nh_(nh), 
        left_undistorter_(nh_priv, "cam0"),
        right_undistorter_(nh_priv, "cam1")
    {
        right_undistorter_.caminfo(&fx_, &fy_, &px_, &py_, &baseline_);
        start_flag_ = 0;
        feature_pos_sub_ = nh_.subscribe<msckf_vio::CameraMeasurement>("feature_pos", 10, &EstimateVelocity::positionCallback, this);
        ROS_INFO("EstimateVelocity initialized");
    }

    void EstimateVelocity::positionCallback(const msckf_vio::CameraMeasurement::ConstPtr& position_msg){
        if(start_flag_ == 0){
            for(int i = 0; i < position_msg->features.size(); i++){
                id_ind_map_[position_msg->features[i].id] = i;
            }
            prev_position_msg_ = *position_msg;
            start_flag_ = 1;
            ROS_INFO("First frame initialized");
            
        } 
        else{
            Eigen::MatrixXf positions(position_msg->features.size(),2);
            Eigen::MatrixXf depths(position_msg->features.size(),1);
            Eigen::MatrixXf flows(position_msg->features.size(),2);
            for(int i = 0; i < position_msg->features.size(); i++){
                if(id_ind_map_.count(position_msg->features[i].id)){
                    positions(i, 0) = position_msg->features[i].u0;
                    positions(i, 1) = position_msg->features[i].v0;
                    
                    float disparity = position_msg->features[i].u1 - position_msg->features[i].u0;
                    //std::cout<<disparity<<std::endl;
                    float depth = -baseline_/(disparity * fx_);
                    depths(i) = depth;

                    
                    flows(i, 0) = 
                        (position_msg->features[i].u0 - 
                         prev_position_msg_.features[id_ind_map_[position_msg->features[i].id]].u0)/
                        (position_msg->header.stamp.toSec() - prev_position_msg_.header.stamp.toSec());
                    flows(i, 1) = 
                        (position_msg->features[i].v0 - 
                         prev_position_msg_.features[id_ind_map_[position_msg->features[i].id]].v0)/
                        (position_msg->header.stamp.toSec() - prev_position_msg_.header.stamp.toSec());
                }
                id_ind_map_[position_msg->features[i].id] = i;
            }        
            
            prev_position_msg_ = *position_msg;
            // std::cout<<"flow"<<flows<<std::endl;
            // std::cout<<"depths"<<depths<<std::endl;
            // std::cout<<"positions"<<positions<<std::endl;

            ///////////////////////TO DO///////////////////////
            velocityFromFlow(flows, depths, positions);
        }
    }

    std::tuple<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf> EstimateVelocity::velocityFromFlow(Eigen::MatrixXf flow, Eigen::MatrixXf depth, Eigen::MatrixXf position) {
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

        Eigen::MatrixXf h_pts(3, num_pts);
        Eigen::MatrixXf new_pos(num_pts, 2);
        h_pts.block(0, 0, 2, num_pts) = position.transpose().block(0, 0, 2, num_pts);
        h_pts.block(2, 0, 1, num_pts) = Eigen::MatrixXf::Ones(1, num_pts).block(0, 0, 1, num_pts);
        h_pts = h_pts.array().rowwise() / h_pts.row(2).array();

        // extract the points in ideal coordinates.
        new_pos = h_pts.transpose().block(0, 0, num_pts, 2);

        // Construct A matrix and b matrix to solve the linear/angular velocity.
        Eigen::MatrixXf ones;
        Eigen::MatrixXf zeros;
        ones = Eigen::MatrixXf::Ones(num_pts, 1);
        zeros = Eigen::MatrixXf::Zero(num_pts, 1);
        Eigen::MatrixXf A_upper(num_pts, 6);
        Eigen::MatrixXf A_lower(num_pts, 6);

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

        Eigen::MatrixXf A(2 * num_pts, 6);
        A.block(0, 0, num_pts, 6) = A_upper.block(0, 0, num_pts, 6);
        A.block(num_pts, 0, num_pts, 6) = A_lower.block(0, 0, num_pts, 6);

        // b matrix
        Eigen::MatrixXf b(2 * num_pts, 1);
        b.block(0, 0, num_pts, 1) = flow.col(0).block(0, 0, num_pts, 1);
        b.block(num_pts, 0, num_pts, 1) = flow.col(1).block(0, 0, num_pts, 1);

        // solve
        int start = clock();
        Eigen::MatrixXf X(6, 1);
        Eigen::ColPivHouseholderQR<Eigen::MatrixXf> dec(A);
        X.col(0) = dec.solve(b);

    //    std::cout << (clock() - start) / double(CLOCKS_PER_SEC) << std::endl;
    //    std::cout << X << std::endl;

        return std::make_tuple(X, A, b);
    }

    std::vector<int> EstimateVelocity::ransac(int iterations, float threshold, Eigen::MatrixXf flow, Eigen::MatrixXf depth, Eigen::MatrixXf position) {
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

        Eigen::MatrixXf X, A, b; // X is not used
        std::tie(X, A, b) = velocityFromFlow(flow, depth, position);

        // RANSAC
        for(int i = 0; i < iterations; i ++){
            // randomly select points
            std::set<int> rand_idx = generateRandomNum(num_rand_pts, num_pts);
            std::set<int>::iterator iterator;

            Eigen::MatrixXf select_postition, select_flow, select_depth;
            select_postition = Eigen::MatrixXf::Zero(num_rand_pts, 2);
            select_flow = Eigen::MatrixXf::Zero(num_rand_pts, 2);
            select_depth = Eigen::MatrixXf::Zero(num_rand_pts, 1);

            int c = 0;
            for(iterator = rand_idx.begin(); iterator != rand_idx.end() && c < num_rand_pts; iterator++){
                int idx = *iterator;
                select_postition.row(c) = position.row(idx);
                select_flow.row(c) = flow.row(idx);
                select_depth.row(c) = depth.row(idx);
                c ++;
            }

            Eigen::MatrixXf vel, A_rand, b_rand;
            std::tie(vel, A_rand, b_rand) = velocityFromFlow(select_flow, select_depth, select_postition);

            Eigen::MatrixXf diff = A * vel - b;
            Eigen::MatrixXf x_diff, y_diff;
            x_diff = Eigen::MatrixXf::Zero(num_pts, 1);
            y_diff = Eigen::MatrixXf::Zero(num_pts, 1);
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

}