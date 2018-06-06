//
// Created by Haumin Zhu on 4/23/18.
//

#include "estimate_velocity/estimate_velocity.h"

using namespace Eigen;
namespace filter_optical_flow{

    EstimateVelocity::EstimateVelocity(ros::NodeHandle nh, ros::NodeHandle nh_priv):
        nh_(nh), 
        left_undistorter_(nh_priv, "cam0"),
        right_undistorter_(nh_priv, "cam1")
    {
        right_undistorter_.caminfo(&fx_, &fy_, &px_, &py_, &baseline_);
        start_flag_ = 0;
        start_gt_flag_ = 0;
        feature_pos_sub_ = nh_.subscribe<msckf_vio::CameraMeasurement>("feature_pos", 10, &EstimateVelocity::positionCallback, this);
        gt_pose_sub_ = nh.subscribe<geometry_msgs::PoseStamped>("pose", 10, &EstimateVelocity::poseCallback, this);
        vel_from_optical_flow_pub_ = nh_.advertise<nav_msgs::Odometry>("vel_from_optical_flow", 10);
        nh_.getParam("do_debug", do_debug_);
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
            Eigen::MatrixXf vel, A, b;
            int iterations = 100;
            float threshold = 0.01;
            vel = RANSAC(iterations, 0.01, flows, depths, positions); 
            //std::cout<<vel<<std::endl;
            Eigen::Vector3f linear_velocity, angular_velocity;
            linear_velocity << vel(0), vel(1), vel(2);
            angular_velocity << vel(3), vel(4), vel(5);

            // std::cout<<linear_velocity<<std::endl;
            // std::cout<<lin_ang_vel_gt_.head(3)<<std::endl;
            if(do_debug_){
                float diff_lin_vel = (linear_velocity - lin_ang_vel_gt_.head(3)).norm();
                avg_lin_diff_vel_ = diff_lin_vel/(lin_vel_count_+1) 
                                              + avg_lin_diff_vel_*lin_vel_count_/(lin_vel_count_+1); 
                lin_vel_count_++;
                ROS_INFO("avg error of eskf:%f", avg_lin_diff_vel_);
            }

            // nav_msgs::Odometry odom;
            // odom.header.stamp = ros::Time::now();
            // tf::vectorEigenToMsg(linear_velocity, odom.twist.twist.linear);
            // tf::vectorEigenToMsg(angular_velocity, odom.twist.twist.angular);
            
            // vel_from_optical_flow_pub_.publish(odom);
        }
    }

    void EstimateVelocity::poseCallback(const geometry_msgs::PoseStampedConstPtr& pose_msg)
    {
        if(start_gt_flag_ == 0){
            prev_gt_pose_ = *pose_msg;    
            start_gt_flag_ = 1;
        } else{
            lin_ang_vel_gt_ = estVelFromPose(prev_gt_pose_,
                                             *pose_msg);
            prev_gt_pose_ = *pose_msg;
        }

    }

    // Uses numerical differentiation to esimate velocity from a set of poses (usually from GT pose).
    Eigen::VectorXf EstimateVelocity::estVelFromPose(const geometry_msgs::PoseStamped& prev_pose_tf,
                                                   const geometry_msgs::PoseStamped& curr_pose_tf) 
    {
        Eigen::VectorXf lin_ang_vel = Eigen::VectorXf::Zero(6);

        Eigen::Isometry3d prev_pose, curr_pose;
        tf::poseMsgToEigen(prev_pose_tf.pose, prev_pose);
        tf::poseMsgToEigen(curr_pose_tf.pose, curr_pose);

        Eigen::Isometry3f pose_prev_curr = (prev_pose.inverse() * curr_pose).cast<float>();
        float dt = (curr_pose_tf.header.stamp - prev_pose_tf.header.stamp).toSec();

        lin_ang_vel.head(3) += pose_prev_curr.translation() / dt;
        Eigen::Matrix3f R = pose_prev_curr.linear();
        Eigen::Matrix3f omega_hat = R.log() / dt;

        lin_ang_vel(3) += omega_hat(2, 1);
        lin_ang_vel(4) += omega_hat(0, 2);
        lin_ang_vel(5) += omega_hat(1, 0);
       
        return lin_ang_vel;
    }

    std::tuple<MatrixXf, MatrixXf, MatrixXf> EstimateVelocity::velocityFromFlow(MatrixXf flow, MatrixXf depth, MatrixXf position) {
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
    //    std::cout << "Col 1: " << h_pts.row(0) << std::endl <<"Col 2: " << std::endl << h_pts.row(1)<< std::endl;

        h_pts = h_pts.array().rowwise() / h_pts.row(2).array();

        // extract the points in ideal coordinates.
        new_pos = h_pts.transpose().block(0, 0, num_pts, 2);

        // convert the flow into ideal coordinates.
    //    flow.col(0) = flow.col(0) / K(0, 0);
    //    flow.col(1) = flow.col(1) / K(1, 1);

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
        X.col(0) = dec.solve(b); // X in camera frame, should be converted to the world frame

    //    std::cout << (clock() - start) / double(CLOCKS_PER_SEC) << std::endl;
    //    std::cout << X << std::endl;
    //    std::cout << "A: " << A << std::endl <<"B: " << std::endl << b << std::endl << "Vel: " << X << std::endl;

        return std::make_tuple(X, A, b);
    }

    MatrixXf EstimateVelocity::RANSAC(int iterations, float threshold, MatrixXf flow, MatrixXf depth, MatrixXf position) {
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
        std::tie(X, A, b) = velocityFromFlow(flow, depth, position); //********ERROR


        // RANSAC
        srand((unsigned) time(NULL));

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
            std::tie(vel, A_rand, b_rand) = velocityFromFlow(select_flow, select_depth, select_postition);

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
    //                std::cout << "Iteration: "<< i + 1 << " inliers number: " << count << std::endl;
                    break;
                };
            }
    //        std::cout << "Iteration: "<< i + 1 << " inliers number: " << count << std::endl;
        }
    //    std::cout << " max inliers number: " << max_inlier << std::endl;

        // estimate velocity
        MatrixXf final_postition, final_flow, final_depth;
        final_postition = MatrixXf::Zero(max_inlier, 2);
        final_flow = MatrixXf::Zero(max_inlier, 2);
        final_depth = MatrixXf::Zero(max_inlier, 1);

        for(int k = 0; k < max_inlier; k ++){
            final_postition.row(k) = position.row(result[k]);
            final_flow.row(k) = flow.row(result[k]);
            final_depth.row(k) = depth.row(result[k]);
        }

        MatrixXf vel_final, A_final, b_final;
        std::tie(vel_final, A_final, b_final) = velocityFromFlow(final_flow, final_depth, final_postition);

        return vel_final;

    }

    std::set<int> EstimateVelocity::generateRandomNum(int num_pts, int range){
        /**
         * num: number of random integer
         * range: range of random number (upper bound)
         */
        std::set<int> output;
    //    srand((unsigned) time(NULL));
        while(output.size() < num_pts){
            int number = rand() % range;
            output.insert(number);
        }
        return output;
    }
}