//
// Created by claude on 4/21/18.
//
#pragma once

#include <eigen3/Eigen/Dense>
#include <tuple>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <set>
#include <tuple>
#include <vector>
#include <math.h>
#include <memory>
#include <ros/ros.h>
#include <unordered_map>
#include <eigen_conversions/eigen_msg.h>
#include <tf_conversions/tf_eigen.h>
#include <nav_msgs/Odometry.h>
#include "estimate_velocity/undistorter.h"
#include "msckf_vio/CameraMeasurement.h"
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
namespace filter_optical_flow{
	class EstimateVelocity {

	Undistorter left_undistorter_, right_undistorter_;
	ros::NodeHandle nh_;
	ros::Subscriber feature_pos_sub_;
	ros::Subscriber gt_pose_sub_;
	ros::Publisher vel_from_optical_flow_pub_;

	std::unordered_map<int, int> id_ind_map_;
	float fx_, fy_, px_, py_, baseline_;
	bool start_flag_;
	bool start_gt_flag_;

	//do debug process or not
    bool do_debug_ = true;

    //debug info
    float avg_lin_diff_vel_ = 0;
    int lin_vel_count_ = 0;
    geometry_msgs::PoseStamped prev_gt_pose_;
	Eigen::VectorXf lin_ang_vel_gt_ = Eigen::VectorXf::Zero(6);

	msckf_vio::CameraMeasurement prev_position_msg_; 

    void positionCallback(const msckf_vio::CameraMeasurement::ConstPtr& position_msg);
    void poseCallback(const geometry_msgs::PoseStampedConstPtr& pose_msg);
    Eigen::VectorXf estVelFromPose(const geometry_msgs::PoseStamped& prev_pose_tf,
                                   const geometry_msgs::PoseStamped& curr_pose_tf);
	std::tuple<MatrixXf, MatrixXf, MatrixXf> velocityFromFlow(MatrixXf flow, MatrixXf depth, MatrixXf position);
	MatrixXf RANSAC(int iterations, float threshold, MatrixXf flow, MatrixXf depth, MatrixXf position);
	std::set<int> generateRandomNum(int num, int range);

	public:

	    EstimateVelocity(ros::NodeHandle nh, ros::NodeHandle nh_priv);

	    ~EstimateVelocity() {
	    }
	};
}
