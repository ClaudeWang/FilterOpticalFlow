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

namespace filter_optical_flow{
	class EstimateVelocity {

	Undistorter left_undistorter_, right_undistorter_;
	ros::NodeHandle nh_;
	ros::Subscriber feature_pos_sub_;
	ros::Publisher vel_from_optical_flow_pub_;

	std::unordered_map<int, int> id_ind_map_;
	float fx_, fy_, px_, py_, baseline_;
	bool start_flag_;

	msckf_vio::CameraMeasurement prev_position_msg_; 

	public:

	    EstimateVelocity(ros::NodeHandle nh, ros::NodeHandle nh_priv);

	    ~EstimateVelocity() {
	    }

	    void positionCallback(const msckf_vio::CameraMeasurement::ConstPtr& position_msg);
	    std::tuple<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf> velocityFromFlow(Eigen::MatrixXf flow, Eigen::MatrixXf depth, Eigen::MatrixXf position);
	    std::vector<int> ransac(int iterations, float threshold, Eigen::MatrixXf flow, Eigen::MatrixXf depth, Eigen::MatrixXf position);
	    std::set<int> generateRandomNum(int num, int range);
	};
}
