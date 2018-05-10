#include "estimate_velocity/estimate_velocity.h"
#include <ros/ros.h>

int main(int argc, char** argv){

    ros::init(argc, argv, "filter_optical_flow");
    ros::NodeHandle nh;
    ros::NodeHandle nh_priv("~");
    filter_optical_flow::EstimateVelocity estimator(nh, nh_priv);

    ros::spin();
}