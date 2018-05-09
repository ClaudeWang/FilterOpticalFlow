#include "estimate_velocity/estimate_velocity.h"
#include <ros/ros.h>

int main(int argc, char** argv){

    ros::init(argc, argv, "filter_optical_flow");
    ros::NodeHandle nh;
    ros::NodeHandle nh_priv("~");
    filter_optical_flow::EstimateVelocity estimator(nh, nh_priv);

    ros::spin();
//     std::tuple<MatrixXf, MatrixXf, MatrixXf> result;
//     result = a.velocityFromFlow(flow, K, depth, pos);
//     MatrixXf X, A, b;
//     std::tie(X, A, b) = result;
// //    std::cout << X << std::endl;
//     std::cout << "Ransac started" <<std::endl;
//     std::vector<int> inlier = a.RANSAC(100,0.3,flow,K,depth, pos);
//     std::cout << "Ransac Finished" << std::endl;
//     return 0;
}