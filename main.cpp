#include <iostream>
//#include <EstimateVelocity.h>
#include "EstimateVelocity.h"

int main() {    
//    EstimateVelocity a;
//    MatrixXf pos(4, 2);
//    MatrixXf flow(4, 2);
//    MatrixXf depth(4, 1);
//    MatrixXf K(3, 3);
//
//    pos << 1, 2, 3, 4, 5, 6, 7, 8;
//    flow << 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4;
//    depth << 1, 1, 1, 1;
//    K << 1, 0, 0, 0, 1, 0, 0, 0, 1;
//
//    std::tuple<MatrixXf, MatrixXf, MatrixXf> result;
//    result = a.velocityFromFlow(flow, K, depth, pos);
//    MatrixXf X, A, b;
//    std::tie(X, A, b) = result;
//    std::cout << X << std::endl;
//    std::cout << "Ransac started" <<std::endl;
//    std::vector<int> inlier = a.RANSAC(100,0.3,flow,K,depth, pos);
//    std::cout << "Ransac Finished" << std::endl;

    std::ifstream dataFile("/home/haumin/Documents/EventCamera/FilterOpticalFlow/data.txt");
//    dataFile.open("sample_data.txt");



    if(!dataFile){
        std::cout << "Unable to open the data file!";
        exit(1);
    }

    std::string line;
    int num_frame = 0;
    int num_pts = 0;
    int count = 0;
    bool first_line = true;
    std::vector<MatrixXf> all_frames;

    while(std::getline(dataFile, line)){
        if(first_line){ // find the number of frame line
            num_frame = std::stoi(line);
            std::cout << "Number of frames: " << num_frame << std::endl;
            first_line = false;
            continue;
        }
        std::size_t frame_signal = line.find(":");
        if(frame_signal != std::string::npos){ // find the number of points in certain frame line
            ++count;
            num_pts = std::stoi(line.substr(0, frame_signal));
            std::cout <<"Frame " << count << " Number of points: " << num_pts << std::endl;
            continue;
        }
        // enter the section read data
        int line_count = 0;
        MatrixXf frame_data(num_pts, 5);

        while(line_count < num_pts){
            int start_pos = 0;
            int pts_count = 0;
            std::size_t pts_signal = line.substr(start_pos).find(" ");

            while(pts_count < 5){
//                std::cout <<"Read Data: " << std::stof(line.substr(start_pos, pts_signal)) << std::endl;
//                std::cout <<"start: " << start_pos  << " signal: " << pts_signal<<std::endl;
                frame_data(line_count, pts_count) = std::stof(line.substr(start_pos, pts_signal));
                start_pos = int(pts_signal) + 1;
                pts_signal = line.find(" ", start_pos);
                ++pts_count;
            }

            ++line_count;
            if(line_count < num_pts){
                std::getline(dataFile, line);
            }
        }
//        std::cout <<"Read Matrix: " << std::endl << frame_data << std::endl;
        all_frames.push_back(frame_data);
    }


    // run velocity estimation on all data.
    int num_frames = all_frames.size();

    EstimateVelocity a;
    MatrixXf K(3, 3);
    K << 311.0520, 0, 201.8724,
            0, 311.3885, 113.6210,
            0, 0, 1;
    int i;
    for (i = 0; i < num_frames; i++) {
        MatrixXf one_frame = all_frames[i];
        int num_points = one_frame.rows();
        MatrixXf flow = one_frame.block(0, 0, num_points, 2);
        MatrixXf position = one_frame.block(0, 2, num_points, 2);
        MatrixXf depth = one_frame.col(4);
//        std::cout << one_frame << std::endl;
//        std::cout <<"Flow: " << flow << std::endl << " Depth: "<< depth << std::endl << " Position: " << position << std::endl;

        MatrixXf vel, A, b;

        std::tie(vel, A, b) = a.velocityFromFlow(flow, K, depth, position);
        std::cout << "A: " << std::endl << A << std::endl;

//        std::cout << "A: " << A << std::endl <<"B: " << std::endl << b << std::endl << "Vel: " << vel << std::endl;

//
        MatrixXf velocity = a.RANSAC(100, 100, flow, K, depth, position);
//        std::cout << "Vel: " << vel << std::endl;
        break;
    }

    return 0;
}