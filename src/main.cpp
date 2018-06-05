#include <iostream>
//#include <EstimateVelocity.h>
#include "include/estimate_velocity/estimate_velocity.h"
#include "include/eskf/eskf.h"

void testRansac();
void testESKF();

int main() {

//    /* main function start**/
//    testRansac();
    testESKF();
//    /**ending**/

    return 0;

}

void testRansac(){

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
            first_line = false;
            continue;
        }
        std::size_t frame_signal = line.find(":");
        if(frame_signal != std::string::npos){ // find the number of points in certain frame line
            ++count;
            num_pts = std::stoi(line.substr(0, frame_signal));
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
                std::string d = line.substr(start_pos, pts_signal - start_pos);
                frame_data(line_count, pts_count) = std::stof(d);
                start_pos = int(pts_signal) + 1;
                pts_signal = line.find(" ", start_pos);
                ++pts_count;
            }

            ++line_count;
            if(line_count < num_pts){
                std::getline(dataFile, line);
            }
        }
        all_frames.push_back(frame_data);
    }


    // run velocity estimation on all data.
    int num_frames = all_frames.size();
    int start = clock();

    EstimateVelocity a;
    MatrixXf K(3, 3);
    K << 311.0520, 0, 201.8724,
            0, 311.3885, 113.6210,
            0, 0, 1;
    int i;
    std::vector<float> vel_x;
    for (i = 0; i < num_frames; i++) {
        MatrixXf one_frame = all_frames[i];
        int num_points = one_frame.rows();
        MatrixXf flow = one_frame.block(0, 0, num_points, 2);
        MatrixXf position = one_frame.block(0, 2, num_points, 2);
        MatrixXf depth = one_frame.col(4);

        MatrixXf velocity = a.velocityFromRANSAC(100, 0.3, flow, K, depth, position);
        float vel_x_curr_frame = velocity(0,0);
        vel_x.push_back(vel_x_curr_frame);

    }
//    std::cout << "XXXXXXX" << std::endl;
//    for(int m = 0; m < vel_x.size(); m++){
//        std::cout << vel_x[m] << std::endl;
//    }

    std::cout << "Time: " << (clock() - start) / double(CLOCKS_PER_SEC) / num_frames << std::endl;
}

void testESKF(){

    std::ifstream dataFile("/home/haumin/Documents/EventCamera/FilterOpticalFlow/data_with_imu.txt");
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
    std::vector<MatrixXf> all_imu;
    std::vector<bool> all_ready;

    while(std::getline(dataFile, line)){
        if(first_line){ // find the number of frame line
            num_frame = std::stoi(line);
            first_line = false;
            continue;
        }
        std::size_t frame_signal = line.find(":");
        if(frame_signal != std::string::npos){ // find the number of points in certain frame line
            ++count;
            num_pts = std::stoi(line.substr(0, frame_signal));
            continue;
        }
        // enter the section read data
        int line_count = 0;

        if(num_pts > 0) {
            MatrixXf frame_data(num_pts, 5);

            while (line_count < num_pts) {
                int start_pos = 0;
                int pts_count = 0;
                std::size_t pts_signal = line.substr(start_pos).find(" ");

                while (pts_count < 5) {
                    frame_data(line_count, pts_count) = std::stof(line.substr(start_pos, pts_signal));
                    start_pos = int(pts_signal) + 1;
                    pts_signal = line.find(" ", start_pos);
                    ++pts_count;
                }

                ++line_count;
                if (line_count < num_pts) {
                    std::getline(dataFile, line);
                }
            }
            all_frames.push_back(frame_data);
            all_ready.push_back(true);
        } else {

            MatrixXf frame_data(1, 5); // if the # of feature points is 0, save a 1 row zero matrix
            frame_data << 0, 0, 0, 0, 0;
            all_frames.push_back(frame_data);
            all_ready.push_back(false);
        }

        std::cout << "Parsing the IMU data......" << std::endl;
        MatrixXf imu_data(1, 6);
        int start_pos = 0;
        int pts_count = 0;
        std::size_t pts_signal = line.substr(start_pos).find(" ");

        while(pts_count < 6){ // ??
            std::string d = line.substr(start_pos, pts_signal - start_pos);
            imu_data(0, pts_count) = std::stof(d);
            start_pos = int(pts_signal) + 1;
            pts_signal = line.find(" ", start_pos);
            ++pts_count;
        }

//        all_frames.push_back(frame_data);
        all_imu.push_back(imu_data);
    }


    // run velocity estimation on all data.
    int num_frames = all_frames.size();
    int num_imu = all_imu.size();
    assert(num_frames == num_imu);

}