#include <iostream>
#include "Kalman.h"
#include "Eigen/Dense"
#include <stdlib.h>
#include <fstream>
#include <string>

using namespace std;

int main()
{   /*
    for(size_t i = 0; i < 20 ; ++i)
        cout << (1 + (rand()%100)/400.0) <<endl;
    */
    ofstream in;
    string tamp;
    in.open("/home/lyw/mycode/kalman_fusion/predict.txt", std::ios::trunc);
    Vector3d x_st;
    Matrix<double, 4, 1> quatss;
    Matrix<double, 3, 3> Cbnss;
    Matrix<double, 3, 1> gvector;
    gvector << 0.0,0.0,9.8113;
    Matrix<double, 3, 1> acc, raw_accelerations, raw_xstate, mxstate, raw_vel;
    acc << 0.0,0.0,0.0;
    int count;
    ifstream imu_file("/home/lyw/mycode/kalman_fusion/imu_data.csv");
    ifstream imu_file_gt("/home/lyw/mycode/kalman_fusion/ground_data.csv");
    string line, line_gt;
    if (!imu_file.good()) {
    cout << "no imu file found at " << "imu0.csv";
    return -1;
    }

        getline(imu_file, line);
        getline(imu_file_gt, line_gt);
        getline(imu_file, line);
        getline(imu_file_gt, line_gt);
        stringstream stream_gt1(line_gt);
        string s_gt1;
        getline(stream_gt1, s_gt1, ',');
        for (int j = 0; j < 3; ++j) {
          getline(stream_gt1, s_gt1, ',');
          raw_xstate(j, 0) = stof(s_gt1);
        }
        for (int j = 0; j < 4; ++j) {
          getline(stream_gt1, s_gt1, ',');
        }
        for (int j = 0; j < 3; ++j) {
          getline(stream_gt1, s_gt1, ',');
          raw_vel(j, 0) = stof(s_gt1);
        }
    Kalman mklm(raw_xstate, raw_vel, acc, 1.0, 0.01);
    while (!imu_file.eof() && !imu_file_gt.eof())
    {
        count ++;
        getline(imu_file, line);
        getline(imu_file_gt, line_gt);
        stringstream stream(line);
        stringstream stream_gt(line_gt);
        string s;
        string s_gt;
        getline(stream, s, ',');
        getline(stream_gt, s_gt, ',');
        tamp = s;


        for (int j = 0; j < 3; ++j) {
          getline(stream, s, ',');
        }
        for (int j = 0; j < 3; ++j) {
          getline(stream_gt, s_gt, ',');
          raw_xstate(j, 0) = stof(s_gt);
        }
        if(!(count % 7)){
            for(int dd = 0; dd <3; ++dd)
            {
                mxstate(dd, 0) = raw_xstate(dd, 0) * (1 + (rand()%100)/400);
            }

            mklm.EKFTransSlam(mxstate);
        }
        for (int j = 0; j < 3; ++j) {
          getline(stream, s, ',');
          raw_accelerations(j,0) = stof(s);
        }
        for (int j = 0; j < 4; ++j) {
          getline(stream_gt, s_gt, ',');
          quatss(j,0) = stof(s_gt);
        }
        if(!(count%2)){
            mklm.quat2rotation(Cbnss, quatss);
            acc = Cbnss * raw_accelerations + gvector;
            mklm.EKFTransIMU(acc);
        }
        mklm.getData(x_st);
        in << tamp  << "\t" << x_st(0) << "\t" << x_st(1) << "\t" << x_st(2) << "\n";

        /*
        for (size_t i = 0; i < 9; i++)
             in << cbn[i] << "\t";
        for (size_t i = 0; i < 3; i++)
             in << pos[i] << "\t";
        in << "\n";
        */
    }

    imu_file.close();
    in.close();

    //mklm.EKFTransIMU(acc);
    //mklm.getData(x_st);

    return 0;
}
