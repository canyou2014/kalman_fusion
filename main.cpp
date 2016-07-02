#include <iostream>
#include "Kalman.h"
#include "Eigen/Dense"
#include <stdlib.h>
#include <fstream>
#include <string>
#include <chrono>
#include <thread>
using namespace std;
using namespace std::chrono;
int main()
{   /*
    for(size_t i = 0; i < 20 ; ++i)
        cout << (1 + (rand()%100)/400.0) <<endl;
    */

    Kalman* mklm = new Kalman();

    double scale = 1.0, mscale;
    ofstream in, scale_in;
    string tamp;
    in.open("predict.txt", std::ios::trunc);
    scale_in.open("predict_scale.txt", std::ios::trunc);

    Vector3d x_st;
    Matrix<double, 3, 1> acc_bias;
    acc_bias << -0.025266,0.136696,0.075593;
    Quaternion<double> quatss;
    Matrix<double, 3, 3> Cbnss;
    Matrix<double, 3, 1> gvector;
    gvector << 0.0,0.0,9.8003;
    Matrix<double, 3, 1> acc, raw_accelerations, raw_xstate, initialstate, mxstate, raw_vel;
    acc << 0.0,0.0,0.0;
    int count;
    ifstream imu_file("imu_data.csv");
    ifstream imu_file_gt("ground_data.csv");
    string line, line_gt;
    bool sa = false;
    if (!imu_file.good()) {
    cout << "no imu file found at " << "imu0.csv";
    return -1;
    }
    /*

    */
        //first line the title
        getline(imu_file, line);
        getline(imu_file_gt, line_gt);
        //second line the initial parameter
        getline(imu_file, line);
        getline(imu_file_gt, line_gt);
        //timestamp
        stringstream stream_gt1(line_gt);
        string s_gt1;
        getline(stream_gt1, s_gt1, ',');
        stringstream stream1(line);
        string s1;
        getline(stream1, s1, ',');
        //
        for (int j = 0; j < 3; ++j) {
          getline(stream_gt1, s_gt1, ',');
          initialstate(j, 0) = stof(s_gt1);
        }

        getline(stream_gt1, s_gt1, ',');
        quatss.w() = stof(s_gt1);
        getline(stream_gt1, s_gt1, ',');
        quatss.x() = stof(s_gt1);
        getline(stream_gt1, s_gt1, ',');
        quatss.y() = stof(s_gt1);
        getline(stream_gt1, s_gt1, ',');
        quatss.z() = stof(s_gt1);
        for (int j = 0; j < 3; ++j) {
          getline(stream_gt1, s_gt1, ',');
          raw_vel(j, 0) = stof(s_gt1);
        }
        for (int j = 0; j < 3; ++j) {
          getline(stream1, s1, ',');
        }
        for (int j = 0; j < 3; ++j) {
          getline(stream1, s1, ',');
          raw_accelerations(j,0) = stof(s1);
        }


        //mklm->quat2rotation(Cbnss, quatss);
        raw_accelerations -= acc_bias;
        //mklm->quat2cbn(quatss, Cbnss);
        Cbnss = quatss.toRotationMatrix();
        acc = (Cbnss)* raw_accelerations - gvector;
    /*

    */
        mklm = new Kalman(initialstate, raw_vel, acc, 2.2, 0.01);

    while (getline(imu_file, line) && getline(imu_file_gt, line_gt))
    {
        count ++;

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

        for (int j = 0; j < 3; ++j) {
          getline(stream, s, ',');
          raw_accelerations(j,0) = stof(s);
        }
        raw_accelerations -= acc_bias;
        getline(stream_gt, s_gt, ',');
        quatss.w() = stof(s_gt);
        getline(stream_gt, s_gt, ',');
        quatss.x() = stof(s_gt);
        getline(stream_gt, s_gt, ',');
        quatss.y() = stof(s_gt);
        getline(stream_gt, s_gt, ',');
        quatss.z() = stof(s_gt);
        if(!(count % 8)){
            raw_xstate -= initialstate;
            scale = (2.0 + (rand()%100 - 50)/2000.0);
            //cout << scale << endl;
            for(int dd = 0; dd <3; ++dd)
            {

                mxstate(dd, 0) = raw_xstate(dd, 0) / scale;
            }
            mklm->EKFTransSlam(mxstate, 0.01);
            sa = true;
       }else if(!(count%2)){
            //mklm->quat2cbn(quatss, Cbnss);
            Cbnss = quatss.toRotationMatrix();
            acc = (Cbnss) * raw_accelerations - gvector;
            if(sa)
                mklm->EKFTransIMU(acc,0.01);
            else
                mklm->EKFTransIMU(acc);
            sa = false;
        }

        mklm->getData(x_st, mscale);
        in << tamp  << "\t" << x_st(0) << "\t" << x_st(1) << "\t" << x_st(2) << "\n";
        scale_in << mscale << "\n";

    }


    this_thread::sleep_for(milliseconds(1));
    imu_file.close();
    in.close();
    scale_in.close();
    //mklm.EKFTransIMU(acc);
    //mklm.getData(x_st);
    delete mklm;
    return 0;
}
