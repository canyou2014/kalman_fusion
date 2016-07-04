/*
 * Kalman.h
 *
 *  Created on: Jun 30, 2016
 *      Author: lyw
 */

#ifndef KALMAN_H_
#define KALMAN_H_
#include "Eigen/Dense"
#include <fstream>
//#include < double, cvd/gl_helpers.h>
using namespace Eigen;





class Kalman {
public:
	Kalman();
	Kalman(Matrix<double, 3, 1>& initPos, Matrix<double, 3, 1>& initVel, Matrix<double, 3, 1>& initAcc, double lambdas, double delta_t);
	void EKFTransSlam(Matrix< double, 3, 1>& slam_measurement, double delta_t);
	void EKFTransSlam(Vector3d& slam_measurement);
	void EKFTransIMU(Matrix< double, 3, 1>& acc_measurement, double delta_t);
	void EKFTransIMU(Matrix< double, 3, 1>& acc_measurement);
	void getData(Matrix< double, 10, 10>& mP, Matrix< double, 10, 1>& x_stat);
	void quat2rotation(Matrix<double, 3, 3>& rotmat, Matrix<double ,4, 1>& q);
	void quat2cbn(Matrix<double, 4, 1>& quaternion, Matrix3d& dcm);
	void getData(Vector3d& x_stat, double& scale_out);
	virtual ~Kalman();
private:
			//F, state transition matrix
			Matrix< double, 10, 10> Ft;
			//states vector and noise vector
	        Matrix< double, 10, 1> x_states, v_noise;
	    	Matrix< double, 3, 10> Hi, Hv;
	    	//Kalman gain
	    	Matrix< double, 10, 3> Ki, Kv, K;
	    	//
	        Matrix< double, 10, 10> P, Q;
	        Matrix< double, 3, 1> aw, aa;

	        Matrix3d Rwc, Rca, Rimu, Rslam;
	        double mlambda, mdelta_t;
	        const double slam_pos_noise = 0.0002, vel_noise  = 0.0001, acc_noise = 0.0001, lambda_noise = 0.008;
	        const double qslam_pos_noise = 0.002, qvel_noise  = 0.08, qacc_noise = 0.0001, qlambda_noise = 0.01;
	        const Matrix< double, 10, 10> I = Matrix< double, 10, 10>::Identity();
	        const Matrix< double, 3, 3> I3 = Matrix< double, 3, 3>::Identity();
	        Matrix<double, 10, 4> G;

};

/* namespace APTAM */

#endif /* KALMAN_H_ */
