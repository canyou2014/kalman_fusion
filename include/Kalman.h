/*
 * Kalman.h
 *
 *  Created on: Jun 30, 2016
 *      Author: lyw
 */

#ifndef KALMAN_H_
#define KALMAN_H_
#include "Eigen/Dense"

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
	void getData(Vector3d& x_stat);
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
	        Matrix< double, 10, 10> P;
	        Matrix< double, 3, 1> aw, aa;

	        Matrix3d Rwc, Rca, Rimu, Rslam;
	        double mlambda, mdelta_t;
	        const double slam_pos_noise = 0.00002, vel_noise  = 0.00001, acc_noise = 0.000001, lambda_noise = 0.0008;
	        const Matrix< double, 10, 10> I = Matrix< double, 10, 10>::Identity();
	        const Matrix< double, 3, 3> I3 = Matrix< double, 3, 3>::Identity();

};

/* namespace APTAM */

#endif /* KALMAN_H_ */
