/*
 * Kalman.cc
 *
 *
 *  Created on: Jun 30, 2016
 *      Author: lyw
 */

#include <Kalman.h>
using namespace Eigen;

Kalman::Kalman(){
		Hv.setZero();
		Hv.block(0, 0, 3, 3) << I3;
		Hi.setZero();
		Hi.block(0, 6, 3, 3) << I3;
		mlambda = 1.0;
		mdelta_t = 0.02;
		//position state
		x_states.block(0, 0, 3, 1) << 0.0,0.0,0.0;;
		//v
		x_states.block(3, 0, 3, 1) << 0.0,0.0,0.0;;
		//a
		x_states.block(6, 0, 3, 1) << 0.0,0.0,0.0;;
		//lambda
		x_states(9, 0) = 1.0;
		//state transition matrix
		Ft.setIdentity();
		Ft.block(0, 3, 3, 3 ) = I3 * mdelta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * mdelta_t * mdelta_t / mlambda;
		Ft.block(3, 6, 3, 3 ) = I3 * mdelta_t;
		P.setIdentity();
		P.block(0, 0, 3, 3) = I3 * slam_pos_noise;
		P.block(3, 3, 3, 3) = I3 * vel_noise;
		P.block(6, 6, 3, 3) = I3 * acc_noise;
		P(9,9) = lambda_noise;
		Rimu.setIdentity();
		Rslam.setIdentity();
		Rimu *= 0.01;
		Rslam *= 0.1;
}
Kalman::Kalman(Matrix<double, 3, 1>& initPos, Matrix<double, 3, 1>& initVel, Matrix<double, 3, 1>& initAcc, double lambdas, double delta_t) {
		Hv.setZero();
		Hv.block(0, 0, 3, 3) << I3;
		Hi.setZero();
		Hi.block(0, 6, 3, 3) << I3;

		x_states.block(0, 0, 3, 1) = initPos;
		x_states.block(3, 0, 3, 1) = initVel;
		x_states.block(6, 0, 3, 1) = initAcc;
		x_states(9, 0) = lambdas;
		mlambda = lambdas;
		mdelta_t = delta_t;
		Ft.setIdentity();
		Ft.block(0, 3, 3, 3 ) = I3 * mdelta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * mdelta_t * mdelta_t / mlambda;
		Ft.block(3, 6, 3, 3 ) = I3 * mdelta_t;
		P.setIdentity();
		P.block(0, 0, 3, 3) = I3 * slam_pos_noise;
		P.block(3, 3, 3, 3) = I3 * vel_noise;
		P.block(6, 6, 3, 3) = I3 * acc_noise;
		P(9,9) = lambda_noise;
		Rimu.setIdentity();
		Rslam.setIdentity();
		Rimu *= 0.001;
		Rslam *= 0.01;

}
void Kalman::EKFTransSlam(Matrix< double, 3, 1>& slam_measurement, double delta_t){

		mlambda = x_states(9, 0);
		Ft.block(0, 3, 3, 3 ) = I3 * delta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * delta_t * delta_t / mlambda;
		Ft.block(3, 6, 3, 3 ) = I3 * delta_t;
		Matrix<double, 10, 10> P_temp;

		P_temp = Ft*P*Ft.transpose();
P_temp = (P_temp + P_temp.transpose()) / 2;
		x_states = Ft*x_states + v_noise;
		Kv = P_temp*Hv.transpose()*(Hv*P_temp*Hv.transpose() + Rslam).transpose();
		x_states = x_states + Kv*(slam_measurement - x_states.block(0, 0, 3, 1));
		P = (I - Kv*Hv)*P_temp;
		P_temp = (P_temp + P_temp.transpose()) / 2;

}
void Kalman::EKFTransSlam(Vector3d& slam_measurement){

		mlambda = x_states(9, 0);
		Ft.block(0, 3, 3, 3 ) = I3 * mdelta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * mdelta_t * mdelta_t / mlambda;
		Ft.block(3, 6, 3, 3 ) = I3 * mdelta_t;
		Matrix<double, 10, 10> P_temp;

		P_temp = Ft*P*Ft.transpose();
 P_temp = (P_temp + P_temp.transpose()) / 2;
		x_states = Ft*x_states + v_noise;
		Kv = P_temp*Hv.transpose()*(Hv*P_temp*Hv.transpose() + Rslam).transpose();
		x_states = x_states + Kv*(slam_measurement - x_states.block(0, 0, 3, 1));
		P = (I - Kv*Hv)*P_temp;
		P_temp = (P_temp + P_temp.transpose()) / 2;

}
void Kalman::EKFTransIMU(Matrix< double, 3, 1>& acc_measurement, double delta_t){

		mlambda = x_states(9, 0);
		Ft.block(0, 3, 3, 3 ) = I3 * delta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * delta_t * delta_t / mlambda;
		Ft.block(3, 6, 3, 3 ) = I3 * delta_t;
		Matrix<double, 10, 10> P_temp;

		P_temp = Ft*P*Ft.transpose();
 P_temp = (P_temp + P_temp.transpose()) / 2;
		x_states = Ft*x_states + v_noise;
		Ki = P_temp*Hi.transpose()*(Hi*P_temp*Hi.transpose() + Rimu).transpose();
		x_states = x_states + Ki*(acc_measurement - x_states.block(6, 0, 3, 1));
		P = (I - Ki*Hi)*P_temp;
		P_temp = (P_temp + P_temp.transpose()) / 2;

}
void Kalman::EKFTransIMU(Matrix< double, 3, 1>& acc_measurement){

		mlambda = x_states(9, 0);
		Ft.block(0, 3, 3, 3 ) = I3 * mdelta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * mdelta_t * mdelta_t / mlambda;
		Ft.block(3, 6, 3, 3 ) = I3 * mdelta_t;
		Matrix<double, 10, 10> P_temp;

		P_temp = Ft*P*Ft.transpose();
 P_temp = (P_temp + P_temp.transpose()) / 2;

		x_states = Ft*x_states + v_noise;
		Ki = P_temp*Hi.transpose()*(Hi*P_temp*Hi.transpose() + Rimu).transpose();
		x_states = x_states + Ki*(acc_measurement - x_states.block(6, 0, 3, 1));
		P = (I - Ki*Hi)*P_temp;
		P_temp = (P_temp + P_temp.transpose()) / 2;

}
void Kalman::getData(Matrix< double, 10, 10>& mP, Matrix< double, 10, 1>& x_stat){
	mP = P;
	x_stat = x_states;
}
void Kalman::getData(Vector3d& x_stat){
	x_stat << x_states(0,0),x_states(1,0),x_states(2,0);

}
void Kalman::quat2rotation(Matrix<double, 3, 3>& rotmat, Matrix<double ,4, 1>& q){

double p[6];

p[0]=q(0, 0)*q(0, 0);
p[1]=q(1, 0)*q(1, 0);
p[2]=q(2, 0)*q(2, 0);
p[3]=q(3, 0)*q(3, 0);

p[4]=p[1]+p[2];

if (p[0]+p[3]+p[4] != 0)
{
	p[5]=2/(p[0]+p[3]+p[4]);
}
else
{
	p[5]=0;
}

rotmat(0,0)=1-p[5]*p[4];
rotmat(1,1)=1-p[5]*(p[0]+p[2]);
rotmat(2,2)=1-p[5]*(p[0]+p[1]);

p[0]=p[5]*q(0,0);
p[1]=p[5]*q(1,0);
p[4]=p[5]*q(2,0)*q(3,0);
p[5]=p[0]*q(1,0);

rotmat(0,1)=p[5]-p[4];
rotmat(1,0)=p[5]+p[4];

p[4]=p[1]*q(3,0);
p[5]=p[0]*q(2,0);

rotmat(0,2)=p[5]+p[4];
rotmat(2,0)=p[5]-p[4];

p[4]=p[0]*q(3,0);
p[5]=p[1]*q(2,0);

rotmat(1,2)=p[5]-p[4];
rotmat(2,1)=p[5]+p[4];


}


Kalman::~Kalman() {
	// TODO Auto-generated destructor stub
}

 /* namespace APTAM */
