/*
 * Kalman.cc
 *
 *
 *  Created on: Jun 30, 2016
 *      Author: lyw
 */

#include <Kalman.h>

namespace APTAM {
Kalman::Kalman(){
		Hv = Zeros;
		Hv.slice<0, 0, 3, 3>() = I3;
		Hi = Zeros;
		Hi.slice<0, 6, 3, 3>() = I3;
		mlambda = 1.0;
		mdelta_t = 0.02;
		//position state
		x_states.slice<0, 0, 3, 1>() = Data(0.0, 0.0, 0.0);
		//v
		x_states.slice<3, 0, 3, 1>() = Data(0.0, 0.0, 0.0);
		//a
		x_states.slice<6, 0, 3, 1>() = Data(0.0, 0.0, 0.0);
		//lambda
		x_states(9, 0) = 1.0;
		//state transition matrix
		Ft = Identity;
		Ft.slice<0, 3, 3, 3 >() = I3 * mdelta_t / mlambda;
		Ft.slice<0, 6, 3, 3 >() = I3 * mdelta_t * mdelta_t / mlambda;
		Ft.slice<3, 6, 3, 3 >() = I3 * mdelta_t;
		P = Identity;
		P.slice<0, 0, 3, 3>() = I3 * slam_pos_noise;
		P.slice<3, 3, 3, 3>() = I3 * vel_noise;
		P.slice<6, 6, 3, 3>() = I3 * acc_noise;
		P(9,9) = lambda_noise;
}
Kalman::Kalman(Matrix<3, 1, double>& initPos, Matrix<3, 1, double>& initVel, Matrix<3, 1, double>& initAcc, double lambdas, double delta_t) {
		x_states.slice<0, 0, 3, 1>() = initPos;
		x_states.slice<3, 0, 3, 1>() = initVel;
		x_states.slice<6, 0, 3, 1>() = initAcc;
		x_states(9, 0) = lambdas;
		mlambda = lambdas;
		mdelta_t = delta_t;
		Ft = Identity;
		Ft.slice<0, 3, 3, 3 >() = I3 * mdelta_t / mlambda;
		Ft.slice<0, 6, 3, 3 >() = I3 * mdelta_t * mdelta_t / mlambda;
		Ft.slice<3, 6, 3, 3 >() = I3 * mdelta_t;
		P = Identity;
		P.slice<0, 0, 3, 3>() = I3 * slam_pos_noise;
		P.slice<3, 3, 3, 3>() = I3 * vel_noise;
		P.slice<6, 6, 3, 3>() = I3 * acc_noise;
		P(9,9) = lambda_noise;

}
void Kalman::EKFTransSlam(Matrix<3, 1, double>& slam_measurement, double delta_t){

		mlambda = x_states(9, 0);
		Ft.slice<0, 3, 3, 3 >() = I3 * delta_t / mlambda;
		Ft.slice<0, 6, 3, 3 >() = I3 * delta_t * delta_t / mlambda;
		Ft.slice<3, 6, 3, 3 >() = I3 * delta_t;
		Matrix<10, 10> P_temp;

		P_temp = Ft*P*Ft.T();

		x_states = Ft*x_states + v_noise;
		Kv = P_temp*Hv.T()*(Hv*P_temp*Hv.T() + Rslam).T();
		x_states = x_states + Kv*(slam_measurement - x_states.slice<0, 0, 3, 1>());
		P = (I - Kv*Hv)*P_temp;
}
void Kalman::EKFTransSlam(Vector<3, double>& slam_measurement){

		mlambda = x_states(9, 0);
		Ft.slice<0, 3, 3, 3 >() = I3 * mdelta_t / mlambda;
		Ft.slice<0, 6, 3, 3 >() = I3 * mdelta_t * mdelta_t / mlambda;
		Ft.slice<3, 6, 3, 3 >() = I3 * mdelta_t;
		Matrix<10, 10> P_temp;

		P_temp = Ft*P*Ft.T();

		x_states = Ft*x_states + v_noise;
		Kv = P_temp*Hv.T()*(Hv*P_temp*Hv.T() + Rslam).T();
		x_states = x_states + Kv*(slam_measurement.as_col() - x_states.slice<0, 0, 3, 1>());
		P = (I - Kv*Hv)*P_temp;
}
void Kalman::EKFTransIMU(Matrix<3, 1, double>& IMU_measurement, double delta_t){

		mlambda = x_states(9, 0);
		Ft.slice<0, 3, 3, 3 >() = I3 * delta_t / mlambda;
		Ft.slice<0, 6, 3, 3 >() = I3 * delta_t * delta_t / mlambda;
		Ft.slice<3, 6, 3, 3 >() = I3 * delta_t;
		Matrix<10, 10> P_temp;

		P_temp = Ft*P*Ft.T();

		x_states = Ft*x_states + v_noise;
		Ki = P_temp*Hi.T()*(Hi*P_temp*Hi.T() + Rimu).T();
		x_states = x_states + Ki*(IMU_measurement - x_states.slice<6, 0, 3, 1>());
		P = (I - Ki*Hi)*P_temp;
}
void Kalman::EKFTransIMU(Matrix<3, 1, double>& IMU_measurement){

		mlambda = x_states(9, 0);
		Ft.slice<0, 3, 3, 3 >() = I3 * mdelta_t / mlambda;
		Ft.slice<0, 6, 3, 3 >() = I3 * mdelta_t * mdelta_t / mlambda;
		Ft.slice<3, 6, 3, 3 >() = I3 * mdelta_t;
		Matrix<10, 10> P_temp;

		P_temp = Ft*P*Ft.T();

		x_states = Ft*x_states + v_noise;
		Ki = P_temp*Hi.T()*(Hi*P_temp*Hi.T() + Rimu).T();
		x_states = x_states + Ki*(IMU_measurement - x_states.slice<6, 0, 3, 1>());
		P = (I - Ki*Hi)*P_temp;
}
void Kalman::getData(Matrix<10, 10, double>& mP, Matrix<10, 1,double>& x_stat){
	mP = P;
	x_stat = x_states;
}
void Kalman::getData(Vector<3,double>& x_stat){
	x_stat = makeVector(x_states(0,0),x_states(1,0),x_states(2,0));

}
Kalman::~Kalman() {
	// TODO Auto-generated destructor stub
}

} /* namespace APTAM */
