/*
 * Kalman.cc
 *
 *
 *  Created on: Jun 30, 2016
 *      Author: lyw
 */
#include <iostream>
#include <Kalman.h>
#include <fstream>
using namespace Eigen;

Kalman::Kalman(){
		Hv.setZero();
		Hv.block(0, 0, 3, 3) = I3;
		Hi.setZero();
		Hi.block(0, 6, 3, 3) = I3;
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
		Rimu *= 0.000001;
		Rslam *= 1;
}
Kalman::Kalman(Matrix<double, 3, 1>& initPos, Matrix<double, 3, 1>& initVel, Matrix<double, 3, 1>& initAcc, double lambdas, double delta_t) {
		Hv.setZero();
		Hv << 1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0;
		Hi.setZero();
		Hi << 0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0;
        v_noise << 0.5,0.5,0.5,0.1,0.1,0.1,0.0005,0.0005,0.0005, 0.0041;
        v_noise *= 0.001;
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
		Ft(0,9) = -(mdelta_t + mdelta_t*mdelta_t/2)/(mlambda*mlambda);

		P.setIdentity();
		P.block(0, 0, 3, 3) = I3 * slam_pos_noise*1;
		P.block(3, 3, 3, 3) = I3 * vel_noise*1;
		P.block(6, 6, 3, 3) = I3 * acc_noise*1;
		P(9,9) = lambda_noise*1;
		Q.setIdentity();
		Q.block(0, 0, 3, 3) = I3 * qslam_pos_noise*0.1;//0.1-----0.1
		Q.block(3, 3, 3, 3) = I3 * qvel_noise*0.004;//0.0005 0.00034-----0.004
		Q.block(6, 6, 3, 3) = I3 * qacc_noise*0.06;//0.002-----0.06
		Q(9,9) = qlambda_noise*0.0012;//0.0012-----0.0012

		P *= 0.005;//0.005----0.005
		Q *= 0.01;//0.01-----0.01
		Rimu.setIdentity();
		Rslam.setIdentity();
		Rimu *= 0.00008;//0.0000001---0.000008
		Rslam *= 0.0004;//0.002------0.004


}
void Kalman::EKFTransSlam(Matrix< double, 3, 1>& slam_measurement, double delta_t){

		mlambda = x_states(9, 0);
		Ft.block(0, 3, 3, 3 ) = I3 * delta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * delta_t * delta_t / (2*mlambda);
		Ft.block(3, 6, 3, 3 ) = I3 * delta_t;
        Ft(0,9) = -(mdelta_t + mdelta_t*mdelta_t/2)/(mlambda*mlambda);



		Matrix<double, 10, 10> P_temp;
        Matrix<double, 10, 1> x_temp, x_gain;

		P_temp = Ft*P*Ft.transpose() + Q*mlambda;
        P_temp = (P_temp + P_temp.transpose()) / 2;

		x_temp = Ft*x_states + v_noise*(rand()%100 - 50)/200.0;

		Kv = P_temp*Hv.transpose()*((Hv*P_temp*Hv.transpose() + Rslam).inverse());
		x_gain = Kv*(slam_measurement - x_temp.block(0, 0, 3, 1));
		x_states = x_temp + x_gain;
		//P = (I - Kv*Hv)*P_temp;
		P = (I - Kv*Hv)*P_temp*((I-Kv*Hv).transpose()) + Kv*Rslam*Kv.transpose();
		P_temp = P.transpose();
		P = (P + P_temp) / 2;
}
void Kalman::EKFTransSlam(Vector3d& slam_measurement){

		mlambda = x_states(9, 0);
		Ft.block(0, 3, 3, 3 ) = I3 * mdelta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * mdelta_t * mdelta_t / (2*mlambda);
		Ft.block(3, 6, 3, 3 ) = I3 * mdelta_t;

		Matrix<double, 10, 10> P_temp;
        Matrix<double, 10, 1> x_temp;

		P_temp = Ft*P*Ft.transpose() + Q;
        P_temp = (P_temp + P_temp.transpose()) / 2;

		x_temp = Ft*x_states + v_noise*(rand()%100 - 50)/200.0;

		Kv = P_temp*Hv.transpose()*((Hv*P_temp*Hv.transpose() + Rslam).inverse());
		x_states = x_temp + Kv*(slam_measurement - x_temp.block(0, 0, 3, 1));
		//P = (I - Kv*Hv)*P_temp;
		P = (I - Kv*Hv)*P_temp*((I-Kv*Hv).transpose()) + Kv*Rslam*Kv.transpose();
		P_temp = P.transpose();
		P = (P + P_temp) / 2;

}
void Kalman::EKFTransIMU(Matrix< double, 3, 1>& acc_measurement, double delta_t){

		mlambda = x_states(9, 0);
		Ft.block(0, 3, 3, 3 ) = I3 * delta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * delta_t * delta_t / (2*mlambda);
		Ft.block(3, 6, 3, 3 ) = I3 * delta_t;
        Ft(0,9) = -(mdelta_t + mdelta_t*mdelta_t/2)/(mlambda*mlambda);
		Matrix<double, 10, 10> P_temp;
        Matrix<double, 10, 1> x_temp, x_gain;

		P_temp = Ft*P*Ft.transpose() + Q*mlambda;
        P_temp = (P_temp + P_temp.transpose()) / 2;

		x_temp = Ft*x_states + v_noise*(rand()%100 - 50)/200.0;

		//Ki = P_temp*Hi.transpose()*((Hi*P_temp*Hi.transpose() + Rimu).inverse());
		Ki = P_temp.block(0, 6, 10, 3)*((P_temp.block(6,6, 3,3) + Rimu).inverse());
		x_gain = Ki*(acc_measurement - x_temp.block(6, 0, 3, 1));
		x_states = x_temp + x_gain;
		//P = (I - Ki*Hi)*P_temp;
		P = (I - Ki*Hi)*P_temp*((I-Ki*Hi).transpose()) + Ki*Rimu*Ki.transpose();
		P_temp = P.transpose();
		P = (P + P.transpose()) / 2;
		if (mlambda < 1.5){
            int d;
            d ++;
		}

}
void Kalman::EKFTransIMU(Matrix< double, 3, 1>& acc_measurement){

		mlambda = x_states(9, 0);

		Ft.block(0, 3, 3, 3 ) = I3 * mdelta_t / mlambda;
		Ft.block(0, 6, 3, 3 ) = I3 * mdelta_t * mdelta_t / (2*mlambda);
		Ft.block(3, 6, 3, 3 ) = I3 * mdelta_t;
        Ft(0,9) = -(mdelta_t + mdelta_t*mdelta_t/2)/(mlambda*mlambda);
    	/*
        Ft(0,3) = mdelta_t / mlambda;
        Ft(1,4) = mdelta_t / mlambda;
        Ft(2,5) = mdelta_t / mlambda;
        Ft(0,3) = mdelta_t / mlambda;
        Ft(0,6) = mdelta_t * mdelta_t / (2*mlambda);
        Ft(1,7) = mdelta_t * mdelta_t / (2*mlambda);
        Ft(2,8) = mdelta_t * mdelta_t / (2*mlambda);
        Ft(3,6) = mdelta_t;
        Ft(4,7) = mdelta_t;
        Ft(5,8) = mdelta_t;
 */
		Matrix<double, 10, 10> P_temp;
        Matrix<double, 10, 1> x_temp, x_gain;

		P_temp = Ft*P*Ft.transpose() + Q*mlambda;
        P_temp = (P_temp + P_temp.transpose()) / 2;

		x_temp = Ft*x_states + v_noise*(rand()%100 - 50)/200.0;

		//Ki = P_temp*Hi.transpose()*((Hi*P_temp*Hi.transpose() + Rimu).inverse());
		Ki = P_temp.block(0, 6, 10, 3)*((P_temp.block(6,6, 3,3) + Rimu).inverse());
		x_gain = Ki*(acc_measurement - x_temp.block(6, 0, 3, 1));
		x_states = x_temp + x_gain;
		P = (I - Ki*Hi)*P_temp*((I-Ki*Hi).transpose()) + Ki*Rimu*Ki.transpose();
		P_temp = P.transpose();
		P = (P + P.transpose()) / 2;

}
void Kalman::getData(Matrix< double, 10, 10>& mP, Matrix< double, 10, 1>& x_stat){
	mP = P;
	x_stat = x_states;
}
void Kalman::getData(Vector3d& x_stat, double& scale_out){
	x_stat << x_states(0,0),x_states(1,0),x_states(2,0);
	std::cout << mlambda << std::endl;
    scale_out = mlambda;
	///std::cout << (rand()%100 - 50)/400.0 << std::endl;
	//x_stat *= mlambda;

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

void Kalman::quat2cbn(Matrix<double, 4, 1>& quaternion, Matrix3d& dcm){

    dcm << (quaternion(0)*quaternion(0)+quaternion(1)*quaternion(1)-quaternion(2)*quaternion(2)-quaternion(3)*quaternion(3)),
            (2*(quaternion(1)*quaternion(2)-quaternion(0)*quaternion(3))),
            (2*(quaternion(1)*quaternion(3)+quaternion(0)*quaternion(2))),
            (2*(quaternion(1)*quaternion(2)+quaternion(0)*quaternion(3))),
            (quaternion(0)*quaternion(0)-quaternion(1)*quaternion(1)+quaternion(2)*quaternion(2)-quaternion(3)*quaternion(3)),
            (2*(quaternion(2)*quaternion(3)-quaternion(0)*quaternion(1))),
            (2*(quaternion(1)*quaternion(3)-quaternion(0)*quaternion(2))),
            (2*(quaternion(2)*quaternion(3)+quaternion(0)*quaternion(1))),
            (quaternion(0)*quaternion(0)-quaternion(1)*quaternion(1)-quaternion(2)*quaternion(2)+quaternion(3)*quaternion(3));

    double normq=quaternion(0)*quaternion(0)+quaternion(1)*quaternion(1)+quaternion(2)*quaternion(2)+quaternion(3)*quaternion(3);
    if (normq != 0) dcm /= sqrt(normq);
    else dcm = I3;

}


Kalman::~Kalman() {

	// TODO Auto-generated destructor stub
}

 /* namespace APTAM */
