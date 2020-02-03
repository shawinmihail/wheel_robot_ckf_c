#include "Ckf.h"
#include <Eigen/Cholesky>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>

Ckf::Ckf() :
	_X(Vector9::Zero()),
	_N(9)
{
	// init P
	Vector9 pDiag;
	pDiag << /*r*/ 1e-2f, 1e-2f, 1e-2f, /*v*/ 1e-3f, 1e-3f, 1e-3f, /*qv*/ 1e-3f, 1e-3f, 1e-3f;
	Eigen::Matrix<float, 9, 9> P = pDiag.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, 9, 9>> cdP(P);
	_sqrtP = cdP.matrixL();

	// init Q
	Vector9 qDiag;
	qDiag << /*r*/ 1e-3f, 1e-3f, 1e-3f, /*v*/ 1e-4f, 1e-4f, 1e-4f, /*qv*/ 1e-4f, 1e-4f, 1e-4f;
	Eigen::Matrix<float, 9, 9> Q = qDiag.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, 9, 9>> cdQ(Q);
	_sqrtQ = cdQ.matrixL();
	
    // init R
	Vector6 rDiag;
	rDiag << /*r*/ 1e-1f, 1e-1f, 1e-1f, /*v*/ 1e-2f, 1e-2f, 1e-2f;
	Eigen::Matrix<float, 6, 6> R = rDiag.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, 6, 6>> cdR(R);
	_sqrtR = cdR.matrixL();
}

Eigen::Matrix<float, 9, 18> Ckf::generateQubaturePoints(const Vector9& X, const Eigen::Matrix<float, 9, 9>& sqrtP)
{
	Eigen::Matrix<float, 9, 18> cubaturePointsX;
	for (int i = 0; i < _N; i++)
	{
		cubaturePointsX.col(i) = X + sqrt((float)_N) * sqrtP.col(i);
	}
	for (int i = _N; i < 2*_N; i++)
	{
		cubaturePointsX.col(i) = X - sqrt((float)_N) * sqrtP.col(i - _N);
	}

	return cubaturePointsX;
}

Eigen::Matrix<float, 9, 9> Ckf::calcPredictedSqrtP(const Eigen::Matrix<float, 9, 18>& hi, const Eigen::Matrix<float, 9, 9>& sqrtQ)
{

	// TODO: optimize decomp

	Eigen::Matrix<float, 9, 18 + 9> hi_sqrtQ;
	hi_sqrtQ << hi, sqrtQ;

	Eigen::Matrix<float, 18 + 9, 9> hi_sqrtQ_transposed = hi_sqrtQ.transpose();
	Eigen::HouseholderQR<Eigen::Matrix<float, 18 + 9, 9>> qrHolder;
	qrHolder.compute(hi_sqrtQ_transposed);
	Eigen::Matrix<float, 18 + 9, 9> R = qrHolder.matrixQR().triangularView<Eigen::Upper>();

	Eigen::Matrix<float, 9, 18 + 9> L = R.transpose();
	Eigen::Matrix<float, 9, 9> sqrtP = L.block(0, 0, 9, 9);

	//Eigen::Matrix<float, 9, 9> sqrtP;
	return sqrtP;
}

Eigen::Matrix<float, 6, 6> Ckf::calcPredictedSqrtTheta(const Eigen::Matrix<float, 6, 18>& theta, const Eigen::Matrix<float, 6, 6>& sqrtR)
{

	// TODO: optimize decomp

	Eigen::Matrix<float, 6, 18 + 6> theta_sqrR;
	theta_sqrR << theta, sqrtR;

	Eigen::Matrix<float, 18 + 6, 6> theta_sqrR_transposed = theta_sqrR.transpose();
	Eigen::HouseholderQR<Eigen::Matrix<float, 18 + 6, 6>> qrHolder;
	qrHolder.compute(theta_sqrR_transposed);
	Eigen::Matrix<float, 18 + 6, 6> R = qrHolder.matrixQR().triangularView<Eigen::Upper>();

	Eigen::Matrix<float, 6, 18 + 6> L = R.transpose();
	Eigen::Matrix<float, 6, 6> sqrtS = L.block(0, 0, 6, 6);

	//Eigen::Matrix<float, 6, 6> sqrtS;
	return sqrtS;
}

Eigen::Matrix<float, 9, 9> Ckf::calcSqrtP(const Eigen::Matrix<float, 9, 18 + 6>& arg)
{

	// TODO: optimize decomp

	Eigen::Matrix<float, 18 + 6, 9> arg_transposed = arg.transpose();
	Eigen::HouseholderQR<Eigen::Matrix<float, 18 + 6, 9>> qrHolder;
	qrHolder.compute(arg_transposed);
	Eigen::Matrix<float, 18 + 6, 9> R = qrHolder.matrixQR().triangularView<Eigen::Upper>();

	Eigen::Matrix<float, 9, 18 + 6> L = R.transpose();
	Eigen::Matrix<float, 9, 9> sqrtP = L.block(0, 0, 9, 9);

	//Eigen::Matrix<float, 9, 9> sqrtP;
	return sqrtP;
}

Vector10 Ckf::imuMesBodyDin(const Vector10& Y, const Vector3& a, const Vector3& w)
{
	// Y = [r v q]
	Vector4 qw(0.0f, w[0], w[1], w[2]);
	const Vector4 qDot = 0.5f * quatMultiply(/*q = Y(6-9)*/Y.segment(6, 4), qw);
	Vector10 yDot; /*y_dot = [v, a, q_dot]*/
	yDot << Y[3], Y[4], Y[5], a[0], a[1], a[2], qDot[0], qDot[1], qDot[2], qDot[3];
	return yDot;
}

Vector6 Ckf::stateToMes(const Vector9& X)
{
	return X.segment(0, 6);
}

Vector9 Ckf::predictWithImu(const Vector9& X0, const Vector3& aMes, const Vector3& wMes, double dt)
{
	// X = [r v qv]
	Vector4 q0 = quatVecToQuat(/*qv = X(6-8)*/X0.segment(6, 3));
	Vector3 a0 = quatRotate(q0, aMes);
	float gz = -10.0f;
	a0[2] = a0[2] + gz;

	Vector3 w0 = wMes; // no preprocess now

	// Y = [r v q]
	Vector10 Y0;
	Y0 << X0[0], X0[1], X0[2], X0[3], X0[4], X0[5], q0[0], q0[1], q0[2], q0[3];

	//// Eul int
	//Vector10 yDot = imuMesBodyDin(Y0, a0, w0);
	//Vector10 Y = Y0 + yDot * dt;

	// RK4 int
	Vector10 yDot1 = imuMesBodyDin(Y0, a0, w0);
	Vector10 yDot2 = imuMesBodyDin(Y0 + 0.5f * yDot1 * dt, a0, w0);
	Vector10 yDot3 = imuMesBodyDin(Y0 + 0.5f * yDot2 * dt, a0, w0);
	Vector10 yDot4 = imuMesBodyDin(Y0 + 1.0f * yDot3 * dt, a0, w0);
	Vector10 Y = Y0 + (1.0f/6.0f)*(yDot1 + 2* yDot2 + 2* yDot3 + yDot4)*dt;

	// Y to X
	Vector4 q = Y.segment(6, 4);
	float qN = q.norm();
	float eps = 1e-6f;
	if (qN > eps)
	{
		q = q / qN;
	}
	else
	{
		q = Vector4(1.0f, 0.0f, 0.0f, 0.0f);
	}
	Vector3 qv = quatToQuatVec(q);

	Vector9 X;
	X << Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], qv[0], qv[1], qv[2];
	return X;
}

void Ckf::updateImu(const Vector3& aMes, const Vector3& wMes, float dt)
{	
	// X = [r v qv]

	//state prediction
	Eigen::Matrix<float, 9, 18> cubaturePointsX = generateQubaturePoints(_X, _sqrtP);
	Eigen::Matrix<float, 9, 18> predictedCubaturePointsX;
	for (int i = 0; i < 2 * _N; i++)
	{
		//printf("%d\n", i);
		predictedCubaturePointsX.col(i) = predictWithImu(cubaturePointsX.col(i), aMes, wMes, dt);
	}

	Vector9 cubaturePointsSum(Vector9::Zero());
	for (int i = 0; i < 2 * _N; i++)
	{
		cubaturePointsSum = cubaturePointsSum + predictedCubaturePointsX.col(i);
	}

	Vector9 predictedX = 1.0f / (2.0f * _N) * cubaturePointsSum;


	//std::cout << _X.transpose() << std::endl << std::endl;
	//std::cout << predictedX.transpose() << std::endl << std::endl;
	//std::cout << "---------" << std::endl << std::endl;

	// state vector covariance prediction
	Eigen::Matrix<float, 9, 18> predictedHi;
	for (int i = 0; i < 2 * _N; i++)
	{
		predictedHi.col(i) = predictedCubaturePointsX.col(i) - predictedX;
	}
	predictedHi = predictedHi / sqrt(2.0f * _N);

	_X = predictedX;
	_sqrtP = calcPredictedSqrtP(predictedHi, _sqrtQ);


	//_X = predictWithImu(_X, aMes, wMes, dt);
}

void Ckf::updateGps(const Vector6& pv)
{
	 //X = [r v qv]
	 //Z = [r v]

	// Evaluate propagated cubature points
	Eigen::Matrix<float, 9, 18> predictedCubaturePointsX = generateQubaturePoints(_X, _sqrtP); // _X is last prediction if imu_rate > gps_rate
	Eigen::Matrix<float, 6, 18> predictedCubaturePointsZ;
	for (int i = 0; i < 2 * _N; i++)
	{
		predictedCubaturePointsZ.col(i) = stateToMes(predictedCubaturePointsX.col(i));
	}

	// Estimate the predicted measurement
	Vector6 cubaturePointsSum(Vector6::Zero());
	for (int i = 0; i < 2 * _N; i++)
	{
		cubaturePointsSum = cubaturePointsSum + predictedCubaturePointsZ.col(i);
	}

	Vector6 predictedZ = 1.0f / (2.0f * _N) * cubaturePointsSum;

	// Estimate the square-root of the innovation covariance matrix
	Eigen::Matrix<float, 6, 18> predictedTheta;
	for (int i = 0; i < 2 * _N; i++)
	{
		predictedTheta.col(i) = predictedCubaturePointsZ.col(i) - predictedZ;
	}
	predictedTheta = predictedTheta / sqrt(2.0f * _N);
	Eigen::Matrix<float, 6, 6> sqrtS = calcPredictedSqrtTheta(predictedTheta, _sqrtR);

	// Estimate the cross - covariance matrix
	Eigen::Matrix<float, 9, 18> predictedHi;
	for (int i = 0; i < 2 * _N; i++)
	{
		predictedHi.col(i) = predictedCubaturePointsX.col(i) - _X;
	}
	predictedHi = predictedHi / sqrt(2.0f * _N);

	Eigen::Matrix<float, 9, 6> Pxz = predictedHi * predictedTheta.transpose();

	//Estimate Kalman gain
	Eigen::Matrix<float, 9, 6> K = Pxz * (sqrtS * sqrtS.transpose()).inverse();

	// Estimate the updated state
	_X = _X + K * (pv - predictedZ);
	Eigen::Matrix<float, 9, 18 + 6> triaArg;
	triaArg << predictedHi - K * predictedTheta, K * _sqrtR;
	_sqrtP = calcSqrtP(triaArg);
}

Vector9 Ckf::getEstState()
{
	return _X;
}


