#include "Ckf.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>

const float CKF_EPS = 1e-6f;

// returns avg rotation
// uses https://stackoverflow.com/questions/12374087/average-of-multiple-quaternions,
// http://www.acsu.buffalo.edu/~johnc/ave_quat07.pdf:
// Markley F. L. et al. Averaging quaternions
// Journal of Guidance, Control, and Dynamics. – 2007. – Ò. 30. – ¹. 4. – Ñ. 1193-1197.
// Q = [q1, ..., qN];
template <int n>
Vector4 avgQuat(const Eigen::Matrix<float, 4, n>& Q)
{
	Eigen::Matrix<float, 4, 4> QQT = Q * Q.transpose();
	Eigen::SelfAdjointEigenSolver <Eigen::Matrix<float, 4, 4>> solver(QQT);
	Vector4 eigenvalues = solver.eigenvalues();

	int k = 0;
	for (int i = 1; i < 4; i++)
	{
		if (eigenvalues[i] > eigenvalues[k])
		{
			k = i;
		}
	}
	return solver.eigenvectors().col(k);
}

/* returns L from decomp A = L * Q, where L is low triang matrix, Q is ortogonal matrix */
// TODO: optimize decomp
template <int a, int b>
Eigen::Matrix<float, a, a> getLowTriang(const Eigen::Matrix<float, a, b>& A)
{ 
	Eigen::Matrix<float, b, a> A_transposed = A.transpose();
	Eigen::HouseholderQR<Eigen::Matrix<float, b, a>> qrHolder;
	qrHolder.compute(A_transposed);
	Eigen::Matrix<float, b, a> R = qrHolder.matrixQR().triangularView<Eigen::Upper>();
	Eigen::Matrix<float, a, b> L = R.transpose();
	Eigen::Matrix<float, a, a> cutedL = L.block(0, 0, a, a);

	return -cutedL;
}

Ckf::Ckf() :
	_X(CkfStateVector::Zero())
{
	// init Q
	CkfStateVector qDiag;
	qDiag << /*r*/ 1e-3f, 1e-3f, 1e-3f, /*v*/ 1e-4f, 1e-4f, 1e-4f, /*qv*/ 1e-8f, 1e-8f, 1e-8f;
	Eigen::Matrix<float, CKF_STATE_DIM, CKF_STATE_DIM> Q = qDiag.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, CKF_STATE_DIM, CKF_STATE_DIM>> cdQ(Q);
	_sqrtQ = cdQ.matrixL();

	// init R
	CkfMesVector rDiag;
	rDiag << /*r*/ 1e-1f, 1e-1f, 1e-1f, /*v*/ 1e-2f, 1e-2f, 1e-2f;
	Eigen::Matrix<float, CKF_MES_DIM, CKF_MES_DIM> R = rDiag.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, CKF_MES_DIM, CKF_MES_DIM>> cdR(R);
	_sqrtR = cdR.matrixL();

	// init P
	CkfStateVector pDiag = 33 * qDiag;
	//pDiag << /*r*/ 1e-2f, 1e-2f, 1e-2f, /*v*/ 1e-3f, 1e-3f, 1e-3f, /*qv*/ 1e-3f, 1e-3f, 1e-3f;
	Eigen::Matrix<float, CKF_STATE_DIM, CKF_STATE_DIM> P = pDiag.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, CKF_STATE_DIM, CKF_STATE_DIM>> cdP(P);
	_sqrtP = cdP.matrixL();
}

Eigen::Matrix<float, CKF_STATE_DIM, 2 * CKF_STATE_DIM> Ckf::generateQubaturePoints(const CkfStateVector& X, const Eigen::Matrix<float, CKF_STATE_DIM, CKF_STATE_DIM>& sqrtP)
{
	Eigen::Matrix<float, CKF_STATE_DIM, 2 * CKF_STATE_DIM> cubaturePointsX;
	for (int i = 0; i < CKF_STATE_DIM; i++)
	{
		cubaturePointsX.col(i) = X + sqrt((float)CKF_STATE_DIM) * sqrtP.col(i);
	}
	for (int i = CKF_STATE_DIM; i < 2 * CKF_STATE_DIM; i++)
	{
		cubaturePointsX.col(i) = X - sqrt((float)CKF_STATE_DIM) * sqrtP.col(i - CKF_STATE_DIM);
	}

	return cubaturePointsX;
}

Vector10 Ckf::imuMesBodyDin(const Vector10& Y, const Vector3& a, const Vector3& w)
{
	// Y = [r v q]
	Vector4 qw(0.0f, w[0], w[1], w[2]);
	Vector4 qDot = 0.5f * quatMultiply(/*q = Y(6-9)*/Y.segment(6, 4), qw);
	Vector10 yDot; /*y_dot = [v, a, q_dot]*/
	yDot << Y[3], Y[4], Y[5], a[0], a[1], a[2], qDot[0], qDot[1], qDot[2], qDot[3];
	return yDot;
}

CkfMesVector Ckf::stateToMes(const CkfStateVector& X)
{
	return X.segment(0, 6);
}

CkfStateVector Ckf::predictWithImu(const CkfStateVector& X0, const Vector3& aMes, const Vector3& wMes, float dt)
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

	// Eul int
	Vector10 yDot = imuMesBodyDin(Y0, a0, w0);
	Vector10 Y = Y0 + yDot * dt;
	if (Y[9] > 1.0f)
	{
		double x = 0;
	}

	//// RK4 int
	//Vector10 yDot1 = imuMesBodyDin(Y0, a0, w0);
	//Vector10 yDot2 = imuMesBodyDin(Y0 + 0.5f * yDot1 * dt, a0, w0);
	//Vector10 yDot3 = imuMesBodyDin(Y0 + 0.5f * yDot2 * dt, a0, w0);
	//Vector10 yDot4 = imuMesBodyDin(Y0 + 1.0f * yDot3 * dt, a0, w0);
	//Vector10 Y = Y0 + (1.0f / 6.0f) * (yDot1 + 2 * yDot2 + 2 * yDot3 + yDot4) * dt;

	// Y to X
	Vector4 q = Y.segment(6, 4);
	float qN = q.norm();
	q = q / qN;

	Vector3 qv = quatToQuatVec(q);

	CkfStateVector X;
	X << Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], qv[0], qv[1], qv[2];
	return X;
}

void Ckf::updateImu(const Vector3& aMes, const Vector3& wMes, float dt)
{
	//X = [r v qv]

	//1. state prediction
	// a. cub points
	// simple addition for r, v
	Eigen::Matrix<float, CKF_STATE_DIM, 2 * CKF_STATE_DIM> cubaturePointsX;
	for (int i = 0; i < CKF_STATE_DIM; i++)
	{
		cubaturePointsX.col(i).segment(0,6) = _X.segment(0, 6) + sqrt((float)CKF_STATE_DIM) * _sqrtP.col(i).segment(0, 6);
	}
	for (int i = CKF_STATE_DIM; i < 2 * CKF_STATE_DIM; i++)
	{
		cubaturePointsX.col(i).segment(0, 6) = _X.segment(0, 6) - sqrt((float)CKF_STATE_DIM) * _sqrtP.col(i - CKF_STATE_DIM).segment(0, 6);
	}
	// quaternion addition for qv
	for (int i = 0; i < CKF_STATE_DIM; i++)
	{
		Vector3 qvAddition = sqrt((float)CKF_STATE_DIM) * _sqrtP.col(i).segment(6, 3);
		Vector4 qAddition = quatVecToQuat(qvAddition);
		Vector3 qv = _X.segment(6, 3);
		Vector4 q = quatVecToQuat(qv);
		cubaturePointsX.col(i).segment(6, 3) = quatToQuatVec(quatMultiply(q, qAddition));
	}
	for (int i = CKF_STATE_DIM; i < 2 * CKF_STATE_DIM; i++)
	{
		Vector3 qvAddition = -sqrt((float)CKF_STATE_DIM) * _sqrtP.col(i - CKF_STATE_DIM).segment(6, 3);
		Vector4 qAddition = quatVecToQuat(qvAddition);
		Vector3 qv = _X.segment(6, 3);
		Vector4 q = quatVecToQuat(qv);
		cubaturePointsX.col(i).segment(6, 3) = quatToQuatVec(quatMultiply(q, qAddition));
	}

	// b. evaluate cub points
	Eigen::Matrix<float, CKF_STATE_DIM, 2 * CKF_STATE_DIM> predictedCubaturePointsX;
	for (int i = 0; i < 2 * CKF_STATE_DIM; i++)
	{
		CkfStateVector predictedState = predictWithImu(cubaturePointsX.col(i), aMes, wMes, dt);
		predictedCubaturePointsX.col(i) = predictedState;
	}

	// c. find state expectation as arithmetic mean of evaluated cub points
	// r v as simple mean
	CkfStateVector predictedX(CkfStateVector::Zero());
	for (int i = 0; i < 2 * CKF_STATE_DIM; i++)
	{
		predictedX.segment(0, 6) = predictedX.segment(0, 6) + predictedCubaturePointsX.col(i).segment(0, 6) / (2.0f * CKF_STATE_DIM);
	}
	// qv as quaternion rotation mean
	Eigen::Matrix<float, 4, 2 * CKF_STATE_DIM> Q;
	for (int i = 0; i < 2 * CKF_STATE_DIM; i++)
	{
		Vector4 qCub = quatVecToQuat(predictedCubaturePointsX.col(i).segment(6, 3));
		Q.col(i) = qCub;
	}
	// find avg rotation
	Vector4 qAvg = avgQuat<2 * CKF_STATE_DIM>(Q);
	predictedX.segment(6, 3) = quatToQuatVec(qAvg);

	// 2. state vector covariance prediction
	// a. find hi
	Eigen::Matrix<float, CKF_STATE_DIM, 2 * CKF_STATE_DIM> predictedHi;
	for (int i = 0; i < 2 * CKF_STATE_DIM; i++)
	{
		// rv simple
		predictedHi.col(i).segment(0, 6) = predictedCubaturePointsX.col(i).segment(0, 6) - predictedX.segment(0, 6);
		//qv quat reduction
		Vector3 qv = predictedCubaturePointsX.col(i).segment(6, 3);
		Vector4 q = quatVecToQuat(qv);
		Vector3 dqv = predictedX.segment(6, 3);
		Vector4 dq = quatVecToQuat(dqv);
		predictedHi.col(i).segment(6, 3) = quatToQuatVec(quatMultiply(q, quatInverse(dq)));
	}
	predictedHi = predictedHi / sqrt(2.0f * CKF_STATE_DIM);

	// b. sqrtP = tria(hi | sqrtQ)
	Eigen::Matrix<float, CKF_STATE_DIM, 3 * CKF_STATE_DIM> hi_sqrtQ;
	hi_sqrtQ << predictedHi, _sqrtQ;
	
	//3. save calcs for correction call
	_sqrtP = getLowTriang<CKF_STATE_DIM, 3 * CKF_STATE_DIM>(hi_sqrtQ);
	_X = predictedX;
	_predictedHi = predictedHi;
	_predictedCubaturePointsX = predictedCubaturePointsX;

	//_X = predictWithImu(_X, aMes, wMes, dt); // for tests
}

void Ckf::updateGps(const CkfMesVector& pv)
{
	// correction works appropriate if prediction was called before
	//X = [r v qv]
	//Z = [r v]
	
    // 1. find expected meas
	// a. translate cub points into meas space
	Eigen::Matrix<float, CKF_MES_DIM, 2 * CKF_STATE_DIM> predictedCubaturePointsZ;
	for (int i = 0; i < 2 * CKF_STATE_DIM; i++)
	{
		predictedCubaturePointsZ.col(i) = stateToMes(_predictedCubaturePointsX.col(i));
	}

	// c. find mes expectation as arithmetic mean of cub points
	// r v as simple mean
	CkfMesVector predictedZ(CkfMesVector::Zero());
	for (int i = 0; i < 2 * CKF_STATE_DIM; i++)
	{
		predictedZ = predictedZ + predictedCubaturePointsZ.col(i) / (2.0f * CKF_STATE_DIM);
	}
	
	// 2. innovation covariance matrix
	Eigen::Matrix<float, CKF_MES_DIM, 2 * CKF_STATE_DIM> predictedTheta;
	for (int i = 0; i < 2 * CKF_STATE_DIM; i++)
	{
		predictedTheta.col(i) = predictedCubaturePointsZ.col(i) - predictedZ;
	}
	predictedTheta = predictedTheta / sqrt(2.0f * CKF_STATE_DIM);

	Eigen::Matrix<float, CKF_MES_DIM, CKF_MES_DIM + 2 * CKF_STATE_DIM> tehta_sqrtQ;
	tehta_sqrtQ << predictedTheta, _sqrtR;
	Eigen::Matrix<float, CKF_MES_DIM, CKF_MES_DIM> sqrtS = getLowTriang<CKF_MES_DIM, CKF_MES_DIM + 2 * CKF_STATE_DIM>(tehta_sqrtQ);

	// 3. Cross cov matrix
	Eigen::Matrix<float, CKF_STATE_DIM, CKF_MES_DIM> Pxz = _predictedHi * predictedTheta.transpose();

	//4. Estimate Kalman gain
	Eigen::Matrix<float, CKF_STATE_DIM, CKF_MES_DIM> K = Pxz * (sqrtS * sqrtS.transpose()).inverse();

	// 5. Estimate the updated state
	_X = _X + 1 * K * (pv - predictedZ);
	Eigen::Matrix<float, CKF_STATE_DIM, 2 * CKF_STATE_DIM + CKF_MES_DIM> triaArg;
	triaArg << _predictedHi - K * predictedTheta, K * _sqrtR;
	_sqrtP = getLowTriang<CKF_STATE_DIM, 2 * CKF_STATE_DIM + CKF_MES_DIM>(triaArg);

	//std::cout << _sqrtP* _sqrtP.transpose() << "\n\n";
	//std::cout << _sqrtP << "\n\n";
}

CkfStateVector Ckf::getEstState()
{
	return _X;
}

