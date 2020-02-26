#pragma once
#include "Definitions.h"
#include "Utils.h"

// state X = [r , v, qv]
// mes Z = [r, v]
#define CKF_STATE_DIM 9
#define CKF_MES_DIM 6
typedef Eigen::Matrix<float, CKF_STATE_DIM, 1> CkfStateVector;
typedef Eigen::Matrix<float, CKF_MES_DIM, 1> CkfMesVector;


class Ckf
{
public:
	Ckf();
	void updateImu(const Vector3& aMes, const Vector3& wMes, float dt);
	void updateGps(const CkfMesVector& pv);
	CkfStateVector getEstState();
private:

	Vector10 imuMesBodyDin(const Vector10& Y, const Vector3& a, const Vector3& w);
	CkfStateVector predictWithImu(const CkfStateVector& X, const Vector3& aMes, const Vector3& wMes, float dt);
	CkfMesVector stateToMes(const CkfStateVector& X);
	Eigen::Matrix<float, CKF_STATE_DIM, 2 * CKF_STATE_DIM > generateQubaturePoints(const CkfStateVector& X, const Eigen::Matrix<float, 9, 9>& sqrtP);

	CkfStateVector _X; // r v qv
	Eigen::Matrix<float, CKF_STATE_DIM, CKF_STATE_DIM> _sqrtQ;
	Eigen::Matrix<float, CKF_MES_DIM, CKF_MES_DIM> _sqrtR;
	Eigen::Matrix<float, CKF_STATE_DIM, CKF_STATE_DIM> _sqrtP;

	Eigen::Matrix<float, CKF_STATE_DIM, 2 * CKF_STATE_DIM> _predictedHi;
	Eigen::Matrix<float, CKF_STATE_DIM, 2 * CKF_STATE_DIM> _predictedCubaturePointsX;
};

