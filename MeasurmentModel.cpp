//#include "MeasurmentModel.h"
//
//MeasurmentModel::MeasurmentModel() :
//_mesState(Mat(15, 1)),
//
//_generator(200),
//_imuAccDist(0.0, 0.05),
//_imuRotVelDist(0.0, 0.05),
//_gpsPosDist(0.0, 0.5),
//_gpsVelDist(0.0, 0.05)
//{
//	_mesState.Zero();
//}
//
//void MeasurmentModel::setData(Mat fullState /* r v a qv w*/)
//{
//	// pose vel in local inertial frame
//
//	_mesState(0, 0) = fullState(0, 0) + _gpsPosDist(_generator);
//	_mesState(1, 0) = fullState(1, 0) + _gpsPosDist(_generator);
//	_mesState(2, 0) = fullState(2, 0) + _gpsPosDist(_generator);
//
//	_mesState(3, 0) = fullState(3, 0) + _gpsVelDist(_generator);
//	_mesState(4, 0) = fullState(4, 0) + _gpsVelDist(_generator);
//	_mesState(5, 0) = fullState(5, 0) + _gpsVelDist(_generator);
//
//	// attitude
//	Mat qv(3, 1);
//	qv(0, 0) = fullState(9, 0);
//	qv(1, 0) = fullState(10, 0);
//	qv(2, 0) = fullState(11, 0);
//	Mat q = quatVecToQuat(qv);
//	_mesState(9, 0) = fullState(9, 0);
//	_mesState(10, 0) = fullState(10, 0);
//	_mesState(11, 0) = fullState(11, 0);
//
//	// imu in body frame
//	double gz = -10;
//	Mat a(3, 1);
//	_mesState(6, 0) = fullState(6, 0) + _imuAccDist(_generator);
//	_mesState(7, 0) = fullState(7, 0) + _imuAccDist(_generator);
//	_mesState(8, 0) = fullState(8, 0) + _imuAccDist(_generator) - gz;
//	Mat quatDual = quatInverse(q);
//	Mat aB = quatRotate(quatDual, a);
//
//	// suppose rot rate the same in I and B frames (only Z rotation)
//	Mat w(3, 1);
//	_mesState(12, 0) = fullState(12, 0) + _imuRotVelDist(_generator);
//	_mesState(13, 0) = fullState(13, 0) + _imuRotVelDist(_generator);
//	_mesState(14, 0) = fullState(14, 0) + _imuRotVelDist(_generator);
//
//}
//
//Mat MeasurmentModel::getMesState()
//{
//	return _mesState;
//}