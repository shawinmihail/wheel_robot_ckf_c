#include "MeasurmentModel.h"

MeasurmentModel::MeasurmentModel() :
_a(Mat(3, 1)),
_w(Mat(3, 1)),
_r(Mat(3, 1)),
_v(Mat(3, 1)),

_generator(200),
_imuAccDist(0.0, 0.05),
_imuRotVelDist(0.0, 0.05),
_gpsPosDist(0.0, 0.5),
_gpsVelDist(0.0, 0.05)
{

}

void MeasurmentModel::setData(Mat fullState /* r v a qv w*/)
{
	// pose vel in local inertial frame
	Mat r(3, 1);
	r(0, 0) = fullState(0, 0) + _gpsPosDist(_generator);
	r(1, 0) = fullState(1, 0) + _gpsPosDist(_generator);
	r(2, 0) = fullState(2, 0) + _gpsPosDist(_generator);
	_r.Move(r);

	Mat v(3, 1);
	v(0, 0) = fullState(3, 0) + _gpsVelDist(_generator);
	v(1, 0) = fullState(4, 0) + _gpsVelDist(_generator);
	v(2, 0) = fullState(5, 0) + _gpsVelDist(_generator);
	_v.Move(v);

	// attitude
	Mat qv(3, 1);
	qv(0, 0) = fullState(9, 0);
	qv(1, 0) = fullState(10, 0);
	qv(2, 0) = fullState(11, 0);
	Mat q = quatVecToQuat(qv);

	// imu in body frame
	double gz = -10;
	Mat a(3, 1);
	a(0, 0) = fullState(6, 0) + _imuAccDist(_generator);
	a(1, 0) = fullState(7, 0) + _imuAccDist(_generator);
	a(2, 0) = fullState(8, 0) + _imuAccDist(_generator) - gz;
	Mat quatDual = quatInverse(q);
	Mat aB = quatRotate(quatDual, a);
	_a.Move(aB);

	// suppose rot rate the same in I and B frames (only Z rotation)
	Mat w(3, 1);
	w(0, 0) = fullState(12, 0) + _imuRotVelDist(_generator);
	w(1, 0) = fullState(13, 0) + _imuRotVelDist(_generator);
	w(2, 0) = fullState(14, 0) + _imuRotVelDist(_generator);
	_w.Move(w);
}