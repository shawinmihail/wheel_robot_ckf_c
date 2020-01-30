#include "QuatMath.h"

Mat quatFromEul(Mat eul)
{
	 Mat quat(4, 1);
	 double roll = eul(0, 0);
	 double pitch = eul(1, 0);
	 double yaw = eul(2, 0);

	 double cy = cos(yaw * 0.5);
	 double sy = sin(yaw * 0.5);
	 double cr = cos(roll * 0.5);
	 double sr = sin(roll * 0.5);
	 double cp = cos(pitch * 0.5);
	 double sp = sin(pitch * 0.5);

	 double q0 = cy * cr * cp + sy * sr * sp;
	 double q1 = cy * sr * cp - sy * cr * sp;
	 double q2 = cy * cr * sp + sy * sr * cp;
	 double q3 = sy * cr * cp - cy * sr * sp;

	 quat(0, 0) = q0;
	 quat(1, 0) = q1;
	 quat(2, 0) = q2;
	 quat(3, 0) = q3;

	 return quat;
}

Mat quatToQuatVec(Mat q) // assume quat scalar part quat[0] > 0;
{
	Mat qv(3,1);
	qv(0, 0) = q(1, 0);
	qv(1, 0) = q(2, 0);
	qv(2, 0) = q(3, 0);

	if (q(0, 0) < 0)
	{
		double sinHalfAlpha = qv.getNorm();
		double eps = 1e-8;
		if (sinHalfAlpha < eps) {
			qv.Zero();
			return qv;
		};
		double scaleFactor = 1.0 / sinHalfAlpha;
		qv.ScaleMat2(&scaleFactor); // pin
		double alpha = 2 * asin(sinHalfAlpha);
		double pi = 3.1415; // use WGS4 PI here;
		double alphaNew = -2 * pi + alpha; // rotate to another dir
		double revertScale = -1.0;
		qv.ScaleMat2(&revertScale); // pin = - pin;
		double sinHalfNewAlpha = sin(alphaNew / 2.0);
		qv.ScaleMat2(&sinHalfNewAlpha); // qv = -pin * sin(alpha_new / 2)
	}
	return qv;
}

Mat quatVecToQuat(Mat qv) // assume quat scalar part quat[0] > 0;
{
	double sq_q0 = 1 - qv(0, 0) * qv(0, 0) - qv(1, 0) * qv(1, 0) - qv(2, 0) * qv(2, 0);
	if (sq_q0 < 0.0) // possible in case of numerical integration error
	{
		sq_q0 = 0.0;
	}
	double q0 = sqrt(sq_q0);

	Mat q(4, 1);
	q(0, 0) = q0;
	q(1, 0) = qv(0, 0);
	q(2, 0) = qv(1, 0);
	q(3, 0) = qv(2, 0);
	return q;
}

Mat quatMultiply(Mat q, Mat r)
{
	Mat p(4, 1);
	p(0, 0) = r(0, 0) * q(0, 0) - r(1, 0) * q(1, 0) - r(2, 0) * q(2, 0) - r(3, 0) * q(3, 0);
	p(1, 0) = r(0, 0) * q(1, 0) + r(1, 0) * q(0, 0) - r(2, 0) * q(3, 0) + r(3, 0) * q(2, 0);
	p(2, 0) = r(0, 0) * q(2, 0) + r(1, 0) * q(3, 0) + r(2, 0) * q(0, 0) - r(3, 0) * q(1, 0);
	p(3, 0) = r(0, 0) * q(3, 0) - r(1, 0) * q(2, 0) + r(2, 0) * q(1, 0) + r(3, 0) * q(0, 0);
	return p;
}

Mat quatInverse(Mat q)
{
	q(1, 0) *= -1.0;
	q(2, 0) *= -1.0;
	q(3, 0) *= -1.0;
	return q;
}

Mat quatRotate(Mat q, Mat v)
{
	Mat qv(4, 1);
	qv(0, 0) = 0.0;
	qv(1, 0) = v(0, 0);
	qv(2, 0) = v(1, 0);
	qv(3, 0) = v(2, 0);

	Mat qDual = quatInverse(q);
	Mat qv1 = quatMultiply(qv, qDual);
	Mat qv2 = quatMultiply(q, qv1);

	v(0, 0) = qv2(1, 0);
	v(1, 0) = qv2(2, 0);
	v(2, 0) = qv2(3, 0);

	return v;
}


