#include "Utils.h"
#include <math.h>

void printVect3(const Vector3& vect)
{
	for (int i = 0; i < 3; i++)
	{
		printf("%.6f, ", vect[i]);
	}
	printf("\n");
}

void printVect6(const Vector6& vect)
{
	for (int i = 0; i < 6; i++)
	{
		printf("%.6f, ", vect[i]);
	}
	printf("\n");
}

void printVect15(const Vector15& vect)
{
	for (int i = 0; i < 15; i++)
	{
		printf("%.6f, ", vect[i]);
	}
	printf("\n");
}

std::string csvStrVect15(const Vector15& vect)
{
	std::ostringstream os;
	for (int i = 0; i < 15; i++)
	{
		os << vect[i];
		os << ",";
	}
	os << "\n";
	std::string str(os.str());
	return str;
}

Vector4 quatFromEul(const Vector3& eul)
{
	 float roll = eul[0];
	 float pitch = eul[1];
	 float yaw = eul[2];

	 float cy = cos(yaw * 0.5f);
	 float sy = sin(yaw * 0.5f);
	 float cr = cos(roll * 0.5f);
	 float sr = sin(roll * 0.5f);
	 float cp = cos(pitch * 0.5f);
	 float sp = sin(pitch * 0.5f);

	 float q0 = cy * cr * cp + sy * sr * sp;
	 float q1 = cy * sr * cp - sy * cr * sp;
	 float q2 = cy * cr * sp + sy * sr * cp;
	 float q3 = sy * cr * cp - cy * sr * sp;

	 return Vector4(q0, q1, q2 ,q3);
}

Vector3 quatToQuatVec(const Vector4& q) // assume quat scalar part quat[0] > 0;
{
	Vector3 qv(q[1], q[2], q[3]);

	if (q[0] < 0.0f)
	{
		float sinHalfAlpha = qv.norm();
		float eps = 1e-8f;
		if (sinHalfAlpha < eps) {
			qv = Vector3(0.0f, 0.0f, 0.0f);
			return qv;
		};
		qv = qv /sinHalfAlpha; // pin
		float alpha = 2.0f * asin(sinHalfAlpha);
		float pi = 3.1415f; // use WGS4 PI here;
		float alphaNew = -2.0f * pi + alpha; // rotate to another dir
		qv = qv * -1.0f; // pin = - pin

		float sinHalfNewAlpha = sin(alphaNew / 2.0f);
		qv = qv * sinHalfNewAlpha;
	}
	return qv;
}

Vector4 quatVecToQuat(const Vector3& qv) // assume quat scalar part quat[0] > 0;
{
	float q0Square = 1 - qv[0] * qv[0] - qv[1] * qv[1] - qv[2] * qv[2];
	if (q0Square < 0.0f) // possible in case of numerical integration error
	{
		q0Square = 0.0f;
	} 
	float q0 = sqrt(q0Square);

	return Vector4(q0, qv[0], qv[1], qv[2]);
}

Vector4 quatMultiply(const Vector4& q, const Vector4& r)
{
	Vector4 p;
	p[0] = r[0] * q[0] - r[1] * q[1] - r[2] * q[2] - r[3] * q[3];
	p[1] = r[0] * q[1] + r[1] * q[0] - r[2] * q[3] + r[3] * q[2];
	p[2] = r[0] * q[2] + r[1] * q[3] + r[2] * q[0] - r[3] * q[1];
	p[3] = r[0] * q[3] - r[1] * q[2] + r[2] * q[1] + r[3] * q[0];
	return p;
}

Vector4 quatInverse(const Vector4& q)
{
	Vector4 qDual (q[0], -q[1], -q[2], -q[3]);
	return qDual;
}

Vector3 quatRotate(const Vector4& q, const Vector3& v)
{
	Vector4 qv(0.0f, v[0], v[1], v[2]);

	Vector4 qDual = quatInverse(q);
	Vector4 qv1 = quatMultiply(qv, qDual);
	Vector4 qv2 = quatMultiply(q, qv1);

	return Vector3(qv2[1], qv2[2], qv2[3]);
}