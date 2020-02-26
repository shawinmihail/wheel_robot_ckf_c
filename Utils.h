#pragma once
#include <stdio.h>
#include "Definitions.h"

void printVect3(const Vector3& vect); // remove with print double*

void printVect6(const Vector6& vect);

void printVect15(const Vector15& vect);

std::string csvStrVect15(const Vector15& vect);

std::string csvStrVect9(const Vector9& vect);

Vector4 quatFromEul(const Vector3& eul);

Vector3 quatToQuatVec(const Vector4& quat); // assume quat scalar part quat[0] > 0;

Vector4 quatVecToQuat(const Vector3& quat); // assume quat scalar part quat[0] > 0;

Vector4 quatMultiply(const Vector4& q1, const Vector4& q2);

Vector4 quatInverse(const Vector4& q);

Vector3 quatRotate(const Vector4& q, const Vector3& v);



