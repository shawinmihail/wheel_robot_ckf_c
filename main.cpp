#pragma once
#include <stdio.h>
#include "DynMat.h"
#include "QuatMath.h"
#include "WheelRobotModel.h"
#include "MeasurmentModel.h"

#include <iostream>
#include <fstream>
#include <sstream>


void printVectN(Mat mat, int N)
{
    for (int i = 0; i < N; i ++)
    {
        printf("%.6f, ", mat(i, 0));
    }
    printf("\n");
}

std::string csvStrVectN(Mat mat, int N)
{
    std::ostringstream os;
    for (int i = 0; i < N; i++)
    {
        os << mat(i, 0);
        os << ",";
    }
    os << "\n";
    std::string str(os.str());
    return str;
}

void main()
{
    //mat v(3, 1);
    //v(0, 0) = 1.0;
    //v(1, 0) = 0.5;
    //v(2, 0) = 1.0;

    //mat q(4, 1);
    //q(0, 0) = 0.71;
    //q(1, 0) = 0.71;
    //q(2, 0) = 0.0;
    //q(3, 0) = 0.0;
    //
    //mat v2 = quatrotate(q, v);
    //printvectn(v2, 3);

    //return;

    // logging init
    std::ofstream myfile;
    myfile.open("wheel_robot_log.csv");

    // init
    WheelRobotModel wheelRobotModel;
    MeasurmentModel mesModel;
    double dt = 1e-2;
    Mat ctrl(2, 1);

    // process
    for (int i = 0; i < 3999; i++)
    {
        // control input
        ctrl(0, 0) = 1.0; // v
        ctrl(1, 0) = 0.5; // u

        // moution
        wheelRobotModel.integrateMoutionEuler(dt, ctrl);
        Mat state = wheelRobotModel.getState();
        Mat fullState = wheelRobotModel.getFullState();

        // mesuarments
        mesModel.setData(fullState);
        
        // estimation
        // TODO

        //printVectN(state, 6);

        // logging
        myfile << csvStrVectN(state, 6);
    }

    myfile.close();
}
