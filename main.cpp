#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#include "Definitions.h"
#include "Utils.h"
#include "WheelRobotModel.h"
#include "MeasurmentModel.h"
#include "Ckf.h"

#include <Eigen/Dense>

void main()
{
    /*Vector4 q1(1, 0, 0, 100);
    Vector4 q2(-5, 0, 0, 100);
    q1 = q1 / q1.norm();
    q2 = q2 / q2.norm();
    std::cout << q1 << "\n\n";
    std::cout << q2 << "\n\n";
    Eigen::Matrix<float, 4, 2> Q;
    Q << q1, q2;
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
    Vector4 qAvg = solver.eigenvectors().col(k);

    std::cout << qAvg << "\n\n";
    std::cout << qAvg.norm() << "\n\n";

    return;*/

    // logging init
    std::ofstream act_state_log;
    act_state_log.open("act_state_log.csv");

    std::ofstream mes_state_log;
    mes_state_log.open("mes_state_log.csv");

    std::ofstream est_state_log;
    est_state_log.open("est_state_log.csv");

    // init
    WheelRobotModel wheelRobotModel;
    MeasurmentModel mesModel;
    Ckf ckf;
    float dt = 1e-2f;
    Vector2 ctrl(0.5f, 0.1f);

    // loop
    for (int i = 0; i < 2000; i++)
    {
        // control input
        ctrl[0]=  0.90f; // v
        ctrl[1] = 0.25f; // u

        // moution
        wheelRobotModel.integrateMoution(ctrl, dt);
        //Vector6 state = wheelRobotModel.getState();
        Vector15 actState = wheelRobotModel.getFullState();

        // mesuarments
        mesModel.setData(actState);
        Vector15 mesState = mesModel.getMesState(); // r v a qv w
        
        // estimation 
        ckf.updateImu(/*a*/mesState.segment(6, 3), /*w*/mesState.segment(12, 3), dt);

        if (i % 10 == 1)
        {
            ckf.updateGps(/*rv*/mesState.segment(0, 6));
        }
        Vector9 estState = ckf.getEstState();
        //std::cout << estState.transpose() << std::endl << std::endl;

        // logging
        act_state_log << csvStrVect15(actState);
        mes_state_log << csvStrVect15(mesState);
        est_state_log << csvStrVect9(estState);

        std::cout << i << std::endl;
    }

    act_state_log.close();
    mes_state_log.close();
    est_state_log.close();
}
