#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#include "Definitions.h"
#include "Utils.h"
#include "WheelRobotModel.h"
#include "MeasurmentModel.h"
#include "Ckf.h"

void main()
{
    //Eigen::HouseholderQR<Eigen::Matrix<float, 3, 2>> qr;
    //Eigen::Matrix<float, 2, 3> A;
    //A << 1, 2, 3, 4, 5, 6;
    //std::cout << A << std::endl << std::endl;

    //Eigen::Matrix<float, 3, 2> AT = A.transpose();
    //qr.compute(AT);

    //Eigen::Matrix<float, 3, 2> R = qr.matrixQR().triangularView<Eigen::Upper>();
    //std::cout <<  R << std::endl << std::endl;

    //Eigen::Matrix<float, 2, 3> L = R.transpose();
    //std::cout << L << std::endl << std::endl;
    //
    //Eigen::Matrix<float, 2, 2> Ld = L.block(0,0,2,2);
    //std::cout << Ld << std::endl << std::endl;

    //return;

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
    for (int i = 0; i < 9999; i++)
    {
        // control input
        ctrl[0]= 0.5f;  // v
        ctrl[1] = 0.1f; // u

        // moution
        wheelRobotModel.integrateMoution(ctrl, dt);
        Vector6 state = wheelRobotModel.getState();
        Vector15 actState = wheelRobotModel.getFullState();

        // mesuarments
        mesModel.setData(actState);
        Vector15 mesState = mesModel.getMesState(); // r v a qv w
        
        // estimation
        ckf.updateImu(/*a*/mesState.segment(6, 3), /*w*/mesState.segment(6, 3), dt);

        if (i % 10 == 9)
        {
            ckf.updateGps(/*rv*/mesState.segment(0, 6));
        }
        Vector9 estState = ckf.getEstState();
        //std::cout << estState.transpose() << std::endl << std::endl;

        // logging
        act_state_log << csvStrVect15(actState);
        mes_state_log << csvStrVect15(mesState);
        est_state_log << csvStrVect9(estState);
    }

    act_state_log.close();
    mes_state_log.close();
    est_state_log.close();

}
