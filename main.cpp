#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#include "Definitions.h"
#include "Utils.h"
#include "WheelRobotModel.h"
#include "MeasurmentModel.h"
#include "Ckf.h"

#include "parse.h"

#include <Eigen/Dense>

void main()
{
    parseCsv();
    return;
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
