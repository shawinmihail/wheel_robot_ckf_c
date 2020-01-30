#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#include "Definitions.h"
#include "Utils.h"
#include "WheelRobotModel.h"
//#include "MeasurmentModel.h"


// 1200 22;

void main()
{
    //Vector15 v;
    //v[14] = 33;
    //printVect15(v);

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
    //MeasurmentModel mesModel;
    float dt = 1e-2f;
    Vector2 ctrl(0.5f, 0.1f);

    // loop
    for (int i = 0; i < 3999; i++)
    {
        // control input
        ctrl[0]= 0.5f;  // v
        ctrl[1] = 0.1f; // u

        // moution
        wheelRobotModel.integrateMoution(ctrl, dt);
        Vector6 state = wheelRobotModel.getState();
        Vector15 actState = wheelRobotModel.getFullState();

        // mesuarments
        //mesModel.setData(actState);
        //Mat mesState = mesModel.getMesState();
        
        // estimation
        // TODO

        //printVectN(state, 6);

        // logging
        act_state_log << csvStrVect15(actState);
        //mes_state_log << csvStrVectN(mesState, 15);
    }

    act_state_log.close();
    mes_state_log.close();
    est_state_log.close();

}
