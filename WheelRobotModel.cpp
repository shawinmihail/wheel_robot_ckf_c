#include "WheelRobotModel.h"
#include "QuatMath.h"
#include <math.h>

WheelRobotModel::WheelRobotModel() :
_state(Mat(6, 1)),
_stateDot(Mat(6, 1)),
_fullState(Mat(15, 1)),
_ctrl(Mat(2, 1))
{
	_state.Zero();
	_stateDot.Zero();
	_fullState.Zero();
	_ctrl.Zero();
}

Mat WheelRobotModel::stateToFullState(Mat state, Mat stateDot)
{
	Mat fullState(15, 1);
	// r
	fullState(0, 0) = state(0, 0); // x
	fullState(1, 0) = state(1, 0); // y
	fullState(2, 0) = 0.0;         // z
	// v
	fullState(3, 0) = state(3, 0); // vx
	fullState(4, 0) = state(4, 0); // vy
	fullState(5, 0) = 0.0;         // vz
	// a
	fullState(6, 0) = stateDot(3, 0); // ax
	fullState(7, 0) = stateDot(4, 0); // ay
	fullState(8, 0) = 0.0;            // az
	// qv
	Mat eul(3, 1);
	eul.Zero();
	eul(2, 0) = state(2, 0); // eul = [0;0;th];
	Mat q = quatFromEul(eul);
	fullState(9, 0) = q(1, 0);  // qx
	fullState(10, 0) = q(2, 0); // qy
	fullState(11, 0) = q(3, 0); // qz
	// w = [0;0;th_dot]
	fullState(12, 0) = 0.0;            // wx
	fullState(13, 0) = 0.0;            // wy
	fullState(14, 0) = stateDot(5, 0); // wz
	

	return fullState;
}


Mat WheelRobotModel::model(Mat ctrl, double dt)
{
	Mat stateDot(6,1);
	stateDot(0, 0) = ctrl(0, 0) * cos(_state(2, 0)); // x_dot = v cos (th)
	stateDot(1, 0) = ctrl(0, 0) * sin(_state(2, 0)); // y_dot = v sin (th)
	stateDot(2, 0) = ctrl(0, 0) * ctrl(1, 0);        // th_dot = v u

	// numerical
	stateDot(3, 0) = (stateDot(0, 0) - _state(3,0))  / dt; // (xdot_k - xdot_k-1) / dt
	stateDot(4, 0) = (stateDot(1, 0) - _state(4, 0)) / dt;
	stateDot(5, 0) = (stateDot(2, 0) - _state(5, 0)) / dt;

	return stateDot;
}

void WheelRobotModel::integrateMoutionEuler(double dt, Mat ctrl)
{
	processCtrlInput(ctrl, dt);
	Mat stateDot = model(_ctrl, dt);
	_stateDot.Move(stateDot);

	
	stateDot.ScaleMat2(&dt); // x_dot * dt
	stateDot.Plus(_state, stateDot); // x + x_dot * dt;

	_state.Move(stateDot);
	_fullState.Move(stateToFullState(_state, _stateDot));
}

void WheelRobotModel::processCtrlInput(Mat ctrl, double dt)
{
	double Kv = 1.0;
	double Ku = 1.0;
	double dv = ctrl(0, 0) - _ctrl(0, 0);
	double du = ctrl(1, 0) - _ctrl(1, 0);
	Mat ctrlSmoothed(2, 1);
	ctrlSmoothed(0, 0) = _ctrl(0, 0) + Kv * dv * dt;
	ctrlSmoothed(1, 0) = _ctrl(1, 0) + Ku * du * dt;

	_ctrl.Move(ctrlSmoothed);
}

Mat WheelRobotModel::getState()
{
	return _state;
}

Mat WheelRobotModel::getFullState()
{
	return _fullState;
}