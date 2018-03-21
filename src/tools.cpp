#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  for (int i = 0; i < estimations.size(); i++) {
    VectorXd diff = estimations[i] - ground_truth[i];
    diff = diff.array()*diff.array(); // element-wise multiplication
    rmse += diff;
  }  
  rmse /= estimations.size(); // mean
  rmse = rmse.array().sqrt(); // element-wise sqrt

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);

  double px, py, vx, vy;
  px = x_state(0);
  py = x_state(1);
  vx = x_state(2);
  vy = x_state(3);

  double pxx_pyy = px*px + py*py;
  double rho = sqrt(pxx_pyy);

  if (pxx_pyy < 0.0001)
    return Hj;

  Hj << (px/rho), (py/rho), 0, 0,
        -(py/pxx_pyy), (px/pxx_pyy), 0, 0,
        py*(vx*py - vy*px)/(pxx_pyy*rho), px*(vy*px - vx*py)/(pxx_pyy*rho), px/rho, py/rho;
  
  return Hj;
}
