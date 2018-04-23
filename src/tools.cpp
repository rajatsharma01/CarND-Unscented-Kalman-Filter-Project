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
  rmse.fill(0);

  if (estimations.size() != ground_truth.size() ||
      estimations.size() == 0) {
      std::cerr << "Invalid estimation or ground truth data" << std::endl;
      return rmse;
  }

  for (unsigned int i=0; i<estimations.size(); i++) {
      VectorXd diff = estimations[i] - ground_truth[i];
      diff = diff.array()*diff.array();
      rmse += diff;
  }

  rmse /= estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}
