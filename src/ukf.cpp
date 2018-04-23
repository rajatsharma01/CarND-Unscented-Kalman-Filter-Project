#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // Set to false, will be initialized on first time initialization
  is_initialized_ = false;
  
  // Initialize dimensions for state and augmented state vectors
  n_x_ = 5;
  n_aug_ = 7;

  // Initialize lambda
  lambda_ = 3 - n_aug_;

  // Initialize NIS values
  nis_ = 0.0;

  // Initialize laser measurement update matrix
  H_laser_ = MatrixXd(2, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  // Allocate sigma point prediction matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Allocate weights vector
  weights_ = VectorXd(2 * n_aug_ + 1);

  // Initialize process covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  // Ignore radar measurement if we are not using it 
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && !use_radar_) {
      return;
  }

  // Ignore lidar measurement if we are not using it 
  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && !use_laser_) {
      return;
  }

  // First measurement? Initialize state variables
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        double rho = meas_package.raw_measurements_[0];
        double phi = meas_package.raw_measurements_[1];
        x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],
              0, 0, 0;
    }
    time_us_ = meas_package.timestamp_;

    // First measurement only initializes state, no prediction or measurement update
    is_initialized_ = true;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
  }

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Step 1: Augment state and covariance matrix

  // create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0);
  x_aug.head(n_x_) = x_;

  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;

  // Step 2: Generate sigma points

  // sqrt(P_aug)
  MatrixXd L = P_aug.llt().matrixL();
  
  // create augmented sigma points  
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  double sqrt_lambda = sqrt(lambda_ + n_aug_);
  for (int i = 0; i < L.cols(); i++) {
      Xsig_aug.col(i+1)         = x_aug + sqrt_lambda * L.col(i);
      Xsig_aug.col(i+1+n_aug_)  = x_aug - sqrt_lambda * L.col(i);
  }

  // Step 3: Predict next state for sigma points generated above

  double delta_t2 = delta_t * delta_t;
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      double px = Xsig_aug(0, i);
      double py = Xsig_aug(1, i);
      double v = Xsig_aug(2, i);
      double yaw = Xsig_aug(3, i);
      double yawd = Xsig_aug(4, i);
      double nu_a = Xsig_aug(5, i);
      double nu_yawdd = Xsig_aug(6, i);

      // predicted state noise
      double px_p, py_p;

      // avoid division by zero
      if (fabs(yawd) > 0.0001) {
          px_p = px + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
          py_p = py + v/yawd * (-cos(yaw + yawd * delta_t) + cos(yaw));
      } else {
          px_p = px + v * delta_t * cos(yaw);
          py_p = py + v * delta_t * sin(yaw);
      }

      double v_p    = v;
      double yaw_p  = yaw + yawd * delta_t;
      double yawd_p = yawd;

      // Add noise
      px_p      += 0.5 * nu_a * delta_t2 * cos(yaw);
      py_p      += 0.5 * nu_a * delta_t2 * sin(yaw);
      v_p       += nu_a * delta_t;
      yaw_p     += 0.5 * nu_yawdd * delta_t2;
      yawd_p    += nu_yawdd * delta_t;

      // Save them to sigma matrix
      Xsig_pred_(0, i) = px_p;
      Xsig_pred_(1, i) = py_p;
      Xsig_pred_(2, i) = v_p;
      Xsig_pred_(3, i) = yaw_p;
      Xsig_pred_(4, i) = yawd_p;
  }

  // Step 4: Predict state mean and covariance

  // setup weights
  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) *= 2 * lambda_;

  // Predict state mean
  x_ = Xsig_pred_ * weights_;

  // predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      VectorXd xdiff = Xsig_pred_.col(i) - x_;
      AngleNormalize(&xdiff(3));
      P_ += weights_(i) * xdiff * xdiff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // Lidar measurements are linear to state vector so no need to use UKF update rules
  // we stick to basic Kalman Filter equation for lidar update
  
  // Predicted state only includes px and py position
  VectorXd z_pred = H_laser_ * x_;
  
  // Measurement error
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  
  // Projection of covariance matrix in measurement space
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht;

  // Add measurement noise to uncentainity matrix
  S(0, 0) += std_laspx_ * std_laspx_;
  S(1, 1) += std_laspy_ * std_laspy_;

  // Kalman Gain
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  // new estimate
  x_ += K*y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;

  // Update NIS value with current measurement
  nis_ = y.transpose() * Si * y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // Step 1: Transform sigma points into measurement space
  
  int n_z = 3; // measurement dimension for radar
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double yaw = Xsig_pred_(3, i);

      double rho = sqrt(px*px + py*py);
      if (fabs(rho) < 0.0001) {
          rho = 0.0001;
      }
      Zsig(0, i) = rho;
      Zsig(1, i) = atan2(py, px);
      Zsig(2, i) = v * (px * cos(yaw) + py * sin(yaw)) / rho;
  }

  // Step 2: Predict next position for sigma points in measurement space

  // mean state prediction
  VectorXd z_pred = Zsig * weights_;

  // predicted measurement covariance. We also calculate cross correlation matrix here
  // to save one extra pass over Zsig matrix.
  MatrixXd S = MatrixXd(n_z, n_z);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  S.fill(0.0);
  Tc.fill(0.0);
  for (int i = 0; i < Zsig.cols(); i++) {
      VectorXd zdiff = Zsig.col(i) - z_pred;
      AngleNormalize(&zdiff(1));

      VectorXd xdiff = Xsig_pred_.col(i) - x_;
      AngleNormalize(&xdiff(3));

      MatrixXd zdiff_trans = zdiff.transpose();
      S += weights_(i) * zdiff * zdiff_trans;
      Tc += weights_(i) * xdiff * zdiff_trans;
  }

  // Add measurement noise to co-variance matrix
  S(0, 0) += std_radr_ * std_radr_;
  S(1, 1) += std_radphi_ * std_radphi_;
  S(2, 2) += std_radrd_ * std_radrd_;

  // Step 3: UKF update state and co-variance matrix

  // Kalman Gain
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc * Si;

  // Finally, update state mean and covariance matrix
  VectorXd zdiff = meas_package.raw_measurements_ - z_pred;
  AngleNormalize(&zdiff(1));
  x_ += K * zdiff;
  P_ -= K * S * K.transpose();

  // Update NIS value with current measurement
  nis_ = zdiff.transpose() * Si * zdiff;
}
