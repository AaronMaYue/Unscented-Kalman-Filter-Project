#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
//small number
#define EPS 0.0001

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
  std_a_ = 6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.45;

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
  is_initialized_ = false;
  //state dimension
  n_x_ = 5;
  // Augmented state dimension
  n_aug_ = 7;
  // number of sigma points
  n_sig_ = 2 * n_aug_ + 1;
  // sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  // predicted sigma points matix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  // weights of sigma points
  weights_ = VectorXd(n_sig_);







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


  if (!is_initialized_)
  {
    x_.fill(0.0);
    P_ = MatrixXd::Identity(5, 5);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho  = meas_package.raw_measurements_[0];
      double phi  = meas_package.raw_measurements_[1];
      //double rhod = meas_package.raw_measurements_[2];

      //while(phi >  M_PI) phi -= 2 * M_PI;
      //while(phi < -M_PI) phi += 2 * M_PI;

      double px = rho * cos(phi);
      double py = rho * sin(phi);

      x_ << px, py, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      while(fabs(x_(0)) < EPS) x_(0) = EPS;
      while(fabs(x_(1)) < EPS) x_(1) = EPS;
    }

    time_us_ = meas_package.timestamp_;

    laser_nis_ = 0.0;

    radar_nis_ = 0.0;

    is_initialized_ = true;
  }

  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_sig_; ++i)
  {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }

  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
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
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;

  //process noise covariance matrix
  MatrixXd Q(2, 2);
  Q << std_a_ * std_a_, 0, 0, std_yawdd_ * std_yawdd_;

  // Augmented P
  MatrixXd P_aug(n_aug_, n_aug_);
  P_aug << P_, MatrixXd::Zero(5, 2), MatrixXd::Zero(2, 5), Q;


  //squar root matrix
  MatrixXd L(n_aug_, n_aug_);
  L = P_aug.llt().matrixL();

  MatrixXd Xsig_aug(n_aug_, n_sig_);
  //cout << Xsig_aug.size() << endl;
  Xsig_aug.col(0) = x_aug;
  //cout << x_aug.size() << endl;

  for (int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i + 1)          = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  for (int i = 0; i < n_sig_; ++i)
  {
    double px       = Xsig_aug(0, i);
    double py       = Xsig_aug(1, i);
    double v        = Xsig_aug(2, i);
    double yaw      = Xsig_aug(3, i);
    double yawd     = Xsig_aug(4, i);
    double nu_a     = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double px_p, py_p;
    if(fabs(yawd) < EPS)
    {
      px_p = px + v * cos(yaw) * delta_t;
      py_p = py + v * sin(yaw) * delta_t;
    }
    else
    {
      px_p = px + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = py + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise process noise
    Xsig_pred_(0, i) = px_p + 0.5 * nu_a * cos(yaw) * delta_t * delta_t;
    Xsig_pred_(1, i) = py_p + 0.5 * nu_a * sin(yaw) * delta_t * delta_t;
    Xsig_pred_(2, i) = v_p + nu_a * delta_t;
    Xsig_pred_(3, i) = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    Xsig_pred_(4, i) = yawd_p + nu_yawdd * delta_t;

  }

  // Predict mean
  x_ = Xsig_pred_ * weights_;

  // Predict covariance
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // normalize angle
    while(x_diff(3) >  M_PI) x_diff(3) -= 2 * M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
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
  int n_z = 2;
  MatrixXd Zsig(n_z, n_sig_);
  Zsig.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
  {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);

    Zsig(0, i) = px;
    Zsig(1, i) = py;
  }

  Update(meas_package, n_z, Zsig);
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
  // Radar sigma points dimesion
  int n_z = 3;
  MatrixXd Zsig(n_z, n_sig_);
  Zsig.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
  {
    double px   = Xsig_pred_(0, i);
    double py   = Xsig_pred_(1, i);
    double v    = Xsig_pred_(2, i);
    double yaw  = Xsig_pred_(3, i);
    //double yawd = Xsig_pred_(4, i);

    Zsig(0, i) = sqrt(px * px + py * py); // rho
    Zsig(1, i) = atan2(py, px); // phi
    Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / sqrt(px * px + py * py); // rho_dot
  }

  Update(meas_package, n_z, Zsig);




}

void UKF::Update(MeasurementPackage meas_package, int n_z, MatrixXd Zsig)
{
  // radar mean
  VectorXd z_pred = Zsig * weights_;
  // radar covariance
  MatrixXd S(n_z, n_z);
  S.fill(0.0);

  for (int i = 0; i < n_sig_; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      while(z_diff(1) >  M_PI) z_diff(1) -= 2 * M_PI;
      while(z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;
    }

    S += weights_(i) * z_diff * z_diff.transpose();
  }


  MatrixXd R(n_z, n_z);
  R.fill(0.0);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radrd_ * std_radrd_;
  }
  else
  {
    R(0, 0) = std_laspx_ * std_laspx_;
    R(1, 1) = std_laspy_ * std_laspy_;
  }

  S += R;

  // Create corss-correlation matrix
  MatrixXd Tc(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while(x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;
    while(x_diff(3) >  M_PI) x_diff(3) -= 2 * M_PI;
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      while(z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;
      while(z_diff(1) >  M_PI) z_diff(1) -= 2 * M_PI;
    }

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    while(z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;
    while(z_diff(1) >  M_PI) z_diff(1) -= 2 * M_PI;
  }

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();


  // caculate NIS
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    radar_nis_ = z_diff.transpose() * S.inverse() * z_diff;
    cout << "radar NIS: " << radar_nis_ << endl;
  }
  else
  {
    laser_nis_ = z_diff.transpose() * S.inverse() * z_diff;
    cout << "laser NIS: " << laser_nis_ << endl;
  }


}

/*void NormalAngle(MeasurementPackage meas_package, VectorXd diff)
{

}*/
