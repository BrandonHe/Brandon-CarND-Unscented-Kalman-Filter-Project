
#include "tools.h"
#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // Initially set to false
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // Set time state, in us
  time_us_ = 0.0;

  // Set state dimension
  n_x_ = x_.size();

  // Set augmented dimension
  n_aug_ = n_x_ + 2;

  // Define Sigma points spreading parameter
  lambda_ = 3 - n_aug_;

  // Define Sigma points number
  n_sig_ = 2 * n_aug_ + 1;

  // Set predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  // Weights of sigma points
  weights_ = VectorXd(n_sig_);

  // NIS for radar
  NIS_radar_ = 0.0;

  // NIS for laser
  NIS_laser_ = 0.0;

  cout << "Finished initializes Unscented Kalman filter!" << endl;
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
  if (!is_initialized_) {
    // initialize covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

     // CTRV(Constant Turn Rate and Velocity Magnitude Model)
     // x_ = [px, py, v, yaw, yaw_dot]     

    if (MeasurementPackage::RADAR == meas_package.sensor_type_) {
      /**
       *Convert radar from polar coordinates to cartesian coordinates
       */      
      float rho = meas_package.raw_measurements_[0];     // range
      float phi = meas_package.raw_measurements_[1];     // bearing
      float rho_dot = meas_package.raw_measurements_[2]; // velocity of rho
      // Coordinates convertion from polar to cartesian 
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      float v = sqrt(vx * vx + vy * vy);
      x_ << px, py, v, 0, 0;
    } else if (MeasurementPackage::LASER == meas_package.sensor_type_) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      if (fabs(x_(0) < 0.001) || fabs(x_(1) < 0.001)) {
        x_(0) = 0.001;
        x_(1) = 0.001;
      }

    }

  // Initialize weights
  double weight_0 = lambda_/(lambda_+ n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i< n_sig_; ++i) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    //cout << "Done initialize!" << endl;

    return;
  }

  /***************************
   * Prediction 
   ***************************/
  // Calculation the timestamp between measurements in seconds
  double dt = meas_package.timestamp_ - time_us_;
  // Convert to seconds
  dt /= 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  /***************************
   * Update
   ***************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    
    cout << "MeasurementPackage:: RADAR update!" << endl;
    //cout << "MeasurementPackage:: RADAR update:" << endl <<
    //  meas_package.raw_measurements_[0] << 
    //  meas_package.raw_measurements_[1] << endl;

    UpdateRadar(meas_package);
  }
  if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {

    cout << "MeasurementPackage:: LASER update!" << endl;
    //cout << "MeasurementPackage:: LASER update:" << endl <<
    //  meas_package.raw_measurements_[0] << 
    //  meas_package.raw_measurements_[1] << endl;

    UpdateLidar(meas_package);
  }

  cout << "Finished UKF::ProcessMeasurement! " << endl;

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

  cout << "Start UKF::Prediction! " << endl;

  /********************************************
   * Generate sigma points
   ********************************************/
  
  // Create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2*n_x_+1);

  // Caculate square root of P
  MatrixXd A = P_.llt().matrixL();

  // Set lambda for non-augmented sigma points
  lambda_ = 3 - n_x_;

  // Create sigma points
  Xsig.col(0) = x_;
  for(int i =0; i < n_x_; ++i) {
    Xsig.col(i+1)       = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    Xsig.col(i+1+n_x_)  = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }

  /********************************************
   * Augment sigma points
   ********************************************/



  // Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // Create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Create sigma points matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  lambda_ = 3 - n_aug_;
  // Create matrix with predicted sigma points as columns
  //MatrixXd Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;


  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  //P_aug(n_x_,n_x_) = std_a_ * std_a_;
  //P_aug(n_x_+1,n_x_+1) = std_yawdd_ * std_yawdd_ ;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_ ;

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    VectorXd sqrt_lambda_a_aug_L = sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1)       = x_aug + sqrt_lambda_a_aug_L;
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt_lambda_a_aug_L;
  }

  /************************************************
   * Predict sigma points
   ************************************************/

  // Predict sigma points
  for (int i = 0; i< n_sig_; ++i)
  {
    // Xxtract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // Predicted state values
    double px_p, py_p;

    // Avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // Add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
 
    // Write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  cout << "Xsig_pred = " << endl << Xsig_pred_ << endl;

  /***************************************************
   * Convert predicted sigma points to mean/covariance
   **************************************************/

  // Predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  cout << "Predicted state:" << endl;
  cout << "x_ =: "<< endl << x_ << endl;

  // Predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) {  //iterate over sigma points

    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    cout << "Before angle normalization x_diff" << endl << x_diff << endl;

    // Angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2. * M_PI;

    cout << "After angle normalization x_diff" << endl << x_diff << endl;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }


  cout << "Predicted covariance matrix:" << endl;
  cout << "P_ =: " << endl << P_ << endl;
  cout << "Finished UKF::Prediction!" << endl;
}


/**
  * Update the state and the state covariance matrix using a measument Package
  * @param meas_package The measurement at k+1
  * @param Zsig The matrix for sigma points in measurement space
  * @param n_z The measurement dimension
  */
void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd Zsig, int n_z) {

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Measurement
  VectorXd z = meas_package.raw_measurements_;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //calculate NIS
  if (MeasurementPackage::RADAR == meas_package.sensor_type_) {
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;  
  } else if (MeasurementPackage::LASER == meas_package.sensor_type_) {
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  }

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
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

  // Set measurement dimension
  int n_z = 2;

  // Create matrix for sigma points in measurement sapce
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  // Transform sigma points into measurement space
  for (int i = 0; i < n_sig_; ++i) {
    // extract values for better readibiality
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurment model
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }

  UpdateUKF(meas_package, Zsig, n_z);
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
  
  // Set measurement dimension, radar can measure rho, phi and rho_dot
  int n_z = 3;

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  // Transform sigma points into measurement space
  for (int i = 0; i < n_sig_; ++i) {
    //extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v_x = v*cos(yaw);
    double v_y = v*sin(yaw);

    // Measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y); //rho
    Zsig(1, i) = atan2(p_y, p_x); //phi
    Zsig(2, i) = (p_x* v_x + p_y*v_y) /  sqrt(p_x*p_x + p_y*p_y); //rho_dot
  }
  UpdateUKF(meas_package, Zsig, n_z);
}
