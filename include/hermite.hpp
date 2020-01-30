#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Dynamic;

#ifndef HERMITE_H
#define HERMITE_H

#define HERMITE_DATA_ROWS 7

class Hermite
{
private:
  int N; // number of particles
  double t; // current simulation time
  double dt; // integration timestep
  double eps; // softening

  Matrix<double, HERMITE_DATA_ROWS, Dynamic> particles;
  Matrix<double, 3, Dynamic> x0; // current positions
  Matrix<double, 3, Dynamic> v0; // current velocities
  Matrix<double, 3, Dynamic> xp; // predicted next positions
  Matrix<double, 3, Dynamic> vp; // predicted next velocities
  Matrix<double, 3, Dynamic> x1; // corrected next positions
  Matrix<double, 3, Dynamic> v1; // corrected next velocities
  Matrix<double, 3, Dynamic> a0; // predicted accelerations
  Matrix<double, 3, Dynamic> a1; // corrected accelerations
  Matrix<double, 3, Dynamic> jerk0; // predicted jerks
  Matrix<double, 3, Dynamic> jerk1; // corrected jerks

  VectorXd energy;
  std::string filename;
  bool lean; // if set to true solver will not write position files because they can get quite large
  int blockSize;
  Matrix<double, 3, Dynamic> totalData; // positions of all particles over all time steps

  void computeEnergy(int step);
  void step();
  MatrixXd computeForces(const MatrixXd mass, const MatrixXd pos, const MatrixXd vel);
public:
  Hermite();
  void enableLean();
  void disableLean();
  void setBlockSize(int size);
  void setSoftening(double newEps);
  void integrate(double dt, int numSteps);
  bool readData(const std::string &filename);
};

#endif // HERMITE_H
