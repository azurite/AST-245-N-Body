#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::MatrixXf;
using Eigen::Dynamic;

#ifndef HERMITE_H
#define HERMITE_H

#define HERMITE_DATA_ROWS 7

class Hermite
{
private:
  int N; // number of particles
  float dt; // integration timestep
  float eps; // softening
  Matrix<float, HERMITE_DATA_ROWS, Dynamic> particles;
  Matrix<float, 3, Dynamic> x0; // predicted next positions
  Matrix<float, 3, Dynamic> x1; // corrected next positions
  Matrix<float, 3, Dynamic> v0; // predicted next velocities
  Matrix<float, 3, Dynamic> v1; // corrected next velocities
  Matrix<float, 3, Dynamic> a0; // predicted accelerations
  Matrix<float, 3, Dynamic> a1; // corrected accelerations
  Matrix<float, 3, Dynamic> jerk0; // predicted jerks
  Matrix<float, 3, Dynamic> jerk1; // corrected jerks
public:
  Hermite();
  void step();
  void integrate(int numSteps);
  bool readData(const std::string &filename);
};

#endif // HERMITE_H
