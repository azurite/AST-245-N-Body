#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::Dynamic;

#ifndef GRAVITYSOLVERS_H
#define GRAVITYSOLVERS_H

/*
  MatrixData[:,j] = [m, x, y, z, vx, vy, vz, fx, fy, fz]
*/
#define MATRIX_DATA_ROWS 10

typedef Matrix<float, MATRIX_DATA_ROWS, Dynamic> MatrixData;

namespace Gravitysolver {
  class Direct
  {
  private:
    float epsilon;
    MatrixData particles;
  public:
    Direct();
    bool readDataOld(const std::string &filename);
    bool readData(const std::string &filename);
    bool writeData(const std::string &filename);
    float softening();
    void setSoftening(float eps);
    void solve();
    const MatrixData &data();
  };

  class PM
  {
  private:
    MatrixData particles;
  public:
    PM();
    bool readDataOld(const std::string &filename);
    bool readData(const std::string &filename);
    bool writeData(const std::string &filename);
    void solve();
    const MatrixData &data();
  };
}

#endif // GRAVITYSOLVERS_H
