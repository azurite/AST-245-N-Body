#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::Dynamic;

#ifndef GRAVITYSOLVERS_H
#define GRAVITYSOLVERS_H

typedef Matrix<float, 10, Dynamic> MatrixData;

namespace Gravitysolver {
  class Direct
  {
  private:
    float epsilon;
    int blockSize;
    MatrixData particles;
  public:
    bool readDataOld(const std::string &filename);
    bool readData(const std::string &filename);
    bool writeData(const std::string &filename);
    void setSoftening(float eps);
    void setBlockSize(int blockSize);
    void solve();
    const MatrixData &data();
  };
}

#endif // GRAVITYSOLVERS_H
