#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Tensor;
using Eigen::VectorXcf;
using Eigen::Vector3i;
using Eigen::Vector3f;

#ifndef GRAVITYSOLVERS_H
#define GRAVITYSOLVERS_H

/*
  MatrixData[:,j] = [m, x, y, z, vx, vy, vz, fx, fy, fz]
*/
#define MATRIX_DATA_ROWS 10
typedef Matrix<float, MATRIX_DATA_ROWS, Dynamic> MatrixData;

#define FIELD_TENSOR_DIMENSIONS 3
typedef Tensor<std::complex<float>, FIELD_TENSOR_DIMENSIONS> FieldTensorCF;
typedef Tensor<float, FIELD_TENSOR_DIMENSIONS> FieldTensorF;

namespace Gravitysolver {
  class DataIO
  {
  protected:
    float epsilon;
    MatrixData particles;
  public:
    DataIO();
    bool readDataOld(const std::string &filename);
    bool readData(const std::string &filename);
    bool writeData(const std::string &filename);
  };

  class Direct : public DataIO
  {
  public:
    Direct();
    float softening();
    void setSoftening(float eps);
    void solve();
    const MatrixData &data();
  };

  class PM : public DataIO
  {
  private:
    int Ng; // number of cells per dimension
    float h; // cell size
    float worldLen; // the size of the smallest enclosing cube of the galaxy
    FieldTensorCF density;
    FieldTensorCF greenFunction;
    FieldTensorCF potential;
    FieldTensorF ax;
    FieldTensorF ay;
    FieldTensorF az;
    Vector3i worldToGrid(float x, float y, float z);
    Vector3f gridToWorld(int i, int j, int k);
    void fft3d(FieldTensorCF &t);
    void ifft3d(FieldTensorCF &t);
    void conv3d(FieldTensorCF &out, FieldTensorCF &in, FieldTensorCF &kernel);
  public:
    PM(int numGridCells);
    void solve();
    const MatrixData &data();
  };
}

#endif // GRAVITYSOLVERS_H
