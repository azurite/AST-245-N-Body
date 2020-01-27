#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <gravitysolvers.hpp>

// ************************************************************************* //
// ******************************** DATA IO ******************************** //
// ************************************************************************* //

Gravitysolver::DataIO::DataIO()
{
  epsilon = 0.0;
}

bool Gravitysolver::DataIO::readDataOld(const std::string &filename)
{
  std::ifstream infile(filename);
  int numParticles, numGasParticles, numStarParticles;

  if(infile.is_open()) {
    infile >> numParticles >> numGasParticles >> numStarParticles;

    particles = Eigen::MatrixXf::Zero(MATRIX_DATA_ROWS, numParticles);

    for(int i = 0; i < numParticles; i++) {
      infile >> particles(0, i); // m
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> particles(1, i); // x
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> particles(2, i); // y
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> particles(3, i); // z
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> particles(4, i); // vx
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> particles(5, i); // vy
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> particles(6, i); // vz
    }

    return true;
  }

  std::cout << "Gravitysolver::DataIO::readDataOld(\"" << filename << "\") " << "could not read file" << std::endl;
  return false;
}

bool Gravitysolver::DataIO::readData(const std::string &filename)
{
  std::ifstream infile(filename);
  int N;

  if(infile.is_open()) {
    infile >> N >> epsilon;
    particles = Eigen::MatrixXf::Zero(MATRIX_DATA_ROWS, N);

    float m, x, y, z, vx, vy, vz, fx, fy, fz;

    for(int j = 0; j < N; j++) {
      infile >> m >> x >> y >> z >> vx >> vy >> vz >> fx >> fy >> fz;
      particles(0, j) = m;
      particles(1, j) = x;
      particles(2, j) = y;
      particles(3, j) = z;
      particles(4, j) = vx;
      particles(5, j) = vy;
      particles(6, j) = vz;
      particles(7, j) = fx;
      particles(8, j) = fy;
      particles(9, j) = fz;
    }

    return true;
  }

  std::cout << "Gravitysolver::DataIO::readData(\"" << filename << "\") " << "could not read file" << std::endl;
  return false;
}

bool Gravitysolver::DataIO::writeData(const std::string &filename)
{
  std::ofstream outfile(filename);
  int N = particles.cols();

  if(outfile.is_open()) {
    outfile << N << " " << epsilon << "\n";

    for(int j = 0; j < N; j++) {
      outfile << particles(0, j) << "\n"; // m
      outfile << particles(1, j) << "\n"; // x
      outfile << particles(2, j) << "\n"; // y
      outfile << particles(3, j) << "\n"; // z
      outfile << particles(4, j) << "\n"; // vx
      outfile << particles(5, j) << "\n"; // vy
      outfile << particles(6, j) << "\n"; // vz
      outfile << particles(7, j) << "\n"; // fx
      outfile << particles(8, j) << "\n"; // fy
      outfile << particles(9, j) << "\n"; // fz
    }

    return true;
  }

  std::cout << "Gravitysolver::DataIO::writeData(\"" << filename << "\") " << "could not write file" << std::endl;
  return false;
}

// ************************************************************************* //
// ***************************** DIRECT SOLVER ***************************** //
// ************************************************************************* //

Gravitysolver::Direct::Direct()
{

}

float Gravitysolver::Direct::softening()
{
  return epsilon;
}

void Gravitysolver::Direct::setSoftening(float eps)
{
  epsilon = eps;
}

void Gravitysolver::Direct::solve()
{
  for(int i = 0; i < particles.cols(); i++) {
    particles(7, i) = .0;
    particles(8, i) = .0;
    particles(9, i) = .0;
  }

  float mi, mj, dx, dy, dz, sqrtInvDist, sqrtInvDist3, fx, fy, fz;

  for(int i = 0; i < particles.cols(); i++) {
    for(int j = i + 1; j < particles.cols(); j++) {
      mi = particles(0, i);
      mj = particles(0, j);
      dx = particles(1, i) - particles(1, j);
      dy = particles(2, i) - particles(2, j);
      dz = particles(3, i) - particles(3, j);

      sqrtInvDist = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz + epsilon * epsilon);
      sqrtInvDist3 = sqrtInvDist * sqrtInvDist * sqrtInvDist;

      // Assumption: G = 1
      fx = -mi * mj * dx * sqrtInvDist3;
      fy = -mi * mj * dy * sqrtInvDist3;
      fz = -mi * mj * dz * sqrtInvDist3;

      particles(7, i) += fx;
      particles(8, i) += fy;
      particles(9, i) += fz;
      particles(7, j) -= fx;
      particles(8, j) -= fy;
      particles(9, j) -= fz;
    }
  }
}

const MatrixData &Gravitysolver::Direct::data()
{
  return particles;
}

// ************************************************************************* //
// ******************************* PM SOLVER ******************************* //
// ************************************************************************* //

Gravitysolver::PM::PM(int numGridCells)
{
  if(numGridCells > 0) {
    Ng = numGridCells;
  }
  else {
    std::cout << "Gravitysolver::PM::PM(" << numGridCells << ") numGridCells must be > 0" << std::endl;
    Ng = 1;
  }

  h = .0;
}

/**
* Returns indices for a mesh cell for the galaxy given world coordinates.
* World coordinates have the center of the galaxy at (0, 0, 0).
* The mapping to grid coordinates is an affine mapping where the galaxy center
* (0, 0, 0) is in the center of the mesh as well.
*/
Vector3i Gravitysolver::PM::worldToGrid(float x, float y, float z)
{
  Vector3i gridCoor;

  gridCoor(0) = std::floor((x + (worldLen / 2)) / worldLen * Ng);
  gridCoor(1) = std::floor((y + (worldLen / 2)) / worldLen * Ng);
  gridCoor(2) = std::floor((z + (worldLen / 2)) / worldLen * Ng);

  return gridCoor;
}

/**
* Returns position at the center of a mesh cell in world coordinates given cell indices
* The position is mapped in a fashin such that (0, 0, 0) is the center of the galaxy
* in world coordinates.
*/
Vector3f Gravitysolver::PM::gridToWorld(int i, int j, int k)
{
  float h = worldLen / Ng; // cell size
  Vector3f worldCoor;

  worldCoor(0) = (i * worldLen / Ng) - (worldLen / 2) + (h / 2);
  worldCoor(1) = (j * worldLen / Ng) - (worldLen / 2) + (h / 2);
  worldCoor(2) = (k * worldLen / Ng) - (worldLen / 2) + (h / 2);

  return worldCoor;
}

void Gravitysolver::PM::fft3d(FieldTensor &t)
{
  const int x = t.dimension(0);
  const int y = t.dimension(1);
  const int z = t.dimension(2);

  Eigen::FFT<float> fft;

  for(int k = 0; k < z; k++) { // for each 2d sheet make a 2d fft
    for(int j = 0; j < y; j++) { // fft in x-dir
      VectorXcf tv(x);
      for(int i = 0; i < x; i++)
          tv(i) = t(i, j, k);

      VectorXcf fv = fft.fwd(tv);
      for(int i = 0; i < x; i++)
          t(i, j, k) = fv(i);
    }

    for(int i = 0; i < x; i++) { // fft in y-dir
      VectorXcf tv(y);
      for(int j = 0; j < y; j++)
          tv(j) = t(i, j, k);

      VectorXcf fv = fft.fwd(tv);
      for(int j = 0; j < y; j++)
          t(i, j, k) = fv(j);
    }
  }

  for(int i = 0; i < x; i++) { // and for each of the x*y spikes pointing upwards in z-dir do a 1D fft
    for(int j = 0; j < y; j++) {

      VectorXcf tv(z);
      for(int k = 0; k < z; k++)
        tv(k) = t(i, j, k);

      VectorXcf fv = fft.fwd(tv);
      for(int k = 0; k < z; k++)
        t(i, j, k) = fv(k);
    }
  }
}

void Gravitysolver::PM::ifft3d(FieldTensor &t)
{
  const int x = t.dimension(0);
  const int y = t.dimension(1);
  const int z = t.dimension(2);

  const float invXYZ = 1.0 / (x*y*z);

  for(int i = 0; i < x; i++) {
    for(int j = 0; j < y; j++) {
      for(int k = 0; k < z; k++) {
        t(i,j,k) = std::conj(t(i,j,k));
      }
    }
  }

  fft3d(t);

  for(int i = 0; i < x; i++) {
    for(int j = 0; j < y; j++) {
      for(int k = 0; k < z; k++) {
        t(i,j,k) = std::conj(t(i,j,k)) * invXYZ;
      }
    }
  }
}

void Gravitysolver::PM::conv3d(FieldTensor &out, FieldTensor &in, FieldTensor &kernel)
{
  int x = in.dimension(0);
  int y = in.dimension(1);
  int z = in.dimension(2);

  fft3d(in);
  fft3d(kernel);

  for(int i = 0; i < x; i++) {
    for(int j = 0; j < y; j++) {
      for(int k = 0; k < z; k++) {
        out(i,j,k) = in(i,j,k) * kernel(i,j,k);
      }
    }
  }

  ifft3d(out);
}

void Gravitysolver::PM::solve()
{
  const int N = particles.cols();
  Vector3f maxPos = particles.block(1, 0, 3, N).rowwise().maxCoeff();
  Vector3f minPosAbs = particles.block(1, 0, 3, N).rowwise().minCoeff().cwiseAbs();

  worldLen = std::max(maxPos.maxCoeff(), minPosAbs.maxCoeff()) * 2;

  // adds one layer of cells in each dimension so particles at the edge will
  // contribute with their entire mass to the density field
  worldLen += (worldLen / Ng);
  h = worldLen / Ng;

  std::cout << "PM solver global parameters" << std::endl;
  std::cout << "---------------------------" << std::endl;
  std::cout << "worldLen:   " << worldLen << std::endl;
  std::cout << "Ng:         " << Ng << std::endl;
  std::cout << "h:          " << h << std::endl;
  std::cout << "---------------------------" << std::endl;
}

const MatrixData &Gravitysolver::PM::data()
{
  return particles;
}
