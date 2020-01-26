#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <gravitysolvers.hpp>

// ************************************************************************* //
// ***************************** DIRECT SOLVER ***************************** //
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

Gravitysolver::PM::PM()
{

}

void Gravitysolver::PM::solve()
{

}

const MatrixData &Gravitysolver::PM::data()
{
  return particles;
}
