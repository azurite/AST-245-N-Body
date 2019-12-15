#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <gravitysolvers.hpp>

bool Gravitysolver::Direct::readDataOld(const std::string &filename)
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

  std::cout << "Gravitysolver::Direct::readDataOld(\"" << filename << "\") " << "could not read file" << std::endl;
  return false;
}

bool Gravitysolver::Direct::readData(const std::string &filename)
{
  std::ifstream infile(filename);
  int N;

  if(infile.is_open()) {
    infile >> N >> blockSize >> epsilon;
    particles = Eigen::MatrixXf::Zero(MATRIX_DATA_ROWS, N);

    for(int j = 0; j < N; j++) {
      infile >> particles(0, j); // m
      infile >> particles(1, j); // x
      infile >> particles(2, j); // y
      infile >> particles(3, j); // z
      infile >> particles(4, j); // vx
      infile >> particles(5, j); // vy
      infile >> particles(6, j); // vz
      infile >> particles(7, j); // ax
      infile >> particles(8, j); // ay
      infile >> particles(9, j); // az
      infile >> particles(10, j); // potential
    }

    return true;
  }

  std::cout << "Gravitysolver::Direct::readData(\"" << filename << "\") " << "could not read file" << std::endl;
  return false;
}

bool Gravitysolver::Direct::writeData(const std::string &filename)
{
  std::ofstream outfile(filename);
  int N = particles.cols();

  if(outfile.is_open()) {
    outfile << N << " " << blockSize << " " << epsilon << "\n";

    for(int j = 0; j < N; j++) {
      outfile << particles(0, j) << "\n"; // m
      outfile << particles(1, j) << "\n"; // x
      outfile << particles(2, j) << "\n"; // y
      outfile << particles(3, j) << "\n"; // z
      outfile << particles(4, j) << "\n"; // vx
      outfile << particles(5, j) << "\n"; // vy
      outfile << particles(6, j) << "\n"; // vz
      outfile << particles(7, j) << "\n"; // ax
      outfile << particles(8, j) << "\n"; // ay
      outfile << particles(9, j) << "\n"; // az
      outfile << particles(10, j) << "\n" // potential
    }

    return true;
  }

  std::cout << "Gravitysolver::Direct::writeData(\"" << filename << "\") " << "could not write file" << std::endl;
  return false;
}

void Gravitysolver::Direct::setSoftening(float eps)
{
  epsilon = eps;
}

void Gravitysolver::Direct::setBlockSize(int newBlockSize)
{
  if(newBlockSize >= 1) {
    blockSize = newBlockSize;
  }
  else {
    std::cout << "Gravitysolver::Direct::setBlockSize(" << newBlockSize << ") " << "newBlockSize must be >= 1" << std::endl;
  }
}

void Gravitysolver::Direct::solve()
{
  for(int i = 0; i < particles.cols(); i++) {
    particles(7, i) = .0;
    particles(8, i) = .0;
    particles(9, i) = .0;
  }

  for(int i = 0; i <= particles.cols() - blockSize; i += blockSize) {
    for(int j = 0; j < particles.cols(); j++) {
      float mj = particles(0, j);
      float dx = particles(1, j) - particles(1, i);
      float dy = particles(2, j) - particles(2, i);
      float dz = particles(3, j) - particles(3, i);

      float d = (dx * dx + dy * dy + dz * dz + epsilon * epsilon);
      float d32 = d * std::sqrt(d);

      // Assumption: G = 1
      particles(7, i) += mj * dx / d32;
      particles(8, i) += mj * dy / d32;
      particles(9, i) += mj * dz / d32;
    }

    float x = particles(1, i);
    float y = particles(2, i);
    float z = particles(3, i);
    float ax = particles(7, i);
    float ay = particles(8, i);
    float az = particles(9, i);

    // project the force vector onto the normalized sphere normal
    particles(10, i) = (x*ax + y*ay + z*az) / std::sqrt(x*x + y*y + z*z);
  }
}

const MatrixData &Gravitysolver::Direct::data()
{
  return particles;
}
