#include <iostream>
#include <fstream>
#include <hermite.hpp>

Hermite::Hermite()
{
  dt = 0.01;
  eps = 0.01; // all files and particles have the same softening
  N = 0;
}

void Hermite::integrate(int numSteps)
{

}

void Hermite::step()
{
  x0 = MatrixXf::Zero(3, N);
  x1 = MatrixXf::Zero(3, N);
  v0 = MatrixXf::Zero(3, N);
  v1 = MatrixXf::Zero(3, N);
  a0 = MatrixXf::Zero(3, N);
  a1 = MatrixXf::Zero(3, N);
  jerk0 = MatrixXf::Zero(3, N);
  jerk1 = MatrixXf::Zero(3, N);

  float mi, mj, dx, dy, dz, dvx, dvy, dvz, sqrtInvDist, sqrtInvDist3, sqrtInvDist5, vdotr, ax, ay, az, jx, jy, jz;

  for(int i = 0; i < N; i++) {
    for(int j = i + 1; j < N; j++) {
      mj = particles(0, j);
      dx = particles(1, i) - particles(1, j);
      dy = particles(2, i) - particles(2, j);
      dz = particles(3, i) - particles(3, j);
      dvx = particles(4, i) - particles(4, j);
      dvy = particles(5, i) - particles(5, j);
      dvz = particles(6, i) - particles(6, j);

      sqrtInvDist = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz + eps * eps);
      sqrtInvDist3 = sqrtInvDist * sqrtInvDist * sqrtInvDist;

      sqrtInvDist5 = sqrtInvDist3 * sqrtInvDist * sqrtInvDist;
      vdotr = dvx * dx + dvy * dy + dvz * dz;

      // Assumption: G = 1
      ax = -mj * dx * sqrtInvDist3;
      ay = -mj * dy * sqrtInvDist3;
      az = -mj * dz * sqrtInvDist3;

      a0(0, i) += ax;
      a0(1, i) += ay;
      a0(2, i) += az;
      a0(0, j) -= ax;
      a0(1, j) -= ay;
      a0(2, j) -= az;

      // Assumption: G = 1
      jx = -mj * ((dvx * sqrtInvDist3) - (3 * vdotr * dx) * sqrtInvDist5);
      jy = -mj * ((dvy * sqrtInvDist3) - (3 * vdotr * dy) * sqrtInvDist5);
      jz = -mj * ((dvz * sqrtInvDist3) - (3 * vdotr * dz) * sqrtInvDist5);

      jerk0(0, i) += jx;
      jerk0(1, i) += jy;
      jerk0(2, i) += jz;
      jerk0(0, j) -= jx;
      jerk0(1, j) -= jy;
      jerk0(2, j) -= jz;
    }
  }
}

bool Hermite::readData(const std::string &filename)
{
  std::ifstream infile(filename);
  int numParticles, numGasParticles, numStarParticles;

  if(infile.is_open()) {
    infile >> numParticles >> numGasParticles >> numStarParticles;

    particles = Eigen::MatrixXf::Zero(HERMITE_DATA_ROWS, numParticles);
    N = numParticles;

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

  std::cout << "Hermite::readData(\"" << filename << "\") " << "could not read file" << std::endl;
  return false;
}
