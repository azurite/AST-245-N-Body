#include <iostream>
#include <fstream>
#include <hermite.hpp>

Hermite::Hermite()
{
  t = .0;
  dt = 0.01;
  eps = 0.01; // all files and particles have the same softening
  N = 0;
  filename = "";
  lean = false;
}

void Hermite::enableLean()
{
    if(!lean)
      lean = true;
}

void Hermite::disableLean()
{
  if(lean)
    lean = false;
}

void Hermite::integrate(double dt, int numSteps)
{
  std::cout << "starting integration of: " << filename << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "N:              " << N << std::endl;
  std::cout << "dt:             " << dt << std::endl;
  std::cout << "softening:      " << eps << std::endl;
  std::cout << "time interval:  [0, " << (numSteps * dt) << "]" << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << std::endl;

  totalData = MatrixXd::Zero(3, numSteps * N);
  energy = VectorXd::Zero(numSteps);

  this->dt = dt;
  t = 0;

  for(int step = 0; step < numSteps; step++) {
    totalData.block(0, step*N, 3, N) = particles.block(1, 0, 3, N);

    this->computeEnergy(step);
    this->step();
    t += dt;
  }

  // we want the relative energy error
  energy = (((energy(0) * VectorXd::Ones(numSteps)) - energy) / energy(0)).cwiseAbs();

  // write particle positions to file
  if(!lean) {
    std::ofstream file_pos(filename + ".output_pos.txt");
    for(int i = 0; i < totalData.rows(); i++) {
      for(int j = 0; j < totalData.cols(); j++) {
        file_pos << totalData(i, j) << " ";
      }
      file_pos << std::endl;
    }
  }

  // write energy to file
  std::ofstream file_energy(filename + ".output_energy.txt");
  for(int i = 0; i < energy.size(); i++) {
    file_energy << energy(i) << " ";
  }
  file_energy << std::endl;

  // write meta data to file
  std::ofstream file_meta(filename + ".output_meta.txt");
  file_meta << dt << " " << (numSteps * dt) << " " << eps << std::endl;
}

void Hermite::computeEnergy(int step)
{
  double etot = .0;

  double mi, mj, dx, dy, dz, vx, vy, vz;

  for(int i = 0; i < N; i++) {
    double subtract = .0;
    mi = particles(0, i);
    vx = particles(4, i);
    vy = particles(5, i);
    vz = particles(6, i);

    for(int j = 0; j < i; j++) {
      mj = particles(0, j);
      dx = particles(1, j) - particles(1, i);
      dy = particles(2, j) - particles(2, i);
      dz = particles(3, j) - particles(3, i);

      // Assumption: G = 1
      subtract += (mj / std::sqrt(dx*dx + dy*dy + dz*dz + eps*eps));
    }

    etot += (mi * (0.5 * (vx*vx + vy*vy + vz*vz) - subtract));
  }

  energy(step) = etot;
}

MatrixXd Hermite::computeForces(const MatrixXd mass, const MatrixXd pos, const MatrixXd vel)
{
  double mi, mj, dx, dy, dz, dvx, dvy, dvz, sqrtInvDist, sqrtInvDist3, sqrtInvDist5, vdotr, ax, ay, az, jx, jy, jz;

  // first 3 rows hold the acceleration vectors, last 3 rows hold the jerk vectors
  MatrixXd acc_and_jerk = MatrixXd::Zero(6, N);

  for(int i = 0; i < N; i++) {
    for(int j = i + 1; j < N; j++) {
      mj = mass(0, j);
      dx = pos(0, i) - pos(0, j);
      dy = pos(1, i) - pos(1, j);
      dz = pos(2, i) - pos(2, j);
      dvx = vel(0, i) - vel(0, j);
      dvy = vel(1, i) - vel(1, j);
      dvz = vel(2, i) - vel(2, j);

      sqrtInvDist = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz + eps * eps);
      sqrtInvDist3 = sqrtInvDist * sqrtInvDist * sqrtInvDist;

      sqrtInvDist5 = sqrtInvDist3 * sqrtInvDist * sqrtInvDist;
      vdotr = dvx * dx + dvy * dy + dvz * dz;

      // Assumption: G = 1
      ax = -mj * dx * sqrtInvDist3;
      ay = -mj * dy * sqrtInvDist3;
      az = -mj * dz * sqrtInvDist3;

      acc_and_jerk(0, i) += ax;
      acc_and_jerk(1, i) += ay;
      acc_and_jerk(2, i) += az;
      acc_and_jerk(0, j) -= ax;
      acc_and_jerk(1, j) -= ay;
      acc_and_jerk(2, j) -= az;

      // Assumption: G = 1
      jx = -mj * ((dvx * sqrtInvDist3) - (3 * vdotr * dx) * sqrtInvDist5);
      jy = -mj * ((dvy * sqrtInvDist3) - (3 * vdotr * dy) * sqrtInvDist5);
      jz = -mj * ((dvz * sqrtInvDist3) - (3 * vdotr * dz) * sqrtInvDist5);

      acc_and_jerk(3, i) += jx;
      acc_and_jerk(4, i) += jy;
      acc_and_jerk(5, i) += jz;
      acc_and_jerk(3, j) -= jx;
      acc_and_jerk(4, j) -= jy;
      acc_and_jerk(5, j) -= jz;
    }
  }

  return acc_and_jerk;
}

void Hermite::step()
{
  MatrixXd masses = particles.topRows(1);
  x0 = particles.block(1, 0, 3, N);
  v0 = particles.block(4, 0, 3, N);

  MatrixXd acc_and_jerk0 = computeForces(masses, x0, v0);
  a0 = acc_and_jerk0.topRows(3);
  jerk0 = acc_and_jerk0.bottomRows(3);

  // predictor step
  xp = (x0 + (v0 * dt) + (a0 * (0.5 * dt * dt)) + (jerk0 * (0.1667 * dt * dt * dt)));
  vp = (v0 + (a0 * dt) + (jerk0 * (0.5 * dt * dt)));

  MatrixXd acc_and_jerk1 = computeForces(masses, xp, vp);
  a1 = acc_and_jerk1.topRows(3);
  jerk1 = acc_and_jerk1.bottomRows(3);

  // corrector step
  v1 = v0 + ((a0 + a1) * 0.5 * dt) + ((jerk0 - jerk1) * 0.0833 * dt * dt);
  x1 = x0 + ((v0 + v1) * 0.5 * dt) + ((a0 - a1) * 0.0833 * dt * dt);

  // https://gravidy.xyz/include/integrator.html
  // MatrixXd a20 = (-6 * (a0 - a1) - dt * (4 * jerk0 + 2 * jerk1));
  // MatrixXd a30 = (-12 * (a0 - a1) - 6 * dt * (jerk0 + jerk1));
  // x1 = xp + (0.041667 * dt*dt * a20) + (0.00833 * dt*dt * a30);
  // v1 = vp + (0.25 * dt * a20) + (0.041667 * dt * a30);

  // assign next positions and velocities to original particle set
  particles.block(1, 0, 3, N) = x1;
  particles.block(4, 0, 3, N) = v1;
}

bool Hermite::readData(const std::string &filename)
{
  std::ifstream infile(filename);
  int numParticles, numGasParticles, numStarParticles;

  if(infile.is_open()) {
    infile >> numParticles >> numGasParticles >> numStarParticles;

    particles = Eigen::MatrixXd::Zero(HERMITE_DATA_ROWS, numParticles);
    N = numParticles;

    this->filename = filename;

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

    infile >> eps;

    return true;
  }

  std::cout << "Hermite::readData(\"" << filename << "\") " << "could not read file" << std::endl;
  return false;
}
