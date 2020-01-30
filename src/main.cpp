#include <iostream>
#include <chrono>
#include <first_task.hpp>
#include <second_task.hpp>

#include <gravitysolvers.hpp>
#include <hermite.hpp>

void precompute_nbody_forces()
{
  float softening = 4.34998; // mean interparticle separation

  Gravitysolver::Direct *solver = new Gravitysolver::Direct();
  solver->readDataOld("data/data.ascii");

  int total = 7;
  int i = 1;

  std::cout << "Starting N-Body-Forces Precomputation" << std::endl;

  for(float eps = softening / 8; eps <= softening * 8; eps *= 2, i++) {
    solver->setSoftening(eps);

    auto start = std::chrono::high_resolution_clock::now();

    solver->solve();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;

    solver->writeData("data/direct-nbody-" + std::to_string(i) + ".txt");

    std::cout << "Done [" << i << "/" << total << "] time: " << time.count() << "s" << std::endl;
  }

  delete solver;
}

int main(int argc, char **argv)
{
  precompute_nbody_forces();
  first_task();
  second_task();
  return 0;
}
