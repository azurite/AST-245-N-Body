#include <iostream>
#include <vector>
#include <second_task.hpp>
#include <hermite.hpp>

void second_task()
{
  const std::vector<std::string> files = {
    "data/data_n2_circular.dat",
    "data/data_n2_elliptic.dat",
    "data/data_n32.dat",
    "data/data_n64.dat",
    "data/data_n128.dat",
    "data/data_n256.dat"
  };

  const double T = 100;
  const std::vector<double> dt = {
    0.001,
    0.001,
    0.00025,
    0.0005,
    0.001,
    0.002
  };

  Hermite *solver = new Hermite();
  solver->setBlockSize(10);

  for(int i = 0; i < files.size(); i++) {

    if(i == 2)
      solver->setBlockSize(250);

    solver->readData(files[i]);
    solver->integrate(dt[i], T / dt[i]);
  }

  delete solver;
}
