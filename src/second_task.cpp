#include <iostream>
#include <second_task.hpp>
#include <hermite.hpp>

void second_task()
{
  Hermite *solver = new Hermite();

  solver->enableLean();
  solver->setBlockSize(100);

  solver->readData("data/data_n2_circular.dat");
  solver->integrate(0.001, 100000);

  solver->readData("data/data_n2_elliptic.dat");
  solver->integrate(0.00001, 10000000);

  solver->readData("data/data_n32.dat");
  solver->integrate(0.0001, 1000000);

  solver->readData("data/data_n64.dat");
  solver->integrate(0.0001, 1000000);

  solver->readData("data/data_n128.dat");
  solver->integrate(0.001, 100000);

  solver->readData("data/data_n256.dat");
  solver->integrate(0.01, 10000);
}
