#include <iostream>
#include <second_task.hpp>
#include <hermite.hpp>

void second_task()
{
  Hermite *solver = new Hermite();
  solver->readData("data/data_n2_circular.dat");
  solver->integrate(0.001, 1000);

  solver->readData("data/data_n2_elliptic.dat");
  solver->integrate(0.00001, 100);

  solver->readData("data/data_n32.dat");
  solver->integrate(0.01, 100);
}
