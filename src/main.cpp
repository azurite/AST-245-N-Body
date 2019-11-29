#include <iostream>
#include <memory>
#include <data.hpp>
#include <particle.hpp>

int main(int argc, char **argv)
{

  std::unique_ptr<const std::vector<Particle>> particles(Data::readFromFile("data.ascii"));

  for(int i = 0; i < 10; i++) {
    std::cout << ((*particles)[i]);
  }

  return 0;
}
