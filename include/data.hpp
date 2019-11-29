#include <fstream>
#include <vector>
#include <particle.hpp>

#ifndef DATA_H
#define DATA_H

class Data
{
public:
  static std::vector<Particle> *readFromFile(const std::string &filename);
};

#endif // DATA_H
