#include <iostream>
#include <data.hpp>

int main(int argc, char **argv)
{
  Data *data = new Data("data.ascii");
  delete data;

  return 0;
}
