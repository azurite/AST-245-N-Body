# AST-245 N-Body

This Repository contains the code for the N-Body project from the [UZH](https://www.uzh.ch/en.html) course on [Computational Astrophysics](http://www.vvz.ethz.ch/Vorlesungsverzeichnis/lerneinheit.view?semkez=2019W&ansicht=LEHRVERANSTALTUNGEN&lerneinheitId=132986&lang=en).

## Dependencies

You'll need [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for the code to compile. If needed change the `EIGEN` variable in the makefile to the path of your eigen source files. You'll need [MathGL 2.0](http://mathgl.sourceforge.net/) or higher to plot the data of the first task and [Python 3](https://www.python.org/download/releases/3.0/) or higher to plot the data of the second task.

## Usage

Simply compile the code with `make` and run it with `./main.out`. The plots for the first task are generated directly. The second task produces output files with data. To plot that data use `python3 plot.py`. The code will perform the necessary tasks for the project but the solvers are fully encapsulated and can be used on their own.

## Documentation

### Gravitysolver::Direct

Computes the n-body forces from a set of initial conditions with the direct summation algorithm in `O(n²)` where n is the number of particles. The initial conditions can be read in two ways

#### bool Gravitysolver::Direct::readDataOld( \<filename\> )

Reads a set of initial conditions in the format according to the data provided by the project. The data file is made of a header and sequential arrays, as follows: Number of particles(N), Number of gas particles(NG), Number of star particles (NS).

```
N NG NS
Masses[i]
x[i]
y[i]
z[i]
Vx[i]
Vy[i]
Vz[i]
```

Note that the data first runns through all N masses then through all N x coordinates, all N y coordinates, all N z coordinates etc. The function returns `true` when successful and `false` otherwise. In the case of an unsuccessful read the console should give some output.

#### bool Gravitysolver::Direct::readData( \<filename\> )

Reads a header in the format `N eps` where N is the number of particles and eps is the softening. It then reads the data particle by particle so `mass, x, y, z, vx, vy, vz` of particle i then `mass, x, y, z, vx, vy, vz` of particle i+1 and so on. The function returns `true` when successful and `false` otherwise. In the case of an unsuccessful read the console should give some output.

#### bool Gravitysolver::Direct::writeData( \<filename\> )

Writes the particles to a file in the same format as being used by the `readData()` function above. The function returns `true` when successfully writing the file and `false` otherwise. In the case of an unsuccessful write the console should give some output.

#### float Gravitysolver::Direct::softening()

Returns the currently used softening for the solver.

#### void Gravitysolver::Direct::setSoftening( \<float\> )

Sets the softening value used by the solver to the function's input.

#### Gravitysolver::Direct::solve()

Computes the n-body forces acting on each particle with the direct summation algorithm.

#### const MatrixData &MatrixGravitysolver::Direct::data()

Returns the current state (mass, position, velocity, force) of all particles in the form of a 10xN Eigen Matrix `Matrix<float, 10, Dynamic>` where the i-th column holds the state of the i-th particle in the form of `[mass, x, y, z, vx, vz, vz, fx, fy, fz]`. Note that until either `readDataOld()` or `readData()` have been called the output of the function remains uninitialized.

### Gravitysolver::PM

* **Note this solver has not been tested for correctness and should therefore not be used!** Right now it uses the unsupported Eigen Tensor for the mesh which is not nearly as efficient as the Eigen matrices.

Computes the n-body forces from a set of initial conditions by constructing a particle mesh whose size in one dimension is specified in the constructor. It computes the forces by solving the poisson equation with a fast fourier transform. It scales as `O(n³log(n))` where n is the number of mesh cells along one dimension. The interface is identical to the direct solver class except that the methods to get and set the softening are missing.

### Hermite

Reads a set of initial conditions, evolves the system for a specified amount of time and writes the state of the system and the relative energy error over the entire time (for each timestep) to ouput files. One file for the particle positions and one file for the energy error.

#### bool Hermite::readData( \<filename\> )

Reads a set of initial conditions from a file with the same format as specified in `Gravitysolver::Direct::readDataOld`. The very last value of the file can hold the softening if you whish to do so.

#### void Hermite::setSoftening( \<double\> )

Sets the softening for the evolution to the function's input value

#### void Hermite::enableLean()

Only writes the energy error to a file. Can be useful if the file for the positions becomes extremely large and one want's to avoid writing it.

#### void Hermite::disableLean()

Inverse of `enableLean()`. If calling the simulation function now the particles positions are written to a file again.

#### void Hermite::setBlockSize( \<int\> )

Only keeps track of every k-th iteration in the simulation where k is the function's input. This allows to save a lot of memory when usign very small step sizes.

#### void Hermite::integrate( \<double\>, \<int\>)

Calling `integrate(dt, numSteps)` evolves the system with timestep `dt` for `numSteps` steps and (depending on the lean mode) writes the result of the simulation to a file.

#### const Matrix<double, 3, Dynamic> &Hermite::data()

Returns the state of the system over the entire simulation time. Each column contains the [x, y, z] position of the particles. The columns 0 to n-1 contain the state of n particles at the first time step, columns n to 2n-1 contain the state of n particles at the next (or a multiple of next depending on the block size) time step and so on.

### Example Usage

```cpp
#include <gravitysolvers.hpp>
#include <hermite.hpp>

int main()
{
    Gravitysolver::Direct *direct = new Gravitysolver::Direct();

    direct->readData("data.txt");
    direct->setSoftening(0.001);
    direct->solve();
    direct->writeData("forces.txt");

    Hermite *hermite = new Hermite();

    hermite->readData("data.txt");
    hermite->setSoftening(0.001);
    hermite->setBlockSize(20);
    hermite->enableLean();
    hermite->integrate(0.01, 1000);

    return 0;
}
```
