
# Installation

The installation instructions here include steps for both ACEtb.jl package and `lib_acetb` API compilation in DFTB+ code.

## ACEtb.jl Installation 

ACEtb.jl is a pure Julia code and can be directly installed using the Julia package following the instructions below.
We recommend Julia with version 1.5 and higher to use with ACEtb.jl package as the code is only tested with version 1.5 and higher (as of now, we also tested it with 1.6).

```
] add https://github.com/ACEsuit/ACEtb.jl.git
```

Here, `]` indicates that REPL is in `pkg` mode at Julia.

## Testing ACEtb.jl installation

ACEtb package includes several unit tests that can be run using Julia `pkg` testing command as follows:

```
] test ACEtb
```

## Installing DFTB+ with ACEtb

DFTB+ integration of ACEtb needs the compilation of `lib_acetb` (currently, within DFTB+) code. This section only covers the installation of DFTB+ with ACEtb module and will not provide information on how DFTB+ is compiled with its other modules. If you are interested in the detailed explanation of how the integration is done, please check DFTB+ Integration section.

### Dependencies of DFTB+ and `lib_acetb`

As well as Julia with version 1.5 and higher is needed and ACEtb is installed in main environment (@v1.5 or equivalent in Julia REPL), the following libraries are also needed to compile DFTB+.

Libraries:
- OpenBLAS
- MPI (Ex.: OpenMPI. Only if parallelisation of DFTB+ is needed.)
- ScaLAPACK (Only needed if DFTB+ is compiled with MPI)

#### Installing Dependencies at macOS

The dependencies may be installed using Homebrew with the following commands:
```
brew install openblas
brew install openmpi
brew install scalapack
```

#### Installing Dependencies at Linux (Ubuntu)
 
The dependencies may be installed using Homebrew with the following commands:
```
apt-get update
apt-get install make build-essential g++ gfortran
apt-get install libblas-dev liblapack-dev libopenmpi-dev libscalapack-mpi-dev openmpi-bin
```

### Installation of DFTB+ with `WITH_ACETB=TRUE` for `lib_acetb`

#### Compiling DFTB+ Serial

In `dftbplus` directory of DFTB+ repository, run the following commands:

```
mkdir build
cd build
cmake -DWITH_ACETB=TRUE ..
make
make install
```

Make sure that the configuration step of CMake run includes the following outputs where Julia and LAPACK libraries are found:

```

```

In case LAPACK libraries cannot be found by CMake, specify the correct library path using `-DLAPACK_LIBRARY_DIRS` as following:
```
cmake -DWITH_ACETB=TRUE -DLAPACK_LIBRARY_DIRS=/usr/local/opt/openblas/lib ..
``` 
Here, `/usr/local/opt/openblas/lib` is the OpenBLAS library path that is installed by Homebrew at macOS. Cange the directory accordingly to match your system's installation path.

#### Compiling DFTB+ Parallel using MPI

Here we assume you use OpenMPI library. In case you would like to use Intel MPI or MPICH, make sure CMake can find the correct library paths.

A typical command to compile DFTB+ with MPI is as follows
```
cmake -DWITH_MPI=TRUE -DWITH_ACETB=TRUE ..
make
make install
```

With the above parameters, CMake try to find the ScaLAPACK library by looking at the default locations and names. In case, the installed ScaLAPACK library may have a different name and path in your system, you may need to explicitly specify the correct name of ScaLAPACK library and path. This can be done using the following flags:

```
cmake -DWITH_MPI=TRUE -DWITH_ACETB=TRUE -DLAPACK_LIBRARY_DIRS=/usr/local/opt/openblas/lib -DSCALAPACK_LIBRARIES="scalapack" -DSCALAPACK_LIBRARY="scalapack" -DSCALAPACK_LIBRARY_DIRS=/usr/local/opt/scalapack/lib ..
```
Here, ScaLAPACK library name is specified with `scalapack` and hence CMake try to find `libscalapack.so` in the directory that is given with `-DSCALAPACK_LIBRARY_DIRS`.

### Setting Python interface

As we did not specify a system-wide installation path above, DFTB+ and its Python iterface `dptools` will only be installed in the `build/install` directory where `build` is the dorectory we run `cmake`. 

To add the path of the installed `dptools` Python module to the Python path, you need to use the following command:

```
export PYTHONPATH=$PYTHONPATH:$HOME/dftbplus/build/install/lib/python3.9/site-packages
```

The above command assumes that you download DFTB+ repository (`dftbplus`) into your home directory and `cmake` found Python 3.9 in your system. To make sure the path is added to Python module paths in all shell instances, add the above line to your `~/.bashrc` and `~/.bash_profile` (for macOS, `~/.profile`) with

```
cat "export PYTHONPATH=$PYTHONPATH:$HOME/dftbplus/build/install/lib/python3.9/site-packages" >> ~/.bashrc
cat "export PYTHONPATH=$PYTHONPATH:$HOME/dftbplus/build/install/lib/python3.9/site-packages" >> ~/.bash_profile
source ~/.bash_profile
```

### Adding DFTB+ to System Path

Running the following lines at command line will add the `bin` directory of installation to the system-wide paths. Once the commands below are issued, executables of DFTB+ such as `dftb+` and `dp_bands` can be run directly at any path.

```
cat "export PATH=$PATH:$HOME/dftbplus/build/install/bin/" >> ~/.bashrc
cat "export PATH=$PATH:$HOME/dftbplus/build/install/bin/" >> ~/.bash_profile
source ~/.bash_profile
```
 
### Running the code (Serial)

Several test configurations are provided in `test` directory of `lib_acetb` distribution. To test the code, enter to `test` directory where `dftb_in.hsd` file resides. To run the code serial, use the following command:

```
dftb+
```

### Running the code with Julia Threads (Parallel)

ACEtb benefits from Julia threads by distributing per-atom calculations to each thread. 

Julia threads can be set by
```
export JULIA_NUM_THREADS=4
```

where `4` is the number of Julia threads that is used by each process. 

Note: ACEtb does not use OpenMP Threads and hence setting `OMP_NUM_THREADS` does not have any effect.

### Running the code with MPI (Parallel)

DFTB+ benefits from MPI paralleisation not only by distributing per-atom calculations to each process but also used BLACS routines at ScaLAPACK. 

Once DFTB+ is compiled with MPI support, you can run the code using the following command:
```
mpirun -np 4 dftb+
```

Note: Keep in mind that setting Julia Threads and using MPI will provide a mixed parallisation scheme and can benefit to speed up in multiple nodes. Example setting:

- For 2 nodes with 4 cores:
```
export JULIA_NUM_THREADS=4
mpirun -np 2 dftb+
```

DFTB+ outputs these settings as following:
```
Julia threads:               4
MPI processes:               2
OpenMP threads:              1
```

See more information on parallelisation of DFTB+ at [Parallel usage of DFTB+](https://dftbplus-recipes.readthedocs.io/en/latest/parallel/compiling.html)
