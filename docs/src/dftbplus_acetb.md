
# DFTB+ Integration Using `lib_acetb`
`lib_acetb` is a module that integrates ACEtb.jl model into DFTB+ using C wrappers of Julia.

## Dependencies

Julia:
- ACEtb

Libraries:
- OpenBLAS
- MPI (Ex.: OpenMPI)

## Installation

```
# Before starting installation make sure you installed Julia >= 1.5
julia
] registry add https://github.com/JuliaMolSim/MolSim.git
] add https://github.com/ACEsuit/ACEtb.jl.git
# CTRL+D
git clone https://github.com/berkonat/lib_acetb
cd lib_acetb/dftbplus
./compile-MPI.sh
cd ../test/
# Serial run
../dftbplus/build/install/bin/dftb+
# Parallel run
mpirun -np 4 ../dftbplus/build/install/bin/dftb+
```

For MacOS, one can install the dependencies using `brew` with
```
brew install openblas
brew install openmpi
brew install scalapack
```

## Python installation of DFTB+

```
export PYTHONPATH=$PYTHONPATH:$HOME/lib_acetb/dftbplus/build/install/lib/python3.9/site-packages
```

## Band Structure Calculation

```
../dftbplus/build/install/bin/dp_bands band.out band_tot.dat
```

For more information on band structure calculation with DFTB+ and how to visulize the result, please check [dftbplus-recipes.readthedocs.io](https://dftbplus-recipes.readthedocs.io/en/latest/basics/bandstruct.html).

## Plotting

Use `Plot-Bands.ipynb` to plot Al FCC results and compare the output of DFTB+ from `band_tot.dat` file.

## Systematic fits

To run fits, run the script `./script_fits.sh` and adjust the parameters that you want to test in the script.
Check also that in the file dftb_in.hsd the following line matches with the file names in `./script_fits.sh`:
```
    Suffix = ".acetb.json"
```

