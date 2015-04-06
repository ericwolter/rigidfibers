# GPU Simulation of Rigid Fibers

This is the code and document for my masters's thesis. It includes a GPU-acclerated implementation of the numerical method "A numerical method for simulations of rigid fiber suspensions" by Anna Karin Tornberg and Katarina Gustavsson. I reduced the required runtime of 24h on a CPU to just 30min on a GPU.

## Prerequisites

- CUDA 6.5 (or later)
- OpenBLAS 0.2.12
- MAGMA 1.6.0
- cmake 2.8.7

## Building

The code is highly optimized thus it will be compiled each time a new simulation is executed. This allows for simulation variables to be compiled into the binary enabling performance increases especially with regard to fixed-size loops etc.

If all the prerequisites are correctly installed the python script *fibers_runner.py* will take care of compiling and running the simulation. Depending on the system configuration it might be necessary to explicitly set the `OPENBLAS_ROOT` environment variable to allow cmake to find the OpenBLAS libary.

  export OPENBLAS_ROOT=/my/custom/path/OpenBLAS

## Running

Running a simulation is done by simply executing *fibers_runner.py* which will take care of compiling and running everything. The script allows to execute both the CPU-based Fortran implementation as well as the GPU-based CUDA implementation if the system supports it.

  ./fibers_runner.py fortran run
                          (--direct | --gmres)
                          (--numerical | --analytical)
                          configuration_file fibers_file
  ./fibers_runner.py fortran validate
                          (--direct | --gmres)
                          (--numerical | --analytical)
                          configuration_file fibers_file
  ./fibers_runner.py fortran benchmark [--max_rse MAX_RSE]
                          (--magma | --gmres | --bicgstab)
                          (--numerical | --analytical)
                          configuration_file

  ./fibers_runner.py cuda run
                          (--magma | --gmres | --bicgstab)
                          (--numerical | --analytical)
                          (--D1 | --D2 | --D3)
                          configuration_file fibers_file
  ./fibers_runner.py cuda validate
                          (--magma | --gmres | --bicgstab)
                          (--numerical | --analytical)
                          (--D1 | --D2 | --D3)
                          configuration_file fibers_file
  ./fibers_runner.py cuda benchmark [--max_rse MAX_RSE]
                          (--magma | --gmres | --bicgstab)
                          (--numerical | --analytical)
                          (--D1 | --D2 | --D3)
                          configuration_file

Examples for the configuration_file and fibers_file format can be found in *tests/*.

## Tools

The directory *tools/* contains various auxiliary tools used for generating fibers files and validating results.

The most important tool is *gen2.f95* which is compiled using the supplied *build.sh*. This tools implements a fast method of generating fiber files with an enforced minimal distance between the fibers

  ./gen (sphere | box - not supported) NUMBER_OF_FIBERS MINIMAL_DISTANCE
