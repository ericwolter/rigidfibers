This documents describes a possible way to setup the prerequisites in order to generate simulation results.

# 1. CUDA

1. Download the latest drivers from NVIDIA for your platform
2. Verify that you can run the sample projects

# 2. OpenBLAS

1. Download source from [http://www.openblas.net/]()
2. Follow instructions on website under *Install from Source*
3. Compile source code

```bash
make
```

4. Install compiled library to user defined directory (i.e. /NOBACKUP/ewolter/openblas)

```bash
make PREFIX=my_openblas_directory install
```

# 3. MAGMA

1. Download source from [http://icl.cs.utk.edu/magma/software/]()
2. Follow instructions in file *README*
3. Rename openblas specific makefile

```bash
mv make.inc.openblas make.inc
```

4. Make sure **OPENBLASDIR** environment variable is set to the user defined directory

```bash
export OPENBLASDIR=my_openblas_directory
```

4. Compile source code

```bash
make
```

5. Install compiled library to user defined directory (i.e /NOBACKUP/ewolter/magma)

```bash
make install prefix=my_magma_directory
```

# 4. Build & Run rigidfibers

1. Make sure cmake will be able to find OpenBLAS by setting the environment variable **OPENBLAS_ROOT**
> This needs to be done for every session, so it might be a good idea to add this line to your bash_profile file.

```bash
export OPENBLAS_ROOT=my_openblas_directory
```

2. Execute *fibers_runner.py* script
> The script will take care of executing *cmake* and building the code correctly. You do not need to manually execute *cmake*

```bash
./fibers_runner.py cuda run --magma --numerical --D2 tests/config.ini tests/XcT_ref100.in
```
