Rocstar
-----

**Rocstar** **M**ulti**P**hysics simulation suite (RocstarMP)

Rocstar is a multiphysics simulation application designed for coupled multiphysics simulations involving fluid-structure interaction (FSI) across moving, reacting interfaces. Rocstar couples multiple domain-specific simulation packages and disparately discretized domains and provides several simulation-supporting services including conservative and accurate data transfer, surface propagation, and parallel I/O. Rocstar is MPI-parallel. Rocstar was originally developed at the University of Illinois Center for Simulation of Advanced Rockets (CSAR) under Department of Energy ASCI program funding.


## Getting Started ##
To acquire RocstarMP, you can download it from Illinois Rocstar's GitHub
or clone it with the following command:
```
$ git clone https://github.com/IllinoisRocstar/Rocstar.git
```
## Build Instructions ##
### Dependencies ###
Make sure to `apt install` following before you start

* build-essential
* cmake
* mpich
* libcgns-dev
* libhdf5-dev
* liblapack-dev
* libblas-dev
* libjpeg-dev

Make sure to compile and install following packages before you start compiling RocstarMP:

* IMPACT (available [here](https://github.com/IllinoisRocstar/IMPACT))
* METIS 4.0.3 (available [here](http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD))

Follow standrad compile instructions for these dependencies. Note that RocstarMP uses an older version of METIS.

### Build RocstarMP ###
**NOTE** Currently RocstarMP is only tested with MPICH compilers. If you have both OpenMPI and MPICH installed make sure `mpicxx`, `mpicc`, and `mpif90` point to the MPICH. In the following we have assumed both MPI libraries are installed.

For the following steps we assume `$ROCSTAR_PROJECT_PATH` is the path to RocstarMP, `$ROCSTAR_INSTALL_PATH` is the desired installation location, `$IMPACT_INSTALL_PATH` is the path to the locatin IMPACT project is installed, `$METIS_LIB_PATH` and `$METIS_INC_PATH` are the paths to the location of library and headers of metis.
Start the build process by:

```
$ cd $ROCSTAR_PROJECT_PATH
$ mkdir build && cd build
$ cmake -DMPI_C_COMPILER=/usr/bin/mpicc.mpich -DMPI_CXX_COMPILER=/usr/bin/mpicxx.mpich -DMPI_Fortran_COMPILER=/usr/bin/mpif90.mpich -DCMAKE_C_COMPILER=mpicc.mpich -DCMAKE_CXX_COMPILER=mpicxx.mpich -DCMAKE_Fortran_COMPILER=mpif90.mpich -DCMAKE_INSTALL_PREFIX=$ROCSTAR_INSTALL_PATH -DENABLE_MPI=ON -DENABLE_CGNS=ON -DIMPACT_HDR=$IMPACT_INSTALL_PATH/include/com.h -DMETIS_INC=$METIS_INC_PATH -DMETIS_LIB=$METIS_LIB_PATH -DSITCOM_LIB=$IMPACT_INSTALL_PATH/lib/libSITCOM.so .. 
$ make -j6
$ make install
```
Executing the commands above will build all libraries and executables.

### Testing RocstarMP ###
**NOTE** The testing is currently heavily dependent on IR's IRAD project. Modern testing is currently under development. In case needed, IRAD testing framework can be revived temporarily. Please contact developers for further instructions.

To perform testing, execute the following in the build directory:
```
$ make test
```
The output of tests are captured in `$ROCSTAR_PROJECT_PATH/testing`. The testing framework also keeps a log of the test outputs in `$ROCSTAR_PROJECT_PATH/Testing` directory. If tests fail seek output log in this directory for more details.

### Manually Build Third Party Libraries ###

#### Building IMPACT ####
The new version of IMPACT project should be obtianed from GitHub [here](https://github.com/IllinoisRocstar/IMPACT).
Follow directions given at the project page to compile and build IMPACT.

#### Building CGNS ####
Obtain latest version of CGNS
```
git clone https://github.com/CGNS/CGNS.git
```
Execute following
```
$ cd CGNS
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=$ROCSTAR_PROJECT_PATH/install/cgns ..
$ make -j
$ make install
```



