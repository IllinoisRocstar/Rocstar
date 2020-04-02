Rocstar
-----

**Rocstar** **M**ulti**P**hysics simulation suite (RocstarMP)

Rocstar is a multiphysics simulation application designed for coupled multiphysics simulations involving fluid-structure interaction (FSI) across moving, reacting interfaces. Rocstar couples multiple domain-specific simulation packages and disparately discretized domains and provides several simulation-supporting services including conservative and accurate data transfer, surface propagation, and parallel I/O. Rocstar is MPI-parallel. Rocstar was originally developed at the University of Illinois Center for Simulation of Advanced Rockets (CSAR) under Department of Energy ASCI program funding.

## Version ##
Version 5.0.0

Rocstar follows semantic versioning. The version will be major.minor.patch.
We will:
* Increase the version for bug fixes, security fixes, and code
documentation. Backwards compatible; no breaking changes.
* Increase the minor version for new features and additions to the library’s
interface. Backwards compatible; no breaking changes.
* Increase the major version for breaking changes to the library’s interface or
breaking changes to behavior.

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
* EITHER mpich OR openmpi
* liblapack-dev
* libblas-dev
* libjpeg-dev

Make sure to compile and install following packages before you start compiling RocstarMP:

* IMPACT 2.1.1 (available [here](https://github.com/IllinoisRocstar/IMPACT))
* METIS 4.0.3 (available [here](http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD))

Follow standard compile instructions for these dependencies. Note that RocstarMP uses an older version of METIS.

### Build RocstarMP ###
**NOTE** Currently RocstarMP is only tested with MPI compilers. Serial build is not recommended.

For the following steps we assume `$ROCSTAR_PROJECT_PATH` is the path to RocstarMP, `$ROCSTAR_INSTALL_PATH` is the desired installation location, `$IMPACT_DIR` is the path to the location the IMPACT project is installed, and `$METIS_LIB` is the location of library of metis.
Start the build process by:

```
$ cd $ROCSTAR_PROJECT_PATH
$ mkdir build && cd build
$ IMPACT_DIR=$IMPACT_DIR cmake -DCMAKE_INSTALL_PREFIX=$ROCSTAR_INSTALL_PATH -DMETIS_LIB=$METIS_LIB .. 
$ make -j$(nproc)
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

### Post-install ###
Add the RocstarMP library directory `$ROCSTAR_INSTALL_PATH/lib` to your environment `$LD_LIBRARY_PATH`.

### Manually Build Third Party Libraries ###

#### Building IMPACT ####
The new version of IMPACT project should be obtianed from GitHub [here](https://github.com/IllinoisRocstar/IMPACT).
Follow directions given at the project page to compile and build IMPACT.

## Run environment for RocstarMP ##
The executables in `bin` require the libraries in `lib` to be available in the runtime environment. When a non-standard installation location is used, this can be done in the shell start-up script or before executing a RocstarMP binary, for example:
```
$ cd $ROCSTAR_INSTALL_PATH/bin
$ LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH ./rocstar
```
