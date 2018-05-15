Rocstar
-----

Rocstar multiphysics simulation application (RocstarMP)

Rocstar is a multiphysics simulation application designed to do fluid-structure interaction (FSI) across moving, reacting interfaces. Rocstar couples multiple domain-specific simulation packages and disparately discretized domains and provides several simulation-supporting services including conservative and accurate data transfer, surface propagation, and parallel I/O. Rocstar is MPI parallel and routinely executes large simulations on massively parallel platforms. Rocstar was originally developed at the University of Illinois Center for Simulation of Advanced Rockets (CSAR) under Department of Energy ASCI program funding. Ongoing development of Rocstar is conducted by Illinois Rocstar LLC with company IR&D and continued DOE SBIR funding.


## Getting Started ##
To acquire RocstarMP, you can download it from Illinois Rocstar's GitHub
or clone it with the following command:
```
$ git clone git@git.illinois.rocstar:rocstar_modern/rocstar.git
```
## Build Instructions ##
### Build Dependencies ###
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

* IMPACT (available at illinoir rocstar git repo)
* METIS 4.0.3 (available http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD)

## Build RocstarMP ##
For the following steps we assume $ROCSTAR_PROJECT_PATH is the path to RocstarMP, $ROCSTAR_INSTALL_PATH is 
the desired installation location, $IMPACT_INSTALL_PATH is the path to the locatin IMPACT project is installed, 
$METIS_LIB_PATH and $METIS_INC_PATH are the paths to the location of library and headers of metis.
Star the build process by executing:

```
$ cd $ROCSTAR_PROJECT_PATH
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=$ROCSTAR_INSTALL_PATH -DENABLE_MPI=ON -DENABLE_CGNS=ON -DIMPACT_HDR=$IMPACT_INSTALL_PATH/include/com.h -DMETIS_INC=$METIS_INC_PATH -DMETIS_LIB=$METIS_LIB_PATH -DSITCOM_LIB=$IMPACT_INSTALL_PATH/lib/libSITCOM.so .. 
$ make -j
$ make install
```

Executing the commands above will build all libraries and executables.

## Testing RocstarMP ##
NOTE: Disregard the testing for now, currently heavily depend on IRAD project. Will be changed as we go.
From the build directory, execute the following command to test the installation:
```
$ make test
```
This will execute several tests in `$ROCSTAR_PROJECT_PATH/testing`. See the testing directories here for more details.

### Manually Build Third Party Libraries ###

#### Building IMPACT ####
A special version of IMPACT (free from IRAD) is used. The project should be obtian from IR's internal git repository.
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



