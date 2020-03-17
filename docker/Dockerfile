FROM centos:7

RUN yum update -y
RUN yum upgrade -y

# Compilers and build tools.
RUN yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm -y
RUN yum groupinstall "Development Tools" -y
RUN yum install cmake3 gcc-c++.x86_64 make mpich-3.2-devel.x86_64 cgnslib-devel hdf5-devel hdf-devel.x86_64 yum install blas-devel lapack-devel environment-modules centos-release-scl scl-utils -y

# Change the default compilers to gcc7
RUN yum install devtoolset-7 -y; 
RUN yum install vim -y

# Finds the module command, loads the mpich module
# RUN echo "source /etc/profile.d/modules.sh" >> /etc/bashrc
# RUN echo "module load mpi/mpich-3.2-x86_64" >> /etc/bashrc
# RUN echo "source scl_source enable devtoolset-7" >> /etc/bashrc
# RUN echo "scl enable devtoolset-7 bash" >> /etc/bashrc
RUN echo "source /opt/rh/devtoolset-7/enable" >> /etc/bashrc
RUN echo "module load mpi/mpich-3.2-x86_64" >> /etc/bashrc
ENV BASH_ENV "/etc/bashrc"

# Copy the TPLs ( metis-4 & IMPACT )
COPY . /home/

WORKDIR /rocstar/
ENTRYPOINT /bin/bash
