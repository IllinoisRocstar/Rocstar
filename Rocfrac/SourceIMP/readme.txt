Implicit Rocfrac and RocTherm Readme

Please note that this implicit version of
Rocfrac is currently only set up to work
on the Turing cluster.  The LAPACK, BLAS,
and BlockSolve95 libraries included may
have to be recompiled for it to work on
different machines.

To use the implicit version of Rocfrac and RocTherm,
make a new directory in rocstar/Codes/Rocfrac named
SourceIMP and SourceTherm.  Copy all of the files from the
SourceIMP and SourceTherm folders here to that folder.  
Then copy all of the files in this Source directory
into /rocstar/Codes/Rocfrac/Source and copy 
the files from this utilities/RocfracPrep 
folder to your 
rocstar/Codes/Rocfrac/utilities/RocfracPrep
folder.  Finally, copy Makefile.basic over 
your current rocstar/Codes/Rocfrac/Makefile.basic.

This process can be expedited by entering the
following command while in this directory:
cp -r * ../.

