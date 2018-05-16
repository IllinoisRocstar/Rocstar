/* *******************************************************************
 * Illinois Open Source License                                      *
 *                                                                   *
 * University of Illinois/NCSA                                       * 
 * Open Source License                                               *
 *                                                                   *
 * Copyright@2008, University of Illinois.  All rights reserved.     *
 *                                                                   *
 *  Developed by:                                                    *
 *                                                                   *
 *     Center for Simulation of Advanced Rockets                     *
 *                                                                   *
 *     University of Illinois                                        *
 *                                                                   *
 *     www.csar.uiuc.edu                                             *
 *                                                                   *
 * Permission is hereby granted, free of charge, to any person       *
 * obtaining a copy of this software and associated documentation    *
 * files (the "Software"), to deal with the Software without         *
 * restriction, including without limitation the rights to use,      *
 * copy, modify, merge, publish, distribute, sublicense, and/or      *
 * sell copies of the Software, and to permit persons to whom the    *
 * Software is furnished to do so, subject to the following          *
 * conditions:                                                       *
 *                                                                   *
 *                                                                   *
 * @ Redistributions of source code must retain the above copyright  * 
 *   notice, this list of conditions and the following disclaimers.  *
 *                                                                   * 
 * @ Redistributions in binary form must reproduce the above         *
 *   copyright notice, this list of conditions and the following     *
 *   disclaimers in the documentation and/or other materials         *
 *   provided with the distribution.                                 *
 *                                                                   *
 * @ Neither the names of the Center for Simulation of Advanced      *
 *   Rockets, the University of Illinois, nor the names of its       *
 *   contributors may be used to endorse or promote products derived * 
 *   from this Software without specific prior written permission.   *
 *                                                                   *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************
 * Please acknowledge The University of Illinois Center for          *
 * Simulation of Advanced Rockets in works and publications          *
 * resulting from this software or its derivatives.                  *
 *********************************************************************/
/*
Reads a Gridgen .grd, .grds, or .grdb input file
and passes it to a blockConsumer.  Returns an 
error string on errors; returns NULL on sucess.

This file is written by Pointwise(tm) Gridgen (ver 15) as:
From the top level menu:
  Input/Output (e) ->
  Grid Pts Export (8) ->
  Block Volumes (4) ->
  Name (k) ->
  (type in name, ending with ".grd", ".grds", or ".grdd") ->
  style == PLOT3D (2), ->
  
For .grd:
  format == ascii (a), precision == double (d), 
For .grds:
  format == binary (^b), precision == single (s), 
For .grdd:
  format == binary (^b), precision == double (d), 

  Done (ent) ->
  Pick All (a), Done (ent)

Orion Sky Lawlor, olawlor@acm.org, 2004/2/18
*/
#include <stdio.h>
#include "gridutil.h"

/* Data types in Gridgen output file */
typedef int grdb_int_t;
typedef double grdb_coord_t;

/* Read an integer from this native-endian binary file */
static const char *readInt(FILE *f,int *dest) {
	grdb_int_t g;
	if (1!=fread(&g,sizeof(g),1,f))
		return "Error reading integer from .grdb file";
	*dest=g;
	return NULL;
}

/* Read a double from this native-endian binary file */
static const char *readDouble(FILE *f,double *dest) {
	grdb_coord_t g;
	if (1!=fread(&g,sizeof(g),1,f))
		return "Error reading coordinate from .grdb file";
	*dest=g;
	return NULL;
}

/// Abstracts the file format of a gridgen .grd file.
class gridgenGrdFormatter {
protected:
	FILE *f;
public:
	void setFile(FILE *f_) {f=f_;}
	virtual const char *readInt(int *dest) =0;
	virtual const char *readDouble(double *dest) =0;
};

/// ASCII format input file:
class gridgenAsciiFormatter : public gridgenGrdFormatter {
public:
	virtual const char *readInt(int *dest) {
		if (1!=fscanf(f,"%d",dest))
			return "Could not read integer from .grd file\n";
		else return NULL;
	}
	virtual const char *readDouble(double *dest) {
		if (1!=fscanf(f,"%lf",dest))
			return "Could not read coordinate from .grd file\n";
		else return NULL;
	}
};

/// Machine binary input file:
template <class coordT>
class gridgenBinaryFormatter : public gridgenGrdFormatter {
public:
	virtual const char *readInt(int *dest) {
		if (1!=fread(dest,sizeof(*dest),1,f))
			return "Could not read integer from binary .grd[sd] file\n";
		else return NULL;
	}
	virtual const char *readDouble(double *dest) {
		coordT in;
		if (1!=fread(&in,sizeof(in),1,f))
			return "Could not read coordinate from binary .grd[sd] file\n";
		*dest=in;
		return NULL;
	}
};

// Work around template instantiation problems on, e.g., Sun CC.
class gridgenBinaryFloatFormatter : public gridgenBinaryFormatter<float> {};
class gridgenBinaryDoubleFormatter : public gridgenBinaryFormatter<double> {};


const char *read_general(const char *gridFile,gridgenGrdFormatter &fmt,blockConsumer &dest)
{
	const char *err;
	//Open the input file
	FILE *f=fopen(gridFile,"r");
	if (f==NULL) return "Couldn't open input .grd file";
	fmt.setFile(f);

	//Read the header
	int nBlocks=0;
	if (NULL!=(err=fmt.readInt(&nBlocks))) return err;

	//Read the block sizes
	blockDim *dims=new blockDim[nBlocks];
	int b;
	for (b=0;b<nBlocks;b++) { 
		int nx,ny,nz;
		if (NULL!=(err=fmt.readInt(&nx))) return err;
		if (NULL!=(err=fmt.readInt(&ny))) return err;
		if (NULL!=(err=fmt.readInt(&nz))) return err;
		dims[b]=blockDim(nx,ny,nz);
	}
	
	//Read and consume each block's coordinates
	for (b=0;b<nBlocks;b++) {
		int nLocs=dims[b].getSize();
		vector3d *locs=dest.allocateBlock(dims[b]);
		for (int c=0;c<3;c++) //*outer* loop is coordinate axis
		  for (int i=0;i<nLocs;i++) {//*inner* loop is location #
			double loc;
			if (NULL!=(err=fmt.readDouble(&loc))) return err;
			locs[i][c]=loc;
		  }
		if (NULL!=(err=dest.consume(dims[b],locs))) return err;
		dest.freeBlock(locs);
	}		
	
	//Finish up
	delete[] dims;
	fclose(f);
	return NULL; //Everything worked
}


const char *read_grd(const char *gridFile,blockConsumer &dest)
{
	gridgenAsciiFormatter fmt;
	return read_general(gridFile,fmt,dest);
}
const char *read_grds(const char *gridFile,blockConsumer &dest)
{
	gridgenBinaryFloatFormatter fmt;
	return read_general(gridFile,fmt,dest);
}
const char *read_grdd(const char *gridFile,blockConsumer &dest)
{
	gridgenBinaryDoubleFormatter fmt;
	return read_general(gridFile,fmt,dest);
}

