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
Determines the format of, and calls the appropriate
routine to read a CSAR fluid block mesh.

Orion Sky Lawlor, olawlor@acm.org, 6/8/2001
*/
#include <cstring>
#include <cstdio>
#include <strings.h>
#include "gridutil.h"
#include "util.h"
//#include "charm-api.h"
#include "FC.h"

static blockConsumer *curConsumer=NULL;

//Fortran numeric data type:
typedef double REAL;

//Fortran-callable mesh input routine
extern "C" void FC_GLOBAL(accept_locations,ACCEPT_LOCATIONS)
	(int *ni,int *nj,int *nk,
	REAL *x,REAL *y,REAL *z)
{
	blockDim dim(*ni,*nj,*nk);
	vector3d *v=curConsumer->allocateBlock(dim);
	blockLoc i;
	for (i[2]=0; i[2]<dim[2]; i[2]++)
	for (i[1]=0; i[1]<dim[1]; i[1]++)
	for (i[0]=0; i[0]<dim[0]; i[0]++) 
	{ //Copy data from Fortran array
		int c_i=dim[i], f_i=dim[i];
		v[c_i].x=x[f_i];
		v[c_i].y=y[f_i];
		v[c_i].z=z[f_i];
	}
	curConsumer->consume(dim,v);
	curConsumer->freeBlock(v);	
}

//Prototypes for read functions:
extern "C" void FC_GLOBAL(read_msh,READ_MSH)(const char *fName);
const char *read_grd(const char *gridFile,blockConsumer &dest);
const char *read_grds(const char *gridFile,blockConsumer &dest);
const char *read_grdd(const char *gridFile,blockConsumer &dest);
const char *read_reg(const char *regFile,blockConsumer &dest);
const char *read_hdf(const char *hdfName,blockConsumer &dest);

//Read blocks from the given file into the given consumer
const char *read_file(const char *file,blockConsumer &dest)
{
	curConsumer=&dest;
	if (endsWith(file,".msh"))
          FC_GLOBAL(read_msh,READ_MSH)(fortranifyString(file));
	else if (endsWith(file,".grd"))
		return read_grd(file,dest);
	else if (endsWith(file,".grds"))
		return read_grds(file,dest);
	else if (endsWith(file,".grdd"))
		return read_grdd(file,dest);
#ifdef USE_HDF
	else if (endsWith(file,".hdf"))
		return read_hdf(file,dest);
#endif
	else if (endsWith(file,".reg"))
		return read_reg(file,dest);
	else
		return "Unrecognized input mesh file extension! ("__FILE__")";
	return NULL;
}

//Read blocks from the given file specification,
//  which may contain several blocks.
const char *read_multiple(const char *file,blockConsumer &dest)
{
	char tmp[1024];
	strcpy(tmp,file);
	do {
		printf("Reading block file '%s'...\n",tmp);
		const char *ret=read_file(tmp,dest);
		if (ret!=NULL) return ret;
		if (!incrementAscii(tmp)) return NULL;
	} while (fileExists(tmp));
	return NULL;
}

/******* BlockConsumer implementation *********/
const char *blockConsumer::read(const char *file)
{
	return read_multiple(file,*this);
}

//Location vector block allocation/deallocation
vector3d *blockConsumer::allocateBlock(blockDim &dim)
{
	return new vector3d[dim.getSize()];
}
void blockConsumer::freeBlock(vector3d *blk)
{
	delete[] blk;
}

blockConsumer::~blockConsumer() { }


	
#ifdef STANDALONE
/*Unit test driver program:
	Prints out a debugging version of the input file.
*/

int printAll=0;

class debugBlockConsumer : public blockConsumer{
	int count;//Number of blocks seen so far

	void print(const vector3d &v,const char *post="") {
		printf("%g\t%g\t%g%s",
			v.x,v.y,v.z,post);
	}

	void checkPoint(int i,int j,int k,
		const blockDim &dim,const vector3d *locs)
	{
		const vector3d &o=locs[dim[blockLoc(i,j,k)]];
		const vector3d &x=locs[dim[blockLoc(i+1,j,k)]];
		const vector3d &y=locs[dim[blockLoc(i,j+1,k)]];
		const vector3d &z=locs[dim[blockLoc(i,j,k+1)]];
		vector3d del(o.dist(x),o.dist(y),o.dist(z));
		printf("\tGrid spacings:");print(del,"\n");
		vector3d nx=(x-o).dir();
		vector3d ny=(y-o).dir();
		vector3d nz=(z-o).dir();
		printf("\tSlant-ness:xy=%.3f\tyz=%.3f\txz=%.3f\n",
			nx.dot(ny),ny.dot(nz),nx.dot(nz));
	}

public:
	debugBlockConsumer() {count=0;}
	
	virtual const char *consume(
		const blockDim &dim,//Dimentions of incoming block
		vector3d *locs) //X,Y,Z coordinates (ni x nj x nk)
	{
		printf("Block #%d:  (%d x %d x %d)\n",++count,
			dim[0],dim[1],dim[2]);
		//Check the inter-point spacing along each axis
		checkPoint(2,3,1, //Check some random point
			dim,locs);
		checkPoint(0,1,2, //Check another random point
			dim,locs);
		
		//Print out a few points
		int l;
		int max=10;
		if (max>dim.getSize() || printAll) max=dim.getSize();
		for (l=0;l<max;l++) {
			printf("\t[%d]=",l);
			print(locs[l],"\n");
		}
		return NULL;
	}
};

int main(int argc,char *argv[]) {
	if (argc<2) {printf("Usage: read <mesh file> [all]\n"); return 1;}
	debugBlockConsumer dest;
	const char *err;
	if (argc>2) printAll=1;
	if (NULL!=(err=dest.read(argv[1])))
		printf("ERROR! %s\n",err);
	else
		printf("Finished successfully.\n");
	return 0;
}
#endif



