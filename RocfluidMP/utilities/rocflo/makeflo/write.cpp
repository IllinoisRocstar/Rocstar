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
routine to write a CSAR fluid block mesh.

Orion Sky Lawlor, olawlor@acm.org, 6/13/2001
*/
#include <stdio.h>
#include <strings.h>
#include "adj.h"
#include "util.h"
//#include "charm-api.h"
#include "FC.h"

//Fortran numeric data type:
typedef double REAL;

extern "C" {
	typedef void (*fn_fortran_write)(int *block,int *ni,int *nj,int *nk,
		REAL *x,REAL *y,REAL *z,const char *fName); // (JB: added block)
}

//Prototypes for write functions:  (JB: added block)
extern "C" {
 void FC_GLOBAL(write_msh,WRITE_MSH)(int *block,int *ni,int *nj,int *nk,
	REAL *x,REAL *y,REAL *z,const char *fName);
 void FC_GLOBAL(write_grda,WRITE_GRDA)(int *block,int *ni,int *nj,int *nk,
	REAL *x,REAL *y,REAL *z,const char *fName);
 void FC_GLOBAL(write_grdb,WRITE_GRDB)(int *block,int *ni,int *nj,int *nk,
	REAL *x,REAL *y,REAL *z,const char *fName);
}

//Send this data to the given fortran routine (JB: added block)
static void write_fortran(int block,blockDim dim,const vector3d *locs,
	const char *fName,fn_fortran_write writeFn)
{
	//Copy the data from C array to fortran arrays
	int len=dim.getSize();
	REAL *x=new REAL[len];
	REAL *y=new REAL[len];
	REAL *z=new REAL[len];
	blockLoc i;
	for (i[2]=0; i[2]<dim[2]; i[2]++)
	for (i[1]=0; i[1]<dim[1]; i[1]++)
	for (i[0]=0; i[0]<dim[0]; i[0]++) 
	{ //Copy data into fortran array
		int c_i=dim[i], f_i=dim[i];
		x[f_i]=locs[c_i].x;
		y[f_i]=locs[c_i].y;
		z[f_i]=locs[c_i].z;
	}

	//Send the data off to fortran (JB: added block)
	(writeFn)(&block,&dim[0],&dim[1],&dim[2],x,y,z,fortranifyString(fName));
	delete[] x; delete[] y; delete[] z;
}

void vecMax(const vector3d &a,vector3d &dest)
{
	if (a.x>dest.x) dest.x=a.x;
	if (a.y>dest.y) dest.y=a.y;
	if (a.z>dest.z) dest.z=a.z;
}
void vecMin(const vector3d &a,vector3d &dest)
{
	if (a.x<dest.x) dest.x=a.x;
	if (a.y<dest.y) dest.y=a.y;
	if (a.z<dest.z) dest.z=a.z;
}

void write_bounds(const blockDim &dim,const vector3d *locs,
	const char *filename_base)
{
	char filename[1024];
	sprintf(filename,"%s.bounds",filename_base);
	FILE *f=fopen(filename,"w");
	if (f==NULL) return;
	//Find the bounds of all locations
	double big=1.0e25;
	vector3d max(-big,-big,-big),min(big,big,big);
	blockLoc i;
	for (i[2]=0; i[2]<dim[2]; i[2]++)
	for (i[1]=0; i[1]<dim[1]; i[1]++)
	for (i[0]=0; i[0]<dim[0]; i[0]++) {
		vecMax(locs[dim[i]],max);
		vecMin(locs[dim[i]],min);
	}
	fprintf(f,"%f %f %f %f %f %f\n",min.x,max.x,min.y,max.y,min.z,max.z);
	fclose(f);
}

const char *write_hdf(const blockDim &dim,const vector3d *locs,
	const char *filename);

//Write blocks to the given file
const char * writeBlocks(vector<block *> &blocks,
                      const char *outMesh)
{
	char oName[1024];
	bool isFortran=false, isHDF=false;
	bool incrementFilenames=true;
	fn_fortran_write writeFn=NULL;
	strcpy(oName,outMesh);
	if (endsWith(outMesh,".msh"))
	{ //Write .msh files:
		isFortran=true;
		writeFn=FC_GLOBAL(write_msh,WRITE_MSH);
	}
	else if (endsWith(outMesh,"1")||
	         endsWith(outMesh,".grda")||
	         endsWith(outMesh,".grdb")) 
	{ //Probably a grda or grdb file:
		isFortran=true;
		incrementFilenames=false; //Same filename every time
		char tmp[1024];
		strcpy(tmp,outMesh);
		//Clip off any trailing numbers:
		int endIdx=strlen(tmp)-1;
		while (endIdx>0&&(isdigit(tmp[endIdx])||(tmp[endIdx]=='_')))
			tmp[endIdx--]=0;
		if (endsWith(tmp,".grda"))
                  writeFn=FC_GLOBAL(write_grda,WRITE_GRDA);
		else if (endsWith(tmp,".grdb")) 
                  writeFn=FC_GLOBAL(write_grdb,WRITE_GRDB);
		else return "Unrecognized output mesh extension! (2)";
	}
#ifdef USE_HDF
	else if (endsWith(outMesh,".hdf"))
	{
		isHDF=true;
	}
#endif
	else if (endsWith(outMesh,".null"))
		{/*Just skip output*/}
	else
		return "Unrecognized output mesh file extension! (1)";
		
	
	for (unsigned int i=0;i<blocks.size();i++) {
		const blockDim &dim=blocks[i]->getDim();
		const vector3d *locs=&blocks[i]->getLoc(blockLoc(0,0,0));

		if (isFortran)
		  if (i+1==blocks.size())
		  { /* Last time around, pass a negative block number
		       so fortran code can close its file. */
		    write_fortran(-(i+1),dim,locs,oName,writeFn);
		  } else {
		    write_fortran(  i+1 ,dim,locs,oName,writeFn);
		  }
               
#ifdef USE_HDF
		else if (isHDF)
			write_hdf(dim,locs,oName);
#endif
		
		printf("."); fflush(stdout);
                if (incrementFilenames) {
		  if (!incrementAscii(oName)) 
			return "Cannot increment output filename!  Use something like 'foo_0001.msh'\n";
                }
	}
	printf("\n");

	return NULL;
}


