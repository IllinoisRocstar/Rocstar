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
Reads a CSAR rocket .hdf input file, an HDF4 format, 
and passes it to a blockConsumer.  Returns an 
error string on errors; returns NULL on sucess.

Orion Sky Lawlor, olawlor@acm.org, 6/8/2001
*/
#include <stdio.h>
#include "gridutil.h"

/*HDF4 include file*/
#include "mfhdf.h"

#define chk(x) checkHDFerr(x,__FILE__,__LINE__)
static int checkHDFerr(int errCode,const char *fName,int lineNo) {
	if (errCode==-1) {
		fprintf(stderr,"HDF I/O Error in (%s:%d)\n",
			fName,lineNo);
		exit(23);
	}
	return errCode;
}

const char *read_hdf(const char *filename,blockConsumer &dest)
{
	int32 i, d;
	//Open the input file
	int32 sd_id = SDstart(filename, DFACC_READ);
	if (sd_id==-1) return "Couldn't open HDF input file!";

	//Select the X, Y, and Z components
	int32 sds_id[3];
	for (i=0;i<3;i++) sds_id[i]=chk(SDselect(sd_id,5+4*i));
	
	//Determine the ni, nj, and nk dimensions
	blockDim dim;
        /* This doesn't work.
	for (d=0;d<3;d++) {
		int32 dim_id = chk(SDgetdimid(sds_id[0],d));
		char dim_name[200];
		int32 dim_size, data_type, n_attrs;
		chk(SDdiminfo(dim_id, dim_name, &dim_size, &data_type, &n_attrs));
		dim[2-d]=dim_size;
	}
	*/
	char name[200];
	int32 rank, size[3], data_type, n_attrs;
	chk(SDgetinfo(sds_id[0], name, &rank, size, &data_type, &n_attrs));
	dim[0]=size[2];
	dim[1]=size[1];
	dim[2]=size[0];
	
	//Allocate storage for the X, Y, and Z output
	vector3d *locs=dest.allocateBlock(dim);

	//Read each of X, Y, and Z
	int nLocs=dim.getSize();
	double *coord=new double[nLocs];
	for (i=0;i<3;i++)
	{
		int32 start[3]={0,0,0};
		int32 edges[3];
		for (d=0;d<3;d++) edges[d]=dim[2-d];
		chk(SDreaddata(sds_id[i], start, NULL, edges, coord));
		//Copy this coordinate into locs
		for (int l=0;l<nLocs;l++) locs[l][i]=coord[l];
	}
	delete[] coord;

	//Pass locations to consumer
	dest.consume(dim,locs);
	dest.freeBlock(locs);
	
	//Close the file and return
	for (i=0;i<3;i++) chk(SDendaccess(sds_id[i]));
	chk(SDend(sd_id));
	return NULL;
}


