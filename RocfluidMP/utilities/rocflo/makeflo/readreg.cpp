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
Reads a proprietary .reg input file, a simple ASCII format, 
and passes it to a blockConsumer.  Returns an 
error string on errors; returns NULL on sucess.

Orion Sky Lawlor, olawlor@acm.org, 4/29/2002
*/
#include <stdio.h>
#include "gridutil.h"

const char *read_reg(const char *regFile,blockConsumer &dest)
{
	//Open the input file
	FILE *in=fopen(regFile,"r");
	if (in==NULL) return "Couldn't open input .reg file";

	//Read the block size
	int dimArr[3];
	vector3d min,max;
	for (int axis=0;axis<3;axis++)
	{
		double m;
		if (1!=fscanf(in,"%d",&dimArr[axis]))
			return "Couldn't read block dimension";	
		if (1!=fscanf(in,"%lf",&m)) return "Couldn't read block min";
		min[axis]=m;
		if (1!=fscanf(in,"%lf",&m)) return "Couldn't read block max";
		max[axis]=m;
	}
	fclose(in);

	//Allocate return grid
	blockDim dim(dimArr[0],dimArr[1],dimArr[2]);
	vector3d *locs=dest.allocateBlock(dim);
	
	//Fill in the grid
	blockLoc i;
	blockSpan s(blockLoc(0,0,0),dim);
	BLOCKSPAN_FOR(i,s) {
		vector3d v;
		for (int axis=0;axis<3;axis++) {
			double scl=i[axis]/(float)(dim[axis]-1);
			v[axis]=min[axis]+scl*(max[axis]-min[axis]);
		}
		locs[dim[i]]=v;
	}
	
	dest.consume(dim,locs);	
	dest.freeBlock(locs);

	return NULL; //Everything worked
}

