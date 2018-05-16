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
Create a single, extremely simple block-- a regular,
axis-aligned box.

Orion Sky Lawlor, olawlor@acm.org, 7/18/2001
*/
#include <stdio.h>
#include "adj.h"
#include "util.h"
#include "makeflo.h"

void printUsage(const char *why)
{
	printf("Usage: simpleblock <output .msh or .hdf file>\n"
	  "     <X start> <X end> <X steps> \n"
	  "     <Y start> <Y end> <Y steps> \n"
	  "     <Z start> <Z end> <Z steps> \n"
	  "  Write out a single axis-aligned block with the given dimensions\n"
	  "and coordinate extents. \n"
	); 
	if (why!=NULL)
		printf("Exiting> %s\n",why);
	exit(1);
}

makefloParam parameters;

int main(int argc,char *argv[]) 
{
//Parse the command line arguments
	if (argc!=11) printUsage("Not enough command-line arguments");
	int curArg=1;
	const char *outMesh=argv[curArg++];
	blockDim dim;
	vector3d start,end;
	int axis;
	static const char *axisNames[]={"X","Y","Z"};
	for (axis=0;axis<3;axis++) {
		double s,e;
		if (1!=sscanf(argv[curArg++],"%lg",&s))
			printUsage("Couldn't parse axis start");
		if (1!=sscanf(argv[curArg++],"%lg",&e))
			printUsage("Couldn't parse axis end");
		start[axis]=s; end[axis]=e;
		if (1!=sscanf(argv[curArg++],"%d",&dim[axis]))
			printUsage("Couldn't parse axis size");
		printf("Along %s axis: (%f - %f) in %d steps\n",axisNames[axis],
			start[axis],end[axis],dim[axis]);
	}
	
//Build a regular 3D block of points
	vector3d *locs=new vector3d[dim.getSize()];
	blockLoc i;
	for (i[2]=0; i[2]<dim[2]; i[2]++)
	for (i[1]=0; i[1]<dim[1]; i[1]++)
	for (i[0]=0; i[0]<dim[0]; i[0]++) { 
		vector3d p;
		for (axis=0;axis<3;axis++)
			p[axis]=start[axis]+((end[axis]-start[axis])*i[axis])/(dim[axis]-1);
		locs[dim[i]]=p;
	}
	printf("Locations built\n");

	vector<block *> blocks;
	blocks.push_back(new block(dim,0,0,locs));
	printf("Block built\n");
	
//Write out the new block
	writeBlocks(blocks,outMesh);
	printf("Program finished successfully\n");
	return 0;
}







