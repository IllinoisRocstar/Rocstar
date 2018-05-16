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
Main routine for grid split-- read the input 
mesh file, write out the .mblock files.

Orion Sky Lawlor, olawlor@acm.org, 6/7/2001
*/
#include <stdio.h>
#include "makeflo.h"

void checkError(const char *errCode) {
  if (errCode==NULL) return; //No problem
  fprintf(stderr,"FATAL ERROR! %s\n",errCode);
  exit(1);
}

void printUsage(const char *why)
{
	printf("Usage: makemblock \n"
	  "     <.grd, or .msh mesh; and .inp and .bc file> \n"
	  "     <# output chunks> <output file prefix>\n"
	  "  Makemblock reads a mesh (first parameter) and boundary condition\n"
	  "(.inp) list, partitions the mesh into the requested number of\n"
	  "chunks, and writes out mblock input files.\n"
	  "  The mesh formats supported are:\n"
	  "<*.grd>  A Gridgen ASCII double-precision mesh description\n"
	  "<*.msh>  A double-precision binary Rocflo mesh file\n"
	  "<*.hdf>  A double-precision 3D HDF mesh file\n"
	  "<*.mblk> A double-precision 3D Mblock mesh file\n"
	  "  When the mesh consists of multiple files, give only the name of\n"
	  "the first file (e.g., 'tstflo_001.hdf').  The number of blocks and\n"
	  "numeric format will be automatically determined.\n"
	  "Part of the Charm++ Tools.  Version " MAKEFLO_VERSION "\n"
	); 
	if (why!=NULL)
		printf("Exiting> %s\n",why);
	exit(1);
}

makefloParam parameters;

int main(int argc,char *argv[]) 
{
//Parse the command line arguments
	if (argc<4) printUsage("Not enough command-line arguments");
	int curArg=1;

	while (argv[curArg][0]=='-') {
		if (0==strcmp(argv[curArg],"-2D")) {
			parameters.skipAxis=2; //Don't bother about the z axis.
			curArg++;
		}
		else if (0==strcmp(argv[curArg],"-top")) {
			parameters.topologyOnly=1; //Only write .bblk
			curArg++;
		}
		else
			printUsage("Unrecongized parameter");
	}

	if (curArg+3!=argc) printUsage("Too many arguments");
	const char *inMesh=argv[curArg+0];
	string inpFile=replaceExtention(argv[curArg+0],".inp");
	int nPieces=atoi(argv[curArg+1]);
	const char *outMesh=argv[curArg+2];
	
	vector<block *> blocks;
	
//Read the blocks of the mesh
	{ //<- scoping block so r's destructor gets called
		blockReader r(blocks);
		checkError(r.read(inMesh));
		printf("Mesh file read successfully.\n");
	}
	
//Read the boundary conditions
	checkError(readBoundaries(blocks,inpFile.c_str()));
	printf("Boundary conditions read successfully.\n");
	
//Split the blocks (for greater parallelism)
	checkError(splitBlocks(blocks,nPieces));
	printf("Split blocks for %d PEs (%d blocks generated).\n",nPieces,blocks.size());	

//Build the block faces & associate all shared nodes
	buildFaces(blocks,false);
	printf("Nodes matched successfully.\n");
	
//Write out the blocks
	checkError(writeMblock(blocks,outMesh));
	printf("Block files written successfully\n");
	return 0;
}







