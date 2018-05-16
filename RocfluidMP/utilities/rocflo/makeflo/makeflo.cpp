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
mesh file, write out the .flo file.

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
	printf("Usage: makeflo [ -splitaxis <-1, 0, 1, or 2)> ]\n"
	  "      [ -splitrcb ]      [ -splithalf ]\n"
	  "      [ -multilevel <maximum level> ]\n"
	  "     <.grd, .grdb, .msh, or .hdf; and .inp and .bc file> \n"
	  "     <# PEs> <.flo or .top output file> <.msh or .hdf output file>\n"
          "\n"
	  "  Makeflo reads a mesh (first parameter) and boundary condition\n"
	  "(.inp) list and description (.bc), partitions the mesh into the\n"
	  "requested number of processors, and writes out a rocflo control\n"
	  "file (.flo or .top) and the partitioned mesh (last parameter).\n"
	  "  The mesh formats supported are:\n"
	  "<*.grd>  A Gridgen ASCII double-precision mesh description\n"
	  "<*.msh>  A double-precision binary Rocflo mesh file\n"
	  "<*.hdf>  A double-precision 3D HDF mesh file\n"
	  "  When the mesh consists of multiple files, give only the name of\n"
	  "the first file (e.g., 'tstflo_001.hdf').  The number of blocks and\n"
	  "numeric format will be automatically determined.\n"
          "\n"
          "Options:\n"
          " -splitaxis <-1, 0, 1, or 2)>: axis to split along. \n"
          " -splitrcb:   split using recursive bisection.\n"
	  " -multilevel <maximum level>: number of multigrid levels to create \n"
          " -splithalf:  split every block in half, <# PEs> is ignored.\n"
          "\n"
	  "Part of the CSAR Tools.  Version " MAKEFLO_VERSION "\n"
	); 
	if (why!=NULL)
		printf("Exiting> %s\n",why);
	exit(1);
}

makefloParam parameters;

int main(int argc,char *argv[]) 
{
//Parse the command line arguments
	if (argc<5) printUsage("Not enough command-line arguments");
	int curArg=1;
	while (argv[curArg][0]=='-') 
	{//Parse a flag argument
		const char *flag=argv[curArg++];
		if (0==strcmp(flag,"-h") ||
		    0==strcmp(flag,"-?") ||
		    0==strcmp(flag,"--help")) 
			printUsage(NULL);
		else if (0==strcmp(flag,"-splitaxis")) {
			if (1!=sscanf(argv[curArg++],"%d",&parameters.splitAxis))
				printUsage("Must pass 0,1, or 2 to -splitaxis");
		}
		else if (0==strcmp(flag,"-splitrcb")) {
			parameters.splitRCB=1;
		}
		else if (0==strcmp(flag,"-splithalf")) {
			parameters.splithalf=1;
		}
		else if (0==strcmp(flag,"-multilevel")) {
			int level=0;
			if ((1!=sscanf(argv[curArg++],"%d",&level))||(level<=0))
				printUsage("Must pass a positive integer to -multilevel");
			parameters.setLevel(level);
		}		
		else printUsage("Unrecognized flag");
	}
	if (curArg+4!=argc) printUsage("Too many arguments");
	const char *inMesh=argv[curArg+0];
	string inpFile=replaceExtention(argv[curArg+0],".inp");
	int nPieces=atoi(argv[curArg+1]);
	const char *floFile=argv[curArg+2];
	const char *outMesh=argv[curArg+3];

	bool forFlo;
	const char *bcExt="?";
	if (endsWith(floFile,".flo")) {
		forFlo=true;
		bcExt=".bc";
	}
	else if (endsWith(floFile,".top")) {
		forFlo=false;
		bcExt=".bcmp";
	}
	else {
		printf("UNRECOGNIZED TOPOLOGY NAME: '%s'\n",floFile);
		exit(1);
	}
	string bcFile=replaceExtention(argv[curArg+0],bcExt);
	
	vector<block *> blocks;
	
//Read the blocks of the mesh
	{ //<- scoping block so r's destructor gets called
		blockReader r(blocks);
		checkError(r.read(inMesh));
		printf("Mesh file read successfully.\n");
	}

// Check the quality of the mesh
	for (unsigned int i=0;i<blocks.size();i++)
		checkQuality(*blocks[i]);
	
//Read the boundary conditions
	checkError(readBoundaries(blocks,inpFile.c_str()));
	printf("Boundary conditions read successfully.\n");
	
	if (parameters.splithalf) {
	  	// split in half for every block
		nPieces = blocks.size()*2;
		printf("Split every block in half for %d pieces.\n", nPieces);
  	}

//Split the blocks (for greater parallelism)
	checkError(splitBlocks(blocks,nPieces));
	printf("Split blocks for %d PEs (%d blocks generated).\n",
		nPieces,(int)blocks.size());	

//Build the block faces & associate all shared nodes
	buildFaces(blocks,forFlo);
	printf("Nodes matched successfully.\n");
	
//Write out the block patches to the topology file
	if (forFlo) { 
		checkError(writeFlo(blocks,nPieces,bcFile.c_str(),floFile));
		printf(".flo file written successfully\n");
	} else /*forTop*/ {
                checkError(writeTop(blocks,bcFile.c_str(),floFile));
                printf(".top file written successfully\n");
        }

//Write out the blocks themselves
	checkError(writeBlocks(blocks,outMesh));
	printf("Program finished successfully\n");
	return 0;
}







