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
Read exterior boundary conditions from a
Gridgen .inp file.

This file is written by Pointwise(tm) Gridgen (ver 15) as:
From the top level menu:
  Input/Output (e) ->
  Export Analysis BCs (5) ->
  Name (k) ->
  (type in name, ending with ".inp")

Orion Sky Lawlor, olawlor@acm.org, 6/7/2001
*/
#include <stdio.h>
#include "makeflo.h"
#include "adj.h"

class bcReader {
	FILE *in;
public:
	bcReader(FILE *in_) 
		:in(in_)
	{
		
	}
	~bcReader() {fclose(in);}
	
	//Read an integer
	bool read(int *dest) {
		int count=fscanf(in,"%d",dest);
		return (count==1);
	}
	//Read a (iMin-iMax) (jMin-jMax) (kMin-kMax) 1-based location
	bool read(blockSpan &ret) {
		blockLoc start,end;
		for (int i=0;i<3;i++) {
			if (!(read(&start[i]) && read(&end[i])))
				return false;
			//Fix up some oddities in the .inp files
			if (end[i]<0) end[i]=-end[i];
			if (start[i]<0) start[i]=-start[i];
			
			//Make zero-based
			end[i]--; start[i]--;

			if (end[i]<start[i]) 
			{ //Start and end are backwards for this axis:
				int tmp=end[i];
				end[i]=start[i];
				start[i]=tmp;
			}
			end[i]++; //Make end into c-style
		}
		ret=blockSpan(start,end);
		return true;
	}
	//Read a dimension
	bool read(blockDim &d) {
		for (int i=0;i<3;i++)
			if (!read(&d[i]))
				return false;
		return true;
	}
	//Skip a line of input
	bool skipLine(void) {
		char str[1024];
		return NULL!=fgets(str,1024,in);
	}
};

void checkSpan(const blockSpan &s,const char *what) {
	for (int axis=0;axis<3;axis++) 
		if ((s.start[axis]&parameters.levelBad)
		  ||((s.end[axis]-1)&parameters.levelBad)) {
			fprintf(stderr,"%s",what);
			parameters.multigridError();
		}
}

#define BC_EXTENTION ".inp"
#define ERR "Error reading " BC_EXTENTION " file: "

const char * readBoundaries(vector<block *> &blocks,
		      const char *inMesh)
{
	FILE *in=fopen(inMesh,"r");
	if (in==NULL) {
		fprintf(stderr,"Cannot open boundary condition file '%s'.\n",inMesh);
		fprintf(stderr,"Continuing without external boundary conditions\n");
		return NULL;//Ignore missing BC file
	}
	bcReader r(in);

	//Check the file header:
	int solver;
	if (!r.read(&solver)) return ERR "Bad file header (line 1)";
	if (solver!=1) return ERR "Must use Gridgen generic solver";
	int nBlocks;
	if (!r.read(&nBlocks)) return ERR "Bad file header (line 2)";
	if (nBlocks!=(int)blocks.size()) return ERR "mesh's block count does not match";
	
	//Read each block
	for (int bn=0;bn<nBlocks;bn++) {
		block *b=blocks[bn];
		//Check the block size
		blockDim dim;
		if (!r.read(dim)) return ERR "Cannot read block size";
		if (dim!=b->getDim()) return ERR "mesh's block size does not match";
		//Skip over the block name line
		r.skipLine();
		r.skipLine();

		//Get the patch count and loop over patches
		int nPatch;
		if (!r.read(&nPatch)) return ERR "Cannot read patch count";
		for (int p=0;p<nPatch;p++) {

			//Get the patch dimensions
			blockSpan span;
			if (!r.read(span)) return ERR "Cannot read patch dimensions";
			//Get the patch type and handle patch
			int type;
			if (!r.read(&type)) return ERR "Cannot read patch type";
			if (type>=0) 
			{//External boundary condition--add to block
				b->addBC(span,type);
				char errBuf[200];
				sprintf(errBuf,"Bad boundary condition for block %d\n",1+bn);
				checkSpan(span,errBuf);
			}
			else
			{//Internal boundary condition-- read and ignore
				blockSpan dest;
				if (!r.read(dest)) return ERR "Cannot read internal patch dimensions";
				int destBlock;
				if (!r.read(&destBlock)) return ERR "Cannot read destination block";
				destBlock--; //Make zero-based
			}
		}
	}
	return NULL;//Everything worked
}



