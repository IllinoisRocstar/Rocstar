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
Routines to write a .top mesh topology file,
given a set of intersecting blocks.

Orion Sky Lawlor, olawlor@acm.org, 2/1/2002
*/
#include <stdio.h>
#include "makeflo.h"
#include "face.h"

/****************** .bc file input *****************/
class bcRecTop {
	int inNo;//gridgen b.c. number
public:	
	int outNo;//RocfloMP b.c. number
	int isCoupled;
	bcRecTop(int in,int out,int coupled)
		:inNo(in),outNo(out),isCoupled(coupled) { }
	bool hasNumber(int n) const { return n==inNo; }
};

class bcListTop {
	vector<bcRecTop *> b;
public:
	bcListTop(FILE *inF);
	bcRecTop *lookup(int bcNo) const;
};


//Read the boundary conditions from the file
bcListTop::bcListTop(FILE *inF) 
{
	int lineCount=0;
	char lineBuf[200];
	while (NULL!=fgets(lineBuf,200,inF)) {
		char *line=lineBuf;
		while (isspace(*line)) line++; //Skip white space
		line[strlen(line)-1]=0; //Erase newline
		lineCount++;
		switch(line[0]) {
		case '#': case '!': break; //Skip comment line
		case '\n': case '\r': case 0: 
			break; //Skip empty lines
		case ' ': case '\t':
		case '0': case '1': case '2': case '3': case '4':
		case '5': case '6': case '7': case '8': case '9':
		{ //A new boundary condition
			int in,out,coupled=0;
			if (sscanf(line,"%d%d%d",&in,&out,&coupled)<2) {
				fprintf(stderr,"ERROR!\n"
					"Can't parse boundary condition number from '%s',\n"
					"found on line %d of the boundary condition file\n",
					line,lineCount);
				exit(1);
			}
			b.push_back(new bcRecTop(in,out,coupled));
			break;
		}
		case '@': //.flo-style boundary condition
			fprintf(stderr,"ERROR!\n"
				"This appears to be a rocflo-style boundary condition file--\n"
				"RocfloMP does not use @ signs.\n");
			exit(1);
		default: //Something bizarre:
			fprintf(stderr,"ERROR!\n"
				"Can't parse line of boundary condition \n"
				"'%s', found on line %d.\n",
				line,lineCount);
			exit(1);
		};
	}
}

//Find a boundary condition and return its string
bcRecTop *bcListTop::lookup(int bcNo) const
{
	unsigned int i;
	for (i=0;i<b.size();i++) 
		if (b[i]->hasNumber(bcNo)) 
			return b[i];
	//If we got here, we couldn't find bcNo
	fprintf(stderr,"ERROR! Can't find boundary condition %d in .bc file!\n",bcNo);
	exit(1);
	return NULL;//<- for whining compilers
}

/****************** .top file output formatting **************/

//Print a face number/l1-range/l2-range
static void write5(FILE *out,int face,int l1lo,int l1hi,int l2lo,int l2hi)
{
	fprintf(out,"%d  %3d %3d  %3d %3d  ",face,l1lo,l1hi,l2lo,l2hi);
}

//Convert a run of nodes (a...b-1) to a run of blocks (o..p)
void node2block(int a,int b,int &o,int &p) {
	if (a<b) 
	{ //Nodes are oriented the right way:
		o=a;
		p=b-2; //Minus 1 for node->block, minus 1 again for excl->incl
	} else {//Nodes run the other way:
		o=a-1;
		p=b-1;
	}
}

//Print a patch extent
static void writeSpan(FILE *out,const blockSpan &s,const blockLoc &sign)
{
	//Map our face numbering (iMin, jMin, kMin, iMax, ...) to 
	// rocfloMP face numbering (iMin,iMax, jMin,jMax, ...)
	const static int faceMap[6]=
		{1,3,5,2,4,6};
	int face=faceMap[s.getFace()];
	// l1 and l2 axes parameterize the face
	int flatAxis=(s.getFace()%3);
	int l1Axis=(flatAxis+1)%3;
	int l2Axis=(flatAxis+2)%3;
	
	int l1S=s.start[l1Axis], l1E=s.end[l1Axis];
	int l2S=s.start[l2Axis], l2E=s.end[l2Axis];	
	node2block(l1S,l1E, l1S,l1E);
	node2block(l2S,l2E, l2S,l2E);
	int l1Sign=sign[l1Axis];
	int l2Sign=sign[l2Axis];
	write5(out,face,
		l1Sign*(1+l1S),l1Sign*(1+l1E),
		l2Sign*(1+l2S),l2Sign*(1+l2E));
}

static void writePatchStart(FILE *out,int patchType) {
	fprintf(out,"%4d   ",patchType);	
}
static void writePatchMiddle(FILE *out,int connBlock) {
	fprintf(out,"%5d    ",connBlock);	
}
static void writePatchEnd(FILE *out,int coupled) {
	fprintf(out,"  %d\n",coupled);	
}

void externalBCpatch::writeTop(FILE *out,const bcListTop &bc)
{
	bcRecTop *r=bc.lookup(bcNo);
	writePatchStart(out,r->outNo);
	writeSpan(out,srcSpan,blockLoc(1,1,1));
	//Destination block and span are zero for external boundaries
	writePatchMiddle(out,0);
	write5(out,0,0,0,0,0);
	writePatchEnd(out,r->isCoupled);
}

void internalBCpatch::writeTop(FILE *out,const bcListTop &bc)
{
	//Compute where the orientation-indicating minus signs should go:
	int flatAxis=(srcSpan.getFace()%3);
	int l1Axis=(flatAxis+1)%3;
	blockLoc srcSign(1,1,1), destSign(1,1,1);
	srcSign[flatAxis]=0; //Shouldn't ever print flat axis
	destSign[orient[flatAxis]]=0; 
	srcSign[l1Axis]=-1; //Flip sign on l1 axis
	destSign[orient[l1Axis]]=-1;

	const int internalBcType=30;
	writePatchStart(out,internalBcType);
	writeSpan(out,srcSpan,srcSign);
	writePatchMiddle(out,1+dest->getBlockNumber());
	writeSpan(out,destSpan,destSign);
	writePatchEnd(out,0); //Internal patches are never coupled
}


//Write a block's data
static void writeBlock(FILE *out,const bcListTop &bc,const block *b) {
	int blockNo=b->getBlockNumber()+1;
	int nGridLevels=parameters.nLevels;
	fprintf(out,
		"%d %d ! ======== BLOCK %d (Split from source block %d) ========\n",
		blockNo,nGridLevels,blockNo,b->getOriginalNumber()+1);
	blockDim d=b->getDim();
	int f,nPatches=0;
	for (f=0;f<block::nFaces;f++)
		nPatches+=b->getFace(f).getPatches().size();
	
	fprintf(out,
		"%d %d %d %d     ! number of patches; block size (ni, nj, nk)\n",
		nPatches,d[0]-1,d[1]-1,d[2]-1);
	
	//Loop over the faces, writing out each patch
	for (f=0;f<block::nFaces;f++) {
		static int order[6]={0,3,1,4,2,5};
		const vector<patch *> &patches=b->getFace(order[f]).getPatches();
		int nPatches=patches.size();
		//Loop over the patches	
		for (int p=0;p<nPatches;p++) {
			patches[p]->writeTop(out,bc);
		}
	}
}

const char * writeTop(vector<block *> &blocks,
		      const char *inBcs,
		      const char *outTop)
{
	FILE *bcs=fopen(inBcs,"r");
	if (bcs==NULL) {
		char *ret=(char *)malloc(sizeof(char)*1000);
		sprintf(ret,"Couldn't open input .bcmp file '%s'!\n",inBcs);
		return ret;
	}
	bcListTop bc(bcs);
	fclose(bcs);

	FILE *top=fopen(outTop,"w");
	if (top==NULL) return "Couldn't open output .top file!\n";

	//Print the .top file header
	int nBlocks=blocks.size();
	fprintf(top,"# RocfloMP topology file, generated by makeflo\n"
	            "# \n"
	            "%d  ! Total number of blocks in this file\n",
	            nBlocks);
  	
	//Print out each block
	for (int bn=0;bn<nBlocks;bn++)
		writeBlock(top,bc,blocks[bn]);
	fclose(top);
	return NULL; //Everything worked!
}













