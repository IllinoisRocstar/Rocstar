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
Routines to write a .mblock mesh description file,
given a set of intersecting blocks.

Orion Sky Lawlor, olawlor@acm.org, 7/17/2001
*/
#include <stdio.h>
#include "makeflo.h"
#include "face.h"

class mblockData {
	int nBlock;
	vector<block *> &blocks;
	//This class stores the patches for each face of a block
	class blockPatches {
	public:
		vector<patch *> patches[block::nFaces];
		void init(block *b) {
			for (int f=0;f<block::nFaces;f++)
				patches[f]=b->getFace(f).getPatches();
		}
		void write(FILE *out,mblockData &d,block *b);
		void swap(int &a,int &b) {int tmp=a;a=b;b=tmp;}
		int getPatchNumber(patch *forPatch,int *retFace=NULL)
		{
			//Find that patch in our list
			int patchNo=0;
			for (int f=0;f<block::nFaces;f++)
			for (unsigned int p=0;p<patches[f].size();p++)
				if (patches[f][p]==forPatch) {
					if (retFace!=NULL) *retFace=f;
					return patchNo;
				}
				else 
					patchNo++;
			fprintf(stderr,"Can't match patch in %s!\n",__FILE__);
			abort();
			return 0;
		}
		int totalPatches(void) const {
			int ret=0;
			for (int f=0;f<block::nFaces;f++)
                          ret+=patches[f].size();
			return ret;
		}
	};
	blockPatches *block2patches;
public:
	mblockData(vector<block *> &blocks_)
		:blocks(blocks_)
	{
		int b;
		nBlock=blocks.size();
		block2patches=new blockPatches[nBlock];
		for (b=0;b<nBlock;b++)
			block2patches[b].init(blocks[b]);
	}
	~mblockData() {
		delete[] block2patches;
	}
	
	void write(block *b,FILE *out) {
		block2patches[b->getBlockNumber()].write(out,*this,b);
	}

	//Get this patch's number in this block's communication list
	int getPatchNumber(block *dest,patch *partner)
	{
		return block2patches[dest->getBlockNumber()].
			getPatchNumber(partner);
	}
	vector<patch *> &getPatchListWith(block *dest,patch *partner,int *hisFace=NULL)
	{
		blockPatches &b=block2patches[dest->getBlockNumber()];
		int face;
		b.getPatchNumber(partner,&face);
		if (hisFace!=NULL) *hisFace=face;
		return b.patches[face];
	}
};

const char * writeMblock(vector<block *> &blocks,
                      const char *outMblock)
{
	mblockData d(blocks);
	for (unsigned int b=0;b<blocks.size();b++) {
		char fName[1024];

		//Write the boundary descriptions
		sprintf(fName,"%s%05d.bblk",outMblock,b);
		FILE *fb=fopen(fName,"w");
		if (fb==NULL) return "Couldn't create .bblk file";
		d.write(blocks[b],fb);
		fclose(fb);
		
		//Write the mesh locations themselves
		if (!parameters.topologyOnly) {
			sprintf(fName,"%s%05d.mblk",outMblock,b);
			FILE *fm=fopen(fName,"wb");
			if (fm==NULL) return "Couldn't create .bblk file";
			fwrite(&blocks[b]->getLoc(blockLoc(0,0,0)),sizeof(vector3d),
				blocks[b]->getDim().getSize(),fm);
			fclose(fm);
		}
	}
	return NULL;
}

internalBCpatch *createSendPatch(internalBCpatch *recv,int faceNo)
{
	internalBCpatch *send=new internalBCpatch(NULL,recv->src,recv->dest,
		recv->srcSpan,recv->destSpan,recv->orient);
	return send;
}

void mblockData::blockPatches::write(FILE *out,mblockData &d,block *b)
{	
	//Write the block header:
	fprintf(out,"# Charm++ Mblock framework block boundary condition file\n");
	fprintf(out,"1.0 # Version number\n");
	fprintf(out,"%d   %d %d %d  # Block number and size\n",b->getBlockNumber(),
		b->getDim()[0],b->getDim()[1],b->getDim()[2]);
	fprintf(out,"%d %d # Number of faces, patches\n",block::nFaces,totalPatches());
	for (int f=0;f<block::nFaces;f++)
	{
		fprintf(out,"\n%d  # Number of patches on this face\n",
			(int)patches[f].size());
		for (unsigned int p=0;p<patches[f].size();p++)
			patches[f][p]->writeMblock(out,d);
	}
}	

//Print a range of grid indices
static void printSpan(FILE *out,const blockSpan &sp)
{
	const blockLoc &s(sp.start), &e(sp.end);
	fprintf(out,"   %d %d    %d %d    %d %d \n",
	       s[0],e[0],s[1],e[1],s[2],e[2]);
}


void externalBCpatch::writeMblock(FILE *out,mblockData &d) 
{
	fprintf(out,"%d ",bcNo);
	printSpan(out,srcSpan);
}

void internalBCpatch::writeMblock(FILE *out,mblockData &d) 
{ 
	fprintf(out,"-1 ");
	printSpan(out,srcSpan);
	fprintf(out,"           %d  %d   ",dest->getBlockNumber(),
		d.getPatchNumber(dest,partner));
	//Figure out the orientation of each axis
	for (int srcAxis=0;srcAxis<3;srcAxis++) {
		int destAxis=orient[srcAxis];
		int isFlipped=0; //Is this axis orientation-reversed?
		if (destSpan.end[destAxis]<destSpan.start[destAxis])
			isFlipped=1;
		fprintf(out,"%c%d ",isFlipped?'-':'+',destAxis+1);
	}
	fprintf(out,"\n");
}


