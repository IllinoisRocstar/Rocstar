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
Split a set of blocks into more blocks, for
parallelism.  The subblocks should have a reasonable
aspect ratio (not too thin/flat) and should all
be roughly the same volume (for load balance).

Orion Sky Lawlor, olawlor@acm.org, 6/13/2001
*/
#include <math.h>
#include "face.h"
#include "makeflo.h"
#include <algorithm>
using std::make_heap;
using std::push_heap;
using std::pop_heap;

void makefloParam::multigridError(void) {
	fprintf(stderr,"To be compatible with multigrid level %d, this\n"
		"must be a multiple of %d.\n",
		nLevels, levelBad+1);
	exit(1);
}

//Force this location in-bounds
blockLoc block::pin(const blockLoc &l) const
{
	blockLoc ret=l;
	for (int i=0;i<3;i++) {
		if (ret[i]<0) ret[i]=0;
		if (ret[i]>dim[i]) ret[i]=dim[i];
	}
	return ret;
}

//Create a new subblock from this block
block *block::subBlock(const blockSpan &span) const
{
	//Create the new node location array
	blockDim ns=span.getDim();
	vector3d *newLocs=new vector3d[ns.getSize()];
	blockLoc i;
	BLOCKSPAN_FOR(i,span) 
		newLocs[ns[i-span.start]]=nodeLocs[dim[i]];
	block *ret=new block(ns,originalNo,-1,newLocs);
	
	//Copy over boundary conditions, remapping locations
	for (int f=0;f<nFaces;f++)
	for (unsigned int c=0;c<BCs[f].size();c++) {
		const externalBCpatch &bc=BCs[f][c];
		blockSpan bs(ret->pin(bc.srcSpan.start-span.start),
		             ret->pin(bc.srcSpan.end  -span.start));
		if (bs.hasArea())
		//Add this boundary to the new block
			ret->addBC(bs,bc.bcNo);
	}
	return ret;
}

//Split this block into n subblocks,
// inserting the subblocks into dest.
void block::split(int nPieces,vector<block *> &dest)
{
	if (nPieces<=0) {
		fprintf(stderr,"Asked to split block into %d pieces!\n",
			nPieces);
		abort();
	}
	else if (nPieces==1) {
		//Reset serial number according to new position
		blockNo=dest.size();
		dest.push_back(this);
	}
	else 
	{ //Must split into several pieces
		int splitAxis=parameters.splitAxis;
		if (splitAxis==-1) //No axis specified--
		{ //Split on longest remaining axis
			int splitLen=0;
			for (int a=0;a<3;a++)
				if (dim[a]>splitLen) {
					splitAxis=a;
					splitLen=dim[a];
				}
		}
		int splitLen=dim[splitAxis];
		
		//Determine the size of each half
		int loPieces=(int)floor(nPieces/2.0);
		int hiPieces=nPieces-loPieces;
		int loSize;
		if (parameters.splitRCB) //Always cut in half
			loSize=(int)(splitLen/2);
		else //Adapt cut location to corresponding number of pieces
			loSize=splitLen*loPieces/nPieces;
		//Round cut to the multigrid:
		loSize=parameters.levelGood&((parameters.levelBad>>1)+loSize);
		
		if (loSize==0) 
		{ /*Asked to split too far*/
			fprintf(stderr,"Couldn't split block %d any further-- you must:\n"
			" -increase the input grid resolution\n"
			" -decrease the number of requested output blocks\n"
			"%s",
			1+getOriginalNumber(),
			(parameters.levelBad)?" -decrease the multigrid resolution\n":"");
			exit(1);
		}
		
		//Split this block along that axis
		blockDim origin(0,0,0);
		blockDim loEnd=dim; loEnd[splitAxis]=loSize+1;
		blockDim hiStart=origin; hiStart[splitAxis]=loSize;
		block *lo=subBlock(blockSpan(origin,loEnd));
		block *hi=subBlock(blockSpan(hiStart,dim));
		
		//Remove ourselves-- we are completely 
		// replaced by our children
		delete this;
		
		//Recursively split the new blocks
		lo->split(loPieces,dest);
		hi->split(hiPieces,dest);
	}
}


//Return the "volume"-- computational effort
// associated with this block
int volume(const block *b) {
	return b->getDim().getSize();
}

//This class used only inside splitBlocks:
class bRec {
public:
	int b;//Block number of input block
	int vol;//Volume of input block
	int n;//Number of pieces to split into
	bRec() {}
	bRec(int b_,int vol_) {
		b=b_;
		vol=vol_;
		n=1;
	}
	//Compare vol/n for these bRecs:
	static bool less(const bRec &a,const bRec &b) {
		return (a.vol*b.n)<(b.vol*a.n);
	}
};
	

//Split this list of blocks into at least n subblocks
const char * splitBlocks(vector<block *> &blocks,
                      int nPieces)
{
	if ((int)blocks.size()>=nPieces)
		return NULL;//Already have enough blocks


//We try to split so that the resulting blocks
// have about the same volume-- distribute pieces
// to blocks by volume
	unsigned int b,nBlocks=blocks.size();
	
	int *nSplit=new int[nBlocks];	

	if (parameters.splithalf == 1) {
		// split every block in half
	  for (b=0;b<nBlocks;b++) {
		nSplit[b]=2;
	  }
	}
	else {   // regular split
	//Create max-heap of blocks, ordered by output volume
	vector<bRec> heap(nBlocks);
	for (b=0;b<nBlocks;b++) {
		heap[b]=bRec(b,volume(blocks[b]));
		nSplit[b]=1;
	}

	std::make_heap(heap.begin(),heap.end(),bRec::less);
	
	//Allocate remaining output blocks to the most needy
	int nRemaining=nPieces-nBlocks;
	while (nRemaining-->0) {
		std::pop_heap(heap.begin(),heap.end(),bRec::less);
		bRec &r=heap[nBlocks-1];
		r.n=++(nSplit[r.b]);
		std::push_heap(heap.begin(),heap.end(),bRec::less);
	}
	}   // end of parameters.splithalf
	
	//Split the blocks
	vector<block *> split;//Accumulates split blocks
	for (b=0;b<nBlocks;b++) {
		printf("Input block %d splits into %d pieces (%d-%d)\n",
			b+1,nSplit[b],(int)(1+split.size()),(int)(1+split.size()+nSplit[b]-1));
		blocks[b]->split(nSplit[b],split);
	}
	delete[] nSplit;

//Check the load balance
	double maxVol=-1e100,minVol=1e100,sumVol=0;
	for (b=0;b<split.size();b++) {
		double vol=volume(split[b]);
		sumVol+=vol;
		if (maxVol<vol) maxVol=vol;
		if (minVol>vol) minVol=vol;
	}
	double avgVol=sumVol/split.size();
	printf("Static load balance results:\n"
		"  %.2f blocks per processor, %.0f nodes per block\n"
		"  Heaviest block is %.02f times the average\n"
		"  Lightest block is %.02f times the average\n",
		((double)split.size())/nPieces,avgVol,
		maxVol/avgVol,minVol/avgVol);
	
//Return split blocks to caller
	blocks=split;
	return NULL;
}


