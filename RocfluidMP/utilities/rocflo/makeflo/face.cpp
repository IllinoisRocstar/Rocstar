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
A face is one side of a block of data, which
may be adjacent to many other faces and/or the 
outside world.
Faces are split up into homogenous, rectangular
"patches" for output.

Orion Sky Lawlor, olawlor@acm.org, 6/8/2001
*/
#include <stdio.h>
#include "makeflo.h"
#include "face.h"


/*********** Orientation-matching utilities *********/

//Return the index of the first nonzero coordinate
static int firstSet(const blockLoc &l) 
{
	for (int i=0;i<3;i++)
		if (l[i]!=0) 
			return i;
	//If we get here, all coordinates are zero-- very bad
	fprintf(stderr,"ERROR! %s:%d> Orientation vector all zeros!\n",
	       __FILE__,__LINE__);
	abort();
	return -1;
}

//Find the one number in {0,1,2} not equal to a or b
static int findOther(int a,int b)
{
	int ret=0;
	while (a==ret || b==ret) 
		ret++; //Keep looking
	return ret; //Found it!
}

//Match up the src and dest axis orientations
// for the given patch using the given directions,
// which point in the face X and Y directions in
// the source and dest blocks.
static void matchOrientations(orient_t &orient,
	const blockLoc &sX,const blockLoc &sY,
	const blockLoc &dX,const blockLoc &dY)
{
	//Compute the block orientations in the face X and Y axes
	int sOX=firstSet(sX),sOY=firstSet(sY);
	int dOX=firstSet(dX),dOY=firstSet(dY);
	
	//Match up the orientations along the face
	orient[sOX]=dOX;
	orient[sOY]=dOY;
	
	//The remaining orientation is found by the
	// process of elimination
	orient[findOther(sOX,sOY)]=findOther(dOX,dOY);
}


node *face::nodeForCoord(const blockLoc &at) 
{
	int dir_z=findOther(dir_x,dir_y);
	if (at[dir_z]!=span.start[dir_z]) return NULL; //Wrong z-plane
	else return nodes[loc2node(at)];	
}

patch *face::patchForCoord(const blockLoc &at)
{
	int i,len=patches.size();
	for (i=0;i<len;i++)
		if (patches[i]->srcSpan.contains(at))
			return patches[i];
	return NULL;
}

/************** Face *************/
//Attach the block to all this face's nodes
face::face(block *source_,nodeMatcher &map,
	vector<externalBCpatch> &extBCs_,int dir_x_,int dir_y_,
	const blockSpan &span_)
	:source(source_), extBCs(extBCs_), dir_x(dir_x_), dir_y(dir_y_), 
	span(span_)
{
	nx=span.end[dir_x]-span.start[dir_x];
	ny=span.end[dir_y]-span.start[dir_y];
	if (findOther(dir_x,dir_y)==parameters.skipAxis) 
	{ //Don't bother about this face-- it's the Z face of a 2D problem
		nodes=NULL;
		return;
	}
	nodes=new node*[nx*ny];
	blockLoc i;
	BLOCKSPAN_FOR(i,span) {
		vector3d nodeLoc=source->getLoc(i);
		node *nu=map.loc2node(nodeLoc);
		nu->addFace(this,i);

		if (nu->getLength()>(3*8)) {
			fprintf(stderr,"WARNING: multiply shared node %p: (%d,%d,%d) for block %d:",
				nu,i[0],i[1],i[2],source->getOriginalNumber()+1);
			nu->print();
		}
		nodes[loc2node(i)]=nu;
	}
}

face::~face() {
	delete[] nodes;
	for (unsigned int i=0;i<patches.size();i++)
		delete patches[i];
}

//Return true if these two spans share no interior points
bool isDisjoint(const blockSpan &a,const blockSpan &b) {
	for (int axis=0;axis<2;axis++) {
		if (a.start[axis]>=b.end[axis]-1) return true;
		if (b.start[axis]>=a.end[axis]-1) return true;
	}
	return false;//No separating axis-- must have a shared point
}


//Represents a face of a single element, bordered by 4 nodes.
// used only inside getPatches
class facet {
	face *dest;//The face we're next to (NULL if an exterior facet)
	bool patched; //Have we been covered by a patch yet?
public:
	//Empty constructor for array allocation
	facet() { } 
	facet(const node **nbor,const face *source) {
		dest=node::intersect(nbor,source);
		patched=false;
	}
	face *getFace(void) const {return dest;}
	bool isPatched(void) const {return patched;}
	void markPatched(void) {patched=true;}
	bool matches(const facet &o) const {
		if (patched) return false;
		return dest==o.dest;
	}
};
//Mark these facets as patched
int face::markPatched(facet *b,const blockSpan &span)
{
	int ret=0;
	int bx=nx-1; /* ,by=ny-1; */ //Size of facet array
	const blockLoc &s=span.start, &e=span.end;
	//Mark all those facets as patched
	for (int y=s[dir_y];y<e[dir_y]-1;y++)
	for (int x=s[dir_x];x<e[dir_x]-1;x++) {
		if (b[x+bx*y].isPatched()) {
			fprintf(stderr,"FATAL ERROR patching block %d> external boundary conditions overlap!\n",
				source->getBlockNumber()+1);
			exit(1);
		}
		b[x+bx*y].markPatched();
		ret++;
	}
	return ret;	
}

//Check if this facet orientation matches the nodes around (x,y)
bool face::nodesMatch(int x,int y,
        const face *dest,const facetOrientation &o) const
{
	if (dest==NULL) return true;
	return nodes[(x  )+nx*(y  )]->hasLoc(dest,o.getDestLoc(x  ,y  ))
	    && nodes[(x+1)+nx*(y+1)]->hasLoc(dest,o.getDestLoc(x+1,y+1));
}

//Find all the patch structures for this face.
//Must be called *after* all faces have been created.
void face::buildPatches(void)
{
#if 0
	//Doublecheck that our nodes agree with our location
	for (int y=0;y<ny;y++)
	for (int x=0;x<nx;x++) {
		if (!nodes[(x  )+nx*(y  )]->hasLoc(this,getSourceLoc(x,y)))
		{
			fprintf(stderr,"ERROR! Face/node location disagreement!\n");
			abort();
		}
	}
#endif	
	if (nodes==NULL) return; //Disabled face

	//Build an array showing the connected face
	// for each *element's* face-- or "facet"
	int x,y; //Loop indices
	int bx=nx-1,by=ny-1; //Size of facet array
	facet *b=new facet[bx*by];//2D facet array (bx x by)
	for (y=0;y<by;y++)
	for (x=0;x<bx;x++) {
		const node *nbor[4];
		fillNeighbors(x,y,nbor);
		b[x+bx*y]=facet(nbor,this);
	}

	//Keep track of the number of facets patched so far
	int nFacets=bx*by;
	int nPatched=0;

	//Apply the existant internal boundary conditions
	for (unsigned int ep=0;ep<patches.size();ep++) {
		nPatched+=markPatched(b,patches[ep]->srcSpan);
	}

	//Apply the external boundary conditions
	for (unsigned int ec=0;ec<extBCs.size();ec++) {
		const externalBCpatch &ep=extBCs[ec];
		patch *p=new externalBCpatch(ep);
		p->setFace(this);
		patches.push_back(p);
		nPatched+=markPatched(b,ep.srcSpan);
	}
	
	//Extract each patch from the facet array
	int itCount=0;//Runaway loop leash
	while (nPatched<nFacets) {
		int startX=-1, startY=-1;
		//Find the first facet not yet patched
		for (y=0;y<by && startY==-1;y++)
		for (x=0;x<bx && startX==-1;x++)
		  if (!b[x+bx*y].isPatched()) { 
		    startX=x;
		    startY=y;
		  }
		//This facet touches the face "dest"
		const facet &destFacet=b[startX+bx*startY];
		face *dest=destFacet.getFace();
		const node *nbor[4];
		fillNeighbors(startX,startY,nbor);
		facetOrientation o(nbor,this,startX,startY);
		
		//Extend the patch as far right as it will go
		x=startX;
		y=startY;
		while (x<bx && b[x+bx*y].matches(destFacet) &&
			nodesMatch(x,y,dest,o))
				x++;
		int endX=x;
		if (startX==endX) {
			blockLoc badLoc=getSourceLoc(startX,startY);
			blockLoc sz=source->getDim();
			fprintf(stderr,"ERROR! Face of block %d is non-conformal at node (%d,%d,%d) of (%d,%d,%d)\n",1+source->getOriginalNumber(),badLoc[0],badLoc[1],badLoc[2],sz[0],sz[1],sz[2]);
			abort();
		}
		//Extend the patch as far down as it will go,
		// marking those facets as patched (face=NULL)
		for (y=startY;y<by;y++) {
			bool rowValid=true;
			for (x=startX;x<endX;x++)
				if (!(b[x+bx*y].matches(destFacet) &&
				      nodesMatch(x,y,dest,o)))
					rowValid=false;
			if (!rowValid) 
				break;//Another block encountered
			//If we get this far, the row is valid-- mark it
			for (x=startX;x<endX;x++)
				b[x+bx*y].markPatched();
			nPatched+=endX-startX;
		}
		int endY=y;
		

		//Build and return a patch structure
		patch *p=NULL;
		blockSpan srcSpan(getSourceLoc(startX,startY),
				  getSourceLoc(endX,endY)+blockLoc(1,1,1)
				);

#if 1		//Make sure new patch is disjoint from all old patches
		for (unsigned int otherP=0;otherP<patches.size();otherP++)
			if (!isDisjoint(srcSpan,patches[otherP]->srcSpan)) {
				fprintf(stderr,"ERROR! Patches not disjoint!\n");
				abort();
			}
#endif

		if (dest==NULL)
		{  //An (unexpected) external boundary condition
			int bcNo=0;//Default boundary condition number
			p=new externalBCpatch(this,source,srcSpan,bcNo);
			/* Warn user that we've fabricated a BC */
			printf("Makeflo WARNING: Missing boundary condition (using zero) for original block %d,\n"
			       "    split block %d, %s face, cells (%d,%d,%d) - (%d,%d,%d)\n",
			       1+source->getOriginalNumber(),
			       1+source->getBlockNumber(),block::face2name[source->getFaceNo(this)],
			       1+srcSpan.start[0],1+srcSpan.start[1],1+srcSpan.start[2], 
			       srcSpan.end[0],srcSpan.end[1],srcSpan.end[2]);
		}
		else 
		{ //An internal boundary-- find other block's location
			//Match up the orientations
			orient_t src2dest;
			matchOrientations(src2dest,
				getSourceLoc(1,0)-getSourceLoc(0,0),
				getSourceLoc(0,1)-getSourceLoc(0,0),
				o.getDestLoc(1,0)-o.getDestLoc(0,0),
				o.getDestLoc(0,1)-o.getDestLoc(0,0)
			);
			blockSpan destSpan(
				  o.getDestLoc(startX,startY),
				  o.getDestLoc(endX,endY)+blockLoc(1,1,1)
				);
			p=new internalBCpatch(
				this,source,dest->getBlock(),
				srcSpan,destSpan,
				src2dest);
			orient_t dest2src=src2dest.inverse();
			destSpan.orient(srcSpan,dest2src);
			dest->patches.push_back(new internalBCpatch(
				dest,dest->getBlock(),source,
				destSpan,srcSpan,
				dest2src
			));
		}
		//Return the new patch
		patches.push_back(p);

		if (itCount++>2000) 
		{ //There can't be 2000 actual patches-- must be an error
			fprintf(stderr,"ERROR! Runaway patching loop (%s:%d)!\n",
				__FILE__,__LINE__);
			abort();
		}
	}
	
	//Free the temporary block array
	delete[] b;
}









