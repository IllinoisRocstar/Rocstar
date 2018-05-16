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
Utilities to help determine which nodes overlap,
where, and with whom.

Orion Sky Lawlor, olawlor@acm.org, 5/30/2001
*/
#include <stdio.h>
#include "makeflo.h"
#include "face.h"

/************* Memory allocation pools **********/
allocPool<adjRec> adjRec::pool;
allocPool<node> node::pool;

/************** AdjRec ***************/
void adjRec::print(void)
{
	printf("%d",b->getBlock()->getBlockNumber()+1);
	if (next!=NULL) {
		printf(" -> "); next->print();
	} else
		printf(".\n");
}

/****************** adjList *******************/

adjList::~adjList() {
	adjRec *cur=next;
	while (cur!=NULL) {
		adjRec *old=cur;			
		cur=cur->getNext();
		delete old;
	}
}
int adjList::getLength(void) const {
	int ret=0;
	for (adjRec *cur=next;cur!=NULL;cur=cur->getNext())
		ret++;
	return ret;
}
bool adjList::hasFace(const face *test) const {
	for (adjRec *cur=next;cur!=NULL;cur=cur->getNext())
		if (test==cur->getFace()) 
			return true;
	return false;
}
bool adjList::hasLoc(const face *test,const blockLoc &l) const {
	for (adjRec *cur=next;cur!=NULL;cur=cur->getNext())
		if (test==cur->getFace() && l==cur->getLoc()) 
			return true;
	return false;
}
void adjList::addFace(face *b,const blockLoc &l)
{
	if (!hasLoc(b,l))
		next=new adjRec(b,l,next);
}
//Find our location in this block
const blockLoc &adjList::getLoc(const face *in) const
{
	for (adjRec *cur=next;cur!=NULL;cur=cur->getNext())
		if (in==cur->getFace()) 
			return cur->getLoc();
	static blockLoc badLoc(-1,-1,-1); 
	return badLoc;
}

void adjList::print(void) const {
	if (next) next->print();
	else printf("Empty.\n");
}

/******************* Node *********************/
/*Return the total number of hops to go from the origin
to the location "dir".
*/
inline int nHops(const blockLoc &dir) {
	return abs(dir[0])+abs(dir[1])+abs(dir[2]);	
}

//Return some face present in all 4 nodes,
// but different from notHim.  Optionally return the	
// destination location and orientations
face *node::intersect(const node** n,
		       const face *notHim,
		       blockLoc *retLoc,
		       blockLoc *oX, blockLoc *oY)
{
/*Nodes in "n" array should be laid out like:
    [0]  [1]
    [2]  [3]
*/
	const int nNodes=4; 
	const int expectedHops[nNodes]={0,1,1,2};
	int l;

	//Check each entry in list 0,
	// and return the first that is present everywhere.
	for (adjRec *test=n[0]->next;test!=NULL;test=test->getNext()) {
		face *tFace=test->getFace();
		
		//Check if it's the "bad" face
		if (tFace==notHim) continue;
		
		//Check if it's in every list at some orientation
		blockLoc orientX,orientY;
		const blockLoc &tLoc=test->getLoc();
		bool isPresent=true;
		for (l=1;l<nNodes && isPresent;l++)
		{ //Check this list for the location
			bool inList=false;
			for (adjRec *c=n[l]->next;c!=NULL;c=c->getNext()) {
				if (c->getFace()!=tFace) 
					continue;//Keep looking
				blockLoc dir=c->getLoc()-tLoc;
				//Make sure it's the proper number of hops to this location:
				if (nHops(dir)==expectedHops[l]) {
					inList=true;
					if (l==1) orientX=dir; //X axis
					if (l==2) orientY=dir; //Y axis
				}
			}
			if (!inList) isPresent=false;
		}
		
		//This location is in every list-- return it
		if (isPresent) {
			if (retLoc!=NULL) *retLoc=tLoc;
			if (oX!=NULL) *oX=orientX;
			if (oY!=NULL) *oY=orientY;
			return tFace;
		}
	}
	//If we get here, the intersection is empty--
	// must be a boundary node
	return NULL;
}


/************** nodeMatcher *************/
/*This value gives the quantization for nodes--
 nodes less distant than this will be considered
 identical.  Because this value is used to round 
 the input coodinates to integers, extremely tiny
 values may cause integer saturation on 32-bit 
 machines.  For any usable coordinate x, we need
    scale*x < 1e9  (since 2^31 is approx. 2e9 and we want + and -)
*/
const double initNodeQuant=1.0e-6;//limits coordinates to < 1000
double hashableVector3d::scale=1.0/initNodeQuant;
static int originMapsTo=(1<<30); //Origin (0) maps here; keeps the sign bit free
double hashableVector3d::offset=originMapsTo+0.5*initNodeQuant;

/*Make sure our quantization is coarse enough to not overflow
if given a vector of this size.
*/
void hashableVector3d::checkVector(const vector3d &v) {
	const double safetyFactor=5; //<- Needed since we only check the corners
	double mag=safetyFactor*v.mag();
	while (mag*scale>1.0e9) {//If our hashableVectors would overflow:
		scale*=0.1; //Use a smaller scale (coarser quantization)
		double newQuant=1.0/scale;
		offset=originMapsTo+0.5*newQuant;
		printf("Encountered large node coordinates: quantizing to %.3g\n",newQuant);
	}
}

hashableVector3d::hashableVector3d(const vector3d &v)
{ 
	x=(int)(scale*v.x+offset); 
	y=(int)(scale*v.y+offset); 
	z=(int)(scale*v.z+offset);
}


//Map this location to a node.
// Creates a new node there if none exists
node *nodeMatcher::loc2node(const vector3d &loc) {
	hashableVector3d hLoc(loc);
	node *n=map.get(hLoc);
	if (n==NULL)
	{ //Node is not yet in table-- add it
		n=new node(loc);
		map.put(hLoc)=n;
	}
	return n;
}

/************** Block ***************/
static int blockBytes=0; //To keep track of memory usage

block::block(const blockDim &dim_,int orig,int blockNo_,vector3d *nodeLocs_)
	:dim(dim_), originalNo(orig),blockNo(blockNo_), nodeLocs(nodeLocs_) 
{
	hashableVector3d::checkVector(nodeLocs[0]);
	hashableVector3d::checkVector(nodeLocs[dim.getSize()-1]);
	blockBytes+=dim.getSize()*sizeof(nodeLocs[0]);
	for (int i=0;i<nFaces;i++)
		faces[i]=0;
}

block::~block()
{ 
	for (int i=0;i<nFaces;i++)
		delete faces[i];
	delete[] nodeLocs; 
}

//Add an external boundary condition at the given locations
void block::addBC(const blockSpan &span,int bcNo)
{
	//Find the face the BC applies to
	int faceNo=span.getFace();
	if (faceNo==-1) {fprintf(stderr,"Boundary condition is volumetric!\n");abort();}
	
	//Apply the BC to that face
	BCs[faceNo].push_back(externalBCpatch(NULL,this,span,bcNo));
}

//Create the 6 faces of our block
void block::buildFaces(nodeMatcher &map) {
	int nx=dim[0];
	int ny=dim[1];
	int nz=dim[2];
	faces[0]=new face(this,map,BCs[0],
		1,2,blockSpan(blockLoc(0,0,0),blockLoc(1,ny,nz)));
	faces[1]=new face(this,map,BCs[1],
		0,2,blockSpan(blockLoc(0,0,0),blockLoc(nx,1,nz)));
	faces[2]=new face(this,map,BCs[2],
		0,1,blockSpan(blockLoc(0,0,0),blockLoc(nx,ny,1)));
	faces[3]=new face(this,map,BCs[3],
		1,2,blockSpan(blockLoc(nx-1,0,0),blockLoc(nx,ny,nz)));	
	faces[4]=new face(this,map,BCs[4],
		0,2,blockSpan(blockLoc(0,ny-1,0),blockLoc(nx,ny,nz)));
	faces[5]=new face(this,map,BCs[5],
		0,1,blockSpan(blockLoc(0,0,nz-1),blockLoc(nx,ny,nz)));
}

//This table lists human-readable names for
// the faces of each block, in storage order.
const char *block::face2name[block::nFaces]={
	"iMin","jMin","kMin",
	"iMax","jMax","kMax"
};

/************** BlockReader **********/
const char *blockReader::consume(
	const blockDim &dim,//Dimentions of incoming block
	vector3d *locs) //X,Y,Z coordinates (ni x nj x nk)
{
	block *b=new block(dim,curBlock,curBlock,locs);
	curBlock++;
	
	for (int axis=0;axis<3;axis++)
		if ((dim[axis]-1)&parameters.levelBad) {
			fprintf(stderr,"ERROR! Block %d's side is %d long.\n",
				curBlock,dim[axis]-1);
			parameters.multigridError();
		}
	
	blocks.push_back(b);
	return NULL;
}
void blockReader::freeBlock(vector3d *locs) 
	{ /*The block will free its locs, not us*/ }


/************** Unit Test Harness **************/

#ifdef STANDALONE
/*Unit test driver program:
	Prints out the patches for each block in a 
human-readable fashion.
*/

template <class T>
void poolStats(const allocPool<T> &p) {
	printf("%.1f MB: (%d blocks serving %d allocd/%d freed)\n",
	       sizeof(T)*p.allocCount*1.0e-6,
	       p.bufferCount,p.allocCount,p.freeCount);
}

void memStats(void) {
	printf("blocks> %.1f MB\n",blockBytes*1.0e-6);
	printf("adjRecs> ");poolStats(adjRec::pool);
	printf("nodes> ");poolStats(node::pool);
}

void printSpan(const blockLoc &a,const blockLoc &b)
{
	printf("(%d %d) (%d %d) (%d %d)\n",
	       1+a[0],1+b[0],1+a[1],1+b[1],1+a[2],1+b[2]);
}

void printPatch(const patch &p)
{
	printf("     %d's ",
	       p.src->getBlockNumber()+1);
	printSpan(p.srcStart,p.srcEnd);
	if (p.dest==NULL)
		printf("  -> exterior\n");
	else {
		printf("  -> %d's ",p.dest->getBlockNumber()+1);
		printSpan(p.destStart,p.destEnd);
		printf("    Orient=(%d %d %d)\n",
			1+p.orient[0],1+p.orient[1],1+p.orient[2]);
	}
}

int main(int argc,char *argv[]) 
{
	if (argc<2) {printf("Usage: test <.grd file>\n"); return 1;}
	vector<block *> blocks;
	blockReader r(blocks);
	const char *err;
	printf("Starting read...\n");
	if (NULL!=(err=r.read(argv[1])))
		printf("ERROR! %s\n",err);
	else
		printf("File read successfully.\n");
	
	nodeMatcher map;	

	//Print out the patches of each block
	for (int bn=0;bn<blocks.size();bn++) {
		memStats();
		block *b=blocks[bn];
		printf("Block %d:\n",b->getBlockNumber()+1);
		b->buildFaces(map);
		for (int f=0;f<block::nFaces;f++) {
			printf("%d> Face %d:\n",b->getBlockNumber()+1,f);
			vector<patch> patches;
			b->faces[f]->getPatches(patches);
			for (int p=0;p<patches.size();p++) 
				printPatch(patches[p]);
		}
	}

	return 0;
}

#endif












