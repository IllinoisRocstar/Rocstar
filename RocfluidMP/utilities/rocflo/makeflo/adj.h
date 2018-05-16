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
Determine which blocks are adjacent, at which
faces.  To do this, we keep a node->face list.

Orion Sky Lawlor, olawlor@acm.org, 5/30/2001
*/
#ifndef __CSAR_ADJ_H
#define __CSAR_ADJ_H
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
using std::vector;
#include "hash_table.h"
#include "gridutil.h"
#include "patch.h"

//A trivial grouped allocator
template <class T>
class allocPool {
	enum {nAlloc=256}; //Number to allocate at once
	T *buf;
	int cur;
	void fillBuffer(void) {
	        bufferCount++;
		buf=(T*)malloc(sizeof(T)*nAlloc);
		cur=0;
	}
public:
	int allocCount,bufferCount,freeCount;
	allocPool() {
	  cur=nAlloc;
	  allocCount=bufferCount=freeCount=0;
	}
	void *alloc(void) { 
	  allocCount++;
	  if (cur>=nAlloc) fillBuffer();
	  return &buf[cur++];
	}
	void free(void *) { 
	  freeCount++; 
	  /*otherwise ignored*/ 
	}
};

class face;
class block;

//Records that the pointing node 
// has the given location in the given face
class adjRec {
	face *b;
	blockLoc l;
	adjRec *next;
public:
	static allocPool<adjRec> pool;
	void *operator new(size_t s) { return pool.alloc(); }
	void operator delete(void *ptr) { pool.free(ptr); }

	adjRec(face *b_,const blockLoc &l_,adjRec *next_=NULL) 
		:b(b_), l(l_), next(next_) { }
	face *getFace(void) {return b;}
	const face *getFace(void) const {return b;}
	const blockLoc &getLoc(void) const {return l;}
	adjRec *getNext(void) {return next;}
	void print(void);
	bool isExternal(void) const;
};

//Maintains a (possibly zero-length) list of adjacent blocks
class adjList {
 protected:
	adjRec *next;//Blocks we're adjacent to
public:
	adjList() {
		next=NULL;
	}
	~adjList();

	int getLength(void) const;

	bool hasFace(const face *test) const;
	bool hasLoc(const face *test,const blockLoc &l) const;

	//Add this block & loc to the list (if they're not already there)
	void addFace(face *b,const blockLoc &l);

	//Find our (first) location in this face
	const blockLoc &getLoc(const face *b) const;

	//Return true if any block lists us with an external BC
	bool isExternal(void) const;

	void print(void) const;
};

//A location in space, owned by (possibly several) faces
class node : public adjList {
	const vector3d &loc; //Coordinates of this node
public:
	static allocPool<node> pool;
	void *operator new(size_t s) { return pool.alloc(); }
	void operator delete(void *ptr) { pool.free(ptr); }
	
	node(const vector3d &l) :loc(l) { }
	
	//Get location in 3D space
	const vector3d &getLoc(void) const { return loc; }

	//Find our index in this face
	const blockLoc &getLoc(const face *b) const { return adjList::getLoc(b); }

	//Return some face present in all 4 nodes,
	// but different from notHim.  Optionally return the
	// destination location and orientations
	static face *intersect(const node** nodes,
			       const face *notHim=NULL,
			       blockLoc *loc=NULL,
			       blockLoc *oX=NULL, blockLoc *oY=NULL);	

};

/*A vector3D, quantized to integers
so roundoff will not affect node matching.
*/
class hashableVector3d {
  int x,y,z;
 public:
  static double scale,offset; //Converts doubles to integers
  static void checkVector(const vector3d &v);
  
  hashableVector3d(const vector3d &v);
};

/*This class accepts node locations read from
the input file and either returns a previously
created node, or creates a new node for the location.
This maps (non-unique) locations onto (unique)
node objects.
*/
class nodeMatcher {
	//This table maps node location to node records
	CkHashtableTslow<hashableVector3d,node *> map;
public:
	//Map this location to a node.  
	// Creates a new node there if none exists
	node *loc2node(const vector3d &loc);
};

//A rectangular 3D grid of nodes
class block {
	blockDim dim;//This block's size
	int originalNo; //The source block's (0-based) serial number
	int blockNo;//This block's (0-based) serial number
	vector3d *nodeLocs;//The locations of all nodes (dim.getSize() vectors)

	//Force this location in-bounds
	blockLoc pin(const blockLoc &l) const;

	//Create a new block as a subregion of this one
	block *subBlock(const blockSpan &span) const;
public:
	block(const blockDim &dim_,int originalNo_,int blockNo_,vector3d *nodeLocs_);
	~block();
	
	//Return our 0-based block number
	int getBlockNumber(void) const { return blockNo; }
	int getOriginalNumber(void) const { return originalNo; }

	const blockDim &getDim(void) const {return dim;}
	const vector3d &getLoc(const blockLoc &l) const 
	  {return nodeLocs[dim[l]];}

	//Split this block into n pieces, which go into dest
	void split(int nPieces,vector<block *> &dest);

	//Describes the 6 faces of the block
	enum {
	  iMin=0,jMin=1,kMin=2,
	  iMax=3,jMax=4,kMax=5,
	  nFaces=6
	};
	//Maps face number to a human-readable name
	static const char *face2name[nFaces];

private: 
	//External boundary conditions per face
	vector<externalBCpatch> BCs[nFaces];
	face *faces[nFaces];
public:
	//Set an external boundary condition
	void addBC(const blockSpan &span,int bcNo);

	//Build/access our face structures
	void buildFaces(nodeMatcher &map);
	face &getFace(int faceNo) {return *faces[faceNo];}
	const face &getFace(int faceNo) const {return *faces[faceNo];}
	int getFaceNo(const face *f) {
		for (int i=0;i<nFaces;i++)
			if (faces[i]==f) return i;
		return nFaces;
	}
};

/************** BlockReader **********/
class blockReader : public blockConsumer {
	int curBlock;//Index of current block
	vector<block *> &blocks;
public:
	blockReader(vector<block *> &dest) :blocks(dest) {curBlock=0;}
	
	virtual const char *consume(
		const blockDim &dim,//Dimentions of incoming block
		vector3d *locs); //X,Y,Z coordinates (ni x nj x nk)

	virtual void freeBlock(vector3d *locs);
};



#endif












