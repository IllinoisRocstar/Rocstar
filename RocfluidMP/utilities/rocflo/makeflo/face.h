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
#ifndef __CSAR_FACE_H
#define __CSAR_FACE_H

#include "adj.h"

//A facet orientation, mapping from face coordinates
// to destination coordinates.
class facetOrientation {
	blockLoc start;//Our starting location on dest
	blockLoc orientX,orientY; //Our orientation on dest	
	int startX,startY;//Offset to start location
public:
	facetOrientation(const node **nbor,const face *source,
		int startX_,int startY_)
		:startX(startX_),startY(startY_)
 	{
		node::intersect(nbor,source,&start,&orientX,&orientY);
	}

	//Map a face index to destination location
	blockLoc getDestLoc(int x,int y) const {
		return start+(x-startX)*orientX+(y-startY)*orientY;
	}
};

class facet;

//A 2D face of a 3D grid block
class face {
	block *source;//The block we front
	vector<externalBCpatch> extBCs;//External boundary conditions 
	vector<patch *> patches; //All boundary conditions for this face

	int dir_x,dir_y;//Direction indices-- the 3D (i,j,k) for our 2D (x,y)
	blockSpan span;//3D location in source block
	node **nodes; //2D node* array (nx x ny)
	int nx,ny;//Width and height of node array
	//Fill this 4-long neighbors array
	void fillNeighbors(int x,int y,const node **nbor) const {
	  nbor[0]=nodes[(x  )+nx*(y  )];
	  nbor[1]=nodes[(x+1)+nx*(y  )];
	  nbor[2]=nodes[(x  )+nx*(y+1)];
	  nbor[3]=nodes[(x+1)+nx*(y+1)];
	}

	//Check if this facet orientation matches the nodes around (x,y)
	bool nodesMatch(int x,int y,
			const face *dest,
			const facetOrientation &o) const;
	
	//Map a grid location to a node array index
	int loc2nodeX(const blockLoc &l) const
	  {return l[dir_x]-span.start[dir_x];}
	int loc2nodeY(const blockLoc &l) const
	  {return l[dir_y]-span.start[dir_y];}
	int loc2node(const blockLoc &l) const 
	  {return loc2nodeX(l)+nx*loc2nodeY(l);}
	
	//Map a node array index to a source grid location
	blockLoc getSourceLoc(int x,int y) const
	{
	  blockLoc ret=span.start;
	  ret[dir_x]+=x;
	  ret[dir_y]+=y;
	  return ret;
	}
	//Mark these facets as patched
	int markPatched(facet *b,const blockSpan &span);
public:
	//Attach the block to all this face's nodes
	face(block *source_,nodeMatcher &map,
	     vector<externalBCpatch> &extBCs_,int dir_x_,int dir_y_,
	     const blockSpan &span_);
	void buildPatches(void);
	node *nodeForCoord(const blockLoc &at);
	patch *patchForCoord(const blockLoc &at);
	~face();

	block *getBlock(void) const {return source;}
	blockDim getDim(void) const {return span.end-span.start;}
	
	//Find all the patch structures for this face.
	//Must be called after all faces have been created.
	const vector<patch *> &getPatches(void) const {return patches;}
};

#endif












