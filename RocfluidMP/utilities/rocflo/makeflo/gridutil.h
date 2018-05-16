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
Simple grid manipulation routines.

Orion Sky Lawlor, olawlor@acm.org, 5/30/2001
*/
#ifndef __CSAR_GRIDUTIL_H
#define __CSAR_GRIDUTIL_H

#include <math.h>
#include <iostream>
using std::ostream;
using std::endl;
using std::cout;
using std::cerr;

//This typedef allows us to easily change the floating-point
// type used in all geometry calculations.
typedef float real;

#include "vector3d.h"


//Permutation of (i,j,k) from source to dest:
//  Source axis i maps to dest. axis orient[i]. (zero-based)
class orient_t {
  int orient[3];
 public:
  int &operator[](int i) {return orient[i];}
  const int &operator[](int i) const {return orient[i];}
  orient_t inverse(void) const {
    orient_t ret;
    for (int i=0;i<3;i++) ret.orient[orient[i]]=i;
    return ret;
  }
};


//An i,j,k location in one 3D grid block
class blockDim;
class blockLoc {
protected:
	int i,j,k;
public:
	blockLoc() { }
	blockLoc(int i_,int j_,int k_) 
		:i(i_), j(j_), k(k_) { }
	blockLoc operator+ (const blockLoc &b) const
	  { return blockLoc(i+b.i,j+b.j,k+b.k); }
	blockLoc &operator+= (const blockLoc &b)
	  { i+=b.i; j+=b.j; k+=b.k; return *this; }
	blockDim operator- (const blockLoc &b) const;
	friend blockLoc operator*(const blockLoc &a,int k) 
	  { return blockLoc(a.i*k,a.j*k,a.k*k); }
	friend blockLoc operator*(int k,const blockLoc &a) 
	  { return a*k; }

	bool operator==(const blockLoc &o) const 
		{return i==o.i && j==o.j && k==o.k; }
	bool operator!=(const blockLoc &o) const 
		{return i!=o.i || j!=o.j || k!=o.k; }

	//Dimension indexing
	int &operator[](int d) {return (&i)[d];}
	int operator[](int d) const {return (&i)[d];}
};

//The i, j, and k dimentions of one block
class blockDim : public blockLoc {
public:
	blockDim() { }
	blockDim(int i_,int j_,int k_) 
		:blockLoc(i_,j_,k_)
		{ }
	int getSize(void) const 
		{ return i*j*k; }
	//Return the (0-based) array index of the (0-based) point (xi,xj,xk)
	int c_index(int xi,int xj,int xk) const 
		{ return xi+i*(xj+j*xk);  }

	//Shorthand for above
	int operator[](const blockLoc &l) const
	  { return c_index(l[0],l[1],l[2]); }

	//Dimension indexing
	int &operator[](int d) {return (&i)[d];}
	int operator[](int d) const {return (&i)[d];}
};

inline blockDim blockLoc::operator- (const blockLoc &b) const
       { return blockDim(i-b.i,j-b.j,k-b.k); }

//Some subset of a block
class blockSpan {
	//Swap the endpoints of the range (start,end+1)
	static void swapSpan(int &start,int &end) {
		end--;
		int tmp=start;
		start=end;
		end=tmp;
		end++;
	}
public:
	blockLoc start; //First included grid location
	blockLoc end; //Last included grid location PLUS 1 ON EACH AXIS

	blockSpan() { }
	blockSpan(const blockLoc &s,const blockLoc &e) 
		:start(s), end(e) { }

	blockDim getDim(void) const { return end-start; }

	//Swap me so start and end are sensible
	void orient(void) {
		for (int axis=0;axis<3;axis++) {
			if (start[axis]>=end[axis]) {
				swapSpan(start[axis],end[axis]);
			}
		}
	}
	//Swap both of us so my start and end is sensible
	void orient(blockSpan &dest,const orient_t &src2dest) {
		for (int axis=0;axis<3;axis++) {
			if (start[axis]>=end[axis]) {
				swapSpan(start[axis],end[axis]);
				int dAxis=src2dest[axis];
				swapSpan(dest.start[dAxis],dest.end[dAxis]);
			}
		}
	}

	//Return the axis we have no thickness along, or -1 if none
	int getFlatAxis(void) const {
		for (int axis=0;axis<3;axis++)
			if (start[axis]+1==end[axis])
				return axis;
		return -1;
	}
	//Return the block::face number we apply to, or -1 if none
	int getFace(void) const {
		int axis=getFlatAxis();
		if (axis==-1) return -1;
		if (start[axis]==0) return axis;
		else return axis+3;
	}

	//Return true if we contain the given location
	bool contains(const blockLoc &l) const
	{
		for (int axis=0;axis<3;axis++)
			if (!(start[axis]<=l[axis] && l[axis]<end[axis]))
				return false;
		return true;
	}

	//Return true if this span represents a planar nonempty block:
	bool hasArea(void) const {
		int nThick=0;
		for (int axis=0;axis<3;axis++) {
			if (start[axis]>=end[axis])
				return false; //A completely empty block
			if (start[axis]<end[axis]-1)
				nThick++; //Some nonzero thickness on this axis
		}
		return 2==nThick;		
	}

	blockSpan operator+(const blockLoc &l) const 
		{return blockSpan(start+l,end+l);}
	blockSpan operator-(const blockLoc &l) const 
		{return blockSpan(start-l,end-l);}

	bool operator==(const blockSpan &o) const
		{return start==o.start && end==o.end;}
	bool operator!=(const blockSpan &o) const
		{return start!=o.start || end!=o.end;}
};

#define BLOCKSPAN_FOR(i,span) \
	blockSpan loop_iter=span; \
	for (i[2]=loop_iter.start[2];i[2]<loop_iter.end[2];i[2]++) \
	for (i[1]=loop_iter.start[1];i[1]<loop_iter.end[1];i[1]++) \
	for (i[0]=loop_iter.start[0];i[0]<loop_iter.end[0];i[0]++)


//Accepts blocks read from a file.
//  Inherit from this abstract class, and call read.
class blockConsumer {
public:
	virtual ~blockConsumer();

	//Pass each block encountered in file to conume.
	// Returns an error string; or NULL on success
	const char *read(const char *gridFile);

	//Location vector block allocation/deallocation
	virtual vector3d *allocateBlock(blockDim &dim);
	virtual void freeBlock(vector3d *blk);

	//Returns error string on error; NULL on success
	virtual const char *consume(
		const blockDim &dim,//Dimentions of incoming block
		vector3d *locs)=0; //X,Y,Z coordinates (ni x nj x nk)

};

#endif









