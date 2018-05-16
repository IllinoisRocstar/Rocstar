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
A patch is a portion of a block's face.

Orion Sky Lawlor, olawlor@acm.org, 6/13/2001
*/
#ifndef __CSAR_PATCH_H
#define __CSAR_PATCH_H

#include "gridutil.h"

#include <string>
using std::string;

class bcList;
class mblockData;
class bcListTop;

//A rectangular portion of a face 
class block;
class face;

class patch {
	face *f;
public:
	//Location in the source block:
	block *src;
	blockSpan srcSpan;
	
	void setFace(face *f_) {f=f_;}
	face *getFace(void) const { return f; }

	virtual void writeFlo(FILE *out,const bcList &bc)=0;
	virtual void writeMblock(FILE *out,mblockData &d)=0;
	virtual void writeTop(FILE *out,const bcListTop &d)=0;
	virtual bool isExternal(void) const;

	patch() {}
	virtual ~patch() {}
protected:
	patch(face *f_,block *src_,const blockSpan &srcSpan_)
	  :f(f_),src(src_),srcSpan(srcSpan_) { }
};

patch *findPatch(const blockSpan &at,const vector<patch *> &from);

//A patch that faces the external world--
// between a block and "outside"
class externalBCpatch : public patch {
	friend class block;
public:
	//Gridgen boundary condition number
	int bcNo;

	virtual void writeFlo(FILE *out,const bcList &bc);
	virtual void writeMblock(FILE *out,mblockData &d);
	virtual void writeTop(FILE *out,const bcListTop &d);
	virtual bool isExternal(void) const;

	externalBCpatch() {}
	externalBCpatch(face *f_,block *src_,
	     const blockSpan &srcSpan_,
	     int bcNo_)
	  : patch(f_,src_,srcSpan_), bcNo(bcNo_) { }
};

//An internal boundary, between blocks
class internalBCpatch : public patch {
public:
	//Location in the destination block:
	block *dest;
	blockSpan destSpan;
	internalBCpatch *partner; //Destination patch
	void setPartner(internalBCpatch *p) {partner=p;}
	int type; //Rocflo boundary condition code
	void setType(int to) {type=to;}
	orient_t orient;

	virtual void writeFlo(FILE *out,const bcList &bc);
	virtual void writeMblock(FILE *out,mblockData &d);
	virtual void writeTop(FILE *out,const bcListTop &d);

	internalBCpatch(face *f_,block *src_,block *dest_,
	     const blockSpan &srcSpan_,const blockSpan &destSpan_,
	     const orient_t &orient_)
	  : patch(f_,src_,srcSpan_), 
	    dest(dest_),destSpan(destSpan_),
	    partner(0),type(1),orient(orient_) { }
};

#endif

