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
Global parameters and declarations for makeflo.

Orion Sky Lawlor, olawlor@acm.org, 2/1/2002
*/

#ifndef __CSAR_MAKEFLO_H
#define __CSAR_MAKEFLO_H
#include <cstring>
#include <cstdlib>
#include "adj.h"
#include "util.h"
#define MAKEFLO_VERSION "1.93"

void checkQuality(const block &b);

const char * readBoundaries(vector<block *> &blocks,
		      const char *inMesh);

const char * splitBlocks(vector<block *> &blocks,
		      int nPieces);

void buildFaces(vector<block *> &blocks,bool buildTypeTwo);

const char * writeFlo(vector<block *> &blocks,int nPEs,
		      const char *inBcs,
		      const char *out);
const char * writeTop(vector<block *> &blocks,
		      const char *inBcs,
		      const char *out);

const char * writeMblock(vector<block *> &blocks,
                      const char *outMblock);

const char * writeBlocks(vector<block *> &blocks,
		      const char *outMesh);

class makefloParam {
 public:
	int skipAxis; //Axis to skip face builds for, or -1 if none
	int topologyOnly; //Don't write out grid locations

	int splitAxis; //Axis to split along, or -1 for any
        int splitRCB; //If 1, split using recursive bisection (optimizes communication)
        int splithalf;  // If 1, always split every block in half
	
	int nLevels; //Number of multigrid levels to create
	unsigned int levelBad; //Cuts and boundaries must not have these bits set
	unsigned int levelGood; //Cuts and boundary can only have these bits set

	makefloParam() {
		skipAxis=-1;
		topologyOnly=0;
		splitAxis=-1;
		splitRCB=0;
		splithalf=0;
		setLevel(1);
	}
	void setLevel(int n) { //Use n multigrid levels
		nLevels=n;
		levelBad=(1<<(n-1))-1;
		levelGood=~levelBad;
	}
	//Print a multigrid-size diagnostic and exit
	void multigridError(void);
};
extern makefloParam parameters;


#endif


