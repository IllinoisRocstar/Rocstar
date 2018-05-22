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
Match up nodelists into actual patches, and
construct all type 2 boundary conditions

Orion Sky Lawlor, olawlor@acm.org, 6/8/2001
*/
#include <stdio.h>
#include "face.h"


bool patch::isExternal(void) const {return false;}
bool externalBCpatch::isExternal(void) const {return true;}

blockLoc dirToCorner(block *b,const blockLoc &dir)
{
	blockLoc ret(0,0,0);
	for (int i=0;i<3;i++)
		if (dir[i])
			ret[i] = b->getDim()[i]-1;
	return ret;
}

//Return true if any block lists us with an external BC
bool adjList::isExternal(void) const
{
	adjRec *cur=next;
	while (cur!=NULL) {
		patch *p=cur->getFace()->patchForCoord(cur->getLoc());
		if (p->isExternal())
			return true;
		cur=cur->getNext();
	}
	return false;
}

inline void swap(int &a,int &b) {int tmp=a;a=b;b=tmp;}

patch *findPatch(const blockSpan &at,const vector<patch *> &from)
{
	int p,len=from.size();
	for (p=0;p<len;p++)
		if (from[p]->srcSpan==at)
			return from[p];
	return NULL;
}

//Find the corresponding patch on the destination block
void findPartner(internalBCpatch *p)
{
	//First orient the start and end locations:
	blockSpan d=p->destSpan;
	d.orient();
	patch *partner=NULL;
	for (int f=0;partner==NULL && f<block::nFaces;f++)
		partner=findPatch(d,p->dest->getFace(f).getPatches());
	if (partner==NULL) {
		fprintf(stderr,"Couldn't find match for internal patch on block %d!\n",
			p->dest->getBlockNumber()+1);
		abort();
	}
	p->setPartner((internalBCpatch *)partner);
}

//Return true if a's face is "smaller" than b's
bool isSmaller(face *a,face *b) {
	int a_area=a->getDim().getSize();
	int b_area=b->getDim().getSize();
	return a_area<=b_area;
}

void buildFaces(vector<block *> &blocks,bool buildTypeTwo)
{
	int f; unsigned int bNo;
	/* Build each block's faces, which matches up the nodes*/
	nodeMatcher map;
	for (bNo=0;bNo<blocks.size();bNo++)
		blocks[bNo]->buildFaces(map);

	/* Build each face's patches, which aggregates matched nodes
	   into rectangular patches.*/
	for (bNo=0;bNo<blocks.size();bNo++)
		for (f=0;f<block::nFaces;f++)
	        	blocks[bNo]->getFace(f).buildPatches();

	/*Match up facing internal patches*/
	for (bNo=0;bNo<blocks.size();bNo++)
		for (f=0;f<block::nFaces;f++) {
			const vector<patch *> &patches=
			    blocks[bNo]->getFace(f).getPatches();
			for (unsigned int pNo=0;pNo<patches.size();pNo++) {
				patch *p=patches[pNo];
				if (p->isExternal()) continue;
				findPartner((internalBCpatch *)p);
			}
		}

	if (!buildTypeTwo) return;
	/*
	Determine which patches should be "type two"--
	patches whose neighbors get pulled around by 
	mesh motion, but don't move themselves.
	Such patches touch a corner node that faces the
	exterior; but none of them are exterior themselves.
	*/
	const int nCorners=8;
	blockLoc cornDirs[nCorners]={
		blockLoc(0,0,0),
		blockLoc(0,0,1),
		blockLoc(0,1,0),
		blockLoc(0,1,1),
		blockLoc(1,0,0),
		blockLoc(1,0,1),
		blockLoc(1,1,0),
		blockLoc(1,1,1)
	};
	/*Loop over blocks*/
	for (bNo=0;bNo<blocks.size();bNo++) {
		block *b=blocks[bNo];
		/*Loop over corners*/
		for (int cNo=0;cNo<nCorners;cNo++) {
			blockLoc c=dirToCorner(b,cornDirs[cNo]);
			/*Determine if this corner faces external world*/
			node *n=map.loc2node(b->getLoc(c));
			if (!n->isExternal())
				continue; /*Nothing to do*/
			/*Corner faces outside-- check if any of our patches do:*/
			bool hasExternalFace=false;
			for (f=0;f<block::nFaces;f++) {
				patch *p=b->getFace(f).patchForCoord(c);
				if (p!=NULL && p->isExternal())
					hasExternalFace=true;
			}
			if (hasExternalFace)
				continue; /*Nothing to do*/
			/*Corner is outside, but none of our faces are:
			mark all patches as type 2.
			*/
			printf("Corner %d of block %d is external> patches on ",
				cNo,b->getBlockNumber()+1);
			for (f=0;f<block::nFaces;f++) {
				patch *p=b->getFace(f).patchForCoord(c);
				if (p!=NULL) {
					printf("%s ",block::face2name[f]);
					((internalBCpatch *)p)->setType(2);
				}
			}
			printf("are all type 2\n");
		}
	}
	/*Loop over patches, eliminating facing pairs of type 2 boundaries*/
	for (bNo=0;bNo<blocks.size();bNo++) {
		block *b=blocks[bNo];
		for (f=0;f<block::nFaces;f++) {
			const vector<patch *> &patches=b->getFace(f).getPatches();
			for (unsigned int pNo=0;pNo<patches.size();pNo++) {
				patch *p=patches[pNo];
				if (!p->isExternal()) {
					internalBCpatch *i=(internalBCpatch *)p;
					if (i->type==2 && i->partner->type==2)
					/*A facing pair of type 2's!*/ 
						i->type=i->partner->type=1;
				}
			}
		}
	}

	/*Loop over patches one final time, adding type 2's to the 
	  smaller face of each pair.*/
	for (bNo=0;bNo<blocks.size();bNo++) {
		block *b=blocks[bNo];
		for (f=0;f<block::nFaces;f++) {
			const vector<patch *> &patches=b->getFace(f).getPatches();
			for (unsigned int pNo=0;pNo<patches.size();pNo++) {
				patch *p=patches[pNo];
				if (!p->isExternal()) {
					internalBCpatch *i=(internalBCpatch *)p;
					if (i->type==1 && i->partner->type==1
						&& isSmaller(i->getFace(),i->partner->getFace()))
						i->type=2;
				}
			}
		}
	}
          
}


