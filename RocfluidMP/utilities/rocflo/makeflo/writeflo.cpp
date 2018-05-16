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
Routines to write a .flo mesh description file,
given a set of intersecting blocks.

Orion Sky Lawlor, olawlor@acm.org, 6/7/2001
*/
#include <stdio.h>
#include "face.h"

/****************** .bc file input *****************/
class bcRec {
	int bcNo;//Gridgen b.c. number
	bool foundLoc; //Have we seen the location marker?
	string pre,post;//Portion before and after location
public:	
	bcRec(int n)
		:bcNo(n),foundLoc(false),pre(""),post("") { }
	void addLine(string l) {
		if (!foundLoc) 
	       	{ //Check the line for a location marker
			int loc=l.find("$");
			if (loc==-1) 
				pre+=l; //no marker
			else { //We found the marker!
				foundLoc=true;
				pre+=l.substr(0,loc);
				post+=l.substr(loc+1);
			}
		} 
		else //Already found the marker
			post+=l;//Just append 
	}
	string getString(string insert) const {
		return pre+insert+post;
	}
	bool hasNumber(int n) const { return n==bcNo; }
};

class bcList {
	vector<bcRec *> b;
public:
	bcList(FILE *inF);
	string lookup(int bcNo,string insert) const;
};

//Read the boundary conditions from the file
bcList::bcList(FILE *inF) 
{
	int lineCount=0;
	char line[200];
	bcRec *curBC=NULL;
	int bcNo;
	while (NULL!=fgets(line,200,inF)) {
		lineCount++;
		switch(line[0]) {
		case '#': break; //Skip comment line
		case '\n': case '\r': case 0: 
			break; //Skip empty lines
		case '@': //Start of new boundary condition
			bcNo=0;
			if (1!=sscanf(&line[1],"%d",&bcNo)) {
				fprintf(stderr,"ERROR!\n"
					"Can't parse boundary condition number from '%s',\n"
					"found on line %d of the boundary condition file\n",
					line,lineCount);
				exit(1);
			}
			//Flush the old bc:
			if (curBC!=NULL) b.push_back(curBC);
			curBC=new bcRec(bcNo);
			break;
		default: //Just a regular description line
			if (curBC!=NULL) curBC->addLine(line);
			break;
		};
	}
	if (curBC!=NULL) b.push_back(curBC);
}

//Find a boundary condition and return its string
string bcList::lookup(int bcNo,string insert) const 
{
	unsigned int i;
	for (i=0;i<b.size();i++) 
		if (b[i]->hasNumber(bcNo)) 
			return b[i]->getString(insert);
	//If we got here, we couldn't find bcNo
	fprintf(stderr,"ERROR! Can't find boundary condition %d in .bc file!\n",bcNo);
	exit(1);
	return "";//<- for whining compilers
}

/****************** .flo file output formatting **************/

//Write the .flo file header
void writeHeader(FILE *out,int nBlocks,int nPEs) {
	fprintf(out,
		"%d %d       ! ndom, nproc\n"
		"\n"
		"PROC_MAP_AUTO\n"
		"\n",
		nBlocks,nPEs);
}


//Print a range of grid indices
string getSpan(const blockSpan &s)
{
	char buf[200];
	sprintf(buf,"   %d %d    %d %d    %d %d",
		1+s.start[0],s.end[0],
		1+s.start[1],s.end[1],
		1+s.start[2],s.end[2]);
	return buf;
}

//An exterior boundary-- lookup the type
void externalBCpatch::writeFlo(FILE *out,const bcList &bc)
{
	string srcS=getSpan(srcSpan);
	string bcDesc=bc.lookup(bcNo,srcS);
	fprintf(out,"%s",bcDesc.c_str());
}

//An interior boundary-- just print it
void internalBCpatch::writeFlo(FILE *out,const bcList &bc)
{
	string srcS=getSpan(srcSpan);
	string destS=getSpan(destSpan);
	int destBlockNo=dest->getBlockNumber()+1;
	const char *selfFlag="";
	if (destBlockNo==src->getBlockNumber()+1) {
		destBlockNo=0;//Self-connecting block
		selfFlag="(self-connecting)";
	}
	fprintf(out,"%d      %s    ! Internal boundary %s%s\n",
		type,srcS.c_str(),(type==2)?"(TYPE TWO)":"",selfFlag);
	fprintf(out,"%-5d  %s\n",
		destBlockNo,destS.c_str());
	fprintf(out,"%d %d %d           ! Orientation\n",
		1+orient[0],1+orient[1],1+orient[2]);
}


//Write a block's data
void writeBlock(FILE *out,const bcList &bc,const block *b) {
	int blockNo=b->getBlockNumber()+1;
	fprintf(out,
		"%d %d 1 ! BLOCK %d (Split from source block %d) ======================\n"
		"\n",
		blockNo,blockNo,blockNo,b->getOriginalNumber()+1);
	blockDim d=b->getDim();
	fprintf(out,
		"%d %d %d     ! Block size (ni, nj, nk)\n",
		d[0],d[1],d[2]);
	//Get the fluid initial conditions
	fprintf(out,"%s\n",bc.lookup(-1,"").c_str());

	//Loop over the faces
	for (int f=0;f<block::nFaces;f++) {
		const vector<patch *> &patches=b->getFace(f).getPatches();
		int nPatches=patches.size();
		
		fprintf(out,
			"%d            ! Block %d, Face (%s) has %d patch\n",
			nPatches,blockNo,block::face2name[f],nPatches);
		//Loop over the patches	
		for (int p=0;p<nPatches;p++) {
			patches[p]->writeFlo(out,bc);
		}
		fprintf(out,"\n");
	}
	fprintf(out,
		"ENDBBC\n"
		"\n");
}

const char * writeFlo(vector<block *> &blocks,
		      int nPEs,
		      const char *inBcs,
		      const char *outFlo)
{
	FILE *bcs=fopen(inBcs,"r");
	if (bcs==NULL) {
		char *ret=(char *)malloc(sizeof(char)*1000);
		sprintf(ret,"Couldn't open input .bc file '%s'!\n",inBcs);
		return ret;
	}
	bcList bc(bcs);
	fclose(bcs);

	FILE *flo=fopen(outFlo,"w");
	if (flo==NULL) return "Couldn't open output .flo file!\n";

	//Print the .flo file header
	int nBlocks=blocks.size();
	writeHeader(flo,nBlocks,nPEs);
  
	//Print out each block
	for (int bn=0;bn<nBlocks;bn++)
		writeBlock(flo,bc,blocks[bn]);
	fclose(flo);
	return NULL; //Everything worked!
}













