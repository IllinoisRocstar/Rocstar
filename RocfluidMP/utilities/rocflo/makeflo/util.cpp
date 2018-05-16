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
Determines the format of, and calls the appropriate
routine to read a CSAR fluid block mesh.

Orion Sky Lawlor, olawlor@acm.org, 6/8/2001
*/
#include <cstdio>
#include <strings.h>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include "util.h"

//Return true if a ends in suffix
bool endsWith(const char *a,const char *suffix)
{
	int al=strlen(a);
	int sl=strlen(suffix);
	if (al<sl) return false;
	const char *a_suff=&a[al-sl];
	return 0==strcmp(a_suff,suffix);
}


//Return true if the indicated file exists and is readable
bool fileExists(const char *fName) {
	FILE *f=fopen(fName,"rb");
	if (f==NULL) 
		return false;
	else {
		fclose(f);
		return true;
	}
}

//Find the numeric part of the file name, and increment it 
// by one.
bool incrementAscii(char *cur)
{
	//Starting from end, work backwards until we reach
	// the first decimal digit.
	int i=strlen(cur)-1;
	while (i>=0 && !isdigit(cur[i])) i--;
	if (i<0) return false; //No digits found
	
	//Carry loop
	while (true) {
		//Increment the current digit
		cur[i]++;
		if (cur[i]!='9'+1) 
			break; //No carry needed-- finished
		//Otherwise, wrap around and carry to next digit
		cur[i]='0';
		i--;
		if (i<0 || !isdigit(cur[i])) return false; //Too many carries
	}
	return true;
}


//Return the filename, padded for fortran use
// Returns a reference to a static internal buffer  
//  (NOT THREADSAFE)
const char *fortranifyString(const char *src)
{
	const int bufLen=500; 
	static char buf[bufLen];
	strcpy(buf,src);
	//Pad with spaces 
	for (int i=strlen(buf);i<bufLen;i++) buf[i]=' ';
	buf[bufLen-1]=0; //Zero-terminate for C use
	return buf;
}

//Replace this string's extention with the given new string
string replaceExtention(const string &a,const string &newSuffix)
{
	string ret(a);
	int dotLoc=ret.rfind(".");
	if (dotLoc==-1) return ret.append(newSuffix);
	return ret.erase(dotLoc).append(newSuffix);
}






