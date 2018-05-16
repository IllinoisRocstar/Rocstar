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
Check the quality of the input mesh, using the
same criteria used by Rocflo.

Inside Rocflo, mesh volumes are computed using the routines:
  RocfluidMP/libfloflu/controlVolume.F90 : VolumeHexa 
     Called by RFLO_calcControlVolumes.F90 and stored in "vol".
  RocfluidMP/libfloflu/faceVector.F90 : FaceVectorQuad  
     Called by RFLO_calcFaceVectors.F90 and stored in "s[ijk]".

These volumes are used by RocfluidMP/rocflo/RFLO_checkMetrics.F90, 
which checks:
  Face area-normal-sum-max against 0.1*(shortest edge squared).
  Cell volume checked against 0.1*(shortest edge cubed).

Orion Sky Lawlor, olawlor@acm.org, 2004/7/6
*/
#include <stdio.h>
#include "adj.h"
#include "ckstatistics.h"

void checkQuality(const block &b) 
{
	int nBadVol=0, nBadFace=0;
	double minVolRatio=1.0;
	blockLoc aBadVol, aBadFace;
	CkSample edgeStats, volumeStats, volRatioStats, faceRatioStats;
	/* loop over all cells in the block */
	blockLoc cellCoor;
	blockSpan cells(blockLoc(0,0,0),b.getDim()-blockLoc(1,1,1));

	enum {nCorners=8,nEdges=12,nFaces=6};
	vector3d c[nCorners];
	int i,j,k,f,e;

	/* Check the quality of each cell */
	BLOCKSPAN_FOR(cellCoor,cells) {

		// Copy out the locations of the cell corners
		for (k=0;k<2;k++) for (j=0;j<2;j++) for (i=0;i<2;i++)
			c[i+2*j+4*k]=b.getLoc(cellCoor+blockLoc(i,j,k));

	/* Find the minimum edge length for this cell.
		FIXME: Rocflo's criterion depends on the minimum edge 
	   length for *any* cell in the (partitioned) block.  This means
	   our per-cell criterion is actually more stringent than Rocflo's.
	 */
		double minEdge=1.0e30;
		
		// Compute the length of the shortest edge
		const static int edges[nEdges][2]={
			{0,1}, {1,3}, {3,2}, {2,0},
			{4,5}, {5,7}, {7,6}, {6,4},
			{0,4}, {1,5}, {3,7}, {2,6}
		};
		for (e=0;e<nEdges;e++) {
			double l=c[edges[e][0]].dist(c[edges[e][1]]);
			if (l<minEdge) minEdge=l;
			edgeStats.add(l);
		}

		// Compute the outward-pointing normals of the faces
		const static int faces[nFaces][4]={
			{0,4,6,2}, {1,3,7,5}, /* imin, imax */
			{0,1,5,4}, {2,6,7,3}, /* jmin, jmax */
			{0,2,3,1}, {4,5,7,6}  /* kmin, kmax */
		};
		vector3d normals[nFaces]; // area-scaled normal vector
		vector3d center[nFaces]; // position of centroid
		for (f=0;f<nFaces;f++) 
		{ // Area and normal come from face diagonals:
			vector3d diagA=c[faces[f][2]]-c[faces[f][0]];
			vector3d diagB=c[faces[f][3]]-c[faces[f][1]];
			normals[f]=0.5*diagA.cross(diagB);
			center[f]=0.25*(c[faces[f][0]]+c[faces[f][1]]+
			                c[faces[f][2]]+c[faces[f][3]]);
		}
		
		// Check face vector sum (closedness, smaller is better)
		vector3d faceVecSum=
		         normals[0]+normals[1]+
			 normals[2]+normals[3]+
			 normals[4]+normals[5];
		double faceRatio=faceVecSum.max()/(minEdge*minEdge);
		if (faceRatio>0.1) { /* non-closed face */
			nBadFace++;
			aBadFace=cellCoor;
		}
		
		// Compute volume of cell relative to smallest edge (fullness, bigger is better)
		double volume=0.0;
		for (f=0;f<nFaces;f++) 
			volume+=(1.0/3.0)*normals[f].dot(center[f]);
		
		double volRatio=volume/(minEdge*minEdge*minEdge);
		if (volRatio<0.1) { /* bad skinny cell */
			nBadVol++;
			if (volRatio<minVolRatio) {
				minVolRatio=volRatio;
				aBadVol=cellCoor;
			}
		}
		
		// Update statistics
		faceRatioStats.add(faceRatio);
		volumeStats.add(volume);
		volRatioStats.add(volRatio);
	}
	
	printf("Original block %d mesh quality:\n",1+b.getOriginalNumber());
	printf("    cell volumes: "); volumeStats.printMinAveMax(stdout);
	printf("    edge lengths: "); edgeStats.printMinAveMax(stdout);
	printf("    volume ratio (fullness): "); volRatioStats.printMinAveMax(stdout);
	printf("    face ratio (skewness): "); faceRatioStats.printMinAveMax(stdout);
	
	if (nBadVol>0) {
	        printf("WARNING: detected %d bad-volume cells in original block %d: \n"
			"     the worst of which is cell (%d,%d,%d), with ratio %.3f\n",
		   nBadVol, b.getOriginalNumber(),
		   1+aBadVol[0],1+aBadVol[1],1+aBadVol[2], minVolRatio);
	}
}


