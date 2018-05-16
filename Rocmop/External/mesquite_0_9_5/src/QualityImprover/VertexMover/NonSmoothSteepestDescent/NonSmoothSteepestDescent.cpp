/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*!  \File NonSmoothSteepestDescent.cpp \brief
  
  Implements the NonSmoothSteepestDescent class member functions.
  
  \author Lori Freitag
  \date 2002-07-20 */

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "NonSmoothSteepestDescent.hpp"
#include "MeshSet.hpp"
#include "MsqTimer.hpp"

using namespace Mesquite;

NonSmoothSteepestDescent::NonSmoothSteepestDescent(ObjectiveFunction* of)
{
  objFunc=of;
  this->set_name("NonSmoothSteepestDescent");
  
  mFunction = (double *)malloc(sizeof(double)*150);
  testFunction = (double *)malloc(sizeof(double)*150);
  originalFunction = (double *)malloc(sizeof(double)*150);
  mGradient = (double **)malloc(sizeof(double *)*150);
  mGS = (double *)malloc(sizeof(double)*150);
  mG = (double **)malloc(sizeof(double *)*150);
  mPDG = (double **)malloc(sizeof(double *)*3);
  prevActiveValues=(double *)malloc(sizeof(double)*150);

  int i;
  for (i=0;i<150;i++) {
    mGradient[i] = (double *)malloc(sizeof(double)*3);
    mG[i] = (double *)malloc(sizeof(double)*150);
  }
  for (i=0;i<3;i++) mPDG[i] = (double *)malloc(sizeof(double)*3);

  mActive = (ActiveSet *)malloc(sizeof(ActiveSet));
  testActive = (ActiveSet *)malloc(sizeof(ActiveSet));
  originalActive = (ActiveSet *)malloc(sizeof(ActiveSet));
  MSQ_DBGOUT(1) << "- Executed NonSmoothSteepestDescent::NonSmoothSteepestDescent()\n";
}  
  
  
void NonSmoothSteepestDescent::initialize(PatchData &/*pd*/, MsqError &err)
{
  this->set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err, 1);
  
  // local parameter initialization
  activeEpsilon = .00003;
  //  activeEpsilon = .000000003;
  minAcceptableImprovement = 1e-6;
  minStepSize = 1e-6;
  MSQ_DBGOUT(1) << "- Executed NonSmoothSteepestDescent::initialize()\n";
}

void NonSmoothSteepestDescent::initialize_mesh_iteration(PatchData &/*pd*/,
                                                         MsqError &/*err*/)
{
}
void NonSmoothSteepestDescent::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
  MSQ_FUNCTION_TIMER( "NonSmoothSteepestDescent" );

  //  cout << "- Executing NonSmoothSteepestDescent::optimize_node_positions()\n";
  /* perform the min max smoothing algorithm */
  MSQ_PRINT(2)("\nInitializing the patch iteration\n");

  numVertices = pd.num_vertices();
  MSQ_PRINT(3)("Number of Vertices: %d\n",numVertices);
  numElements = pd.num_elements();
  MSQ_PRINT(3)("Number of Elements: %d\n",numElements);
    //Michael: Note: is this a reliable way to get the dimension?
  mDimension = get_mesh_set()->space_dim();
  MSQ_PRINT(3)("Spatial Dimension: %d\n",mDimension);

  numFree=pd.num_free_vertices(err); MSQ_ERRRTN(err);
  MSQ_PRINT(3)("Num Free = %d\n",numFree);

  MsqFreeVertexIndexIterator free_iter(&pd, err); MSQ_ERRRTN(err);
  free_iter.reset();
  free_iter.next(); 
  freeVertexIndex = free_iter.value();
  MSQ_PRINT(3)("Free Vertex Index = %d\n",freeVertexIndex);

  mCoords = pd.get_vertex_array(err); MSQ_ERRRTN(err);

  if (MSQ_DBG(3)) {
    for (int i99=0;i99<numVertices;i99++) {
      MSQ_PRINT(3)("coords: %g %g\n",mCoords[i99][0],mCoords[i99][1]);
    }
  }

  mConnectivity = pd.get_element_array(err); MSQ_ERRRTN(err);
  if (MSQ_DBG(3)) {
    msq_std::vector<size_t> indices;
    for (int i99=0;i99<numElements;i99++) {
      mConnectivity[i99].get_vertex_indices(indices);
      MSQ_PRINT(3)("connectivity: %d %d %d\n",indices[0],
	      indices[1],indices[2]);
    }
  }

  // TODO - need to switch to validity via metric evaluations should
  // be associated with the compute_function somehow
  /* check for an invalid mesh; if it's invalid return and ask the user 
     to use untangle */
  if (this->validity_check(err)!=1) {
      MSQ_PRINT(1)("ERROR: Invalid mesh\n");
      MSQ_SETERR(err)("Invalid Mesh: Use untangle to create a valid "
                      "triangulation", MsqError::INVALID_MESH);
      return;
  }

  /* assumes one function value per element */
  // TODO - need to include vertex metrics 
  numFunctionValues = numElements;

  /* initialize the optimization data up to numFunctionValues */
  this->init_opt(err); MSQ_ERRRTN(err);
  this->init_max_step_length(err);  MSQ_ERRRTN(err);
  MSQ_PRINT(3)("Done initializing optimization\n");

  /* compute the initial function values */
  //TODO this should return a bool with the validity
  this->compute_function(&pd, originalFunction, err);  MSQ_ERRRTN(err);
 
  // find the initial active set
  this->find_active_set(originalFunction, mActive, err);   MSQ_ERRRTN(err);

  this->minmax_opt(pd,err);  MSQ_ERRRTN(err);
}


void NonSmoothSteepestDescent::terminate_mesh_iteration(PatchData &/*pd*/,
                                                        MsqError &/*err*/)
{
}
  
void NonSmoothSteepestDescent::cleanup()
{
  MSQ_DBGOUT(1) << "- Executing NonSmoothSteepestDescent::cleanup()\n";
  int i;
  for (i=0;i<150;i++) {
    free(mGradient[i]);
    free(mG[i]);
  }
  for (i=0;i<3;i++) free(mPDG[i]);

  free(mFunction);
  free(testFunction);
  free(originalFunction);
  free(mGradient);
  free(mGS);
  free(mG);
  free(mPDG);
  free(prevActiveValues);
  free(mActive);
  free(testActive);
  free(originalActive);
  MSQ_DBGOUT(1) << "- Done with NonSmoothSteepestDescent::cleanup()\n";
}



int NonSmoothSteepestDescent::improvement_check(MsqError &/*err*/)
{
  int improved = 1;
  
  /* check to see that the mesh didn't get worse */
  if (originalValue < mActive->true_active_value) {
     MSQ_PRINT(2)("The local mesh got worse; initial value %f; final value %f\n",
	       originalValue,  mActive->true_active_value );
       improved = 0;
   }

  return(improved);

}




void NonSmoothSteepestDescent::find_plane_points(int dir1, int dir2,
                                                 double **vec, int num_vec,
                                                 double *pt1, double *pt2,
                                                 double* /*pt3*/, int *status,
                                                 MsqError &/*err*/)
{
    int i;
    int ind[50], num_min, num_max;
    int rotate=MSQ_CW;
    int num_rotated=0;
    double pt_1, pt_2;
    double min, inv_slope;
    double min_inv_slope=0.;
    double max; 
    double max_inv_slope=0;
    double inv_origin_slope=0;

    *status = MSQ_CHECK_BOTTOM_UP;
    /* find the minimum points in dir1 starting at -1 */
    num_min = 0; ind[0]=-1; ind[1]=-1; ind[2]=-1; min=1.0;
    for (i=0;i<num_vec;i++) {
      if (vec[i][dir1]<min) {
	min = vec[i][dir1]; ind[0] = i; num_min = 1;
      } else if (fabs(vec[i][dir1] - min) < MSQ_MACHINE_EPS) {
	ind[num_min++] = i;
      }
    }
    if (min >= 0) *status = MSQ_NO_EQUIL;
 
    if (*status != MSQ_NO_EQUIL) {
      switch(num_min) {
      case 1: /* rotate to find the next point */
	MSQ_COPY_VECTOR(pt1,vec[ind[0]],3);
	pt_1 = pt1[dir1]; pt_2 = pt1[dir2];
	if (pt1[dir2] <= 0){rotate=MSQ_CCW; max_inv_slope=MSQ_BIG_NEG_NMBR;}
	if (pt1[dir2] > 0){rotate=MSQ_CW; min_inv_slope=MSQ_BIG_POS_NMBR;}
	switch(rotate) {
	case MSQ_CCW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if ((inv_slope>max_inv_slope) &&  
		  (fabs(inv_slope - max_inv_slope) > MSQ_MACHINE_EPS)) {
		ind[1] = i; max_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - max_inv_slope) < MSQ_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	  break;
	case MSQ_CW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if ((inv_slope<min_inv_slope) && 
		  (fabs(inv_slope - max_inv_slope) > MSQ_MACHINE_EPS)){
		ind[1] = i; min_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - min_inv_slope) < MSQ_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	}
	switch(num_rotated) {
	case 0:
	  MSQ_PRINT(3)("No points in the rotation ... odd\n");
	    *status = MSQ_HULL_TEST_ERROR;
	  break;
	case 1:
	  MSQ_PRINT(3)("Found a line in the convex hull\n");
	  MSQ_COPY_VECTOR(pt2,vec[ind[1]],3); *status = MSQ_TWO_PT_PLANE;
	  break;
	default:
	  MSQ_PRINT(3)("Found 2 or more points in the rotation\n");
	    if (fabs(pt_1) > MSQ_MACHINE_EPS) inv_origin_slope = pt_2/pt_1;
	  switch(rotate) {
	  case MSQ_CCW:
	    if (inv_origin_slope >= max_inv_slope) *status=MSQ_NO_EQUIL;
	    else *status=MSQ_CHECK_TOP_DOWN;
	    break;
	  case MSQ_CW:
	    if (inv_origin_slope <= min_inv_slope) *status=MSQ_NO_EQUIL;
	    else *status=MSQ_CHECK_TOP_DOWN;
	  }
	}
	break;
      case 2: /* use these two points to define the plane */
	MSQ_PRINT(3)("Found two minimum points to define the plane\n");
                MSQ_COPY_VECTOR(pt1,vec[ind[0]],3);
	MSQ_COPY_VECTOR(pt2,vec[ind[1]],3);
	*status = MSQ_TWO_PT_PLANE;
	break;
      default: /* check to see if all > 0 */
	MSQ_PRINT(3)("Found 3 or more points in min plane %f\n",min);
	  if (vec[ind[0]][dir1] >= 0) *status = MSQ_NO_EQUIL;
	  else *status = MSQ_CHECK_TOP_DOWN;
    }
    }

    /***************************/
    /*  failed to find any information, checking top/down this coord*/
    /***************************/

    if (*status == MSQ_CHECK_TOP_DOWN) {
    /* find the maximum points in dir1 starting at 1 */
    num_max = 0; ind[0]=-1; ind[1]=-1; ind[2]=-1; max=-1.0;
    for (i=0;i<num_vec;i++) {
      if (vec[i][dir1] > max) {
	max = vec[i][dir1]; ind[0] = i; num_max = 1;
      } else if (fabs(vec[i][dir1] - max) < MSQ_MACHINE_EPS) {
	ind[num_max++] = i;
      }
    }
    if (max <= 0) *status = MSQ_NO_EQUIL;
 
    if (*status != MSQ_NO_EQUIL) {
      switch(num_max) {
      case 1: /* rotate to find the next point */
	MSQ_COPY_VECTOR(pt1,vec[ind[0]],3);
	pt_1 = pt1[dir1];  pt_2 = pt1[dir2];
	if (pt1[dir2] < 0){rotate=MSQ_CW; min_inv_slope=1E300;}
	if (pt1[dir2] >= 0){rotate=MSQ_CCW; max_inv_slope=-1E300;}
	switch(rotate) {
	case MSQ_CCW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if (inv_slope>max_inv_slope) {
		ind[1] = i; max_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - max_inv_slope) < MSQ_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	  break;
	case MSQ_CW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if (inv_slope<min_inv_slope) {
		ind[1] = i; min_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - min_inv_slope) < MSQ_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	}
	switch(num_rotated) {
	case 0:
	  MSQ_PRINT(3)("No points in the rotation ... odd\n");
	  *status = MSQ_HULL_TEST_ERROR;
	  break;
	case 1:
	  MSQ_PRINT(3)("Found a line in the convex hull\n");
          MSQ_COPY_VECTOR(pt2,vec[ind[1]],3);
	  *status = MSQ_TWO_PT_PLANE;
	  break;
	default:
	  MSQ_PRINT(3)("Found 2 or more points in the rotation\n");
	    /* check to see if rotation got past origin */
	  inv_origin_slope = pt_2/pt_1;
	  switch(rotate) {
	  case MSQ_CCW:
	    if (inv_origin_slope >= max_inv_slope) *status=MSQ_NO_EQUIL;
	    else if (dir1 == 2) *status=MSQ_CHECK_Y_COORD_DIRECTION;
	    else if (dir1 == 1) *status=MSQ_CHECK_X_COORD_DIRECTION;
	    else *status=MSQ_EQUIL;
	    break;
	  case MSQ_CW:
	    if (inv_origin_slope <= min_inv_slope) *status=MSQ_NO_EQUIL;
	    else if (dir1 == 2) *status=MSQ_CHECK_Y_COORD_DIRECTION;
	    else if (dir1 == 1) *status=MSQ_CHECK_X_COORD_DIRECTION;
	    else *status=MSQ_EQUIL;
	  }
	}
	break;
      case 2: /* use these two points to define the plane */
	MSQ_COPY_VECTOR(pt1,vec[ind[0]],3);
	MSQ_COPY_VECTOR(pt2,vec[ind[1]],3);
	*status = MSQ_TWO_PT_PLANE;
	break;
      default: /* check to see if all > 0 */
	MSQ_PRINT(3)("Found 3 in max plane %f\n",max);
	if (vec[ind[0]][dir1] <= 0) *status = MSQ_NO_EQUIL;
	else if (dir1==2) *status=MSQ_CHECK_Y_COORD_DIRECTION;
	else if (dir1==1) *status=MSQ_CHECK_X_COORD_DIRECTION;
	else *status = MSQ_EQUIL;
      }
    }
  }

}

void NonSmoothSteepestDescent::search_direction(PatchData &/*pd*/,
                                                MsqError &err)
{
   int        i;
   int        viable;
   int        singular;
   double     a, b, c, denom;
   double     **dir;
   double     R0, R1;
   double     **P, *x;
   double     search_mag;

   int num_active = mActive->num_active;

   //TODO This might be o.k. actually - i don't see any dependence
   // on the element geometry here... try it and see if it works.
   // if not, try taking all of the gradients in the active set
   // and let the search direction be the average of those.
//   MSQ_FUNCTION_TIMER( "Search Direction" );

   MSQ_PRINT(2)("\nIn Search Direction\n");
   this->print_active_set(mActive, mFunction, err);  MSQ_ERRRTN(err);
   
   if (num_active==0) {
       MSQ_SETERR(err)("No active values in search",MsqError::INVALID_STATE);
       return;
    }

    switch(num_active) {
    case 1: 
        MSQ_COPY_VECTOR(mSearch,mGradient[mActive->active_ind[0]],mDimension);
        mSteepest = mActive->active_ind[0];
        break;
    case 2:
        /* if there are two active points, move in the direction of the
	   intersection of the planes.  This is the steepest descent
           direction found by analytically solving the QP */
        
        /* set up the active gradient directions */
        this->get_active_directions(mGradient,&dir,err);  
        if (MSQ_CHKERR(err)) { free(dir); return; }

        /* form the grammian */
        this->form_grammian(dir,err); 
        if (MSQ_CHKERR(err)) { free(dir); return; }
        this->form_PD_grammian(err); 
        if (MSQ_CHKERR(err)) { free(dir); return; }

        denom = (mG[0][0] + mG[1][1] - 2*mG[0][1]);
        viable = 1;
        if (fabs(denom) > MSQ_MACHINE_EPS) {
	  /* gradients are LI, move along their intersection */
           b = (mG[0][0] - mG[0][1])/denom;  
           a = 1 - b;
           if ((b < 0) || (b > 1)) viable=0;  /* 0 < b < 1 */
           if (viable) {
             for (i=0;i<mDimension;i++) {
               mSearch[i] = a*dir[0][i] + b*dir[1][i];
             }
           } else {
             /* the gradients are dependent, move along one face */
             MSQ_COPY_VECTOR(mSearch,dir[0],mDimension);
           }
        } else {
	   /* the gradients are dependent, move along one face */
           MSQ_COPY_VECTOR(mSearch,dir[0],mDimension);
        }
        mSteepest = mActive->active_ind[0];

        for (i=0;i<num_active;i++) free(dir[i]);
	free(dir);

        break;
    default:
        /* as in case 2: solve the QP problem to find the steepest
           descent direction.  This can be done analytically - as
           is done in Gill, Murray and Wright 
             for 3 active points in 3 directions - test PD of G
             otherwise we know it's SP SD so search edges and faces */

        /* get the active gradient directions */
        this->get_active_directions(mGradient,&dir,err);  MSQ_ERRRTN(err);

        /* form the entries of the grammian matrix */
        this->form_grammian(dir,err);  
        if (MSQ_CHKERR(err)) { free(dir); return; }
        this->form_PD_grammian(err);  
        if (MSQ_CHKERR(err)) { free(dir); return; }

        switch(mDimension) {
        case 2:
  	    this->search_edges_faces(dir,err);  
            if (MSQ_CHKERR(err)) { free(dir); return; }
            break;
        case 3:
	  if (num_active == 3) {
              this->singular_test(num_active,mG,&singular,err);   
              if (MSQ_CHKERR(err)) { free(dir); return; }
              if (!singular) {
	        /* form the entries of P=Z^T G Z where Z = [-1...-1; I ] */
                this->form_reduced_matrix(&P,err);   
                if (MSQ_CHKERR(err)) { free(dir); return; }
                /* form  the RHS and solve the system for the coeffs */
                R0 = mG[0][0] - mG[1][0];  R1 = mG[0][0] - mG[2][0];
                this->solve2x2(P[0][0],P[0][1],P[1][0],P[1][1],R0,R1,&x,err);
                if (MSQ_CHKERR(err)) { free(dir); return; }
                if (x!=NULL) {
                	a = 1 - x[0] - x[1];  b = x[0];  c = x[1];
                	for (i=0;i<mDimension;i++) {
                    	  mSearch[i] = a*dir[0][i] + b*dir[1][i] + 
                       	                      c*dir[2][i];
                	}
                        mSteepest = mActive->active_ind[0];
                	for (i=0;i<num_active-1;i++)  free(P[i]);  
                	free(P);  free(x);
                } else { 
                  	this->search_edges_faces(dir, err);
                        if (MSQ_CHKERR(err)) { free(dir); return; }
                }
	      } else {
                 this->search_edges_faces(dir, err); 
                 if (MSQ_CHKERR(err)) { free(dir); return; }
	      }
            } else {
              this->search_edges_faces(dir, err);
              if (MSQ_CHKERR(err)) { free(dir); return; }
            }
            break;
        }
        for (i=0;i<num_active;i++) free(dir[i]);
	free(dir);
    }

    /* if the search direction is essentially zero, equilibrium pt */
    MSQ_DOT(search_mag,mSearch,mSearch,mDimension);
    MSQ_PRINT(3)("  Search Magnitude %g \n",search_mag);

    if (fabs(search_mag)<1E-13) optStatus = MSQ_ZERO_SEARCH;
    else MSQ_NORMALIZE(mSearch,mDimension);
    MSQ_PRINT(3)("  Search Direction %g %g  Steepest %d\n",mSearch[0],mSearch[1],mSteepest);
}

void NonSmoothSteepestDescent::minmax_opt(PatchData &pd, MsqError &err)
{
//      int valid;
      MSQ_FUNCTION_TIMER( "Minmax Opt" );
      MSQ_PRINT(2)("In minmax_opt\n");

      MSQ_COPY_VECTOR(mFunction,originalFunction,numFunctionValues);
      originalValue = mActive->true_active_value;

      iterCount = 0;
      optIterCount = 0;

      MSQ_PRINT(3)("Done copying original function to function\n");

      this->find_active_set(mFunction, mActive, err); MSQ_ERRRTN(err);
      prevActiveValues[0] = mActive->true_active_value;

     /* check for equilibrium point */
     /* compute the gradient */
     mGradient = this->compute_gradient(&pd, err); MSQ_ERRRTN(err);
     
     if (mActive->num_active >= 2) {
	MSQ_PRINT(3)("Testing for an equilibrium point \n");
	this->check_equilibrium(&equilibriumPt, &optStatus, err); MSQ_ERRRTN(err);

	if (MSQ_DBG(2) && equilibriumPt ) 
	  MSQ_PRINT(2)("Optimization Exiting: An equilibrium point \n");
     }

    /* terminate if we have found an equilibrium point or if the step is
       too small to be worthwhile continuing */
    while ((optStatus != MSQ_EQUILIBRIUM) && 
	   (optStatus != MSQ_STEP_TOO_SMALL) &&
	   (optStatus != MSQ_IMP_TOO_SMALL) &&
	   (optStatus != MSQ_FLAT_NO_IMP) &&
           (optStatus != MSQ_ZERO_SEARCH) &&
	   (optStatus != MSQ_MAX_ITER_EXCEEDED)) {

	/* increase the iteration count by one */
        /* smooth_param->iter_count += 1; */
        iterCount += 1;
        optIterCount += 1;
        if (iterCount > MSQ_MAX_OPT_ITER) optStatus = MSQ_MAX_ITER_EXCEEDED;

	MSQ_PRINT(3)("\nITERATION %d \n",iterCount);
	    
	/* compute the gradient */
	mGradient = this->compute_gradient(&pd, err); MSQ_ERRRTN(err);
        
	MSQ_PRINT(3)("Computing the search direction \n");
	this->search_direction(pd, err); MSQ_ERRRTN(err);

	/* if there are viable directions to search */
	if ((optStatus != MSQ_ZERO_SEARCH) &&
            (optStatus != MSQ_MAX_ITER_EXCEEDED)) {

	    MSQ_PRINT(3)("Computing the projections of the gradients \n");
	    this->get_gradient_projections(err); MSQ_ERRRTN(err);

	    MSQ_PRINT(3)("Computing the initial step size \n");
	    this->compute_alpha(err); MSQ_ERRRTN(err);

	    MSQ_PRINT(3)("Testing whether to accept this step \n");
	    this->step_acceptance(pd, err); MSQ_ERRRTN(err);
            MSQ_PRINT(3)("The new free vertex position is %f %f %f\n",
              mCoords[freeVertexIndex][0],mCoords[freeVertexIndex][1],mCoords[freeVertexIndex][2]);

	    if (MSQ_DBG(3)) {
     		/* Print the active set */
	     	this->print_active_set(mActive, mFunction, err);
                MSQ_ERRRTN(err);
	    }

	    /* check for equilibrium point */
	    if (mActive->num_active >= 2) {
		MSQ_PRINT(3)("Testing for an equilibrium point \n");
                this->check_equilibrium(&equilibriumPt, &optStatus, err); 
                       MSQ_ERRRTN(err);

		if (MSQ_DBG(2) && equilibriumPt) 
		    MSQ_PRINT(2)("Optimization Exiting: An equilibrium point \n");
	    }

	    /* record the values */
	    prevActiveValues[iterCount] = mActive->true_active_value;

	} else {
	    /* decrease the iteration count by one */
	    /* smooth_param->iter_count -= 1; */
	    iterCount -= 1;
	    if (MSQ_DBG(2)) {
		MSQ_PRINT(2)("Optimization Exiting: No viable directions; equilibrium point \n");
		/* Print the old active set */
		this->print_active_set(mActive,mFunction,err); MSQ_ERRRTN(err);
	    }
	}
      }

      MSQ_PRINT(2)("Checking the validity of the mesh\n");
      if (!this->validity_check(err)) MSQ_PRINT(2)("The final mesh is not valid\n");
      MSQ_ERRRTN(err);

      MSQ_PRINT(2)("Number of optimization iterations %d\n", iterCount);
 
      switch(optStatus) {
	case MSQ_EQUILIBRIUM:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Equilibrium\n"); break;
	case MSQ_STEP_TOO_SMALL:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Step Too Small\n"); break;
	case MSQ_IMP_TOO_SMALL:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Improvement Too Small\n"); break;
	case MSQ_FLAT_NO_IMP:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Flat No Improvement\n"); break;
	case MSQ_ZERO_SEARCH:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Zero Search\n"); break;
	case MSQ_MAX_ITER_EXCEEDED:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Max Iter Exceeded\n"); break;
      }
}


void NonSmoothSteepestDescent::step_acceptance(PatchData &pd, MsqError &err)
{
//  int        ierr;
  int        i;
  int        num_values, num_steps;
  int        valid = 1, step_status;
  int        accept_alpha;
  double     estimated_improvement;
  double     current_improvement = 1E300;
  double     previous_improvement = 1E300;
  double     current_percent_diff = 1E300;
  double     original_point[3];

//  MSQ_FUNCTION_TIMER( "Step Acceptance" );
  num_values = numFunctionValues;

  step_status = MSQ_STEP_NOT_DONE;
  optStatus = 0;
  num_steps = 0;

  if (mAlpha < minStepSize) {
      optStatus = MSQ_IMP_TOO_SMALL;
      step_status = MSQ_STEP_DONE;
      MSQ_PRINT(3)("Alpha starts too small, no improvement\n");
  }

  /* save the original function and active set */
  MSQ_COPY_VECTOR(original_point,mCoords[freeVertexIndex],mDimension);
  MSQ_COPY_VECTOR(originalFunction, mFunction, num_values);
  this->copy_active(mActive, originalActive, err); MSQ_ERRRTN(err);

  while (step_status == MSQ_STEP_NOT_DONE) {

    num_steps++;  if (num_steps >= 100) step_status = MSQ_STEP_DONE;

    accept_alpha = MSQ_FALSE;

    while (!accept_alpha && mAlpha>minStepSize) {

      /* make the step */
      for (i=0;i<mDimension;i++) {
         mCoords[freeVertexIndex][i] -= mAlpha*mSearch[i];
      }
        //pd.set_coords_array_element(mCoords[freeVertexIndex],0,err);

      MSQ_PRINT(2)("search direction %f %f \n",mSearch[0],mSearch[1]); 
      MSQ_PRINT(2)("new vertex position %f %f \n",mCoords[freeVertexIndex][0],mCoords[freeVertexIndex][1]); 

      /* assume alpha is acceptable */
      accept_alpha=MSQ_TRUE;

      /* never take a step that makes a valid mesh invalid or worsens the quality */
      // TODO Validity check revision -- do the compute function up here
      // and then the rest based on validity
      valid = validity_check(err); MSQ_ERRRTN(err);
      if (valid) valid=improvement_check(err); MSQ_ERRRTN(err);
      if (!valid) {
          accept_alpha=MSQ_FALSE;
          for (i=0;i<mDimension;i++) {
             mCoords[freeVertexIndex][i] += mAlpha*mSearch[i];
          }
            //pd.set_coords_array_element(mCoords[freeVertexIndex],0,err);
          mAlpha = mAlpha/2;
          MSQ_PRINT(2)("Step not accepted, the new alpha %f\n",mAlpha); 

          if (mAlpha < minStepSize) {
 	        optStatus = MSQ_STEP_TOO_SMALL;
                step_status = MSQ_STEP_DONE;
                MSQ_PRINT(2)("Step too small\n");
 	        /* get back the original point, mFunction, and active set */
                MSQ_COPY_VECTOR(mCoords[freeVertexIndex],original_point,mDimension);
                  //pd.set_coords_array_element(mCoords[freeVertexIndex],0,err);
	        MSQ_COPY_VECTOR(mFunction,originalFunction,num_values);
	        this->copy_active(originalActive, mActive, err); 
	  }
       }
    } 
         
    if (valid  && (mAlpha > minStepSize)) {
      /* compute the new function and active set */
      this->compute_function(&pd, mFunction, err); MSQ_ERRRTN(err);
      this->find_active_set(mFunction, mActive, err); MSQ_ERRRTN(err);
	
      /* estimate the minimum improvement by taking this step */
      this->get_min_estimate(&estimated_improvement, err); MSQ_ERRRTN(err);
      MSQ_PRINT(2)("The estimated improvement for this step: %f\n",
		   estimated_improvement); 
	
      /* calculate the actual increase */
      current_improvement = mActive->true_active_value - prevActiveValues[iterCount-1];

      MSQ_PRINT(3)("Actual improvement %f\n",current_improvement);

      /* calculate the percent difference from estimated increase */
      current_percent_diff = fabs(current_improvement-estimated_improvement)/
	fabs(estimated_improvement);

      /* determine whether to accept a step */
      if ((fabs(previous_improvement) > fabs(current_improvement)) && 
	  (previous_improvement < 0)) {
	/* accept the previous step - it was better */
	     MSQ_PRINT(2)("Accepting the previous step\n");
 
	/* subtract alpha in again (previous step) */
	for (i=0;i<mDimension;i++) {
	  mCoords[freeVertexIndex][i] -= mAlpha*mSearch[i];
	}
            //pd.set_coords_array_element(mCoords[freeVertexIndex],0,err);

	/* does this make an invalid mesh valid? */
   //TODO Validity check revisison
        valid = 1;
        valid = validity_check(err); MSQ_ERRRTN(err);
        if (valid) valid=improvement_check(err); MSQ_ERRRTN(err);

	/* copy test function and active set */
	MSQ_COPY_VECTOR(mFunction,testFunction,numFunctionValues);
	this->copy_active(testActive, mActive, err); MSQ_ERRRTN(err);
 
	optStatus = MSQ_STEP_ACCEPTED;  step_status = MSQ_STEP_DONE;
            
	/* check to see that we're still making good improvements */
	if (fabs(previous_improvement) < minAcceptableImprovement) {
	  optStatus = MSQ_IMP_TOO_SMALL; step_status = MSQ_STEP_DONE;
	  MSQ_PRINT(2)("Optimization Exiting: Improvement too small\n");
	}

      } else if (((fabs(current_improvement) > fabs(estimated_improvement)) ||
		  (current_percent_diff < .1)) && (current_improvement<0)) {
	/* accept this step, exceeded estimated increase or was close */
	optStatus = MSQ_STEP_ACCEPTED;  step_status = MSQ_STEP_DONE;

	/* check to see that we're still making good improvements */
	if (fabs(current_improvement) < minAcceptableImprovement) {
	  MSQ_PRINT(2)("Optimization Exiting: Improvement too small\n");
	  optStatus = MSQ_IMP_TOO_SMALL; step_status = MSQ_STEP_DONE;
	}

      } else if ((current_improvement > 0) && (previous_improvement > 0) &&
		 (fabs(current_improvement) < minAcceptableImprovement) &&
		 (fabs(previous_improvement) < minAcceptableImprovement)) {

	/* we are making no progress, quit */
	optStatus = MSQ_FLAT_NO_IMP; step_status = MSQ_STEP_DONE;
	MSQ_PRINT(2)("Opimization Exiting: Flat no improvement\n");
           
	/* get back the original point, function, and active set */
	MSQ_COPY_VECTOR(mCoords[freeVertexIndex],original_point,mDimension);
            //pd.set_coords_array_element(mCoords[freeVertexIndex],0,err);
	MSQ_COPY_VECTOR(mFunction,originalFunction,numFunctionValues);
	this->copy_active(originalActive, mActive, err); MSQ_CHKERR(err);

      }
      else
      {
	/* halve alpha and try again */
	/* add out the old step */
	for (i=0;i<mDimension;i++) mCoords[freeVertexIndex][i] += mAlpha*mSearch[i];
            //pd.set_coords_array_element(mCoords[freeVertexIndex],0,err);

	/* halve step size */
	mAlpha = mAlpha/2; 
	MSQ_PRINT(3)("Step not accepted, the new alpha %f\n",mAlpha);

	if (mAlpha < minStepSize)
          {
	  /* get back the original point, function, and active set */
	  MSQ_PRINT(2)("Optimization Exiting: Step too small\n");
	  MSQ_COPY_VECTOR(mCoords[freeVertexIndex],original_point,mDimension);
              //pd.set_coords_array_element(mCoords[freeVertexIndex],0,err);
	  MSQ_COPY_VECTOR(mFunction,originalFunction,numFunctionValues);
	  this->copy_active(originalActive, mActive, err); MSQ_ERRRTN(err);
	  optStatus = MSQ_STEP_TOO_SMALL;  step_status = MSQ_STEP_DONE;
	}
          else
          {
	  MSQ_COPY_VECTOR(testFunction, mFunction, numFunctionValues);
	  this->copy_active(mActive, testActive, err); MSQ_ERRRTN(err);
	  previous_improvement = current_improvement;
	}
      }
    }
  }
  if (current_improvement>0 && optStatus==MSQ_STEP_ACCEPTED) {
    MSQ_PRINT(2)("Accepted a negative step %f \n",current_improvement);
  }

}
