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
/*!
  \file   NonSmoothSteepestDescent.hpp
  \brief  

  The NonSmoothSteepestDescent Class implements the steepest descent algorythm in
  order to move a free vertex to an optimal position given an
  ObjectiveFunction object and a QaulityMetric object.

  \author Thomas Leurent
  \date   2002-06-13
*/

#ifndef Mesquite_NonSmoothSteepestDescent_hpp 
#define Mesquite_NonSmoothSteepestDescent_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "ObjectiveFunction.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqDebug.hpp"

namespace Mesquite
{
#define MSQ_XDIR 0
#define MSQ_YDIR 1
#define MSQ_ZDIR 2
#define MSQ_BIG_POS_NMBR    1E300
#define MSQ_BIG_NEG_NMBR   -1E300
#define MSQ_MAX_OPT_ITER    20

#define MSQ_CCW               1
#define MSQ_CW                0
#define MSQ_NO_EQUIL 101
#define MSQ_CHECK_TOP_DOWN 102
#define MSQ_CHECK_BOTTOM_UP 103
#define MSQ_TWO_PT_PLANE 104
#define MSQ_THREE_PT_PLANE 105
#define MSQ_CHECK_Y_COORD_DIRECTION 106
#define MSQ_CHECK_X_COORD_DIRECTION 107
#define MSQ_CHECK_Z_COORD_DIRECTION 108
#define MSQ_EQUIL 109
#define MSQ_HULL_TEST_ERROR 110

#define MSQ_STEP_ACCEPTED     100
#define MSQ_IMP_TOO_SMALL     101
#define MSQ_FLAT_NO_IMP       102
#define MSQ_STEP_TOO_SMALL    103
#define MSQ_EQUILIBRIUM       104
#define MSQ_ZERO_SEARCH       105
#define MSQ_MAX_ITER_EXCEEDED 106

#define MSQ_STEP_DONE        101
#define MSQ_STEP_NOT_DONE    102

#define MAX_NUM_ELEMENTS 150
#define MAX_FUNC_PER_ELEMENT 6
#define MSQ_MACHINE_EPS     1E-16
#define MSQ_TRUE  1
#define MSQ_FALSE 0
#define MSQ_MAX(a,b) (a > b ? a : b)
#define MSQ_MIN(a,b) (a < b ? a : b)
#define MSQ_LESS_THAN_MACHINE_EPS(x)   ( ((fabs(x)+1.0) > 1.0) ? 0 : 1 )

#define MSQ_DOT(c,a,b,n) {\
  int i99; \
  if (n==2) c = a[0]*b[0] + a[1]*b[1]; \
  else if (n==3) c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];\
  else { \
    c = 0; \
    for (i99=0;i99<n;i99++) c += a[i99]*b[i99]; \
  } \
}

#define MSQ_NORMALIZE(v,n) {\
    int i99; \
    double mag99; \
    if (n==2){ \
       mag99 = sqrt(v[0]*v[0] + v[1]*v[1]) ; \
       if (mag99 != 0) { \
          v[0] = v[0]/mag99; \
          v[1] = v[1]/mag99; \
       } \
    } else if (n==3) {\
     mag99 = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) ; \
     if (mag99 != 0) { \
         v[0] = v[0]/mag99; \
         v[1] = v[1]/mag99; \
         v[2] = v[2]/mag99; \
     } \
   } else { \
     mag99 = 0; \
     for (i99=0;i99<n;i99++) mag99+=v[i99]+v[i99]; \
     if (mag99 != 0) { \
       for (i99=0;i99<n;i99++) v[i99] = v[i99]/mag99;\
     } \
   }\
}

#define MSQ_COPY_VECTOR(a,b,n) { \
  int i99; \
  if (n==2) { \
     a[0] = b[0];  a[1] = b[1];  \
  } else if (n==3) {\
     a[0] = b[0];  a[1] = b[1];  a[2] = b[2]; \
  } else { \
     for (i99=0;i99<n;i99++) a[i99] = b[i99]; \
  } \
}


  struct ActiveSet
  {
    int num_active;
    int num_equal;
    double true_active_value;
    int active_ind[150];  // need a better way of setting max number of active values
  };

  class NonSmoothSteepestDescent : public VertexMover 
  {
  public:
    NonSmoothSteepestDescent(ObjectiveFunction* of);

    virtual ~NonSmoothSteepestDescent() { }
    
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();
    
  private:
      /* local copy of patch data */
      //    PatchData patch_data;
    int mDimension;
    int numVertices;
    int numElements;
    MsqVertex* mCoords;
    MsqMeshEntity* mConnectivity;
    int numFree;
    int freeVertexIndex;
    
      /* smoothing parameters */
    double activeEpsilon;
    double minAcceptableImprovement;
    double minStepSize;
    
      /* optimization data */
    double originalValue;
    int iterCount;
    int optIterCount;
    int numFunctionValues;
    double *mFunction;
    double *testFunction;
    double *originalFunction;
    double **mGradient;
    int optStatus;
    int equilibriumPt;
    int mSteepest;
    double mSearch[3];
    double mAlpha;
    double maxAlpha;
    double *mGS;
    double *prevActiveValues;
    double **mG;
    double **mPDG;
    int pdgInd[3];
    ActiveSet *mActive;
    ActiveSet *testActive;
    ActiveSet *originalActive;
    
      /* functions */
    void init_opt(MsqError &err);
    void init_max_step_length(MsqError &err);
    
      /* optimize */
    void minmax_opt(PatchData &pd, MsqError &err);
    void step_acceptance(PatchData &pd, MsqError &err); 
    void get_min_estimate(double *final_est, MsqError &err);
    void get_gradient_projections(MsqError &err);
    void compute_alpha(MsqError &err);
    void copy_active(ActiveSet *active1, ActiveSet *active2, MsqError &err);
    
      /* function/gradient/active set computations */
    bool compute_function(PatchData *pd, double *function, MsqError &err);
    double** compute_gradient(PatchData *pd, MsqError &err);
    void find_active_set(double *function, ActiveSet *active_set, MsqError &err);
    void print_active_set(ActiveSet *active_set, double *func, MsqError &err);
    
      /* checking validity/improvement */
    int improvement_check(MsqError &err);
    int validity_check(MsqError &err);
    
      /* checking equilibrium routines */
    void check_equilibrium(int *equil, int *opt_status, MsqError &err);
    int convex_hull_test(double **vec, int num_vec, MsqError &err);
    int check_vector_dots(double **vec, int num_vec, double *normal, MsqError &err);
    void find_plane_normal(double pt1[3], double pt2[3], double pt3[3], 
                           double *cross, MsqError &err);
    void find_plane_points(int dir1, int dir2, double **vec, int num_vec, double *pt1,
                           double *pt2, double*pt3, int *opt_status, MsqError &err);
    
      /* from the matrix file */
    void form_grammian(double **vec, MsqError &err);
    void form_PD_grammian(MsqError &err);
    void singular_test(int n, double **A, int *singular, MsqError &err);
    void condition3x3(double **A, double *cond, MsqError &err);
    void solve2x2(double a11, double a12, double a21, double a22, 
                  double b1, double b2, double **x,MsqError &err);
    void form_reduced_matrix(double ***P, MsqError &err);
    
      /* search direction */
    void search_direction(PatchData &pd, MsqError &err);
    void search_edges_faces(double **dir, MsqError &err);
    void get_active_directions(double **gradient,
                               double ***dir, MsqError &err);
  };
  

inline bool NonSmoothSteepestDescent::compute_function(PatchData *patch_data, double *func, MsqError &err)
{
  // ASSUMES ONE VALUE PER ELEMENT; ALSO NEED 1.0/FUNCTION WHICH IS ONLY
  // TRUE OF CONDITION NUMBER

  //  MSQ_DEBUG_PRINT(2,"Computing Function\n");
//  FUNCTION_TIMER_START("Compute Function");

  //TODO need to switch this to element or vertex metric evaluations
  //TODO need to include boolean testing for validity
  int i;
  bool valid_bool=true;

  for (i=0;i<numElements;i++) func[i]=0.0;
  QualityMetric* currentQM=objFunc->get_quality_metric();
  if(currentQM==NULL){
    currentQM = objFunc->get_quality_metric_list().front();
  }
  
  for (i=0;i<numElements;i++) {
    valid_bool = currentQM->evaluate_element(*patch_data,
                                &(patch_data->element_by_index(i)),
                                func[i], err); MSQ_ERRZERO(err);
    //    MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  Function value[%d]=%g\n",i,func[i]);});
  }
//  FUNCTION_TIMER_END();
  return(valid_bool);
}


inline double** NonSmoothSteepestDescent::compute_gradient(PatchData *patch_data, MsqError &err)
{
//  FUNCTION_TIMER_START("Compute Gradient");

  MSQ_DBGOUT(2) << "Computing Gradient\n";

  double delta = 10e-6;

  for (int i=0;i<numElements;i++) {
    for (int j=0;j<3;j++) mGradient[i][j] = 0.0;
  }
  QualityMetric* currentQM=objFunc->get_quality_metric();
  if(currentQM==NULL)
    currentQM = objFunc->get_quality_metric_list().front();

  double *func, *fdelta;
  func = (double *)malloc(sizeof(double)*150);
  fdelta = (double *)malloc(sizeof(double)*150);

  this->compute_function(patch_data, func, err); 
  if (MSQ_CHKERR(err)) {
    free(func);
    free(fdelta);
    return 0;
  }

  /* gradient in the x, y, z direction */

  for (int j=0;j<3;j++) {

    // perturb the coordinates of the free vertex in the j direction by delta
    mCoords[freeVertexIndex][j] += delta;

    //compute the function at the perturbed point location
    this->compute_function(patch_data, fdelta, err);  
    if (MSQ_CHKERR(err)) {
      free(func);
      free(fdelta);
      return 0;
    }

    //compute the numerical gradient
    for (int i=0;i<numFunctionValues;i++) {
       mGradient[i][j] = (fdelta[i] - func[i])/delta;
       // MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  Gradient value[%d][%d]=%g\n",i,j,mGradient[i][j]);});
    }

    // put the coordinates back where they belong
    mCoords[freeVertexIndex][j] -= delta;
  }
  
  free(func);
  free(fdelta);
//  FUNCTION_TIMER_END();
  return(mGradient);
}

inline void NonSmoothSteepestDescent::find_active_set(double *function,
                                                      ActiveSet *active_set,
                                                      MsqError & /*err*/ )
{ 
    int         i, ind;
    double      function_val;
    double      active_value0;
    double      temp;

//    FUNCTION_TIMER_START("Find Active Set");
    MSQ_DBGOUT(2) << "\nFinding the active set\n";

    // initialize the active set indices to zero
    for (i=0;i<numFunctionValues;i++) active_set->active_ind[i] = 0; 

    /* the first function value is our initial active value */
    active_set->num_active = 1;
    active_set->num_equal = 0;
    active_set->active_ind[0] = 0;
    active_set->true_active_value = function[0];
    //    MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  function value[0]: %g\n",function[0]);});

    /* first sort out the active set... 
       all vals within active_epsilon of largest val */

    for (i=1;i<numFunctionValues;i++) {
	function_val = function[i];
        active_set->true_active_value = MSQ_MAX(function_val,active_set->true_active_value);
	active_value0 = function[active_set->active_ind[0]];
	temp = fabs(function_val - active_value0);
	//        MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  function value[%d]: %g\n",i,function[i]);});
	if ( function_val > active_value0 ) {
	    if ( temp > activeEpsilon) {
		active_set->num_active = 1;
		active_set->num_equal = 0;
		active_set->active_ind[0] = i;
	    } else if ( temp < activeEpsilon) {
		active_set->num_active += 1;
		ind = active_set->num_active - 1;
		active_set->active_ind[ind] = i;
		if (fabs(function_val - active_value0) < MSQ_MACHINE_EPS) {
		    active_set->num_equal += 1;
		}
	    }
	} else {
	    if (temp < activeEpsilon) {
		active_set->num_active += 1;
		ind = active_set->num_active - 1;
		active_set->active_ind[ind] = i;
		if (fabs(function_val - active_value0) < MSQ_MACHINE_EPS) {
		    active_set->num_equal += 1;
		}
	    }
	}
    }

}


inline int NonSmoothSteepestDescent::validity_check(MsqError &/*err*/)
        
{
//  FUNCTION_TIMER_START("Validity Check");

  // ONLY FOR SIMPLICIAL MESHES - THERE SHOULD BE A VALIDITY CHECKER ASSOCIATED
  // WITH MSQ ELEMENTS
  
  /* check that the simplicial mesh is still valid, based on right handedness. 
       Returns a 1 or a 0 */

  // TODO as a first step we can switch this over to the function
  // evaluation and use the rest of the code as is
  int valid = 1;
  double dEps = 1.e-13;

  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;

  if (mDimension == 2)
  {
    for (int i=0;i<numElements;i++)
    {
      double dummy = 0;
      mCoords[mConnectivity[i].get_vertex_index(0)].get_coordinates(x1, y1, dummy);
      mCoords[mConnectivity[i].get_vertex_index(1)].get_coordinates(x2, y2, dummy);
      mCoords[mConnectivity[i].get_vertex_index(2)].get_coordinates(x3, y3, dummy);
      
      double a = x2*y3 - x3*y2;
      double b = y2 - y3;
      double c = x3 - x2;
      
      if (.5*(a+b*x1+c*y1) < .01*MSQ_MACHINE_EPS)
        valid=0;
    }
  }

  if (mDimension == 3)
  {
    for (int i=0;i<numElements;i++)
    {
      mCoords[mConnectivity[i].get_vertex_index(0)].get_coordinates(x1, y1, z1);
      mCoords[mConnectivity[i].get_vertex_index(1)].get_coordinates(x2, y2, z2);
      mCoords[mConnectivity[i].get_vertex_index(2)].get_coordinates(x3, y3, z3);
      mCoords[mConnectivity[i].get_vertex_index(3)].get_coordinates(x4, y4, z4);
      
      double dDX2 = x2 - x1;
      double dDX3 = x3 - x1;
      double dDX4 = x4 - x1;
      
      double dDY2 = y2 - y1;
      double dDY3 = y3 - y1;
      double dDY4 = y4 - y1;

      double dDZ2 = z2 - z1;
      double dDZ3 = z3 - z1;
      double dDZ4 = z4 - z1;
      
        /* dDet is proportional to the cell volume */
      double dDet = dDX2*dDY3*dDZ4 + dDX3*dDY4*dDZ2 + dDX4*dDY2*dDZ3
        - dDZ2*dDY3*dDX4 - dDZ3*dDY4*dDX2 - dDZ4*dDY2*dDX3 ;

        /* Compute a length scale based on edge lengths. */
      double dScale = ( sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) +
                             (z1-z2)*(z1-z2)) +
                        sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) +
                             (z1-z3)*(z1-z3)) +
                        sqrt((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4) +
                             (z1-z4)*(z1-z4)) +
                        sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) +
                             (z2-z3)*(z2-z3)) +
                        sqrt((x2-x4)*(x2-x4) + (y2-y4)*(y2-y4) +
                             (z2-z4)*(z2-z4)) +
                        sqrt((x3-x4)*(x3-x4) + (y3-y4)*(y3-y4) +
                             (z3-z4)*(z3-z4)) ) / 6.;
      
        /* Use the length scale to get a better idea if the tet is flat or
           just really small. */
      if (fabs(dScale) < MSQ_MACHINE_EPS)
      {
        return(valid = 0);
      }
      else
      {
        dDet /= (dScale*dScale*dScale);
      }
      
      if (dDet > dEps)
      {
        valid = 1;
      }
      else if (dDet < -dEps)
      {
        valid = -1;
      }
      else
      {
        valid = 0;
      }
    }  // end for i=1,numElements
  }  // end mDimension==3
  
  //  MSQ_DEBUG_ACTION(2,{fprintf(stdout,"Mesh Validity is: %d \n",valid);});
  
//  FUNCTION_TIMER_END();
  return(valid);
}


inline void NonSmoothSteepestDescent::get_active_directions(double **gradient, 
                                                            double ***dir,
                                                            MsqError &/*err*/)
{
    int i;
    int num_active = mActive->num_active;

    (*dir) =(double **)malloc(sizeof(double *)*num_active);
    for (i=0;i<num_active;i++) {
        (*dir)[i] =(double *)malloc(sizeof(double)*mDimension);
        MSQ_COPY_VECTOR((*dir)[i],gradient[mActive->active_ind[i]],mDimension);
    }
}


inline int NonSmoothSteepestDescent::check_vector_dots(double **vec,
                                                       int num_vec, 
                                                       double *normal,
                                                       MsqError &/*err*/)
{
    int equil;
    int i, ind;
    double test_dot, dot;

    equil = MSQ_FALSE;
    MSQ_DOT(test_dot,vec[0],normal,3);
    ind = 1;
    while ((fabs(test_dot) < MSQ_MACHINE_EPS) && (ind<num_vec)) {
      MSQ_DOT(test_dot,vec[ind],normal,3);
      ind++;
    }
      
    for (i=ind;i<num_vec;i++) {
       MSQ_DOT(dot,vec[i],normal,3);
       if ( ((dot>0 && test_dot<0) || (dot<0 && test_dot>0)) &&
            (fabs(dot)>MSQ_MACHINE_EPS)) {
          return(equil = MSQ_TRUE);

       }
    }
    return(equil);
}


inline void NonSmoothSteepestDescent::find_plane_normal(double pt1[3],
                                                        double pt2[3],
                                                        double pt3[3],
                                                        double *cross,
                                                        MsqError &/*err*/)
{
  int i;
  double vecA[3], vecB[3];

  for (i=0;i<3;i++) {
    vecA[i] = pt2[i] - pt1[i];
    vecB[i] = pt3[i] - pt1[i];
  }
  cross[0] = vecA[1]*vecB[2] - vecA[2]*vecB[1];
  cross[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2];
  cross[2] = vecA[0]*vecB[1] - vecA[1]*vecB[0];
  MSQ_NORMALIZE(cross, 3);
}


inline int NonSmoothSteepestDescent::convex_hull_test(double **vec, int num_vec, MsqError &err)
{
//    int ierr;
    int equil;
    int status, dir_done;
    double pt1[3], pt2[3], pt3[3];
    double normal[3];

//    FUNCTION_TIMER_START("Convex Hull Test");
    /* tries to determine equilibrium for the 3D case */
    equil = 0;
    status = MSQ_CHECK_Z_COORD_DIRECTION;
    dir_done = -1;

    if (num_vec <= 2) status = MSQ_NO_EQUIL;

    while ((status != MSQ_EQUIL) && (status != MSQ_NO_EQUIL) && 
           (status != MSQ_HULL_TEST_ERROR)) {
       if (status == MSQ_CHECK_Z_COORD_DIRECTION) {
          this->find_plane_points(MSQ_ZDIR, MSQ_YDIR, 
                          vec, num_vec, pt1, pt2, pt3, &status, err);
          dir_done = 2;
       }
       if (status == MSQ_CHECK_Y_COORD_DIRECTION) {
          this->find_plane_points(MSQ_YDIR, MSQ_XDIR, 
                          vec, num_vec, pt1, pt2, pt3, &status, err);
          dir_done = 1;
       }
       if (status == MSQ_CHECK_X_COORD_DIRECTION) {
          this->find_plane_points(MSQ_XDIR, MSQ_ZDIR, 
                          vec, num_vec, pt1, pt2, pt3, &status, err);
          dir_done = 0;
       }
       if (status == MSQ_TWO_PT_PLANE) {
          pt3[0]=0.; pt3[1]=0.; pt3[2]=0.;
       }
       if ((status == MSQ_TWO_PT_PLANE) || (status == MSQ_THREE_PT_PLANE)){
           this->find_plane_normal(pt1,pt2,pt3,normal,err); 
           equil = this->check_vector_dots(vec,num_vec,normal,err); 
           if (equil == 1) {
             switch(dir_done){
             case 2:
               equil = 0; status = MSQ_CHECK_Y_COORD_DIRECTION;
               break;
             case 1:
               equil = 0; status = MSQ_CHECK_X_COORD_DIRECTION;
               break;
             case 0:
               equil = 1; status = MSQ_EQUIL;
             }
           } else if (equil == 0) {
               status = MSQ_NO_EQUIL;
           } else {
              MSQ_SETERR(err)("equil flag not set to 0 or 1",MsqError::INVALID_STATE);
           }
       }
    }
    switch (status){
    case MSQ_NO_EQUIL:
      MSQ_PRINT(3)("Not an equilibrium point\n");
      equil = 0; 
      break;
    case MSQ_EQUIL:
      MSQ_PRINT(3)("An equilibrium point\n");
      equil = 1;
      break;
    default:
      MSQ_PRINT(3)("Failed to determine equil or not; status = %d\n",status);
    }
//    FUNCTION_TIMER_END();
    return (equil);
}

inline void NonSmoothSteepestDescent::form_grammian(double **vec, MsqError &err)
{
   int i, j;
   int num_active = mActive->num_active;

   if (num_active > 150) {
      MSQ_SETERR(err)("Exceeded maximum allowed active values",MsqError::INVALID_STATE);
      return;
   }
   /* form the grammian with the dot products of the gradients */
   for (i=0; i<num_active; i++) {
      for (j=i; j<num_active; j++) {
         mG[i][j] = 0.;
	 MSQ_DOT(mG[i][j],vec[i],vec[j],mDimension);
	 mG[j][i] = mG[i][j];	 
      }
   }
}

inline void NonSmoothSteepestDescent::check_equilibrium(int *equil, int *status, MsqError &err)
{
//    int  ierr;
    int  i,j;
    int  ind1, ind2;  
    double min;
    double **dir;
    double mid_vec[3], mid_cos, test_cos;

    //TODO - this subroutine is no longer clear to me... is it still
    // appropriate for quads and hexes?  I think it might be in 2D, but
    // 3D is less clear.  Is there a more general algorithm to use?
    // ask Todd/check in numerical optimization

    *equil = MSQ_FALSE;
    ind1 = ind2 = -1;

    int num_active = mActive->num_active;

    if (num_active==numFunctionValues)
    {
         *equil = 1; *status = MSQ_EQUILIBRIUM;
         MSQ_PRINT(3)("All the function values are in the active set\n"); 
    }

    /* set up the active mGradient directions */
    this->get_active_directions(mGradient,&dir,err);

    /* normalize the active directions */
    for (j=0;j<num_active;j++) MSQ_NORMALIZE(dir[j],mDimension);

    if (mDimension == 2) {
      /* form the grammian */
      this->form_grammian(dir,err);  

    /* find the minimum element in the upper triangular portion of G*/
    min = 1;
    for (i=0;i<num_active;i++) {
      for (j=i+1;j<num_active;j++) {
        if ( fabs(-1 - mG[i][j]) < 1E-08 ) {
           *equil = 1; *status = MSQ_EQUILIBRIUM;
           MSQ_PRINT(3)("The gradients are antiparallel, eq. pt\n"); 
         }
         if (mG[i][j]  < min) { 
           ind1 = i; ind2 = j;
           min = mG[i][j];
        }
      }
    }

    if ((ind1 != -1) && (ind2 != -1)) {
      /* find the diagonal of the parallelepiped */
      for (j=0;j<mDimension;j++) {
       mid_vec[j]=.5*(dir[ind1][j]+dir[ind2][j]);
      }
      MSQ_NORMALIZE(mid_vec,mDimension);
      MSQ_DOT(mid_cos,dir[ind1],mid_vec,mDimension);

      /* test the other vectors to be sure they lie in the cone */
      for (i=0;i<num_active;i++) {
         if ((i != ind1) && (i != ind2)) {
            MSQ_DOT(test_cos,dir[i],mid_vec,mDimension);
            if ((test_cos < mid_cos)  &&  fabs(test_cos-mid_cos) > MSQ_MACHINE_EPS) {
              MSQ_PRINT(3)("An equilibrium point \n");
              *equil = 1; *status = MSQ_EQUILIBRIUM;
            }
         }
       }
     }
    }
    if (mDimension == 3) {
       *equil = this->convex_hull_test(dir,num_active,err);
       if (*equil == 1) *status = MSQ_EQUILIBRIUM;
    }
}



inline void NonSmoothSteepestDescent::condition3x3(double **A, double *cond,
                                                   MsqError &/*err*/) 
{
//   int ierr;
   double a11, a12, a13;
   double a21, a22, a23;
   double a31, a32, a33;
//   double s1, s2, s4, s3, t0;
   double s1, s2, s3;
   double denom;
//   double one = 1.0;
   double temp;
   int zero_denom = MSQ_TRUE;

   a11 = A[0][0]; a12=A[0][1]; a13=A[0][2];
   a21 = A[1][0]; a22=A[1][1]; a23=A[1][2];
   a31 = A[2][0]; a32=A[2][1]; a33=A[2][2];

   denom = -a11*a22*a33+a11*a23*a32+a21*a12*a33-a21*a13*a32-
            a31*a12*a23+a31*a13*a22;

   if ( (fabs(a11) > MSQ_MACHINE_EPS) && 
        (fabs(denom/a11) > MSQ_MACHINE_EPS)) {
         zero_denom = MSQ_FALSE;
   }
   if ( (fabs(a22) > MSQ_MACHINE_EPS) && 
        (fabs(denom/a22) > MSQ_MACHINE_EPS)) {
         zero_denom = MSQ_FALSE;
   }       
   if ( (fabs(a33) > MSQ_MACHINE_EPS) && 
        (fabs(denom/a33) > MSQ_MACHINE_EPS)) {
         zero_denom = MSQ_FALSE;
   }

   if (zero_denom) {
     (*cond) = 1E300;
   } else {
     s1 = sqrt(a11*a11 + a12*a12 + a13*a13 + 
               a21*a21 + a22*a22 + a23*a23 + 
               a31*a31 + a32*a32 + a33*a33);

     temp = (-a22*a33+a23*a32)/denom;
     s3 = temp*temp;
     temp =(a12*a33-a13*a32)/denom;
     s3 += temp*temp;
     temp = (a12*a23-a13*a22)/denom;
     s3 += temp*temp;
     temp = (a21*a33-a23*a31)/denom;
     s3 += temp*temp;
     temp = (a11*a33-a13*a31)/denom;
     s3 += temp*temp;
     temp = (a11*a23-a13*a21)/denom;
     s3 += temp*temp;
     temp = (a21*a32-a22*a31)/denom;
     s3 += temp*temp;
     temp = (-a11*a32+a12*a31)/denom;
     s3 += temp*temp;
     temp = (-a11*a22+a12*a21)/denom;
     s3 += temp*temp;

     s2 = sqrt(s3);
     (*cond) = s1*s2;
   }
}

inline void NonSmoothSteepestDescent::singular_test(int n, double **A, int *singular, MsqError &err) 
{
//    int test;
//    double determinant;
    double cond;

    if ((n>3) || (n<1)) {
      MSQ_SETERR(err)("Singular test works only for n=1 to n=3",MsqError::INVALID_ARG);
      return;
    }

    (*singular)=MSQ_TRUE;
    switch(n) {
    case 1:
        if (A[0][0] > 0) (*singular) = MSQ_FALSE;
        break;
    case 2:
        if (fabs(A[0][0]*A[1][1] - A[0][1]*A[1][0]) > MSQ_MACHINE_EPS)           
            (*singular) = MSQ_FALSE;
        break;
    case 3:
       /* calculate the condition number */
        this->condition3x3(A, &cond, err); 
        if (cond < 1E14) (*singular)=MSQ_FALSE;
        break;
    }
}


inline void NonSmoothSteepestDescent::form_PD_grammian(MsqError &err)
{
    int  i,j,k;
    int  g_ind_1;
    int  singular;

    int num_active = mActive->num_active;
        
    /* this assumes the grammian has been formed */
    for (i=0;i<num_active;i++) {
      for (j=0;j<num_active;j++) {
        if (mG[i][j]==-1) {
          MSQ_SETERR(err)("Grammian not computed properly",MsqError::INVALID_STATE);
          return;
        }
      }
    }

    /* use the first gradient in the active set */
    g_ind_1 = 0;
    mPDG[0][0] = mG[0][0];
    pdgInd[0] = mActive->active_ind[0];

    /* test the rest and add them as appropriate */
    k = 1; i = 1;
    while( (k<mDimension) && (i < num_active) ) {
        mPDG[0][k] = mPDG[k][0] = mG[0][i];
        mPDG[k][k] = mG[i][i];
        if ( k == 2) { /* add the dot product of g1 and g2 */
           mPDG[1][k] = mPDG[k][1] = mG[g_ind_1][i];
        }
        this->singular_test(k+1,mPDG,&singular,err);
        if (!singular) {
           pdgInd[k] = mActive->active_ind[i];
           if (k==1) g_ind_1 = i;
           k++;
        }
        i++;
    }
}


inline void NonSmoothSteepestDescent::search_edges_faces(double **dir, MsqError &err)
{
//    int ierr;
    int i,j,k;
    int viable;
    double a,b,denom;
    double g_bar[3];
    double temp_search[3];
    double projection, min_projection;

    int num_active = mActive->num_active;

    if ( (mDimension != 2) && (mDimension != 3)) {
       MSQ_SETERR(err)("Dimension must be 2 or 3", MsqError::INVALID_MESH);
    }

    /* initialize the search direction to 0,0 */
    for (i=0;i<mDimension;i++) temp_search[i] = 0;

    /* Check for viable faces */
    min_projection = 1E300;
    for (i=0; i<num_active; i++) {
        /* FACE I */
        viable = 1;

        /* test the viability */
        for (j=0;j<num_active;j++) {       /* lagrange multipliers>0 */
             if (mG[j][i] < 0) viable = 0;
        }
       
        /* find the minimum of viable directions */
        if ((viable) && (mG[i][i] < min_projection)) {
            min_projection = mG[i][i];
            MSQ_COPY_VECTOR(temp_search,dir[i],mDimension);
            mSteepest = mActive->active_ind[i];
        }
    
       /* INTERSECTION IJ */
       for (j=i+1; j<num_active; j++) {
          viable = 1;

          /* find the coefficients of the intersection 
             and test the viability */
          denom = 2*mG[i][j] - mG[i][i] - mG[j][j];
          a = b = 0;
          if (fabs(denom) > MSQ_MACHINE_EPS) {
             b = (mG[i][j] - mG[i][i])/denom;
             a = 1 - b;
             if ((b < 0) || (b > 1)) viable=0;  /* 0 < b < 1 */
	     for (k=0;k<num_active;k++) {       /* lagrange multipliers>0 */
                 if ((a*mG[k][i] + b*mG[k][j]) <= 0) viable=0;
             }
          } else {
             viable = 0;                        /* Linearly dependent */
          }

          /* find the minimum of viable directions */
          if (viable) {
             for (k=0;k<mDimension;k++) {
                g_bar[k] = a * dir[i][k] + b * dir[j][k];
             }
             MSQ_DOT(projection,g_bar,g_bar,mDimension);
             if (projection < min_projection) {
	        min_projection = projection;
                MSQ_COPY_VECTOR(temp_search,g_bar,mDimension);
                mSteepest = mActive->active_ind[i];
             }
          }
       }
    }
    if (optStatus != MSQ_EQUILIBRIUM) {
        MSQ_COPY_VECTOR(mSearch,temp_search,mDimension);
    }
}         


inline void NonSmoothSteepestDescent::solve2x2(double a11, double a12,
                                               double a21, double a22, 
                                               double b1, double b2,
                                               double **x, MsqError &/*err*/)
{
    double factor;

    /* if the system is not singular, solve it */
    if (fabs(a11*a22 - a21*a12) > MSQ_MACHINE_EPS) {
	(*x)=(double *)malloc(sizeof(double)*2);
	if (fabs(a11) > MSQ_MACHINE_EPS) {
	    factor = (a21/a11);
	    (*x)[1] = (b2 - factor*b1)/(a22 - factor*a12);
	    (*x)[0] = (b1 - a12*(*x)[1])/a11;
	} else if (fabs(a21) > MSQ_MACHINE_EPS) {
	    factor = (a11/a21);
	    (*x)[1] = (b1 - factor*b2)/(a12 - factor*a22);
	    (*x)[0] = (b2 - a22*(*x)[1])/a21;
	}
    } else {
	(*x) = NULL;
    }
}


inline void NonSmoothSteepestDescent::form_reduced_matrix(double ***P,
                                                          MsqError &/*err*/)
{
    int i,j;
    int num_active = mActive->num_active;

    (*P)=(double **)malloc(sizeof(double *)*(num_active-1));
    for (i=0; i<num_active-1; i++) 
        (*P)[i]=(double *)malloc(sizeof(double)*(num_active-1));

    for (i=0;i<num_active-1;i++) {
        (*P)[i][i] = mG[0][0] - 2*mG[0][i+1] + mG[i+1][i+1];
        for (j=i+1;j<num_active-1;j++) {
            (*P)[i][j] = mG[0][0] - mG[0][j+1] - mG[i+1][0] + mG[i+1][j+1];
            (*P)[j][i] = (*P)[i][j];
        }
    }
}


inline void  NonSmoothSteepestDescent::get_min_estimate(double *final_est,
                                                        MsqError &/*err*/)
{
    int    i;
    double est_imp;

    *final_est = -1E300;
    for (i=0;i<mActive->num_active;i++) {
	est_imp = -mAlpha*mGS[mActive->active_ind[i]];
        if (est_imp>*final_est) *final_est = est_imp;
    }
    if (*final_est == 0) {
	*final_est = -1E300;
	for (i=0;i<numFunctionValues;i++) {
	    est_imp = -mAlpha*mGS[i];
	    if ((est_imp>*final_est) && (fabs(est_imp) > MSQ_MACHINE_EPS)) {
		*final_est = est_imp;
	    }
	}
    }
}


inline void NonSmoothSteepestDescent::get_gradient_projections(MsqError &/*err*/)
{
    for (int i=0;i<numFunctionValues;i++) 
	MSQ_DOT(mGS[i],mGradient[i],mSearch,mDimension);

    MSQ_PRINT(3)("steepest in get_gradient_projections %d\n",mSteepest);
}


inline void NonSmoothSteepestDescent::compute_alpha(MsqError &/*err*/)
{
//    int       ierr;
//    int       j;
    int       i;
//    int       ind;
    int       num_values;
    double    steepest_function;
    double    steepest_grad;
    double    alpha_i;
    double    min_positive_value=1E300;

//    FUNCTION_TIMER_START("Compute Alpha");

    MSQ_PRINT(2)("In compute alpha\n");

    num_values = numFunctionValues;
    mAlpha = 1E300;

    steepest_function = mFunction[mSteepest];
    steepest_grad = mGS[mSteepest];
    for (i=0;i<num_values;i++)
    {
        /* if it's not active */
      if (i!=mSteepest)
      {
	  alpha_i = steepest_function - mFunction[i];
	   
	  if (fabs(mGS[mSteepest] - mGS[i])>1E-13) {
	     /* compute line intersection */
	     alpha_i = alpha_i/(steepest_grad - mGS[i]);
	  } else {
	     /* the lines don't intersect - it's not under consideration*/
	     alpha_i = 0;
	  }
	  if ((alpha_i > minStepSize ) && (fabs(alpha_i) < fabs(mAlpha))) {
	    mAlpha = fabs(alpha_i); 
            MSQ_PRINT(3)("Setting alpha %d %g\n",i,alpha_i);
	  }
          if ((alpha_i > 0) && (alpha_i < min_positive_value)) {
            min_positive_value = alpha_i;
          }
       }
    }

    if ((mAlpha == 1E300) && (min_positive_value != 1E300)) {
      mAlpha = min_positive_value;
    }

    /* if it never gets set, set it to the default */
    if (mAlpha == 1E300) {
      mAlpha = maxAlpha;
      MSQ_PRINT(3)("Setting alpha to the maximum step length\n");
    }

    MSQ_PRINT(3)("  The initial step size: %f\n",mAlpha);

//    FUNCTION_TIMER_END();
}


inline void NonSmoothSteepestDescent::copy_active(ActiveSet *active1, ActiveSet *active2, 
                                          MsqError &err)
{
    if (active1==NULL || active2==NULL) {
       MSQ_SETERR(err)("Null memory in copy_active\n",MsqError::INVALID_ARG);
       return;
    }

    active2->num_active = active1->num_active;
    active2->num_equal  = active1->num_equal;
    active2->true_active_value = active1->true_active_value;
    for (int i=0;i<active1->num_active;i++) {
	active2->active_ind[i] = active1->active_ind[i];
    }
}


inline void NonSmoothSteepestDescent::print_active_set(ActiveSet *active_set, 
                                                       double * func,
                                                       MsqError &err)
{
    if (active_set==0) {
      MSQ_SETERR(err)("Null ActiveSet", MsqError::INVALID_ARG);
      return;
    }
 
    if (active_set->num_active == 0) MSQ_DBGOUT(3)<< "No active values\n";
    /* print the active set */
    for (int i=0;i<active_set->num_active;i++) {
     MSQ_PRINT(3)("Active value %d:   %f \n",
	             i+1,func[active_set->active_ind[i]]); 
    }
}


inline void NonSmoothSteepestDescent::init_opt(MsqError &err)
{
    int        i, j;

    MSQ_PRINT(2)("\nInitializing Optimization \n");
    if (numFunctionValues > 150) {
      MSQ_SETERR(err)("num_values exceeds 150", MsqError::INVALID_STATE);
    }
    /* for the purposes of initialization will be set to zero after */
    equilibriumPt = 0;
    optStatus = 0;
    iterCount = 0;
    optIterCount = 0;
    mSteepest = 0;
    mAlpha = 0;
    maxAlpha = 0;

    MSQ_PRINT(3)("  Initialized Constants \n");
    for (i=0;i<3;i++) {
      mSearch[i] = 0;
      pdgInd[i] = -1;
      for (j=0;j<3;j++) mPDG[i][j] = 0;
    }

    MSQ_PRINT(3)("  Initialized search and PDG \n");
    for (i=0;i<numFunctionValues;i++) {
       mFunction[i] = 0;
       testFunction[i] = 0;
       originalFunction[i] = 0;
       mGS[i] = 0;
       for (j=0;j<3;j++) {
           mGradient[i][j] = 0;
       }
    }
    MSQ_PRINT(3)("  Initialized function/gradient \n");
    if (numFunctionValues > 150) {
      for (i=0;i<150;i++) {
       for (j=0;j<150;j++) mG[i][j] = -1;
      }
    } else {
      for (i=0;i<numFunctionValues;i++) {
       for (j=0;j<numFunctionValues;j++) mG[i][j] = -1;
      }
    }
    MSQ_PRINT(3)("  Initialized G\n");
 
    for (i=0;i<100;i++) prevActiveValues[i] = 0;
    MSQ_PRINT(3)("  Initialized prevActiveValues\n");
}


inline void NonSmoothSteepestDescent::init_max_step_length(MsqError &err)
{
  int i, j, k;
  double max_diff = 0;
  double diff=0;

  MSQ_PRINT(2)("In init_max_step_length\n");

  /* check that the input data is correct */
  if (numElements==0) {
    MSQ_SETERR(err)("Num incident vtx = 0\n",MsqError::INVALID_MESH);
    return;
  }
  if ((mDimension!=2) && (mDimension!=3)) {
     MSQ_SETERR(err)("Problem dimension is incorrect", MsqError::INVALID_MESH);
     return;
  }

  /* find the maximum distance between two incident vertex locations */
  for (i=1;i<numVertices;i++) {
    for (j=i;j<numVertices+1;j++) {
      diff=0;
      for (k=0;k<mDimension;k++) {
        diff += (mCoords[i][k]-mCoords[j][k])*(mCoords[i][k]-mCoords[j][k]);
      }
      if (max_diff < diff) max_diff=diff;
    } 
  }
  max_diff = sqrt(max_diff);
  if (max_diff==0) {
     MSQ_SETERR(err)("Maximum distance between incident vertices = 0\n",
                     MsqError::INVALID_MESH);
    return;
  }
  maxAlpha = max_diff/100;

  MSQ_PRINT(3)("  Maximum step is %g\n",maxAlpha);
}


} // namespace

#endif  // Mesquite_NonSmoothSteepestDescent_hpp 
