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
  \file   ObjectiveFunction.cpp
  \brief  

  \author Michael Brewer
  \author Thomas Leurent
  
  \date   2002-08-02
*/

#include "ObjectiveFunction.hpp"
#include "MsqVertex.hpp"
#include "MsqDebug.hpp"
#include "MsqFreeVertexIndexIterator.hpp"

namespace Mesquite {

/*! 
  Numerically Calculates the gradient of the ObjectiveFunction for the
  free vertices in the patch.  Returns 'false' if the patch is outside
  of a required feasible region, returns 'ture' otherwise.
  The behavior of the function depends on the value of the boolean
  useLocalGradient.  If useLocalGradient is set to
  'true', compute_numerical_gradient creates a sub-patch around a free
  vertex, and then perturbs that vertex in one of the coordinate directions.
  Only the ObjectiveFunction value on the local sub-patch is used in the
  computation of the gradient.  Therefore, useLocalGradient should only
  be set to 'true' for ObjectiveFunctions which can use this method.  Unless
  the concrete ObjectiveFunction sets useLocalGradient to 'true' in its
  constructor, the value will be 'false'.  In this case, the objective
  function value for the entire patch is used in the calculation of the
  gradient.  This is computationally expensive, but it is numerically
  correct for all (C_1) functions.
  \param pd  PatchData on which the gradient is taken.
  \param grad  Array of Vector3D of length the number of vertices used to store gradient.
  \param OF_val will be set to the objective function value.
  \param array_size Either the length of grad or 0.
 */
bool ObjectiveFunction::compute_numerical_gradient(Mesquite::PatchData &pd,
                                                   Vector3D *const &grad,
                                                   double &OF_val,
                                                   MsqError &err,
                                                   size_t array_size)
{
  size_t num_vtx=pd.num_vertices();
  if(num_vtx!=array_size && array_size>0)
    MSQ_DBGOUT(1) << "\nArray size not equal to the number of vertices.\n";

  OF_val = 0.; // in case of return false. 
  MsqVertex* vertices=pd.get_vertex_array(err);
  double flocal=0;
  double flocald=0;
  double eps=0;
  size_t m=0;
  short j;
  
  if(useLocalGradient){
    //********************useLocalGradient***************************
    //if useLocalGradient is turned on, do more efficient computation
    PatchData sub_patch;
    for (m=0; m<num_vtx; ++m) {
      if (vertices[m].is_free_vertex()) {
        pd.get_subpatch(m, sub_patch, err); MSQ_ERRZERO(err);
        //If sub_patch is not in the feasible region, do not
        //calculate anything.  Just return false.
        bool b = evaluate(sub_patch,flocal,err);
        if(MSQ_CHKERR(err) || !b) {
          return false;
        }

        //loop over the three coords x,y,z
        for(j=0;j<3;++j){
          eps=get_eps(sub_patch, flocald, j, (&vertices[m]), err); MSQ_ERRZERO(err);
          //PRINT_INFO("\nin obj num grad j=%i, eps=%20.19f",j,eps);
          if(eps==0){
            MSQ_SETERR(err)("Dividing by zero in Objective Functions numerical grad",
                            MsqError::INVALID_STATE);
            return false;
          }
          grad[m][j]=(flocald-flocal)/eps;
        }
      }
      else {
        for(j=0;j<3;++j)
          grad[m][j] = 0.0;
      }
    }
    evaluate(pd, OF_val, err);  MSQ_ERRZERO(err);
  }
  else {
    //********************DO NOT useLocalGradient********************
    //if useLocalGradient is turned off, we do inefficient computation
    for (m=0; m<num_vtx; ++m) {

      if (vertices[m].is_free_vertex()) {
        //If pd is not in the feasible region, do not calculate anything.
        //Just return false.
        bool b = evaluate(pd,flocal,err);
        if(MSQ_CHKERR(err) || !b) {
          return false;
        }
        OF_val = flocal;
        //loop over the three coords x,y,z
        for(j=0;j<3;++j){
          eps=get_eps(pd, flocald, j, (&vertices[m]), err); MSQ_ERRZERO(err);
          //PRINT_INFO("\nin obj num grad j=%i, eps=%20.19f",j,eps);
          if(eps==0){
            MSQ_SETERR(err)("Dividing by zero in Objective Functions numerical grad",
                            MsqError::INVALID_STATE);
            return false;
          }
          grad[m][j]=(flocald-flocal)/eps;
        }
      }
      else {
        for(j=0;j<3;++j)
          grad[m][j] = 0.0;
      }
      //PRINT_INFO("  gradx = %f, grady = %f, gradz = %f\n",grad[m][0],grad[m][1],grad[m][2]);   
    }//end loop over all vertices
  }
  //*****************END of DO NOT useLocalGradient*****************
                        return true;
}

bool ObjectiveFunction::compute_analytical_gradient(PatchData &patch,
                                             Vector3D *const &grad,
                                             double &OF_val,
                                             MsqError &err, size_t array_size){
      set_gradient_type(NUMERICAL_GRADIENT);
      bool result = compute_numerical_gradient(patch, grad, OF_val, err, array_size);
      return !MSQ_CHKERR(err) && result;
    }

bool ObjectiveFunction::compute_analytical_hessian(PatchData &/*patch*/,
                                            MsqHessian &/*hessian*/,
                                            Vector3D *const &/*grad*/,
                                            double &/*OF_val*/,
                                            MsqError &err) {
      MSQ_SETERR(err)("Analytic hessian not implemented for this Objective "
                    "Function. Feasible Newton algorythm cannot be used.\n",
                    MsqError::INVALID_STATE);
      return false;
    }



} // namespace Mesquite

