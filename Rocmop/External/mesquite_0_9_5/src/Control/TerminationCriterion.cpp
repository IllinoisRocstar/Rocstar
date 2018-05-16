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
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
// DESCRIPTION:
// ============
/*! \file TerminationCriterion.cpp
  
    \brief  Member functions of the Mesquite::TerminationCriterion class

    \author Michael Brewer
    \author Thomas Leurent
    \date   Feb. 14, 2003
 */

#include "TerminationCriterion.hpp"
#include "MeshSet.hpp"
#include "MsqVertex.hpp"
#include "MsqInterrupt.hpp"
#include "ObjectiveFunction.hpp"
#include "MsqError.hpp"
#include "MsqDebug.hpp"

namespace Mesquite {


const unsigned long GRAD_FLAGS = TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE |
                                 TerminationCriterion::GRADIENT_INF_NORM_ABSOLUTE |
                                 TerminationCriterion::GRADIENT_L2_NORM_RELATIVE |
                                 TerminationCriterion::GRADIENT_INF_NORM_RELATIVE;
const unsigned long OF_FLAGS   = TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE |
                                 TerminationCriterion::QUALITY_IMPROVEMENT_RELATIVE |
                                 TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_ABSOLUTE |
                                 TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_RELATIVE;

/*!Constructor initializes all of the data members which are not
  necessarily automatically initialized in their constructors.*/
TerminationCriterion::TerminationCriterion() 
  : mGrad(8),
    initialVerticesMemento(0),
    previousVerticesMemento(0),
    debugLevel(2)
{
  terminationCriterionFlag=NONE;
  cullingMethodFlag=NONE;
  cullingEps=0.0;
  initialOFValue=0.0;
  previousOFValue=0.0;
  currentOFValue = 0.0;
  lowerOFBound=0.0;
  initialGradL2Norm=0.0;
  initialGradInfNorm=0.0;
    //initial size of the gradient array
  gradL2NormAbsoluteEps=0.0;
  gradL2NormRelativeEps=0.0;
  gradInfNormAbsoluteEps=0.0;
  gradInfNormRelativeEps=0.0;
  qualityImprovementAbsoluteEps=0.0;
  qualityImprovementRelativeEps=0.0;
  iterationBound=0;
  iterationCounter=0;
  timeBound=0.0;
  vertexMovementAbsoluteEps=0.0;
  vertexMovementRelativeEps=0.0;
  successiveImprovementsAbsoluteEps=0.0;
  successiveImprovementsRelativeEps=0.0;
  boundedVertexMovementEps=0.0;
  
}



/*! Function to add a type of termination criterion to this object.  It is
  only valid if the specified criterion type requires a single double
  value.*/
void TerminationCriterion::add_criterion_type_with_double(TCType tc_type,
                                                       double eps,
                                                       MsqError &err)
{
  switch(tc_type){
    case GRADIENT_L2_NORM_ABSOLUTE:
       terminationCriterionFlag|=GRADIENT_L2_NORM_ABSOLUTE;
       gradL2NormAbsoluteEps=eps;
       break; 
    case GRADIENT_INF_NORM_ABSOLUTE:
       terminationCriterionFlag|=GRADIENT_INF_NORM_ABSOLUTE;
       gradInfNormAbsoluteEps=eps;
       break;
    case GRADIENT_L2_NORM_RELATIVE:
       terminationCriterionFlag|=GRADIENT_L2_NORM_RELATIVE;
       gradL2NormRelativeEps=eps;
       break;  
    case GRADIENT_INF_NORM_RELATIVE:
       terminationCriterionFlag|=GRADIENT_INF_NORM_RELATIVE;
       gradInfNormRelativeEps=eps;
       break;  
    case QUALITY_IMPROVEMENT_ABSOLUTE:
       terminationCriterionFlag|=QUALITY_IMPROVEMENT_ABSOLUTE;
       qualityImprovementAbsoluteEps=eps;
       break;
    case QUALITY_IMPROVEMENT_RELATIVE:
       terminationCriterionFlag|=QUALITY_IMPROVEMENT_RELATIVE;
       qualityImprovementRelativeEps=eps;
       break;
    case CPU_TIME:
       terminationCriterionFlag|=CPU_TIME;
       timeBound=eps;
       break;
    case VERTEX_MOVEMENT_ABSOLUTE:
       terminationCriterionFlag|=VERTEX_MOVEMENT_ABSOLUTE;
         //we actually compare squared movement to squared epsilon
       vertexMovementAbsoluteEps=(eps*eps);
       break;
    case VERTEX_MOVEMENT_RELATIVE:
       terminationCriterionFlag|=VERTEX_MOVEMENT_RELATIVE;
         //we actually compare squared movement to squared epsilon
       vertexMovementRelativeEps=(eps*eps);
       break;
    case SUCCESSIVE_IMPROVEMENTS_ABSOLUTE:
       terminationCriterionFlag|=SUCCESSIVE_IMPROVEMENTS_ABSOLUTE;
       successiveImprovementsAbsoluteEps=eps;
       break;
    case SUCCESSIVE_IMPROVEMENTS_RELATIVE:
       terminationCriterionFlag|=SUCCESSIVE_IMPROVEMENTS_RELATIVE;
       successiveImprovementsRelativeEps=eps;
       break;
    case BOUNDED_VERTEX_MOVEMENT:
       terminationCriterionFlag|=BOUNDED_VERTEX_MOVEMENT;
       boundedVertexMovementEps=eps;
       break;
    default:
       MSQ_SETERR(err)("TCType not valid for this function.",MsqError::INVALID_ARG);
  };
}

/*! Function to add a type of termination criterion to this object.  It is
  only valid if the specified criterion type requires a single integer
  value.*/
void TerminationCriterion::add_criterion_type_with_int(TCType tc_type,
                                                    int bound,
                                                    MsqError &err)
{
  switch(tc_type){
    case NUMBER_OF_ITERATES:
       terminationCriterionFlag|=NUMBER_OF_ITERATES;
       iterationBound=bound;
       break;
    default:
       MSQ_SETERR(err)("TCType not valid for this function.",MsqError::INVALID_ARG);
  };
}
/*! Function to remove a previously set criterion type.*/
void TerminationCriterion::remove_criterion_type(TCType tc_type,
                                                 MsqError &/*err*/)
{
  terminationCriterionFlag&=(~tc_type);
}


/*! Function to add a type of termination criterion to this object which
  will be used for culling purposes only.  As part of the global
  termination criterion, the outer criterion will also check to make
  sure that there are still free vertices in the mesh.*/
void TerminationCriterion::set_culling_type(TCType tc_type, double eps,
                                         MsqError &err)
{
  switch(tc_type){
    case QUALITY_IMPROVEMENT_ABSOLUTE:
       cullingMethodFlag=QUALITY_IMPROVEMENT_ABSOLUTE;
       break;
    case QUALITY_IMPROVEMENT_RELATIVE:
       cullingMethodFlag=QUALITY_IMPROVEMENT_RELATIVE;
       break;
    case VERTEX_MOVEMENT_ABSOLUTE:
       cullingMethodFlag=VERTEX_MOVEMENT_ABSOLUTE;
       break;
    case VERTEX_MOVEMENT_RELATIVE:
       cullingMethodFlag=VERTEX_MOVEMENT_RELATIVE;
       break;   
    case SUCCESSIVE_IMPROVEMENTS_ABSOLUTE:
       cullingMethodFlag=SUCCESSIVE_IMPROVEMENTS_ABSOLUTE;
       break;
     case SUCCESSIVE_IMPROVEMENTS_RELATIVE:
       cullingMethodFlag=SUCCESSIVE_IMPROVEMENTS_RELATIVE;
       break;  
    default:
       MSQ_SETERR(err)("TCType not valid for this function.",MsqError::INVALID_ARG);
  };
  cullingEps=eps;
}

/*!Sets the culling type to be NONE.*/
void TerminationCriterion::remove_culling(MsqError &/*err*/)
{
  cullingMethodFlag=NONE;
}

  

/*!This version of reset is called using a MeshSet, which implies
  it is only called when this criterion is used as the 'outer' termination
  criterion.  
 */
void TerminationCriterion::reset_outer(MeshSet &ms, ObjectiveFunction* obj_ptr,
                                       MsqError &err)
{
  const unsigned long totalFlag = terminationCriterionFlag | cullingMethodFlag;
  PatchData global_patch;
  
    //if we need to fill out the global patch data object.
  if (totalFlag & (GRAD_FLAGS | OF_FLAGS | VERTEX_MOVEMENT_RELATIVE))
  {
    PatchDataParameters global_patch_params;
    global_patch_params.set_patch_type( PatchData::GLOBAL_PATCH, err, 0, 0 );
       MSQ_ERRRTN(err);
    ms.get_next_patch( global_patch, global_patch_params, err ); 
       MSQ_ERRRTN(err);
  }

    //now call the other reset
  reset_inner( global_patch, obj_ptr, err ); MSQ_ERRRTN(err);
}
    
/*!Reset function using using a PatchData object.  This function is
  called for the inner-stopping criterion directly from the
  loop over mesh function in VertexMover.  For outer criterion,
  it is called from the reset function which takes a MeshSet object.
  This function prepares the object to be used by setting the initial
  values of some of the data members.  As examples, if needed, it resets
  the cpu timer to zero, the iteration counter to zero, and the
  initial and previous objective function values to the current
  objective function value for this patch.
  The return value for this function is similar to that of terminate().
  The function returns false if the checked criteria have not been
  satisfied, and true if they have been.  reset() only checks the
  GRADIENT_INF_NORM_ABSOLUTE, GRADIENT_L2_NORM_ABSOLUTE, and the
  QUALITY_IMPROVEMENT_ABSOLUTE criteria.  Checking these criteria
  allows the QualityImprover to skip the entire optimization if
  the initial mesh satisfies the appropriate conditions.
 */
void TerminationCriterion::reset_inner(PatchData &pd, ObjectiveFunction* obj_ptr,
                                    MsqError &err)
{
  const unsigned long totalFlag = terminationCriterionFlag | cullingMethodFlag;
  OFPtr = obj_ptr;
  
    // clear flag for BOUNDED_VERTEX_MOVEMENT
  vertexMovementExceedsBound = 0;
  
    // Use -1 to denote that this isn't initialized yet.
    // As all valid values must be >= 0.0, a negative
    // value indicates that it is uninitialized and is
    // always less than any valid value.
  maxSquaredMovement = -1;
  
    // Clear the iteration count.
  iterationCounter = 0;
  
    //reset the inner timer if needed
  if(totalFlag & CPU_TIME){
    mTimer.reset();
  }
   
    //GRADIENT
  if(totalFlag & GRAD_FLAGS)
  {
    int num_vertices=pd.num_vertices();
    mGrad.resize( num_vertices );

      //get gradient and make sure it is valid
    bool b = obj_ptr->compute_gradient(pd, &mGrad[0] , currentOFValue, 
                                       err, num_vertices); MSQ_ERRRTN(err);
    if (!b) {
      MSQ_SETERR(err)("Initial patch is invalid for gradient computation.", 
                      MsqError::INVALID_STATE);
      return;
    } 

      //get the gradient norms
    if (totalFlag & (GRADIENT_INF_NORM_ABSOLUTE|GRADIENT_INF_NORM_RELATIVE))
    {
      currentGradInfNorm = initialGradInfNorm = Linf(&mGrad[0], num_vertices);
      MSQ_DBGOUT(debugLevel) << "  o Initial gradient Inf norm: " 
        << initialGradInfNorm << msq_stdio::endl;
    }  
    if (totalFlag & (GRADIENT_L2_NORM_ABSOLUTE|GRADIENT_L2_NORM_RELATIVE))
    {
      currentGradL2Norm = initialGradL2Norm = length(&mGrad[0], num_vertices);
      MSQ_DBGOUT(debugLevel) << "  o Initial gradient L2 norm: " 
        << initialGradL2Norm << msq_stdio::endl;
    }  
      //the OFvalue comes for free, so save it
    previousOFValue=currentOFValue;
    initialOFValue=currentOFValue;
  }
  //find the initial objective function value if needed and not already
  //computed.  If we needed the gradient, we have the OF value for free.
  else if (totalFlag & OF_FLAGS)
  {
      //ensure the obj_ptr is not null
    if(obj_ptr==NULL){
      MSQ_SETERR(err)("Error termination criteria set which uses objective "
                      "functions, but no objective function is available.",
                      MsqError::INVALID_STATE);
      return;
    }
    
    bool b = obj_ptr->evaluate(pd, currentOFValue, err); MSQ_ERRRTN(err);
    if (!b){
      MSQ_SETERR(err)("Initial patch is invalid for evaluation.",MsqError::INVALID_STATE);
      return;
    }
      //std::cout<<"\nReseting initial of value = "<<initialOFValue;
    previousOFValue=currentOFValue;
    initialOFValue=currentOFValue;
  }
  
  if (totalFlag & (GRAD_FLAGS|OF_FLAGS))
    MSQ_DBGOUT(debugLevel) << "  o Initial OF value: " << initialOFValue << msq_stdio::endl;
  
    // Store current vertex locations now, because we'll
    // need them later to compare the current movement with.
  if (totalFlag & VERTEX_MOVEMENT_RELATIVE)
  {
    if (initialVerticesMemento)
    {
      pd.recreate_vertices_memento( initialVerticesMemento, err );
    }
    else
    {
      initialVerticesMemento = pd.create_vertices_memento( err );
    }
    MSQ_ERRRTN(err);
    maxSquaredInitialMovement = DBL_MAX;
  }
}

void TerminationCriterion::reset_patch(PatchData &pd, MsqError &err)
{
  const unsigned long totalFlag = terminationCriterionFlag | cullingMethodFlag;
  if (totalFlag & (VERTEX_MOVEMENT_ABSOLUTE | VERTEX_MOVEMENT_RELATIVE))
  {
    if (previousVerticesMemento)
      pd.recreate_vertices_memento(previousVerticesMemento,err); 
    else
      previousVerticesMemento = pd.create_vertices_memento(err);
    MSQ_ERRRTN(err);
  }
}

void TerminationCriterion::accumulate_inner( PatchData& pd, MsqError& err )
{
  double of_value = 0;
  
  if (terminationCriterionFlag & GRAD_FLAGS)
  {
    mGrad.resize( pd.num_vertices() );
    bool b = OFPtr->compute_gradient(pd, &mGrad[0], of_value, err, pd.num_vertices());
    MSQ_ERRRTN(err);
    if (!b) {
      MSQ_SETERR(err)("Initial patch is invalid for gradient compuation.",
                      MsqError::INVALID_MESH);
      return;
    }
  }
  else if (terminationCriterionFlag & OF_FLAGS)
  {
    bool b = OFPtr->evaluate(pd, of_value, err); MSQ_ERRRTN(err);
    if (!b) {
      MSQ_SETERR(err)("Invalid patch passed to TerminationCriterion.",
                      MsqError::INVALID_MESH);
      return;
    }
  }

  accumulate_inner( pd, of_value, &mGrad[0], err );  MSQ_CHKERR(err);
}


void TerminationCriterion::accumulate_inner( PatchData& pd, 
                                             double of_value,
                                             Vector3D* grad_array,
                                             MsqError& err )
{
  //if terminating on the norm of the gradient
  currentGradL2Norm = 10e6;
  if (terminationCriterionFlag & (GRADIENT_L2_NORM_ABSOLUTE | GRADIENT_L2_NORM_RELATIVE)) 
  {
    currentGradL2Norm = length(grad_array, pd.num_vertices()); // get the L2 norm
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- gradient L2 norm: " 
      << currentGradL2Norm << msq_stdio::endl;
  }
  currentGradInfNorm = 10e6;
  if (terminationCriterionFlag & (GRADIENT_INF_NORM_ABSOLUTE | GRADIENT_INF_NORM_RELATIVE)) 
  {
    currentGradInfNorm = length(grad_array, pd.num_vertices()); // get the Linf norm
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- gradient Inf norm: " 
      << currentGradInfNorm << msq_stdio::endl;
  } 
  
  if (terminationCriterionFlag & VERTEX_MOVEMENT_RELATIVE)
  {
    maxSquaredInitialMovement = pd.get_max_vertex_movement_squared(
                               initialVerticesMemento, err );  MSQ_ERRRTN(err);
  }
  
  previousOFValue = currentOFValue;
  currentOFValue = of_value;
  if (terminationCriterionFlag & OF_FLAGS)
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- OF Value: " << of_value << msq_stdio::endl;
  else if (grad_array)
    MSQ_DBGOUT(debugLevel) << "  o OF Value: " << of_value << msq_stdio::endl;

  ++iterationCounter;
}


void TerminationCriterion::accumulate_outer(MeshSet &ms, MsqError &err)
{
  PatchData global_patch;
  
    //if we need to fill out the global patch data object.
  if (terminationCriterionFlag & (GRAD_FLAGS|OF_FLAGS|VERTEX_MOVEMENT_RELATIVE))
  {
    PatchDataParameters global_params;
    global_params.set_patch_type( PatchData::GLOBAL_PATCH, err, 0, 0 );MSQ_ERRRTN(err);
    ms.get_next_patch( global_patch, global_params, err );             MSQ_ERRRTN(err);
  }
  
  accumulate_inner( global_patch, err );                             MSQ_ERRRTN(err);
}


void TerminationCriterion::accumulate_patch( PatchData& pd, MsqError& err )
{
  if (terminationCriterionFlag & (VERTEX_MOVEMENT_ABSOLUTE|VERTEX_MOVEMENT_RELATIVE))
  {
    double patch_max_dist = pd.get_max_vertex_movement_squared( previousVerticesMemento, err );
    if (patch_max_dist > maxSquaredMovement)
      maxSquaredMovement = patch_max_dist;
    pd.recreate_vertices_memento( previousVerticesMemento, err );  MSQ_ERRRTN(err);
  }
    
    //if terminating on bounded vertex movement (a bounding box for the mesh)
  if(terminationCriterionFlag & BOUNDED_VERTEX_MOVEMENT)
  {
    MsqVertex* vert = pd.get_vertex_array(err);
    int num_vert = pd.num_vertices();
    int i=0;
      //for each vertex
    for(i=0;i<num_vert;++i)
    {
        //if any of the coordinates are greater than eps
      if( (vert[i][0]>boundedVertexMovementEps) ||
          (vert[i][1]>boundedVertexMovementEps) ||
          (vert[i][2]>boundedVertexMovementEps) )
      {
        ++vertexMovementExceedsBound;
      }
    }
  }
}


/*!  This function evaluates the needed information and then evaluates
  the termination criteria.  If any of the selected criteria are satisfied,
  the function returns true.  Otherwise, the function returns false.
 */
bool TerminationCriterion::terminate( )
{
  bool return_flag = false;
  //  cout<<"\nInside terminate(pd,of,err):  flag = "<<terminationCriterionFlag << endl;

    //First check for an interrupt signal
  if (MsqInterrupt::interrupt())
  {
     MSQ_DBGOUT(debugLevel) << "  o TermCrit -- INTERRUPTED" << msq_stdio::endl;
    return true;
  }
  
    //if terminating on numbering of inner iterations
  if (NUMBER_OF_ITERATES & terminationCriterionFlag
    && iterationCounter >= iterationBound)
  {
    return_flag = true;
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- Reached " << iterationBound << " iterations." << msq_stdio::endl;
  }
  
  if (CPU_TIME & terminationCriterionFlag && mTimer.since_birth()>=timeBound)
  {
    return_flag=true;
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- Exceeded CPU time." << msq_stdio::endl;
  }
  
  
  if ((VERTEX_MOVEMENT_ABSOLUTE|VERTEX_MOVEMENT_RELATIVE) & terminationCriterionFlag
      && maxSquaredMovement >= 0.0)
  {
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- Maximuim vertex movement: "
    << sqrt(maxSquaredMovement) << msq_stdio::endl;

    if (VERTEX_MOVEMENT_ABSOLUTE & terminationCriterionFlag 
        && maxSquaredMovement <= vertexMovementAbsoluteEps)
    {
      return_flag = true;
    }
  
    if (VERTEX_MOVEMENT_RELATIVE & terminationCriterionFlag
        && maxSquaredMovement <= vertexMovementRelativeEps*maxSquaredInitialMovement)
    {
      return_flag = true;
    }
    
      // Clear this value at the end of each iteration.
    maxSquaredMovement = -1.0;
  }

  if (GRADIENT_L2_NORM_ABSOLUTE & terminationCriterionFlag &&
      currentGradL2Norm <= gradL2NormAbsoluteEps)
  {
    return_flag = true;
  }
  
  if (GRADIENT_INF_NORM_ABSOLUTE & terminationCriterionFlag &&
      currentGradInfNorm <= gradInfNormAbsoluteEps)
  {
    return_flag = true;
  }
  
  if (GRADIENT_L2_NORM_RELATIVE & terminationCriterionFlag &&
      currentGradL2Norm <= (gradL2NormRelativeEps * initialGradL2Norm))
  {
    return_flag = true;
  }
  
  if (GRADIENT_INF_NORM_RELATIVE & terminationCriterionFlag &&
      currentGradInfNorm <= (gradInfNormRelativeEps * initialGradInfNorm))
  {
    return_flag = true;
  }
  
  if (QUALITY_IMPROVEMENT_ABSOLUTE & terminationCriterionFlag &&
      currentOFValue <= qualityImprovementAbsoluteEps)
  {
    return_flag = true;
  }
  
  if (QUALITY_IMPROVEMENT_RELATIVE & terminationCriterionFlag &&
      (currentOFValue - lowerOFBound) <= 
        qualityImprovementRelativeEps * (initialOFValue - lowerOFBound))
  {
    return_flag = true;
  }
  
  if (SUCCESSIVE_IMPROVEMENTS_ABSOLUTE & terminationCriterionFlag &&
     (previousOFValue - currentOFValue) <= successiveImprovementsRelativeEps)
  {
    return_flag = true;
  }
  
  if (SUCCESSIVE_IMPROVEMENTS_RELATIVE & terminationCriterionFlag &&
      (previousOFValue - currentOFValue) <= 
        successiveImprovementsRelativeEps * (initialOFValue - currentOFValue))
  {
    return_flag = true;
  }
  
  if (BOUNDED_VERTEX_MOVEMENT & terminationCriterionFlag && vertexMovementExceedsBound)
  {
    return_flag = true;
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- " << vertexMovementExceedsBound
                           << " vertices out of bounds." << msq_stdio::endl;
  }
  
    // clear this value at the end of each iteration
  vertexMovementExceedsBound = 0;

    //if none of the criteria were satisfied
  return return_flag;
}


  

/*!This function checks the culling method criterion supplied to the object
  by the user.  If the user does not supply a culling method criterion,
  the default criterion is NONE, and in that case, no culling is performed.
  If the culling method criterion is satisfied, the interior vertices
  of the given patch are flagged as soft_fixed.  Otherwise, the soft_fixed
  flag is removed from each of the vertices in the patch (interior and
  boundary vertices).  Also, if the criterion was satisfied, then the
  function returns true.  Otherwise, the function returns false.
 */
bool TerminationCriterion::cull_vertices(PatchData &pd,
                                      ObjectiveFunction* obj_ptr,
                                      MsqError &err)
{
    //PRINT_INFO("CULLING_METHOD FLAG = %i",cullingMethodFlag);
  
    //cull_bool will be changed to true if the criterion is satisfied
  bool b, cull_bool=false;
  double prev_m, init_m;
  switch(cullingMethodFlag){
      //if no culling is requested, always return false
    case NONE:
       return cull_bool;
         //if culling on quality improvement absolute
    case QUALITY_IMPROVEMENT_ABSOLUTE:
         //get objective function value
       b = obj_ptr->evaluate(pd, currentOFValue, err);
       if (MSQ_CHKERR(err)) return false;
       if (!b) {
         MSQ_SETERR(err)(MsqError::INVALID_MESH);
         return false;
       }
         //if the improvement was enough, cull
       if(currentOFValue <= cullingEps)
       {
         cull_bool=true;  
       }
         //PRINT_INFO("\ncurrentOFValue = %f, bool = %i\n",currentOFValue,cull_bool);
       
       break;
         //if culing on quality improvement relative
    case QUALITY_IMPROVEMENT_RELATIVE:
         //get objective function value
       b = obj_ptr->evaluate(pd, currentOFValue, err);
       if (MSQ_CHKERR(err)) return false;
       if(!b){
         MSQ_SETERR(err)(MsqError::INVALID_MESH);
         return false;
       }
         //if the improvement was enough, cull
       if((currentOFValue-lowerOFBound)<=
          (cullingEps*(initialOFValue-lowerOFBound)))
       {
         cull_bool=true;  
       }
       break;
         //if culling on vertex movement absolute
    case VERTEX_MOVEMENT_ABSOLUTE:
         //if movement was enough, cull
       prev_m = pd.get_max_vertex_movement_squared(previousVerticesMemento,err);
       MSQ_ERRZERO(err);
       if(prev_m <= cullingEps){
         cull_bool=true;  
       }
       
       break;
         //if culling on vertex movement relative
    case VERTEX_MOVEMENT_RELATIVE:
         //if movement was small enough, cull
       prev_m = pd.get_max_vertex_movement_squared(previousVerticesMemento,err);
       MSQ_ERRZERO(err);
       init_m = pd.get_max_vertex_movement_squared(initialVerticesMemento,err);
       MSQ_ERRZERO(err);
       if(prev_m <= (cullingEps * init_m)){
         cull_bool=true;  
       }
       break;
    default:
       MSQ_SETERR(err)("Requested culling method not yet implemented.",
                       MsqError::NOT_IMPLEMENTED);
       return false;
  };
    //Now actually have patch data cull vertices
  if(cull_bool)
  {
    pd.set_free_vertices_soft_fixed(err); MSQ_ERRZERO(err);
  }
  else
  {
    pd.set_all_vertices_soft_free(err); MSQ_ERRZERO(err);
  }
  return cull_bool;
}

/*!
  Currently this only deletes the memento of the vertex positions and the
  mGrad vector if neccessary.
  When culling, we remove the soft fixed flags from all of the vertices.
 */
void TerminationCriterion::cleanup(MeshSet &ms, MsqError &err)
{
  delete previousVerticesMemento;
  delete initialVerticesMemento;
  previousVerticesMemento = 0;
  initialVerticesMemento = 0;

  if(cullingMethodFlag){
    ms.clear_all_soft_fixed_flags(err); MSQ_ERRRTN(err);
  }
}

} //namespace Mesquite
