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

/*! \file TerminationCriterion.hpp

Header file for the TerminationCriterion classes.

  \author Michael Brewer
  \author Thomas Leurent
  \date   Feb. 14, 2003
 */


#ifndef TerminationCriterion_hpp
#define TerminationCriterion_hpp

#include "Mesquite.hpp"
#include "PatchDataUser.hpp"
#include "MsqTimer.hpp"

#include <string>

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
#endif

namespace Mesquite
{
   class MeshSet;
   class MsqError;
   class ObjectiveFunction;
   
  /*! \class TerminationCriterion

      \brief The TerminationCriterion class contains functionality to
      terminate the VertexMover's optimization.

      The TerminationCriterion class has three roles.  It
      is used to terminate the optimization on a single patch; it
      is used to terminate the iterations over all patches in the
      mesh; and it is used to cull vertices frm the optimization
      processes.  Thus, for each optimzation, two TerminationCriterion
      objects are used.  The class contains five important member
      functions used in the VertexMover:  initialize(), reset(),
      terminate(), cull_vertices(), and cleanup().  These functions
      are each explained in detail below.  In general, the only one
      of these functions called directly from a concrete VertexMover
      is terminate() which allows the concrete VertexMover to determine
      when to stop producing new iterates on a given patch.  All other
      functionality is handled from the base VertexMover base class.

      There are several different types of termination criteria
      available.  These types are listed in teh enumberation
      TCType.  Multiple criteria types can be set on a given
      TermiantionCriterion object, and when this occurs, the
      optimization process will terminate whenever any of the
      criteria have been satisfied.
      
      The following is a brief description of how TerminationCriterion
      is used within Mesquite.  Functions called during QualityImprovement
      can be devided into three groups:
        reset_*      - Initialize data for an iteration
        accumulate_* - Update TC for changed data during iteration
        terminate    - Check if the termination criterion has been met.
      There are three different forms of the reset_* and accumulate_*
      functions which are called on the inner, outer, or both 
      TerminationCriterion classes:
        *_outer      - Called on outer termination criterion.
        *_inner      - Called on inner termination criterion.
        *_patch      - Called on outer termination criterion for
                       each patch and on inner termination criterion
                       for each inner iteration.
      
      If implementing a new TerminationCriterion, the following rules
      should be followed.  If the value must be calculated on a global
      patch for the outer TC, then:
        o The functionality should be added to *_inner (yes, INNER) 
        o The *_outer methods should be updated to call the *_inner 
            with a global patch when your TC is requested.
        o The internal data for any such TC should be initialized 
          in the reset_inner method.  
      If the value for the outer criterion can be calculated from each 
      local patch when iterating over the mesh with local patches, then:
        o The functionality should be added to *_patch
        o Any state values pertaining to the entire iteration must be 
           initialized in reset_inner(..) and cleared in terminate()
        o Any patch-specific data should be initialized in reset_patch
        o Care should be taken that terminate() does not check 
          uninitialized data if called before the first call to
          accumulate_patch()
  */
  class TerminationCriterion
  {
  public:
    
     /*! \enum TCType  defines the termination criterion */
    enum TCType {
       NONE    = 0,
       //! checks the gradient \f$\nabla f \f$ of objective function 
       //! \f$f : I\!\!R^{3N} \rightarrow I\!\!R \f$ against a double \f$d\f$  
       //! and stops when \f$\sqrt{\sum_{i=1}^{3N}\nabla f_i^2}<d\f$  
       GRADIENT_L2_NORM_ABSOLUTE = 1,  
       //! checks the gradient \f$\nabla f \f$ of objective function 
       //! \f$f : I\!\!R^{3N} \rightarrow I\!\!R \f$ against a double \f$d\f$  
       //! and stops when \f$ \max_{i=1}^{3N} \nabla f_i < d \f$  
       GRADIENT_INF_NORM_ABSOLUTE = 2,
         //!terminates on the j_th iteration when
         //! \f$\sqrt{\sum_{i=1}^{3N}\nabla f_{i,j}^2}<d\sqrt{\sum_{i=1}^{3N}\nabla f_{i,0}^2}\f$
         //!  That is, terminates when the norm of the gradient is small
         //! than some scaling factor times the norm of the original gradient. 
       GRADIENT_L2_NORM_RELATIVE = 4,
       //!terminates on the j_th iteration when
         //! \f$\max_{i=1 \cdots 3N}\nabla f_{i,j}<d \max_{i=1 \cdots 3N}\nabla f_{i,0}\f$
         //!  That is, terminates when the norm of the gradient is small
         //! than some scaling factor times the norm of the original gradient.
         //! (Using the infinity norm.)
       GRADIENT_INF_NORM_RELATIVE = 8,
         //! Not yet implemented.
       KKT  = 16,
         //!Terminates when the objective function value is smaller than
         //! the given scalar value.
       QUALITY_IMPROVEMENT_ABSOLUTE = 32,
         //!Terminates when the objective function value is smaller than
         //! the given scalar value times the original objective function
         //! value.
       QUALITY_IMPROVEMENT_RELATIVE = 64,
         //!Terminates when the number of iterations exceeds a given integer.
       NUMBER_OF_ITERATES = 128,
         //!Terminates when the algorithm exceeds an allotted time limit
         //! (given in seconds).
       CPU_TIME  = 256,
         //!Terminates when a the maximum distance moved by any vertex
         //! during the previous iteration is below the given value.
       VERTEX_MOVEMENT_ABSOLUTE  = 512,
         //!Terminates when a the maximum distance moved by any vertex
         //! during the previous iteration is below the given value
         //! times the maximum distance moved by any vertex over the
         //! entire course of the optimization.
       VERTEX_MOVEMENT_RELATIVE  = 1024,
         //!Terminates when the decrease in the objective function value since
         //! the previous iteration is below the given value.
       SUCCESSIVE_IMPROVEMENTS_ABSOLUTE = 2048,
         //!Terminates when the decrease in the objective function value since
         //! the previous iteration is below the given value times the
         //! decrease in the objective function value since the beginning
         //! of this optimization process.
       SUCCESSIVE_IMPROVEMENTS_RELATIVE = 4096,
         //!Terminates when any vertex leaves the bounding box, defined
         //! by the given value, d.  That is, when the absolute value of
         //! a single coordinate of vertex's position exceeds d.
       BOUNDED_VERTEX_MOVEMENT = 8192
    };

      //!Constructor which does not take any arguements
    TerminationCriterion();
    
      //!Destructor
    ~TerminationCriterion(){};

      //Functions with which the user can specify the criteria to be used
      //!Sets the criterion by specifing the TCType and the eps value
    void add_criterion_type_with_double(TCType tc_type, double eps,
                                        MsqError &err);
      //!Sets the criterion by specifing the TCType and the integer value
    void add_criterion_type_with_int(TCType tc_type, int bound,
                                     MsqError &err);
      //!Removes the criterion by specifing just the TCType.
    void remove_criterion_type(TCType tc_type, MsqError &err);
    
      //!Sets the type of criterion that the user would like to
      //! use for culling purposes (along with the associated tolerance.
    void set_culling_type(TCType tc_type, double eps, MsqError &err);
      //!Removes any previously set culling types (sets the culling
      //! type to be NONE).
    void remove_culling(MsqError &err);
    
      //! Clear any data accumulated during an outer iteration
    void reset_outer( MeshSet& ms, ObjectiveFunction* of, MsqError& err );
    
      //! Clear any data accumulated during an inner iteration
    void reset_inner( PatchData& pd, ObjectiveFunction* of, MsqError& err );
    
      //! Shared inner and outer initialization during inner loop
    void reset_patch( PatchData& pd, MsqError& err );
    
      //! Accumulate data during inner iteration
    void accumulate_inner( PatchData& pd, MsqError& err );
    
      //! Accumulate data during inner iteration
    void accumulate_inner( PatchData& pd, double of_value, Vector3D* of_grads, 
                           MsqError& err );
    
      //! Common code for both inner and outer termination 
      //! criteria during inner iteration.                       
    void accumulate_patch( PatchData& pd, MsqError& err );
    
    void accumulate_outer( MeshSet& ms, MsqError& err );
    
      //! Check if termination criterion has been met
    bool terminate();
    
    
      //!Function which determines whether this patch should be 'culled'
    bool cull_vertices(PatchData &pd, ObjectiveFunction* obj_ptr, MsqError &err);
      //!Cleans up after the TerminationCriterion is finished.
    void cleanup(MeshSet &ms, MsqError &err);

      //!This function returns the current function value.
      /*! \todo Michael:  this function is not reliable.  It
        needs to be more robust.  How do we know whether
        currentOFValue got updated or not?  We may want to
        make sure that all the criteria get checked.*/
    double get_current_function_value()
       {return currentOFValue;}
       
    void set_debug_output_level( int i )
      { debugLevel = i; }
    
 protected:
    
 private:
    //PRIVATE DATA MEMBERS
    long unsigned int terminationCriterionFlag;//!<Bit flag of termination crit
    long unsigned int cullingMethodFlag;/*!<Bit flag of criterion for culling*/
      //epsiloon used in culling methods.
    double cullingEps;

      // ObjectiveFunction pointer
    ObjectiveFunction* OFPtr;

      //Data not specific to a single criterion
    double initialOFValue;
    double previousOFValue;
    double currentOFValue;
    double lowerOFBound;

      //Data specific to termination criterion 1 (gradient bounds)
    msq_std::vector<Vector3D> mGrad;
    double initialGradL2Norm;
    double currentGradL2Norm;
    double gradL2NormAbsoluteEps;
    double gradL2NormRelativeEps;
    double initialGradInfNorm;
    double currentGradInfNorm;
    double gradInfNormAbsoluteEps;
    double gradInfNormRelativeEps;
      //Data specific to termination criterion 2 (KKT)
      //???????????????????????????????????????????
      //Data specific to termination criterion 3 (Quality Improvement)
    double qualityImprovementAbsoluteEps;
    double qualityImprovementRelativeEps;
      //Data specific to termination criterion 4 (inner iterations)
    int iterationBound;
    int iterationCounter;
      //Data specific to termination criterion 5 (cpu time)
    Timer mTimer;
    double timeBound;
      //Data specific to termination criterion 6 (vertex movement)
    PatchDataVerticesMemento* initialVerticesMemento;
    PatchDataVerticesMemento* previousVerticesMemento;//if we want relative
    double vertexMovementAbsoluteEps;
    double vertexMovementRelativeEps;
    double maxSquaredInitialMovement;
    double maxSquaredMovement;
    
      //Data specific to termination criterion 7 (successive improvement to F)
    double successiveImprovementsAbsoluteEps;
    double successiveImprovementsRelativeEps;
      //crit 8
    double boundedVertexMovementEps;
    int vertexMovementExceedsBound;
    
    int debugLevel;
    
  };

} //namespace


#endif // TerminationCriterion_hpp
