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
  \file   QualityImprover.hpp
  \brief  

  The Quality Improver Class is the base class for all the algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_QualityImprover_hpp 
#define Mesquite_QualityImprover_hpp

#include <string>

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "TerminationCriterion.hpp"
#include "PatchDataUser.hpp"

namespace Mesquite
{

  class MeshSet;
  
  /*! \class QualityImprover
    \brief Base class for all quality improvers.
    Mote that the PatchData settings are inherited from the PathDataUser class. 

  */ 
  class QualityImprover : public PatchDataUser
  {
  public:

    // Constructor is protected ... see below.
    
     // virtual destructor ensures use of polymorphism during destruction
    virtual ~QualityImprover() { };
    
    virtual double loop_over_mesh(MeshSet &ms, MsqError &err) = 0;

    //! provides a name to the QualityImprover (use it in constructor).
    void set_name(msq_std::string name)
      {
        qualityImproverName = name;
      };
    
    //! retrieves the QualityImprover name. A default name should be set in the constructor.
    virtual msq_std::string get_name() { return qualityImproverName; }
    virtual AlgorithmType get_algorithm_type() { return QUALITY_IMPROVER; }

      //!Sets in the termination criterion for the concrete solver's
      //! optimization.
    void set_inner_termination_criterion(TerminationCriterion* crit)
      {
        innerTerminationCriterion=crit;
      }
      //!Sets in the termination criterion for the outer loop over 
      //! patches.
    void set_outer_termination_criterion(TerminationCriterion* crit)
      {
        outerTerminationCriterion=crit;
      }

  protected:

    /*! The default constructor initialises a few member variables
        to default values.
        This can be reused by concrete class constructor. */    
    QualityImprover() : mMeshSet(0), qualityImproverName("noname")
      {
          //Temporary solution to not having an err object
        MsqError temp_err;
        defaultOuterCriterion.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,temp_err);
        outerTerminationCriterion = & defaultOuterCriterion;
        innerTerminationCriterion = & defaultInnerCriterion;
      }
    
    friend class MeshSet;
      //friend double QualityMetric::evaluate_element(MsqMeshEntity* element, MsqError &err);
      //friend double QualityMetric::evaluate_node(MsqNode* node, MsqError &err);
      //will not be needed when we remove stopping criterion
    const MeshSet* get_mesh_set() const
      { return mMeshSet; }
    MeshSet* get_mesh_set()
      { return mMeshSet; }
    
    
    void set_mesh_set(MeshSet *ms)
      {
        mMeshSet=ms;
      }
      //!return the outer termination criterion pointer 
    TerminationCriterion* get_outer_termination_criterion()
      { return outerTerminationCriterion; }
      //!return the inner termination criterion pointer       
    TerminationCriterion* get_inner_termination_criterion()
      { return innerTerminationCriterion; } 
    
  private:
    MeshSet* mMeshSet;
    msq_std::string qualityImproverName;
    int patchDepth;
    
    TerminationCriterion* innerTerminationCriterion;
    TerminationCriterion* outerTerminationCriterion;
      //default TerminationCriterion for outer loop will be set in constructor
    TerminationCriterion defaultOuterCriterion;
      //default TerminationCriterion for inner loop set by concrete improver
    TerminationCriterion defaultInnerCriterion;
  };

}

#endif
