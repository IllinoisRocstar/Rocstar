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
 
/*! \file LVQDTargetCalculators.hpp

Header file for the Mesquite::LVQDTargetCalculator class

  \author Thomas Leurent
  \date   2004-09-31
 */


#ifndef LVQDTargetCalculators_hpp
#define LVQDTargetCalculators_hpp

#include "Mesquite.hpp"
#include "TargetCalculator.hpp"

namespace Mesquite
{
  /*! \class LVQDTargetCalculator
    \brief This is an intermediary class.
    Concrete classes will simply instantiate the various guide enums in their constructor. 
  */
  class LVQDTargetCalculator : public TargetCalculator
  {
  protected:
    //! the constructor is protected since this is not intended to be a concrete class.
    //! Concrete classes will instantiate the various guide enums in their constructor. 
    LVQDTargetCalculator()
    { }

  public:
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~LVQDTargetCalculator()
      {};

    virtual void compute_target_matrices(PatchData& pd, MsqError& err);
      //! Function called by
      //! compute_target_matrices(PatchData& pd, MsqError& err) to compute
      //! the matrices after the reference PatchData has been created.
      //! \todo Michael:  this function should be protected, but it is used in
      //!  QualityMetricTest to avoid creating the MeshSet in the unit tests.
    void compute_target_matrices(PatchData& pd, PatchData& ref_pd,
                                 MsqError& err);


  protected:
    enum Lambda_type lambdaBase; 
    enum guide_type guideLambda;
    enum guide_type guideV;
    enum guide_type guideQ;
    enum guide_type guideDelta;

  };
  

} //namespace


#endif // LVQDTargetCalculator_hpp
