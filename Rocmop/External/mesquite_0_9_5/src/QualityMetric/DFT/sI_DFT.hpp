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

/*! \file sI_DFT.hpp

Header file for the Mesquite::sI_DFT class

  \author Thomas Leurent
  \date   2004-09-29
 */


#ifndef sI_DFT_hpp
#define sI_DFT_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "DistanceFromTarget.hpp"

namespace Mesquite
{
  
  /*! \class sI_DFT
    \brief Class containing the target corner matrices for the context based smoothing. 
  */
  class sI_DFT : public DistanceFromTarget
  {
  public:
    
    sI_DFT()
    {
      MsqError err;
      set_averaging_method(LINEAR, err); 
      set_metric_type(ELEMENT_BASED);
      set_gradient_type(NUMERICAL_GRADIENT);
      set_hessian_type(NUMERICAL_HESSIAN);
    }
    
    //! virtual destructor ensures use of polymorphism during destruction
    virtual ~sI_DFT()
       {};

    virtual bool evaluate_element(PatchData& pd,
                            MsqMeshEntity* element,
                            double& value, MsqError &err);
    
  protected:
 
  private:
    
  };

} //namespace


#endif // sI_DFT_hpp
