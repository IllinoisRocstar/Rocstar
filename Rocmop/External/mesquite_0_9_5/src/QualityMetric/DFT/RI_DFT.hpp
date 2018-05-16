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

/*! \file RI_DFT.hpp

Header file for the Mesquite::RI_DFT class

  \author Thomas Leurent
  \date   2004-09-29
 */


#ifndef RI_DFT_hpp
#define RI_DFT_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "DistanceFromTarget.hpp"
#include "Exponent.hpp"

namespace Mesquite
{
  
  /*! \class RI_DFT
    \brief Class containing the target corner matrices for the context based smoothing. 
  */
  class RI_DFT : public DistanceFromTarget
  {
  public:

    RI_DFT()
      : a(pow(2.0, MSQ_ONE_THIRD)), b(1.0), c(-4.0/3.0)
    {
      MsqError err;
      set_averaging_method(LINEAR, err); MSQ_CHKERR(err);
      set_metric_type(ELEMENT_BASED);
      set_gradient_type(NUMERICAL_GRADIENT);
      set_hessian_type(NUMERICAL_HESSIAN);
    }
    
    //! virtual destructor ensures use of polymorphism during destruction
    virtual ~RI_DFT()
       {};

    bool evaluate_element(PatchData &pd,
			  MsqMeshEntity *e,
			  double &m, MsqError &err);
    
    bool compute_element_analytical_gradient(PatchData &pd,
					     MsqMeshEntity *e,
					     MsqVertex *fv[], 
					     Vector3D g[],
					     int nfv, 
					     double &m,
					     MsqError &err);

    bool compute_element_analytical_hessian(PatchData &pd,
					    MsqMeshEntity *e,
					    MsqVertex *fv[], 
					    Vector3D g[],
					    Matrix3D h[],
					    int nfv, 
					    double &m,
					    MsqError &err);

  protected:
 
  private:
    // variables used in the definition of the metric (2d and 3d)
    double a;
    Exponent b;
    Exponent c;

    // variables used during the analytic gradient calculations
    Vector3D mNormal;		// Normal vector for merit function
    Vector3D mCoords[4]; 	// Vertex coordinates
    Vector3D mGrads[4];		// Gradients for element
    Vector3D mAccGrads[8];	// Accumulated gradients
    Matrix3D mHessians[10];	// Hessian values for element
    Matrix3D invW;		// Inverse matrix
  };

  
} //namespace


#endif // RI_DFT_hpp
