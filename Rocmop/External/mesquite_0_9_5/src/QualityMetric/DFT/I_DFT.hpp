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

/*! \file I_DFT.hpp

Header file for the Mesquite::I_DFT class

  \author Thomas Leurent
  \date   2004-09-29
 */


#ifndef I_DFT_hpp
#define I_DFT_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "DistanceFromTarget.hpp"
#include "Exponent.hpp"

namespace Mesquite
{
  
  /*! \class I_DFT
    \brief Class containing the target corner matrices for the context based smoothing.

The form of this metric is as follows (taken from I_DFTFamilyFunctions.hpp,
see that file for more detail):

                  mAlpha * || A*inv(W) - mBeta * I ||_F^2                    \n
   ------------------------------------------------------------------        \n
   0.5^(mGamma)*(det(A*inv(W)) + sqrt(det(A*inv(W))^2 + 4*delta^2))^(mGamma) \n

The default for data members (corresponding to the variables above):

      mAlpha = 1/2.0; \n
      mBeta = 1.0; \n
      mGamma = MSQ_TWO_THIRDS; \n

delta, above, is calculated using PatchData::get_barrier_delta(MsqError &err),
if useBarrierDelta == true.  Otherwise, delta is zero.

*/
  
  class I_DFT : public DistanceFromTarget
  {
  public:
    
    I_DFT()
    {
      MsqError err;
      set_averaging_method(LINEAR, err);
      set_metric_type(ELEMENT_BASED);
      set_gradient_type(ANALYTICAL_GRADIENT);
      set_hessian_type(ANALYTICAL_HESSIAN);
      set_name("I_DFT");
      mAlpha = 1/2.0;
      mBeta = 1.0;
      mGamma = MSQ_TWO_THIRDS;
      useBarrierDelta = true;
   }
    
    //! virtual destructor ensures use of polymorphism during destruction
    virtual ~I_DFT()
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
    
      //! access function to set mAlpha
    void p_set_alpha(double alpha)
      {mAlpha = alpha;}
      //! access function to get mAlpha  
    double p_get_alpha()
      {return mAlpha;}
      //! access function to set mBeta
    void p_set_beta(double beta)
      {mBeta = beta;}
      //! access function to get mBeta  
    double p_get_beta()
      {return mBeta;}
      //! access function to set mGamma
    void p_set_gamma(double gamma)
      {mGamma = gamma;}
      //! access function to get mGamma  
    double p_get_gamma()
      {return mGamma;}
      //! access function to set useBarrierDelta
    void p_set_use_barrier_delta(bool use_delta)
      {useBarrierDelta = use_delta;}
      //! access function to get useBarrierDelta  
    bool p_get_use_barrier_delta()
      {return useBarrierDelta;}

    
  private:
    // variables used in the definition of the metric (2d and 3d)
    double mAlpha;
    double mBeta;
    Exponent mGamma;
    bool useBarrierDelta;
    // variables used during the analytic gradient calculations
    Vector3D mNormals[4];		// Normal vector for merit function
    Vector3D mCoords[4]; 	// Vertex coordinates
    Vector3D mGrads[4];		// Gradients for element
    Vector3D mAccGrads[8];	// Accumulated gradients
    Matrix3D mHessians[10];	// Hessian values for element
    Matrix3D mR;		// R (in QR factorization of W)
    Matrix3D invR;		// Inverse matrix of R (in QR factorization)
    Matrix3D mQ;		// Q (in QR factorization of W)
  };
  
} //namespace


#endif // I_DFT_hpp
