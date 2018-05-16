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

/*! \file UntangleBetaQualityMetric.hpp

Header file for the Mesquite::UntangleBetaQualityMetric class

  \author Michael Brewer
  \date   2002-09-10
 */

#ifndef UNTANGLE_BETA_QUALITY_METRIC_HPP
#define UNTANGLE_BETA_QUALITY_METRIC_HPP

#include "Mesquite.hpp"
#include "UntangleQualityMetric.hpp"
#include "PatchData.hpp"
namespace Mesquite
{
     /*! \class UntangleBetaQualityMetric
       \brief The untangle beta quality metric.
       
       Given a scalar value beta and local signed element volume alpha_i,
       define delta_i to be alpha_i minus beta.  The Untangle beta value
       is then defined as square root of the sum over sample points
       of the absolute value of delta_i minus delta_i, difference squared.
       That is, the root mean square of the difference, abs(delta_i) minus
       delta_i.

       The constructor defaults to RMS AveragingMethod and
       ELEMENT_VERTICES evaluationMode.  The default beta value is
       .05.
     */
   
   class UntangleBetaQualityMetric : public UntangleQualityMetric
   {
   public:
     
     UntangleBetaQualityMetric(double bet=0.05);

       // virtual destructor ensures use of polymorphism during destruction
     virtual ~UntangleBetaQualityMetric()
        {}
       /*!Evaluate the Untangle Beta metric value for an element.
         \todo This function needs to be modified so that it no longer
         uses compute_weighted_jacobian.  It also needs to set an error
         whenever sent a 2D element and the surface normal information
         is not available.*/
     bool evaluate_element(PatchData &pd,MsqMeshEntity* element,
                           double &fval,MsqError &err);
     
       /*!Function to allow users to set the beta value after the
         metric has already been created. */
     void set_beta(double beta_in)
       {mBeta = beta_in;}
       /*!Function to allow the user to check the value of beta.*/
     double get_beta()
       {return mBeta;}
     
   protected:
     inline void untangle_function_2d(Vector3D temp_vec[],size_t e_ind,
                                      PatchData &pd, double &fval,
                                      MsqError &err);
     
     inline void untangle_function_3d(Vector3D temp_vec[],double &fval,
                                      MsqError &err);
     
   private:
     double mBeta;
   };
     //************BEGIN INLINE FUNCTIONS**************************
   
   inline void UntangleBetaQualityMetric::untangle_function_2d(Vector3D temp_vec[],size_t e_ind,PatchData &pd, double &fval, MsqError &err)
   {
     Vector3D surface_normal;
     pd.get_domain_normal_at_element(e_ind,surface_normal,err); MSQ_ERRRTN(err);
     Vector3D cross_vec=temp_vec[0]*temp_vec[1];
       //cout<<"\nsurface_normal "<<surface_normal;
       //cout<<"\cross_vec "<<cross_vec;
     double temp_var=cross_vec.length();
     if(cross_vec%surface_normal<0.0){
       temp_var*=-1;
     }
     temp_var -= mBeta;
       //cout<<"temp_var == "<<temp_var;
     fval=0.0;
     if(temp_var<0.0){
        fval=fabs(temp_var)-temp_var;
     }
       //  cout<<"\nfval == "<<fval<<"  e_ind "<<e_ind;
   }
   
   inline void UntangleBetaQualityMetric::untangle_function_3d(Vector3D temp_vec[],double &fval, MsqError &/*err*/)
   {
     double temp_var=temp_vec[0]%(temp_vec[1]*temp_vec[2]);
     temp_var-=mBeta;
     fval=0.0;
     if(temp_var<0.0){
       fval=fabs(temp_var)-temp_var;
     }
   }
   
} //namespace


#endif // UntangleBetaQualityMetric_hpp


