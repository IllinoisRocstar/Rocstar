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

/*! \file ObjectiveFunction.hpp

Header file for the Mesquite::ObjectiveFunction class

  \author Michael Brewer
  \author Thomas Leurent
  \date   2002-05-23
 */


#ifndef OBJECTIVE_FUNCTION_HPP
#define OBJECTIVE_FUNCTION_HPP

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#else
#  include <list>
#endif

namespace Mesquite
{
   class PatchData;
   class MsqHessian;
   class QualityMetric;
   
  /*! \class ObjectiveFunction
       \brief Base class for concrete Objective Functions
       ObjectiveFunction contains a pointer to a QualityMetric.  If
       the ObjectiveFunction is associated with more than one
       QualityMetric (i.e., the Objective is a composite, and the
       composed ObjectiveFunctions are associated with different
       QualityMetrics), then the QualityMetric pointer is set
       to NULL..
      */
  class ObjectiveFunction
  {
   public:
    ObjectiveFunction ():
        useLocalGradient(false)
       {}
    
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~ObjectiveFunction()
       {};
    
      /*! 
        Evaluate the objective function on a given patch.
      */
    virtual bool concrete_evaluate(PatchData &patch, double &fval,
                                   MsqError &err)=0;

      /*!
        Computes the value of the objective funciton as fval*negateFlag,
        where fval is computed in concrete_evaluate(patch, fval, err)
        and negateFlag is either 1 or -1 depending on whether the
        function needs to be minimized or maximized, respectively.
        Returns the bool given as a return value from concrete_evaluate.
        If the bool is 'false', then the patch is not within the feasible
        region required by the associated QualityMetric(s).
      */
    bool evaluate(PatchData &patch, double &fval, MsqError &err)
       {
         bool return_bool = concrete_evaluate(patch, fval, err);
         fval *= negateFlag;
         return return_bool;
       }
    
    //! \enum GRADIENT_TYPE 
    enum GRADIENT_TYPE
    {
       NUMERICAL_GRADIENT,  //!< can be very slow. Should be for tests only. 
       ANALYTICAL_GRADIENT  //!< every differentiable function should have
                            //!< an analytical gradient implemented.
    };

    //! Set gradType to either NUMERICAL_GRADIENT or ANALYTICAL_GRADIENT.
    void set_gradient_type(GRADIENT_TYPE grad)
    { gradType = grad; }


      /*!\brief
        Calls either compute_numerical_gradient or compute_analytical_gradient
        depending on the value of gradType.
        Function returns 'false' if the patch is not within a required
        feasible regeion.  Otherwise, it returns 'true'.
      */
    bool compute_gradient(PatchData &patch, Vector3D *const &grad,
                          double &OF_val, MsqError &err, size_t array_size=0);          
              
      /*!\brief
        Calls compute_analytical_hessian. 
        Function returns 'false' if the patch is not within a required
        feasible regeion.  Otherwise, it returns 'true'.
      */
    bool compute_hessian(PatchData &patch, MsqHessian &hessian,
                         Vector3D *const &grad,
                         double &OF_val,
                         MsqError &err);          
              
  /*! 
        Return the quality metric associated with this objective function.
        Returns null for composite functions which have multiple associated
        quality metrics.  Use get_quality_metric_list() to retrieve all
        metrics.
      */
    QualityMetric*  get_quality_metric(){
         return qMetric;
       }
    
      /*! 
        returns a list of all associated metrics;
      */
    virtual msq_std::list<QualityMetric*> get_quality_metric_list()
       {
         msq_std::list<QualityMetric*> temp_list;
         temp_list.push_front(qMetric);
         return temp_list;
       }

      //!Set the value of qMetric.
    void set_quality_metric(QualityMetric* qm)
       {
         qMetric=qm;
       }
    
      /*! 
        \brief Set the value of ObjectiveFunction's negateFlag.  Unless
        composite, concrete ObjectiveFunctions should set this flag to
        to the value of the associated QualityMetric's negateFLag.
      */
    void set_negate_flag(int neg)
       {
         negateFlag=neg;
       }

      //!Returns negateFlag
    int get_negate_flag()
       {
         return negateFlag;
       }
    
   protected:

      /*! \brief Non-virtual function which numerically computes
        the gradient of the Objective Function.
         Function returns 'false' if the patch is not within a required
        feasible regeion.  Otherwise, it returns 'true'.
      */
    bool compute_numerical_gradient(PatchData &patch, Vector3D *const &grad,
                                    double &OF_val,
                                    MsqError &err, size_t array_size);

     /*! 
        Fills an array of Vector3D, grad, with the gradient of
        the objective function computed using the gradient of the
        quality metric.  If the function has not been over-riden
        in the concrete Objective Function, the base class implementation
        prints a warning and then defaults to numerical gradient.
        Function returns 'false' if the patch is not within a required
        feasible regeion.  Otherwise, it returns 'true'.

        \param patch The PatchData object for which the objective function
        gradient is computed.
        \param grad An array of Vector3D, at least the size of the number
        of vertices in the patch.
        \param OF_val is set to the value of the objective function.
        \param array_size is the size of the grad Vector3D[] array and
        must correspond to the number of vertices in the patch.
     */
    virtual bool compute_analytical_gradient(PatchData &patch,
                                             Vector3D *const &grad,
                                             double &OF_val,
                                             MsqError &err, size_t array_size);
    
     /*! 
        Fills a MsqHessian object with the Hessian of
        the objective function computed using the hessian of the
        quality metric.  If the function has not been over-riden
        in the concrete Objective Function, the base class implementation
        prints a warning and returns false.
        Function returns 'false' if the patch is not within a required
        feasible regeion.  Otherwise, it returns 'true'.
      */
    virtual bool compute_analytical_hessian(PatchData &/*patch*/,
                                            MsqHessian &/*hessian*/,
                                            Vector3D *const &/*grad*/,
                                            double &/*OF_val*/,
                                            MsqError &/*err*/);
    
      //!Returns eps used in the numerical gradient calculation.
    inline double get_eps(PatchData &pd, double &local_val,
                          int k,MsqVertex* vertex, MsqError &err);
    
      //!Sets useLocalGradient
      //!This variable determines whether compute_numercial_gradient
      //!can use the most efficient gradient calculation.
    void set_use_local_gradient(bool new_bool)
       {
         useLocalGradient=new_bool;
       }
       
    
 private:
      
    enum GRADIENT_TYPE gradType;//!Flag for numerical or analytical gradient.
      
    QualityMetric* qMetric;//!Pointer to associated QualityMetric.
      
    int negateFlag; /*!Equals one if ObjectiveFunction needs to
        be minimized; equals negative one if ObjectiveFunction needs
        to be maximized.*/
    bool useLocalGradient;/*!Set to true if we should use the more efficient
                            sub-patch method of computing the numerical
                            gradient.  Otherwise, it is set to false.*/
    
  };

//BEGIN INLINE

     /*!  
       Calls either compute_numerical_gradient or compute_analytical_gradient
       depending on the value of gradType.           
      */
   inline bool ObjectiveFunction::compute_gradient(PatchData &patch,
                                                   Vector3D *const &grad,
                                                   double &OF_val,
                                                   MsqError &err,
                                                   size_t array_size)
   {
     bool obj_bool = false;
     switch(gradType){
       case NUMERICAL_GRADIENT:
          obj_bool=compute_numerical_gradient(patch, grad, OF_val,
                                              err, array_size);
          break;
       case ANALYTICAL_GRADIENT:
          obj_bool=compute_analytical_gradient(patch, grad, OF_val,
                                               err, array_size);
          break;
     }
     return !MSQ_CHKERR(err) && obj_bool;
   }


     /*!  
       Calls compute_analytical_hessian.
       Numerical objective function hessians are only used for test purposes. 
     */
   inline bool ObjectiveFunction::compute_hessian(PatchData &patch,
                                                  MsqHessian &hessian,
                                                  Vector3D *const &grad,
                                                  double &OF_val,
                                                  MsqError &err)
   {
     bool result = compute_analytical_hessian(patch, hessian,
                                       grad, OF_val, err);
     return !MSQ_CHKERR(err) && result;
   }


  /*!Returns an appropiate value (eps) to use as a delta step for
    MsqVertex vertex in dimension k (i.e. k=0 -> x, k=1 -> y, k=2 -> z).
    The objective function value at the perturbed vertex position is given
    in local_val.
  */
  inline double ObjectiveFunction::get_eps(PatchData &pd, double &local_val,
                                           int k,MsqVertex* vertex, MsqError& err)
  {
    double eps = 1.e-07;
    //  double rho=.5;
    int imax=20;
    int i=0;
    bool feasible=false;
    double tmp_var=0.0;
    while (i<imax && !feasible)
      {
        i++;
        //perturb kth coord val and check feas if needed
        tmp_var=(*vertex)[k];
        (*vertex)[k]+=eps;
        feasible = evaluate(pd,local_val,err); MSQ_ERRZERO(err);
        //if step was too big, shorten it         
        if(!feasible)
          eps*=0.5;
        //revert kth coord val
        (*vertex)[k]=tmp_var;
      }//end while looking for feasible eps
    return eps;
  }//end function get_eps

  
} //namespace


#endif // ObjectiveFunction_hpp


