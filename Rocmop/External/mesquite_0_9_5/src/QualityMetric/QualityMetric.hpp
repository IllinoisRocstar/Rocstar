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

/*! \file QualityMetric.hpp
    \brief
Header file for the Mesquite::QualityMetric class

  \author Thomas Leurent
  \author Michael Brewer
  \date   2002-05-01
 */

#ifndef QualityMetric_hpp
#define QualityMetric_hpp

#ifndef MSQ_USE_OLD_C_HEADERS
#include <cmath>
#include <cstring>
#else
#include <math.h>
#include <string.h>
#endif


#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"

namespace Mesquite
{
   
     /*! \class QualityMetric
       \brief Base class for concrete quality metrics.
     */
   class MsqVertex;
   class MsqMeshEntity;
   class PatchData;
   
   class QualityMetric
   {
   protected:
     /*!Constructor defaults concrete QualityMetric's to
       gradType=NUMERCIAL_GRADIENT and negateFlag=1.
       Concrete QualityMetric constructors over-write these defaults
       when appropriate.
       */
     QualityMetric() :
       evalMode(EEM_UNDEFINED),
       mType(MT_UNDEFINED),
       gradType(NUMERICAL_GRADIENT),
       hessianType(NUMERICAL_HESSIAN),
       negateFlag(1)
     {}

   public:
       // This is defined in each concrete class.  It isn't virtual, so
       // it doesn't exist in the base class.
       //   static void create_new() = 0;
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~QualityMetric()
        {};
     
     
       /*! \enum MetricType
       is a property of the metric. It should be set correctly in the constructor
       of the concrete QualityMetric.
       An example of a (mediocre) VERTEX_BASED metric is the smallest edge
       connected to a vertex.
       An example of a (mediocre) ELEMENT_BASED metric is the aspect ratio of an element.
       */
     enum MetricType
     {
        MT_UNDEFINED,
        VERTEX_BASED,
        ELEMENT_BASED,
        VERTEX_BASED_FREE_ONLY
     };
     
     MetricType get_metric_type() { return mType; }

     /*!\enum ElementEvaluationMode
       is for metrics of type ELEMENT_BASED only.
       It allows you to indicate whether we are evaluating
       the metric based on element vertices, or element gauss points. */
     enum ElementEvaluationMode
     {
       EEM_UNDEFINED,
       ELEMENT_VERTICES,
       LINEAR_GAUSS_POINTS,
       QUADRATIC_GAUSS_POINTS,
       CUBIC_GAUSS_POINTS
     };
     
     //!Sets the evaluation mode for the ELEMENT_BASED metrics.
     void set_element_evaluation_mode(ElementEvaluationMode mode, MsqError &err);
     
     //!Returns the evaluation mode for the metric
     inline ElementEvaluationMode  get_element_evaluation_mode()
        { return evalMode; }
     
       /*!AveragingMethod allows you to set how the quality metric values
         attained at each sample point will be averaged together to produce
         a single metric value for an element.
       */
     enum AveragingMethod
     {
        NONE,
        LINEAR,
        RMS,
        HMS,
        MINIMUM,
        MAXIMUM,
        HARMONIC,
        GEOMETRIC,
        SUM,
        SUM_SQUARED,
        GENERALIZED_MEAN,
        STANDARD_DEVIATION,
        MAX_OVER_MIN,
        MAX_MINUS_MIN,
        SUM_OF_RATIOS_SQUARED
     };
     
       /*!Set the averaging method for the quality metric. Current
         options are
         NONE: the values are not averaged,
         GEOMETRIC:  the geometric average,
         HARMONIC:  the harmonic average,
         LINEAR:  the linear average,
         MAXIMUM:  the maximum value,
         MINIMUM:  the minimum value,
         RMS:  the root-mean-squared average,
         HMS:  the harmonic-mean-squared average,
         SUM:  the sum of the values,
         SUM_SQUARED:  the sum of the squares of the values,
         GENERALIZED_MEAN: self explainatory,
         STANDARD_DEVIATION:  the standard deviation squared of the values,
         MAX_MINUS_MIN:  the maximum value minus the minum value,
         MAX_OVER_MIN:  the maximum value divided by the minimum value,
         SUM_OF_RATIOS_SQUARED:  (1/(N^2))*(SUM (SUM (v_i/v_j)^2))
       */
     inline void set_averaging_method(AveragingMethod method, MsqError &err);
     
       /*! Set feasible flag (i.e., does this metric have a feasible region
         that the mesh must maintain.)
       */
     inline void set_feasible_constraint(int alpha)
        { feasible=alpha; }
     
       //!Returns the feasible flag for this metric
     inline int get_feasible_constraint()
        { return feasible; }
     
       //!Sets the name of this metric
     inline void set_name(msq_std::string st)
        { metricName=st; }
     
       //!Returns the name of this metric (as a string).
     inline msq_std::string get_name()
        { return metricName; }

       //!Escobar Barrier Function for Shape and Other Metrics
       // det = signed determinant of Jacobian Matrix at a Vertex
       // delta = scaling parameter
     inline double vertex_barrier_function(double det, double delta) 
            { return 0.5*(det+sqrt(det*det+4*delta*delta)); }
     
       //!Evaluate the metric for a vertex
     virtual bool evaluate_vertex(PatchData& /*pd*/, MsqVertex* /*vertex*/,
                                  double& /*value*/, MsqError &err);
     
       //!Evaluate the metric for an element
     virtual bool evaluate_element(PatchData& /*pd*/,
                                   MsqMeshEntity* /*element*/,
                                   double& /*value*/, MsqError &err);
     
       /*!\enum GRADIENT_TYPE Sets to either NUMERICAL_GRADIENT or
         ANALYTICAL_GRADIENT*/
     enum GRADIENT_TYPE
     {
        NUMERICAL_GRADIENT,
        ANALYTICAL_GRADIENT
     };
     
       //!Sets gradType for this metric.
     void set_gradient_type(GRADIENT_TYPE grad)
        { gradType=grad; }
     
       /*!\enum HESSIAN_TYPE Sets to either NUMERICAL_HESSIAN or
         ANALYTICAL_HESSIAN*/
     enum HESSIAN_TYPE
     {
        NUMERICAL_HESSIAN,
        ANALYTICAL_HESSIAN
     };
     
       //!Sets hessianType for this metric.
     void set_hessian_type(HESSIAN_TYPE ht)
        { hessianType=ht; }
     
       /*!For MetricType == VERTEX_BASED.
         Calls either compute_vertex_numerical_gradient or
         compute_vertex_analytical_gradient for gradType equal
         NUMERICAL_GRADIENT or ANALYTICAL_GRADIENT, respectively.

         \return true if the element is valid, false otherwise. 
       */
     bool compute_vertex_gradient(PatchData &pd,MsqVertex &vertex,
                                  MsqVertex* vertices[],Vector3D grad_vec[],
                                  int num_vtx, double &metric_value,
                                  MsqError &err);
     
       /*! \brief For MetricType == ELEMENT_BASED.
         Calls either compute_element_numerical_gradient() or
         compute_element_analytical_gradient() for gradType equal
         NUMERICAL_GRADIENT or ANALYTICAL_GRADIENT, respectively.
       */
     bool compute_element_gradient(PatchData &pd, MsqMeshEntity* element,
                                   MsqVertex* free_vtces[], Vector3D grad_vec[],
                                   int num_free_vtx, double &metric_value, MsqError &err);

     /*! same as compute_element_gradient(), but fills fixed vertices spots with
       zeros instead of not returning values for fixed vertices. Also, the vertices
       are now ordered according to the element vertices array.
       */
     bool compute_element_gradient_expanded(PatchData &pd, MsqMeshEntity* element,
                                   MsqVertex* free_vtces[], Vector3D grad_vec[],
                                   int num_free_vtx, double &metric_value, MsqError &err);
     
       /*! \brief For MetricType == ELEMENT_BASED.
         Calls either compute_element_numerical_hessian() or
         compute_element_analytical_hessian() for hessianType equal
         NUMERICAL_HESSIAN or ANALYTICAL_HESSIAN, respectively.
       */
     bool compute_element_hessian(PatchData &pd, MsqMeshEntity* element,
                                  MsqVertex* free_vtces[], Vector3D grad_vec[],
                                  Matrix3D hessian[],
                                  int num_free_vtx, double &metric_value, MsqError &err);
     
       /*! Set the value of QualityMetric's negateFlag.  Concrete
         QualityMetrics should set this flag to -1 if the QualityMetric
         needs to be maximized.
       */
     void set_negate_flag(int neg)
        { negateFlag=neg; }
     
       //!Returns negateFlag.
     int get_negate_flag()
        { return negateFlag; }
       /*! This function is user accessible and virtual.  The base
         class implementation sets an error, because many metrics
         will only be defined as Element_based or Vertex_based, and
         this function will not be needed.  Some concrete metrics
         will have both Element_based and Vertex_based definintions,
         and those metrics will re-implement this function to the
         MetricType to be changed to either QualityMetric::VERTEX_BASED
         or QualityMetric::ELEMENT_BASED.*/
     virtual void change_metric_type(MetricType t, MsqError &err);
     
     
  protected:
     
     //! This function should be used in the constructor of every concrete
     //! quality metric. Errors will result if type is left to MT_UNDEFINED.
     void set_metric_type(MetricType t) { mType = t; }
     
     //! average_metrics takes an array of length num_values and averages the
     //! contents using averaging method data member avgMethod .
     double average_metrics(const double metric_values[], const int& num_values,
                            MsqError &err);
                            
     //! Given a list of metric values, calculate the average metric
     //! valude according to the current avgMethod and write into
     //! the passed metric_values array the the value weight/count to
     //! use when averaging gradient vectors for the metric.
     //!\param metric_values : As input, a set of quality metric values
     //!                       to average.  As output, the fraction of
     //!                       the corresponding gradient vector that
     //!                       contributes to the average gradient.
     //!\param num_metric_values The number of values in the passed array.
     double average_metric_and_weights( double metric_values[],
                                        int num_metric_values,
                                        MsqError& err );
     
     //! takes an array of coefficients and an array of metrics (both of length num_value)
     //! and averages the contents using averaging method 'method'.
     double weighted_average_metrics(const double coef[],
                                    const double metric_values[],
                                    const int& num_values, MsqError &err);
     
       /*!Non-virtual function which numerically computes the gradient
         of a QualityMetric of a given free vertex. This is used by metric
         which mType is VERTEX_BASED. 
         \return true if the element is valid, false otherwise. */
     bool compute_vertex_numerical_gradient(PatchData &pd,
                                            MsqVertex &vertex,
                                            MsqVertex* vertices[],
                                            Vector3D grad_vec[],
                                            int num_vtx,
                                            double &metric_value,
                                            MsqError &err);
     
     
       /*!\brief Non-virtual function which numerically computes the gradient
         of a QualityMetric of a given element for a given set of free vertices
         on that element.
         This is used by metric which mType is ELEMENT_BASED.
         For parameters, see compute_element_gradient() . */
     bool compute_element_numerical_gradient(PatchData &pd, MsqMeshEntity* element,
                                             MsqVertex* free_vtces[], Vector3D grad_vec[],
                                             int num_free_vtx, double &metric_value,
                                             MsqError &err);

     /*! \brief Virtual function that computes the gradient of the QualityMetric
         analytically.  The base class implementation of this function
         simply prints a warning and calls compute_numerical_gradient
         to calculate the gradient. This is used by metric
         which mType is VERTEX_BASED. */
     virtual bool compute_vertex_analytical_gradient(PatchData &pd,
                                                     MsqVertex &vertex,
                                                     MsqVertex* vertices[],
                                                     Vector3D grad_vec[],
                                                     int num_vtx,
                                                     double &metric_value,
                                                     MsqError &err);
     
     
     /*! \brief Virtual function that computes the gradient of the QualityMetric
         analytically.  The base class implementation of this function
         simply prints a warning and calls compute_numerical_gradient
         to calculate the gradient. This is used by metric
         which mType is ELEMENT_BASED.
         For parameters, see compute_element_gradient() . */
     virtual bool compute_element_analytical_gradient(PatchData &pd,
                                                      MsqMeshEntity* element,
                                                      MsqVertex* free_vtces[],
                                                      Vector3D grad_vec[],
                                                      int num_free_vtx,
                                                      double &metric_value,
                                                      MsqError &err);


     bool compute_element_numerical_hessian(PatchData &pd,
                                            MsqMeshEntity* element,
                                            MsqVertex* free_vtces[],
                                            Vector3D grad_vec[],
                                            Matrix3D hessian[],
                                            int num_free_vtx,
                                            double &metric_value,
                                            MsqError &err);

     virtual bool compute_element_analytical_hessian(PatchData &pd,
                                            MsqMeshEntity* element,
                                            MsqVertex* free_vtces[],
                                            Vector3D grad_vec[],
                                            Matrix3D hessian[],
                                            int num_free_vtx,
                                            double &metric_value,
                                            MsqError &err);

     friend class MsqMeshEntity;

     // TODO : pass this private and write protected access fucntions.
     AveragingMethod avgMethod;
     int feasible;
     msq_std::string metricName;
  private:
     ElementEvaluationMode evalMode;
     MetricType mType;
     GRADIENT_TYPE gradType;
     HESSIAN_TYPE hessianType;
     int negateFlag;
   };

  
  inline void QualityMetric::set_element_evaluation_mode(ElementEvaluationMode mode, MsqError &err)
  {
    if (mType == VERTEX_BASED) {
      MSQ_SETERR(err)("function must only be used for ELEMENT_BASED metrics.", MsqError::INVALID_STATE);
      return;
    }
    
    switch(mode)
      {
      case(ELEMENT_VERTICES):
      case(LINEAR_GAUSS_POINTS):
      case(QUADRATIC_GAUSS_POINTS):
      case(CUBIC_GAUSS_POINTS):
        evalMode=mode;
        break;
      default:
        MSQ_SETERR(err)("Requested mode not implemented", MsqError::NOT_IMPLEMENTED);
      }
    return;
  }

  
  inline void  QualityMetric::set_averaging_method(AveragingMethod method, MsqError &err)
  {
    switch(method)
    {
      case(NONE):
      case(GEOMETRIC):
      case(HARMONIC):
      case(LINEAR):
      case(MAXIMUM):
      case(MINIMUM):
      case(RMS):
      case(HMS):
      case(STANDARD_DEVIATION):
      case(SUM):
      case(SUM_SQUARED):
      case(MAX_OVER_MIN):
      case(MAX_MINUS_MIN):
      case(SUM_OF_RATIOS_SQUARED):
        avgMethod=method;
        break;
      default:
       MSQ_SETERR(err)("Requested Averaging Method Not Implemented", MsqError::NOT_IMPLEMENTED);
      };
    return;
  }
  

/*! 
  \brief Calls compute_vertex_numerical_gradient if gradType equals
  NUMERCIAL_GRADIENT.  Calls compute_vertex_analytical_gradient if 
  gradType equals ANALYTICAL_GRADIENT;
*/
   inline bool QualityMetric::compute_vertex_gradient(PatchData &pd,
                                                      MsqVertex &vertex,
                                                      MsqVertex* vertices[],
                                                      Vector3D grad_vec[],
                                                      int num_vtx,
                                                      double &metric_value,
                                                      MsqError &err)
   {
     bool ret=false;;
     switch(gradType)
     {
       case NUMERICAL_GRADIENT:
          ret = compute_vertex_numerical_gradient(pd, vertex, vertices,
                                                  grad_vec, num_vtx,
                                                  metric_value, err);
          MSQ_CHKERR(err);
          break;
       case ANALYTICAL_GRADIENT:
          ret = compute_vertex_analytical_gradient(pd, vertex, vertices,
                                                   grad_vec,num_vtx,
                                                   metric_value, err);
          MSQ_CHKERR(err);
          break;
     }
     return ret;
   }
   

/*! 
    \param free_vtces base address of an array of pointers to the element vertices which
    are considered free for purposes of computing the gradient. The quality metric
    gradient relative to each of those vertices is computed and stored in grad_vec.
    \param grad_vec base address of an array of Vector3D where the gradient is stored,
    in the order specified by the free_vtces array.
    \param num_free_vtx This is the size of the vertices and gradient arrays. 
    \param metric_value Since the metric is computed, we return it. 
    \return true if the element is valid, false otherwise.
*/
   inline bool QualityMetric::compute_element_gradient(PatchData &pd,
                                                       MsqMeshEntity* el,
                                                       MsqVertex* free_vtces[],
                                                       Vector3D grad_vec[],
                                                       int num_free_vtx,
                                                       double &metric_value,
                                                       MsqError &err)
   {
     bool ret=false;
     switch(gradType)
     {
       case NUMERICAL_GRADIENT:
          ret = compute_element_numerical_gradient(pd, el, free_vtces, grad_vec,
                                                  num_free_vtx, metric_value, err);
          MSQ_CHKERR(err);
          break;
       case ANALYTICAL_GRADIENT:
          ret = compute_element_analytical_gradient(pd, el, free_vtces, grad_vec,
                                                   num_free_vtx, metric_value, err);
          MSQ_CHKERR(err);
          break;
     }
     return ret;
   }
   

     /*! 
       average_metrics takes an array of length num_value and averages the
       contents using averaging method 'method'.
     */
   inline double QualityMetric::average_metrics(const double metric_values[],
                                           const int& num_values, MsqError &err)
   {
       //MSQ_MAX needs to be made global?
     //double MSQ_MAX=1e10;
     double total_value=0.0;
     double temp_value=0.0;
     int i=0;
     int j=0;
       //if no values, return zero
     if (num_values<=0){
       return 0.0;
     }
     
     switch(avgMethod){
       case GEOMETRIC:
          total_value=1.0;
          for (i=0;i<num_values;++i){
            total_value*=(metric_values[i]);
          }
          total_value=pow(total_value, (1/((double) num_values)));
          break;
          
       case HARMONIC:
            //ensure no divide by zero, return zero
          for (i=0;i<num_values;++i){
            if(metric_values[i]<MSQ_MIN){
              if(metric_values[i]>MSQ_MIN){
                return 0.0;
              }
            }
            total_value+=(1/metric_values[i]);
          }
            //ensure no divide by zero, return MSQ_MAX_CAP
          if(total_value<MSQ_MIN){
            if(total_value>MSQ_MIN){
              return MSQ_MAX_CAP;
            }
          }
          total_value=num_values/total_value;
          break;
          
       case LINEAR:
          for (i=0;i<num_values;++i){
            total_value+=metric_values[i];
          }
          total_value/= (double) num_values;
          break;
          
       case MAXIMUM:
          total_value = metric_values[0];
          for (i = 1; i < num_values; ++i){
            if (metric_values[i] > total_value){
              total_value = metric_values[i];
            }
          }
          break;
          
       case MINIMUM:
          total_value = metric_values[0];
          for (i = 1; i < num_values; ++i){
            if (metric_values[i] < total_value) {
              total_value = metric_values[i];
            }
          }
          break;
          
       case NONE:
          MSQ_SETERR(err)("Averaging method set to NONE", MsqError::INVALID_ARG);
          break;
          
       case RMS:
          for (i=0;i<num_values;++i){
            total_value+=(metric_values[i]*metric_values[i]);
          }
          total_value/= (double) num_values;
          total_value=sqrt(total_value);
          break;
          
       case HMS:
          //ensure no divide by zero, return zero
          for (i=0;i<num_values;++i){
            if (metric_values[i]*metric_values[i] < MSQ_MIN) {
              return 0.0;
            }
            total_value += (1.0/(metric_values[i]*metric_values[i]));
          }

          //ensure no divide by zero, return MSQ_MAX_CAP
          if (total_value < MSQ_MIN) {
            return MSQ_MAX_CAP;
          }
          total_value = sqrt(num_values/total_value);
          break;

       case STANDARD_DEVIATION:
          total_value=0;
          temp_value=0;
          for (i=0;i<num_values;++i){
            temp_value+=metric_values[i];
            total_value+=(metric_values[i]*metric_values[i]);
          }
          temp_value/= (double) num_values;
          temp_value*=temp_value;
          total_value/= (double) num_values;
          total_value-=temp_value;
          break;

       case SUM:
          for (i=0;i<num_values;++i){
            total_value+=metric_values[i];
          }
          break;

       case SUM_SQUARED:
          for (i=0;i<num_values;++i){
            total_value+= (metric_values[i]*metric_values[i]);
          }
          break;

       case MAX_MINUS_MIN:
          //total_value used to store the maximum
            //temp_value used to store the minimum
          temp_value=MSQ_MAX_CAP;
          for (i=0;i<num_values;++i){
            if(metric_values[i]<temp_value){
              temp_value=metric_values[i];
            }
            if(metric_values[i]>total_value){
              total_value=metric_values[i];
            }
          }
          
            //ensure no divide by zero, return MSQ_MAX_CAP
          if (temp_value < MSQ_MIN) {
            return MSQ_MAX_CAP;
          }
          total_value-=temp_value;
          break;

       case MAX_OVER_MIN:
            //total_value used to store the maximum
            //temp_value used to store the minimum
          temp_value=MSQ_MAX_CAP;
          for (i=0;i<num_values;++i){
            if(metric_values[i]<temp_value){
              temp_value=metric_values[i];
            }
            if(metric_values[i]>total_value){
              total_value=metric_values[i];
            }
          }
          
            //ensure no divide by zero, return MSQ_MAX_CAP
          if (temp_value < MSQ_MIN) {
            return MSQ_MAX_CAP;
          }
          total_value/=temp_value;
          break;
          
       case SUM_OF_RATIOS_SQUARED:
          for (j=0;j<num_values;++j){
            //ensure no divide by zero, return MSQ_MAX_CAP
            if (metric_values[j] < MSQ_MIN) {
              return MSQ_MAX_CAP;
            }
            for (i=0;i<num_values;++i){
              total_value+=((metric_values[i]/metric_values[j])*
                            (metric_values[i]/metric_values[j]));
            }
          }
          total_value/=(double)(num_values*num_values);
          break;
          
       default:
            //Return error saying Averaging Method mode not implemented
          MSQ_SETERR(err)("Requested Averaging Method Not Implemented", MsqError::NOT_IMPLEMENTED);
          return 0;
     }
     return total_value;
   }
   
   

  inline double QualityMetric::weighted_average_metrics(const double coef[],
                                                        const double metric_values[],
                                                        const int& num_values, MsqError &err)
  {
    //MSQ_MAX needs to be made global?
    //double MSQ_MAX=1e10;
    double total_value=0.0;
    int i=0;
    //if no values, return zero
    if (num_values<=0){
      return 0.0;
    }
     
    switch(avgMethod){

    case LINEAR:
      for (i=0;i<num_values;++i){
        total_value += coef[i]*metric_values[i];
      }
      total_value /= (double) num_values;
      break;
          
    default:
      //Return error saying Averaging Method mode not implemented
      MSQ_SETERR(err)("Requested Averaging Method Not Implemented",MsqError::NOT_IMPLEMENTED);
      return 0;
    }
    return total_value;
  }
   
   
} //namespace


#endif // QualityMetric_hpp
