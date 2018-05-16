/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************/
// $Id: Algebraic_Metrics_3.h,v 1.5 2008/12/06 08:45:23 mtcampbe Exp $

/*! \file Algebraic_Metrics_3.h 
    \brief 3d algebraic quality %Metric declarations.
*/

#ifndef __ALGEBRAIC_METRICS_3_H__
#define __ALGEBRAIC_METRICS_3_H__

#include "Metric.h"

MOP_BEGIN_NAMESPACE

//! 3D Algebraic %Metric Base Class
/** 
 * This class is the base class for the implementation of
 * the specific 3D algebraic metrics defined by Patrick Knupp in
 * "Aglebraic Mesh Quality Metrics for Unstructured Initial Meshes"
 */ 
class Alg_Metric_Base_3 : public Metric {
  friend std::ostream & operator<<(std::ostream &, 
				   const Alg_Metric_Base_3 &);
public:

  //! Constructor
  Alg_Metric_Base_3() {}

  //! Virtual Destructor
  /**
   *  Although the destructor does nothing, this member function must be declared
   *  since the class has other virtual functions.
   */  
  virtual ~Alg_Metric_Base_3() {}  
  
  //! Initialize a 3D Algebraic Metric
  /**
   * \param n[] Ordered set of nodal coordinates.
   * \param type Element type: TET or QUAD
   */
  virtual void initialize(Vector_3<double> n[], int type);

  virtual void initialize(Element_node_enumerator & ene);

  //! The maximum value for this metric.
  virtual double maxValue() const { return 1.0; }

  //! The minimum value for this metric.
  virtual double minValue() const { return 0.0; }

protected:

  //! Compute the size metric.
  /**
   * \param ref_area the ideal area for this element.
   */
  double compute_size(double ref_voume=1.) const;

  //! Compute the shape metric.
  double compute_shape() const;

  //! Compute the skew metric.
  double compute_skew() const;

protected:
  double   alpha[8];
  J_Matrix A[8];
  Matrix   L[8];
  int type_;
};

//! 3D Shape %Metric Class
class Shape_Metric_3 : public Alg_Metric_Base_3 {
public:

  //! Construct a 3D Shape metric
  Shape_Metric_3() {}

  //! Calculate the metric value on this element.
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;
};

//! 3D Size %Metric Class
class Size_Metric_3 : public Alg_Metric_Base_3 {
public:
  
  //! Construct a 3D Size metric
  /**
   * \param r_area The reference area for a 2D metric.
   */
  explicit Size_Metric_3(double r_vol) { ref_vol=r_vol; }

  //! Calculate the metric value on this element.
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;

private:
  double ref_vol;
};

//! 3D Size-Shape %Metric Class
class Size_Shape_Metric_3 : public Alg_Metric_Base_3 {
public:

  //! Construct a 3D Size-Shape metric
  /**
   * \param r_area The reference volume for this element.
   */
  explicit Size_Shape_Metric_3(double r_vol) { ref_vol=r_vol; }

  //! Calculate the metric value on this element.
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;

private:
  double ref_vol;
};

//! 3D Skew %Metric Class
class Skew_Metric_3 : public Alg_Metric_Base_3 {
public:

  //! Construc a 3D Skew Metric
  Skew_Metric_3() {}

  //! Calculate the metric value on this element.
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;
};

//! 3D Size-Skew %Metric Class
class Size_Skew_Metric_3 : public Alg_Metric_Base_3 {
public:

  //! Construct a 3D Size-Skew metric
  /**
   * \param r_area The reference volume for this element.
   */
  explicit Size_Skew_Metric_3(double r_vol) { ref_vol=r_vol; }

  //! Calculate the metric value on this element.
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;

private:
  double ref_vol;
};

MOP_END_NAMESPACE

#endif






