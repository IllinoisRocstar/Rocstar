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
// $Id: Algebraic_Metrics_2.h,v 1.5 2008/12/06 08:45:23 mtcampbe Exp $

/*! \file Algebraic_Metrics_2.h 
    \brief 2D algebraic quality %Metric declaration..
*/

#ifndef __ALGEBRAIC_METRICS_2_H__
#define __ALGEBRAIC_METRICS_2_H__

#include "Metric.h"

MOP_BEGIN_NAMESPACE

//! 2D Algebraic %Metric Base Class
/** 
 * This class is the base class for the implementation of
 * the specific 2D algebraic metrics defined by Patrick Knupp in
 * "Aglebraic Mesh Quality Metrics for Unstructured Initial Meshes"
 */ 
class Alg_Metric_Base_2 : public Metric {
  friend std::ostream & operator<<(std::ostream &, 
				   const Alg_Metric_Base_2 &);
public:
  
  //! Constructor
  Alg_Metric_Base_2() {}
  
  //! Virtual Destructor
  /**
   *  Although the destructor does nothing, this member function must be declared
   *  since the class has other virtual functions.
   */
  virtual ~Alg_Metric_Base_2() {}  

  //! Initialize a 2D Algebraic Metric
  /**
   * \param n[] Ordered set of nodal coordinates.
   * \param type Element type: TRI or QUAD
   */
  virtual void initialize(Vector_3<double> n[], int type);

  virtual void initialize(Element_node_enumerator &ene);

  //! The maximum value for this metric.
  virtual double maxValue() const { return 1.0; }

  //! The minimum value for this metric.
  virtual double minValue() const { return 0.0; }

protected:

  //! Compute the size metric.
  /**
   * \param ref_area the ideal area for this element.
   */
  double compute_size(double ref_area=1.) const;

  //! Compute the shape metric.
  double compute_shape() const;

  //! Compute the skew metric.
  double compute_skew() const;

protected:
  double   alpha[4];
  J_Matrix A[4];
  Matrix   L[4];
  int type_;
};

//! 2D Shape %Metric Class
class Shape_Metric_2 : public Alg_Metric_Base_2 {
public:

  //! Initialize a 2D Shape metric.
  Shape_Metric_2() {}

  //! Calculate the shape metric value
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;
};

//! 2D Size %Metric Class
class Size_Metric_2 : public Alg_Metric_Base_2 {
public:

  //! Construct a 2D Size metric
  /**
   * \param r_area The reference area for this element..
   */
  explicit Size_Metric_2(double r_area) { ref_area=r_area; }

  //! Calculate the shape metric value
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;

private:
  double ref_area;
};

//! 2D Size-Shape %Metric Class
class Size_Shape_Metric_2 : public Alg_Metric_Base_2 {
public:

  //! Construct a Size_Shape_Metric_2
  /**
   * \param r_area the reference area for this element.
   */
  explicit Size_Shape_Metric_2(double r_area) { ref_area=r_area; }

  //! Calculate the metric value
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;

  //! modify the reference area for this element.
  void set_ref_area(double r_area) { ref_area = r_area; }

private:
  double ref_area;
};

//! 2D Skew %Metric Class
class Skew_Metric_2 : public Alg_Metric_Base_2 {
public:

  // Construct a Skew_Metric_2
  Skew_Metric_2() {}

  //! Calculate the metric value
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;
};

//! 2D Size-Skew %Metric Class
class Size_Skew_Metric_2 : public Alg_Metric_Base_2 {
public:

  //! Construct Size_Skew_Metric_2
  /**
   * \param r_area the reference area for this element.
   */
  explicit Size_Skew_Metric_2(double r_area) { ref_area=r_area; }
 
  //! Calculate the metric value
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;

private:
  double ref_area;
};

MOP_END_NAMESPACE

#endif






