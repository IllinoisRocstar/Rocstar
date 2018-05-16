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
// $Id: Geometric_Metrics_2.h,v 1.5 2008/12/06 08:45:23 mtcampbe Exp $

/*! \file Geometric_Metrics_2.h 
    \brief 2D geometric quality %Metric declarations.
*/

#ifndef __GEOMETRIC_METRICS_2_H__
#define __GEOMETRIC_METRICS_2_H__

#include "Metric.h"

MOP_BEGIN_NAMESPACE

//! 2D Geometric %Metric Base Class
/** 
 * This class is the base class for the implementation of
 * various 2D geometric element quality metrics.
 */ 
class Geo_Metric_Base_2 : public Metric {

public:
  
  //! Constructor
  Geo_Metric_Base_2() {}
  
  //! Initialize a 2D Geometric Metric
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

  //! Compute the min and max angles.
  void compute_angles(double& min, double& max) const;

  //! Compute the circumradius, inradius, and shortest edge length.
  void compute_aspects(double& R, double& r, double& l) const;

protected:
  Vector_3<double> v[4];
  int type_;
};

//! 2D Max and Min Angle %Metric Class
class Angle_Metric_2 : public Geo_Metric_Base_2 {
public:

  //! Initialize a 2D Angle Metric.
  Angle_Metric_2() {}

  //! Calculate max and min angles
  /**  
   * \param atts[] Two computed values placed here.
   */
  virtual void compute(double atts[]) const;

  //! The maximum value for this metric
  double maxValue() const;

  //! The minimum value for this metric
  double minValue() const;
};

//! 2D Aspect Ratios %Metric Class
class Aspect_Metric_2 : public Geo_Metric_Base_2 {
public:

  //! Initialize an Aspect Ratio Metric.
  Aspect_Metric_2() {}

  //! Calculate scaled R/r and R/l.
  /**  
   * \param atts[] R/r placed in atts[0], R/l in atts[1]
   *
   * Note that R/r and R/l are only implemented for simplical elements.
   * This method will return -1 for both metrics if used on a 
   * quadrilateral.
   */
  virtual void compute(double atts[]) const;

  //! The maximum value for this metric.
  double maxValue() const;

  //! The minimum value for this metric.
  double minValue() const;
};

MOP_END_NAMESPACE

#endif






