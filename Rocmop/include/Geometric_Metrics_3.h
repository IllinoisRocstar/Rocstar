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
// $Id: Geometric_Metrics_3.h,v 1.7 2008/12/06 08:45:23 mtcampbe Exp $

/*! \file Geometric_Metrics_3.h 
    \brief 3D geometric quality %Metric declarations.
*/

#ifndef __GEOMETRIC_METRICS_3_H__
#define __GEOMETRIC_METRICS_3_H__

#include "Metric.h"

#include "Pane.hpp"

MOP_BEGIN_NAMESPACE

//! 3D Geometric %Metric Base Class
/** 
 * This class is the base class for the implementation of
 * various 3D geometric element quality metrics.
 */ 
class Geo_Metric_Base_3 : public Metric {

public:

  // Constructor
  Geo_Metric_Base_3() : _ene_pane(NULL), _ene_i(0){;}

  virtual ~Geo_Metric_Base_3(){;}

  //! Initialize a 3D Geometric Metric by nodal coords and element type
  /**
   * \param n[] Ordered set of nodal coordinates.
   * \param type Element type: TET or HEX
   */
  virtual void initialize(Vector_3<double> n[], int type);
  
  //! Initialize a 3D Geometric Metric by Element_node_enumerator
  /**
   * \param ene Element_node_enumerator
   */
  virtual void initialize(Element_node_enumerator &ene);

  //! The maximum value for this metric.
  virtual double maxValue() const { return 1.0; }

  //! The minimum value for this metric.
  virtual double minValue() const { return 0.0; }

protected:

  //! Compute min and max dihedral angles
  void compute_angles(double& min, double& max) const;

  //! Compute circumradius, inradius, and shortest edge length
  void compute_aspects(double& R, double& r, double& l) const;
  
protected:
  std::vector<Vector_3<double> > v;
  int _type;
  const COM::Pane * _ene_pane;
  int _ene_i;
};

//! 3D Max and Min Angle %Metric Class
class Angle_Metric_3 : public Geo_Metric_Base_3 {
public:

  //! Initialize a 3D Angle Metric.
  Angle_Metric_3() {}

  //! Calculate the metric value on this element.
  /*!  
   *  \param atts[] Computed value(s) placed here.
   */
  virtual void compute(double atts[]) const;

  //! Get the maximum value for this metric.
  double maxValue() const;

  //! Get the minimum value for this metric
  double minValue() const;
};

//! 3D Aspect Ratios %Metric Class
class Aspect_Metric_3 : public Geo_Metric_Base_3 {
public:

  //! Initialize a 3D Aspect Metric.
  Aspect_Metric_3() {}

  //! Calculate the metric value on this element.
  /*!  
   *  \param atts[] Computed value(s) placed here.
   *
   * Note that R/r and R/l are only implemented for simplical elements.
   * This method will return -1 for both metrics if used on a 
   * hexahedron.
   */  
  virtual void compute(double atts[]) const;

  //! Get the maximum value for this metric.
  double maxValue() const;

  //! Get the minimum value for this metric
  double minValue() const;

  //! Get the geometric aspects
  void getAspects(double& R, double& r, double& l){
    compute_aspects(R,r,l);
  }
};

MOP_END_NAMESPACE

#endif






