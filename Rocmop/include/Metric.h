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
// $Id: Metric.h,v 1.6 2008/12/06 08:45:23 mtcampbe Exp $


/** \file Metric.h 
 *
 * The %Metric Base Class, from which algebraic and
 * geometric metric clases are derived.
 *
*/

#ifndef __METRIC_H__
#define __METRIC_H__

#include "mopbasic.h"
#include "Matrix.h"
#include "Connectivity.hpp"
#include "Element_accessors.hpp"

MOP_BEGIN_NAMESPACE
using COM::Element_node_enumerator;
using COM::Element_node_vectors_k_const;

//! Metric Base Class
/*! 
 */ 
class Metric {
public:

  Metric() {}

  virtual ~Metric(){;}

  virtual void initialize(Vector_3<double> n[], int type)=0;
  virtual void initialize(Element_node_enumerator &ene)=0;
  
  //! Return maximum value for this metric.
  virtual double maxValue() const =0;

  //! Return minimum value for this metric.
  virtual double minValue() const =0;

  //! Compute the value of this metric.
  virtual void compute(double atts[]) const=0;
};

MOP_END_NAMESPACE

#endif






