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
// $Id: PointPropagate.h,v 1.4 2008/12/06 08:45:28 mtcampbe Exp $

#ifndef _POINT_PROPAGATE_H_
#define _POINT_PROPAGATE_H_

#include "speeds.h"
#include "Element_accessors.hpp"

PROP_BEGIN_NAMESPACE

class PointPropagate {
public:
  // Propagate the faces of a given window.
  static void propagate_faces( const std::string &wname, const Speed &s,
			       double t, double dt, const std::string &vname);

  // Propagate the nodes of a given window.
  static void propagate_nodes( const std::string &wname, const Speed &s,
			       double t, double dt, const std::string &vname);

protected:
  // Integrate the motion of a given point
  static Vector_3 time_integrate( const Point_3 &p, const Speed &s, 
				  double t, double dt, int order);
  
  // Integrate the motion of the Gauss points of a given face
  static void propagate_gp_of_element(const Point_3 *pnts, 
				      COM::Element_node_enumerator &ene, 
				      const Speed &s, double t, double dt, 
				      Vector_3 *qp_vels);

};

PROP_END_NAMESPACE

#endif






