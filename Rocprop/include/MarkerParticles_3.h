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
// $Id: MarkerParticles_3.h,v 1.7 2008/12/06 08:45:27 mtcampbe Exp $

/*******************************************************************
 * This file implements the marker particle method for 
 * interface propagation. It inherits some basic primitives
 * from the base class Propagation_3.
 *******************************************************************/

#ifndef __MARKER_PARTICLES_3_H_
#define __MARKER_PARTICLES_3_H_

#include "Propagation_3.h"

PROP_BEGIN_NAMESPACE

class MarkerParticles_3 : public Propagation_3 {
public:
  /// Construct an object from a window manifold
  MarkerParticles_3( Manifold *wm, COM::Window *buf) 
    : Propagation_3( wm, buf) {}

  /// Main entry of the algorithm
  virtual double time_stepping( const COM::DataItem *spd, double dt,
				COM::DataItem *disp, int *smoothed=NULL);

protected:
  // Multiply attribute a by nodal normals and save into attribute b
  void multiply_nodal_normals( const COM::DataItem *a,
			       COM::DataItem *b);
};

PROP_END_NAMESPACE

#endif






