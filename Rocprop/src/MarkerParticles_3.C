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
// $Id: MarkerParticles_3.C,v 1.8 2008/12/06 08:45:27 mtcampbe Exp $

#include "MarkerParticles_3.h"
#include "Rocblas.h"

PROP_BEGIN_NAMESPACE

double
MarkerParticles_3::time_stepping( const COM::DataItem *spd, double dt,
				  COM::DataItem *disp, int *smoothed) {
  COM_assertion( spd->size_of_components()==1);
  COM_assertion( disp->size_of_components()==3 && disp->is_nodal());

  multiply_nodal_normals( spd, disp);
  Rocblas::mul_scalar( disp, &dt, disp);
  
  if ( smoothed) *smoothed = false;
  return dt;
}

// Multiply dataitem speed by nodal normals and save into attribute disp
void 
MarkerParticles_3::multiply_nodal_normals( const COM::DataItem *spd,
					   COM::DataItem *disp) {
  assert( spd->size_of_components()==1);

  if ( spd->is_nodal()) { 
    // compute normals
    _surf->compute_normals( disp);
    
    // If spd is not elemental, simply need call Rocblas
    Rocblas::mul( disp, spd, disp);
  }
  else {
    // Create a buffer window to store the element normals
    COM::Window *win = disp->window();
    COM::Window buf( win->name()+"-marker", win->get_communicator());
    buf.inherit( win->dataitem( COM::COM_MESH), "", false, true, NULL, 0);
    DataItem *node_disps = buf.inherit( disp, "disps", false, true, NULL, 0);
    buf.inherit(const_cast<COM::DataItem*>(spd), "spd", false, true, NULL, 0);
    DataItem *elem_normals = 
      buf.new_dataitem( "elem_normals", 'e', COM_DOUBLE, 3, "");
    buf.resize_array( "elem_normals", 0, NULL);
    buf.init_done();

    // First, compute elemental normals
    _surf->compute_normals( elem_normals);

    // Multiply elemental normals by speed
    Rocblas::mul( elem_normals, spd, elem_normals);

    // Convert elemental motion to nodal motion
    _surf->elements_to_nodes( elem_normals, node_disps, SURF::E2N_ANGLE);
  }
}


PROP_END_NAMESPACE






