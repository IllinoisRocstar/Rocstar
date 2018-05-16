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
// $Id: PointPropagate.C,v 1.4 2008/12/06 08:45:28 mtcampbe Exp $

#include "PointPropagate.h"
#include "Generic_element_2.h"

PROP_BEGIN_NAMESPACE

// Perform time integration at a specific point.
Vector_3 PointPropagate::
time_integrate( const Point_3 &p, const Speed &s, 
		double t, double dt, int order) {
  Vector_3 disp1 = dt*s.get_velocity( p, t);

  // If first order, using Euler's method
  if ( order <= 1) { return disp1; }

  // Otherwise, use fourth-order Runge-Kutta.
  Vector_3 disp2 = dt*s.get_velocity( p+0.5*disp1, t+0.5*dt);
  Vector_3 disp3 = dt*s.get_velocity( p+0.5*disp2, t+0.5*dt);
  Vector_3 disp4 = dt*s.get_velocity( p+disp3, t+dt);
  
  return (disp1+disp4)/6+(disp2+disp3)/3;
}

// Propagate the Gauss points of a given element.
void PointPropagate::
propagate_gp_of_element(const Point_3 *pnts, COM::Element_node_enumerator &ene,
			const Speed &s, double t, double dt, Vector_3 *disps) {

  // Only support triangular elements now
  COM_assertion( ene.size_of_edges()==3);

  // Evaluate current face normal 
  COM::Element_node_vectors_k_const<Point_3> ps; ps.set( pnts, ene, 1);
  SURF::Generic_element_2 e(ene.size_of_edges(), ene.size_of_nodes());

  // Look through the Gauss points and compute their displacements.
  for ( int i=0; i<3; ++i) {
    // Obtain the positions of the Gauss point
    Vector_2 nc; e.get_gp_nat_coor( i, nc, 2);
    Point_3  p;  e.interpolate( ps, nc, &p);

    // Obtain the displacement at the Gauss point using 4th order RK
    disps[ i] = time_integrate( p, s, t, dt, 4);
  }
}

// Propagate the Gauss points of all faces of a given mesh
void PointPropagate::
propagate_faces( const std::string &wname, const Speed &s, 
		 double t, double dt, const std::string &vname) {
  COM::Window *win = COM_get_com()->get_window_object(wname);
  COM::DataItem *velos = win->dataitem( vname);

  COM_assertion_msg( velos->is_elemental() && velos->size_of_components()==9,
		     "Velos must be elemental data associated with quadrature points");

  // Loop through all the panes and then all the faces
  std::vector<COM::Pane*> panes; win->panes( panes);

  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = panes.begin();
  for ( int i=0, local_npanes = panes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;
    const Point_3 *pnts = reinterpret_cast<Point_3*>(pane->coordinates());
    Vector_3 *velos_ptr = reinterpret_cast<Vector_3*>
      (pane->dataitem( velos->id())->pointer());
    
    // Loop through real elements of the current pane
    COM::Element_node_enumerator ene( pane, 1); 
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      // First compute the displacements at Gauss points using high-order
      // time integration
      propagate_gp_of_element( pnts, ene, s, t, dt, &velos_ptr[3*j]);

      // Convert displacements into velocities by dividing by time step.
      velos_ptr[3*j] /= dt; velos_ptr[3*j+1] /= dt; velos_ptr[3*j+2] /= dt;
    }
  }
}

// Propagate the nodes of a given mesh
void PointPropagate::
propagate_nodes( const std::string &wname, const Speed &s, 
		 double t, double dt, const std::string &vname) {
  COM::Window *win = COM_get_com()->get_window_object(wname);
  COM::DataItem *velos = win->dataitem( vname);

  COM_assertion_msg( velos->is_nodal() && velos->size_of_components()==3,
		     "Velos must be nodal data associated");

  double a = 1./dt;
  
  // Loop through all the panes and then all the vertices
  std::vector<COM::Pane*> panes; win->panes( panes);

  std::vector< COM::Pane*>::iterator it = panes.begin();
  for ( int i=0, local_npanes = panes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;
    const Point_3 *pnts = reinterpret_cast<Point_3*>(pane->coordinates());
    Vector_3 *velos_ptr = reinterpret_cast<Vector_3*>
      (pane->dataitem( velos->id())->pointer());
    
    // Loop through real vertices of the current pane
    for ( int j=0, nj=pane->size_of_real_nodes(); j<nj; ++j) {
      velos_ptr[j] = a*time_integrate( pnts[j], s, t, dt, 4);
    }
  }
}

PROP_END_NAMESPACE






