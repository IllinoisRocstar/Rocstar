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
// $Id: AngleBasedSmoothing.C,v 1.4 2008/12/06 08:45:27 mtcampbe Exp $

#include "FaceOffset_3.h"
#include "../Rocblas/include/Rocblas.h"
#include "../Rocsurf/include/Generic_element_2.h"
#include "../Rocsurf/include/Rocsurf.h"

PROP_BEGIN_NAMESPACE

void FaceOffset_3::
compute_angle_based_vertex_centers() {
  COM_assertion_msg( COMMPI_Comm_size( _buf->get_communicator())==1, 
		     "Angle-based smoothing supports only serial runs now");

  double zero = 0., eps = 1.e-100;
  Rocblas::copy_scalar( &zero, _vcenters);
  Rocblas::copy_scalar( &eps, _weights);

  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  for ( int i=0, local_npanes = _panes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    
    const Point_3 *pnts = reinterpret_cast<Point_3*>
      (pane->dataitem(COM_NC)->pointer());
    const char *tranks = reinterpret_cast<const char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    const Vector_3 *es = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_eigvecs->id())->pointer());
    const Vector_3 *vnrm = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_vnormals->id())->pointer());

    Vector_3 *vcnts = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_vcenters->id())->pointer());
    double   *ws = reinterpret_cast<double*>
      ( pane->dataitem(_weights->id())->pointer());

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1);
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne = ene.size_of_edges();

      // Loop through all Halfedge of current face.
      for ( int k=0; k<ne; ++k) {
	Halfedge h( &*pm_it, Edge_ID( j+1,k), SURF::REAL_PANE);
	int vindex=ene[k]-1; // Destination of the halfedge 
	
	Vector_3 dirs[2] = { es[3*vindex+1], es[3*vindex+2] };
	// Skip corners.
	switch ( tranks[vindex]) {
	case 0:
	  continue;
	case 1:	  
	  // Special handling for ridges.
	  dirs[0] = Vector_3::cross_product(dirs[1], vnrm[vindex]);
	default: {
	  Point_3  pv = h.origin().point();
	  Vector_3 ds[3] = { h.opposite().prev().origin().point()-pv,
			     h.destination().point() - pv,
			     h.next().destination().point() - pv};
	
	  Vector_2 p0 = proj( ds[0], dirs[0], dirs[1]);
	  Vector_2 p1 = proj( ds[1], dirs[0], dirs[1]);
	  Vector_2 p2 = proj( ds[2], dirs[0], dirs[1]);
	  
	  Vector_2 d1 = (-p1).normalize();
	  Vector_2 d2 = (p2-p1).normalize();
	  Vector_2 d0 = (p0-p1).normalize();

	  double alpha1 = std::acos( d1*d2);
	  double alpha2 = std::acos( d0*d1);
	  double w = 1/(alpha1+alpha2)/(alpha1+alpha2);

	  ws[vindex] += w;

	  double beta = 0.5*(alpha2-alpha1);
	  double cosb=std::cos(beta), sinb=std::sin(beta);
	  Vector_2 pnt( p1[0]-p1[0]*cosb+p1[1]*sinb, 
			p1[1]-p1[0]*sinb-p1[1]*cosb);
	  
	  vcnts[vindex] += w*(pnt[0]*dirs[0]+pnt[1]*dirs[1]);
	}
	}
      }
    }
  }

  _surf->reduce_on_shared_nodes( _vcenters, Manifold::OP_SUM);
  _surf->reduce_on_shared_nodes( _weights, Manifold::OP_SUM);
  Rocblas::div( _vcenters, _weights, _vcenters);
}

PROP_END_NAMESPACE






