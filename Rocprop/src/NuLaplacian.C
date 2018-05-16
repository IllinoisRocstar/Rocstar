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
// $Id: NuLaplacian.C,v 1.8 2008/12/06 08:45:27 mtcampbe Exp $

#include "FaceOffset_3.h"
#include "Rocblas.h"
#include "Generic_element_2.h"
#include "Rocsurf.h"

PROP_BEGIN_NAMESPACE

// At input, vcenters contains displacements due to Laplacian smoothing,
// disps is used as buffer space during smoothing.
// At output, vcenters contains the tangential motion.
void FaceOffset_3::
nulaplacian_smooth( const COM::DataItem *vert_normals,
		    const COM::DataItem *tangranks,
		    const COM::DataItem *disps,
		    COM::DataItem *vcenters,
		    COM::DataItem *buf,
		    COM::DataItem *vert_weights_buf,
		    const COM::DataItem *edge_weights) {
  const double alpha = 0.03;

  // First, sum up weights for each vertex.
  double eps = 1.e-100;
  Rocblas::copy_scalar( &eps, vert_weights_buf);

  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  for ( int i=0, local_npanes = _panes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it; 
    
    // Obtain address for buffer space for vertex weights
    double *vws = reinterpret_cast<double*>
      (pane->dataitem(vert_weights_buf->id())->pointer());
    // Obtain address for edge weights
    const double *ews = edge_weights ? reinterpret_cast<const double*>
      (pane->dataitem(edge_weights->id())->pointer()) : NULL;

    Element_node_enumerator ene( pane, 1); 
    // Loop through the real elements of the pane
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne = ene.size_of_edges();     

      COM_assertion_msg( ne==4, "NuLapacian supports only quadrilateral meshes");

      // Loop through the vertices of the element to compute weights
      for ( int k=0; k<ne; ++k) {
	vws[ene[k]-1] += ews ? ews[j*4+k] + ews[j*4+(k+3)%4] : 1;
      }
    } 
  } 
  _surf->reduce_on_shared_nodes( vert_weights_buf, Manifold::OP_SUM);

  // Zero out buf
  double zero=0;
  Rocblas::copy_scalar( &zero, buf);

  // Loop through the panes and its real faces
  it = _panes.begin(); 
  for ( int i=0, local_npanes = _panes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it; 
    
    // Obtain nodal coordinates of current pane, assuming contiguous layout
    const Vector_3 *pnts = reinterpret_cast<Vector_3*>
      (pane->dataitem(COM_NC)->pointer());
    const Vector_3 *ds = reinterpret_cast<Vector_3*>
      (pane->dataitem(disps->id())->pointer());
    const Vector_3 *nrms = reinterpret_cast<Vector_3*>
      (pane->dataitem(vert_normals->id())->pointer());
    const Vector_3 *vcnt = reinterpret_cast<Vector_3*>
      (pane->dataitem(vcenters->id())->pointer());
    
    // Obtain address for output displacements
    Vector_3 *dbuf = reinterpret_cast<Vector_3*>
      (pane->dataitem(buf->id())->pointer());
    // Obtain address for buffer space for vertex weights
    const double *vws = reinterpret_cast<const double*>
      (pane->dataitem(vert_weights_buf->id())->pointer());
    // Obtain address for edge weights
    const double *ews = edge_weights ? reinterpret_cast<const double*>
      (pane->dataitem(edge_weights->id())->pointer()) : NULL;

    // pnts_buf stores difference between center and perturbed position.
    Vector_3 pnts_buf[4];
    int      is_regular[4];

    Element_node_enumerator ene( pane, 1); 
    // Loop through the real elements of the pane
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne = ene.size_of_edges();     

      COM_assertion_msg( ne==4, 
			 "NuLapacian now supports only quadrilateral meshes");

      // Compute face center of the face perturbed along normal direction.
      Vector_3 cnt(0,0,0);
      for ( int k=0; k<ne; ++k) {
	int uindex = ene[k]-1;

	is_regular[k] = std::abs(vws[uindex]-4)<1.e-6;
	pnts_buf[k] = pnts[uindex] + ds[uindex];

	cnt += pnts_buf[k];
      }
      cnt *= 0.25; // Divide cnt by ne

      // pnts_buf contains the displacements from the face center to the vertex
      for ( int k=0; k<ne; ++k) pnts_buf[k] -= cnt;
      
      // Loop through the edges of the element
      int uindex=ene[ne-1]-1, vindex=ene[0]-1;
      for ( int k=0; k<ne; ++k, uindex=vindex, vindex=ene[(k==ne)?0:k]-1) {
	// Compute normal of the perturbed edge as the cross product 
	// of the vertex normal and edge tangent
	Vector_3 tng = pnts[uindex]-pnts[vindex]+ vcnt[uindex]-vcnt[vindex];
	double w = ews ? ews[j*4+k] : 0.5;

	int k_1 = (k+3)%4;
	if ( is_regular[k_1]) {
	  // Accumulate buf and vert_weights for vertex uindex
	  Vector_3 nrm_e=Vector_3::cross_product(tng, nrms[uindex]).normalize();
	  
	  dbuf[uindex] -= w * ( pnts_buf[k_1] * nrm_e * nrm_e +
				alpha * pnts_buf[k_1]);
	}
	
	if ( is_regular[k]) {
	  // Accumulate buf and vert_weights for vertex vindex
	  Vector_3 nrm_e=Vector_3::cross_product(tng, nrms[vindex]).normalize();

	  dbuf[vindex] -= w * ( pnts_buf[k] * nrm_e * nrm_e +
				alpha * pnts_buf[k]);
	}
      }
    } 
  }

  // Reduce on shared nodes for buf and vert_weights_buf
  _surf->reduce_on_shared_nodes( buf, Manifold::OP_SUM);

  // Loop through all real vertices to update displacement vector
  it = _panes.begin();
  for (int i=0, local_npanes = _panes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;

    const char *tranks = reinterpret_cast<const char*>
      ( pane->dataitem( tangranks->id())->pointer());
    const double *vws = reinterpret_cast<const double*>
      ( pane->dataitem( vert_weights_buf->id())->pointer());

    Vector_3 *vcnt = reinterpret_cast<Vector_3*>
      ( pane->dataitem( vcenters->id())->pointer());
    Vector_3 *dbuf = reinterpret_cast<Vector_3*>
      ( pane->dataitem( buf->id())->pointer());

    const Vector_3 *es = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_eigvecs->id())->pointer());

    // Loop through all real nodes of the pane
    for ( int j=0, nj=pane->size_of_real_nodes(); j<nj; ++j) {
      // If vertex is unmasked and is regular
      if ( std::abs(vws[j]-4) < 1.e-6) {
	// Assign the displacement to be the weighted average of projections
	vcnt[j] = dbuf[j] / vws[j];
	
	if ( tranks[j] == 1) // Project onto the direction
	  vcnt[j] = (vcnt[j]*es[3*j+2])*es[3*j+2];
	else if ( tranks[j] == 2) // Project onto the plane
	  vcnt[j] -= (vcnt[j]*es[3*j])*es[3*j];
	else
	  vcnt[j] = Vector_3(0,0,0);
      }
    }
  }
}

PROP_END_NAMESPACE






