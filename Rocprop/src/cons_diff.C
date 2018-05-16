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
// $Id: cons_diff.C,v 1.4 2008/12/06 08:45:28 mtcampbe Exp $

#include "FaceOffset_3.h"
#include "Rocblas.h"
#include "Generic_element_2.h"
#include "Rocsurf.h"

PROP_BEGIN_NAMESPACE

void FaceOffset_3::
balance_mass() {
  COM::DataItem *disps = 
    _buf->new_dataitem( "disps_buf_bm", 'n', COM_DOUBLE, 3, "");
  _buf->resize_array( disps, 0);
  COM::DataItem *pos = 
    _buf->new_dataitem( "pos_buf_bm", 'n', COM_DOUBLE, 3, "");
  _buf->resize_array( pos, 0);
  COM::DataItem *fvol_buf = 
    _buf->new_dataitem( "vfol_buf_bm", 'e', COM_DOUBLE, 1, "");
  _buf->resize_array( fvol_buf, 0);
  _buf->init_done( false);

  // Use face-area as buffer for facial-volume
  COM::DataItem *fvol = _faceareas;

  // Compute volume change at each vertex.
  const COM::DataItem *coors = _buf->dataitem( COM::COM_NC);
  Rocblas::add( coors, _vcenters, pos);

  // Obtain density-vector for vertices.
  int normalize=false;
  SURF::Rocsurf::compute_element_normals( _facenormals, &normalize, pos);
  _surf->elements_to_nodes( _facenormals, disps, SURF::E2N_ONE, NULL, NULL, 1);
  
  // Use _scales for vortex area and _weights for vertex-volume.
  Rocblas::dot( disps, _vnormals, _scales);
  double zero=0;
  Rocblas::copy_scalar( &zero, disps);

  for ( int k=0;k<2;++k) {
    if ( k) {
      SURF::Rocsurf::compute_swept_volumes( pos, disps, fvol);

      Rocblas::sub( fvol_buf, fvol, fvol);
    }
    else {
      SURF::Rocsurf::compute_bounded_volumes( pos, coors, fvol);
      Rocblas::copy( fvol, fvol_buf);
    }
    distribute_volume_e2n( fvol, _tangranks, _weights);


    std::vector< COM::Pane*>::iterator it = _panes.begin();
    // Loop through the panes and its real vertices
    for (int i=0, local_npanes = _panes.size(); i<local_npanes; ++i, ++it) { 
      COM::Pane *pane = *it;
      
      const double *vs = reinterpret_cast<const double*>
	( pane->dataitem(_weights->id())->pointer());
      const double *as = reinterpret_cast<const double*>
	( pane->dataitem(_scales->id())->pointer());
      const Vector_3 *ds_m = reinterpret_cast<const Vector_3*>
	( pane->dataitem(_vnormals->id())->pointer());
      Vector_3 *ds = reinterpret_cast<Vector_3*>
	( pane->dataitem(disps->id())->pointer());

      // Loop through all real nodes of the pane
      for ( int j=0, jn=pane->size_of_real_nodes(); j<jn; ++j) {
	ds[j] += vs[j]/as[j]*ds_m[j]; // Omit the factor 3 here
      }
    }
  }

  Rocblas::add( _vcenters, disps, _vcenters);
  
  _buf->delete_dataitem( fvol_buf->name());
  _buf->delete_dataitem( pos->name());
  _buf->delete_dataitem( disps->name());
  _buf->init_done( false);
}

void FaceOffset_3::
distribute_volume_e2n( const COM::DataItem *fvol,
		       const COM::DataItem *tranks,
		       COM::DataItem *vvol) {

  double zero=0;
  Rocblas::copy_scalar( &zero, vvol);

  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  for ( int i=0, local_npanes = _panes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;

    const double *fv = reinterpret_cast<const double*>
      ( pane->dataitem(fvol->id())->pointer());
    double   *vv = reinterpret_cast<double*>
      ( pane->dataitem(vvol->id())->pointer());

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1); 
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne = ene.size_of_edges();

      for (int k=0; k<ne; ++k) {
	int index = ene[k]-1;
	vv[index] += fv[j]; // Omit the factor 1/3 here
      }
    }
  }
  
  // Reduce on shared nodes
  _surf->reduce_on_shared_nodes( vvol, Manifold::OP_SUM);
 
}

PROP_END_NAMESPACE






