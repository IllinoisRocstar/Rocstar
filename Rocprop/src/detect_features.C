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
// $Id: detect_features.C,v 1.7 2008/12/06 08:45:28 mtcampbe Exp $

#include "FaceOffset_3.h"
#include "Rocblas.h"
#include "Rocsurf.h"
#include <algorithm>

PROP_BEGIN_NAMESPACE

void FaceOffset_3::
filter_and_identify_ridge_edges( bool filter_curves) 
{
  // Calculate boundary normals
  _surf->update_bd_normals( _facenormals, false);

  COM::DataItem *maxdianglev = 
    _buf->new_dataitem( "FO_maxdianglev", 'n', COM_DOUBLE, 1, "");
  _buf->resize_array( maxdianglev, 0);

  double t=-pi;
  Rocblas::copy_scalar( &t, maxdianglev);

  COM::DataItem *mindianglev = 
    _buf->new_dataitem( "FO_mindianglev", 'n', COM_DOUBLE, 1, "");
  _buf->resize_array( mindianglev, 0);
  t=pi;
  Rocblas::copy_scalar( &t, mindianglev);

  COM::DataItem *maxtrnangv = 
    _buf->new_dataitem( "FO_maxtrnangv", 'n', COM_DOUBLE, 1, "");
  _buf->resize_array( maxtrnangv, 0);

  t=-2;
  Rocblas::copy_scalar( &t, maxtrnangv);

  COM::DataItem *mintrnangv = 
    _buf->new_dataitem( "FO_mintrnangv", 'n', COM_DOUBLE, 1, "");
  _buf->resize_array( mintrnangv, 0);
  t=2;
  Rocblas::copy_scalar( &t, mintrnangv);

  COM::DataItem *neighbor_features = 
    _buf->new_dataitem( "FO_neighbor_features", 'n', COM_CHAR, 1, "");
  _buf->resize_array( neighbor_features, 0);

  COM::DataItem *qsedges = 
    _buf->new_dataitem( "FO_qsedges", 'e', COM_CHAR, 1, "");
  _buf->resize_array( qsedges, 0);
  _buf->init_done( false);

  int zero=0;
  Rocblas::copy_scalar( &zero, _strong);

  // Clear ridge edges
  _edges.clear(); _edges.resize( _panes.size());

  bool to_filter = filter_curves && _panes.size() == 1 && !COMMPI_Initialized();

  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();

  int nedges=0, nvertices=0;
  std::vector<std::map<Edge_ID, double> > edge_maps( _panes.size());
  std::vector<std::map<Edge_ID, double> > diangle_maps( _panes.size());

  for ( int i=0, local_npanes = _panes.size();
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    std::map< Edge_ID, double> &emap_pn = edge_maps[i];
    std::map< Edge_ID, double> &diangle_pn = diangle_maps[i];

    nvertices += pane->size_of_real_nodes();

    // Obtain nodal coordinates of current pane, assuming contiguous layout
    const Point_3 *pnts = reinterpret_cast<Point_3*>
      (pane->dataitem(COM_NC)->pointer());

    const Vector_3 *es = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_eigvecs->id())->pointer());
    char *tranks = reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    char *strong = reinterpret_cast<char*>
      ( pane->dataitem(_strong->id())->pointer());
    double *minvs = reinterpret_cast<double*>
      ( pane->dataitem(mindianglev->id())->pointer());
    double *maxvs = reinterpret_cast<double*>
      ( pane->dataitem(maxdianglev->id())->pointer());
    double *minta = reinterpret_cast<double*>
      ( pane->dataitem(mintrnangv->id())->pointer());
    double *maxta = reinterpret_cast<double*>
      ( pane->dataitem(maxtrnangv->id())->pointer());
    Vector_3 *fnrms = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_facenormals->id())->pointer());

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1);
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne=ene.size_of_edges(), uindex=ene[ne-1]-1, vindex=ene[0]-1, windex=0;
      nedges += ne;
	
      for ( int k=0; k<ne; ++k, uindex=vindex, vindex=windex) {
	windex = ene[ (k+1)%ne]-1;

	Edge_ID eid(j+1,k), eid_opp = (*pm_it)->get_opposite_real_edge( eid);

	// Consider border edges as flat and skip them when identifying ridges
	if ( (*pm_it)->is_physical_border_edge( eid_opp)) 
	{ ++nedges; continue; }

	Vector_3 tng = pnts[ windex] - pnts[vindex];
	tng.normalize();
	double feg = es[3*vindex+2]* tng;

	// Safeguard strong angles.
	// Important: We assume this computation gives identical results 
	// for eid and eid_opp near the thresholds.
	const Vector_3 &n1 = fnrms[j];
	const Vector_3 &n2 = eid_opp.is_border() ? 
	  (*pm_it)->get_bd_normal( eid_opp) : fnrms[eid_opp.eid()-1];
	double costheta = n1*n2;
	if ( costheta>1) costheta = 1;
	else if (costheta<-1) costheta = -1;

	double fangle = std::acos( costheta);
	if ( fangle >  _tol_angle_strong) {
	  // Upgrade vertices.
	  if ( tranks[vindex]==2) tranks[vindex] = 1;
	  if ( tranks[windex]==2) tranks[windex] = 1;

	  if ( fangle>pi*0.505 && fangle !=pi)
	    strong[vindex] = strong[windex] = 2;
	  else
	    strong[vindex] = strong[windex] = 1;

	  std::map< Edge_ID, double>::iterator rit = diangle_pn.find(eid);
	  // Make sure the fangle is used consistently.
	  if ( rit == diangle_pn.end()) diangle_pn[eid] = fangle;
	  else fangle = rit->second;
 
	  if ( feg>0) {
	    if ( fangle > maxvs[vindex]) maxvs[vindex] = fangle; 
	    if ( feg > maxta[vindex]) maxta[vindex] = feg;
	    emap_pn[ eid] = feg;
	  }
	  else {
	    if ( -fangle<minvs[vindex]) minvs[vindex] = -fangle;   
	    if ( feg < minta[vindex]) minta[vindex] = feg;

	    emap_pn[ eid] = feg;
	  }
	}
	else if ( fangle > _tol_angle_weak) {
	  if ( feg > 0) {
	    if ( fangle > maxvs[vindex]) maxvs[vindex] = fangle;
	    if ( feg > maxta[vindex])    maxta[vindex] = feg;
	  }
	  else {
	    if ( -fangle < minvs[vindex]) minvs[vindex] = -fangle;
	    if ( feg < minta[vindex])    minta[vindex] = feg;
	  }

 	  // Save the current edge for comparison
	  diangle_pn[eid] = fangle; emap_pn[ eid] = feg;
	}
      } // for k (edges)
    } // for j (elements)
  } // for i (panes)

  if ( to_filter)
    std::cout << "There are " << nedges/2 << " edges and " 
	      << nvertices << " vertices " << std::endl;
  
  // Reduce on shared nodes for mindianglev and maxdianglev
  _surf->reduce_on_shared_nodes( mindianglev, Manifold::OP_MIN);
  _surf->reduce_on_shared_nodes( maxdianglev, Manifold::OP_MAX);

  _surf->reduce_on_shared_nodes( mintrnangv, Manifold::OP_MIN);
  _surf->reduce_on_shared_nodes( maxtrnangv, Manifold::OP_MAX);
  _surf->reduce_on_shared_nodes( _strong, Manifold::OP_MAX);

  // Count the number of incident strongly attached edges.
  it = _panes.begin();
  pm_it=_surf->pm_begin();
  Rocblas::copy_scalar( &zero, neighbor_features);

  for ( int i=0, local_npanes = _panes.size(); i<local_npanes; 
	++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    std::map< Edge_ID, double> &emap_pn = edge_maps[i];
    std::map< Edge_ID, double> &diangle_pn = diangle_maps[i];

    const double *minvs = reinterpret_cast<double*>
      ( pane->dataitem(mindianglev->id())->pointer());
    const double *maxvs = reinterpret_cast<double*>
      ( pane->dataitem(maxdianglev->id())->pointer());
    const double *minta = reinterpret_cast<double*>
      ( pane->dataitem(mintrnangv->id())->pointer());
    double *maxta = reinterpret_cast<double*>
      ( pane->dataitem(maxtrnangv->id())->pointer());

    char *tranks = reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    char *nnf = reinterpret_cast<char*>
      ( pane->dataitem(neighbor_features->id())->pointer());

    COM_assertion( emap_pn.size() == diangle_pn.size());
    // Loop through ebuf and put ridges edges onto emap_pn
    for ( std::map< Edge_ID, double>::iterator eit=emap_pn.begin(), 
	    rit=diangle_pn.begin(), eend=emap_pn.end(); eit!=eend; ++eit, ++rit) {
      Edge_ID eid = eit->first;
      double feg = eit->second, fangle = rit->second;

      COM_assertion(! eid.is_border());
      Element_node_enumerator ene( pane, eid.eid()); 
      int vi = eid.lid(), vindex = ene[vi]-1;

      // Determine the number of strongly attached edges
      if ( !tranks[vindex] || tranks[vindex] == char(-1) || fangle > _tol_angle_strong ||
	   ( feg < 0 && fangle == -minvs[vindex] || 
	     feg > 0 && fangle == maxvs[vindex]) &&
	   ( feg < -_tol_cos_osta && feg == minta[vindex] || 
	     feg > _tol_cos_osta && feg == maxta[vindex])) {
	++nnf[vindex];
      }
    } // for eit
  } // for i
   
  _surf->reduce_on_shared_nodes( neighbor_features, Manifold::OP_SUM);

  it = _panes.begin();
  pm_it=_surf->pm_begin();

  for ( int i=0, local_npanes = _panes.size(); i<local_npanes; 
	++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    std::map< Edge_ID, double> &emap_pn = edge_maps[i];
    std::map< Edge_ID, double> &diangle_pn = diangle_maps[i];
    std::set< Edge_ID> &eset_pn = _edges[i];

    const double *minvs = reinterpret_cast<double*>
      ( pane->dataitem(mindianglev->id())->pointer());
    const double *maxvs = reinterpret_cast<double*>
      ( pane->dataitem(maxdianglev->id())->pointer());
    const double *minta = reinterpret_cast<double*>
      ( pane->dataitem(mintrnangv->id())->pointer());
    double *maxta = reinterpret_cast<double*>
      ( pane->dataitem(maxtrnangv->id())->pointer());
    
    char *strong = reinterpret_cast<char*>
      ( pane->dataitem(_strong->id())->pointer());
    char *tranks = reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    char *qses = reinterpret_cast<char*>
      ( pane->dataitem(qsedges->id())->pointer());
    char *nnf = reinterpret_cast<char*>
      ( pane->dataitem(neighbor_features->id())->pointer());
 
    COM_assertion( emap_pn.size() == diangle_pn.size());
    // Loop through ebuf and put ridges edges onto emap_pn
    for ( std::map< Edge_ID, double>::const_iterator eit=emap_pn.begin(), 
	    rit=diangle_pn.begin(), eend=emap_pn.end(); eit!=eend; ++eit, ++rit) {
      Edge_ID eid = eit->first, eid_opp = (*pm_it)->get_opposite_real_edge( eid);
      COM_assertion( !eid.is_border());
      double feg = eit->second, fangle = rit->second, feg_opp;

      Element_node_enumerator ene( pane, eid.eid());
      int vi=eid.lid(), ne = ene.size_of_edges();
      int vindex = ene[vi]-1, windex = ene[ (vi+1)%ne]-1;

      if ( (nnf[vindex] && nnf[windex]) &&
	   ( strong[vindex] || !tranks[vindex] || tranks[vindex] == char(-1) ||
	     feg < 0 && fangle == -minvs[vindex] || 
	     feg > 0 && fangle == maxvs[vindex] ||
	     feg < -_tol_cos_osta && feg == minta[vindex] ||
	     feg > _tol_cos_osta && feg == maxta[vindex]) ||
	   nnf[vindex]==2 && nnf[windex]==2 && 
	   ( (feg_opp=emap_pn.find( eid_opp)->second)< 0 && 
	     fangle == -minvs[windex] ||
	     feg_opp > 0 && fangle == maxvs[windex]) &&
	   ( feg_opp < -_tol_cos_osta && feg_opp == minta[windex] ||
	     feg_opp > _tol_cos_osta && feg_opp == maxta[windex])) {
	eset_pn.insert( eid);
	qses[ eid.eid()-1] |= (1<<eid.lid());
      }
    } // for eit
  } // for i

  // Update border-edge count
  _surf->update_bdedge_bitmap( qsedges);

  // Select the edges whose both half-edges are quasi-strong.
  it = _panes.begin();
  pm_it=_surf->pm_begin();

  for ( int i=0, local_npanes = _panes.size(); i<local_npanes; 
	++i, ++it, ++pm_it) {
    std::set< Edge_ID> &eset_pn = _edges[i];

    // Select quasi-strong edges
    for ( std::set< Edge_ID>::iterator eit=eset_pn.begin(); 
	  eit!=eset_pn.end();) {
      Edge_ID eid = *eit, eid_opp = (*pm_it)->get_opposite_real_edge( eid);

      if ( !eid_opp.is_border() && eset_pn.find(eid_opp)==eset_pn.end() ||
	   eid_opp.is_border() && !(*pm_it)->get_bdedge_flag( eid_opp)) {
	std::set< Edge_ID>::iterator to_delete = eit;
	++eit;
	eset_pn.erase( to_delete);
      }
      else 
	++eit;
    }
  }

  // Note: neighbor_features is used as buffers in reclassify_ridge_vertices
  reclassify_ridge_vertices( 0, 1, neighbor_features, 
			     _tangranks, to_filter);
  
  // Currently the filtering step does not yet work in parallel.
  // Filter out obscure-ended ridges. Currently supports only sequential meshes
  if ( to_filter) for ( ;;) {
    std::vector< ObscendSet > obscends;

    int nobscended=build_ridge_neighbor_list( maxtrnangv,mintrnangv, edge_maps,
					      maxdianglev,mindianglev, 
					      diangle_maps, obscends);
    if ( nobscended==0) break;

    bool dropped = filter_obscended_ridge( edge_maps, diangle_maps, obscends);
    if ( !dropped) break;

    reclassify_ridge_vertices( 0, 0, neighbor_features, _tangranks, 1);
  }
  // Reclassify ridge edges to upgrade potentially missing corners.
  reclassify_ridge_vertices( 1, 0, neighbor_features, _tangranks, to_filter);

  Rocblas::copy( _tangranks, _ctangranks);
  // Update vertices on constraint boundary.
  if ( insert_boundary_edges( _ctangranks))
    // Reclassify ridge edges to upgrade potentially missing corners.
    reclassify_ridge_vertices(1, 1, neighbor_features, _ctangranks, to_filter);

  int nredges=0; // Number of ridge edges
  // Export ridge edges into list
  it = _panes.begin();
  pm_it=_surf->pm_begin();
  for ( int i=0, local_npanes = _panes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    std::set< Edge_ID> &eset_pn = _edges[i];

    // Allocate memory for _ridges and save them. 
    COM::DataItem *att_rdg = pane->dataitem(_ridges->id());

    int n=eset_pn.size();
    att_rdg->allocate( 2, std::max(n,att_rdg->size_of_items()), 0);
    int *rdgs = reinterpret_cast<int*>
      ( pane->dataitem(_ridges->id())->pointer()), ii=0;

    for ( std::set< Edge_ID>::const_iterator eit=eset_pn.begin();
	  eit!=eset_pn.end(); ++eit) {
      // Make sure that eid_opp is not border
      Edge_ID eid = *eit, eid_opp = (*pm_it)->get_opposite_real_edge( eid);

      // Obtain vertex indices.
      Element_node_enumerator ene( pane, eid.eid()); 
      int vindex = eid.lid(), neighbor = (vindex+1)%ene.size_of_edges();
      vindex = ene[vindex]; neighbor=ene[neighbor];
     
      if ( vindex<neighbor || eid_opp.is_border())
      { rdgs[ii++] = vindex; rdgs[ii++] = neighbor; }
    }
    
    att_rdg->set_size( ii/2);
    nredges += att_rdg->size_of_items();
  }

  if ( to_filter)
    std::cout << "Found " << nredges << " ridge edges " << std::endl;

  // Deallocate mindianglev and maxdianglev
  _buf->delete_dataitem( qsedges->name());
  _buf->delete_dataitem( neighbor_features->name());
  _buf->delete_dataitem( mintrnangv->name());
  _buf->delete_dataitem( maxtrnangv->name());
  _buf->delete_dataitem( mindianglev->name());
  _buf->delete_dataitem( maxdianglev->name());
  _buf->init_done( false);
}

int FaceOffset_3::
insert_boundary_edges( COM::DataItem *tr_attr) {
  std::vector< COM::Pane*>::iterator it = _panes.begin(), iend=_panes.end();
  Manifold::PM_iterator pm_it=_surf->pm_begin();

  int inserted=0;
  // Loop through the panes and its real elements
  for ( int i=0; it!=iend; ++i, ++it, ++pm_it) {
    COM::Pane *pane = *it;
    std::set< Edge_ID> &eset_pn = _edges[i];

    char *tranks = reinterpret_cast<char*>
      ( pane->dataitem(tr_attr->id())->pointer());
    const char *val_bndry=(const char*)pane->dataitem
      (_cnstr_bndry_edges->id())->pointer();

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1);
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      for ( int k=0, ne=ene.size_of_edges(); k<ne; ++k) {
	Edge_ID eid(j+1,k), eid_opp = (*pm_it)->get_opposite_real_edge( eid);
	
	// If a vertex has mixed geometric and topological edges, 
	// set it to be a corner
	if ( (*pm_it)->is_physical_border_edge( eid_opp) || 
	     (val_bndry && (val_bndry[j] & (1<<k)))) {
	  int vindex = ene[k]-1;

	  if (tranks[ vindex]==2 || eset_pn.find(eid)==eset_pn.end()) { 
	    if ( tranks[ vindex]==1) tranks[ vindex] = 0;
	    eset_pn.insert(eid); ++inserted; 
	  }
	}
      } // for k
    } // for j
  } // for it

  int inserted_total = inserted;
  if ( COMMPI_Initialized()) {
    MPI_Allreduce( &inserted, &inserted_total, 1, MPI_INT, MPI_SUM, 
		   _buf->get_communicator());
  }

  return inserted_total;
}

void FaceOffset_3::
reclassify_ridge_vertices( const bool upgrade_corners, 
			   const bool upgrade_ridge,
			   COM::DataItem *neighbor_feas, 
			   COM::DataItem *tr_attr, 
			   bool to_filter) {

  COM::DataItem *vecsums_buf = NULL;
  if ( upgrade_corners) {
    vecsums_buf = _buf->new_dataitem( "FO_vecsums_buf", 'n', COM_DOUBLE, 3, "");
    _buf->resize_array( vecsums_buf, 0);
  }

  int zero=0;
  Rocblas::copy_scalar( &zero, neighbor_feas);

  Vector_3 *pnts=NULL, *vss=NULL;
  // Loop through the panes and the candidate edges to count incident edges
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  for ( int ii=0, local_npanes = _panes.size(); 
	ii<local_npanes; ++ii, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    char *nnf = reinterpret_cast<char*>
      ( pane->dataitem(neighbor_feas->id())->pointer());
    if ( upgrade_corners) {
      pnts = reinterpret_cast<Vector_3*>
	(pane->dataitem(COM_NC)->pointer());
      vss = reinterpret_cast<Vector_3*>
	(pane->dataitem(vecsums_buf->id())->pointer());
    }
 
    std::set< Edge_ID> &eset_pn = _edges[ii];
    for ( std::set< Edge_ID>::const_iterator eit=eset_pn.begin(),
	    eend = eset_pn.end();  eit!=eend; ++eit) {
      Edge_ID eid = *eit;  COM_assertion(!eid.is_border());
      Element_node_enumerator ene( pane, eid.eid());
      int lid = eid.lid(), vindex = ene[lid]-1;
      ++nnf[vindex];

      Edge_ID eid_opp = (*pm_it)->get_opposite_real_edge( eid);
      bool is_border_edge = (*pm_it)->is_physical_border_edge( eid_opp);
      if ( is_border_edge)
	++nnf[ ene[(lid+1)%ene.size_of_edges()]-1];

      if ( upgrade_corners) {
	int windex = ene[(lid+1)%ene.size_of_edges()]-1;
	Vector_3 diff = pnts[ windex] - pnts[vindex];
	vss[vindex] += diff.normalize();
      
	if ( is_border_edge) vss[windex] -= diff;
      }
    }
  }
  _surf->reduce_on_shared_nodes( neighbor_feas, Manifold::OP_SUM);

  if ( upgrade_corners)
    _surf->reduce_on_shared_nodes( vecsums_buf, Manifold::OP_SUM);

  double tol_sqsin = 4*( 1 - _tol_cos_osta*_tol_cos_osta);
  int ncrns = 0; // Number of corners
  // Loop through all vertices to upgrade all the vertices incident on
  // one or more than two ridge edges.  
  it = _panes.begin(); 
  for ( int ii=0, local_npanes = _panes.size(); ii<local_npanes; ++ii, ++it) {
    COM::Pane *pane = *it;
    char *nnf = reinterpret_cast<char*>
      ( pane->dataitem(neighbor_feas->id())->pointer());
    char *tranks = reinterpret_cast<char*>
      ( pane->dataitem(tr_attr->id())->pointer());
    if ( upgrade_corners)
      vss = reinterpret_cast<Vector_3*>
	(pane->dataitem(vecsums_buf->id())->pointer());

    for ( int j=0, nj=pane->size_of_real_nodes(); j<nj; ++j) {
      if ( upgrade_corners && ( nnf[j]>=3 || nnf[j] == 2 &&
				vss[j].squared_norm() >= tol_sqsin))
	tranks[j] = 0;   // Upgrade joints to corners
      else if ( nnf[j] == 0 && std::abs(tranks[j]) == 1)
	tranks[j] = 2; 	// Downgrade ridge vertices without neighbors
      else if ( upgrade_ridge && tranks[j] == 2 && nnf[j] >= 2 || 
		upgrade_corners && tranks[j] == -1)
	tranks[j] = 1; // Upgrade ridge vertices with two or more neighbors

      ncrns += !tranks[j];
    }
  }

  // Reduce tangranks to ensure shared nodes have the same value across panes
  _surf->reduce_on_shared_nodes( tr_attr, Manifold::OP_MIN);
  
  if ( upgrade_corners) {
    _buf->delete_dataitem( vecsums_buf->name());
    _buf->init_done( false);

    if ( to_filter)
      std::cout << "Found " << ncrns << " corner vertices" << std::endl;
  }
}

struct RidgeNeighbor {
  RidgeNeighbor( int v, Edge_ID e) : vid(v), eid(e) {}

  int       vid;
  Edge_ID   eid;
};

// Construct the neighbor vertices of ridges.
// This routine only supports single-pane meshes.
int FaceOffset_3::
build_ridge_neighbor_list( const COM::DataItem *maxtrnangv,
			   const COM::DataItem *mintrnangv,
			   const std::vector<std::map<Edge_ID, double> > &edge_maps,
			   const COM::DataItem *maxdianglev, 
			   const COM::DataItem *mindianglev,
			   const std::vector<std::map<Edge_ID, double> > &diangle_maps,
			   std::vector< ObscendSet > &obscends) {

  COM_assertion( _panes.size()==1); // Precondition
  
  // Loop through the panes and ridge edges
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  COM::Pane *pane = _panes[0];

  // Build ACH list
  typedef std::map< int, std::list< RidgeNeighbor> >  ACH_Lists;
  ACH_Lists ach; // ACH list

  const std::set< Edge_ID> &eset_pn = _edges[0];  // List of ridge edges
  for ( std::set< Edge_ID>::const_iterator eit=eset_pn.begin(),
	  eend=eset_pn.end(); eit!=eend; ++eit) {
    // Make sure that eid_opp is not border
    Edge_ID eid = *eit, eid_opp = (*pm_it)->get_opposite_real_edge( eid);
    COM_assertion( eset_pn.find(eid_opp)!=eset_pn.end() &&
		   !eid.is_border());
    Element_node_enumerator ene( pane, eid.eid()); 
    int lid=eid.lid(), neighbor = (lid+1)%ene.size_of_edges();

    // Fill in ridge-neighbor into _ridgeneighbors.
    int vid = ene[lid], neighbor_id=ene[neighbor];
    ach[vid].push_back( RidgeNeighbor( neighbor_id, eid));
  }

  // initialize storage
  int zero=0;
  Rocblas::copy_scalar( &zero, _ridgeneighbors);
  obscends.clear(); obscends.resize( _panes.size());

  // Check ACH list and build rdgngbs
  RidgeNeighbor *rdgngbs = reinterpret_cast<RidgeNeighbor*>
    ( pane->dataitem(_ridgeneighbors->id())->pointer());
  const Point_3 *pnts = reinterpret_cast<Point_3*>
    (pane->dataitem(COM_NC)->pointer());
  const char *tranks = reinterpret_cast<char*>
    ( pane->dataitem(_tangranks->id())->pointer());
  const double *minvs = reinterpret_cast<double*>
    ( pane->dataitem(mindianglev->id())->pointer());
  const double *maxvs = reinterpret_cast<double*>
    ( pane->dataitem(maxdianglev->id())->pointer());
  const double *minta = reinterpret_cast<double*>
    (pane->dataitem(mintrnangv->id())->pointer());
  const double *maxta = reinterpret_cast<double*>
    ( pane->dataitem(maxtrnangv->id())->pointer());
  char *strong = reinterpret_cast<char*>
    ( pane->dataitem(_strong->id())->pointer());
  const std::map< Edge_ID, double> &diangle_pn = diangle_maps[0];
  const std::map< Edge_ID, double> &emap_pn = edge_maps[0];
  ObscendSet &obscend_pn = obscends[0];

  for ( ACH_Lists::iterator ach_it = ach.begin(); ach_it!=ach.end(); ++ach_it) {
    // Check angle defect and turning angle
    int vid= ach_it->first, index = (vid-1)*2;
    int nedges = ach_it->second.size();
    std::list< RidgeNeighbor>::iterator ait = ach_it->second.begin();

    if (nedges>2) {
      // Insert into rdgngbs
      for (; ait!=ach_it->second.end(); ) {
	// Remove non-joint halfedges
	Edge_ID eid = ait->eid;
	double fangle = diangle_pn.find( eid)->second;
	double feg = emap_pn.find( eid)->second;

	// If disjoint, then add edge into obscend_pn.
	if ( fangle < _tol_eangle &&
	     ( strong[vid-1] || tranks[vid-1]==1 && 
	       ( (-feg < _tol_cos_osta || feg != minta[vid-1]) &&
		 ( feg < _tol_cos_osta || feg != maxta[vid-1]) ||
		 ( feg<0 || fangle != maxvs[vid-1]) && 
		 ( feg>0 || fangle != -minvs[vid-1]))) ||
	     fangle < _tol_angle_strong && 
	     ( strong[vid-1]==2 || tranks[vid-1]==1 &&
	       ( std::abs(feg) < _tol_cos_osta ||
		 feg != minta[vid-1] && feg != maxta[vid-1]) && 
	       ( feg<0 || fangle != maxvs[vid-1]) && 
	       ( feg>0 || fangle != -minvs[vid-1]))) {
	  obscend_pn[ ait->eid] = std::make_pair( vid, ait->vid);
	  ait = ach_it->second.erase( ait);
        }
        else
	  ++ait;
      }

      rdgngbs[index].vid = 0; rdgngbs[index+1].vid = -1; 
    }

    // Classify dangling and semi-joint halfedges
    ait = ach_it->second.begin();
    if ( nedges == 1) {
      // Insert into obscure list to indicate dangling edges
      rdgngbs[index] = *ait; rdgngbs[index+1].vid=0;
      Edge_ID eid = ait->eid;
      double fangle = diangle_pn.find( eid)->second;

      if ( fangle < _tol_angle_strong || tranks[vid-1])
	obscend_pn[ eid] = std::make_pair( ach_it->first, ait->vid);
    }
    else if ( nedges == 2) {
      std::list< RidgeNeighbor>::iterator ait2 = ait; ++ait2;
      int v1 = ait->vid, v2 = ait2->vid;
      double fangle1 = diangle_pn.find( ait->eid)->second;
      double fangle2 = diangle_pn.find( ait2->eid)->second;

      if ( (fangle1<_tol_angle_strong || fangle2<_tol_angle_strong) &&
	   ( !tranks[vid-1] || eval_angle
	     ( pnts[v1-1]-pnts[vid-1], pnts[vid-1]-pnts[v2-1]) > _tol_turn_angle)) {
	// Do not insert the edges into rdgngbs but insert into obscend_pn
	obscend_pn[ ait->eid] = std::make_pair( vid, v1);
	obscend_pn[ ait2->eid] = std::make_pair( vid, v2);
      }
      else if (!rdgngbs[index+1].vid)
      { rdgngbs[index] = *ait; rdgngbs[index+1] = *ait2; }
    }
  }
    
  int  n_obscended = obscends[0].size();
  std::cerr << "There are " << n_obscended 
	    << " obscure-ended edges" << std::endl;
  return n_obscended;
}

bool FaceOffset_3::
filter_obscended_ridge( const std::vector<std::map<Edge_ID, double> > &edge_maps,
			const std::vector<std::map<Edge_ID, double> > &diangle_maps,
			const std::vector< ObscendSet > &obscends) {
  // Loop through the panes and ridge edges
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  int local_npanes = _panes.size();
  std::vector< ObscendSet > todelete(local_npanes);

  for ( int i=0; i<local_npanes; ++i, ++it, ++pm_it) { 
    const ObscendSet &obscend_pn = obscends[i];
    ObscendSet &todelete_pn = todelete[i];
    if ( obscend_pn.empty()) continue;

    COM::Pane *pane = *it;
    const char *tranks = reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    char *strong = reinterpret_cast<char*>
      ( pane->dataitem(_strong->id())->pointer());
    const RidgeNeighbor *rdgngbs = reinterpret_cast<RidgeNeighbor*>
      ( pane->dataitem(_ridgeneighbors->id())->pointer());
    const std::map< Edge_ID, double> &diangle_pn = diangle_maps[i];

    for ( ObscendSet::const_iterator rit=obscend_pn.begin(),
	    rend=obscend_pn.end(); rit!=rend; ++rit) {
      // Loop through edges to count
      int cur=rit->second.first, next=rit->second.second, start=cur;

      Edge_ID eid = rit->first;
      bool is_dangling = ( rdgngbs[ (cur-1)*2].vid==next && 
			   !rdgngbs[ (cur-1)*2+1].vid ||
			   !rdgngbs[ (cur-1)*2].vid && 
			   !rdgngbs[ (cur-1)*2+1].vid);
      bool is_not_joint1 = ( rdgngbs[ (cur-1)*2].vid!=next && 
			     rdgngbs[ (cur-1)*2+1].vid!=next), is_not_joint2=false;

      if ( !is_dangling && !is_not_joint1) continue;
      bool is_corner=false, is_strong1=strong[cur-1] || !tranks[cur-1];
      bool is_strong2=strong[next-1] || !tranks[next-1];

      int count_ks = 0;
      for ( ;!(is_corner = !tranks[next-1]) && next && next!=start; ) {
	count_ks += diangle_pn.find( eid)->second>_tol_kangle;
	
	int  index = (next-1)*2;
	if ( rdgngbs[ index].vid == cur) { 
	  cur = next; next = rdgngbs[index+1].vid; eid = rdgngbs[index+1].eid;
	  is_strong2=strong[next-1] || !tranks[next-1];
	  is_not_joint2 = false; 
	}
	else {
	  Edge_ID eid_opp = (*pm_it)->get_opposite_real_edge( eid);

	  is_not_joint2 = (rdgngbs[ index+1].vid!=cur) && 
	    (obscend_pn.find(eid_opp) != obscend_pn.end());

	  if ( rdgngbs[ index+1].vid == cur) { 
	    cur = next; next = rdgngbs[index].vid; eid = rdgngbs[ index].eid; 
	    is_strong2=strong[next-1] || !tranks[next-1];
	  }
	  else
	  { cur = next; next = 0; }
	}
      }

      // Check whether the edge is quasi-obscure.
      if ( (is_dangling || is_not_joint1 && is_not_joint2) && count_ks<=_tol_kstrong ||
	   is_strong1 && is_strong2 && count_ks==0) {
	todelete_pn.insert( *rit);
      }
    }
  }
  
  return remove_obscure_curves( todelete);
}

// Remove edges in quasi-obscure curves from candidate-edge list.
int FaceOffset_3::
remove_obscure_curves( const std::vector< ObscendSet > &obscends) {
  // Loop through the panes and ridge edges
  Manifold::PM_iterator pm_it=_surf->pm_begin();

  // Remove quasi-obscure edges.
  int dropped = 0;
  for ( int i=0, local_npanes = _panes.size(); 
	i<local_npanes; ++i, ++pm_it) { 
    std::set< Edge_ID> &eset_pn = _edges[i];
    const ObscendSet &obscend_pn = obscends[i];
    COM::Pane *pane = _panes[i];

    char *tranks = reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    const RidgeNeighbor *rdgngbs = reinterpret_cast<RidgeNeighbor*>
      ( pane->dataitem(_ridgeneighbors->id())->pointer());

    std::cout << "There are " << obscend_pn.size() 
	      << " obscure curves. " << std::endl;

    for ( ObscendSet::const_iterator rit=obscend_pn.begin(),
	    rend=obscend_pn.end(); rit!=rend; ++rit) {
      dropped += eset_pn.erase( rit->first);
      eset_pn.erase( (*pm_it)->get_opposite_real_edge( rit->first));

      int cur = rit->second.first, next = rit->second.second, start=cur;
      for ( ; tranks[next-1]; ) {
	int  index = (next-1)*2;
	Edge_ID eid;

	if ( rdgngbs[ index].vid == cur) 
	{ cur = next; next = rdgngbs[index+1].vid; eid = rdgngbs[ index+1].eid; }
	else if ( rdgngbs[ index+1].vid == cur) 
	{ cur = next; next = rdgngbs[index].vid; eid = rdgngbs[ index].eid; }
	else
	{ cur = next; next = 0; }
	if ( next==0 || next==start) break;

	dropped += eset_pn.erase( eid);
	eset_pn.erase( (*pm_it)->get_opposite_real_edge( eid));
      }
    }
  }

  if ( dropped)
    std::cerr << "Dropped " << dropped << " false candidate edges. " << std::endl;

  return dropped;
}

PROP_END_NAMESPACE






