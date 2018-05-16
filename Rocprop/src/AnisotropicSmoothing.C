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
// $Id: AnisotropicSmoothing.C,v 1.6 2008/12/06 08:45:27 mtcampbe Exp $

#include "FaceOffset_3.h"
#include "Rocblas.h"
#include "Generic_element_2.h"
#include "Rocsurf.h"

PROP_BEGIN_NAMESPACE

#define ANISOTROPIC            0
#define MIN_ANGLE_WEIGHTED     1

bool FaceOffset_3::
is_strong( const Vector_3 *es, const Vector_3 &lambdas,
	   const Vector_3 &mn, int nrank) const {
  if ( lambdas[1]>_dir_thres_strong*lambdas[0]) return true;
  
  return is_acute_ridge( es, lambdas, mn, nrank) || 
    is_acute_corner( es, lambdas, mn, nrank);
}


// Obtain weight for nonuniform mesh fairing.
double FaceOffset_3::
get_nonuniform_weight( const Vector_3 &lambdas, int trank) {

  const double eps=1.e-6;
  double w = 1/lambdas[0];
  
  if (trank == 2) return w;
  else if ( trank==1 ) return eps*w;
  else return eps*eps*w;
}

void FaceOffset_3::
get_tangents( const Vector_3 *es, const Vector_3 &lambdas,
	      const Vector_3 &mn, int nrank,
	      Vector_3 vs[2], double scales[2], int *is_strong) const {
  vs[0] = es[1];
  vs[1] = es[2];

  // Bound from two ends appear to be the most robust.
  const double max_scale=0.07, min_scale=0.005;

  if ( scales) {
    scales[0] = std::min( max_scale, std::max(lambdas[1]/lambdas[0], min_scale));
    scales[1] = std::min( max_scale, std::max(lambdas[2]/lambdas[0], min_scale));

    scales[0] = std::sqrt( scales[0]);
    scales[1] = std::sqrt( scales[1]);
    
    if ( is_strong) {
      bool acute_ridge = FaceOffset_3::is_acute_ridge( es, lambdas, mn, nrank);
      bool acute_corner = FaceOffset_3::is_acute_corner( es, lambdas, mn, nrank);

      *is_strong = acute_ridge || acute_corner || 
	lambdas[1]>_dir_thres_strong*lambdas[0];
    }
  }
}

void FaceOffset_3::
compute_anisotropic_vertex_centers( const COM::DataItem *nodal_disps) {

  const double zero = 0., eps = 1.e-100;

  Rocblas::copy_scalar( &zero, _vcenters);
  Rocblas::copy_scalar( &eps, _weights);
  
  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  for ( int i=0, local_npanes = _panes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;

    const Point_3 *pnts = reinterpret_cast<const Point_3*>
      (pane->dataitem(COM_NC)->pointer());
    const Vector_3 *disps = nodal_disps?reinterpret_cast<const Vector_3*>
      (pane->dataitem(nodal_disps->id())->pointer()):NULL;
    const Point_3 *fcnts = reinterpret_cast<const Point_3*>
      ( pane->dataitem(_facecenters->id())->pointer());
    Vector_3 *vcnts = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_vcenters->id())->pointer());
    double   *ws = reinterpret_cast<double*>
      ( pane->dataitem(_weights->id())->pointer());

#if ANISOTROPIC    
    const Vector_3 *As = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_As->id())->pointer());
    const Vector_3 *bs = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_bmedial->id())->pointer());

    const char *tranks = reinterpret_cast<const char*>
      ( pane->dataitem(_tangranks->id())->pointer());

    const char *strong = reinterpret_cast<char*>
      ( pane->dataitem(_strong->id())->pointer());
#endif

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1);
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne = ene.size_of_edges();
      const Point_3  &cnt = fcnts[j];
      Vector_3 anipnt;

#if ANISOTROPIC
      bool has_feature=false;
      Vector_3 Ae[3], mn;   // Matrix at an edge

      // First, compute medial quadric for the face
      Ae[0] = Ae[1] = Ae[2] = mn = Vector_3(0,0,0);

      for ( int k=0; k<ne; ++k) {
	int vindex = ene[k]-1;
	has_feature |= tranks[vindex]!=2;
	
	// Obtain the quadric at the face 
	const Vector_3 *As_v=&As[3*vindex];
	Ae[0] += As_v[0];
	Ae[1] += As_v[1];
	Ae[2] += As_v[2];
	    
	mn += bs[vindex];
      }

      // Perform eigen-analysis
      Vector_3 lambdas = mn;
      int nrank = eigen_analyze_vertex( Ae, lambdas);
      
      Vector_3 vs_e[2];
      double   scales_e[2];

      mn = Ae[0]; nrank = 1;
      get_tangents( Ae, lambdas, mn, nrank, vs_e, scales_e);

      Vector_3 mn_face;
      Vector_3 vs_f[2];
      double   scales_f[2];
      if ( has_feature) {
	// Compute eigen-vector of weighted-sum of the tensors.
	Ae[0] = Ae[1] = Ae[2] = Vector_3(0,0,0);
	for ( int k=0; k<ne; ++k) {
	  int vindex = ene[k]-1;

	  double w=1;
	  switch (tranks[ vindex]) {
	  case 0: {
	    w = 1.e-12;
	    break;
	  }
	  case 1: {
	    if ( strong[vindex]) { 
	      w = 1.e-6;
	      break;
	    }
	  }
	  default: w = 1.;
	  }
	  
	  // Obtain the quadric at the face 
	  const Vector_3 *As_v=&As[3*vindex];
	  Ae[0] += w*As_v[0];
	  Ae[1] += w*As_v[1];
	  Ae[2] += w*As_v[2];
	}
	
	// Perform eigen-analysis
	eigen_analyze_vertex( Ae, lambdas);
	mn_face = Ae[0];
	
	get_tangents( Ae, lambdas, mn_face, 1, vs_f, scales_f);
      }
      else {
	mn_face = Ae[0];

	scales_f[0] = scales_e[0]; scales_f[1] = scales_e[1]; 
	vs_f[0] = vs_e[0]; vs_f[1] = vs_e[1]; 
      }
#endif

      int uindex=ene[ne-1]-1, vindex=ene[0]-1;
      // Loop through all vertices of current face.
      for ( int k=0; k<ne; ++k) {
	int windex = ene[(k+1==ne)?0:k+1]-1;

	Point_3 ps[3] = { pnts[uindex], pnts[vindex], pnts[windex] };
	if ( disps) { 
	  ps[0] += disps[uindex]; ps[1] += disps[vindex]; 
	  ps[2] += disps[windex]; 
	}

#if ANISOTROPIC	
	double scales_v[2];
	Vector_3 vs_v[2];

	// Note: is_border node.
	if ( _feature_layer && 
	     ( tranks[uindex]==2 || tranks[windex]==2)) {
	  scales_v[0] = scales_e[0]; scales_v[1] = scales_e[1]; 
	  vs_v[0] = vs_e[0]; vs_v[1] = vs_e[1];
	}
	else  {
	  scales_v[0] = scales_f[0]; scales_v[1] = scales_f[1]; 
	  vs_v[0] = vs_f[0]; vs_v[1] = vs_f[1];
	}
	
	const Vector_3 diff = cnt-ps[1];
	anipnt = scales_v[0]*(diff*vs_v[0])*vs_v[0]+
	  scales_v[1]*(diff*vs_v[1])*vs_v[1];
	if ( disps) anipnt += disps[vindex];
#else
	anipnt = cnt-ps[1];
#endif

#if MIN_ANGLE_WEIGHTED
      	double w = std::min( eval_angle( ps[1]-ps[0], ps[2]-ps[0]),
			     eval_angle( ps[1]-ps[2], ps[0]-ps[2]));
#else
	double  w=1;
#endif

	vcnts[vindex] += w*anipnt;
	ws[vindex] += w;

	uindex=vindex; vindex=windex;
      }
    }
  }

  _surf->reduce_on_shared_nodes( _vcenters, Manifold::OP_SUM);
  _surf->reduce_on_shared_nodes( _weights, Manifold::OP_SUM);
  Rocblas::div( _vcenters, _weights, _vcenters);
}

void FaceOffset_3::
denoise_surface( const COM::DataItem *nodal_disps, 
		 COM::DataItem *normal_motion) {
  
  const double zero = 0., eps = 1.e-100;

  Rocblas::copy_scalar( &zero, normal_motion);
  Rocblas::copy_scalar( &eps, _weights);

  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  for ( int i=0, local_npanes = _panes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    
    const Vector_3 *As = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_As->id())->pointer());
    const Vector_3 *vnrms = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_vnormals->id())->pointer());

    const Point_3 *pnts = reinterpret_cast<const Point_3*>
      (pane->dataitem(COM_NC)->pointer());
    const Vector_3 *disps = nodal_disps ? reinterpret_cast<const Vector_3*>
      (pane->dataitem(nodal_disps->id())->pointer()) : NULL;
    const Point_3 *fcnts = reinterpret_cast<const Point_3*>
      ( pane->dataitem(_facecenters->id())->pointer());
    const char *tranks = reinterpret_cast<const char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    const Vector_3 *lambdas = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_eigvalues->id())->pointer());

    Vector_3 *nmotion = normal_motion ? reinterpret_cast<Vector_3*>
      ( pane->dataitem(normal_motion->id())->pointer()) : NULL;
    double *ws = normal_motion ? reinterpret_cast<double*>
      ( pane->dataitem(_weights->id())->pointer()) : NULL;
    const double *farea = normal_motion ? reinterpret_cast<double*>
      ( pane->dataitem(_faceareas->id())->pointer()) : NULL;
    const char *weak = reinterpret_cast<const char*>
      ( pane->dataitem(_weak->id())->pointer());
    const int *cnstrs = _cnstr_nodes ? reinterpret_cast<const int*>
      ( pane->dataitem(_cnstr_nodes->id())->pointer()) : NULL;

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1);
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne = ene.size_of_edges();
      const Point_3  &cnt = fcnts[j];

      Vector_3 Ae[3];   // Matrix at an edge
      Ae[0] = Ae[1] = Ae[2] = Vector_3(0,0,0);
      bool is_weak = false;

      for ( int k=0; k<ne; ++k) {
	int vindex = ene[k]-1;
	if ( weak[vindex]) { is_weak = true; break; }
      }
      if ( !is_weak) continue;


      for ( int k=0; k<ne; ++k) {
	int vindex = ene[k]-1;
	
	// Obtain the quadric at the face 
	const Vector_3 *As_v=&As[3*vindex];

	double w = get_nonuniform_weight
	  ( lambdas[vindex], tranks[vindex]-(cnstrs && cnstrs[vindex]));

	Ae[0] += w*As_v[0];
	Ae[1] += w*As_v[1];
	Ae[2] += w*As_v[2];
      }

      // Perform eigen-analysis to obtain diffused face normal
      Vector_3 ls;
      compute_eigenvectors( Ae, ls);
      Vector_3 mn_face = Ae[0];

      for ( int k=0; k<ne; ++k) {
	int vindex = ene[k]-1;
	
	if ( weak[vindex]) {
	  Vector_3 diff = (cnt-pnts[vindex]);
	  if ( disps) diff -= disps[vindex];

	  double offset = diff*mn_face, cos_a=mn_face*vnrms[vindex];
	  if (cos_a<0) offset = -offset;

	  // Use eigenvalue to control amount of dissipation
	  double w = farea[j];
	  double a = lambdas[vindex][1]/lambdas[vindex][0];
	  nmotion[vindex] += a*(w*offset)*vnrms[vindex];

	  ws[vindex] += w;
	}
      }
    }
  }

  _surf->reduce_on_shared_nodes( normal_motion, Manifold::OP_SUM);
  _surf->reduce_on_shared_nodes( _weights, Manifold::OP_SUM);
  Rocblas::div( normal_motion, _weights, normal_motion);
}

PROP_END_NAMESPACE
  






