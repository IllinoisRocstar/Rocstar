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
// $Id: Remesher_Simmetrix.C,v 1.6 2008/12/06 08:45:28 mtcampbe Exp $

/** This file provides a wrapper to invoke GeomSim from Simmetrix to
 *  remesh a surface mesh */

#include "Remesher_Simmetrix.h"
#include "MeshSim.h"
#include "SimModel.h"
#include "SimDiscrete.h"

PROP_BEGIN_NAMESPACE

int Remesher_Simmetrix::instances = 0;

/** Initialize the remesher. In particular, initialize MeshSim.
 */
Remesher_Simmetrix::Remesher_Simmetrix( const char *logfile)
  : _size_type(0), _size_val(0.), _fangle(0), _shrtRatio(5) {
  
  if ( logfile) Sim_logOn( logfile);

  if ( instances == 0) {
    MS_init(); // start up MeshSim
    Sim_registerKey("cen002 surface 20060331 0 mTcaTqnWjLIUId26P29r3w=="); // for surface
    Sim_registerKey("cen002 discrete 20060331 0 8VnUzbSAnldsPdfPbpWFVw=="); // for discrete
    
    SimDiscrete_start(0);  
  }
  ++instances;
}

/** Finalize the remesher. In particular, finalize MeshSim.
 */
Remesher_Simmetrix::~Remesher_Simmetrix() {
  --instances;

  if ( instances==0) {
    SimDiscrete_stop(0);
    MS_exit();
    Sim_logOff();
    Sim_unregisterAllKeys(); 
  }
}

void *Remesher_Simmetrix::
window_to_simmesh( const COM::Window *outwin) {
  // Obtain the list of panes
  std::vector<const COM::Pane*> panes; 
  outwin->panes(panes);
  COM_assertion_msg( panes.size() == 1, "Window must be a serial one");

  // Get the sizes
  int numVerts = panes[0]->size_of_nodes();
  int numElems = panes[0]->size_of_elements();

  // Get the pointers to 
  std::vector<const COM::Connectivity*> elems;
  panes[0]->elements( elems);
  COM_assertion_msg( elems.size() == 1 && 
		     elems[0]->element_type()==COM::Connectivity::TRI3,
		     "Window must be a serial one");

  pMesh mesh = M_new(0,0);
  // Initialize pMesh.
  std::vector<int> elemTypes(numElems, 5 /* Triangles */);
  /* offset nodal IDs to start from 0 */
  int *tris = const_cast<int*>(elems[0]->pointer());
  std::vector<int> tris_vec( numElems*3,0);
  for ( int i=0; i<numElems*3; ++i)
    tris_vec[i] = tris[i]-1;
  
  double *coors = const_cast<double*>(panes[0]->coordinates());
  M_importFromData( mesh, numVerts, coors, numElems, &elemTypes[0], 
		    &tris_vec[0], NULL, NULL);
  
  return mesh;
}

void Remesher_Simmetrix::
simmesh_to_window( void *m, COM::Window *outwin) {
  pMesh mesh = (pMesh)m;

  // Obtain the panes
  std::vector<COM::Pane*> panes; 
  outwin->panes(panes);
  COM_assertion_msg( panes.size() == 1, "Window must be a serial one");
  int pid = panes[0]->id();

  // Obtain coordinates
  COM::DataItem *coor = panes[0]->dataitem( COM::COM_NC);
  outwin->set_size( coor->name(), pid, M_numVertices(mesh));

  double *coors;
  outwin->resize_array( coor, (void**)&coors);

  VIter vIter = M_vertexIter(mesh);
  pVertex vertex;

  for (int i=0; (vertex=VIter_next(vIter)); ++i) {
    EN_setID((pEntity)vertex, i+1);
    V_coord(vertex, &coors[3*i]);
  }
  VIter_delete(vIter);

  // Obtain connectivity
  std::vector<COM::Connectivity*> elems;
  panes[0]->elements( elems);

  outwin->set_size( elems[0]->name(), pid, M_numFaces(mesh));

  int *tris;
  outwin->resize_array( elems[0], (void**)&tris);

  FIter fIter = M_faceIter(mesh);
  pFace face;
  for (int i=0; (face = FIter_next(fIter)); ++i){
    pPList fVerts = F_vertices(face,1);
    int j;
    for(j=0; j < 3; j++){
      vertex = (pVertex)PList_item( fVerts, j);
      tris[3*i+j] = EN_id((pEntity)vertex);
    }
  }
  FIter_delete(fIter);  
}

void Remesher_Simmetrix::
remesh_serial( Manifold *wm, COM::DataItem *mesh_out,
	       double lave, double fangle) {

  std::cout << "Setting average edge length to " << lave << std::endl;
  std::cout << "Setting feature angle to " << fangle << std::endl;

  // Set options
  set_global_mesh_size( 1, lave);
  set_fangle( fangle);
  set_short_ratio( 5);


  COM::Window *outwin = mesh_out->window();

  // Serialize input mesh and put into output mesh
  wm->serialize_window( outwin);

  pMesh mesh = (pMesh)window_to_simmesh( outwin);
  pDiscreteModel model = DM_createFromMesh( mesh, true);

  if ( _fangle>0) {
    // Perform feature detection
    double angleRadian = _fangle/180*3.14159265358979;

    DM_findEdgesByFaceNormals( model, angleRadian);
    DM_eliminateDanglingEdges( model);
  }
  DM_completeTopology( model);
  
  // Obtain a new copy of the mesh for alteration
  pMesh newMesh = DM_getMeshCopy( model);
  pACase newCase = MS_newMeshCase(model); 
	
  // set the desired mesh size 
  if ( _size_type != 0) 
  {
     pModelMember mdomain = GM_domain(model);	
     MS_setMeshSize(newCase,mdomain, _size_type, _size_val,0);
  }

  // Modify the surface mesh so that it conforms to the new size
  DM_modifySurfaceMesh( newMesh, newCase ,0,0,0,0, /* swap, split, min, and max angles*/
			0, 0, /* Collapes face and edge angles */
			0, _shrtRatio, 0);
  
  simmesh_to_window( newMesh, outwin);
  
  M_release(newMesh);
  GM_release((SGModel *) model);
}

PROP_END_NAMESPACE






