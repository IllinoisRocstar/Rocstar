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
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <list>
#include <cassert>

#include "TRAIL_UnixUtils.H"

//#ifdef _TRAIL_MPI_
#include "mpi.h"
//#endif
#include "GEM.H"
#include "TRAIL.H"
//#ifdef _ROCSTAR_X_
#include "com.h"
COM_EXTERN_MODULE(SimOUT);
//#endif

using namespace std;


void
GEM_Partition::report()
{
  if(_out){
    *_out  << "GEM_Partition ID: " << _id << endl
	   << "Mesh Entities: Total Real Ghosts" << endl
	   << "Number of Nodes:    (" << _nc.size()/3 << "," 
	   << _nc.size()/3 - _ngnodes << "," << _ngnodes << ")" 
	   << endl
	   << "Number of Elements: (" << nelem() << ","
	   << nelem() - _ngtet - _ngpyr - _ngpris - _nghex << ","
	   << _ngtet + _ngpyr + _ngpris + _nghex << ")" << endl
	   << "Cell Mapping: (" << _cell_ordering[0] << "," 
	   << _cell_ordering[1] << "," << _cell_ordering[2] << ","
	   << _cell_ordering[3] << ")" << endl;
    if(_tetconn.size() > 0) 
      *_out << "Tets:     (" << _tetconn.size()/4  << "," 
	    << _tetconn.size()/4 - _ngtet << "," << _ngtet  << ")" << endl;
    if(_pyrconn.size() > 0)
      *_out << "Pyramids: (" << _pyrconn.size()/5  << "," 
	    << _pyrconn.size()/5 - _ngpyr << "," << _ngpyr  << ")" << endl;
    if(_prisconn.size() > 0)
      *_out   << "Prisms:   (" << _prisconn.size()/6 << "," 
	      << _prisconn.size()/6 - _ngpris << "," << _ngpris << ")" << endl;
    if(_hexconn.size() > 0)
      *_out   << "Hexes:    (" << _hexconn.size()/8  << "," 
	      << _hexconn.size()/8 - _nghex << "," << _nghex  << ")" << endl;
    if(_debug){
      debug(false);
      *_out << "Cell debugging: " << endl;
      unsigned int nel,el;
      unsigned int nlin = 5;
      if(_tetconn.size() > 0) {
	*_out << "Tet Cell IDs: " << endl;
	nel = _tetconn.size()/4;
	nlin = 5;
	el = 0;
	while(el < nel){
	  *_out << setw(12) << Elem2Cell(std::make_pair((unsigned int)1,++el));
	  if(!(el%nlin))
	    *_out << endl;
	}
	if(el%nlin) *_out << endl;
      }
      if(_pyrconn.size() > 0){
	*_out << "Pyr Cell IDs: " << endl;
	nel = _pyrconn.size()/5;
	el = 0;
	while(el < nel){
	  *_out << setw(12) << Elem2Cell(std::make_pair((unsigned int)2,++el));
	  if(!(el%nlin))
	    *_out << endl;
	}
	if(el%nlin) *_out << endl;
      }
      if(_prisconn.size() > 0){
	*_out << "Pris Cell IDs:" << endl;
	nel = _prisconn.size()/6;
	el = 0;
	while(el < nel){
	  *_out << setw(12) << Elem2Cell(std::make_pair((unsigned int)3,++el));
	  if(!(el%nlin))
	    *_out << endl;
	}
	if(el%nlin) *_out << endl;
      }
      if(_hexconn.size() > 0){
	*_out << "Hex Cell IDs: " << endl;
	nel = _hexconn.size()/8;
	el = 0;
	while(el < nel){
	  *_out << setw(12) << Elem2Cell(std::make_pair((unsigned int)4,++el));
	  if(!(el%nlin))
	    *_out << endl;
	}
	if(el%nlin) *_out << endl;
      }
      debug(true);
    }
  }
}

void
GEM_Partition::report_partition_boundaries()
{
  if(_out){
    *_out << "==================================================" << endl  
	  << "Number of partition boundaries: " << _pb.size() << endl; 
    unsigned int border = 0;
    while(border < _pb.size()){
      GEM_PartitionBoundary &fb = _pb[border++];
      *_out << "------------------------------------------------" << endl  
	    << " Partition Boundary #" << border+1 << endl;
      fb.report();
      *_out << "------------------------------------------------" << endl;
    }
    *_out << "==================================================" << endl;
  }
}


void
GEM_PartitionBoundary::report()
{
  if(_out){
    *_out  << "Remote Partition: " << _rpart << endl
	   << "Nodes: (" << _sharenodes.size() << "," 
	   << _sendnodes.size() << ","  << _recvnodes.size() 
	   << ")" << endl
	   << "Cells: (" << _sendcells.size() << "," 
	   << _recvcells.size() << ")" << endl;
  }
}

void 
GEM_Partition::report_domain_boundaries()
{
  if(_out){
    *_out << "==================================================" << endl  
	  << "Number of Domain Boundaries: " << _db.size() << endl;
    unsigned int patch = 0;
    while(patch < _db.size()){
      GEM_DomainBoundary &fp = _db[patch++];
      *_out << "------------------------------------------------" << endl;
      fp.report();
      *_out << "------------------------------------------------" << endl;
    }
    *_out << "==================================================" << endl;  
  }
}


void
GEM_DomainBoundary::report()
{
  if(_out){
    *_out << "Domain Boundary  ID: " << _id << endl
	  << " Triangles:  (" << _triconn.size()/3 << "," << _ngtri << ")" 
	  << endl
	  << " Quads:  (" << _quadconn.size()/4 << "," << _ngquad 
	  << ")" << endl;
  }
}

void 
GEM_PartitionBoundary::populate(int rpid, int nnshared,int nnsend,int nnrecv,
				int ncsend,int ncrecv, int *sharedn, 
				int *sendn,int *recvn, int *sendc,int *recvc)
{
  int indy = 0;
  _rpart = rpid;
  _sendcells.resize(ncsend);
  while(indy < ncsend){
    _sendcells[indy] = sendc[indy];
    indy++;
  }
  _recvcells.resize(ncrecv);
  indy = 0;
  while(indy < ncrecv){
    _recvcells[indy] = recvc[indy];
    indy++;
  }
  indy = 0;
  _sendnodes.resize(nnsend);
  while(indy < nnsend){
    _sendnodes[indy] = sendn[indy];
    indy++;
  }
  indy = 0;
  _recvnodes.resize(nnrecv);
  while(indy < nnrecv){
    _recvnodes[indy] = recvn[indy];
    indy++;
  }
  indy = 0;
  _sharenodes.resize(nnshared);
  while(indy < nnshared){
    _sharenodes[indy] = sharedn[indy];
    indy++;
  }
}

 

void 
GEM_DomainBoundary::PopulateSurfaceArrays(const std::vector<double> &ic,
					  unsigned int ngnodes)
{
  if(_debug && _out)
    *_out  <<"GEM_DomainBoundary::PopulateSurfaceArrays(" << _id 
	   << "): Enter\n";
  map<unsigned int,unsigned int> s2v_imap; // idex map surface NC to volume NC 
  list<unsigned int> nodelist;
  std::vector<unsigned int>::iterator ci = _triconn.begin();
  unsigned int indy = 0;
  while(ci != _triconn.end()){
    if(*ci == 0){
      if(_out)
	*_out << "GEM_DomainBoundary::PopulateSurfaceArrays: Error: Found 0 in"
	      << " triangle conn. for boundary_id = " << _id  << "\n"
	      << " triangle/node (" << indy/3 << "," << (indy-((indy/3)*3)+1) 
	      << ")\n"
	      << " Number of total triangles = " << _triconn.size()/3 << "\n"
	      << " Number of ghost tri       = " << _ngtri << "\n";
      exit(1);
    }
    nodelist.push_back(*ci++);
    indy++;
  }
  ci = _quadconn.begin();
  while(ci != _quadconn.end()){
    assert(*ci != 0);
    nodelist.push_back(*ci++);
  }
  if(_debug && _out)
    *_out  << "GEM_DomainBoundary::PopulateSurfaceArrays(" << _id << "): "
	   << "nnodes before sort/unique: " << nodelist.size() << "\n";
  nodelist.sort();
  nodelist.unique();
  if(_debug && _out)
    *_out  << "GEM_DomainBoundary::PopulateSurfaceArrays)" << _id << "): "
	   << "nnodes after sort/unique: " << nodelist.size() << "\n";

  unsigned int nreal_volume_nodes = ic.size()/3 - ngnodes;
  unsigned int nreal_tris = _triconn.size()/3 - _ngtri;
  unsigned int nreal_quads = _quadconn.size()/3 - _ngquad;
  // Run a debugging check on the surface nodes
  if(_debug){
    bool fail = false;
    list<unsigned int>::iterator nli = nodelist.begin();
    while(nli != nodelist.end()){
      // Make sure NO ghost nodes live in real elements
      if(*nli > nreal_volume_nodes){
	unsigned int nrel = 0;
	// Check real triangles
	vector<unsigned int>::iterator ti = _triconn.begin();
	while(nrel < nreal_tris){
	  int ein = 0;
	  while(ein < 3){
	    if((unsigned int)*nli == (unsigned int)*ti++){
	      if(_out)
		*_out << "GEM_DomainBoundary::PopulateSurfaceArrays:"
		      << " Found ghost node, volume_id(" 
		      << *nli 
		      << ") in real tri(" << nrel+1 
		      << ") of surface " << _id << ". Aborting." << endl;
	      fail = true;
	    }
	    ein++;
	  }
	  nrel++;
	}
	// Check real quads
	nrel = 0;
	ti = _quadconn.begin();
	while(nrel < nreal_quads){
	  int ein = 0;
	  while(ein < 4){
	    if(*nli == *ti++){
	      if(_out)
		*_out << "GEM_DomainBoundary::PopulateSurfaceArrays:"
		      << " Found ghost node, volume_id(" 
		      << *nli << ") in real quad(" << nrel+1 
		      << ") of surface " << _id << ". Aborting." << endl;
	      fail = true;
	    }
	    ein++;
	  }
	  nrel++;
	}	
      }
      // It's a real node, warn if it's isolated
      else { 
	unsigned int nrel = 0;
	// Check real triangles
	vector<unsigned int>::iterator ti = _triconn.begin();
	bool found = false;
	while(nrel < nreal_tris && !found){
	  int ein = 0;
	  while(ein < 3 && !found){
	    if((unsigned int)*nli == (unsigned int)*ti++)
	      found = true;
	    ein++;
	  }
	  nrel++;
	}
	// Check real quads
	nrel = 0;
	ti = _quadconn.begin();
	while(nrel < nreal_quads && !found){
	  int ein = 0;
	  while(ein < 4 && !found){
	    if(*nli == *ti++)
	      found = true;
	    ein++;
	  }
	  nrel++;
	}
	if(!found && _out)
	  *_out << "GEM_DomainBoundary::PopulateSurfaceArrays: WARNING: "
		<< " Found isolated real node(" << *nli << ") on surface."
		<< std::endl;
      }
      nli++;
    }
    if(fail){
      if(_out)
	*_out << "GEM_DomainBoundary::PopulateSurfaceArrays: "
	      << "Aborting due to previous errors." << std::endl;
      exit(1);
    }
  }
  surface_coordinates.resize(3*nodelist.size());
  unsigned int nghosts = 0;
  unsigned int nvol_real = ic.size()/3 - ngnodes;
  unsigned int real_node = 0;
  unsigned int ghost_node = 0;
  list<unsigned int>::iterator nli = nodelist.begin();
  while(nli != nodelist.end())
    if(*nli++ > nvol_real)
      nghosts++;
  nli = nodelist.begin();
  unsigned int nreal_nodes = nodelist.size() - nghosts;
  if(_debug && _out){
    *_out  << "GEM_DomainBoundary::PopulateSurfaceArrays(" << _id 
	   << "): Nodes("
	   << nreal_nodes << "," << nghosts << ")\n";
  }
  while(nli != nodelist.end()){
    unsigned int node = *nli++; // the volume node index
    if(node > nvol_real){
      surface_coordinates[(nreal_nodes+ghost_node)*3]   = ic[(node -1)*3];
      surface_coordinates[(nreal_nodes+ghost_node)*3+1] = ic[(node -1)*3+1];
      surface_coordinates[(nreal_nodes+ghost_node)*3+2] = ic[(node -1)*3+2];
      ghost_node++;
      s2v_imap.insert(make_pair(node,nreal_nodes+ghost_node));
    }
    else{
      surface_coordinates[real_node*3]   = ic[(node-1)*3];
      surface_coordinates[real_node*3+1] = ic[(node-1)*3+1];
      surface_coordinates[real_node*3+2] = ic[(node-1)*3+2];
      real_node++;
      s2v_imap.insert(make_pair(node,real_node));
    }
  }
  surface_ngnodes = ghost_node;
  ci = _triconn.begin();
  surface_tri.resize(_triconn.size());
  surface_quad.resize(_quadconn.size());
  unsigned int ind = 0;
  while(ci != _triconn.end())
    surface_tri[ind++] = s2v_imap[*ci++];
  ind = 0;
  ci = _quadconn.begin();
  while(ci != _quadconn.end())
    surface_quad[ind++] = s2v_imap[*ci++];
  if(_debug && _out)
    *_out  << "GEM_DomainBoundary::PopulateSurfaceArrays(" << _id 
	   << "): Exit\n";
}

//#ifdef _ROCSTAR_X_

// utility for registering straightforward volume solution fields
void
GEM_Partition::Create_com_volsoln(const string &fname,
				  std::vector<double> &fvec,
				  unsigned int ncomp,const string &unit)
{
  string aname(volume_window+"."+fname);
  int ahndl = COM_get_dataitem_handle(aname);
  if(ahndl<=0)
    COM_new_dataitem(aname,'e',COM_DOUBLE,ncomp,unit.c_str());
  if(fvec.size() > 0)
    COM_set_array(aname,pane_id,&fvec[0],ncomp);
}

void
GEM_DomainBoundary::Create_com_surfsoln(const string &wname,
					const string &fname,
					std::vector<double> &fvec,
					unsigned int ncomp,const string &unit)
{
  string aname(wname+"."+fname);
  int ahndl = COM_get_dataitem_handle(aname);
  if(!ahndl)
    COM_new_dataitem(aname,'e',COM_DOUBLE,ncomp,unit.c_str());
  if(fvec.size() > 0)
    COM_set_array(aname,pane_id,&fvec[0],ncomp);
}

// Adds a section to the Pconn (shared, send, recv, etc)
void
AddPconnSection(std::vector<unsigned int> rpids,
		std::vector<std::vector<unsigned int> > &indices,
		unsigned int &section_size,
		std::vector<int> &pconn)
{
  section_size = 0;
  unsigned int nrp = rpids.size();
  pconn.push_back(nrp);
  section_size++;
  std::vector<unsigned int>::iterator rpi = rpids.begin();
  unsigned int rpindex = 0;
  while(rpi != rpids.end()){
    pconn.push_back(*rpi++);
    section_size++;
    unsigned int ne = indices[rpindex].size();
    pconn.push_back(ne);
    section_size += (ne+1);
    std::vector<unsigned int>::iterator ei = indices[rpindex].begin();
    while(ei != indices[rpindex].end())
      pconn.push_back(*ei++);
    rpindex++;
  }
}

//
// Build the pane connectivity array
// <nremotepanes> <rpaneid> <nshnode> <shnode indices>  <rpaneid> <nshnode> 
// <shnode indices> <...> 
// <nremotepanes> <rpaneid> <nnsend> <sendnode indices> <rpaneid> <nnsend> 
// <sendnode indices>  <...>
// <same for recv nodes>
// <same for send cells>
// <same for recv cells>
//
void
GEM_Partition::Create_com_pconn(std::vector<unsigned int> rpids,
				std::vector<std::vector<
				std::vector<unsigned int> > > &index_vectors,
				unsigned int &nreal,unsigned int &ng,
				std::vector<int> &pc)
{
  
  unsigned int section_size = 0;
  unsigned int i = 0;
  ng = 0;
  if(index_vectors[i].size() == 0 && _out){
    *_out 
      << "Roccom_create_pconn::Attempt to create empty pconn. Aborting.\n"; 
    exit(1);
  }
  pc.resize(0);
  AddPconnSection(rpids,index_vectors[i++],section_size,pc);
  nreal = section_size;
  while(i < 5 && index_vectors[i].size() != 0){
    section_size = 0;
    AddPconnSection(rpids,index_vectors[i++],section_size,pc);
    ng += section_size;
  }
}

// utility for registering built in Roccom volume connectivities
void
GEM_Partition::Register_com_volconn(const string &wname,int paneid,
				    unsigned int nel,unsigned int ngel,
				    std::vector<unsigned int> &conn,
				    unsigned int esize,bool ghost_part)
{
  if(nel == 0 || esize == 0) return;
  unsigned int nreal = conn.size()/esize - ngel;
  string mesh_type;
  switch(esize){
  case 4: // Tets
    mesh_type = ":T4";
    break;
  case 5:
    mesh_type = ":P5";
    break;
  case 6:
    mesh_type = ":P6";
    break;
  case 8:
    mesh_type = ":H8";
    break;
  default:
    if(_out)
      *_out << "Roccom_register_mesh::Error: Unknown mesh type, aborting.\n"; 
    exit(1);
  }
  mesh_type = wname + "." + mesh_type + ":";
  string entity;
  if(nreal > 0 && !ghost_part){
    entity = mesh_type + "real";
    COM_set_size(entity,paneid,nreal);
    COM_set_array(entity,paneid,&(conn[0]),esize);
  }
  if(ngel > 0 && ghost_part){
    entity = mesh_type + "virtual";
    COM_set_size(entity,paneid,ngel,ngel);
    COM_set_array(entity,paneid,&(conn[nreal*esize]),esize);
  } 
}


bool
GEM_DomainBoundary::Register_com_surfmesh(const string &wname)
{
  // Set the surface mesh entity sizes and register the arrays
  COM_set_size((wname+".nc"),pane_id,surface_coordinates.size()/3,
	       surface_ngnodes);
  COM_set_array((wname+".nc"),pane_id,&surface_coordinates[0],3);
  unsigned int nreal = _triconn.size()/3 - _ngtri;
  if(nreal > 0){
    COM_set_size((wname+".:t3:real"),pane_id,nreal);
    COM_set_array((wname+".:t3:real"),pane_id,&(surface_tri[0]),3);
  }
  nreal = _quadconn.size()/4 - _ngquad;
  if(nreal > 0){
    COM_set_size((wname+".:q4:real"),pane_id,nreal);
    COM_set_array((wname+".:q4:real"),pane_id,&(surface_quad[0]),4);
  }
  if(_ngtri > 0){
    nreal = _triconn.size()/3-_ngtri;
    COM_set_size((wname+".:t3:virtual"),pane_id,_ngtri,_ngtri);
    COM_set_array((wname+".:t3:virtual"),pane_id,&(surface_tri[nreal*3]),3);
  }
  if(_ngquad > 0){
    nreal = _quadconn.size()/4-_ngquad;
    COM_set_size((wname+".:q4:virtual"),pane_id,_ngquad,_ngquad);
    COM_set_array((wname+".:q4:virtual"),pane_id,&(surface_quad[nreal*4]),4);
  }
  return(true);
}

bool 
GEM_Partition::WindowInitDone()
{
  if(volume_window.empty() || surface_window.empty())
    return(false);
  COM_window_init_done(volume_window);
  COM_window_init_done(surface_window);
  return(true);
}

bool GEM_Partition::DestroyWindows()
{
  if(volume_window.empty() || surface_window.empty())
    return(false);
  COM_delete_window(volume_window);
  COM_delete_window(surface_window);
  //  volume_window.erase();
  //  surface_window.erase();
  return(true);
}
bool
GEM_Partition::ReadRocstar(const std::string &prefix,double t)
{
  return(false);
}

bool
GEM_Partition::WriteRocstar(const std::string &prefix,double t)
{
  string pre(prefix);
  if(pre.empty())
    pre = ".";
  if(volume_window.empty() || surface_window.empty())
    return false;
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id << ")::WriteRocstar: Writing volume window"
	  << " in " << prefix << ". CWD = " << TRAIL_CWD() << std::endl;
  if(!TRAIL_WriteWindow(volume_window,pre,volume_window,pre,t,_id,_comm))
    return false;
//   COM_LOAD_MODULE_STATIC_DYNAMIC(Rocout,"Rocout");
//   int OUT_set_option = COM_get_function_handle( "Rocout.set_option");
//   string rankstr("0");
//   COM_call_function( OUT_set_option, "rankwidth", rankstr.c_str());
  
//   // Build a filename and write the volume window
//   string timestring(TRAIL_TimeString(t));
//   ostringstream Ostr;
//   Ostr << pre << "/" << volume_window << "_"
//        << timestring << "_" << setw(5) << setfill('0')
//        << _id;
//   int whand = COM_get_function_handle("Rocout.write_dataitem");
//   int all = COM_get_dataitem_handle((volume_window+".all"));
//   if(_debug && _out)
//     *_out << "GEM_Partition(" << _id 
// 	  << ")::WriteRocstar: Writing volume window\n";
//   COM_call_function(whand,Ostr.str().c_str(),&all,volume_window.c_str(),
// 		    timestring.c_str());
//   // Write Rocin control file
//   std::vector<int> pane_ids;
//   string controlfilename;
//   COM_get_panes(volume_window.c_str(),pane_ids);
//   ofstream Ouf;
//   controlfilename = Ostr.str() + "_in.txt";
//   Ouf.open(controlfilename.c_str());
//   Ouf << "@Proc: " << _id - 1 << endl
//       << "@Files: " << volume_window << "_" << timestring << "_" 
//       << setw(5) << setfill('0') << _id << ".hdf" << endl;
//   Ouf.clear();
//   Ouf << "@Panes: " << pane_ids[0] << endl;
//   Ouf.close();

//   Ostr.clear();
//   pane_ids.resize(0);
//   Ostr.str("");
//   Ostr << pre << "/" << surface_window << "_"
//        << timestring << "_" << setw(5) << setfill('0')
//        << _id;
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id  
	  << ")::WriteRocstar: Writing surface window\n";
  if(!TRAIL_WriteWindow(surface_window,pre,surface_window,pre,t,_id,_comm))
    return false;
//   all = COM_get_dataitem_handle((surface_window+".all"));
//   COM_call_function(whand,Ostr.str().c_str(),&all,surface_window.c_str(),
// 		    timestring.c_str());
//   // Write Rocin control file
//   controlfilename = Ostr.str() + "_in.txt";
//   COM_get_panes(surface_window.c_str(),pane_ids);
//   Ouf.open(controlfilename.c_str());
//   Ouf << "@Proc: " << _id - 1 << endl
//       << "@Files: " << surface_window << "_" << timestring << "_" 
//       << setw(5) << setfill('0') << _id << ".hdf" << endl;
//   Ouf.clear();
//   Ouf << "@Panes: ";
//   std::vector<int>::iterator pii = pane_ids.begin();
//   while(pii != pane_ids.end())
//     Ouf << *pii++ << " ";
//   Ouf << endl;
//   Ouf.close();
//   if(_debug && _out)
//     *_out << "GEM_Partition(" << _id << ")::WriteRocstar: Unloading Rocout\n";
//   COM_UNLOAD_MODULE_STATIC_DYNAMIC(Rocout,"Rocout");
  
  return(true);
}

// This function is an all-in-one window creation/registration utility
// Note that your data structures need to be consistent after this is 
// called until you have finished using the Windows this function
// has created.
bool 
GEM_Partition::InitRoccomWindows(const string &wname)
{
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id 
	  << ")::InitRoccomWindows: Creating windows" << endl;
  volume_window = wname+"_vol";
  surface_window = wname+"_surf";
  COM_new_window(volume_window);
  COM_new_window(surface_window);
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id 
	  << ")::InitRoccomWindows: Populate Volume Window" 
	  << endl;
  if(!PopulateVolumeWindow(volume_window))
    return(false);
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id 
	  << ")::InitRoccomWindows: Populate Surface Window" 
	  << endl;
  if(!PopulateSurfaceWindow(surface_window))
    return(false);
  return(true);
}

bool 
GEM_Partition::CreatePconn(const string &wname)
{
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id << ")::CreatePconn: enter" << endl;
  pconn_nghost = 0;
  unsigned int nrp = _pb.size(); // number of remote partitions
  std::vector<unsigned int> rpid_v; // remote pane id's
  std::vector<std::vector<std::vector<unsigned int> > > indices_v(5);
  rpid_v.resize(nrp);
  indices_v[0].resize(nrp); // shared nodes
  indices_v[1].resize(nrp); // sent nodes
  indices_v[2].resize(nrp); // recv nodes
  indices_v[3].resize(nrp); // send elements
  indices_v[4].resize(nrp); // recv elements
  unsigned int rpin = 0;
  while(rpin < nrp){
    GEM_PartitionBoundary &fb = _pb[rpin];
    unsigned int rpid = fb._rpart * 100 + 1;
    rpid_v[rpin] = rpid;
    indices_v[0][rpin] = fb._sharenodes;
    indices_v[1][rpin] = fb._sendnodes;
    indices_v[2][rpin] = fb._recvnodes;
    indices_v[3][rpin] = fb._sendcells;
    indices_v[4][rpin] = fb._recvcells;
    rpin++;
  }
  unsigned int nreal = 0;
  Create_com_pconn(rpid_v,indices_v,nreal,pconn_nghost,pconn);
  assert(pconn.size() - pconn_nghost  == nreal);
  COM_set_size((wname+".pconn"),pane_id,pconn.size(),pconn_nghost);
  COM_set_array((wname+".pconn"),pane_id,&pconn[0],1);
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id << ")::CreatePconn: exit" << endl;
  return(true);
}
  
bool 
GEM_Partition::PopulateVolumeWindow(const string &wname)
{
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id 
	  << ")::PopulateVolumeWindow: enter" << endl;
  pane_id = _id * 100 + 1; 
  // Set the volume mesh entity sizes and register the arrays
  unsigned int nnodes = _nc.size()/3;
  COM_set_size((wname+".nc"),pane_id,nnodes,_ngnodes);
  COM_set_array((wname+".nc"),pane_id,&_nc[0],3);
  unsigned int el_type = 0;
  while(el_type < 4){
    unsigned int nreal;
    switch(_cell_ordering[el_type]){
    case 1:
      nreal = _tetconn.size()/4 - _ngtet;
      if(nreal > 0)
	Register_com_volconn(wname,pane_id,nreal,_ngtet,_tetconn,4,false);
      break;
    case 2:
      nreal = _pyrconn.size()/5 - _ngpyr;
      if(nreal > 0)
	Register_com_volconn(wname,pane_id,nreal,_ngpyr,_pyrconn,5,false);
      break;
    case 3:
      nreal = _prisconn.size()/6 - _ngpris;
      if(nreal > 0)
	Register_com_volconn(wname,pane_id,nreal,_ngpris,_prisconn,6,false);
      break;
    case 4:
      nreal = _hexconn.size()/8 - _nghex;
      if(nreal > 0)
	Register_com_volconn(wname,pane_id,nreal,_nghex,_hexconn,8,false);
      break;
    }
    el_type++;
  }
  el_type = 0;
  while(el_type < 4){
    unsigned int nreal;
    switch(_cell_ordering[el_type]){
    case 1:
      nreal = _tetconn.size()/4 - _ngtet;
      if(_ngtet > 0)
	Register_com_volconn(wname,pane_id,_ngtet,_ngtet,_tetconn,4,true);
      break;
    case 2:
      nreal = _pyrconn.size()/5 - _ngpyr;
      if(_ngpyr > 0)
	Register_com_volconn(wname,pane_id,_ngpyr,_ngpyr,_pyrconn,5,true);
      break;
    case 3:
      nreal = _prisconn.size()/6 - _ngpris;
      if(_ngpris > 0)
	Register_com_volconn(wname,pane_id,_ngpris,_ngpris,_prisconn,6,true);
      break;
    case 4:
      nreal = _hexconn.size()/8 - _nghex;
      if(_nghex > 0)
	Register_com_volconn(wname,pane_id,_nghex,_nghex,_hexconn,8,true);
      break;
    }
    el_type++;
  }
  CreatePconn(wname);
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id 
	  << ")::PopulateVolumeWindow: exit" << endl;
  return(true);
}

bool 
GEM_Partition::PopulateSurfaceWindow(const string &wname)
{
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id 
	  << ")::PopulateSurfaceWindow: enter" << endl;
  // For each domain boundary
  int npatches = _db.size();
  int patch = 0;
  while(patch < npatches){
    GEM_DomainBoundary &fp = _db[patch];
    fp.pane_id = _id * 100 + (patch+1) + 1;
    fp.PopulateSurfaceArrays(_nc,_ngnodes);
    fp.Register_com_surfmesh(wname);
    patch++;
  }
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id 
	  << ")::PopulateSurfaceWindow: exit" << endl;
  return(true);
}

//#endif

// Utilities
bool 
flip_elements(std::vector<unsigned int> &conn,unsigned int es)
{
  unsigned int nel = conn.size()/es;
  unsigned int temp;
  unsigned int el = 0;
  switch(es){
  case 2:
  case 3:
  case 4:
    while(el < nel){
      temp = conn[el];
      conn[el] = conn[el+1];
      el = el+1;
      conn[el] = temp;
      el = el+1;
      if(es > 2)
	el++;
      if(es > 3)
	el++;
    }
    return(true);
    break;
  default:
    return(false);
    break;
  }
  return(false);
}

// Return the element type and index for a given cell number
pair<unsigned int,unsigned int>
GEM_Partition::Cell2Elem(unsigned int cell)
{
  unsigned int elem_type = 0;
  unsigned int offset = 0;
  unsigned int nreal_elems;
  while(elem_type < 4){
    switch(_cell_ordering[elem_type]){
    case 1:
      nreal_elems = _tetconn.size()/4 - _ngtet;
      if(cell <= (offset+nreal_elems))
	return(make_pair(1,cell-offset));
      else
	offset += nreal_elems;
      break;
    case 2:
      nreal_elems = _pyrconn.size()/5 - _ngpyr;
      if(cell <= (offset+nreal_elems))
	return(make_pair(2,cell-offset));
      else
	offset += nreal_elems;
      break;
    case 3:
      nreal_elems = _prisconn.size()/6 - _ngpris;
      if(cell <= (offset+nreal_elems))
	return(make_pair(3,cell-offset));
      else
	offset += nreal_elems;
      break;
    case 4:
      nreal_elems = _hexconn.size()/8 - _nghex;
      if(cell <= (offset+nreal_elems))
	return(make_pair(4,cell-offset));
      else
	offset += nreal_elems;
      break;
    }
    elem_type++;
  }
  elem_type = 0;
  while(elem_type < 4){
    switch(_cell_ordering[elem_type]){
    case 1:
      nreal_elems = _tetconn.size()/4 - _ngtet;
      if(cell <= (offset+_ngtet))
	return(make_pair(1,cell-offset+nreal_elems));
      else
	offset += _ngtet;
      break;
    case 2:
      nreal_elems = _pyrconn.size()/5 - _ngpyr;
      if(cell <= (offset+_ngpyr))
	return(make_pair(2,cell-offset+nreal_elems));
      else
	offset += _ngpyr;
      break;
    case 3:
      nreal_elems = _prisconn.size()/6 - _ngpris;
      if(cell <= (offset+_ngpris))
	return(make_pair(3,cell-offset+nreal_elems));
      else
	offset += _ngpris;
      break;
    case 4:
      nreal_elems = _hexconn.size()/8 - _nghex;
      if(cell <= (offset+_nghex))
	return(make_pair(4,cell-offset+nreal_elems));
      else
	offset += _nghex;
      break;
    }
    elem_type++;
  }
  if(_out)
    *_out << "GEM_Partition(" << _id 
	  << ")::Cell2Elem: Fatal error - Could not find cell " 
	  << cell << ", dying.\n";
  exit(1);
}

// Return cell number for a given element type and index
unsigned int
GEM_Partition::Elem2Cell(pair<unsigned int,unsigned int> ti)
{
  unsigned int elem_type = 0;
  unsigned int offset = 0;
  unsigned int nreal_elem;
  while(elem_type < 4){
    switch(_cell_ordering[elem_type]){
    case 1:
      nreal_elem = _tetconn.size()/4 - _ngtet;
      if(ti.first != 1 || ti.second > nreal_elem)
	offset += nreal_elem;
      else{
	return(offset+ti.second);
      }
      break;
    case 2:
      nreal_elem = _pyrconn.size()/5 - _ngpyr;
      if(ti.first != 2 || ti.second > nreal_elem)
	offset += nreal_elem;
      else{
	return(offset+ti.second);
      }
      break;
    case 3:
      nreal_elem = _prisconn.size()/6 - _ngpris;
      if(ti.first != 3 || ti.second > nreal_elem)
	offset += nreal_elem;
      else{
	return(offset+ti.second);
      }
      break;
    case 4:
      nreal_elem = _hexconn.size()/8 - _nghex;
      if(ti.first != 4 || ti.second > nreal_elem)
	offset += nreal_elem;
      else{
	return(offset+ti.second);
      }
      break;
    }
    elem_type++;
  }
  elem_type = 0;
  while(elem_type < 4){
    switch(_cell_ordering[elem_type]){
    case 1:
      nreal_elem = _tetconn.size()/4 - _ngtet;
      if(ti.first != 1)
	offset += _ngtet;
      else{
	return(offset+ti.second-nreal_elem);
      }
      break;
    case 2:
      nreal_elem = _pyrconn.size()/5 - _ngpyr;
      if(ti.first != 2)
	offset += _ngpyr;
      else{
	return(offset+ti.second-nreal_elem);
      }
      break;
    case 3:
      nreal_elem = _prisconn.size()/6 - _ngpris;
      if(ti.first != 3)
	offset += _ngpris;
      else{
	return(offset+ti.second-nreal_elem);
      }
      break;
    case 4:
      nreal_elem = _hexconn.size()/8 - _nghex;
      if(ti.first != 4)
	offset += _nghex;
      else{
	return(offset+ti.second-nreal_elem);
      }
      break;
    }
    elem_type++;
  }
  if(_out)
    *_out << "GEM_Partition(" << _id 
	  << ")::Elem2Cell: Fatal error.  Could not find element "
	  << "(" << ti.first << "," << ti.second << "). Dying.\n";
  exit(1);
}


bool
GEM_Partition::SetSolverDataBlock(const string &wname,double *cell_data,
				  int nval_cells,double *node_data,
				  int nval_nodes)
{
  _solver_data._string_data.push_back(wname);
  _solver_data._field_data.resize(2);
  unsigned int ncells = _tetconn.size()/4 + _hexconn.size()/8 +
    _prisconn.size()/6 + _pyrconn.size()/5;
  unsigned int nnodes = _nc.size()/3;
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id << ")::SetSolverDataBlock: "
	  << "Receiving data for " << nval_cells << " doubles on " << ncells 
	  << " cells and " << nval_nodes <<  " doubles on " << nnodes 
	  << " nodes." << std::endl; 
  _solver_data._field_data[0].resize(nval_cells*ncells);
  _solver_data._field_data[1].resize(nval_nodes*nnodes);
  _solver_data._stride_field.resize(2);
  _solver_data._stride_field[0] = nval_cells;
  _solver_data._stride_field[1] = nval_nodes;
  memcpy(&_solver_data._field_data[0][0],cell_data,
	 sizeof(double)*nval_cells*ncells);
  memcpy(&_solver_data._field_data[1][0],node_data,
	 sizeof(double)*nval_nodes*nnodes);
  return(true);
}

bool
GEM_Partition::AddSolverDataBlock(const string &wname,double *cell_data,
				  int nval_cells,double *node_data,
				  int nval_nodes)
{
  _solver_data._string_data.push_back(wname);
  unsigned int current_size = _solver_data._field_data.size();
  if(current_size == 0)
    _solver_data._field_data.resize(2);
  else{
    std::vector<double> temp;
    _solver_data._field_data.push_back(temp);
    _solver_data._field_data.push_back(temp);
    _solver_data._stride_field.push_back(0);
    _solver_data._stride_field.push_back(0);
  }
  unsigned int ncells = _tetconn.size()/4 + _hexconn.size()/8 +
    _prisconn.size()/6 + _pyrconn.size()/5;
  unsigned int nnodes = _nc.size()/3;
  _solver_data._field_data[current_size].resize(nval_cells*ncells);
  _solver_data._field_data[current_size+1].resize(nval_nodes*nnodes);
  _solver_data._stride_field.resize(2);
  _solver_data._stride_field[current_size] = nval_cells;
  _solver_data._stride_field[current_size+1] = nval_nodes;
  memcpy(&_solver_data._field_data[current_size][0],cell_data,
	 sizeof(double)*nval_cells*ncells);
  memcpy(&_solver_data._field_data[current_size+1][0],node_data,
	 sizeof(double)*nval_nodes*nnodes);
  return(true);
}

unsigned int
GEM_DomainBoundary::NNodes()
{
  if(_nnodes <= 0){
    list<unsigned int> snlist;
    vector<unsigned int>::iterator ci = _triconn.begin();
    while(ci != _triconn.end())
      snlist.push_back(*ci++);
    ci = _quadconn.begin();
    while(ci != _quadconn.end())
      snlist.push_back(*ci++);
    snlist.sort();
    snlist.unique();
    _nnodes = snlist.size();
  }
  return(_nnodes);
}

bool
GEM_DomainBoundary::SetSolverDataBlock(const string &wname,
				       double *cell_data,
				       int nval_cells, 
				       double *node_data,int nval_nodes)
{
  _solver_data._string_data.push_back(wname);
  _solver_data._field_data.resize(2);
  unsigned int ncells = _triconn.size()/3 + _quadconn.size()/4;
  unsigned int nnodes = NNodes();
  if(_debug && _out)
    *_out << "GEM_DomainBoundary(" << _id << ")::SetSolverDataBlock: "
	  << "Receiving data for " << nval_cells << " doubles on " << ncells 
	  << " cells and " << nval_nodes <<  " doubles on " << nnodes 
	  << " nodes." << std::endl; 
  _solver_data._field_data[0].resize(nval_cells*ncells);
  _solver_data._field_data[1].resize(nval_nodes*nnodes);
  _solver_data._stride_field.resize(2);
  _solver_data._stride_field[0] = nval_cells;
  _solver_data._stride_field[1] = nval_nodes;
  memcpy(&_solver_data._field_data[0][0],cell_data,
	 sizeof(double)*nval_cells*ncells);
  memcpy(&_solver_data._field_data[1][0],node_data,
	 sizeof(double)*nval_nodes*nnodes);
  return(true);
}


bool
GEM_Partition::validate_comm_list(int ncsend,int ncrecv,int *csend,int *crecv)
{
  int index = 0;
  int nreal_cell = _tetconn.size()/4 + _prisconn.size()/6 +
    _pyrconn.size()/5 + _hexconn.size()/8 - (_ngtet + _ngpris +
					     _ngpyr + _nghex);
  bool rval = true;
  while(index < ncsend){
    int ind = index++;
    if(!(csend[ind] <= nreal_cell)){
      if(_out)
	*_out << "SEND CELL " << index << " is a ghost cell!!" << endl;
      rval = false;
    }
    if(!(csend[ind] > 0)){
      if(_out)
	*_out << "SEND CELL " << index << " is zero or negative!" << endl;
      rval = false;
    }
  }
  index = 0;
  list<int> recvcell_list;
  while(index < ncrecv) {
    int ind = index++;
    if(!(crecv[ind] > nreal_cell)){
      if(_out)
	*_out << "RECV CELL " << index << " is a real cell!!" << endl;
      rval = false;
    }
    if(!(crecv[ind] > 0)){
      if(_out)
	*_out << "RECV CELL " << index << " is zero or negative!" << endl;
      rval = false;
    }
    bool duped = false;
    list<int>::iterator rci = recvcell_list.begin();
    while(rci != recvcell_list.end() && !duped){
      if(crecv[ind] == *rci++){
	if(_out)
	  *_out << "RECV CELL " << index 
		<< " is duplicated in the receive list!" 
		<< endl;
	duped = true;
      }
    }
    if(!duped)
      recvcell_list.push_back(crecv[ind]);
  }
  return(rval);
}

void 
GEM_Partition::AddParitionBoundary(int rpid,int nnshare, int nnsend,
				   int nnrecv,int ncsend,int ncrecv,
				   int *nshared,int *nsend,int *nrecv,
				   int *csend,int *crecv)
{
  assert(rpid > 0);
  if(_debug && _out)
    *_out << "GEM_Mesh(" << _id << ")::AddPartitionBoundary: "
	  << "Adding Border with"
	  << " partition " << rpid << "." << endl;
  GEM_PartitionBoundary new_pb;
  new_pb._out = _out;
  new_pb._debug = _debug;
  if(!validate_comm_list(ncsend,ncrecv,csend,crecv)){
    if(_out)
      *_out << "GEM_Mesh(" << _id << ")::AddPartitionBoundary"
	    << ": Validation of "
	    << "communication arrays failed, aborting." << endl;
    exit(-1);
  }
  new_pb.populate(rpid,nnshare,nnsend,nnrecv,ncsend,ncrecv,
		  nshared,nsend,nrecv,csend,crecv);
  _pb.push_back(new_pb);
}

void 
GEM_Partition::AddDomainBoundary(int db_id,int ntri, int ngtri, int *tris,
				 int nquad,int ngquad, int *quads)
{
  assert(ntri >= ngtri && nquad >= ngquad);
  if(_debug && _out)
    *_out << "GEM_Mesh(" << _id << ")::AddDomainBoundary: "
	  << "Adding domain boundary with"
	  << " id " << db_id << "." << endl;
  GEM_DomainBoundary new_db;
  new_db._id = db_id;
  new_db._ngtri = ngtri;
  new_db._ngquad = ngquad;
  int indy = 0;
  new_db._triconn.resize(3*ntri);
  new_db._quadconn.resize(4*nquad);
  new_db._out = _out;
  new_db._debug = _debug;
  while(indy < 3*ntri){
    assert(tris[indy] != 0);
    new_db._triconn[indy] = tris[indy];
    indy++;
  }
  indy = 0;
  while(indy < 4*nquad){
    assert(quads[indy] != 0);
    new_db._quadconn[indy] = quads[indy];
    indy++;
  }
  new_db._nnodes = new_db.NNodes();
  unsigned int csize = _db.size();
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id << ")::AddDomainBoundary: "
	  << "DomainBoundary " << csize << " has (nodes,tri,gtri,quad,gquad)" 
	  << " = (" << new_db._nnodes << "," << ntri << "," << ngtri << "," 
	  << nquad << "," << ngquad << ")" << std::endl;
  _db.push_back(new_db);
}

bool 
GEM_Partition::debug(bool s)
{
  _debug = s;
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id << ")::Debugging turned on\n";
  return(_debug);
}

// Simple interface utilities
void 
GEM_Partition::SetNodalCoordinates(double *data,int nn,int ng)
{
  if(_debug && _out){
    ostringstream Ostr;
    Ostr << "GEM_Partition(" << _id << ")::SetNodalCoordinates "
	 << " total nodes = " << nn << ", " << ng << " ghosts.\n";
    *_out << Ostr.str();
  }
  assert(nn >= ng);
  _ngnodes = ng;
  _nc.resize(3*nn);
  memcpy((void *)&_nc[0],(void *)data,3*nn*sizeof(double));
}

void 
GEM_Partition::SetVolumeElements(int *data,int ncells,int ng,int npe)
{
  if(_debug && _out){
    *_out << "GEM_Partition(" << _id << ")::SetVolumeElements "
	  << " total cells = " << ncells << " of size " << npe
	  << " of which " << ng << " are ghosts." << endl;
  }
  unsigned int datasize = ncells*npe*sizeof(int);
  assert(ncells >= ng);
  int ndatp = npe*ncells;
  for(int ndat =0;ndat < ndatp;ndat++)
    assert(data[ndat] != 0);
  void *dest;
  switch(npe){
  case 4:
    _tetconn.resize(4*ncells);
    dest = &_tetconn[0];
    _ngtet = ng;
    break;
  case 5:
    _pyrconn.resize(5*ncells);
    dest = &_pyrconn[0];
    _ngpyr = ng;
    break;
  case 6:
    _prisconn.resize(6*ncells);
    dest = &_prisconn[0];
    _ngpris = ng;
    break;
  case 8:
    _hexconn.resize(8*ncells);
    dest = &_hexconn[0];
    _nghex = ng;
    break;
  default:
    if(_out)
      *_out << "GEM_Partition::Unknown volume element type.  Aborting." 
	    << endl;
    exit(1);
  }
  memcpy(dest,data,datasize);
}

// bool
// GEM_Partition::ValidateMesh()
// {
//   bool error = false;
//   //  unsigned int nvolnodes = _nc.size()/3;
//   //  unsigned int nvolnodes_real = nvolnodes - _ngnodes;
//   unsigned int ntet  = _tetconn.size()/4;
//   unsigned int ntet_real = ntet - _ngtet;
//   unsigned int npyr  = _pyrconn.size()/5;
//   unsigned int npyr_real = npyr - _ngpyr;
//   unsigned int npris = _prisconn.size()/6;
//   unsigned int npris_real = npris - _ngpris;
//   unsigned int nhex  = _hexconn.size()/8;
//   unsigned int nhex_real = nhex - _nghex;
//   unsigned int ncells = nhex + npris + npyr + ntet;
//   //  unsigned int ngcells = ncells - 
//   //  (ntet_real+nhex_real+npyr_real+npris_real);
//   // Make sure no ghost nodes live in a real element
//   // Check for isolated mesh entities
//   //  - real nodes not belonging to any real element
//   //  - dangling faces and edges
//   return (!error);
// }
bool 
GEM_Partition::PopulatePartitionBoundaries(std::vector<GEM_PartitionBoundary> 
					   &pb)
{
  unsigned int nborders = pb.size();
  _pb.resize(nborders);
  //    if(_solver_data._int_data.empty())
  //      _solver_data._int_data.resize(2);
  //    _solver_data._int_data[0].resize(nborders);
  unsigned int border = 0;
  while(border < nborders){
    _pb[border]._rpart      = pb[border]._rpart;
    _pb[border]._sendcells  = pb[border]._sendcells;
    _pb[border]._recvcells  = pb[border]._recvcells;
    _pb[border]._sharenodes = pb[border]._sharenodes;
    _pb[border]._sendnodes  = pb[border]._sendnodes;
    _pb[border]._recvnodes  = pb[border]._recvnodes;
    //      _pb[border]._rbid       = pb[border]._rbid;
    _pb[border]._out        = pb[border]._out;
    border++;
  }
  return(true);
}

// Maps the domain boundaries on one mesh representation to domain
// boundaries on another.  Useful when mapping BC's from GridGen to
// application dependent BC's as often multiple BC's from GridGen 
// map to a common surface with the same application specific BC.
void GEM_Partition::MapDomainBoundaries(map<unsigned int,unsigned int> &bcmap)
{
  std::vector<GEM_DomainBoundary> indb(_db);
  std::vector<GEM_DomainBoundary> &outdb = _db;
  // Determine how many unique target BC's we have on this 
  // source DB.
  list<unsigned int> local_patches;
  unsigned int partpatch = 0;
  unsigned int npartpatch = indb.size();
  while(partpatch < npartpatch)
    local_patches.push_back(bcmap[indb[partpatch++]._id]);
  local_patches.sort();
  local_patches.unique();
  unsigned int nlocal_patches = local_patches.size();
  if(_debug && _out){
    *_out << "GEM_Partition(" << _id 
	  << ")::MapDomainBoundaries: Local patches: " << endl
	  << "GEM_Partition(" << _id << ")::MapDomainBoundaries: ";
    list<unsigned int>::iterator pli = local_patches.begin();
    while(pli != local_patches.end())
      *_out << *pli++ << " ";
    *_out << endl;
  }
  // resize the storage for the local db's to the number of unique 
  // bc types on this partition
  outdb.resize(nlocal_patches);
  std::vector< list<unsigned int> > ppatch_list;
  ppatch_list.resize(nlocal_patches);
  
  // For every bc type on the source db, assign the target bc type to the 
  // local storage array.  Also construct an index mapping so we
  // can map bctype ----> local_patch_storage_array_index
  unsigned int local_patch = 0;
  list<unsigned int>::iterator li = local_patches.begin();
  while(li != local_patches.end()){  
    GEM_DomainBoundary &fp = outdb[local_patch];
    fp._id     = *li;
    fp._ngtri  = 0;
    fp._ngquad = 0;
    if(_debug) 
      fp.debug();
    unsigned int local_tri_size = 0;
    unsigned int local_quad_size = 0;
    // Now - loop through every source bc and find out which ones 
    // contribute to this particular target bc so we can determine the total
    // size of the connectivities for pre-allocation.
    partpatch = 0;
    while(partpatch < npartpatch){
      const GEM_DomainBoundary &pp = indb[partpatch];
      unsigned int ggbpid = pp._id;
      unsigned int tdbid = bcmap[ggbpid];
      if(tdbid == *li){
	local_tri_size  += pp._triconn.size()/3;
	fp._ngtri  += pp._ngtri;
	local_quad_size  += pp._quadconn.size()/4;
	fp._ngquad += pp._ngquad;
	ppatch_list[local_patch].push_back(partpatch);
      }
      partpatch++;
    }
    fp._triconn.resize(3*local_tri_size);
    fp._quadconn.resize(4*local_quad_size);
    li++;
    local_patch++;
  }
  
  // Step through all the target db's and check the list that tells us
  // which source db's to get our surface element connectivities from.
  // This is the second pass to actually populate the connectivities. 
  local_patch = 0;
  while(local_patch < nlocal_patches){
    GEM_DomainBoundary &fp = outdb[local_patch];
    fp._out = indb[0]._out;
    unsigned int tri       = 0;
    unsigned int gtri      = 0;
    unsigned int quad      = 0;
    unsigned int gquad     = 0;
    unsigned int realtri   = fp._triconn.size()/3 - fp._ngtri;
    unsigned int realquads = fp._quadconn.size()/4 - fp._ngquad;
    li = ppatch_list[local_patch].begin();
    while(li != ppatch_list[local_patch].end()){
      const GEM_DomainBoundary &pp = indb[*li++];
      unsigned int pptri = 0;
      unsigned int ppquad = 0;
      unsigned int ntri = pp._triconn.size()/3;
      unsigned int nrealtri = ntri - pp._ngtri;
      unsigned int nquad = pp._quadconn.size()/4;
      unsigned int nrealquad = nquad - pp._ngquad;
      while(ppquad < nquad){
	if(ppquad < nrealquad){
	  fp._quadconn[4*quad]   = pp._quadconn[4*ppquad];
	  fp._quadconn[4*quad+1] = pp._quadconn[4*ppquad+1];
	  fp._quadconn[4*quad+2] = pp._quadconn[4*ppquad+2];
	  fp._quadconn[4*quad+3] = pp._quadconn[4*ppquad+3];
	  quad++;
	  ppquad++;
	}
	else{
	  fp._quadconn[(4*realquads)+(4*gquad)]   = pp._quadconn[4*ppquad];
	  fp._quadconn[(4*realquads)+(4*gquad)+1] = pp._quadconn[4*ppquad+1];
	  fp._quadconn[(4*realquads)+(4*gquad)+2] = pp._quadconn[4*ppquad+2];
	  fp._quadconn[(4*realquads)+(4*gquad)+3] = pp._quadconn[4*ppquad+3];
	  gquad++;
	  ppquad++;
	}
      }
      while(pptri < ntri){
	if(pptri < nrealtri){
	  fp._triconn[3*tri]   = pp._triconn[3*pptri];
	  fp._triconn[3*tri+1] = pp._triconn[3*pptri+1];
	  fp._triconn[3*tri+2] = pp._triconn[3*pptri+2];
	  tri++;
	  pptri++;
	}
	else{
	  fp._triconn[(3*realtri)+(3*gtri)]   = pp._triconn[3*pptri];
	  fp._triconn[(3*realtri)+(3*gtri)+1] = pp._triconn[3*pptri+1];
	  fp._triconn[(3*realtri)+(3*gtri)+2] = pp._triconn[3*pptri+2];
	  gtri++;
	  pptri++;
	}
      }
    }
    local_patch++;
  }
  return;
}
// The procedure for using this function is to first copy the source
// representation's partition boundary structures.  Then arrange the
// _cell_ordering[] array appropriately for your own cell mapping. Then
// call this function, passing in the source's Partition.  The code below
// is quite obvious.
//
// If you are calling this function, then your own cells are in the "wrong"
// order already.  This function will rearrange the cell mapping from the 
// source partition (sp) representation to the local one defined by the
// _cell_ordering[] array.
//template<typename PB,typename TP>
void 
GEM_Partition::ResolveCellMapping(GEM_Partition &sp)
{
  if(_debug && _out)
    *_out << "GEM_Partition::ResolveCellMapping: enter";
  unsigned int npb = _pb.size();
  unsigned int p = 0;
  while(p < npb){
    GEM_PartitionBoundary &pbi = _pb[p++];
    unsigned int nsend = pbi._sendcells.size();
    unsigned int cell = 0;
    while(cell < nsend){
      pbi._sendcells[cell] = Elem2Cell(sp.Cell2Elem(pbi._sendcells[cell]));
      cell++;
    }
    unsigned int nrecv = pbi._recvcells.size();
    cell = 0;
    while(cell < nrecv){
      pbi._recvcells[cell] = Elem2Cell(sp.Cell2Elem(pbi._recvcells[cell]));
      cell++;
    }
  }
  if(_debug && _out)
    *_out << "GEM_Partition(" << _id << ")::ResolveCellMapping: exit";
} 
//   void AddData(const string &name,int *data,int stride,int nitems)
//   {
//     // Search through existing data to see if we already have a dataset
//     // with this name, if so, use it - if not, create one.
//     std::vector<GEM_UserData>::iterator gdi = _data.begin();
//     while(gdi != _data.end() && gdi->_name != name)
//       gdi++;
//     if(gdi == _data.end()){
//       GEM_UserData newdata;
//       newdata._name = name;
//       newdata._int_data.resize(1);
//       newdata._int_data[0].resize(stride*nitems);
//       memcpy(&newdata._int_data[0],data,sizeof(int)*stride*nitems);
//       newdata._stride_int.resize(1);
//       newdata._stride_int[0] = stride;
//       _data.push_back(newdata);
//     }
//     else {
//       std::vector<int> newdata;
//       newdata.resize(stride*nitems);
//       memcpy(&newdata[0],data,sizeof(int)*stride*nitems);
//       gdi->_stride_int.push_back(stride);
//       gdi->_int_data.push_back(newdata);
//     }
//   };
//   void AddData(const string &name,double *data,int stride,int nitems)
//   {
//     // Search through existing data to see if we already have a dataset
//     // with this name, if so, use it - if not, create one.
//     std::vector<GEM_UserData>::iterator gdi = _data.begin();
//     while(gdi != _data.end() && gdi->_name != name)
//       gdi++;
//     if(gdi == _data.end()){
//       GEM_UserData newdata;
//       newdata._name = name;
//       newdata._field_data.resize(1);
//       newdata._field_data[0].resize(stride*nitems);
//       memcpy(&newdata._field_data[0],data,sizeof(double)*stride*nitems);
//       newdata._stride_field.resize(1);
//       newdata._stride_field[0] = stride;
//       _data.push_back(newdata);
//     }
//     else {
//       std::vector<double> newdata;
//       newdata.resize(stride*nitems);
//       memcpy(&newdata[0],data,sizeof(double)*stride*nitems);
//       gdi->_stride_field.push_back(stride);
//       gdi->_field_data.push_back(newdata);
//     }
//   };












