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
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <map>
#include <list>
#include <cmath>
#include <cassert>

using namespace std;

#include "roccom.h"

#include "GEM.H"
#include "FluRegion.H"

unsigned int 
FluRegion::nelem(){
  return(_ntet + _npyr + _npris + _nhex);
}

void 
FluRegion::report()
{
  if(_out){
    *_out  << "********************************************************" 
	   << "\n"
	   << "Region ID: " << _id << "\n"
	   << "Number of Nodes: (" << _nnodes << "," << _ngnodes << ")" << "\n"
	   << "Number of Elements: " << nelem() << "\n"
	   << "Tets:     (" << _ntet  << "," << _ngtet  << ")" << "\n"
	   << "Pyramids: (" << _npyr  << "," << _ngpyr  << ")" << "\n"
	   << "Prisms:   (" << _npris << "," << _ngpris << ")" << "\n"
	   << "Hexes:    (" << _nhex  << "," << _nghex  << ")" << "\n"
	   << "Region Borders: " << "\n"
	   << "========================================================" 
	   << "\n" ;
    unsigned int border = 0;
    while(border < _borders.size()){
      *_out  << " Local Border #" << border+1 << "\n" ;
      FluBorder *fb = &_borders[border++];
      *_out << "-----------------------------------------------------" << "\n"
	    << "   Region/Border: (" << fb->_rpart << "," << fb->_rbid 
	    << ")" << "\n"
	    << "   Nodes: (" << fb->_sharenodes.size() << "," 
	    << fb->_sendnodes.size() 
	    << ","  << fb->_recvnodes.size() << ")" << "\n"
	    << "   Cells: (" << fb->_sendcells.size() << "," 
	    << fb->_recvcells.size()
	    << ")" << "\n"
	    << "-----------------------------------------------------" << "\n";
    }
    *_out   
      << "========================================================" << "\n"
      << "Domain Boundaries: " << "\n"
      << "========================================================" << "\n";
    unsigned int patch = 0;
    while(patch < _patches.size()){
      *_out << " Local patch #" << patch+1 << "\n"
	    << "-----------------------------------------------------" << "\n";
      FluPatch *fp = &_patches[patch++];
      *_out 
	<< "  Patch ID: " << fp->_id << "\n"
	<< "   Triangles: (" << fp->_ntri << "," << fp->_ngtri << ")" << "\n"
	<< "   Quads:     (" << fp->_nquad << "," << fp->_ngquad << ")" 
	<< "\n"
	<< "-----------------------------------------------------" << "\n";
    }
    *_out   
      << "========================================================" << "\n"
      << "\n"
      << "********************************************************" << "\n";
  }
}

bool
FluRegion::ReadControlFile()
{
  ifstream Inf;
  Inf.open("Rocflu/RocfluControl.txt");
  if(!Inf)
    Inf.open("RocfluControl.txt");
  if(!Inf)
    return(false);
  Inf >> _casename >> _casepath;
  Inf.close();
  return(true);
}
bool 
FluRegion::WriteFluNative(const string &prefix)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casepath.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  WriteFluCMP(pre);
  WriteFluDIM(pre);
  WriteFluCOM(pre);
  return(true);
} 

bool
FluRegion::WriteFluCOM(const string &prefix)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename + ".";
  }
  else
    pre += ".";
  //  if(pre.empty())
  //    pre = "./";
  ofstream Ouf;
  ostringstream Ostr;
  Ostr << pre << "com_";
  Ostr << setw(5) << setfill('0');
  Ostr << _id;
  Ouf.open(Ostr.str().c_str());
  if(!Ouf)
    return(false);
  Ouf << "# ROCFLU communication lists file" << endl
      << "# Dimensions"  << endl
      << setw(8) << _borders.size() << endl
      << "# Information" << endl;
  vector<FluBorder>::iterator fbi = _borders.begin();
  while(fbi != _borders.end()){
    Ouf << setw(8) << fbi->_rpart << setw(8) << fbi->_rbid << endl;
    fbi++;
  }
  Ouf << "# Cells" << endl;
  fbi = _borders.begin();
  while(fbi != _borders.end()){
    Ouf << setw(8) << fbi->_sendcells.size() << setw(8) 
	<< fbi->_recvcells.size() << endl;
    int line = 0;
    vector<unsigned int>::iterator sci = fbi->_sendcells.begin();
    while(sci != fbi->_sendcells.end()){
      Ouf << setw(8) << *sci++;
      line++;
      if(line == 10){ 
	Ouf << endl;
	line = 0;
      }
    }
    if(line){
      Ouf << endl;
      line = 0;
    }
    sci = fbi->_recvcells.begin();
    while(sci != fbi->_recvcells.end()){
      Ouf << setw(8) << *sci++;
      line++;
      if(line == 10){ 
	Ouf << endl;
	line = 0;
      }
    }
    if(line){
      Ouf << endl;
      line = 0;
    }
    fbi++;
  }    
  Ouf << "# Vertices" << endl;
  fbi = _borders.begin();
  while(fbi != _borders.end()){
    Ouf << setw(8) << fbi->_sendnodes.size() << setw(8) 
	<< fbi->_recvnodes.size() << setw(8) << fbi->_sharenodes.size() 
	<< endl;
    vector<unsigned int>::iterator ni = fbi->_sendnodes.begin();
    int line = 0;
    while(ni != fbi->_sendnodes.end()){
      Ouf << setw(8) << *ni++;
      line++;
      if(line == 10){
	Ouf << endl;
	line = 0;
      }
    } 
    if(line || fbi->_sendnodes.empty()){
      Ouf << endl;
      line = 0;
    }
    ni = fbi->_recvnodes.begin();
    line = 0;
    while(ni != fbi->_recvnodes.end()){
      Ouf << setw(8) << *ni++;
      line++;
      if(line == 10){
	Ouf << endl;
	line = 0;
      }
    } 
    if(line || fbi->_recvnodes.empty()){
      Ouf << endl;
      line = 0;
    }
    ni = fbi->_sharenodes.begin();
    line = 0;
    while(ni != fbi->_sharenodes.end()){
      Ouf << setw(8) << *ni++;
      line++;
      if(line == 10){
	Ouf << endl;
	line = 0;
      }
    } 
    if(line || fbi->_sharenodes.empty()){
      Ouf << endl;
      line = 0;
    }
    fbi++;
  }
  Ouf << "# End" << endl;
  return(true);
}

bool
FluRegion::ReadFluCOM(const string &prefix)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  //  if(pre.empty())
  //    pre = "./";
  ifstream Inf;
  ostringstream Ostr;
  Ostr << pre << ".com_";
  Ostr << setw(5) << setfill('0');
  Ostr << _id;
  Inf.open(Ostr.str().c_str());
  if(!Inf)
    return(false);
  string line;
  SkipLines(Inf,2);
  unsigned int nborders = 0;
  Inf >> nborders;
  if(_debug && _out)
    *_out << "FluRegion::ReadFluCOM(" << _id << "): nborders: " 
	 << nborders << "\n";
  assert(nborders == _borders.size());
  SkipLines(Inf,2);
  vector<FluBorder>::iterator fbi = _borders.begin();
  while(fbi != _borders.end()){
    unsigned int rpart,rbid;
    Inf >> rpart >> rbid;
    assert(fbi->_rpart == rpart && fbi->_rbid == rbid);
    fbi++;
  }
  SkipLines(Inf,2);
  fbi = _borders.begin();
  while(fbi != _borders.end()){
    unsigned int rcells,scells;
    Inf >> scells >> rcells;
    assert(scells == fbi->_sendcells.size() &&
	   rcells == fbi->_recvcells.size());
    int cell = 0;
    int ncells = fbi->_sendcells.size();
    while(cell < ncells)
      Inf >> fbi->_sendcells[cell++];
    cell = 0;
    ncells = fbi->_recvcells.size();
    while(cell < ncells)
      Inf >> fbi->_recvcells[cell++];
    fbi++;
  }    
  SkipLines(Inf,2);
  fbi = _borders.begin();
  while(fbi != _borders.end()){
    unsigned int nshare,nsend,nrecv;
    Inf >> nsend >> nrecv >> nshare;
    assert(nsend  == fbi->_sendnodes.size() &&
	   nrecv  == fbi->_recvnodes.size() &&
	   nshare == fbi->_sharenodes.size());
    unsigned int node = 0;
    unsigned int nnodes = fbi->_sendnodes.size();
    while(node < nnodes)
      Inf >> fbi->_sendnodes[node++];
    node = 0;
    nnodes = fbi->_recvnodes.size();
    while(node < nnodes)
      Inf >> fbi->_recvnodes[node++];
    nnodes = fbi->_sharenodes.size();
    node = 0;
    while(node < nnodes)
      Inf >> fbi->_sharenodes[node++];
    fbi++;
  }
  Inf.close();
  return(true);
}

bool
FluRegion::CreateRegionMapFile(const string &prefix,unsigned int nproc,
			       unsigned int nregions)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename + ".";
  }
  else
    pre += ".";
  //  if(pre.empty())
  //    pre = "./";
  ofstream Ouf;
  ostringstream Ostr;
  Ostr << pre << "map";
  Ouf.open(Ostr.str().c_str());
  if(!Ouf)
    return(false);
  Ouf << "# ROCFLU region mapping file" << endl
      << "# Number of regions" << endl
      << setw(8) << nregions << endl
      << "# Number of processes" << endl
      << setw(8) << nproc << endl;
  unsigned int regpproc = nregions/nproc;
  if(regpproc == 0)
    return(false);
  unsigned int left     = nregions%nproc;
  unsigned int proc = 0;
  unsigned int reg = 0;
  while(proc < nproc){
    unsigned int nreg = regpproc;
    if(left > 0){
      nreg++;
      left--;
    }
    Ouf << "# Process " << setw(6) << setfill('0') << proc+1 << endl
	<< setfill(' ') << setw(8) << nreg << endl;
    unsigned int uplimit = reg + nreg;
    while(reg < uplimit)
      Ouf << setw(8) << (1+reg++) << endl;
    proc++;
  }
  Ouf << "# End" << endl;
  Ouf.close();
  return(true);
}

bool
FluRegion::WriteFluDIM(const string &prefix,double t,bool unsteady)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename + ".";
  }
  else
    pre += ".";
  //  if(pre.empty())
  //    pre = "./";
  ofstream Ouf;
  ostringstream Ostr;
  Ostr << pre << "dim_";
  Ostr << setw(5) << setfill('0');
  Ostr << _id;
  string filebase(Ostr.str());
  if(unsteady){
    Ostr.clear();
    Ostr.str("");
    Ostr << "_" << scientific << setprecision(5) << t;
    string timestring(Ostr.str());
    timestring.replace(8,1,"E");
    filebase+=timestring;
  }
  Ouf.open(filebase.c_str());
  if(!Ouf)
    return(false);
  _ntet = _tetconn.size()/4;
  _nnodes = _nc.size()/3;
  _nhex = _hexconn.size()/8;
  _npris = _prisconn.size()/6;
  _npyr = _pyrconn.size()/5;
  Ouf << "# ROCFLU dimensions file" << endl
      << "# Vertices" << endl
      << setw(8) << _nnodes - _ngnodes << setw(8) << _nnodes
      << setw(8) << _nnodes+(_nnodes/5) << endl
      << "# Cells" << endl
      << setw(8) << (nelem()-_ngtet-_ngpyr-_ngpris-_nghex)
      << setw(8) << nelem() 
      << setw(8) << (nelem()+(nelem()/5)) << endl
      << "# Tetrahedra" << endl
      << setw(8) << _ntet - _ngtet << setw(8) << _ntet << setw(8) 
      << (_ntet+(_ntet/5))
      << endl
      << "# Hexahedra" << endl
      << setw(8) << _nhex - _nghex << setw(8) << _nhex << setw(8) 
      << (_nhex+(_nhex/5))
      << endl
      << "# Prisms" << endl
      << setw(8) << _npris - _ngpris << setw(8) << _npris << setw(8) 
      << _npris+(_npris/5)
      << endl
      << "# Pyramids" << endl
      << setw(8) << _npyr - _ngpyr << setw(8) << _npyr 
      << setw(8) << _npyr+(_npyr/5)
      << endl
      << "# Patches" << endl
      << setw(8) << _patches.size() << setw(8) << _npatches_total << endl;
  vector<FluPatch>::iterator fpi = _patches.begin();
  while(fpi != _patches.end()){
    Ouf << setw(8) << fpi->_id << setw(8) << fpi->_ntri - fpi->_ngtri
	<< setw(8) << fpi->_ntri << setw(8) << fpi->_nquad - fpi->_ngquad
	<< setw(8) << fpi->_nquad << endl;
    fpi++;
  }
  Ouf << "# Borders" << endl
      << setw(8) << _borders.size() << endl;
  vector<FluBorder>::iterator fbi = _borders.begin();
  while(fbi != _borders.end()){
    Ouf << setw(8) << fbi->_rpart << setw(8) << fbi->_rbid << setw(8) 
	<< fbi->_sendcells.size() << setw(8) << fbi->_recvcells.size() 
	<< setw(8) 
	<< fbi->_sendnodes.size() << setw(8) << fbi->_recvnodes.size() 
	<< setw(8) 
	<< fbi->_sharenodes.size() << endl;
    fbi++;
  }
  Ouf << "# End" << endl;
  return(true);
}

bool
FluRegion::ReadFluDIM(const string &prefix,double t,bool unsteady)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  //  if(pre.empty())
  //    pre = "./";
  if(_debug && _out)
    *_out << "FluRegion::ReadFluDIM: Entry" << "\n";
  _soln._current_time = t;
  _soln._unsteady = unsteady;
  ifstream Inf;
  ostringstream Ostr;
  Ostr << pre << ".dim_";
  Ostr << setw(5) << setfill('0');
  Ostr << _id;
  string filebase(Ostr.str());
  if(unsteady){
    Ostr.clear();
    Ostr.str("");
    Ostr << "_" << scientific << setprecision(5) << t;
    string timestring(Ostr.str());
    timestring.replace(8,1,"E");
    filebase+=timestring;
  }
  Inf.open(filebase.c_str());
  if(!Inf){
    if(_debug && _out)
      *_out << "FluRegion::ReadFluDIM: Unable to open " << filebase 
	   << " for reading." << endl;
    return(false);
  }
  string line;
  unsigned int nreal_cells = 0;
  unsigned int ncells = 0;
  unsigned int npatches = 0;
  unsigned int nreal_tets = 0;
  unsigned int nreal_hex = 0;
  unsigned int nreal_pris = 0;
  unsigned int nreal_pyr = 0;
  unsigned int nreal_nodes = 0;
  SkipLines(Inf,2);
  Inf >> nreal_nodes >> _nnodes;
  SkipLines(Inf,2);
  Inf >> nreal_cells >> ncells;
  SkipLines(Inf,2);
  Inf >> nreal_tets >> _ntet;
  SkipLines(Inf,2);
  Inf >> nreal_hex >> _nhex;
  SkipLines(Inf,2);
  Inf >> nreal_pris >> _npris;
  SkipLines(Inf,2);
  Inf >> nreal_pyr >> _npyr;
  SkipLines(Inf,2);
  Inf >> npatches >> _npatches_total;
  _ngnodes = _nnodes - nreal_nodes;
  _ngtet = _ntet - nreal_tets;
  _nghex = _nhex - nreal_hex;
  _ngpyr = _npyr - nreal_pyr;
  _ngpris = _npris - nreal_pris;
  assert(nelem() == ncells);
  if(_debug && _out)
    *_out << "FluRegion::ReadFluDIM:: Elems: (" << nelem() << "," 
	 << _ngtet+_nghex+_ngpyr+_ngpris << ")" << endl
	 << "FluRegion::ReadFluDIM:: Nodes: (" << _nnodes << "," 
	 << _ngnodes << ")" << endl
	 << "FluRegion::ReadFluDIM:: Patches: (" << npatches << "," 
	 << _npatches_total << ")" << endl;  
  _patches.resize(npatches);
  _tetconn.resize(4*_ntet);
  _hexconn.resize(8*_nhex);
  _prisconn.resize(6*_npris);
  _pyrconn.resize(5*_npyr);
  _nc.resize(3*_nnodes);
  vector<FluPatch>::iterator fpi = _patches.begin();
  while(fpi != _patches.end()){
    unsigned int nreal_tri = 0;
    unsigned int nreal_quad = 0;
    Inf >> fpi->_id >> nreal_tri >> fpi->_ntri >> nreal_quad >> fpi->_nquad;
    fpi->_ngtri  = fpi->_ntri - nreal_tri;
    fpi->_ngquad = fpi->_nquad - nreal_quad;
    fpi->_triconn.resize(fpi->_ntri*3);
    fpi->_quadconn.resize(fpi->_nquad*4);
    fpi++;
  }
  unsigned int nborders = 0;
  SkipLines(Inf,2);
  Inf >> nborders;
  if(_debug && _out)
    *_out << "FluRegion::ReadFluDIM:: Number of borders: " << nborders << endl;
  _borders.resize(nborders);
  vector<FluBorder>::iterator fbi = _borders.begin();
  while(fbi != _borders.end()){
    unsigned int nsend,nrecv,nshared,csend,crecv;
    Inf >> fbi->_rpart >> fbi->_rbid >> csend >> crecv
	>> nsend >> nrecv >> nshared;
    if(_debug && _out)
      *_out << "FluRegion::ReadFluDIM:: Border: (" << fbi->_rpart << "," 
	   << fbi->_rbid << "," << csend << "," << crecv << "," << nsend 
	   << "," << nrecv << "," << nshared << ")" << endl;
    fbi->_sharenodes.resize(nshared);
    fbi->_sendnodes.resize(nsend);
    fbi->_recvnodes.resize(nrecv);
    fbi->_sendcells.resize(csend);
    fbi->_recvcells.resize(crecv);
    fbi++;
  }
  Inf.close();
  if(_debug && _out)
    *_out << "FluRegion::ReadFluDIM: Exit" << endl;
  return(true);
}

bool
FluRegion::PopRemBordIndFILE(const string &prefix,double t,bool unsteady)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  vector<FluBorder>::iterator fbi = _borders.begin();
  while(fbi != _borders.end()){
    ifstream Inf;
    ostringstream Ostr;
    Ostr << pre << ".dim_";
    Ostr << setw(5) << setfill('0');
    Ostr << fbi->_rpart;
    string filebase(Ostr.str());
    if(unsteady){
      Ostr.clear();
      Ostr.str("");
      Ostr << "_" << scientific << setprecision(5) << t;
      string timestring(Ostr.str());
      timestring.replace(8,1,"E");
      filebase+=timestring;
    }
    Inf.open(filebase.c_str());
    if(!Inf){
      if(_debug && _out)
	*_out 
	     << "FluRegion(" << _id << ")::PopRemBordIndFILE: Unable to open "
	     << filebase << " for reading.\n";
      return(false);
    }
    string line;
    unsigned int fakeout = 0;
    unsigned int nreal_cells = 0;
    unsigned int nnodes = 0;
    unsigned int nhex = 0;
    unsigned int ntet = 0;
    unsigned int npris = 0;
    unsigned int npyr = 0;
    unsigned int ncells = 0;
    unsigned int npatches = 0;
    unsigned int nreal_tets = 0;
    unsigned int nreal_hex = 0;
    unsigned int nreal_pris = 0;
    unsigned int nreal_pyr = 0;
    unsigned int nreal_nodes = 0;
    unsigned int npatches_total = 0;
    unsigned int ngnodes = 0;
    unsigned int ngtet = 0;
    unsigned int nghex = 0;
    unsigned int ngpris = 0;
    unsigned int ngpyr = 0;
    SkipLines(Inf,2);
    Inf >> nreal_nodes >> nnodes;
    SkipLines(Inf,2);
    Inf >> nreal_cells >> ncells;
    SkipLines(Inf,2);
    Inf >> nreal_tets >> ntet;
    SkipLines(Inf,2);
    Inf >> nreal_hex >> nhex;
    SkipLines(Inf,2);
    Inf >> nreal_pris >> npris;
    SkipLines(Inf,2);
    Inf >> nreal_pyr >> npyr;
    SkipLines(Inf,2);
    Inf >> npatches >> npatches_total;
    ngnodes = nnodes - nreal_nodes;
    ngtet = ntet - nreal_tets;
    nghex = nhex - nreal_hex;
    ngpyr = npyr - nreal_pyr;
    ngpris = npris - nreal_pris;
    if(_debug && _out)
      *_out  << "FluRegion(" << _id << ")::PopRemBordIndFILE::Elems: (" 
	    << ntet+npris+nhex+npyr
	    << "," << ngtet+nghex+ngpyr+ngpris << ")\n"
	    << "FluRegion(" << _id << ")::PopRemBordIndFILE" << "::Nodes: (" 
	    << nnodes << "," 
	    << ngnodes << ")\n"
	    << "FluRegion(" << _id << ")::PopRemBordIndFILE:: Patches: (" 
	    << npatches << "," 
	    << npatches_total << ")\n";
    for(unsigned int pind = 0;pind < npatches;pind++){
      unsigned int nreal_tri = 0;
      unsigned int nreal_quad = 0;
      Inf >> fakeout >> nreal_tri >> fakeout >> nreal_quad >> fakeout;
    }
    unsigned int nborders = 0;
    SkipLines(Inf,2);
    Inf >> nborders;
    if(_debug && _out)
      *_out  << "FluRegion(" << _id << ")::PopRemBordIndFILE:: " 
	   << "Number of borders: " << nborders 
	   << "\n";
    bool done = false;
    for(unsigned int bind = 0;(bind < nborders && !done);bind++){
      unsigned int rpart;
      Inf >> rpart;
      if(rpart == this->_id){
	fbi->_rbid = bind+1;
	done = true;
      }
      else
	Inf >> fakeout >> fakeout >> fakeout >> fakeout >> fakeout >> fakeout;
    }
    if(!done)
      return(false);
    Inf.close();
    //    this->WriteFluDIM(prefix,t,unsteady);
    fbi++;
  }
  return(true);
}

bool 
FluRegion::WriteFluCMP(const string &prefix)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename + ".";
  }
  else
    pre += ".";
  ofstream Ouf;
  ostringstream Ostr;
  Ostr << pre << "cmp_";
  Ostr << setw(5) << setfill('0');
  Ostr << _id;
  Ouf.open(Ostr.str().c_str());
  if(!Ouf)
    return(false);
  _ntet = _tetconn.size()/4;
  _nhex = _hexconn.size()/8;
  _npris = _prisconn.size()/6;
  _npyr = _pyrconn.size()/5;
  //  BuildCellMapping();
  Ouf << "# ROCFLU cell mapping file" << endl
      << "# Dimensions" << endl
      << setw(8) << _ntet << setw(8) << _nhex << setw(8) << _npris 
      << setw(8) << _npyr << endl;
  vector<pair<unsigned int,unsigned int> >::iterator ci = _cellmap.begin();
  //  int line = 0;
  unsigned int el = 0;
  unsigned int npl = 10;
  if(_ntet > 0){
    el = 0;
    Ouf << "# Tetrahedra" << endl;
    while(el < _ntet){
      Ouf << setw(8) << Elem2Cell(make_pair(1,++el));
      if(!(el%npl))
	Ouf << endl;
    }
    if(el%npl) Ouf << endl;
  }
  if(_nhex > 0){
    el = 0;
    Ouf << "# Hexahedra" << endl;
    while(el < _nhex){
      Ouf << setw(8) << Elem2Cell(make_pair(4,++el));
      if(!(el%npl))
	Ouf << endl;
    }
    if(el%npl) Ouf << endl;
  }
  if(_npris > 0){
    el = 0;
    Ouf << "# Prisms" << endl;
    while(el < _npris){
      Ouf << setw(8) << Elem2Cell(make_pair(3,++el));
      if(!(el%npl))
	Ouf << endl;
    }
    if(el%npl) Ouf << endl;
  }
  if(_npyr > 0){
    el = 0;
    Ouf << "# Pyramids" << endl;
    while(el<_npyr){
      Ouf << setw(8) << Elem2Cell(make_pair(2,++el));
      if(!(el%npl))
	Ouf << endl;
    }
    if(el%npl) Ouf << endl;
  }
  Ouf << "# End" << endl;
  Ouf.close();
  return(true);
}


bool
FluRegion::InitVolumeSoln(const string &prefix)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  string fname(pre+".inp");
  ifstream Inf;
  Inf.open(fname.c_str());
  if(!Inf)
    return false;
  string line;
  double cp    = -1.0;
  double gamma = -1.0;
  double velx  =  0.0;
  double vely  =  0.0;
  double velz  =  0.0;
  double press = -1.0;
  double rho   = -1.0;
  while(getline(Inf,line)){
    if(line[0] != '#'){
      string::size_type xtime  = line.find("STARTTIME");
      string::size_type ftype  = line.find("FLOWTYPE");
      string::size_type xgamma = line.find("GAMMA");
      string::size_type xpress = line.find("PRESS");
      string::size_type x_dens = line.find("DENS");
      string::size_type x_velx = line.find("VELX");
      string::size_type x_vely = line.find("VELY");
      string::size_type x_velz = line.find("VELZ");
      string::size_type x_cp   = line.find("CP");
      string tstr;
      istringstream Istr(line);
      if(xtime  != string::npos)
	Istr >> tstr >> _soln._current_time;
      if(ftype  != string::npos)
	Istr >> tstr >> _soln._unsteady;
      if(xgamma != string::npos)
	Istr >> tstr >> gamma;
      if(xpress != string::npos)
	Istr >> tstr >> press;
      if(x_dens != string::npos)
	Istr >> tstr >> rho;
      if(x_velx != string::npos)
	Istr >> tstr >> velx;
      if(x_vely != string::npos)
	Istr >> tstr >> vely;
      if(x_velz != string::npos)
	Istr >> tstr >> velz;
      if(x_cp   != string::npos)
	Istr >> tstr >> cp;
    }
  }
  if(gamma < 0 || press < 0 || rho < 0 || cp < 0)
    return(false);
  _ntet = _tetconn.size()/4;
  _npyr = _pyrconn.size()/5;
  _nhex = _hexconn.size()/8;
  _npris = _prisconn.size()/6;
  unsigned int nelem = _ntet + _npyr + _npris + _nhex;
  _soln._rhof.resize(nelem,rho);
  double E = press/(rho*(gamma - 1.0)) + 
    (velx*velx + vely*vely + velz*velz)/2.0;
  _soln._rhovf.resize(3*nelem,0.0);
  if(velx != 0.0 || vely != 0.0 || velz != 0.0){
    unsigned int index = 0;
    for(unsigned int el = 0;el < nelem;el++){
      _soln._rhovf[index++] = rho*velx;
      _soln._rhovf[index++] = rho*vely;
      _soln._rhovf[index++] = rho*velz;
    }
  }
  double R = (cp*(gamma-1.0))/gamma;
  double Ti = press/(rho*R);
  double C  = sqrt(gamma*press/rho);
  _soln._rhoEf.resize(nelem,rho*E);
  _soln._pf.resize(nelem,press); 
  _soln._Tf.resize(nelem,Ti);
  _soln._af.resize(nelem,C);
  _soln._gsp.resize(_nvface,0.0);
  return(true);
}

bool
FluRegion::InitSurfaceSoln(const string &prefix)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  vector<FluPatch>::iterator fpi = _patches.begin();
  while(fpi != _patches.end()){
    FluPatch &fp = *fpi++;
    if(!fp.InitSurfaceSoln(pre,true))
      return(false);
  }
  return(true);
}


bool 
FluPatch::InitSurfaceSoln(const string &prefix,bool all)
{
  string pre(prefix);
  if(_debug && _out)
    *_out << "FluPatch(" << _id << ")::InitSurfaceSoln Enter" << endl;
  if(pre.empty())
    pre = "./";
  if(all){
    _ntri = _triconn.size()/3;
    _nquad = _quadconn.size()/4;
    //#ifdef _ROCSTAR_X_
    _soln.Resize(_ntri+_nquad,surface_coordinates.size()/3);
    //#else
    //    _soln.gsp.resize(_ntri+_nquad - (_ngtri+_ngquad),0.0);
    //#endif
  }
  //#ifdef _ROCSTAR_X_
  int constraint_map[] = { 2 , 120, 121, 122, -122, -121, -120, 0};
  // Populate the Rocstar data living on the surface
  string fname(pre+".bc");
  ifstream Inf;
  Inf.open(fname.c_str());
  if(!Inf){
    if(_out){
      *_out << "FluPatch::InitSurfaceSoln: Error: could not locate bc file."
	   << endl;
    }
    return false;
  }
  string line;
  bool found_patch = false;
  while(getline(Inf,line) && !found_patch){
    if(line[0] != '#'){
      string::size_type patchx = line.find("PATCH");
      if( patchx != string::npos){
	string tstr;
	unsigned int ubound,lbound;
	istringstream Istr(line);
	Istr >> tstr >> lbound >> ubound;
	if(_id >= lbound && _id <= ubound)
	  found_patch = true;
      }
    }
  }
  if(!found_patch){
    if(_out){
      *_out << "FluPatch::InitSurfaceSoln: Error: Did not find patch " 
	   << _id << " in bc file." << endl;
    }
    return(false);
  }
  // Now we want to search forward for the relevant data until we encounter #
  bool end_section = false;
  _constr_type = 2;
  _bcflag      = 2;
  while(getline(Inf,line) && !end_section){
    if(line[0] == '#') end_section = true;
    else{
      string::size_type movx  = line.find("MOVEDIR");
      string::size_type coupx = line.find("COUPLED");
      string::size_type burn  = line.find("BFLAG");
      if(movx != string::npos){
	string tstr;
	int mvdir;
	istringstream Istr(line);
	Istr >> tstr >> mvdir;
	_constr_type = constraint_map[mvdir];
      }
      if(coupx != string::npos){
	string tstr;
	istringstream Istr(line);
	Istr >> tstr >> _bcflag;
      }
      if(burn != string::npos){
	string tstr;
	unsigned int burnf = 0;
	istringstream Istr(line);
	Istr >> tstr >> burnf;
	_soln._bflag.resize(_ntri+_nquad,burnf);
      }
    }
  }
  //#endif
  Inf.close();
  if(_debug && _out)
    *_out << "FluPatch(" << _id << ")::InitSurfaceSoln Exit" << endl;
  return true;
}

// Populates our connectivites from the ASCII grid files
// pre should contain the path + casename
// Ex: /my/directory/StarSlice
bool
FluRegion::ReadFluGridASCII(const string &prefix,double t,bool unsteady)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  ostringstream Ostr;
  Ostr << pre << ".grda_" << setw(5) << setfill('0') << _id;
  string filebase(Ostr.str());
  Ostr.clear();
  Ostr.str("");
  Ostr << "_" << scientific << setprecision(5) << t;
  string timestring(Ostr.str());
  timestring.replace(8,1,"E"); 
  if(unsteady) filebase+=timestring;
  ifstream Inf;
  Inf.open(filebase.c_str());
  if(!Inf){
    if(_debug && _out)
      *_out << "FluRegion::ReadFluGridASCII: Unable to open " << filebase 
	   << " for reading.\n";
    return(false);
  }
  string line;
  SkipLines(Inf,6);
  unsigned int nnode,ntet,nhex,npris,npyr;
  Inf >> nnode >> ntet >> nhex >> npris >> npyr;
  assert(nnode == _nnodes && ntet == _ntet && 
	 nhex == _nhex    && npris == _npris &&
	 npyr == _npyr    && _nc.size() == _nnodes*3);
  SkipLines(Inf,2);
  for(int i = 0; i < 3;i++){
    unsigned int node = 0;
    while(node < _nnodes)
      Inf >> _nc[(node++)*3+i];
  }
  if(_ntet > 0){
    _tetconn.resize(_ntet*4);
    SkipLines(Inf,2);
    for(int i = 0;i < 4;i++){
      unsigned int el = 0;
      while(el < _ntet)
	Inf >> _tetconn[(el++)*4+i];
    }
  }
  if(_nhex > 0){
    _hexconn.resize(_nhex*8);
    SkipLines(Inf,2);
    for(int i = 0;i < 8;i++){
      unsigned int el = 0;
      while(el < _nhex)
	Inf >> _hexconn[(el++)*8+i];
    }
  }
  if(_npris > 0){
    _prisconn.resize(_npris*6);
    SkipLines(Inf,2);
    for(int i = 0;i < 6;i++){
      unsigned int el = 0;
      while(el < _npris)
	Inf >> _prisconn[(el++)*6+i];
    }
  }
  if(_npyr > 0){
    _pyrconn.resize(_npyr*5);
    SkipLines(Inf,2);
    for(int i = 0;i < 5;i++){
      unsigned int el = 0;
      while(el < _npyr)
	Inf >> _pyrconn[(el++)*5+i];
    }
  }
  SkipLines(Inf,2);
  unsigned int npatches = 0;
  unsigned int patch    = 0;
  Inf >> npatches;
  assert(npatches == _patches.size());
  while(patch < npatches){
    FluPatch &fp = _patches[patch++];
    unsigned int ntri,nquad;
    Inf >> ntri >> nquad;
    assert(ntri == fp._ntri && nquad == fp._nquad);
    if(fp._ntri > 0){
      fp._triconn.resize(3*fp._ntri);
      for(int i = 0;i < 3;i++){
	unsigned int tri = 0;
	while(tri < fp._ntri)
	  Inf >> fp._triconn[(tri++)*3+i];
      }
    }
    if(fp._nquad > 0){
      fp._quadconn.resize(4*fp._nquad);
      for(int i = 0;i < 4;i++){
	unsigned int quad = 0;
	while(quad < fp._nquad)
	  Inf >> fp._quadconn[(quad++)*4+i];
      }
    }
  }
  Inf.close();
  return true;
}

// Populates our connectivites from the ASCII grid files
// pre should contain the path + casename
// Ex: /my/directory/StarSlice
bool
FluRegion::WriteFluGridASCII(const string &prefix,double t,bool unsteady)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename + ".";
  }
  else
    pre += ".";
  ostringstream Ostr;
  Ostr << pre << "grda_" << setw(5) << setfill('0') << _id;
  string filebase(Ostr.str());
  if(unsteady){
    Ostr.clear();
    Ostr.str("");
    Ostr << "_" << scientific << setprecision(5) << t;
    string timestring(Ostr.str());
    timestring.replace(8,1,"E");
    filebase+=timestring;
  }
  ofstream Ouf;
  Ouf.open(filebase.c_str());
  if(!Ouf)
    return(false);
  _nnodes = _nc.size()/3;
  _ntet = _tetconn.size()/4;
  _nhex = _hexconn.size()/8;
  _npyr = _pyrconn.size()/5;
  _npris = _prisconn.size()/6;
  Ouf << "# ROCFLU grid file" << endl
      << "# Precision and range" << endl
      << setw(8) << 15 << setw(8) << 307 << endl
      << "# Physical time" << endl
      << scientific << setprecision(16) << setw(23)
      << _soln._current_time << endl
      << "# Dimensions" << endl
      << setw(8) << _nnodes << setw(8) << _ntet << setw(8)
      << _nhex << setw(8) << _npris << setw(8) << _npyr << endl
      << "# Coordinates" << endl;
  unsigned int ll = 5;
  for(int i = 0; i < 3;i++){
    unsigned int node = 0;
    while(node < _nnodes){
      Ouf << setw(23) << scientific << setprecision(16) 
	  << _nc[(node++)*3+i];
      if(!(node%ll))
	Ouf << endl;
    }
    if(node%ll)
      Ouf << endl;
  }
  ll = 10;
  if(_ntet > 0){
    Ouf << "# Tetrahedra" << endl;
    for(int i = 0;i < 4;i++){
      unsigned int el = 0;
      while(el < _ntet){
	Ouf << setw(8) << _tetconn[(el++)*4+i];
	if(!(el%ll))
	  Ouf << endl;
      }
      if(el%ll)
	Ouf << endl;
    }
  }
  if(_nhex > 0){
    Ouf << "# Hexahedra" << endl;
    for(int i = 0;i < 8;i++){
      unsigned int el = 0;
      while(el < _nhex){
	Ouf << setw(8) << _hexconn[(el++)*8+i];
	if(!(el%ll))
	  Ouf << endl;
      }
      if(el%ll)
	Ouf << endl;
    }
  }
  if(_npris > 0){
    Ouf << "# Prisms" << endl;
    for(int i = 0;i < 6;i++){
      unsigned int el = 0;
      while(el < _nhex){
	Ouf << setw(8) << _prisconn[(el++)*6+i];
	if(!(el%ll))
	  Ouf << endl;
      }
      if(el%ll)
	Ouf << endl;
    }
  }
  if(_npyr > 0){
    Ouf << "# Pyramids" << endl;
    for(int i = 0;i < 5;i++){
      unsigned int el = 0;
      while(el < _npyr){
	Ouf << setw(8) << _pyrconn[(el++)*5+i];
	if(!(el%ll))
	  Ouf << endl;
      }
      if(el%ll)
	Ouf << endl;
    }
  }
  Ouf << "# Boundaries" << endl
      << setw(8) << _patches.size() << endl;
  unsigned int npatches = _patches.size();
  unsigned int patch = 0;
  while(patch < npatches){
    FluPatch &fp = _patches[patch++];
    Ouf << setw(8) <<  fp._ntri << setw(8) <<  fp._nquad << endl;
    if(fp._ntri > 0){
      for(int i = 0;i < 3;i++){
	unsigned int tri = 0;
	while(tri < fp._ntri){
	  Ouf << setw(8) <<  fp._triconn[(tri++)*3+i];
	  if(!(tri%ll))
	    Ouf << endl;
	}
	if(tri%ll)
	  Ouf << endl;
      }
    }
    if(fp._nquad > 0){
      for(int i = 0;i < 4;i++){
	unsigned int quad = 0;
	while(quad < fp._nquad){
	  Ouf << setw(8) <<  fp._quadconn[(quad++)*4+i];
	  if(!(quad%ll))
	    Ouf << endl;
	}
	if(quad%ll)
	  Ouf << endl;
      }
    }
  }
  Ouf << "# End" << endl;
  Ouf.close();
  return true;
}

bool
FluRegion::ReadFluSolnASCII(const string &prefix,unsigned int n,
			    double t,bool unsteady)
{
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  ostringstream Ostr;
  Ostr << pre << ".mixt.cva_" << setw(5) << setfill('0') << _id << "_";
  string filebase(Ostr.str());
  Ostr.str("");
  Ostr << setw(6) << setfill('0') << n;
  string iterstring(Ostr.str());
  Ostr.clear();
  Ostr.str("");
  Ostr << scientific << setprecision(5) << t;
  string timestring(Ostr.str());
  timestring.replace(7,1,"E");
  ifstream Inf;
  if(unsteady) filebase+=timestring;
  else filebase+=iterstring;
  Inf.open(filebase.c_str());
  if(!Inf){
    if(_debug && _out)
      *_out << "FluRegion::ReadFluSolnASCII: Unable to open " << filebase 
	   << " for reading.\n";
    return(false);
  }
  string line;
  SkipLines(Inf,4);
  Inf >> _soln._resid;
  SkipLines(Inf,2);
  Inf >> _soln._current_time;
  SkipLines(Inf,2);
  unsigned int ncells = 0;
  Inf >> ncells;
  assert(ncells == nelem());
  SkipLines(Inf,2);
  unsigned int cell = 0;
  _soln._rhof.resize(ncells,0.0);
  _soln._rhovf.resize(3*ncells,0.0);
  _soln._rhoEf.resize(ncells,0.0);
  while(cell < ncells)
    Inf >> _soln._rhof[cell++];
  SkipLines(Inf,2);  
  for(int i = 0;i < 3;i++){
    cell = 0;
    while(cell < ncells)
      Inf >> _soln._rhovf[(cell++)*3+i];
    SkipLines(Inf,2);
  }
  cell = 0;
  while(cell < ncells)
    Inf >> _soln._rhoEf[cell++];
  if(!ReadGSPASCII(pre,n,t,unsteady)){
    if(_debug && _out)
      *_out 
	<< "FluRegion::ReadFluSolnASCII: Unable to process grid speed file."
	<< "\n";
    return(false);
  }
  Inf.close();
  return(true);
}

bool
FluRegion::WriteFluSolnASCII(const string &prefix,unsigned int n,
			     double t,bool unsteady)
{
  if(_debug && _out)
    *_out << "FluRegion::WriteFluSolnASCII: Enter\n";
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename + ".";
  }
  else
    pre += ".";
  ostringstream Ostr;
  Ostr << pre << "mixt.cva_" << setw(5) << setfill('0') << _id << "_";
  string filebase(Ostr.str());
  Ostr.str("");
  Ostr << setw(6) << setfill('0') << n;
  string iterstring(Ostr.str());
  Ostr.clear();
  Ostr.str("");
  Ostr << scientific << setprecision(5) << t;
  string timestring(Ostr.str());
  timestring.replace(7,1,"E");
  ofstream Ouf;
  if(!unsteady)
    filebase+=iterstring;
  else
    filebase+=timestring;
  Ouf.open(filebase.c_str());
  if(!Ouf)
    return false;
  unsigned int ncells = nelem();
  unsigned int n_cvar = 5;
  Ouf << "# ROCFLU flow file" << endl
      << "# Precision and range" << endl
      << setw(8) << 15 << setw(8) << 307 << endl
      << "# Initial residual" << endl
      << setw(23) << scientific << setprecision(16)
      << _soln._resid << endl
      << "# Physical time" << endl
      << setw(23) << scientific << setprecision(16)
      << _soln._current_time << endl
      << "# Dimensions" << endl
      << setw(8) << ncells << setw(8) << n_cvar << endl
      << "# Mixture density" << endl;
  unsigned int cell = 0;
  unsigned int ll = 5;
  while(cell < ncells){
    Ouf <<  setw(23) << scientific << setprecision(16)
	<< _soln._rhof[cell++];
    if(!(cell%ll))
      Ouf << endl;
  }
  if(cell%ll)
    Ouf << endl;
  for(int i = 0;i < 3;i++){
    Ouf << "# Mixture " 
	<< (i == 0 ? "x-momentum" :
	    i == 1 ? "y-momentum" : "z-momentum") 
	<< endl; 
    cell = 0;
    while(cell < ncells){
      Ouf << setw(23) << scientific << setprecision(16) 
	  << _soln._rhovf[(cell++)*3+i];
      if(!(cell%ll))
	Ouf << endl;
    }
    if(cell%ll)
      Ouf << endl;
  }
  Ouf << "# Mixture total internal energy" << endl;
  cell = 0;
  while(cell < ncells){
    Ouf << setw(23) << scientific << setprecision(16) 
	<< _soln._rhoEf[cell++];
    if(!(cell%ll))
      Ouf << endl;
  }
  if(cell%ll)
    Ouf << endl;
  Ouf << "# End" << endl;
  Ouf.close();
  WriteGSPASCII(pre,n,t);
  if(_debug && _out)
    *_out << "FluRegion::WriteFluSolnASCII: Exit\n";
  return(true);
}

bool
FluRegion::WriteGSPASCII(const string &prefix,unsigned int n,
			 double t,bool unsteady)
{
  if(_debug && _out)
    *_out << "FluRegion::WriteGSPASCII: Enter\n";
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename + ".";
  }
  else
    pre += ".";
  ostringstream Ostr;
  Ostr << pre << "gspa_" << setw(5) << setfill('0') << _id << "_";
  string filebase(Ostr.str());
  Ostr.str("");
  Ostr << setw(6) << setfill('0') << n;
  string iterstring(Ostr.str());
  Ostr.clear();
  Ostr.str("");
  Ostr << scientific << setprecision(5) << t;
  string timestring(Ostr.str());
  timestring.replace(7,1,"E");
  ofstream Ouf;
  if(!unsteady) filebase+=iterstring;
  else filebase+=timestring;
  Ouf.open(filebase.c_str());
  if(!Ouf)
    return false;
  Ouf << "# ROCFLU grid speeds file" << endl
      << "# Precision and range" << endl
      << setw(8) << 15 << setw(8) << 307 << endl
      << "# Physical time" << endl
      << setw(23) << scientific << setprecision(16) 
      << t << endl
      << "# Dimensions" << endl
      << setw(8) << _nvface << endl
      << "# Grid speeds" << endl;
  unsigned int face = 0;
  unsigned int ll = 5;
  bool gsp = (_soln._gsp.size() != 0);
  if(_debug && _out){
    *_out << "FluRegion::WriteGSPASCII: nvface = " << _nvface << "\n";
    if(!gsp) 
      *_out 
	<< "FluRegion::WriteGSPASCII:no gridspeeds for the volume faces.\n";
  }
  while(face < _nvface){
    Ouf << setw(23) << scientific << setprecision(16)
	<< (gsp ? _soln._gsp[face] : 0);
    face++;
    if(!(face%ll))
      Ouf << endl;
  }
  if(face%ll)
    Ouf << endl;
  if(_debug && _out)
    *_out << "FluRegion::WriteGSPASCII: npatches = " << _patches.size() 
	  << "\n";
  vector<FluPatch>::iterator fpi = _patches.begin();
  while(fpi != _patches.end()){
    FluPatch &fp = *fpi++;
    gsp = (fp._soln._gsp.size() != 0);
    unsigned int nfaces = fp._ntri+fp._nquad-(fp._ngtri+fp._ngquad);
    if(_debug && _out){
    *_out << "FluRegion::WriteGSPASCII: npatch_faces = " << nfaces << "\n";
      if(!gsp) 
	*_out << "FluRegion::WriteGSPASCII: no gridspeeds for patch faces." 
	     << "\n";
    }
    face = 0;
    while(face < nfaces){
      Ouf << setw(23) << scientific << setprecision(16)
	  << (gsp ? fp._soln._gsp[face] : 0);
      face++;
      if(!(face%ll))
	Ouf << endl;
    }
    if(face%ll)
      Ouf << endl;
  }
  Ouf << "# End" << endl;
  Ouf.close();  
  if(_debug && _out)
    *_out << "FluRegion::WriteGSPASCII: Exit\n";
  return true;
}

bool
FluRegion::ReadGSPASCII(const string &prefix,unsigned int n,
			double t,bool unsteady)
{
  if(_debug && _out)
    *_out << "FluRegion::ReadGSPASCII: Enter\n";
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  ostringstream Ostr;
  Ostr << pre << ".gspa_" << setw(5) << setfill('0') << _id << "_";
  string filebase(Ostr.str());
  Ostr.str("");
  Ostr << setw(6) << setfill('0') << n;
  string iterstring(Ostr.str());
  Ostr.clear();
  Ostr.str("");
  Ostr << scientific << setprecision(5) << t;
  string timestring(Ostr.str());
  timestring.replace(7,1,"E");
  ifstream Inf;
  if(!unsteady)
    Inf.open((filebase+iterstring).c_str());
  else
    Inf.open((filebase+timestring).c_str());
  if(!Inf)
    return false;
  string line;
  SkipLines(Inf,6);
  Inf >> _nvface;
  _soln._gsp.resize(_nvface,0.0);
  SkipLines(Inf,2);
  unsigned int face = 0;
  while(face < _nvface)
    Inf >> _soln._gsp[face++];
  vector<FluPatch>::iterator fpi = _patches.begin();
  while(fpi != _patches.end()){
    FluPatch &fp = *fpi++;
    unsigned int nfaces = fp._ntri+fp._nquad-(fp._ngtri+fp._ngquad);
    fp._soln._gsp.resize(nfaces,0.0);
    face = 0;
    while(face < nfaces)
      Inf >> fp._soln._gsp[face++];
    fp.InitSurfaceSoln(pre,false);
  }
  Inf.close();
  return true;
}


bool 
FluRegion::BuildPatchMapping(map<unsigned int,unsigned int> &patch_mapping,
			     const string &prefix)
{
  // Build the patch mapping from the hand produced patch mapping file
  ostringstream Ostr;
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  Ostr << pre << ".cgi";
  ifstream bcfile;
  bcfile.open(Ostr.str().c_str());
  if(!bcfile)
    return(false);
  //  map< unsigned int,unsigned int > patch_mapping;
  unsigned int nggpatch_mapping_lines;
  bcfile >> _npatches_total >> nggpatch_mapping_lines;
  unsigned int patch_line = 0;
  while(patch_line < nggpatch_mapping_lines){
    unsigned int low, high, flupatch;
    bcfile >> low >> high >> flupatch;
    unsigned int p = low;
    while(p <= high)
      patch_mapping[p++] = flupatch;
    patch_line++;
  }
  if(_debug && _out){
    *_out << "FluRegion::BuildPatchMapping: Patch Mapping:\n";
    map<unsigned int,unsigned int>::iterator mi = patch_mapping.begin();
    while(mi != patch_mapping.end()){
      *_out << "FluRegion::BuildPatchMapping: (" << mi->first << "," 
	   << mi->second << ")\n";
      mi++;
    }
  }
  bcfile.close();
  return(true);
}  

bool
FluRegion::ReadRegionASCII(const string &prefix,unsigned int region_id,
			   unsigned int n,double t,bool unsteady)
{
  if(_debug && _out)
    *_out << "FluRegion::ReadRegionASCII:: Reading Rocflu Region "
	 << "from Rocflu native files...\n";
  string pre(prefix);
  if(pre.empty()){
    if(_casename.empty())
      pre = "./";
    else
      pre = _casepath + "/" + _casename;
  }
  _id = region_id;
  if(!ReadFluDIM(pre,t,unsteady))
    return(false);
  if(!ReadFluCOM(pre)){
    if(_debug && _out)
      *_out << "FluRegion::ReadRegionASCII: ReadFluCOM failed.\n";
    return(false);
  }
  if(!ReadFluGridASCII(pre,t,unsteady)){
    if(_debug && _out)
      *_out << "FluRegion::ReadRegionASCII: ReadFluGridASCII failed.\n";
    return(false);
  }
  if(!ReadFluSolnASCII(pre,n,t,unsteady)){
    if(_debug && _out)
      *_out << "FluRegion::ReadRegionASCII: ReadFluSolnASCII failed.\n";
    return(false);
  }
  if(_debug && _out)
    *_out << "FluRegion::ReadRegionASCII:: Done reading Rocflu Region "
	 << "from Rocflu native files.\n";
  return(true);
}

void
FluRegion::PopRemBordIndMPI()
{
  // Assumes rank = partition - 1
  int *local_index = NULL;
  int *remote_index = NULL;
  MPI_Request *req = NULL;
  int ind = 0;
  MPI_Status *status = NULL;
  if(_borders.size() > 0){
    local_index = new int [_borders.size()];
    req = new MPI_Request [2*_borders.size()];
    remote_index = new int [_borders.size()];
    status = new MPI_Status [_borders.size()*2];
  }
  vector<FluBorder>::iterator fbi = _borders.begin();
  while(fbi != _borders.end()){
    int remote_proc = fbi->_rpart - 1;
    local_index[ind] = ind + 1;
    MPI_Isend(&local_index[ind],1,MPI_INT,remote_proc,remote_proc,
	      MPI_COMM_WORLD,&req[ind]);
    fbi++;
    ind++;
  }
  int ind2 = 0;
  fbi = _borders.begin();
  while(fbi != _borders.end()){
    int remote_proc = fbi->_rpart - 1;
    MPI_Irecv(&remote_index[ind2++],1,MPI_INT,remote_proc,_id-1,
	      MPI_COMM_WORLD,&req[ind]);
    fbi++;
    ind++;
  }
  if(_borders.size() > 0)
    MPI_Waitall(_borders.size()*2,req,status);
}

//#ifdef _ROCSTAR_X_
bool
FluRegion::ComRemeshInitData(const string &wname,double *cell_data,
			     int nval_cells, double *node_data,int nval_nodes)
{
  _ntet = _tetconn.size()/4;
  _nhex = _hexconn.size()/8;
  _npris = _prisconn.size()/6;
  _npyr = _pyrconn.size()/5;
  _nnodes = _nc.size()/3;
  unsigned int n_elem = nelem();
  _soln.Resize(n_elem,0,_nnodes);
  _soln._unsteady = true;
  unsigned int ind = 0;
  unsigned int cdind = 0;
  unsigned int vind = 0;
  while(ind < n_elem){
    _soln._rhof[ind] = cell_data[cdind++];
    _soln._rhovf[vind++] = cell_data[cdind++];
    _soln._rhovf[vind++] = cell_data[cdind++];
    _soln._rhovf[vind++] = cell_data[cdind++];
    _soln._rhoEf[ind] = cell_data[cdind++];
    _soln._pf[ind] = cell_data[cdind++];
    _soln._Tf[ind] = cell_data[cdind++];
    _soln._af[ind] = cell_data[cdind++];
    ind++;
  }
  ind = 0;
  while(ind < _nnodes){
    _soln._disp[ind*3] = node_data[ind*3];
    _soln._disp[ind*3+1] = node_data[ind*3+1];
    _soln._disp[ind*3+2] = node_data[ind*3+2];
    ind++;
  }
  return(true);
}

bool
FluRegion::RegisterVolumeSoln(bool all)
{
  Create_com_volsoln("rhof",_soln._rhof,1,"kg/m^3");
  Create_com_volsoln("rhovf",_soln._rhovf,3,"kg/(m^2s)");
  Create_com_volsoln("rhoEf",_soln._rhoEf,1,"J/m^3");
  COM_new_attribute((volume_window+".gs"),'p',COM_DOUBLE,1,"m/s");
  COM_set_size((volume_window+".gs"),pane_id,_nvface);
  if(_soln._gsp.size() != 0)
    COM_set_array((volume_window+".gs"),pane_id,&(_soln._gsp[0]));
  Create_com_volsoln("pf",_soln._pf,1,"Pa");
  Create_com_volsoln("Tf",_soln._Tf,1,"K");
  Create_com_volsoln("af",_soln._af,1,"m/s");
  return(true);
}

bool
FluPatch::RegisterSoln(const string &wname,bool all)
{
  if(_debug && _out)
    *_out << "FluPatch(" << _id << ")::RegisterSoln: "
	  << "Registering pane " << pane_id << "." << endl;
  if(_soln._gsp.size() != 0){
    COM_set_array((wname+".gs"),pane_id,&(_soln._gsp[0]));
  }
  bool do_t0 = (t0_coords.size() != 0);
  COM_set_size((wname+".patchNo"),pane_id,1);
  COM_set_array((wname+".patchNo"),pane_id,&_id);
  if(do_t0){
    COM_set_size((wname+"0"+".bcflag"),pane_id,1);
    COM_set_array((wname+"0"+".bcflag"),pane_id,&_bcflag);
  }
  COM_set_size((wname+".bcflag"),pane_id,1);
  COM_set_array((wname+".bcflag"),pane_id,&_bcflag);
  COM_set_size((wname+".constr_type"),pane_id,1);
  COM_set_array((wname+".constr_type"),pane_id,&_constr_type);
  if(all){
    COM_set_array((wname+".du_alp"),pane_id,&_soln._du_alp[0]);
    Create_com_surfsoln(wname,"rhofvf_alp",_soln._rhofvf_alp,3,"kg/(m^2s)");
    Create_com_surfsoln(wname,"nf_alp",_soln._nf_alp,3,"");
    Create_com_surfsoln(wname,"rhof_alp",_soln._rhof_alp,1,"kg/m^3");
    Create_com_surfsoln(wname,"pf",_soln._pf,1,"Pa");
    Create_com_surfsoln(wname,"qc",_soln._qc,1,"W/m^2");
    Create_com_surfsoln(wname,"qr",_soln._qr,1,"W/m^2");
    Create_com_surfsoln(wname,"tf",_soln._tf,3,"Pa");
    Create_com_surfsoln(wname,"Tb_alp",_soln._Tb_alp,1,"K");
    Create_com_surfsoln(wname,"mdot_alp",_soln._mdot_alp,1,"kg/(m^2s)");
    Create_com_surfsoln(wname,"Tflm_alp",_soln._Tflm_alp,1,"K");
    COM_set_array((wname+".bflag"),pane_id,&(_soln._bflag[0]));
  }
  return(true);
}

bool
FluRegion::RegisterSurfaceSoln(bool all)
{
 
  if(_debug && _out)
    *_out << "FluRegion(" << _id << ")::RegisterSurfaceSoln Enter" << endl;
  COM_new_attribute((surface_window+".gs"),'e',COM_DOUBLE,1,"m/s");
  if(_t0_coords.size() != 0)
    COM_new_attribute((surface_window+"0"+".bcflag"),'p',COM_INT,1,"");
  COM_new_attribute((surface_window+".bcflag"),'p',COM_INT,1,"");
  COM_new_attribute((surface_window+".constr_type"),'p',COM_INT,1,"");
  COM_new_attribute((surface_window+".patchNo"),'p',COM_INT,1,"");
  if(all){
    COM_new_attribute((surface_window+".du_alp"),'n',COM_DOUBLE,3,"m");
    COM_new_attribute((surface_window+".bflag"),'e',COM_INT,1,"");
    COM_new_attribute(surface_window+".rhofvf_alp",'e',
		      COM_DOUBLE,3,"kg/(m^2s)");
    COM_new_attribute(surface_window+".nf_alp",'e',COM_DOUBLE,3,"");
    COM_new_attribute(surface_window+".rhof_alp",'e',COM_DOUBLE,1,"kg/m^3");
    COM_new_attribute(surface_window+".pf",'e',COM_DOUBLE,1,"Pa");
    COM_new_attribute(surface_window+".qc",'e',COM_DOUBLE,1,"W/m^2");
    COM_new_attribute(surface_window+".qr",'e',COM_DOUBLE,1,"W/m^2");
    COM_new_attribute(surface_window+".tf",'e',COM_DOUBLE,3,"Pa");
    COM_new_attribute(surface_window+".Tb_alp",'e',COM_DOUBLE,1,"K");
    COM_new_attribute(surface_window+".mdot_alp",'e',COM_DOUBLE,1,"kg/(m^2s)");
    COM_new_attribute(surface_window+".Tflm_alp",'e',COM_DOUBLE,1,"K");
  }
  vector<FluPatch>::iterator fpi = _patches.begin();
  if(_debug && _out)
    *_out << "FluRegion(" << _id << ")::RegisterSurfaceSoln: "
	  << "Registering patches..." << endl;
  while(fpi != _patches.end()){
    FluPatch &fp = *fpi++;
    fp.RegisterSoln(surface_window,all);
  }
  if(_debug && _out)
    *_out << "FluRegion(" << _id << ")::RegisterSurfaceSoln: "
	  << "Registering patche done." << endl;
  return(true);
}

bool
FluRegion::RegisterFluSurfaceMesh()
{
  COM_new_attribute((surface_window+".t3g:real"),'p',COM_INT,3,"");
  COM_new_attribute((surface_window+".t3g:virtual"),'p',COM_INT,3,"");      
  COM_new_attribute((surface_window+".q4g:real"),'p',COM_INT,4,"");
  COM_new_attribute((surface_window+".q4g:virtual"),'p',COM_INT,4,"");      

  vector<FluPatch>::iterator fpi = _patches.begin();
  while(fpi != _patches.end()){
    FluPatch &fp = *fpi++;
    unsigned int nreal = fp._triconn.size()/3 - fp._ngtri;
    if(nreal > 0){
      COM_set_size((surface_window+".t3g:real"),fp.pane_id,nreal);    
      COM_set_array((surface_window+".t3g:real"),fp.pane_id,
		    &(fp._triconn[0]),3);
    }
    nreal = fp._quadconn.size()/4 - fp._ngquad;
    if(nreal > 0){
      COM_set_size((surface_window+".q4g:real"),fp.pane_id,nreal);
      COM_set_array((surface_window+".q4g:real"),fp.pane_id,
		    &(fp._quadconn[0]),4);
    }
    if(fp._ngtri > 0){
      nreal = fp._triconn.size()/3 - fp._ngtri;
      COM_set_size((surface_window+".t3g:virtual"),fp.pane_id,
		   fp._ngtri,fp._ngtri);
      COM_set_array((surface_window+".t3g:virtual"),fp.pane_id,
		    &(fp._triconn[nreal*3]),3);
    }
    if(fp._ngquad > 0){
      nreal = fp._quadconn.size()/4 - fp._ngquad;
      COM_set_size((surface_window+".q4g:virtual"),fp.pane_id,
		   fp._ngquad,fp._ngquad);
      COM_set_array((surface_window+".q4g:virtual"),fp.pane_id,
		    &(fp._quadconn[nreal*4]),4);
    }
  }
  return(true);
}
//#endif

istream &
SkipLines(istream &ioo,unsigned int n)
{
  string line;
  unsigned int l = 0;
  while(l++ < n)
    getline(ioo,line);
  return(ioo);
}






