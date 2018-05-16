#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <algorithm>

#include "Mesh.H"
#include "Profiler.H"
#include "Global.H"
#include "primitive_utilities.H"

class T3DGeomEnt : public std::pair<unsigned int,unsigned int>
{
public:
  T3DGeomEnt() : std::pair<unsigned int,unsigned int>()
  {};
  T3DGeomEnt(const T3DGeomEnt &ing)
  {
    this->first = ing.first;
    this->second = ing.second;
  };
  const T3DGeomEnt &operator=(const T3DGeomEnt &ing)
  {
    this->first = ing.first;
    this->second = ing.second;
    return(*this);
  };
  bool operator<(const T3DGeomEnt &ge) const { return (first < ge.first ? true  :
						       first > ge.first ? false : second < ge.second); };
  bool operator==(const T3DGeomEnt &ge) const { return (first == ge.first ? second==ge.second : false); };
};


int
main(int argc,char *argv[])
{
  std::string line;
  std::istringstream Istr;
  unsigned int mesh_type = 0;
  unsigned int elem_order = 0;
  unsigned int nnodes  = 0;
  unsigned int nedges  = 0;
  unsigned int ntris   = 0;
  unsigned int nquads  = 0;
  unsigned int ntets   = 0;
  unsigned int nhexs   = 0;
  unsigned int npyrs   = 0;
  unsigned int npriss  = 0;


  // Read Block 1
  std::getline(std::cin,line);
  Istr.str(line);
  Istr >> mesh_type >> elem_order;
  std::getline(std::cin,line);
  Istr.str(line);
  Istr >> nnodes >> nedges;
  switch(mesh_type){
  case 3:
    Istr >> ntris >> ntets;
    break;
  case 4:
    Istr >> nquads >> nhexs;
    break;
  case 7:
    Istr >> ntris >> nquads >> ntets >> npyrs >> npriss >> nhexs;
    break;
  default:
    std::cerr << "Fatal Error: Unsupported mesh type, "
	      << mesh_type << "." << std::endl;
    return(1);
    break;

  }
  std::getline(std::cin,line);

  // Calculate some mesh parameters
  unsigned int nelems;
  unsigned int n_2d_elem = ntris + nquads;
  unsigned int n_3d_elem = ntets + npyrs + npriss + nhexs;
  if(n_2d_elem > 0 && n_3d_elem > 0){
    std::cerr << "Fatal Error: Mixing 2d and 3d elements not supported at this time."
	      << std::endl;
    return(1);
  }
  nelems = n_2d_elem > 0 ? n_2d_elem : n_3d_elem;

  std::cout << nnodes << std::endl;
  // Read Block 2 (nodes)
  // Open a file to store properties
  std::ofstream PropOut;
  std::list<T3DGeomEnt> geomlist;
  std::list<int> nodal_properties;
  std::list<int> elemental_properties;
  std::vector<T3DGeomEnt> nodal_geoms(nnodes);
  std::vector<T3DGeomEnt> elem_geoms(nelems);
  for(unsigned int i = 0;i < nnodes;i++){
    T3DGeomEnt ge;
    std::getline(std::cin,line);
    int superfluous_id;
    Istr.clear();
    Istr.str(line);
    Istr >> superfluous_id;
    double x,y,z;
    Istr >> x >> y >> z;
    std::cout << x << "\t" << y << "\t" << z << std::endl;
    int property = 0;
    Istr >> nodal_geoms[i].first >> nodal_geoms[i].second >> property;
    //    std::getline(Istr,line);
    //    std::string::size_type xpos = line.find_first_not_of(" ");
    //    line = line.substr(xpos);
    //    PropOut << line << std::endl;
    nodal_properties.push_back(property);
  }
  // add a blank line between nodal and elemental properties
  //  PropOut << std::endl;
  std::getline(std::cin,line);
  std::cout << nelems << std::endl;

  // Skip Block 3 (edges)
  if(nedges > 0){
    std::cerr << "Warning: skipping edges." << std::endl;
    for(unsigned int i = 0;i < nedges;i++)
      std::getline(std::cin,line);
    std::getline(std::cin,line);
  }

  // Read Block 4 (triangle elements)
  if(mesh_type != 4 && ntris > 0){
    for(unsigned int i = 0;i < ntris;i++){
      std::getline(std::cin,line);
      std::string::size_type xpos = line.find_first_not_of(" ");
      line = line.substr(xpos);
      Istr.clear();
      Istr.str(line);
      unsigned int superfluous_id;
      Istr >> superfluous_id;
      unsigned int node_id;
      for(unsigned int j = 0;j < elem_order*3;j++){
	Istr >> node_id;
	std::cout << node_id;
	if(j != (elem_order*3 - 1)) std::cout << "\t";
      }
      std::cout << std::endl;
      int property = 0;
      Istr >> elem_geoms[i].first >> elem_geoms[i].second >> property;
      elemental_properties.push_back(property);
      //      std::getline(Istr,line);
      //      std::string::size_type x = line.find_first_not_of(" ");
      //      line = line.substr(x);
      //      PropOut << line << std::endl;
    }
    std::getline(std::cin,line);
  }
  Mesh::IndexType egoffset = ntris;

  // Read Block 5 (quad elements)
  if(mesh_type != 3 && nquads > 0){
    for(unsigned int i = 0;i < nquads;i++){
      std::getline(std::cin,line);
      std::string::size_type xpos = line.find_first_not_of(" ");
      line = line.substr(xpos);
      Istr.clear();
      Istr.str(line);
      unsigned int superfluous_id;
      Istr >> superfluous_id;
      unsigned int node_id;
      for(unsigned int j = 0;j < elem_order*4;j++){
	Istr >> node_id;
	std::cout << node_id;
	if(j != (elem_order*4 - 1)) std::cout << "\t";
      }
      std::cout << std::endl;
      int property = 0;
      Istr >> elem_geoms[i+egoffset].first >> elem_geoms[i+egoffset].second >> property;
      elemental_properties.push_back(property);
    }
    std::getline(std::cin,line);
  }
  egoffset += nquads;
  // Read Block 6 (tet elements)
  if(mesh_type != 4 && ntets > 0){
    unsigned int elem_size = elem_order==1 ? 4 : 10;
    for(unsigned int i = 0;i < ntets;i++){
      std::getline(std::cin,line);
      std::string::size_type xpos = line.find_first_not_of(" ");
      line = line.substr(xpos);
      Istr.clear();
      Istr.str(line);
      unsigned int superfluous_id;
      Istr >> superfluous_id;
      unsigned int node_id;
      for(unsigned int j = 0;j < elem_size;j++){
	Istr >> node_id;
	std::cout << node_id;
	if(j != (elem_size - 1)) std::cout << "\t";
      }
      std::cout << std::endl;
      int property = 0;
      Istr >> elem_geoms[i+egoffset].first >> elem_geoms[i+egoffset].second >> property;
      elemental_properties.push_back(property);
    }
    std::getline(std::cin,line);
  }
  egoffset += ntets;
  // Read Block 7 (pyr elements)
  if(mesh_type == 7 && npyrs > 0){
    unsigned int elem_size = elem_order==1 ? 5 : 8;
    for(unsigned int i = 0;i < npyrs;i++){
      std::getline(std::cin,line);
      std::string::size_type xpos = line.find_first_not_of(" ");
      line = line.substr(xpos);
      Istr.clear();
      Istr.str(line);
      unsigned int superfluous_id;
      Istr >> superfluous_id;
      unsigned int node_id;
      for(unsigned int j = 0;j < elem_size;j++){
	Istr >> node_id;
	std::cout << node_id;
	if(j != (elem_size - 1)) std::cout << "\t";
      }
      std::cout << std::endl;
      int property = 0;
      Istr >> elem_geoms[i+egoffset].first >> elem_geoms[i+egoffset].second >> property;
      elemental_properties.push_back(property);
    }
    std::getline(std::cin,line);
  }

  egoffset += npyrs;

  // Read Block 8 (pris elements)
  if(mesh_type == 7 && npriss > 0){
    unsigned int elem_size = elem_order==1 ? 6 : 15;
    for(unsigned int i = 0;i < npriss;i++){
      std::getline(std::cin,line);
      std::string::size_type xpos = line.find_first_not_of(" ");
      line = line.substr(xpos);
      Istr.clear();
      Istr.str(line);
      unsigned int superfluous_id;
      Istr >> superfluous_id;
      unsigned int node_id;
      for(unsigned int j = 0;j < elem_size;j++){
	Istr >> node_id;
	std::cout << node_id;
	if(j != (elem_size - 1)) std::cout << "\t";
      }
      std::cout << std::endl;
      int property = 0;
      Istr >> elem_geoms[i+egoffset].first >> elem_geoms[i+egoffset].second >> property;
      elemental_properties.push_back(property);
    }
    std::getline(std::cin,line);
  }
  egoffset += npriss;
  // Read Block 9 (hex elements)
  if(mesh_type != 3 && nhexs > 0){
    unsigned int elem_size = elem_order==1 ? 8 : 20;
    for(unsigned int i = 0;i < nhexs;i++){
      std::getline(std::cin,line);
      std::string::size_type xpos = line.find_first_not_of(" ");
      line = line.substr(xpos);
      Istr.clear();
      Istr.str(line);
      unsigned int superfluous_id;
      Istr >> superfluous_id;
      unsigned int node_id;
      for(unsigned int j = 0;j < elem_size;j++){
	Istr >> node_id;
	std::cout << node_id;
	if(j != (elem_size - 1)) std::cout << "\t";
      }
      std::cout << std::endl;
      int property = 0;
      Istr >> elem_geoms[i+egoffset].first >> elem_geoms[i+egoffset].second >> property;
      elemental_properties.push_back(property);
    }
    std::getline(std::cin,line);
  }
  //  PropOut.close();
  geomlist.resize(0);
  // Now sort out the geometry information - real slow
  std::vector<T3DGeomEnt>::iterator gi = nodal_geoms.begin();
  while(gi != nodal_geoms.end()){
    std::list<T3DGeomEnt>::iterator gli = std::find(geomlist.begin(),geomlist.end(),*gi);
    if(gli == geomlist.end()){
      geomlist.push_back(*gi);
    }
    gi++;
  }
  gi = elem_geoms.begin();
  while(gi != elem_geoms.end())
    {
      std::list<T3DGeomEnt>::iterator gli = std::find(geomlist.begin(),geomlist.end(),*gi);
      if(gli == geomlist.end())
	geomlist.push_back(*gi);
      gi++;
    }
  geomlist.sort();
  std::list<T3DGeomEnt>::iterator gli = geomlist.begin();
  std::cerr << "GEOMETRIES:" << std::endl;
  while(gli != geomlist.end())
    {
      std::cerr << "Dimension: " << gli->first << ", Id = " << gli->second << std::endl;
      gli++;
    }
  Mesh::IndexType ngeoms = geomlist.size();
  Mesh::IndexType nvertex = 0;
  Mesh::IndexType ncurve  = 0;
  Mesh::IndexType nsurf   = 0;
  Mesh::IndexType nregion = 0;
  Mesh::IndexType npatch  = 0;
  //  Mesh::IndexType nshell  = 0;
  //  Mesh::IndexType ninter  = 0;
  std::vector<Mesh::GeometricEntity> geometries(ngeoms);
  std::vector<T3DGeomEnt> t3d_geovec(ngeoms);
  gli = geomlist.begin();
  while(gli != geomlist.end()){
    ubyte dim = (gli->first == 1 ? 0 :
		 gli->first == 2 ? 1 :
		 gli->first == 3 ? 2 :
		 gli->first == 4 ? 3 :
		 gli->first == 5 ? 2 :
		 gli->first == 6 ? 2 :
		 gli->first);
    if(dim > 3){
      std::cerr << "Error: Cannot deal with interface mesh type"
		<< std::endl;
      return(1);
    }
    std::ostringstream Ostr;
    if(dim == 0)
      Ostr << "Vertex-" << ++nvertex;
    if(dim == 1)
      Ostr << "Curve-" << ++ncurve;
    if(dim == 2)
      Ostr << "Surface-" << ++nsurf;
    if(dim == 3)
      Ostr << "Region-" << ++nregion;
    Mesh::IndexType index = nvertex+ncurve+nsurf+nregion-1;
    t3d_geovec[index] = *gli;
    geometries[index].first = Ostr.str();
    geometries[index].second = dim;
    geometries[index]._collections.resize((unsigned int)dim+1);
    gli++;
  }
  geomlist.resize(0);
  std::cerr << "Geometries created, populating..." << std::endl;
  gi = nodal_geoms.begin();
  Mesh::IndexType nid = 1;
  std::cerr << "     populating nodal geoms..." << std::endl;
  PropOut.open("t3d_properties");
  std::list<int>::iterator pi = nodal_properties.begin();
  while(gi != nodal_geoms.end()){
    Mesh::IndexType gindex = std::find(t3d_geovec.begin(),t3d_geovec.end(),*gi++) -
      t3d_geovec.begin();
    geometries[gindex]._collections[0].push_back(nid++);
    PropOut << gindex+1 << "\t" << *pi++ << std::endl;
  }
  // blank line between nodal and elemental properties
  PropOut << std::endl;
  gi = elem_geoms.begin();
  pi = elemental_properties.begin();
  std::cerr << "     populating elemental geoms..." << std::endl;
  Mesh::IndexType eid = 1;
  while(gi != elem_geoms.end()){
    Mesh::IndexType gindex = std::find(t3d_geovec.begin(),t3d_geovec.end(),*gi++) -
      t3d_geovec.begin();
    geometries[gindex]._collections[geometries[gindex].second].push_back(eid++);
    PropOut << gindex+1 << "\t" << *pi++ << std::endl;
  }
  PropOut.close();
  std::ofstream GeomOut;
  GeomOut.open("t3d_geometries");
  GeomOut << geometries.size() << std::endl;
  for(Mesh::IndexType gindex = 0;gindex < geometries.size();gindex++)
    GeomOut << geometries[gindex].first << "\t"
	    << (uint)geometries[gindex].second << std::endl;
  for(Mesh::IndexType gindex = 0;gindex < geometries.size();gindex++)
    GeomOut << geometries[gindex] << std::endl;
  GeomOut.close();
  return(0);
}
