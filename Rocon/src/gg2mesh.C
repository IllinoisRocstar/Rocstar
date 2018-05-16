#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Mesh.H"
#include "Profiler.H"
#include "Global.H"
#include "primitive_utilities.H"

int
main(int argc,char *argv[])
{
  std::string line;
  bool done = false;
  while(!done){
    std::getline(std::cin,line);
    std::string::size_type x = line.find("Nodes");
    if(x != std::string::npos)
      done = true;
  }
  Mesh::NodalCoordinates nc;
  std::cin >> nc;
  std::cout << nc << std::endl;
  done = false;
  while(!done){
    std::getline(std::cin,line);
    std::string::size_type x = line.find("Boundary");
    if(x != std::string::npos)
      done = true;
  }
  Mesh::IndexType nfaces;
  std::getline(std::cin,line);
  std::istringstream Dum(line);
  Dum >> nfaces;
  for(Mesh::IndexType i = 0;i<nfaces;i++)
    getline(std::cin,line);
  done = false;
  while(!done){
    std::getline(std::cin,line);
    std::string::size_type x = line.find("Elements");
    if(x != std::string::npos)
      done = true;
  }
  done = false;
  Mesh::Connectivity ec;
  while(!done){
    std::getline(std::cin,line);
    std::string::size_type x = line.find("Variables");
    if(x == std::string::npos){
      std::istringstream Istr(line);
      Mesh::IndexType dummy;
      Istr >> dummy >> dummy;
      std::vector<Mesh::IndexType> nodes;
      while(Istr >> dummy)
	nodes.push_back(dummy);
      ec.AddElement(nodes);
    }
    else
      done = true;
  }
  std::cout << ec << std::endl;
  return(0);
}
