#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "TRAIL.H"
#include "roccom.h"

void CheckCoordinates(double *coords,int number_of_points,double tol)
{
  for(int i = 0; i < number_of_points*3;i++)
    {
      if(std::fabs(coords[i]) < tol && (coords[i] != 0.0)){
	std::cout << "WRLSURF2HDF> WARNING: node " << (i/3)+1 << "/" 
		  << number_of_points << " "   
		  << (((i%3)==0 ? "x " : ((i%3) < 2) ? "y " : "z "))
		  << "coordinate = " << coords[i] << std::endl;
	coords[i] = 0.0;
      }
    }
}

int
main(int argc,char *argv[])
{
  COM_init(&argc,&argv);
  std::string line;
  std::getline(std::cin,line);
  std::string::size_type x = 0;
  Mesh::NodalCoordinates nc;
  std::vector<double> coords;
  //  nc.resize(0);
  x = line.find("point [");
  while(x == std::string::npos  && !std::cin.eof()){
    std::getline(std::cin,line);
    x = line.find("point [");
  }
  if(std::cin.eof()){
    std::cerr << "File error 1." << std::endl;
    exit(1);
  }
  std::getline(std::cin,line);
  x = line.find("]");
  unsigned int count = 1;
  while(x == std::string::npos && !std::cin.eof()){
    std::istringstream Istr(line);
    double xc, yc, zc;
    Istr >> xc >> yc >> zc;
    
    coords.push_back(xc);
    coords.push_back(yc);
    coords.push_back(zc);
    std::getline(std::cin,line);
    x = line.find("]");
  }
  if(std::cin.eof()){
    std::cerr << "File error 2." << std::endl;
    exit(1);
  }
  unsigned int number_of_nodes = coords.size()/3;
  CheckCoordinates(&(coords[0]),number_of_nodes,1.0e-9);
  nc.init(number_of_nodes,&(coords[0]));
  getline(std::cin,line);
  x = line.find("coordIndex [");
  while(x == std::string::npos && !std::cin.eof()){
    getline(std::cin,line);
    x = line.find("coordIndex [");
  }
  if(std::cin.eof()){
    std::cerr << "File error 3." << std::endl;
    exit(1);
  }
  Mesh::Connectivity ec;
  std::vector<Mesh::IndexType> element(3);
  std::getline(std::cin,line);
  x = line.find("]");
  while(x == std::string::npos && !std::cin.eof()){
    std::istringstream Istr(line);
    std::string delim;
    int dummy;
    Istr >> element[0] >> delim >> element[1] >> delim >> element[2] >> delim >> dummy;
    element[0]+=1;
    element[1]+=1;
    element[2]+=1;
    ec.AddElement(element);
    std::getline(std::cin,line);
    x = line.find("]");
  }
  if(std::cin.eof()){
    std::cerr << "File error 4." << std::endl;
    exit(1);
  }
  ec.Sync();
  unsigned int nvert = nc.Size();
  unsigned int ntri  = ec.Nelem();
  std::cout << "Number of vertices:  " << nvert << std::endl
	    << "Number of triangles: "  << ntri << std::endl;
  TRAIL_SurfaceMesh2Window("tempwin",1,nc,ec);
  TRAIL_WriteWindow("tempwin",".","","",0,0,MPI_COMM_NULL,NULL);
  return(0);
}
