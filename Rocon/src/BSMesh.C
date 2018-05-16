///
/// \file
/// \ingroup support
/// \brief Mesh stuff implementation
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include "BSMesh.H"

namespace Mesh {
  int GenerateCartesianGrid(Mesh::NodalCoordinates &nc,
			    Mesh::BSExtent<Mesh::IndexType> &gridextent,
			    std::vector<Mesh::IndexType> &gridsizes,
			    GeoPrim::CBox &box)
  {
    if(gridsizes.size() != 3)
      return 1;
    // sizes is ncells
    Mesh::IndexType ncells_x = gridsizes[0]; 
    Mesh::IndexType ncells_y = gridsizes[1];
    Mesh::IndexType ncells_z = gridsizes[2];
    Mesh::IndexType nnodes_x = ncells_x + 1;
    Mesh::IndexType nnodes_y = ncells_y + 1;
    Mesh::IndexType nnodes_z = ncells_z + 1;
    std::vector<Mesh::IndexType> flatextent(6);
    flatextent[0] = 1;
    flatextent[1] = nnodes_x;
    flatextent[2] = 1;
    flatextent[3] = nnodes_y;
    flatextent[4] = 1;
    flatextent[5] = nnodes_z;
    gridextent.Init(flatextent);
    unsigned int npoints = flatextent[1]*flatextent[3]*flatextent[5];
    nc.init(npoints);
    std::vector<double> xcoords(ncells_x+1);
    std::vector<double> ycoords(ncells_y+1);
    GeoPrim::CPoint p(box.P2() - box.P1());
    double h_x = (ncells_x > 0 ? p.x()/ncells_x : 0);
    double h_y = (ncells_y > 0 ? p.y()/ncells_y : 0);
    double h_z = (ncells_z > 0 ? p.z()/ncells_z : 0);
    for(unsigned int i = 0;i < nnodes_y;i++)
      ycoords[i] = i*h_y + box.P1().y();
    for(unsigned int i = 0;i < nnodes_x;i++)
      xcoords[i] = i*h_x + box.P1().x();
    unsigned int point_id = 1;
    for(unsigned int k = 0;k < nnodes_z;k++){
      double z = k*h_z + box.P1().z();
      for(unsigned int j = 0;j < nnodes_y;j++){
	for(unsigned int i = 0;i < nnodes_x;i++){
	  GeoPrim::C3Point gridpoint(nc[point_id++]);
	  gridpoint.x() = xcoords[i];
	  gridpoint.y() = ycoords[j];
	  gridpoint.z() = z;
	}
      }
    }
  }
}
