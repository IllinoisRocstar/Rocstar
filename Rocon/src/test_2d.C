#include "Mesh.H"

int main(int argc,char *argv[])
{
  Mesh::NodalCoordinates nc(8);
  Mesh::Connectivity     con(6);
  Mesh::Connectivity     dual_con;
  std::cout << "Setting nodes." << std::endl;
  nc.init_node(1,GeoPrim::CPoint(0,0,0));
  nc.init_node(2,GeoPrim::CPoint(1,0,0));
  nc.init_node(3,GeoPrim::CPoint(1,1,0));
  nc.init_node(4,GeoPrim::CPoint(0,1,0));
  nc.init_node(5,GeoPrim::CPoint(0,0,1));
  nc.init_node(6,GeoPrim::CPoint(1,0,1));
  nc.init_node(7,GeoPrim::CPoint(1,1,1));
  nc.init_node(8,GeoPrim::CPoint(0,1,1));
  std::cout << "Building Mesh." << std::endl;
  con[0].push_back(1);  con[0].push_back(2);
  con[0].push_back(3);  con[0].push_back(4);
  con[1].push_back(5);  con[1].push_back(6);
  con[1].push_back(7);  con[1].push_back(8);
  con[2].push_back(1);  con[2].push_back(2);
  con[2].push_back(6);  con[2].push_back(5);
  con[3].push_back(2);  con[3].push_back(3);
  con[3].push_back(7);  con[3].push_back(6);
  con[4].push_back(3);  con[4].push_back(4);
  con[4].push_back(8);  con[4].push_back(7);
  con[5].push_back(4);  con[5].push_back(1);
  con[5].push_back(5);  con[5].push_back(8);
  con.Sync();
  std::cout << "Getting Inverse." << std::endl;
  con.Inverse(dual_con,8);
  std::cout << "Creating points." << std::endl;
  std::vector<GeoPrim::CPoint> points;
  GeoPrim::CPoint point1(.5,.5,0);
  GeoPrim::CPoint point2(.5,0,0);
  GeoPrim::CPoint point3(.5,1,0);
  GeoPrim::CPoint point4(0,.5,0);
  GeoPrim::CPoint point5(1,.5,0);
  GeoPrim::CPoint point6(.5,.5,1);
  GeoPrim::CPoint point7(.5,.5,2);
  GeoPrim::CPoint point8(2,.5,1);
  points.push_back(point1); points.push_back(point2);
  points.push_back(point3); points.push_back(point4);
  points.push_back(point5); points.push_back(point6);
  points.push_back(point7); points.push_back(point8);
  std::cout << "Building box." << std::endl;
  GeoPrim::CBox box(GeoPrim::CPoint(0,0,0),GeoPrim::CPoint(1,1,1));
  GeoPrim::CVector natc;
  std::cout << "Searching for points" << std::endl;
  std::vector<GeoPrim::CPoint>::iterator pi = points.begin();
  while(pi != points.end()){
    Mesh::IndexType point_number = pi - points.begin() + 1;
    Mesh::IndexType element_id = FindPointInMesh_2(*pi++,nc,con,dual_con,box,natc);
    if(!element_id)
      std::cout << "Could not find point " << point_number 
		<< " in the mesh" << std::endl;
    else {
      std::cout << "Found point " << point_number << " in element " 
		<< element_id << " with natc = " << natc << std::endl;
      GeoPrim::CVector nodal_coords[4];
      nodal_coords[0].init(nc[con[element_id-1][0]]);
      nodal_coords[1].init(nc[con[element_id-1][1]]);
      nodal_coords[2].init(nc[con[element_id-1][2]]);
      nodal_coords[3].init(nc[con[element_id-1][3]]);
      Mesh::GenericCell_2 cell(4);
      GeoPrim::CVector point_coords;
      cell.interpolate(nodal_coords,natc,point_coords);
      std::cout << "The interpolated coordinates are: " << point_coords << std::endl;
    }
  }
  return(0);
}
