/// 
/// VTK utilities for the Mesh data structures
/// Author: George Zagaris (gzagaris@illinois.edu)
///
#ifndef _MESH_VTK_H_
#define _MESH_VTK_H_

#include "Mesh.H"

namespace Mesh {

/**
 * @brief Writes the merged grid and its data in VTK UNSTRUCTURED_GRID file format.
 * @pre Program.globalNodeList.size( ) 	  > 0.
 * @pre Program.globalElementList.size( ) > 0.
 * @pre Program.smdv.size( ) > 0.
 * @pre Program.globalSolution.size( ) == Program.smdv.size( ).
 * @pre Program.solutionignore.size( ) == Program.smdv.size( ).
 * @pre Program.outputfile != ""
 * @post A VTK formatted file is written in Program.outputfile.vtk
 * @note The extension in the outputfile will be appended.
 */
  void writeVtkData( );

/**
 * @brief Writes the set of meshes into separate VTK files that can be
 * visualized with Paraview. The format of the files is {outputFile}.paneId.vtk.
 * @note This method is mainly used for debugging.
 * @param m the set of meshes.
 */
void writeVtkFiles( std::vector< Mesh::UnstructuredMesh > &m );

/**
 * @brief This method prints the global grid in VTK file format.
 * @note Mostly used for debugging.
 * @param vlist the global vertex list.
 * @param elist the element list.
 */
void printVtk( std::vector< double > &vlist, std::vector< std::vector< int > > &elist );

}

#endif
