

#include "HexElement.hpp"
#include "QuadFace.hpp"
#include "Mesh.hpp"

HexElement::HexElement() :
Element( e_hex ){

  d_nodes = new Node*[getNumNodes()];
  d_faces = new Face*[getNumFaces()];

}


HexElement::HexElement(Node ** thenodes) :
Element( e_hex ){

  d_nodes = new Node*[getNumNodes()];
  d_faces = new Face*[getNumFaces()];
  setFromNodes( thenodes );

}

HexElement::~HexElement() {

  int numn = getNumNodes();
  int i;
  for( i = 0; i < numn; i++ ){
    d_nodes[i]->removeElement( this );
  }
  int numf = getNumFaces();
  for( i = 0; i < numf; i++ ){
    d_faces[i]->removeElement( this );
  }
}

int HexElement::getNumNodes() const
{ return 8; }

int HexElement::getNumFaces() const
{ return 6; }

void HexElement::setFaceFromNodes(int num, Node** nodes){

  // find if face exist (if not creates it)
  Face* face = nodes[0]->sharedFace( nodes[1], d_nodes[2], d_nodes[3] );
  if( !face ){
    face = new QuadFace( nodes );
    s_mesh->addFace( face );
  }
  setFace( num, face );
}


void HexElement::setFromMyNodes(){

  // !!! check ordering

  Node *nodes[4];
  nodes[0] = d_nodes[7];
  nodes[1] = d_nodes[6];
  nodes[2] = d_nodes[5];
  nodes[3] = d_nodes[4];
  setFaceFromNodes( 0, nodes );

  nodes[0] = d_nodes[0];
  nodes[1] = d_nodes[1];
  nodes[2] = d_nodes[2];
  nodes[3] = d_nodes[3];
  setFaceFromNodes( 1, nodes );

  nodes[0] = d_nodes[3];
  nodes[1] = d_nodes[7];
  nodes[2] = d_nodes[4];
  nodes[3] = d_nodes[0];
  setFaceFromNodes( 2, nodes );

  nodes[0] = d_nodes[1];
  nodes[1] = d_nodes[5];
  nodes[2] = d_nodes[6];
  nodes[3] = d_nodes[2];
  setFaceFromNodes( 3, nodes );

  nodes[0] = d_nodes[0];
  nodes[1] = d_nodes[4];
  nodes[2] = d_nodes[5];
  nodes[3] = d_nodes[1];
  setFaceFromNodes( 4, nodes );

  nodes[0] = d_nodes[2];
  nodes[1] = d_nodes[6];
  nodes[2] = d_nodes[7];
  nodes[3] = d_nodes[3];
  setFaceFromNodes( 5, nodes );
}
