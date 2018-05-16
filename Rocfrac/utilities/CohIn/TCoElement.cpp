

#include "TCoElement.hpp"
#include "TriFace.hpp"
#include "Mesh.hpp"

TCoElement::TCoElement() :
Element( e_tri_cohesive ){

  d_nodes = new Node*[getNumNodes()];
  d_faces = new Face*[getNumFaces()];

}


TCoElement::TCoElement(Node ** thenodes) :
Element( e_tri_cohesive ){

  d_nodes = new Node*[getNumNodes()];
  d_faces = new Face*[getNumFaces()];
  setFromNodes( thenodes );

}

TCoElement::~TCoElement() {

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

int TCoElement::getNumNodes() const
{ return 6; }

int TCoElement::getNumFaces() const
{ return 2; }


void TCoElement::replaceFaceNode (Node *node, Node *new_node, Face* face ){

  int start = ( face == d_faces[0] ? 0 : 3 );
  int i;
  for( i = start; i < start+3; i++ ){
    if( d_nodes[i] == node ){
      d_nodes[i] = new_node;
      node->removeElement( this );
      new_node->addElement( this );
    }
  }
  face->replaceNode( node, new_node );
}

void TCoElement::setFaceFromNodes(int num, Node** nodes){

  // find if face exist (if not creates it)
  Face* face = nodes[0]->sharedFace( nodes[1], nodes[2] );
  if( !face ){
    face = new TriFace( nodes );
    s_mesh->addFace( face );
  }
  setFace( num, face );
}


void TCoElement::setFromMyNodes(){

  // !!! check ordering

  Node *nodes[3];
  nodes[0] = d_nodes[0];
  nodes[1] = d_nodes[1];
  nodes[2] = d_nodes[2];
  setFaceFromNodes( 0, nodes );

  nodes[0] = d_nodes[5];
  nodes[1] = d_nodes[4];
  nodes[2] = d_nodes[3];
  setFaceFromNodes( 1, nodes );

  // here separate node loops
  int i;
  for( i = 0; i < 3; i++ ){
    if( d_nodes[i] != d_nodes[ 3 + i] ){
      d_nodes[i]->addNextLink( d_nodes[ 3 + i ] );
    }
  }
}

