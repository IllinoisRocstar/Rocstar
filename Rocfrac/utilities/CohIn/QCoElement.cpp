

#include "QCoElement.hpp"
#include "QuadFace.hpp"
#include "Mesh.hpp"

QCoElement::QCoElement() :
Element( e_quad_cohesive ){

  d_nodes = new Node*[getNumNodes()];
  d_faces = new Face*[getNumFaces()];

}


QCoElement::QCoElement(Node ** thenodes) :
Element( e_quad_cohesive ){

  d_nodes = new Node*[getNumNodes()];
  d_faces = new Face*[getNumFaces()];
  setFromNodes( thenodes );
}

QCoElement::~QCoElement() {

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

int QCoElement::getNumNodes() const
{ return 8; }

int QCoElement::getNumFaces() const
{ return 2; }


void QCoElement::replaceFaceNode (Node *node, Node *new_node, Face* face ){

  int start = ( face == d_faces[0] ? 0 : 4 );
  int i;
  for( i = start; i < start+4; i++ ){
    if( d_nodes[i] == node ){
      d_nodes[i] = new_node;
      node->removeElement( this );
      new_node->addElement( this );
    }
  }
  face->replaceNode( node, new_node );
}

void QCoElement::setFaceFromNodes(int num, Node** nodes){

  // find if face exist (if not creates it)
  Face* face = nodes[0]->sharedFace( nodes[1], d_nodes[2], nodes[3] );
  if( !face ){
    face = new QuadFace( nodes );
    s_mesh->addFace( face );
  }
  setFace( num, face );
}


void QCoElement::setFromMyNodes(){

  // !!! check ordering

  Node *nodes[4];
  nodes[0] = d_nodes[0];
  nodes[1] = d_nodes[1];
  nodes[2] = d_nodes[2];
  nodes[3] = d_nodes[3];
  setFaceFromNodes( 0, nodes );

  nodes[0] = d_nodes[7];
  nodes[1] = d_nodes[6];
  nodes[2] = d_nodes[5];
  nodes[3] = d_nodes[4];
  setFaceFromNodes( 1, nodes );

  // here separate node loops
  int i;
  for( i = 0; i < 4; i++ ){
    if( d_nodes[i] != d_nodes[ 4 + i] ){
      d_nodes[i]->addNextLink( d_nodes[ 4 + i ] );
    }
  }

}

