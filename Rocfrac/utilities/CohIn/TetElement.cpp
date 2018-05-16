

#include "TetElement.hpp"
#include "TriFace.hpp"
#include "Mesh.hpp"

TetElement::TetElement() :
Element( e_tet ){

  d_nodes = new Node*[getNumNodes()];
  d_faces = new Face*[getNumFaces()];

}


TetElement::TetElement(Node ** thenodes) :
Element( e_tet ){

  d_nodes = new Node*[getNumNodes()];
  d_faces = new Face*[getNumFaces()];
  setFromNodes( thenodes );

}

TetElement::~TetElement() {

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

int TetElement::getNumNodes() const
{ return 4; }

int TetElement::getNumFaces() const
{ return 4; }


void TetElement::setFaceFromNodes(int num, Node** nodes){

  // find if face exist (if not creates it)
  Face* face = nodes[0]->sharedFace( nodes[1], nodes[2] );
  if( !face ){
    face = new TriFace( nodes );
    s_mesh->addFace( face );
  }
  setFace( num, face );
}


void TetElement::setFromMyNodes(){

  // !!! check ordering

  Node *nodes[3];
  nodes[0] = d_nodes[1];
  nodes[1] = d_nodes[3];
  nodes[2] = d_nodes[2];
  setFaceFromNodes( 0, nodes );

  nodes[0] = d_nodes[2];
  nodes[1] = d_nodes[3];
  nodes[2] = d_nodes[0];
  setFaceFromNodes( 1, nodes );

  nodes[0] = d_nodes[0];
  nodes[1] = d_nodes[3];
  nodes[2] = d_nodes[1];
  setFaceFromNodes( 2, nodes );

  nodes[0] = d_nodes[0];
  nodes[1] = d_nodes[1];
  nodes[2] = d_nodes[2];
  setFaceFromNodes( 3, nodes );

}

double TetElement::getMinEdgeLength(){

   double minEdge = 1E+25;
   double currentLength;
   MVec pos[4];
   
   //there are six edge lengths to check:
   //0 -> 1, 0 -> 2, 0 -> 3
   //1 -> 2, 1 -> 3, 2 -> 3
   // for nodes numbered 0, 1, 2, 3
   
   for (int i=0; i<4; i++) {
   	pos[i] = d_nodes[i]->getPosition();
   }
   
   currentLength = pos[0].distance_between(pos[1]);
   if (currentLength < minEdge) minEdge = currentLength;
   currentLength = pos[0].distance_between(pos[2]);
   if (currentLength < minEdge) minEdge = currentLength;
   currentLength = pos[0].distance_between(pos[3]);
   if (currentLength < minEdge) minEdge = currentLength;
   currentLength = pos[1].distance_between(pos[2]);
   if (currentLength < minEdge) minEdge = currentLength;
   currentLength = pos[1].distance_between(pos[2]);
   if (currentLength < minEdge) minEdge = currentLength;
   currentLength = pos[2].distance_between(pos[3]);
   if (currentLength < minEdge) minEdge = currentLength;
   
   assert(minEdge < 1E+24);
   
   return(minEdge); 

}

