#include <iomanip>
#include "DriverProgram.H"

namespace GridConversion{ namespace DriverProgram{

  int SerialProgram::HigherOrderTets(SolverUtils::Mesh::Connectivity elemEdgeToElems,
                                     SolverUtils::Mesh::Connectivity nodesToDomains){

      std::ostringstream Ostr;
      std::stringstream ss;

      FunctionEntry("HigherOrderTets");

      //new numNodesPerElem
      int newNumNPE=10;

      //a vector for the new elements
      std::vector<unsigned int> newElems(numElems*newNumNPE,0);

      //Loop over every element
      for(int i=0; i < numElems; i++){
       std::cout << "elem " << i+1 << std::endl;
       //loop over all the nodes for the element
       for(int j=0; j < numNodesPerElem; j++){
         //first populate the new element with the original nodes
         newElems[i*newNumNPE + j] = elems[i*numNodesPerElem + j];
       }//loop over nodes for the element
      }//loop over all the elements

      //Loop over every edge (we want to add a node in the middle
      //of each edge
      double x, y, z;
      unsigned int n1, n2;
      for(int i=0; i < numElemEdges; i++){

        std::cout << "edge " << i+1 << std::endl;
        //create the new node
        numNodes++;
        n1 = elemEdges[i*numNodesPerElemEdge]-1;
        n2 = elemEdges[i*numNodesPerElemEdge + 1]-1;

        x = (nodes[3*n1 + 0] + nodes[3*n2 + 0])/2.0;
        y = (nodes[3*n1 + 1] + nodes[3*n2 + 1])/2.0;
        z = (nodes[3*n1 + 2] + nodes[3*n2 + 2])/2.0;

        nodes.push_back(x);
        nodes.push_back(y);
        nodes.push_back(z);

        std::cout << " new node " << numNodes << ": "
                  << x << " " << y << " " << z << std::endl;

        //add the new node to the appropriate domain
        //if its on one
        //loop over the domains for the first node
        for(int j=0; j < nodesToDomains[n1].size(); j++){
          //loop over the domains for the second node
          for(int k=0; k < nodesToDomains[n2].size(); k++){
            //If they share a domain this new node is also on it
            if(nodesToDomains[n1][j] == nodesToDomains[n2][k]){
              int domain = nodesToDomains[n1][j]-1;
              domains[domain].push_back(numNodes);
            }
          }
        } 

        //now add the node to the appropriate places
        //in the element connectivity

        //loop over every element for the edge
        for(int j=0; j < elemEdgeToElems[i].size(); j++){
          int elem = elemEdgeToElems[i][j];
          std::cout << "  element " << elem << std::endl;
          //find which edge this edge is for the element
          //(1st, 2nd, ... 6th edge)
          //loop over the edges for the element
          for(int k=0; k < numEdgesPerElem; k++){
            std::cout << "    edge " << elemToElemEdges[numEdgesPerElem*(elem-1) + k] << std::endl;
            //if we found the edge number add the node
            //to the element in the appropriate place
            //(if it hasn't already been added)
            if(elemToElemEdges[numEdgesPerElem*(elem-1) + k] == i+1){
              switch (k) {
                case 0:
                  //Add the node if it hasn't already been added
                  if(newElems[newNumNPE*(elem-1) + 4] == 0){
                    newElems[newNumNPE*(elem-1) + 4] = numNodes;
                  }
                  break;
                case 1:
                  //Add the node if it hasn't already been added
                  if(newElems[newNumNPE*(elem-1) + 6] == 0){
                    newElems[newNumNPE*(elem-1) + 6] = numNodes;
                  }
                  break;
                case 2:
                  //Add the node if it hasn't already been added
                  if(newElems[newNumNPE*(elem-1) + 7] == 0){
                    newElems[newNumNPE*(elem-1) + 7] = numNodes;
                  }
                  break;
                case 3:
                  //Add the node if it hasn't already been added
                  if(newElems[newNumNPE*(elem-1) + 5] == 0){
                    newElems[newNumNPE*(elem-1) + 5] = numNodes;
                  }
                  break;
                case 4:
                  //Add the node if it hasn't already been added
                  if(newElems[newNumNPE*(elem-1) + 8] == 0){
                    newElems[newNumNPE*(elem-1) + 8] = numNodes;
                  }
                  break;
                case 5:
                  //Add the node if it hasn't already been added
                  if(newElems[newNumNPE*(elem-1) + 9] == 0){
                    newElems[newNumNPE*(elem-1) + 9] = numNodes;
                  }
                  break;
              }//switch over k
            }//if we found the edge
          }//loop over edges for the element
        }//loop over elements for the edge
      }//loop over every edge

      std::cout << "line " << __LINE__ << std::endl;

      //print to check
      if(verblevel > 3){
        Ostr << "new elements: " << std::endl;
        for(int i=0; i < numElems; i++){
          Ostr << i+1 << ":" << std::endl;
          for(int j=0; j < newNumNPE; j++){
            Ostr << newElems[i*newNumNPE + j] << " ";
          }
          Ostr << std::endl;
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");
 
      std::cout << "line " << __LINE__ << std::endl;
      //save the new elements arrray to the old one
      elems.resize(numElems*newNumNPE);
      elems = newElems;
      numNodesPerElem = newNumNPE;

      StdOut(Ostr.str());
      Ostr.str("");

      std::cout << "line " << __LINE__ << std::endl;
      FunctionExit("HigherOrderTets");

      std::cout << "line " << __LINE__ << std::endl;
    return 0;
 
  } //HigherOrderTets function

}; //DriverProgram namespace
}; //GridConversion namespace
