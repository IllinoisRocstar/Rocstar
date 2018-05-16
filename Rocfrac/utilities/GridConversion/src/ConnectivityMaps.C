#include <iomanip>
#include "DriverProgram.H"

namespace GridConversion{ namespace DriverProgram{

  int SerialProgram::ConnectivityMaps(SolverUtils::Mesh::Connectivity nodesToElems){

      std::ostringstream Ostr;
      std::stringstream ss;

      FunctionEntry("ConnectivityMaps");

      std::cout << "In ConnectivityMaps function" << std::endl;

      //resize nodeToNode map
      nodeToNode.resize(numNodes);

      int elem, node;

      //Generate node to node map
      //loop over every node
      for(int i=0; i < numNodes; i++){
        std::cout << "node " << i+1 << ":" << std::endl;
        //loop over every element for that node
        for(int j=0; j < nodesToElems[i].size(); j++){
          std::cout << "element " << nodesToElems[i][j] << ":" << std::endl;
          elem = nodesToElems[i][j]-1;
          //loop over every node for that element
          for(int k=0; k < numNodesPerElem; k++){
            node = elems[numNodesPerElem*elem + k];
            //don't add the node to its own list
            if(node == i+1)
              continue;
            //loop over all the entries for the nodeToNode map and
            //see if this node has been added yet
            bool newValue=true;
            for(int l=0; l < nodeToNode[i].size(); l++){
              if(node == nodeToNode[i][l]){
                newValue=false;
                break;
              }
            }
            //if its a new value add it
            if(newValue){
              nodeToNode[i].push_back(node);
            }
          }
        }//element around node loop
      }//node lope

      //print nodeToNode to check
      if(verblevel > 3){
        Ostr << "Node to node map:" << std::endl;
        for(int i=0; i < nodeToNode.size(); i++){
          Ostr << i+1 << ":" << std::endl;
          for(int j=0; j < nodeToNode[i].size(); j++){
            Ostr << nodeToNode[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }

      //Generate list of element edges from node to node map
      //loop over every node
      numNodesPerElemEdge = 2;
      numElemEdges = 0;
      for(int i=0; i < numNodes; i++){
        //loop over every node connected to that node
        for(int j=0; j < nodeToNode[i].size(); j++){
          node = nodeToNode[i][j];
          //if the node # is greater than the current node then it
          //hasn't been processed yet so its a new edge
          if(node > i+1){
            elemEdges.push_back(i+1);
            elemEdges.push_back(node);
            numElemEdges++;
          }
        }//nodes around node loop
      }//node loop      

      //print nodeToNode to check
      if(verblevel > 3){
        Ostr << "Element edges:" << std::endl;
        for(int i=0; i < numElemEdges; i++){
          Ostr << i+1 << ":" << std::endl;
          for(int j=0; j < numNodesPerElemEdge; j++){
            Ostr << elemEdges[i*numNodesPerElemEdge + j] << " ";
          }
          Ostr << std::endl;
        }
      }

      //Genereate a map from elements to element edges
      numEdgesPerElem = 6; //Everything is pretty much hard coded
                           //for linear tets right now
      //loop over every element
      for(int i=0; i < numElems; i++){
        //loop over every node for the element
        for(int j=0; j < numNodesPerElem; j++){
          int n1, n2;
          n1 = elems[i*numNodesPerElem + j];
          //loop over every other node for the element
          for(int k=j+1; k < numNodesPerElem; k++){   
            n2 = elems[i*numNodesPerElem + k];
            //loop over every element edge to find what edge 
            //these two nodes are
            for(int l=0; l < numElemEdges; l++){
              if(n1 == elemEdges[l*numNodesPerElemEdge] &&
                 n2 == elemEdges[l*numNodesPerElemEdge + 1])
                elemToElemEdges.push_back(l+1);         
              else if(n2 == elemEdges[l*numNodesPerElemEdge] &&
                 n1 == elemEdges[l*numNodesPerElemEdge + 1])
                elemToElemEdges.push_back(l+1);         
            }//loop over every element edge
          }//loop over every other node for the element
        }//loop over every node for the element
      }//loop over every element

      //print elemToElemEdges to check
      if(verblevel > 3){
        Ostr << "Element to element edges:" << std::endl;
        for(int i=0; i < numElems; i++){
          Ostr << i+1 << ":" << std::endl;
          for(int j=0; j < numEdgesPerElem; j++){
            Ostr << elemToElemEdges[i*numEdgesPerElem + j] << " ";
          }
          Ostr << std::endl;
        }
      }

      StdOut(Ostr.str());
      Ostr.str("");

      FunctionExit("ConnectivityMaps");

    return 0;
 
  } //HigherOrderTets function

}; //DriverProgram namespace
}; //GridConversion namespace
