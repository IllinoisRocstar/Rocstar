/// 
/// @file
/// @ingroup gridconversion_group
/// @brief Driver program for Grid Conversion.
/// @author Jessica Kress (jkress@illinoisrocstar.com)
/// @date 11/3/2015
/// 

#include "DriverProgram.H"
#include "Mesh.H"

namespace GridConversion {

  namespace DriverProgram {

    int SerialProgram::Run()
    {
      // FunctionEntry("NAME"): Updates the user-defined stack and 
      // the program profiles with timing information. The placement
      // of FunctionEntry and FunctionExit calls is at the developer's
      // discretion.  
      FunctionEntry("Run");

      // ---------- The Program -----------------
      // Drives the grid conversion process

      std::ostringstream Ostr;
      int error=0;        

      // Parse the input file
      error = ReadStanfordInput();
      if(error == 1){
        std::ostringstream ErrOstr;
        ErrOstr << "Error reading input file." << std::endl;
        StdOut(Ostr.str());
        Ostr.str("");
        ErrOut(ErrOstr.str());
        ErrOstr.str("");
        exit(1);
      }
  
      //Print to check nodes
      if(verblevel > 3){
        Ostr << "nodes: " << std::endl;
        for(int i=0; i < numNodes; i++){
          Ostr << i+1 << ":" << std::endl;
          for(int j=0; j < 3; j++){
            Ostr << nodes[i*3 + j] << " ";
          }
          Ostr << std::endl;
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");

      //Print to check elems
      if(verblevel > 3){
        Ostr << "elements: " << std::endl;
        for(int i=0; i < numElems; i++){
          Ostr << i+1 << ":" << std::endl;
          for(int j=0; j < numNodesPerElem; j++){
            Ostr << elems[i*numNodesPerElem + j] << " ";
          }
          Ostr << std::endl;
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");

      //Print to check domains
      if(verblevel > 3){
        Ostr << "domains: " << std::endl;
        for(int i=0; i < domains.size(); i++){
          Ostr << i+1 << ":" << std::endl;
          for(int j=0; j < domains[i].size(); j++){
            Ostr << domains[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }

     StdOut(Ostr.str());
     Ostr.str("");
     std::cout << "line " << __LINE__ << std::endl;

      // Populate element connectivty with elems vector
      SolverUtils::Mesh::Connectivity utilsElems;
      utilsElems.AddElements(numElems, numNodesPerElem, elems);
     std::cout << "line " << __LINE__ << std::endl;

      // Get a map of nodes to elements
      SolverUtils::Mesh::Connectivity nodesToElems;
      utilsElems.Sync();
      utilsElems.Inverse(nodesToElems,numNodes);
     std::cout << "line " << __LINE__ << std::endl;

      // Populate domain connectivity with domain vectors
      SolverUtils::Mesh::Connectivity utilsDomains;
      for(int i=0; i < domains.size(); i++){
        utilsDomains.AddElement(domains[i]);
      }
     std::cout << "line " << __LINE__ << std::endl;


     StdOut(Ostr.str());
     Ostr.str("");
     std::cout << "line " << __LINE__ << std::endl;

      //Print to check domain connectivity
      if(verblevel > 3){
        Ostr << "Domain connectivity: " << std::endl;
        for(int i=0; i < utilsDomains.size(); i++){
          Ostr << i << ":" << std::endl;
          for(int j=0; j < utilsDomains[i].size(); j++){
            Ostr << utilsDomains[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");

      // Get a map of nodes to domains
      utilsDomains.Sync();
      utilsDomains.Inverse(nodesToDomains,numNodes);

      // Print to check the nodes to elements
      if(verblevel > 3){
        Ostr << "Nodes to elements: " << std::endl;
        for(int i=0; i < nodesToElems.size(); i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < nodesToElems[i].size(); j++){
            Ostr << nodesToElems[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }

      // Print to check the nodes to domains
      if(verblevel > 3){
        Ostr << "Nodes to domains: " << std::endl;
        for(int i=0; i < nodesToDomains.size(); i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < nodesToDomains[i].size(); j++){
            Ostr << nodesToDomains[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");

      //If we have tet elements and the user specifies higher
      //order then we need to make the tets quadratic     
      //NOTE: This functionality is not complete!!!!!!!!!!! 
      if(quadratic){
        if(numNodesPerElem != 4){
          std::ostringstream ErrOstr;
          ErrOstr << "Only 4 node tetrahedra can be made into" << std::endl
                  << "quadratic elements!" << std::endl;
          StdOut(Ostr.str());
          Ostr.str("");
          ErrOut(ErrOstr.str());
          ErrOstr.str("");
          exit(1);
        }

        //Create node to node map
        ConnectivityMaps(nodesToElems);
  
        // Populate element to element edge for Solver Utils
        SolverUtils::Mesh::Connectivity utilsElemToElemEdges;
        utilsElemToElemEdges.AddElements(numElems, numEdgesPerElem, elemToElemEdges);
        std::cout << "line " << __LINE__ << std::endl;

        // Get a map of nodes to elements
        SolverUtils::Mesh::Connectivity elemEdgeToElems;
        utilsElemToElemEdges.Sync();
        utilsElemToElemEdges.Inverse(elemEdgeToElems,numElemEdges);

        //Print to check elemEdgesToElems
        if(verblevel > 3){
          Ostr << "Element edges to elements: " << std::endl;
          for(int i=0; i < elemEdgeToElems.size(); i++){
            Ostr << i+1 << ":" << std::endl;
            for(int j=0; j < elemEdgeToElems[i].size(); j++){
              Ostr << elemEdgeToElems[i][j] << " ";
            }
            Ostr << std::endl;
          }
          std::cout << "line " << __LINE__ << std::endl;
        }
 
        //Create higher order tets
        HigherOrderTets(elemEdgeToElems,nodesToDomains);

        //Print to check nodes
        if(verblevel > 3){
          Ostr << "(after higher order) nodes: " << std::endl;
          for(int i=0; i < numNodes; i++){
            Ostr << i+1 << ":" << std::endl;
            for(int j=0; j < 3; j++){
              Ostr << nodes[i*3 + j] << " ";
            }
            Ostr << std::endl;
          }
        }
        StdOut(Ostr.str());
        Ostr.str("");

        //Print to check elems
        if(verblevel > 3){
          Ostr << "(after higher order) elements: " << std::endl;
          for(int i=0; i < numElems; i++){
            Ostr << i+1 << ":" << std::endl;
            for(int j=0; j < numNodesPerElem; j++){
              Ostr << elems[i*numNodesPerElem + j] << " ";
            }
            Ostr << std::endl;
          }
          std::cout << "line " << __LINE__ << std::endl;
        }
        StdOut(Ostr.str());
        Ostr.str("");
        std::cout << "line " << __LINE__ << std::endl;
       
        //Print to check domains
        if(verblevel > 3){
          Ostr << "(after higher order) domains: " << std::endl;
          for(int i=0; i < domains.size(); i++){
            Ostr << i+1 << ":" << std::endl;
            for(int j=0; j < domains[i].size(); j++){
              Ostr << domains[i][j] << " ";
            }
            Ostr << std::endl;
          }
        }

      // Populate element connectivty with elems vector
      SolverUtils::Mesh::Connectivity newUtilsElems;
      newUtilsElems.AddElements(numElems, numNodesPerElem, elems);
     std::cout << "line " << __LINE__ << std::endl;

      // Get a map of nodes to elements
      nodesToElems.clear();
      newUtilsElems.Sync();
      newUtilsElems.Inverse(nodesToElems,numNodes);
     std::cout << "line " << __LINE__ << std::endl;

      // Populate domain connectivity with domain vectors
      SolverUtils::Mesh::Connectivity newUtilsDomains;
      for(int i=0; i < domains.size(); i++){
        newUtilsDomains.AddElement(domains[i]);
      }
     std::cout << "line " << __LINE__ << std::endl;


     StdOut(Ostr.str());
     Ostr.str("");
     std::cout << "line " << __LINE__ << std::endl;

      // Get a map of nodes to domains
      nodesToDomains.clear();
      newUtilsDomains.Sync();
      newUtilsDomains.Inverse(nodesToDomains,numNodes);

      // Print to check the nodes to elements
      if(verblevel > 3){
        Ostr << "Nodes to elements: " << std::endl;
        for(int i=0; i < nodesToElems.size(); i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < nodesToElems[i].size(); j++){
            Ostr << nodesToElems[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }

      // Print to check the nodes to domains
      if(verblevel > 3){
        Ostr << "Nodes to domains: " << std::endl;
        for(int i=0; i < nodesToDomains.size(); i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < nodesToDomains[i].size(); j++){
            Ostr << nodesToDomains[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");


      }//creating higher order tets 

      //Populate elems to domains vector using nodes to domains
      //loop over all the elements
      elemsToDomains.resize(numElems);
      

      for(int i=0; i < numElems; i++){
        std::map<int, int>domCount;
        std::cout << "elem " << i+1 << std::endl;
        //loop over the nodes for the element
        for(int j=0; j < numNodesPerElem; j++){
          int node = elems[numNodesPerElem*i + j];
          std::cout << "node " << j+1 << " (" << node << ")" << std::endl;
          node--;
          //loop over the domains for that node
          for(int k=0; k < nodesToDomains[node].size(); k++){
            std::cout << "domain " << k+1 << " (" 
                      << nodesToDomains[node][k] << ")" << std::endl;
            int domain = nodesToDomains[node][k];
            //We only want to add the domain if 3 or more of the element's
            //nodes reside on it (i.e., an entire face of the element)
            if( domCount.find(domain) == domCount.end() )
              domCount[domain] = 1;
            else
              domCount[domain]++;

            std::cout << "domCount = " << domCount[domain] << std::endl;
 
            bool newValue=true;
            //See if the domain has been added yet
            for(int l=0; l < elemsToDomains[i].size(); l++){
              if(elemsToDomains[i][l] == nodesToDomains[node][k]){
                newValue=false;
                break;
              }
            }//loop over domains for elem
            //add the domain to the elem
            if(newValue && domCount[domain] > 2){
              std::cout << "Adding domain" << std::endl;
              elemsToDomains[i].push_back(nodesToDomains[node][k]);
            }
          }//loop over domains for node
        }// loop over nodes for elem
      }//loop over elements

     // Print to check elems to domains
      if(verblevel > 3){
        Ostr << "Elems to domains: " << std::endl;
        for(int i=0; i < elemsToDomains.size(); i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < elemsToDomains[i].size(); j++){
            Ostr << elemsToDomains[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }

      // Print the domain bcs to check
      if(verblevel > 3){
        Ostr << "domain bcs: " << std::endl;
        for(int i=0; i < domainBCs.size(); i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < domainBCs[i].size(); j++){
            Ostr << domainBCs[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }

      // Print the domain bc values to check
      if(verblevel > 3){
        Ostr << "domain bc values: " << std::endl;
        for(int i=0; i < domainBCValues.size(); i++){
          Ostr << i+1 << ": " << std::endl;
          for(int j=0; j < domainBCValues[i].size(); j++){
            for(int k=0; k < domainBCValues[i][j].size(); k++){
              Ostr << "    " <<  domainBCValues[i][j][k] << " ";
            }
            Ostr << std::endl;
          }
        }
      }


      // Print the edges to check
      if(verblevel > 3){
        Ostr << "edges: " << std::endl;
        for(int i=0; i < edges.size(); i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < edges[i].size(); j++){
            Ostr << edges[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }

      // Print the edge bcs to check
      if(verblevel > 3){
        Ostr << "edge bcs: " << std::endl;
        for(int i=0; i < edgeBCs.size(); i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < edgeBCs[i].size(); j++){
            Ostr << edgeBCs[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }

      // Print the edge bc values to check
      if(verblevel > 3){
        Ostr << "edge bc values: " << std::endl;
        for(int i=0; i < edgeBCValues.size(); i++){
          Ostr << i+1 << ": " << std::endl;
          for(int j=0; j < edgeBCValues[i].size(); j++){
            for(int k=0; k < edgeBCValues[i][j].size(); k++){
              Ostr << "    " <<  edgeBCValues[i][j][k] << " ";
            }
            Ostr << std::endl;
          }
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");

      std::cout << "line " << __LINE__ << std::endl; 
      // Populate the nodal bc vectors from the domain bc vectors
      nodeBCs.resize(numNodes);
      nodeBCValues.resize(numNodes);
      //loop over all domains
      for(int i=0; i < domains.size(); i++){
        // loop over all nodes in the domain
        for(int j=0; j < domains[i].size(); j++){
          int node;
          node = domains[i][j]-1;
          // they are on
          // loop over all the bc flags for this domain
          for(int k=0; k < domainBCs[i].size(); k++){
            bool newValue = true;
            int type8BC=-1;
            // loop over all the bc flags for this node
            for(int l=0; l < nodeBCs[node].size(); l++){
              // if the node already has this bc & value skip it
              if(nodeBCs[node][l] == domainBCs[i][k]){
                type8BC=l;
                for(int m=0; m < domainBCValues[i][k].size(); m++){
                  if(domainBCValues[i][k][m] != nodeBCValues[node][l][m])
                    break;
                  else if(m == domainBCValues[i][k].size()-1)
                    newValue = false;
                }
              }
            }// node bc flags
            // if it is a new flag add it and its values to 
            // the node's vectors
            if(newValue){
              //nodes can have as many type 6 bcs as needed
              if(domainBCs[i][k] == 6){
                nodeBCs[node].push_back(domainBCs[i][k]);
                nodeBCValues[node].push_back(domainBCValues[i][k]); 
              }
              //nodes can only have one type 8 bc which is determined
              //by a "priority" number that is the last bc value number
              if(domainBCs[i][k] == 8){
                //the node doesn't have any type 8 bcs yet
                int last = domainBCValues[i][k].size()-1;
                if(type8BC == -1){
                  nodeBCs[node].push_back(domainBCs[i][k]);
                  nodeBCValues[node].push_back(domainBCValues[i][k]); 
                }
                //it does have a type 8 bc already & we should keep the one
                //with the highest priority (here a higher priority means a
                //lower number in the last value position)
                else if(nodeBCValues[node][type8BC][last] > domainBCValues[i][k][last]){
                  for(int l=0; l < nodeBCValues[node][type8BC].size(); l++)
                    nodeBCValues[node][type8BC][l] = domainBCValues[i][k][l];
                }
              }//new type 8 bc
            }//it was a new bc flag for the node
          }//domain bc flag loop
        }// loop over nodes in domain
      } //loop over all domains

      std::cout << "line " << __LINE__ << std::endl; 
      // Populate the nodal bc vectors from the edge bc vectors
      //loop over all edges
      for(int i=0; i < edges.size(); i++){
        //loop over all nodes
        for(int j=0; j < nodesToDomains.size(); j++){
          int onEdge=0;
          int node=j;

          //check to see if the node is on the domains of the edge
          for(int k=0; k < nodesToDomains[node].size(); k++){
            if(nodesToDomains[node][k] == edges[i][0] 
              || nodesToDomains[node][k] == edges[i][1])
              onEdge++;
          }

          // If the node is on the edge add the edge bcs to the node's bcs
          if(onEdge == 2){
            bool newValue = true;
            //loop over all the bc flags for this edge
            for(int k=0; k < edgeBCs[i].size(); k++){
              int type8BC=-1;
              // loop over all the bc flags for this node
              for(int l=0; l < nodeBCs[node].size(); l++){
                // if the node already has this bc & value skip it
                if(nodeBCs[node][l] == edgeBCs[i][k]){
                  type8BC=l;
                  for(int m=0; m < edgeBCValues[i][k].size(); m++){
                    if(edgeBCValues[i][k][m] != nodeBCValues[node][l][m])
                      break;
                    else if(m == edgeBCValues[i][k].size()-1)
                      newValue = false;
                  }
                }
              }// node bc flags
              // if it is a new flag add it and its values to 
              // the node's vectors
              if(newValue){
                //nodes can have as many type 6 bcs as needed
                if(edgeBCs[i][k] == 6){
                  nodeBCs[node].push_back(edgeBCs[i][k]);
                  nodeBCValues[node].push_back(edgeBCValues[i][k]); 
                }
                //nodes can only have one type 8 bc which is determined
                //by a "priority" number that is the last bc value number
                if(edgeBCs[i][k] == 8){
                  int last = edgeBCValues[i][k].size()-1;
                  //the node doesn't have any type 8 bcs yet
                  if(type8BC == -1){
                    nodeBCs[node].push_back(edgeBCs[i][k]);
                    nodeBCValues[node].push_back(edgeBCValues[i][k]); 
                  }
                  //it does have a type 8 bc already & we should keep the one
                  //with the highest priority (here a higher priority means a
                  //lower number in the last value position)
                  else if(nodeBCValues[node][type8BC][last] > edgeBCValues[i][k][last]){
                    for(int l=0; l < nodeBCValues[node][type8BC].size(); l++)
                      nodeBCValues[node][type8BC][l] = edgeBCValues[i][k][l];
                  }
                }//new type 8 bc
              }//it was a new bc flag for the node
            }//loop over the bcs for the edge 
          } //if the node is on the edge
        } //loop over all the nodes
      } //loop over all edges

      std::cout << "line " << __LINE__ << std::endl; 
      // Print the nodal bcs to check
      if(verblevel > 3){
        Ostr << "nodal bcs: " << std::endl;
        for(int i=0; i < numNodes; i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < nodeBCs[i].size(); j++){
            Ostr << nodeBCs[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }

      std::cout << "line " << __LINE__ << std::endl; 
      // Print the nodal bc values to check
      if(verblevel > 3){
        Ostr << "nodal bc values: " << std::endl;
        for(int i=0; i < numNodes; i++){
          Ostr << i+1 << ": " << std::endl;
          for(int j=0; j < nodeBCValues[i].size(); j++){
            for(int k=0; k < nodeBCValues[i][j].size(); k++){
              Ostr << "    " << nodeBCValues[i][j][k] << " ";
            }
            Ostr << std::endl;
          }
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");


      std::cout << "line " << __LINE__ << std::endl; 
      // Populate the element bc vectors from the nodal bc vectors
      elemBCs.resize(numElems);
      elemBCValues.resize(numElems);
      //loop over all nodes
      for(int i=0; i < numNodes; i++){
        std::cout << "node " << i+1 << std::endl;
        // loop over all elements for the node
        for(int j=0; j < nodesToElems[i].size(); j++){
          int elem;
          elem = nodesToElems[i][j]-1;
          std::cout << "elem " << elem+1 << std::endl;
          // elems will get all bc flags for all the nodes
          // they contain
          // loop over all the bc flags for this node
          std::cout << "line " << __LINE__ << std::endl; 
          for(int k=0; k < nodeBCs[i].size(); k++){
            bool newValue = true;
            std::cout << "nodeBC " << nodeBCs[i][k] << std::endl;
            // loop over all the bc flags for this elem
            for(int l=0; l < elemBCs[elem].size(); l++){
              std::cout << "elemBC " << elemBCs[elem][l] << std::endl;
              // if the elem already has this bc & value skip it
              if(elemBCs[elem][l] == nodeBCs[i][k]){
                for(int m=0; m < nodeBCValues[i][k].size(); m++){
                  if(nodeBCValues[i][k][m] != elemBCValues[elem][l][m])
                    break;
                  else if(m == nodeBCValues[i][k].size()-1)
                    newValue = false;
                }
              }
            }// elem bc flags
            // if it is a new flag add it and its values to 
            // the elem's vectors
            if(newValue){
              std::cout << "newValue" << std::endl;
              elemBCs[elem].push_back(nodeBCs[i][k]);
              elemBCValues[elem].push_back(nodeBCValues[i][k]); 
            }//it was a new bc flag for the elem
          }//node bc flag loop
          std::cout << "line " << __LINE__ << std::endl; 
        }// loop over elements for the node
      } //loop over all nodes


      std::cout << "line " << __LINE__ << std::endl; 
      // Print the element bcs for a check
      if(verblevel > 3){
        Ostr << "elemBCs:" << std::endl;
        for(int i=0; i < numElems; i++){
          Ostr << i+1 << ": ";
          for(int j=0; j < elemBCs[i].size(); j++){
            Ostr << elemBCs[i][j] << " ";
          }
          Ostr << std::endl;
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");
      std::cout << "line " << __LINE__ << std::endl; 

      // Print the element bc values for a check
      if(verblevel > 3){
        Ostr << "elem BC values:" << std::endl;
        for(int i=0; i < numElems; i++){
          Ostr << i+1 << ": " << std::endl;
          for(int j=0; j < elemBCValues[i].size(); j++){
            for(int k=0; k < elemBCValues[i][j].size(); k++){
              Ostr << "    " <<  elemBCValues[i][j][k] << " ";
            }
            Ostr << std::endl;
          }
        }
      }
      StdOut(Ostr.str());
      Ostr.str("");

      std::cout << "line " << __LINE__ << std::endl; 
      //Write the output file
      error = WriteOutput();
      std::cout << "line " << __LINE__ << std::endl; 
      if(error == 1){
        std::ostringstream ErrOstr;
        ErrOstr << "Error writing output file." << std::endl;
        StdOut(Ostr.str());
        Ostr.str("");
        ErrOut(ErrOstr.str());
        ErrOstr.str("");
        exit(1);
      }
      StdOut(Ostr.str());
      Ostr.str("");
      //
      // ---------- Program End -----------------
      

      // Update the stacker/profiler that we are exiting 
      // this function.
      FunctionExit("Run");
      std::cout << "line " << __LINE__ << std::endl; 
      // return 0 for success
      return(0);
    };
  };
};
