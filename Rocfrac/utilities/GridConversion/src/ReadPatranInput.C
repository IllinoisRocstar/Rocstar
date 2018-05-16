#include "DriverProgram.H"

namespace GridConversion{ namespace DriverProgram{

  int SerialProgram::ReadPatranInput(){

      std::ostringstream Ostr;

      FunctionEntry("ReadPatranInput");
      // Open the specified input file for reading
      Inf.open(input_name.c_str());
      if(!Inf){
        // If the input file failed to open, notify to
        // the error stream and return non-zero
        std::ostringstream ErrOstr;
        ErrOstr << "Could not open input file, '" 
             << input_name << "'.\n";
        StdOut(Ostr.str());
        Ostr.str("");
        ErrOut(ErrOstr.str());
        ErrOstr.str("");
        // don't forget to tell the profiler/stacker the
        // function is exiting.
        FunctionExit("ReadPatranInput");
        return(1);
      }

      std::string line;
      double valueD;
      unsigned int valueI;
      std::string valueS;
      std::stringstream ss;

      //Read Input file in Como (Patran) format from Gridgen
      std::getline(Inf,line);
      std::getline(Inf,line);
      std::getline(Inf,line);
      ss << line;
      ss >> valueD >> valueD >> valueD >> valueD >> numNodes >> numElems;
      std::getline(Inf,line);

      //Resize nodes vector
      nodes.resize(3*numNodes);

      //Print for check
      if(verblevel > 1){
        Ostr << "Number of nodes = " << numNodes << std::endl;
        Ostr << "Number of elements = " << numElems << std::endl;
      }

      for(int i=0; i < numNodes; i++){
        std::getline(Inf,line);
        std::getline(Inf,line);
        ss.clear();
        ss.str("");
        ss << line;
        ss >> nodes[3*i + 0] >> nodes[3*i + 1] >> nodes[3*i + 2];
        std::getline(Inf,line);
      }

      for(int i=0; i < numElems; i++){
        std::getline(Inf,line);
        ss.clear();
        ss.str("");
        ss << line;
        ss >> valueI >> valueI >> valueI;

        //Check that valueI is a valid element shape
        if(valueI <= 1 || valueI == 6 || valueI > 9){
          std::ostringstream ErrOstr;
          ErrOstr << "element shape value of " << valueI << " is invalid." << std::endl
                    << "Valid shape values:" << std::endl
                    << "  2: bar" << std::endl
                    << "  3: triangle" << std::endl
                    << "  4: quadrilateral" << std::endl
                    << "  5: tetrahedron" << std::endl
                    << "  7: triangular prism" << std::endl
                    << "  8: hexahedron" << std::endl
                    << "  9: pyramid" << std::endl;
          StdOut(Ostr.str());
          Ostr.str("");
          ErrOut(ErrOstr.str());
          ErrOstr.str("");
          exit(1);
        }
        
        //Make sure the elements all have the same shape
        //For now we only support meshes with one shape of element
        //(this is what Rocfrac requires as well).
        if(i == 0){
          elemShape = valueI;
          //Print for check
          if(verblevel > 1)
            Ostr << "Element shape value " << elemShape << ": " << shapes[elemShape] << std::endl;
        }
        else{
          if(valueI != elemShape){
            std::ostringstream ErrOstr;
            ErrOstr << "Meshes must have only one element shape!" << std::endl
                 << "Element " << i << " has shape value " << valueI << " but "
                 << "previous elements have shape value "  << elemShape << std::endl; 
            StdOut(Ostr.str());
            Ostr.str("");
            ErrOut(ErrOstr.str());
            ErrOstr.str("");
            exit(1);
          }
        }

        std::getline(Inf,line);
        ss.clear();
        ss.str("");
        ss << line;
        ss >> valueI;    

        //Check that valueI is a valid number of nodes per element
        if(valueI < 2 || valueI == 7 || valueI == 9 || valueI > 10){
          std::ostringstream ErrOstr;
          ErrOstr << "Error: " << valueI << " number of nodes per element is invalid." << std::endl
                    << "Valid number of nodes:" << std::endl
                    << "  2: bar" << std::endl
                    << "  3: triangle" << std::endl
                    << "  4: quadrilateral or tetrahedron" << std::endl
                    << "  5: pyramid" << std::endl
                    << "  6: triangular prism" << std::endl
                    << "  8: hexahedron" << std::endl
                    << " 10: tehtrahedron (higher order)" << std::endl;
          StdOut(Ostr.str());
          Ostr.str("");
          ErrOut(ErrOstr.str());
          ErrOstr.str("");
          exit(1);
        }
 
        //Make sure the elements are all of the same type (number of nodes)
        //For now we only support meshes with one type of element
        //(this is what Rocfrac requires as well).
        if(i == 0){
          numNodesPerElem = valueI;
          //Print for check
          if(verblevel > 1)
            Ostr << "Number of nodes per element " << numNodesPerElem << "." << std::endl;
        }
        else{
          if(valueI != numNodesPerElem){
            std::ostringstream ErrOstr;
            ErrOstr << "Meshes must have elements with the same number of nodes!" << std::endl
                 << "Element " << i << " has " << valueI << " nodes but "
                 << "previous elements have "  << numNodesPerElem << " nodes." <<  std::endl; 
            StdOut(Ostr.str());
            Ostr.str("");
            ErrOut(ErrOstr.str());
            ErrOstr.str("");
            exit(1);
          }
        }

        //Resize the elems array
        elems.resize(numElems*numNodesPerElem);

        //Read in the nodes for the element
        std::getline(Inf,line);
        ss.clear();
        ss.str("");
        ss << line;
        for(int j=0; j < numNodesPerElem; j++){
          ss >> elems[i*numNodesPerElem + j];
        }
      }

      // Read the boundary domain info
      while(getline(Inf,line)){
        int packet, lines=0, domainNodes;
        std::vector<unsigned int> oneDomain;

        ss.clear();
        ss.str("");
        ss << line;
        ss >> packet >> valueI >> domainNodes >> lines;
        if(packet == 99)//Patran exit criterion
          break;       
 
        // Read the nodes that are on the domain
        std::getline(Inf,line);
        for(int i=0; i < lines-1; i++){
          std::getline(Inf,line);
          ss.clear();
          ss.str("");
          ss << line;
          int j=0;
          while(ss >> valueI){
            if(j%2 == 1)
              oneDomain.push_back(valueI);
            j++;
          }
        }
        domains.push_back(oneDomain);
      } 

      //Close input file
      Inf.close();

      // Open the bc input file for reading
      Inf.open(bc_input_name.c_str());
      if(!Inf){
        // If the input file failed to open, notify to
        // the error stream and return non-zero
        std::ostringstream ErrOstr;
        ErrOstr << "Could not open input file, '" 
             << bc_input_name << "'.\n";
        StdOut(Ostr.str());
        Ostr.str("");
        ErrOut(ErrOstr.str());
        ErrOstr.str("");
        // don't forget to tell the profiler/stacker the
        // function is exiting.
        FunctionExit("ReadPatranInput");
        return(1);
      }

      domainBCs.resize(domains.size());
      domainBCValues.resize(domains.size());
      while(std::getline(Inf,line)){
        int numFlags=0, domainNum;
        std::string edgeOrDomain;
        bool domain = true;
        ss.clear();
        ss.str("");
        ss << line;
        ss >> edgeOrDomain >> domainNum;
        domainNum--;
        std::cout << "edgeOrDomain = " << edgeOrDomain << std::endl;
        //bcs for an edge
        if(edgeOrDomain == "edge" || edgeOrDomain == "Edge" 
           || edgeOrDomain == "EDGE"){
          domain = false;
          ss >> valueI;
          edges[domainNum].push_back(valueI);
          ss >> valueI;
          edges[domainNum].push_back(valueI);
        }
        //tell the total number of edges with a bc
        else if(edgeOrDomain == "edges" || edgeOrDomain == "Edges"
                || edgeOrDomain == "EDGES"){
          edges.resize(domainNum+1);
          edgeBCs.resize(domainNum+1);
          edgeBCValues.resize(domainNum+1);
          continue;
        }
        //bcs for a domain
        else if(edgeOrDomain != "domain" && edgeOrDomain != "Domain"
           && edgeOrDomain != "DOMAIN"){
          std::ostringstream ErrOstr;
          ErrOstr << "First word of every line in bc file must be 'domain,'" 
                  << std::endl << "'edge', or 'edges.'" << std::endl;
          StdOut(Ostr.str());
          Ostr.str("");
          ErrOut(ErrOstr.str());
          ErrOstr.str("");
          // don't forget to tell the profiler/stacker the
          // function is exiting.
          FunctionExit("ReadPatranInput");
          return(1);
        }
      
        std::cout << "edges.size() = " << edges.size() << std::endl; 
        std::cout << "edgeBCs.size() = " << edgeBCs.size() << std::endl; 
        std::cout << "edgeBCValues.size() = " << edgeBCValues.size() << std::endl; 
        std::cout << "domainNum = " << domainNum << std::endl;
        ss >> numFlags;
        if(domain) 
          domainBCValues[domainNum].resize(numFlags);
        else 
          edgeBCValues[domainNum].resize(numFlags);
       
        std::cout << "line " << __LINE__ << std::endl;
 
        for(int j=0; j < numFlags; j++){
          ss >> valueI;
          if(domain)
            domainBCs[domainNum].push_back(valueI);
          else
            edgeBCs[domainNum].push_back(valueI);
          
        std::cout << "line " << __LINE__ << std::endl;
          if(valueI == 8){
            for(int i=0; i < 6; i++){
              ss >> valueD;
              if(domain)
                domainBCValues[domainNum][j].push_back(valueD);            
              else
                edgeBCValues[domainNum][j].push_back(valueD);            
            } 
          }
          else if(valueI == 6){
            ss >> valueD;
            if(domain)
              domainBCValues[domainNum][j].push_back(valueD);
            else
              edgeBCValues[domainNum][j].push_back(valueD);
          }
          else{
            // Error for invalid bc type
            std::ostringstream ErrOstr;
            ErrOstr << "Invalid boundary condition type of " << valueI 
                    << " for boundary " << domainBCs.size() << "." << std::endl;
            StdOut(Ostr.str());
            Ostr.str("");
            ErrOut(ErrOstr.str());
            ErrOstr.str("");
            // don't forget to tell the profiler/stacker the
            // function is exiting.
            FunctionExit("ReadPatranInput");
            return(1);
          }
        std::cout << "line " << __LINE__ << std::endl;
        }
      }

        std::cout << "line " << __LINE__ << std::endl;
      //Close the bc input file
      Inf.close();

      StdOut(Ostr.str());
      Ostr.str("");
      FunctionExit("ReadPatranInput");

    return 0;
 
  } //ReadInput function

}; //DriverProgram namespace
}; //GridConversion namespace
