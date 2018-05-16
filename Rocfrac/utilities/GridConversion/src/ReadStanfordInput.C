#include "DriverProgram.H"

namespace GridConversion{ namespace DriverProgram{

  int SerialProgram::ReadStanfordInput(){

      FunctionEntry("ReadStanfordInput");
      std::ostringstream Ostr;

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
        FunctionExit("ReadStanfordInput");
        return(1);
      }

      std::string line;
      double valueD;
      unsigned int valueI;
      std::string valueS;
      std::stringstream ss;

      //Read Input file in Stanford format from Poinwise
      for(int i=0; i < 8; i++)
        std::getline(Inf,line);
      ss.clear();
      ss.str("");
      ss << line;
      ss >> valueS >> numElems;

      //Print for check
      if(verblevel > 1){
        Ostr << "Number of nodes = " << numNodes << std::endl;
        Ostr << "Number of elements = " << numElems << std::endl;
      }
      for(int i=0; i < numElems; i++){
        std::getline(Inf,line);
        ss.clear();
        ss.str("");
        ss << line;
        ss >> valueI;

        //Make sure the elements all have the same shape
        //For now we only support meshes with one shape of element
        //(this is what Rocfrac requires as well).
        if(i == 0){
          elemShape = valueI;
          //Print for check
          if(verblevel > 1)
            Ostr << "Element shape value: " << elemShape << std::endl;
          if(elemShape == 10)
            numNodesPerElem = 4;
          else if(elemShape == 12)
            numNodesPerElem = 8;
          else{
            std::ostringstream ErrOstr;
            ErrOstr << "For Stanford format input only element shapes " << std::endl
                    << "10 and 12 are supported! Input file has shape " 
                    << elemShape << "!" << std::endl;
            StdOut(Ostr.str());
            Ostr.str("");
            ErrOut(ErrOstr.str());
            ErrOstr.str("");
            exit(1);
          }
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

        //Read in 8 nodes for shape 12 (hexahedral) or 4 nodes for shape 10 (tetrahedral) 
        for(int j=0; j < numNodesPerElem; j++){
          ss >> valueI;
          //We are starting are node count at 1
          valueI++;
          elems.push_back(valueI);
        }
      }

      for(int i=0; i < 4; i++)
        std::getline(Inf,line);
      ss.clear();
      ss.str("");
      ss << line;
      ss >> valueS >> numNodes;

      //Resize nodes vector
      nodes.resize(3*numNodes);

      //Print for check
      if(verblevel > 1){
        Ostr << "Number of nodes = " << numNodes << std::endl;
        Ostr << "Number of elements = " << numElems << std::endl;
      }

      for(int i=0; i < numNodes; i++){
        std::getline(Inf,line);
        ss.clear();
        ss.str("");
        ss << line;
        ss >> nodes[3*i + 0] >> nodes[3*i + 1] >> nodes[3*i + 2];
      }


      // Read the boundary domain info
      for(int i=0; i < 4; i++)
        getline(Inf,line);
      ss.clear();
      ss.str("");
      ss << line;
   
      ss >> valueS >> numDomains;

      // loop over the domains
      for(int i=0; i < numDomains; i++){
        int packet, lines=0, domainNodes;
        std::vector<unsigned int> oneDomain;
        int numNodesPerFaceElem;

        getline(Inf,line);
        getline(Inf,line);
        ss.clear();
        ss.str("");
        ss << line;

        ss >> valueS >> lines; 
      
        //loop over number of faces on domain
        for(int j=0; j < lines; j++){
          getline(Inf,line);
          ss.clear();
          ss.str("");
          ss << line;

          ss >> valueI;
          if(valueI == 9)
            numNodesPerFaceElem = 4;
          else if(valueI == 5)
            numNodesPerFaceElem = 3;
          else{
            std::ostringstream ErrOstr;
            ErrOstr << "Face elements must be of type 9 (quad) or 5 (tri)!" << std::endl
                 << "Face type " << valueI << " was read!" << std::endl;
            StdOut(Ostr.str());
            Ostr.str("");
            ErrOut(ErrOstr.str());
            ErrOstr.str("");
            exit(1);
          }

          //loop over the number of nodes for the face
          for(int k=0; k < numNodesPerFaceElem; k++){
            ss >> valueI;
            //We are starting our node count at 1
            valueI++;
            //check to see if the node has been added to
            //the domain already
            bool newValue=true;
            for(int l=0; l < oneDomain.size(); l++){
              if(oneDomain[l] == valueI){
                newValue=false;
                break;
              }
            }
            //Add the value to the domain if its new
            if(newValue){
              oneDomain.push_back(valueI);
            }
          }
        }//faces on domain loop
        //Add the one domain we read in to the vector of all the domains
        domains.push_back(oneDomain);
      }//domain loop

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
        FunctionExit("ReadStanfordInput");
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
        //std::cout << "edgeOrDomain = " << edgeOrDomain << std::endl;
        //std::cout << "domainNum = " << domainNum << std::endl;
        domainNum--;
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
          FunctionExit("ReadStanfordInput");
          return(1);
        }
      
        ss >> numFlags;
        //std::cout << "numFlags = " << numFlags << std::endl;
        if(domain) 
          domainBCValues[domainNum].resize(numFlags);
        else 
          edgeBCValues[domainNum].resize(numFlags);
       
        for(int j=0; j < numFlags; j++){
          ss >> valueI;
          //std::cout << "domain type = " << valueI << std::endl;
          if(domain)
            domainBCs[domainNum].push_back(valueI);
          else
            edgeBCs[domainNum].push_back(valueI);
         
          //std::cout << "values = " << std::endl; 
          if(valueI == 8){
            for(int i=0; i < 7; i++){
              ss >> valueD;
              //std::cout << valueD << " ";
              if(domain)
                domainBCValues[domainNum][j].push_back(valueD);            
              else
                edgeBCValues[domainNum][j].push_back(valueD);            
            } 
          }
          else if(valueI == 6){
            ss >> valueD;
            //std::cout << valueD << " ";
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
            FunctionExit("ReadStanfordInput");
            return(1);
          }
          //std::cout << std::endl;
        }
      }

      //Close the bc input file
      Inf.close();

      StdOut(Ostr.str());
      Ostr.str("");

      FunctionExit("ReadStanfordInput");
    return 0;
 
  } //ReadInput function

}; //DriverProgram namespace
}; //GridConversion namespace
