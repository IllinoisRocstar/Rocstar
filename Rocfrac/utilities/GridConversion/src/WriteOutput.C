#include <iomanip>
#include "DriverProgram.H"

namespace GridConversion{ namespace DriverProgram{

  int SerialProgram::WriteOutput(){

      std::ostringstream Ostr;
      std::stringstream ss;


      // Open the specified output file for writing
      bool use_outfile = false;
      if(!output_name.empty()){
        use_outfile = true;
        Ouf.open(output_name.c_str());
        if(!Ouf){
          // If the output file failed to open, notify
          // to error stream and return non-zero
          std::ostringstream ErrOstr;
          ErrOstr << "Error: Unable to open output file, " << output_name << ".";
          StdOut(Ostr.str());
          Ostr.str("");
          ErrOut(ErrOstr.str());
          ErrOstr.str("");
          // don't forget to tell the profiler/stacker the
          // function is exiting.
          FunctionExit("Run");
          return(1);
        }
      }

      // Write mesh info to output file or screen
      ss.clear();
      ss.str("");

      //Write Patran packet 25 Intro (unused) 
      ss << "unused line" << std::endl;
      ss << "unused line" << std::endl;
      //Write Patran packet 26 
      ss << "26 1 1 1 " << numNodes << " " << numElems << " 1 1 1" << std::endl;
      ss << "unused line" << std::endl;

      unsigned int ePos, expPos; 
      std::stringstream numSS;
      std::string output, exponent;
      //Write Patran packet 1 node coordinates
      for(int i=0; i < numNodes; i++){
        ss << "1 " << i+1 << std::endl;
      
        if(fabs(nodes[3*i+0]) > 1.0e-10 || nodes[3*i+0] == 0.0){
        numSS << std::setprecision(9) <<  std::scientific << nodes[3*i + 0]; 
        ePos = numSS.str().find("e");
        expPos = numSS.str().find_last_of("0");
        //std::cout << "numSS.str() = " << numSS.str() << std::endl;
        //std::cout << "expPos = " << expPos << std::endl;
        //std::cout << "numSS.str().size() = " << numSS.str().size() << std::endl;
        if(expPos == numSS.str().size()-1)
          exponent = "0";
        else
          exponent = numSS.str().substr(expPos+1,numSS.str().size());       
        output = numSS.str().substr(0,ePos+2) + exponent;
        }
        else{//we have to be careful with numbers that have two digit exponents
        //in scientific notation because Rocfrac expects all numbers to take a  
        //specific number of columns.
          if(nodes[3*i+0] < 0.0)
            numSS << std::setprecision(8) << std::scientific << nodes[3*i + 0];
          else
            numSS << std::setprecision(9) << std::scientific << nodes[3*i + 0];
         
          output = numSS.str();
        }
        ss << std::setw(16) << output;
        numSS.str("");
        numSS.clear();

        if(fabs(nodes[3*i+1]) > 1.0e-10 || nodes[3*i+1] == 0.0){
        numSS << std::setprecision(9) <<  std::scientific << nodes[3*i + 1]; 
        ePos = numSS.str().find("e");
        expPos = numSS.str().find_last_of("0");
        //std::cout << "numSS.str() = " << numSS.str() << std::endl;
        //std::cout << "expPos = " << expPos << std::endl;
        //std::cout << "numSS.str().size() = " << numSS.str().size() << std::endl;
        if(expPos == numSS.str().size()-1)
          exponent = "0";
        else
          exponent = numSS.str().substr(expPos+1,numSS.str().size());       
        output = numSS.str().substr(0,ePos+2) + exponent;
        }
        else{//we have to be careful with numbers that have two digit exponents
        //in scientific notation because Rocfrac expects all numbers to take a  
        //specific number of columns.
          if(nodes[3*i+1] < 0.0)
            numSS << std::setprecision(8) << std::scientific << nodes[3*i + 1];
          else
            numSS << std::setprecision(9) << std::scientific << nodes[3*i + 1];
         
          output = numSS.str();
        }
        ss << std::setw(16) << output;
        numSS.str("");
        numSS.clear();

        if(fabs(nodes[3*i+2]) > 1.0e-10 || nodes[3*i+2] == 0.0){
        numSS << std::setprecision(9) <<  std::scientific << nodes[3*i + 2]; 
        ePos = numSS.str().find("e");
        expPos = numSS.str().find_last_of("0");
        //std::cout << "numSS.str() = " << numSS.str() << std::endl;
        //std::cout << "expPos = " << expPos << std::endl;
        //std::cout << "numSS.str().size() = " << numSS.str().size() << std::endl;
        if(expPos == numSS.str().size()-1)
          exponent = "0";
        else
          exponent = numSS.str().substr(expPos+1,numSS.str().size());       
        output = numSS.str().substr(0,ePos+2) + exponent;
        }
        else{//we have to be careful with numbers that have two digit exponents
        //in scientific notation because Rocfrac expects all numbers to take a  
        //specific number of columns.
          if(nodes[3*i+2] < 0.0)
            numSS << std::setprecision(8) << std::scientific << nodes[3*i + 2];
          else
            numSS << std::setprecision(9) << std::scientific << nodes[3*i + 2];
         
          output = numSS.str();
        }
        ss << std::setw(16) << output << std::endl;
        numSS.str("");
        numSS.clear();

        ss << "unused line" << std::endl;
      }
    
      //Write Patran packet 2 element connectivities 
      for(int i=0; i < numElems; i++){
        ss << "2 " << i+1 << " " << elemShape << std::endl
           << "       " <<  numNodesPerElem << "       0"
           << "       0" << std::endl;
        for(int j=0; j < numNodesPerElem; j++){
          ss << elems[i*numNodesPerElem + j] << " ";
        }
        ss << std::endl;
      }

      //Write Patran packet 8 (all nodes with this boundary type)
      for(int i=0; i < numNodes; i++){
        for(int k=0; k < nodeBCs[i].size(); k++){
          if(nodeBCs[i][k] == 8){
            ss << "8 " << i+1 << " 1 2" << std::endl;
            ss << "       0";
            for(int j=0; j < 3; j++)
              ss << int(nodeBCValues[i][k][j]);
            ss << "000" <<  std::endl;
            for(int j=3; j < 6; j++){
              if(nodeBCValues[i][k][j-3] != 0)
                ss << " " << nodeBCValues[i][k][j];
            }
            ss << std::endl;
          }
        }
      }
   

      //Write Patran packet 6 (all elements with this boundary type)
      //loop over all the elements
      //std::cout << "bcs of type 6: " << std::endl;
      for(int i=0; i < numElems; i++){
        //std::cout << "elem " << i+1 << ": " << std::endl;
        //loop over the domains for that element (we only want
        //to print one domain's bcs for the element at a time
        for(int doms=0; doms < elemsToDomains[i].size(); doms++){
          int domain=elemsToDomains[i][doms]-1;
          //std::cout << "domain " << domain+1 << std::endl;
          //loop over the boundary flags for that domain
          for(int j=0; j < domainBCs[domain].size(); j++){ 
            int bcCount=0;
            if(domainBCs[domain][j] == 6){
              //std::cout << "bc " << j+1 << ": " << std::endl;
              ss << "6 " << i+1 << " 1 2 0 0 0 0 0" << std::endl;
              ss << "111100000";
              //loop over the nodes for that element
              for(int k=0; k < numNodesPerElem; k++){
                bool onBC=false;
                int node = elems[i*numNodesPerElem + k]-1;
                //std::cout << "node " << k+1 << " (" << node+1 << ")" << std::endl; 
                //loop over the bc flags for that node
                for(int l=0; l < nodeBCs[node].size(); l++){
                  if(nodeBCs[node][l] == 6){
                    //check if the bcs have the same values
                    for(int m=0; m < nodeBCValues[node][l].size(); m++){
                      if(nodeBCValues[node][l][m] != domainBCValues[domain][j][m])
                        break;
                      else if(m == nodeBCValues[node][l].size()-1){
                        onBC=true;
                        bcCount++;
                        //std::cout << "has bc" << std::endl;
                      }
                    }
                  }
                }//loop over node bcs
                if(onBC){
                  //now check that the node is on that domain
                  int domCount;
                  for(domCount=0; domCount < nodesToDomains[node].size(); domCount++){
                    if(nodesToDomains[node][domCount] == domain+1){
                      //std::cout << "node on domain (writing it)" << std::endl;
                      ss << "1";
                      break;
                    }
                  }
                  if(domCount == nodesToDomains[node].size())
                    ss << "0";
                }
                else
                  ss << "0";
              }//loop over nodes for the element
              //Fill in with 0's at the end if we are using tet elements
              if(numNodesPerElem == 4)
                ss << "0000";
              ss << std::endl;
              for(int m=0; m < domainBCValues[domain][j].size(); m++)
                ss << domainBCValues[domain][j][m] << " ";
              ss << std::endl;
            }//if it is bc type 6
          }//loop over bc flags for domain
        }//loop over the domains for the element
      }//loop over all elements
      ss << "99 0 0 1";

      if(use_outfile){
        Ouf << ss.str();
      }  
      else
        StdOut(ss.str());
      // Close output file
      Ouf.close();
 
      StdOut(Ostr.str());
      Ostr.str("");

    return 0;
 
  } //ReadOutput function

}; //DriverProgram namespace
}; //GridConversion namespace
