

/* A "C" main */

#include <stdio.h>
#include <stdlib.h>
#include <strstream.h>
#include <fstream.h>
#include "Mesh.hpp"
#include "Element.hpp"
#include "Face.hpp"

int main( int argc, char** argv){

  if( argc < 4 || argc > 5){
    cerr << "Wrong format. Should be:\n"
	 << argv[0] << " [file name] [type1] [type2] <type3>\n";
    exit(-1);
  }
  Mesh mesh;
  ifstream ifs( argv[1] );
  ifs >> mesh;
  int type1 = atoi( argv[2] );
  int type2 = atoi( argv[3] );
  int type3 = -1;
  if( argc == 5 ){
    type3 = atoi( argv[4] );
  }
  
  mesh.addCohesive( type1, type2, type3);

  char line[80];
  ostrstream str(line,sizeof(line));
  str << argv[1] <<"." << "coh"<< '\0';
  ofstream ocohesive( str.str() );
  ocohesive << mesh;  

  char line1[80];
  ostrstream str1(line1,sizeof(line1));
  str1 << argv[1] <<"." << "bound"<< '\0';
  ofstream obound( str1.str() );
  mesh.write_boundary( obound );  
  
  return(0);
}



