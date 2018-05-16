/// 
/// @file
/// @ingroup gridconversion_group
/// @brief Main for example serial program.
/// @author Mike Campbell (mtcampbe@illinois.edu)
/// @date 
/// 
#include "ExampleProgram.H"

typedef GridConversion::ExampleProgram::SEProgramType ProgramType;

int main(int argc,char *argv[])
{
  return(GridConversion::ExampleProgram::Driver<ProgramType>(argc,argv));
}
