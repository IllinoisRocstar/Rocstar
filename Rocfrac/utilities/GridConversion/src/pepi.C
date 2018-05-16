/// 
/// @file
/// @ingroup gridconversion_group
/// @brief Main for example parallel program.
/// @author Mike Campbell (mtcampbe@illinois.edu)
/// @date 
/// 
#include "ExampleProgram.H"

typedef GridConversion::ExampleProgram::PEProgramType ProgramType;

int main(int argc,char *argv[])
{
  return(GridConversion::ExampleProgram::Driver<ProgramType>(argc,argv));
}
