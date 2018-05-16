/// 
/// @file
/// @ingroup gridconversion_group
/// @brief Main for example driver program.
/// @author Jessica Kress (jkress@illinoisrocstar.com)
/// @date 
/// 
#include "DriverProgram.H"

typedef GridConversion::DriverProgram::SEProgramType ProgramType;

int main(int argc,char *argv[])
{
  return(GridConversion::DriverProgram::Driver<ProgramType>(argc,argv));
}
