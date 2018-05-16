#include "ProbeMon.H"

std::string BuildProbeFileName(const std::string &caseName,int probeNumber)
{
  std::ostringstream outString;
  outString << caseName << ".prb_"
            << (probeNumber >= 1000 ? "" :
                (probeNumber >= 100 ? "0" :
                 (probeNumber >= 10 ? "00" :
                  "000"))) << probeNumber;
  return(outString.str());
}

int 
main(int argc,char *argv[])
{
  if(argc < 2) return(1);
  std::string caseName(argv[1]);

  std::string configFileName(caseName+".probe_config");
  std::ifstream configIn;
  configIn.open(configFileName.c_str());
  if(!configIn){
    return(1);
  }
  int nProbes = 0;
  configIn >> nProbes;
  std::vector<double> probeLocations;
  double location = 0.0;
  while(configIn >> location)
    probeLocations.push_back(location);
  configIn.close();
  int nProbeCoordinates = probeLocations.size();
  assert(nProbeCoordinates == nProbes*3);
  std::vector<double>::iterator coordinateIterator = probeLocations.begin();

  std::ostringstream htmlOut;
  htmlOut << "<html>" << std::endl
	  << "  <head>" << std::endl
	  << "    <title>" 
	  << caseName << " Probe Monitor" 
	  << "</title>" << std::endl
	  << "    <meta http-equiv=\"refresh\" content=\"300\" />" 
	  << std::endl
	  << "  </head>" << std::endl
	  << "  <body>" << std::endl
	  << "    <p>" << std::endl
	  << "      <h1><center>" << caseName 
	  << " Probe Monitor</center></h1>" << std::endl
	  << "    </p>" << std::endl << std::endl
	  << "    <hr>" << std::endl << std::endl;

  
  std::vector<IR::flowprobe> probesData;
  for(int iProbe = 1;iProbe <= nProbes;iProbe++){

    double probeLocationX = *coordinateIterator++;
    double probeLocationY = *coordinateIterator++;
    double probeLocationZ = *coordinateIterator++;

    htmlOut << "    <p>" << std::endl
	    << "      <h2><center>Probe " 
	    << iProbe << "</center></h2>" << std::endl
	    << "      <img src=\"velocity_plot_" << iProbe << ".png\" "
	    << "style=\"float: left; width: 50%\">" << std::endl 
	    << "      <img src=\"pandt_plot_" << iProbe << ".png\" "
	    << "style=\"float: left; width: 50%\">" << std::endl
	    << "    </p>" << std::endl << std::endl 
	    << "    <hr>" << std::endl
	    << std::endl;
    
    std::string probeFileName(BuildProbeFileName(caseName,iProbe));
    std::fstream probeFile;
    probeFile.open(probeFileName.c_str());
    if(!probeFile){
      std::cout << "probemon: Unable to open probe file: " 
		<< probeFileName << std::endl;
      return(1);
    }

    IR::flowprobe flowProbe;
    unsigned int nTimes = flowProbe.ReadProbeData(probeFile);
    probeFile.close();
    flowProbe.SetLocation(probeLocationX,probeLocationY,probeLocationZ);
    probesData.push_back(flowProbe);
    std::cout << "probemon: Got " << nTimes << " probe values from " 
              << probeFileName << "." << std::endl;


    std::ofstream plotFile;
    std::ostringstream dataStreamOut;
    dataStreamOut << "plot_data_" << iProbe << ".txt";
    std::string plotFileName(dataStreamOut.str());
    plotFile.open(plotFileName.c_str());
    if(!plotFile){
      std::cout << "probemon: Unable to open plot data file: " << plotFileName << std::endl;
      return(1);
    }
    flowProbe.WriteData(plotFile,true);
    plotFile.close();

    std::ostringstream plotCommands;
    plotCommands << "set title 'Probe " << iProbe << " Velocity'" << std::endl
                 << "set ylabel 'Velocity (m/s)" << std::endl
                 << "set term png" << std::endl
                 << "set key top right title ''" << std::endl
                 << "set output 'velocity_plot_"<< iProbe << ".png'" << std::endl
                 << "plot '" << plotFileName << "' using 4:6 w l t 'X-velocity','"
                 << plotFileName << "' "
                 << "using 4:7 w l t 'Y-velocity','" << plotFileName 
                 << "' using 4:8 w l t 'Z-velocity'"
                 << std::endl;
    IR::GNUPlot(plotCommands.str());

    plotCommands.clear();
    plotCommands.str("");
    plotCommands << "set title 'Probe " << iProbe << " Pressure and Temperature'" << std::endl
                 << "set ylabel 'Pressure (Pa)" << std::endl
                 << "set y2label 'Temperature (K)'" << std::endl
                 << "set term png" << std::endl
                 << "set key bottom right title ''" << std::endl
                 << "set ytics nomirror" << std::endl
                 << "set y2tics" << std::endl
                 << "set output 'pandt_plot_" << iProbe << ".png'" << std::endl
                 << "plot '"<< plotFileName << "' using 4:9 w l t 'Pressure'"
                 << ",'" << plotFileName << "' using 4:10 w l t 'Temperature' axes x1y2" 
                 << std::endl;
    IR::GNUPlot(plotCommands.str());
  }

  configFileName.assign(caseName+".derived_config");
  configIn.open(configFileName.c_str());
  if(configIn){
    // Process derived quantities
    std::vector<IR::derivedprobe> derivedProbes;
    IR::derivedprobe derivedProbe;
//     while(derivedProbe << configIn)
//       derivedProbes.push_back(derivedProbe);
    
  }
  htmlOut << "  </body>" << std::endl
	  << "</html>" << std::endl;
  std::ofstream htmlFile;
  htmlFile.open("probemonitor.html");
  htmlFile << htmlOut.str();
  htmlFile.close();
  return(0);
}
