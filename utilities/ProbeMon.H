#ifndef _IR_PLOT_UTILS_H_
#define _IR_PLOT_UTILS_H_
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstdio>

namespace IR {
  class flowprobe : public std::vector<std::vector<double> > 
  {
  public:
    enum probedata {TIME=0,RHO,U,V,W,P,T,NPROBEDATA};
    flowprobe()
    {
      this->resize(NPROBEDATA);
    };
    flowprobe(const flowprobe &inFlowProbe)
    {
      this->resize(NPROBEDATA);
      (*this)[TIME] = inFlowProbe[TIME];
      (*this)[RHO]  = inFlowProbe[RHO];
      (*this)[U]    = inFlowProbe[U];
      (*this)[V]    = inFlowProbe[V];
      (*this)[W]    = inFlowProbe[W];
      (*this)[P]    = inFlowProbe[P];
      (*this)[T]    = inFlowProbe[T];
    };
    void Clear()
    {
      this->resize(0);
      this->resize(NPROBEDATA);
    };
    unsigned int ReadProbeData(std::istream &inStream)
    {
      Clear();
      std::string probeLine;
      while(std::getline(inStream,probeLine)){
        if(probeLine[0] == '#')
          continue;
        std::istringstream Instr(probeLine);
        double time  = 0;
        double rho   = 0;
        double rhou  = 0;
        double rhov  = 0;
        double rhow  = 0;
        double press = 0;
        double temp  = 0;
        Instr >> time >> rho >> rhou >> rhov >> rhow >> press >> temp;
        (*this)[TIME].push_back(time);
        (*this)[RHO].push_back(rho);
        (*this)[U].push_back(rhou/rho);
        (*this)[V].push_back(rhov/rho);
        (*this)[W].push_back(rhow/rho);
        (*this)[P].push_back(press);
        (*this)[T].push_back(temp);
      }
      return((*this)[TIME].size());
    };
    void WriteData(std::ostream &outStream,bool writeLoc = false)
    {
      unsigned int nTimes = (*this)[TIME].size();
      for(unsigned int i = 0;i < nTimes;i++){
        if(writeLoc){
          outStream << locX << " " << locY << " " << locZ << " ";
        }
        for(unsigned int j = 0;j < NPROBEDATA;j++){
          outStream << (*this)[j][i] << " ";
        }
        outStream << std::endl;
      }
    };
    void SetLocation(double xLoc,double yLoc,double zLoc)
    {
      locX = xLoc;
      locY = yLoc;
      locZ = zLoc;
    };
  protected:
    double locX;
    double locY;
    double locZ;
  };
  class derivedprobe : public std::vector<std::vector<double> >
  {

  };
  inline int GNUPlot(const std::string &commands)
  {
    FILE *fp = popen("gnuplot -persist","w");
    fprintf(fp,commands.c_str());
    pclose(fp);
    return(0);
  }
}
#endif
