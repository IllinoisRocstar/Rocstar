/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************/
// This program compares two Rocflu/Rocflo probe files in some time interval 
// by reporting the L2 norm of the pressure difference.
//  
// The L2 norm of the error is computed as:
// SQRT(INTEGRAL((P1 - P2)^2))
//
// and average energy normalization is computed as:
// SQRT((INTEGRAL(P1))^2 + (INTEGRAL(P2))^2)
//
// Compile it with c++ -O3 probDiff.C -o probediff
//
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>


// This routine determines and returns the value of the probe data
// at time t.  It also returns the indices of the interval, (n0,n1),
// where time t lives. It starts at index 'start' to keep it from 
// having to search the entire curve with every call.
double
determine_probe_value(std::vector<std::pair<double, double> > &probe,
		      int &n0, int &n1, double t, int start)
{
  double TOL = 1e-20;
  int n = start;
  if(probe[n].first > t && !(std::fabs(probe[n].first - t) <= TOL)){
    std::cerr << "Interval error!" << std::endl;
    exit(1);
  }
  while(n < probe.size() && probe[n].first < t && 
  	!(std::fabs(probe[n].first - t) <= TOL))
    n++;
  if(n >= probe.size()){
    std::cerr << "Ran off the end of the probe data!" << std::endl;
    exit(1);
  }
  if(std::fabs(probe[n].first - t) <= TOL){
    n0 = n;
    n1 = n + 1;
  }
  else{
    n1 = n;
    n0 = n - 1;
  }

  // This little check will account for probe hiccups that
  // sometimes happen where we will get multiple entries
  // for a single time.
  double dt = probe[n1].first - probe[n0].first;
  if(dt < 0 || (std::fabs(dt) <= TOL))
    n1++;

  // Return the linearly interpolated probe value at time t
  double slope = (probe[n1].second - probe[n0].second)/
    (probe[n1].first - probe[n0].first);
  return(probe[n0].second + slope*(t - probe[n0].first));
}

int
main(int argc,char *argv[])
{

  // This is all the necessary argument parsing and file reading
  std::vector<std::string> args;
  int nn = 0;
  while(argv[nn])
    args.push_back(std::string(argv[nn++]));
  std::string program_name(args[0]);
  
  if(argc < 3){
    std::cerr << "Usage: " << program_name << " [-t0 time0] [-t1 time1] "
	      << "<probe file 1> <probe file 2>" << std::endl;
    exit(0);
  }

  bool minset = false;
  bool maxset = false;
  double t0 = 0.0;
  double t1 = 10000000.0;
  
  std::string timeopt;
  std::vector<std::string>::iterator ai = args.begin();
  while(ai != args.end() && !minset){
    if(*ai == "-t0"){
      ai++;
      if(ai == args.end()){
	std::cerr << program_name << ": Expected argument after -t0 option."
		  << std::endl;
	exit(1);
      }
      timeopt = *ai;
      std::istringstream Istr(timeopt);
      // This removes whitespace from the argument
      Istr >> timeopt;
      Istr.clear();
      Istr.str(timeopt);
      Istr >> t0;
      if(t0 < 0.0){
	std::cerr << program_name << ": Invalid value for t0 option, " 
		  << t0 << "." << std::endl;
	exit(1);
      }
      minset = true;
    }
    else
      ai++;
  }
  ai = args.begin();
  while(ai != args.end() && !maxset){
    if(*ai == "-t1"){
      ai++;
      if(ai == args.end()){
	std::cerr << program_name << ": Expected argument after -t1 option."
		  << std::endl;
	exit(1);
      }
      timeopt = *ai;
      // This removes whitespace from the argument
      std::istringstream Istr(timeopt);
      Istr >> timeopt;
      Istr.clear();
      Istr.str(timeopt);
      Istr >> t1;
      if(t1 < 0.0){
	std::cerr << program_name << ": Invalid value for t1 option, " 
		  << t1 << "." << std::endl;
	exit(1);
      }
      maxset = true;
    }
    else
      ai++;
  }
  if((t1 - t0) < 0.0){
    std::cerr << program_name << ": Invalid time interval, " << t0 << " - "
	      << t1 << std::endl;
    exit(1);
  }
  int acount = 1;
  while(acount < argc && args[acount][0] == '-'){
    acount++;
    acount++;
  }
  if(acount >= argc){
    std::cerr << "Usage: " << program_name << " [-t0 time0] [-t1 time1] "
	      << "<probe file 1> <probe file 2>" << std::endl;
    exit(1);
  }
  std::string probe1_file(args[acount++]);
  if(acount >= argc){
    std::cerr << "Usage: " << program_name << " [-t0 time0] [-t1 time1] "
	      << "<probe file 1> <probe file 2>" << std::endl;
    exit(1);
  }
  std::string probe2_file(args[acount]);

  std::ifstream Inf1;
  std::ifstream Inf2;
  
  Inf1.open(probe1_file.c_str());
  if(!Inf1){
    std::cerr << program_name << ": Unable to open probe file, " << probe1_file
	      << "." << std::endl;
    exit(1);
  }
  Inf2.open(probe2_file.c_str());
  if(!Inf2){
    std::cerr << program_name << ": Unable to open probe file, " << probe2_file
	      << "." << std::endl;
    exit(1);
  }
  std::vector<std::pair<double,double> > probe1;
  std::vector<std::pair<double,double> > probe2;
  int n = 0;
  double dum;
  std::cout << program_name << ": Reading probe 1 ..." << std::flush;
  while(!Inf1.eof()){
    std::pair<double,double> inpair;
    Inf1 >> inpair.first >> dum >> dum >> dum >> dum >> inpair.second
	 >> dum;
    probe1.push_back(inpair);
    n++;
  }
  n--;
  probe1.pop_back();
  Inf1.close();
  std::cout << "... done" << std::endl
	    << program_name << ": Probe1(" << n << "): (" << probe1[0].first 
	    << ", " << probe1[0].second << ") - (" << probe1[n-1].first 
	    << ", " << probe1[n-1].second << ")" << std::endl;
  n = 0;
  std::cout << program_name << ": Reading probe 2 ...";
  while(!Inf2.eof()){
    std::pair<double,double> inpair;
    Inf2 >> inpair.first >> dum >> dum >> dum >> dum >> inpair.second
	 >> dum;
    probe2.push_back(inpair);
    n++;
  }

  // The stupid thing reads too many - not sure why it does this,
  // but the next two lines fix the problem
  n--;
  probe2.pop_back();
  Inf2.close();

  int n1 = probe1.size();
  int n2 = probe2.size();
  double probe_tmin = (probe1[0].first < probe2[0].first ? 
		       probe2[0].first : probe1[0].first);
  double probe_tmax = (probe1[n1-1].first > probe2[n2-1].first ?
		       probe2[n2-1].first : probe1[n1-1].first);
  std::cout << "... done" << std::endl 
	    << program_name << ": Probe2(" << n << "): (" << probe2[0].first 
	    << ", " << probe2[0].second << ") - (" << probe2[n-1].first 
	    << ", " << probe2[n-1].second << ")" << std::endl;
  if(minset && (probe1[0].first > t0 || probe2[0].first > t0))
    std::cout << program_name << ": t0(" << t0 << ") is less than the "
	      << " earliest common probe time: " << probe_tmin 
	      << ", resetting." << std::endl;
  if(t0 < probe_tmin)
    t0 = probe_tmin;
  if(maxset && (probe1[n1-1].first < t1 || probe2[n2-1].first < t1))
    std::cout << program_name << ": t1(" << t1 << ") is greater than the "
	      << " latest common probe time: " << probe_tmax << ", resetting." 
	      << std::endl;
  if(t1 > probe_tmax)
    t1 = probe_tmax;
  std::cout << program_name << ": Comparing probes in time interval (" << t0 
	    << ", " << t1 << ")" << std::endl;

  // Now we have two probe data arrays.  Each array entry is a pair of doubles
  // with the time and pressure respectively.  Time to compare them.
  //
  // Loop from t0 to t1
  int n10 = 0;
  int n11 = 0;
  int n20 = 0;
  int n21 = 0;
  double integral = 0.0;
  double normal1 = 0.0;
  double normal2 = 0.0;
  double ti = t0;
  double infnorm = 0.0;

  // Find both probe values by linear interpolation at t0
  double pp10 = determine_probe_value(probe1,n10,n11,ti,0);
  double pp20 = determine_probe_value(probe2,n20,n21,ti,0); 
  double tf = ti;

  while(tf < t1){
    
    // Set tf to be the next time we have a known value on 
    // either curve.
    tf = ((probe1[n11].first < probe2[n21].first) &&
	  (probe1[n11].first > ti) ? probe1[n11].first :
	  probe2[n21].first);

    // Big problem if tf = ti Note: this can actually be okay in 
    // exceptional cases, but we are dealing with those in the 
    // determine_probe_value routine so that when we get to this
    // point, if tf = ti, we have to fail.  
    if(tf == ti){
      std::cerr << program_name << ": Encountered 0 interval. "
		<< "Cannot proceed." << std::endl
		<< "Here's some debugging information:" << std::endl
		<< "Getting probe 1 data. (" << n10 << "," << n11 
		<< "," << tf << "," << probe1[n10].first << ")" 
		<< std::flush << std::endl;
      double pp11 = determine_probe_value(probe1,n10,n11,tf,n10);
      std::cerr << "Got probe 1 data. (" << n10 << "," << n11 
		<< "," << tf << "," << probe1[n10].first << ")" 
		<< std::flush << std::endl
		<< "Getting probe 2 data. (" << n20 << "," << n21 
		<< "," << tf << "," << probe2[n20].first << ")" 
		<< std::flush << std::endl;
      double pp21 = determine_probe_value(probe2,n20,n21,tf,n20);
      std::cerr << "Got probe 2 data. (" << n20 << "," << n21 
		<< "," << tf << "," << probe2[n20].first << ")" 
		<< std::flush << std::endl
		<< "probe1(" << n10 << "," << n11 << ")["
		<< probe1[n11].first << "]" << std::endl
		<< "probe2(" << n20 << "," << n21 << ")["
		<< probe2[n21].first << "]" << std::endl;
      exit(1);
    }

    // Reset the tf if it exceeds the maximum time in our interval
    if(tf > t1) tf = t1;

    // Determine the probe values at tf by linear interpolation
    double pp11 = determine_probe_value(probe1,n10,n11,tf,n10);
    double pp21 = determine_probe_value(probe2,n20,n21,tf,n20);

    // Evaluate the integrated functions exactly as integrated lines
    // on this dt interval
    double dt = tf - ti;
    double dp1 = pp10 - pp11;
    double dp2 = pp20 - pp21;
    double A1 = dp1/dt;
    double A2 = dp2/dt;
    double p12 = pp10 - pp20;
    double A12 = A1 - A2;

    // Integral(f - g) for infinity norm
    double infterm = (p12*dt + A12*dt*dt/2.0);
    if(std::fabs(infterm) > std::fabs(infnorm))
      infnorm = infterm;
    
    // Integral(f) and Integral(g) for normalization
    normal1 += (pp10*dt + A1*dt*dt/2.0);
    normal2 += (pp20*dt + A2*dt*dt/2.0);

    // Integral((f-g)^2) for L2
    integral += (p12*p12*dt + p12*A12*dt*dt + A12*A12*dt*dt*dt/3.0);

    // Update our place along the curves before continuing
    ti = tf;
    pp10 = pp11;
    pp20 = pp21;
  }

  // L2 norm integral
  integral = std::sqrt(integral);

  // Normalization (not sure if this is the way to go)
  double normal = std::sqrt(normal1*normal1+normal2*normal2);

  std::cout.precision(12);
  std::cout << program_name << ": Error Integral = "  
	    << std::scientific << integral << std::endl
	    << program_name << ": Probe 1 Area = " 
	    << std::scientific << normal1 << std::endl
	    << program_name << ": Probe 2 Area = " 
	    << std::scientific << normal2 << std::endl
	    << program_name << ": Max Error = " 
	    << std::scientific << infnorm << std::endl
	    << program_name << ": Normalized Errors:" << std::endl
	    << program_name << ": L2 = " 
	    << std::scientific << integral/normal << std::endl
	    << program_name << ": LINF = " 
	    << std::scientific << infnorm/normal << std::endl;
  exit(0);
}






