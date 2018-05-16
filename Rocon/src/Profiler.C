///
/// @file
/// @ingroup irad_group
/// @brief Performance Profiling implementation
///
#include <cmath>
#include <iomanip>

#include "Profiler.H"
#include "primitive_utilities.H"

namespace IRAD {
  namespace Profiler {

    ///
    /// faster timeval math 
    ///
    inline struct timeval 
    operator-(const struct timeval &t1, const struct timeval &t2)
    {
      struct timeval rtv;
      rtv.tv_usec = t1.tv_usec - t2.tv_usec;
      rtv.tv_sec  = t1.tv_sec - t2.tv_sec;
      return(rtv);
    }
  
    ///
    /// faster timeval math 
    ///
    inline struct timeval &
    operator-=(struct timeval &t1,const struct timeval &t2)
    {
      t1.tv_usec -= t2.tv_usec;
      t1.tv_sec -= t2.tv_sec;
      return(t1);
    }
  
    ///
    /// faster timeval math 
    ///
    inline struct timeval
    operator+(const struct timeval &t1,const struct timeval &t2)
    {
      struct timeval rtv;
      rtv.tv_usec = t1.tv_usec + t2.tv_usec;
      rtv.tv_sec = t1.tv_sec + t2.tv_sec;
      return(rtv);
    }
  
    ///
    /// faster timeval math 
    ///
    inline struct timeval &
    operator+=(struct timeval &t1,const struct timeval &t2)
    {
      t1.tv_usec += t2.tv_usec;
      t1.tv_sec  += t2.tv_sec;
      return(t1);
    }
  
    ///
    /// faster timeval math 
    ///
    inline bool
    operator==(const struct timeval &t1,const struct timeval &t2)
    {
      return((t1.tv_sec + (t1.tv_usec/1000000.0)) == 
	     (t2.tv_sec + (t2.tv_usec/1000000.0)));
    }
  
    ///
    /// faster timeval math 
    ///
    inline bool
    operator!=(const struct timeval &t1,const struct timeval &t2)
    {
      return(!(t1==t2));
    }
  
    ///
    /// faster timeval math 
    ///
    inline bool 
    operator<(const timeval &t1,const timeval &t2)
    {
      return((t1.tv_sec + (t1.tv_usec/1000000.0)) < 
	     (t2.tv_sec + (t2.tv_usec/1000000.0)));
    }
  
    ///
    /// faster timeval math 
    ///
    inline bool
    operator>(const struct timeval &t1,const struct timeval &t2)
    {
      return((t1.tv_sec + (t1.tv_usec/1000000.0)) > 
	     (t2.tv_sec + (t2.tv_usec/1000000.0)));
    }
  
    ///
    /// faster timeval math 
    ///
    inline bool
    operator<=(const struct timeval &t1,const struct timeval &t2)
    {
      return((t1 < t2) || (t1 == t2));
    }
  
    ///
    /// faster timeval math 
    ///
    inline bool
    operator>=(const struct timeval &t1,const struct timeval &t2)
    {
      return((t1 > t2) || (t1 == t2));
    }
  

    ///
    /// Event math
    ///
    bool
    operator<(const std::pair<unsigned int,std::list<Event> > &p1,
	      const std::pair<unsigned int,std::list<Event> > &p2)
    {
      return (p1.first < p2.first);
    };

    /// Output function for Event
    std::ostream &
    operator<<(std::ostream &ost,const Event &e)
    {
      ost << e._id << " " << e._timestamp << " " 
	  << e._inclusive << " " << e._exclusive;
      return(ost);
    }

    /// Input function for Event
    std::istream &
    operator>>(std::istream &ist,Event &e)
    {
      ist >> e._id >> e._timestamp >> e._inclusive >> e._exclusive;
      return(ist);
    }

    // TYPE DEFINITIONS
    typedef std::map<unsigned int,cumulative_stats> StatMap;
    typedef std::list<std::pair<unsigned int,StatMap> > PStatList;
    typedef std::map<unsigned int,parallel_stats> PStatMap;
    typedef std::map<unsigned int,PStatMap> ScalaMap;
    typedef std::list<std::pair<unsigned int,std::list<Event> > > PEventList;
    typedef std::map<unsigned int,scalability_stats> ScalaStatMap;
  
    ProfilerObj::ProfilerObj(){
      profiler_rank = 0;
      verblevel = 0;
      Out = NULL;
      Err = NULL;
      //    function_map["Application"] = 0;
      //    configmap[0] = "Application";
      nfunc = 0;
      _initd = false;
      _finalized = false;
    };
  
    int ProfilerObj::Init(int id){
      if(_initd){
	std::cerr << "ProfilerObj::Init: Error: already initialized." << std::endl;
	return(1);
      }
      profiler_rank = (unsigned int)id;
      Event begin(0);
      time0 = Time();
      begin.timestamp(time0);
      open_event_list.push_back(begin);
#ifdef WITH_HPM_TOOLKIT
      hpmInit((int)id,(char *)configmap[0].c_str());
#endif
#ifdef WITH_PAPI
      // Init PAPI Library
      //  if(!PAPIX_library_init())
      //    exit(1);
      // Find max number of counters and set num_hwc
      // Read config file or environment variables for possible user settings
      // Default for any remaining unset settings
      // Populate EventSet
#endif
      if(function_map.empty()){
	function_map["Application"] = 0;
	configmap[0] = "Application";
      }
      _initd = true;
      return(0);
    };
  
    int ProfilerObj::Init(const std::string &name,int id){
      if(_initd){
	std::cerr << "ProfilerObj::Init: Error: already initialized. Tried to reinit with " 
		  << name << "." << std::endl;
	return(1);
      }      //    ProfileName = name.c_str();
      //    function_map.insert(make_pair(name,0));
      function_map[name] = 0;
      configmap[0] = name;
      //	event_file.assign(name);
      return(Init(id));
    };
  
    int ProfilerObj::FunctionEntry(const std::string &name){
      unsigned int id = function_map[name];
      if(id == 0){
	id = ++nfunc;
	function_map[name] = id;
	configmap[id] = name;
      }
      Event e(id);
      double t = Time() - time0;
      //  assert(t > 0);
      e.timestamp(t);
      e.exclusive(0.0);
#ifdef WITH_HPM_TOOLKIT
      hpmStart(id,(char *)name.c_str());
#endif
#ifdef WITH_PAPI
      // Read HWC into e.hwc_vals,e.hwc_tree to 0
#endif
      open_event_list.push_front(e);
      return(0);
    };
    int ProfilerObj::FunctionEntry(int id){
      assert(!((unsigned int)id <= 0));
      Event e((unsigned int)id);
      double t = Time() - time0;
      e.timestamp(t);
      e.exclusive(0.0);
#ifdef WITH_HPM_TOOLKIT
      std::string name = configmap[id];
      if(name.empty()){
	std::ostringstream Ostr;
	Ostr << "Function" << id;
	name = Ostr.str();
	function_map[name] = id;
	configmap[id] = name;
      }
      hpmStart(id,(char *)name.c_str());
#endif
#ifdef WITH_PAPI
      // Read HWC into e.hwc_vals,e.hwc_tree to 0
#endif
      open_event_list.push_front(e);
      return(0);
    };
  
    int ProfilerObj::FunctionExit(const std::string &name){
      unsigned int id = function_map[name];
      std::list<Event>::iterator ei = open_event_list.begin();
      // This means unmatched function name
      if(id == 0)
	std::cerr << "Mismatched(" << profiler_rank
		  << "):" << name 
		  << ", expected " << configmap[ei->id()]
		  << std::endl;
      assert(id != 0);
      if(id != ei->id())
	std::cerr << "Mismatched(" << profiler_rank 
		  << "):" << name << ", expected "
		  << configmap[ei->id()] << std::endl;
      assert(id == ei->id());
      double t = Time() - time0;
      double inclusive = t - ei->timestamp();
      double exclusive = inclusive - ei->exclusive();
      //  assert(inclusive > 0 && exclusive > 0);
      ei->inclusive(inclusive);
      ei->exclusive(exclusive);
#ifdef WITH_HPM_TOOLKIT
      hpmStop(id);
#endif
#ifdef WITH_PAPI
      // Read HWC
      // ei->hwc_vals = hwc_vals - ei->hwc_vals;
      // ei->hwc_tree = ei->hwc_vals - ei->hwc_tree;
      // hwc_vals = ei->hwc_vals;
#endif
      event_list.push_back(*ei);
      ei++;
      t = ei->exclusive();
      ei->exclusive(t + inclusive);
#ifdef WITH_PAPI
      // ei->hwc_tree += hwc_vals;
#endif
      open_event_list.pop_front();
      return(0);
    };
    int ProfilerObj::FunctionExit(int id){
#ifdef WITH_HPM_TOOLKIT
      hpmStop((int)id);
#endif
      std::list<Event>::iterator ei = open_event_list.begin();
      assert((unsigned int)id == ei->id());
      double t = Time() - time0;
      double inclusive = t - ei->timestamp();
      double exclusive = inclusive - ei->exclusive();
      //  assert(inclusive > 0 && exclusive > 0);
      ei->inclusive(inclusive);
      ei->exclusive(exclusive);
#ifdef WITH_PAPI
      // Read HWC
      // ei->hwc_vals = hwc_vals - ei->hwc_vals;
      // ei->hwc_tree = ei->hwc_vals - ei->hwc_tree;
      // hwc_vals = ei->hwc_vals;
#endif
      event_list.push_back(*ei);
      ei++;
      t = ei->exclusive();
      ei->exclusive(t + inclusive);
#ifdef WITH_PAPI
      // ei->hwc_tree += hwc_vals;
#endif
      open_event_list.pop_front();
      return(0);
    };
    
    /// Close all preparing for some emergency exit probably.
    int ProfilerObj::FunctionExitAll(){
      std::list<Event>::iterator ei = open_event_list.begin();
      while(ei != open_event_list.end()){
        unsigned int id = ei->id();        
#ifdef WITH_HPM_TOOLKIT
        hpmStop((int)id);
#endif
        double t = Time() - time0;
        double inclusive = t - ei->timestamp();
        double exclusive = inclusive - ei->exclusive();
        //  assert(inclusive > 0 && exclusive > 0);
        ei->inclusive(inclusive);
        ei->exclusive(exclusive);
#ifdef WITH_PAPI
        // Read HWC
        // ei->hwc_vals = hwc_vals - ei->hwc_vals;
        // ei->hwc_tree = ei->hwc_vals - ei->hwc_tree;
        // hwc_vals = ei->hwc_vals;
#endif
        event_list.push_back(*ei);
        ei++;
        if(ei != open_event_list.end()){
          t = ei->exclusive();
          ei->exclusive(t + inclusive);
        }
        open_event_list.pop_front();
      }
#ifdef WITH_PAPI
      // ei->hwc_tree += hwc_vals;
#endif
      return(0);
    };

    /// dump closed events to Ostr, clear mem
    int ProfilerObj::Dump(std::ostream &Ostr){return(0);};
    int ProfilerObj::ReadConfig(const std::string &cfname)
    {
      if(cfname.empty()){
	if(Err)
	  *Err << "ProfilerObj::Error: Config file name empty."
	       << std::endl;
	return(1);
      }
      std::ifstream conf_file;
      unsigned int cid;
      std::string routine;
      conf_file.open(cfname.c_str());
      if(!conf_file){
	if(Err)
	  *Err << "ProfilerObj::Error: Could not open config file, "
	       << cfname << "." << std::endl;
	return(1);
      }
      std::string configline;
      while(std::getline(conf_file,configline)){
	std::istringstream Istr(configline);
	Istr >> cid;
	std::getline(Istr,routine);
	configmap[cid] = routine;
	function_map[routine] = cid;
      }
      conf_file.close();
      return(0);
    }
    void ProfilerObj::SummarizeSerialExecution(std::ostream &Ostr)
    {  
      std::string application_name = "Application";
      std::map<unsigned int,cumulative_stats> statmap;
      std::list<Event>::iterator ei = event_list.begin();
      ei++;
      while(ei != event_list.end()){
	std::map<unsigned int,cumulative_stats>::iterator si;
	si = statmap.find(ei->id());
	if(si == statmap.end()){
	  cumulative_stats cs;
	  cs.incl = ei->inclusive();
	  cs.excl = ei->exclusive();
	  cs.ncalls = 1;
	  cs.incl_dev = cs.incl * cs.incl;
	  cs.excl_dev = cs.excl * cs.excl;
	  statmap.insert(std::make_pair(ei->id(),cs));
	}
	else{
	  si->second.incl += ei->inclusive();
	  si->second.excl += ei->exclusive();
	  si->second.ncalls++;
	  si->second.incl_dev += (ei->inclusive() * ei->inclusive());
	  si->second.excl_dev += (ei->exclusive() * ei->exclusive());
	}
	ei++;
      }
      std::map<unsigned int,cumulative_stats>::iterator si = statmap.begin();
      std::map<unsigned int,std::string>::iterator cmi = configmap.find(0);
      if(cmi != configmap.end())
	application_name = cmi->second;
      ei = event_list.begin();
      Ostr << "#Statistics for " << application_name << ":" << std::endl
	   << std::endl
	   << "#Total Execution Time: " << ei->inclusive() << std::endl
	   << "#------------------------------------------" 
	   << "Breakdown by Routine" 
	   << "------------------------------------------" << std::endl
	   << "#Routine Name                     Inclusive     Exclusive"
	   << "  Count          I-Mean" 
	   << "       I-Std        E-Mean        E-Std   " << std::endl
	   << "#------------------------------- ------------ ------------"
	   << " -----       ------------"
	   << " ------------ ------------ ------------" << std::endl;
      while(si != statmap.end()){
	std::string routine_name = "Unknown";
	cmi = configmap.find(si->first);
	Ostr << std::setiosflags(std::ios::left);
	if(cmi != configmap.end())
	  routine_name = cmi->second;
	else{
	  std::ostringstream Ostr;
	  Ostr << routine_name << " (" << si->first << ")";
	  routine_name = Ostr.str();
	}
	double imean = si->second.incl/(double)si->second.ncalls;
	double emean = si->second.excl/(double)si->second.ncalls;
	double imean2 = si->second.incl_dev/(double)si->second.ncalls;
	double emean2 = si->second.excl_dev/(double)si->second.ncalls;
	Ostr << std::setw(32) << routine_name << " ";
	Ostr << std::setprecision(5)
	     << std::setw(12) << si->second.incl << " " 
	     << std::setw(12) << si->second.excl << " " 
	     << std::setw(12) << si->second.ncalls << " "
	     << std::setw(12) << imean << " " 
	     << std::setw(12) 
	     << (si->second.ncalls == 1 ? 0 : 
		 sqrt(std::fabs(imean2-(imean*imean)))) 
	     << " " 
	     << std::setw(12) << emean << " " 
	     << std::setw(12) 
	     << (si->second.ncalls == 1 ? 0 : 
		 sqrt(std::fabs(emean2-(emean*emean)))) 
	     << std::endl;
	si++;
      }
    };
    void ProfilerObj::DumpEvents(std::ostream &Ostr)
    {
      Ostr << profiler_rank << std::endl;
      Util::DumpContents(Ostr,event_list);
      Ostr << std::endl;
    }
    void ProfilerObj::WriteEventFile()
    {
      std::ofstream eventfile;
      std::ostringstream Ostr;
      Ostr << configmap[0] << ".prof_";
      if(!(profiler_rank/10000))
	Ostr << "0";
      if(!(profiler_rank/1000))
	Ostr << "0";
      if(!(profiler_rank/100))
	Ostr << "0";
      if(!(profiler_rank/10))
	Ostr << "0";
      Ostr << profiler_rank;
      eventfile.open(Ostr.str().c_str());
      DumpEvents(eventfile);
      eventfile.close();
    };
    int ProfilerObj::Finalize(){
      if(_finalized)
	return(0);
      double t = Time();
      // This means there are unclosed events
      assert(open_event_list.size() == 1);
      std::list<Event>::iterator ei = open_event_list.begin();
      ei->inclusive(t - ei->timestamp());
      ei->exclusive(ei->inclusive()-ei->exclusive());
      ei->timestamp(0.);
#ifdef WITH_PAPI
      // Read the HWC
      // Update the overall counter (begin - from init):
      //   - begin.hwc_vals = hwc_vals - begin.hwc_vals;
      //   - begin..hwc_tree = begin.hwc_vals - begin.hwc_tree;
#endif
      event_list.push_front(*ei);
      event_list.sort();
      if(profiler_rank == 0){
	if(!configmap[0].empty()){
	  std::ofstream configfile;
	  std::ostringstream Bfn;
	  Bfn << configmap[0] << ".rpconfig";
	  configfile.open(Bfn.str().c_str());
	  FunctionMap::iterator fmi = function_map.begin();
	  while(fmi != function_map.end()){
	    configfile << fmi->second << " " << fmi->first << std::endl;
	    fmi++;
	  }
	  configfile.close();
	}
      } 
      WriteEventFile();
      //      if(summary && profiler_rank==0)
      //	summarize_execution();
#ifdef WITH_HPM_TOOLKIT
      hpmTerminate(profiler_rank);
#endif
#ifdef WITH_PAPI
      // Shutdown PAPI lib?  Maybe this is not needed.
#endif
      _finalized = true;
      return(0);
    };

    int ProfilerObj::ReadEventsFromFile(const std::string &filename){
      std::ifstream datafile;
      datafile.open(filename.c_str());
      if(!datafile){
	if(Err)
	  *Err << "ProfilerObj::ReadEventsFromFile: Error: Could not"
	       << " open datafile, " << filename << "." << std::endl;
	return(1);
      }
      Profiler::Event e;
      datafile >> profiler_rank;
      while(datafile >> e)
	event_list.push_back(e);
      datafile.close();
      return(0);
    }
     
    int 
    ProfilerObj::ReadParallelEventFiles(const std::vector<std::string> &ifiles,
					Profiler::PEventList &par_event_list)
    {
      if(ifiles.empty()){
	if(Err)
	  *Err << "ProfilerObj::ReadParallelEventFiles:Error: No input files."
	       << std::endl;
	return(1);
      }
      //    unsigned int number_of_processors = ifiles.size();
      std::vector<std::string>::const_iterator ifi = ifiles.begin();
      while(ifi != ifiles.end()){
	std::ifstream Inf;
	Inf.open(ifi->c_str());
	if(!Inf){
	  if(Err)
	    *Err << "ProfilerObj::ReadParallelEventFiles:Error: Unable to open" 
		 << " event file, " << *ifi << "." << std::endl;
	  return(1);
	}
	event_list.clear();
	Inf >> profiler_rank;
	Event e;
	while(Inf >> e)
	  event_list.push_back(e);
	Inf.close();
	event_list.sort();
	par_event_list.push_back(std::make_pair(profiler_rank,event_list));
	ifi++;
      }
      par_event_list.sort();
      return(0);
    }


    int ProfilerObj::SummarizeParallelExecution(std::ostream &Ostr,
						std::ostream &Ouf,
						PEventList &parallel_event_list)
    {
      if(parallel_event_list.empty())
	return(1);
      std::string application_name = "Application";
      std::map<unsigned int,std::string>::iterator cmi = configmap.find(0);
      if(cmi != configmap.end())
	application_name = cmi->second;
      PEventList::reverse_iterator peri = parallel_event_list.rbegin();
      // Assuming ranks of 0 to nproc-1
      unsigned int number_of_processors = peri->first + 1;
      Ouf << number_of_processors << std::endl;
      PStatMap pstat_map;
      PStatList parallel_cstat_list;
      std::map<unsigned int,cumulative_stats> statmap;
      PEventList::iterator peli = parallel_event_list.begin();
      while(peli != parallel_event_list.end()){
	profiler_rank = peli->first;
	std::list<Event>::const_iterator eli = peli->second.begin();
	while(eli != peli->second.end()){
	  std::map<unsigned int,cumulative_stats>::iterator si;
	  si = statmap.find(eli->id());
	  if(si == statmap.end()){
	    cumulative_stats cs;
	    cs.incl = eli->inclusive();
	    cs.excl = eli->exclusive();
	    cs.ncalls = 1;
	    cs.incl_dev = cs.incl * cs.incl;
	    cs.excl_dev = cs.excl * cs.excl;
	    statmap.insert(std::make_pair(eli->id(),cs));
	  }
	  else{
	    if(eli->inclusive() < 0){
	      if(Err)
		*Err << "ProfilerObj::SummarizeParallelExecution:Error: "
		     << "Read existing negative inclusive value:" 
		     << eli->inclusive() << std::endl;
	      return(1);
	    }
	    if(eli->exclusive() < 0){
	      if(Err)
		*Err << "ProfilerObj::SummarizeParallelExecution:Error: "
		     << "Read existing negative exclusive value:" 
		     << eli->inclusive() << std::endl;
	      return(1);
	    }
	    si->second.incl += eli->inclusive();
	    si->second.excl += eli->exclusive();
	    si->second.ncalls++;
	    si->second.incl_dev += (eli->inclusive() * eli->inclusive());
	    si->second.excl_dev += (eli->exclusive() * eli->exclusive());
	  }
	  eli++;
	}
	parallel_cstat_list.push_back(std::make_pair(profiler_rank,statmap));
	statmap.clear();
	peli++;
      }
      // Loop over the cumulative stats on each processor, and build the map 
      // containing the Min, Max, Mean, StdDev, of the cumulative times and 
      // the number of calls for each routine.
      PStatList::iterator psli = parallel_cstat_list.begin();
      // Loop over processors
      while(psli != parallel_cstat_list.end()){
	profiler_rank = psli->first;
	// Loop over routines
	std::map<unsigned int,cumulative_stats>::iterator si = 
	  psli->second.begin();
	while(si != psli->second.end()){
	  std::map<unsigned int,parallel_stats>::iterator psmi;
	  // Find the entry in the pstat_map that is for this routine
	  psmi = pstat_map.find(si->first);
	  // Process the very first entry for this routine.
	  if(psmi == pstat_map.end()){
	    parallel_stats ps;
	    ps.incl_min = si->second.incl;
	    ps.incl_max = si->second.incl;
	    ps.incl_mean = si->second.incl;
	    ps.incl_stdev = si->second.incl * si->second.incl;
	    ps.incl_minrank = profiler_rank;
	    ps.incl_maxrank = profiler_rank;
	    ps.excl_min = si->second.excl;
	    ps.excl_max = si->second.excl;
	    ps.excl_mean = si->second.excl;
	    ps.excl_stdev = si->second.excl * si->second.excl;
	    ps.excl_minrank = profiler_rank;
	    ps.excl_maxrank = profiler_rank;
	    ps.call_max = si->second.ncalls;
	    ps.call_min = si->second.ncalls;
	    ps.call_mean = si->second.ncalls;
	    ps.call_stdev = si->second.ncalls * si->second.ncalls;
	    ps.call_minrank = profiler_rank;
	    ps.call_maxrank = profiler_rank;
	    pstat_map.insert(std::make_pair(si->first,ps));
	  }
	  // Accumulate data for existing entry
	  else {
	    if(psmi->second.incl_min > si->second.incl){
	      psmi->second.incl_min = si->second.incl;
	      psmi->second.incl_minrank = profiler_rank;
	    }
	    if(psmi->second.incl_max < si->second.incl){
	      psmi->second.incl_max = si->second.incl;
	      psmi->second.incl_maxrank = profiler_rank;
	    }
	    psmi->second.incl_mean += si->second.incl;
	    psmi->second.incl_stdev += si->second.incl * si->second.incl;
	  
	    if(psmi->second.excl_min > si->second.excl){
	      psmi->second.excl_min = si->second.excl;
	      psmi->second.excl_minrank = profiler_rank;
	    }
	    if(psmi->second.excl_max < si->second.excl){
	      psmi->second.excl_max = si->second.excl;
	      psmi->second.excl_maxrank = profiler_rank;
	    }
	    psmi->second.excl_mean += si->second.excl;
	    psmi->second.excl_stdev += si->second.excl * si->second.excl;
	  
	    if(psmi->second.call_max < si->second.ncalls){
	      psmi->second.call_max = si->second.ncalls;
	      psmi->second.call_maxrank = profiler_rank;
	    }
	    if(psmi->second.call_min > si->second.ncalls){
	      psmi->second.call_min = si->second.ncalls;
	      psmi->second.call_minrank = profiler_rank;
	    }
	    psmi->second.call_mean += si->second.ncalls;
	    psmi->second.call_stdev += si->second.ncalls * si->second.ncalls;
	  }
	  si++;
	}
	psli++;
      }      
      std::map<unsigned int,parallel_stats>::iterator si = pstat_map.begin();
      Ostr << "#Statistics for " << application_name << " ("
	   << number_of_processors << " procs):" << std::endl << std::endl
	   << "#----------------------------------Inclusive Statistics" 
	   << "------------------------------" << std::endl
	   << "#                                   Min Inc     Min    "
	   << "Max Inc     Max" 
	   << "    Mean Inc               " << std::endl
	   << "#Routine Name                      Duration    Rank   "
	   << "Duration    Rank" 
	   << "   Duration     Std Dev   " << std::endl
	   << "#--------------------            ------------ ----- "
	   << "------------ -----"
	   << " ------------ ------------" << std::endl;
      while(si != pstat_map.end()){
	std::string routine_name = "Unknown";
	cmi = configmap.find(si->first);
	Ostr << std::setiosflags(std::ios::left);
	if(cmi != configmap.end())
	  routine_name = cmi->second;
	else{
	  std::ostringstream Ostr2;
	  Ostr2 << routine_name << " (" << si->first << ")";
	  routine_name = Ostr2.str();
	}
	double imean = si->second.incl_mean/(double)number_of_processors;
	double imean2 = (si->second.incl_min == si->second.incl_max ? 0 :
			 (sqrt(fabs((si->second.incl_stdev/
				     (double)number_of_processors) 
				    - (imean * imean)))));
	Ostr << std::setw(32) << routine_name << " "
	     << std::setw(12) << si->second.incl_min << " "
	     << std::setw(5)  << si->second.incl_minrank << " "
	     << std::setw(12) << si->second.incl_max << " " 
	     << std::setw(5)  << si->second.incl_maxrank << " "
	     << std::setw(12) << imean  << " " 
	     << std::setw(12) << imean2 << " "
	     << std::endl;
	Ouf << si->first << " " << si->second.incl_min << " "
	    << si->second.incl_minrank << " " 
	    << si->second.incl_max << " " << si->second.incl_maxrank
	    << " " << imean << " " << imean2 << std::endl;
	si++;
      }
      si = pstat_map.begin();
      Ostr << "#----------------------------------Exclusive Statistics"
	   << "------------------------------" << std::endl
	   << "#                                   Min Exc     Min    "
	   << "Max Exc     Max"
	   << "    Mean Exc               " << std::endl
	   << "#Routine Name                      Duration    Rank   "
	   << "Duration    Rank"
	   << "   Duration     Std Dev   " << std::endl
	   << "#--------------------            ------------ ----- "
	   << "------------ -----"
	   << " ------------ ------------" << std::endl;
      while(si != pstat_map.end()){
	std::string routine_name = "Unknown";
	cmi = configmap.find(si->first);
	Ostr << std::setiosflags(std::ios::left);
	if(cmi != configmap.end())
	  routine_name = cmi->second;
	else{
	  std::ostringstream Ostr2;
	  Ostr2 << routine_name << " (" << si->first << ")";
	  routine_name = Ostr2.str();
	}
	double emean = si->second.excl_mean/(double)number_of_processors;
	double emean2 = (si->second.excl_max == si->second.excl_min ? 0 :
			 (sqrt(fabs((si->second.excl_stdev/
				     (double)number_of_processors)
				    - (emean * emean)))));
	Ostr << std::setw(32) << routine_name << " "
	     << std::setw(12) << si->second.excl_min << " "
	     << std::setw(5)  << si->second.excl_minrank << " "
	     << std::setw(12) << si->second.excl_max << " " 
	     << std::setw(5)  << si->second.excl_maxrank << " "
	     << std::setw(12) << emean  << " " 
	     << std::setw(12) << emean2 << " "
	     << std::endl;
	Ouf << si->first << " " << si->second.excl_min << " "
	    << si->second.excl_minrank << " " 
	    << si->second.excl_max << " " << si->second.excl_maxrank
	    << " " << emean << " " << emean2 << std::endl;
	si++;
      }
      return(0);
    }

    int
    ProfilerObj::ReadSummaryFiles(const std::vector<std::string> &input_files,
				  ScalaMap &scala_map)
    {
    
      // I guess it's okay to have only a single summary file - it's only a 
      // problem if it's intended use is a scalability summary
      int number_of_runs = input_files.size();
      if(number_of_runs < 1){
	if(Err)
	  *Err << "ProfilerObj::ReadSummaryFiles:Error: No input files."
	       << std::endl;
	return(1);
      }
      std::vector<int> problem_sizes;
      problem_sizes.resize(number_of_runs);
      number_of_runs = 0;
      unsigned int runsize = 0;
      std::vector<std::string>::const_iterator ifi = input_files.begin();
      while(ifi != input_files.end()){
	PStatMap pstats;
	std::ifstream Inf;
	Inf.open(ifi->c_str());
	if(!Inf){
	  if(Err)
	    *Err << "ProfilerObj::ReadSummaryFiles:Error: Cannot open summary "
		 << "file, " << *ifi << "." << std::endl;
	  return(1);
	}
	unsigned int id, minrank, maxrank;
	double min, max, mean, stddev;
	Inf >> runsize;
	while(Inf >> id >> min >> minrank >> max >> maxrank >> mean >> stddev){
	  // if this function's id cannot be found in (PStatMap)pstats, then it's
	  // the inclusive section, populate a new parallel_stats object and add
	  // it and the id to the pstats object.
	  PStatMap::iterator psi = pstats.find(id);
	  if(psi == pstats.end()){
	    parallel_stats pso;
	    pso.incl_min = min;
	    pso.incl_minrank = minrank;
	    pso.incl_max = max;
	    pso.incl_maxrank = maxrank;
	    pso.incl_mean = mean;
	    pso.incl_stdev = stddev;
	    pstats.insert(std::make_pair(id,pso));
	  }
	  // - else -
	  // it's the exclusive section, retrieve the existing parallel_stats 
	  // from the pstats and populate the exclusive items.
	  else {
	    psi->second.excl_min = min;
	    psi->second.excl_minrank = minrank;
	    psi->second.excl_max = max;
	    psi->second.excl_maxrank = maxrank;
	    psi->second.excl_mean = mean;
	    psi->second.excl_stdev = stddev;
	  }
	}
	Inf.close();
	// Add the (PStatMap)pstats and the number of processors (runsize) to the
	// ScalaMap object
	//
	// First, if a run of this size exists already, then fail as we have no
	// idea how to handle that (yet).
	ScalaMap::iterator smi = scala_map.find(runsize);
	if(smi != scala_map.end()){
	  if(Err)
	    *Err << "ProfilerObj::ReadSummaryFiles:Error: Cannot process "
		 << "multiple runs of the same size." << std::endl;
	  return(1);
	}
	problem_sizes[number_of_runs++] = runsize;
	//    vout << "Processed problem size: " << runsize << endl;
	scala_map.insert(std::make_pair(runsize,pstats));
	ifi++;
      }
      return(0);
    }
  
    enum {IMIN_T=0,IMAX_T,IMEAN_T,IMIN_E,IMAX_E,IMEAN_E,
	  IMIN_S,IMAX_S,IMEAN_S,EMIN_T,EMAX_T,EMEAN_T,
	  EMIN_E,EMAX_E,EMEAN_E,EMIN_S,EMAX_S,EMEAN_S};
  
    int ProfilerObj::PopulateScalaMap(ScalaMap &scala_map,
				      ScalaStatMap &scala_statmap,
				      bool is_scaled)
    {
    
      if(scala_map.empty() || (scala_map.size() == 1)){
	if(Err)
	  *Err << "ProfilerObj::Error: Invalid number of runs for scalability analysis." 
	       << std::endl;
	return(1);
      }
  
      // When we loop thru problem sizes, make sure it's ascending order
      //  sort(problem_sizes.begin(),problem_sizes.end());
      //  unsigned int number_of_runs = problem_sizes.size();
      //  tagit(DEBUGV,"populate_scalamap");
      //  vout << "Number of runs: " << number_of_runs << endl;
      // For every problem size
      ScalaMap::iterator scalamap_i = scala_map.begin();
      //  for(unsigned int ps = 0;ps < number_of_runs;ps++){
      while(scalamap_i != scala_map.end()){
	int number_of_processors = scalamap_i->first;
	//    tagit(DEBUGV,"populate_scalamap");
	//    vout << "Processing run: " << number_of_processors << endl;
	//      if(scalamap_i == scala_map.end()){
	//	tagit(1,"populate_scalamap");
	//	vout << "Error: No data for " << number_of_processors << " processors." 
	//	     << endl;
	//	return(false);
	//      }
	// For every Fid in the ScalaMap's PStat object
	PStatMap::iterator psm_i = scalamap_i->second.begin();
	while(psm_i != scalamap_i->second.end()){
	  // Look for this function's id in the scalastatmap
	  ScalaStatMap::iterator ssmi = scala_statmap.find(psm_i->first);
	  if(ssmi == scala_statmap.end()){  // Not found
	    // Make a new entry for this Fid
	    scalability_stats ss;
	    ss.nprocs.push_back(number_of_processors);
	    unsigned int cpos = ss.nprocs.size() - 1;
	    double probsize = (is_scaled ? number_of_processors : 1);
	    double probsize0 = (is_scaled ? ss.nprocs[0] : 1); 
	    ss.sstats[IMIN_T].push_back(psm_i->second.incl_min);
	    ss.sstats[IMAX_T].push_back(psm_i->second.incl_max);
	    ss.sstats[IMEAN_T].push_back(psm_i->second.incl_mean);
	    ss.sstats[IMIN_E].push_back(((ss.sstats[IMIN_T][0]*ss.nprocs[0])/
					 probsize0)/((ss.sstats[IMIN_T][cpos]*
						      ss.nprocs[cpos])/probsize));
	    ss.sstats[IMAX_E].push_back(((ss.sstats[IMAX_T][0]*ss.nprocs[0])/
					 probsize0)/((ss.sstats[IMAX_T][cpos]*
						      ss.nprocs[cpos])/probsize));
	    ss.sstats[IMEAN_E].push_back(((ss.sstats[IMEAN_T][0]*ss.nprocs[0])/
					  probsize0)/((ss.sstats[IMEAN_T][cpos]*
						       ss.nprocs[cpos])/probsize));
	    ss.sstats[IMIN_S].push_back(ss.sstats[IMIN_E][cpos]*
					number_of_processors);
	    ss.sstats[IMAX_S].push_back(ss.sstats[IMAX_E][cpos]*
					number_of_processors);
	    ss.sstats[IMEAN_S].push_back(ss.sstats[IMEAN_E][cpos]*
					 number_of_processors);
	  
	    ss.sstats[EMIN_T].push_back(psm_i->second.excl_min);
	    ss.sstats[EMAX_T].push_back(psm_i->second.excl_max);
	    ss.sstats[EMEAN_T].push_back(psm_i->second.excl_mean);
	    ss.sstats[EMIN_E].push_back(((ss.sstats[EMIN_T][0]*ss.nprocs[0])/
					 probsize0)/((ss.sstats[EMIN_T][cpos]*
						      ss.nprocs[cpos])/probsize));
	    ss.sstats[EMAX_E].push_back(((ss.sstats[EMAX_T][0]*ss.nprocs[0])/
					 probsize0)/((ss.sstats[EMAX_T][cpos]*
						      ss.nprocs[cpos])/probsize));
	    ss.sstats[EMEAN_E].push_back(((ss.sstats[EMEAN_T][0]*ss.nprocs[0])/
					  probsize0)/((ss.sstats[EMEAN_T][cpos]*
						       ss.nprocs[cpos])/probsize));
	    ss.sstats[EMIN_S].push_back(ss.sstats[EMIN_E][cpos]*
					number_of_processors);
	    ss.sstats[EMAX_S].push_back(ss.sstats[EMAX_E][cpos]*
					number_of_processors);
	    ss.sstats[EMEAN_S].push_back(ss.sstats[EMEAN_E][cpos]*
					 number_of_processors);
	    // Add an entry for this function into the scala_statmap
	    scala_statmap.insert(std::make_pair(psm_i->first,ss));
	  }
	  // There was an entry, this function had statistics at some other
	  // runsize.  Add to it's scalastats entries for this problem size.
	  else { 
	    // Make a new entry for this Fid
	    ssmi->second.nprocs.push_back(number_of_processors);
	    unsigned int cpos = ssmi->second.nprocs.size() - 1;
	    double probsize = (is_scaled ? number_of_processors : 1);
	    double probsize0 = (is_scaled ? ssmi->second.nprocs[0] : 1); 
	    ssmi->second.sstats[IMIN_T].push_back(psm_i->second.incl_min);
	    ssmi->second.sstats[IMAX_T].push_back(psm_i->second.incl_max);
	    ssmi->second.sstats[IMEAN_T].push_back(psm_i->second.incl_mean);
	    ssmi->second.sstats[IMIN_E].push_back
	      (((ssmi->second.sstats[IMIN_T][0]*ssmi->second.nprocs[0])/probsize0)/
	       ((ssmi->second.sstats[IMIN_T][cpos]*ssmi->second.nprocs[cpos])/
		probsize));
	    ssmi->second.sstats[IMAX_E].push_back
	      (((ssmi->second.sstats[IMAX_T][0]*ssmi->second.nprocs[0])/probsize0)/
	       ((ssmi->second.sstats[IMAX_T][cpos]*ssmi->second.nprocs[cpos])/
		probsize));
	    ssmi->second.sstats[IMEAN_E].push_back
	      (((ssmi->second.sstats[IMEAN_T][0]*ssmi->second.nprocs[0])/
		probsize0)/((ssmi->second.sstats[IMEAN_T][cpos]*
			     ssmi->second.nprocs[cpos])/probsize));
	    ssmi->second.sstats[IMIN_S].push_back
	      (ssmi->second.sstats[IMIN_E][cpos]*number_of_processors);
	    ssmi->second.sstats[IMAX_S].push_back
	      (ssmi->second.sstats[IMAX_E][cpos]*number_of_processors);
	    ssmi->second.sstats[IMEAN_S].push_back
	      (ssmi->second.sstats[IMEAN_E][cpos]*number_of_processors);
	  
	    ssmi->second.sstats[EMIN_T].push_back(psm_i->second.excl_min);
	    ssmi->second.sstats[EMAX_T].push_back(psm_i->second.excl_max);
	    ssmi->second.sstats[EMEAN_T].push_back(psm_i->second.excl_mean);
	    ssmi->second.sstats[EMIN_E].push_back
	      (((ssmi->second.sstats[EMIN_T][0]*ssmi->second.nprocs[0])/
		probsize0)/((ssmi->second.sstats[EMIN_T][cpos]*
			     ssmi->second.nprocs[cpos])/probsize));
	    ssmi->second.sstats[EMAX_E].push_back
	      (((ssmi->second.sstats[EMAX_T][0]*ssmi->second.nprocs[0])/
		probsize0)/((ssmi->second.sstats[EMAX_T][cpos]*
			     ssmi->second.nprocs[cpos])/probsize));
	    ssmi->second.sstats[EMEAN_E].push_back
	      (((ssmi->second.sstats[EMEAN_T][0]*ssmi->second.nprocs[0])/
		probsize0)/((ssmi->second.sstats[EMEAN_T][cpos]*
			     ssmi->second.nprocs[cpos])/probsize));
	    ssmi->second.sstats[EMIN_S].push_back
	      (ssmi->second.sstats[EMIN_E][cpos]*number_of_processors);
	    ssmi->second.sstats[EMAX_S].push_back
	      (ssmi->second.sstats[EMAX_E][cpos]*number_of_processors);
	    ssmi->second.sstats[EMEAN_S].push_back
	      (ssmi->second.sstats[EMEAN_E][cpos]*number_of_processors);
	  }
	  psm_i++;
	}
	scalamap_i++;
      }
      return(0);
    }

    int ProfilerObj::ScalabilitySummary(ScalaStatMap &scala_statmap,std::ostream &Out)
    {
      std::string appname = "Unknown";
      std::map<unsigned int,std::string>::iterator cfi = configmap.find(0);
      if(cfi != configmap.end())
	appname.assign(cfi->second);
      // For every routine in the scalamap, print out the scalability information
      ScalaStatMap::iterator ssm_i = scala_statmap.begin();
      Out << "############# Scalability Summary for " << appname << " "
	  << "#############" << std::endl; 
      while(ssm_i != scala_statmap.end()){
	std::string routine_name;
	std::map<unsigned int,std::string>::iterator cfi = configmap.find(ssm_i->first);
	if(cfi != configmap.end())
	  routine_name.assign(cfi->second);
	else {
	  std::ostringstream Ostr;
	  Ostr << "Unknown(" << ssm_i->first << ")";
	  routine_name.assign(Ostr.str());
	}
	scalability_stats ss = ssm_i->second;
	unsigned int nps = ss.nprocs.size();
	//      tagit(DEBUGV,"scalability_summary");
	//      vout << "Number of runs to process = " << nps << endl;
	Out << "# " << routine_name << ":" << std::endl
	    << "#-------------------------------------------------"
	    << "------------------------------------------------------------"
	    << "-----------------------------------------------" << std::endl
	    << "#          Inclusive Max            Inclusive Mean   "   
	    << "        Inclusive Min             Exclusive Max        "
	    << "   Exclusive Mean"
	    << "            Exclusive Min" << std::endl
	    << "# NProc  Time(Eff)(Speedup)       Time(Eff)(Speedup)"
	    << "       Time(Eff)(Speedup)       Time(Eff)(Speedup)       "
	    << "Time(Eff)(Speedup)       Time(Eff)(Speedup)" << std::endl
	    << "#-------------------------------------------------"
	    << "------------------------------------------------------------"
	    << "-----------------------------------------------" << std::endl;
	for(unsigned int a = 0;a < nps;a++){
	  Out << std::resetiosflags(std::ios::floatfield)
	      << std::setiosflags(std::ios::right) << std::setprecision(0)
	      << std::setw(5) << ss.nprocs[a]  << " " 
	  
	      << std::setprecision(3) << std::showpoint << std::fixed
	      << std::setw(9) << std::setiosflags(std::ios::right) << ss.sstats[IMAX_T][a] 
	      << " " << std::setiosflags(std::ios::left) << std::setprecision(2) 
	      << std::setw(4) << ss.sstats[IMAX_E][a] << " " 
	      << std::setw(7) << std::setprecision(1) << ss.sstats[IMAX_S][a] << "   "
	
	      << std::setprecision(3) << std::showpoint << std::fixed
	      << std::setw(9) << std::setiosflags(std::ios::right) << ss.sstats[IMEAN_T][a] 
	      << " " << std::setiosflags(std::ios::left) << std::setprecision(2) << std::setw(4) 
	      << ss.sstats[IMEAN_E][a] << " " 
	      << std::setw(7) << std::setprecision(1) << ss.sstats[IMEAN_S][a] << "   "
	
	      << std::setprecision(3) << std::showpoint << std::fixed
	      << std::setw(9) << std::setiosflags(std::ios::right) << ss.sstats[IMIN_T][a] 
	      << " " << std::setiosflags(std::ios::left) << std::setprecision(2) << std::setw(4) 
	      << ss.sstats[IMIN_E][a] << " " 
	      << std::setw(7) << std::setprecision(1) << ss.sstats[IMIN_S][a] << "   "
	
	      << std::setprecision(3) << std::showpoint << std::fixed
	      << std::setw(9) << std::setiosflags(std::ios::right) << ss.sstats[EMAX_T][a] 
	      << " " << std::setiosflags(std::ios::left) << std::setprecision(2) << std::setw(4) 
	      << ss.sstats[EMAX_E][a] << " " 
	      << std::setw(7) << std::setprecision(1) << ss.sstats[EMAX_S][a] << "   "
	  
	      << std::setprecision(3) << std::showpoint << std::fixed
	      << std::setw(9) << std::setiosflags(std::ios::right) << ss.sstats[EMEAN_T][a] 
	      << " " << std::setiosflags(std::ios::left) << std::setprecision(2) << std::setw(4) 
	      << ss.sstats[EMEAN_E][a] << " " 
	      << std::setw(7) << std::setprecision(1) << ss.sstats[EMEAN_S][a] << "   "
	
	      << std::setprecision(3) << std::showpoint << std::fixed
	      << std::setw(9) << std::setiosflags(std::ios::right) << ss.sstats[EMIN_T][a] 
	      << " " << std::setiosflags(std::ios::left) << std::setprecision(2) << std::setw(4) 
	      << ss.sstats[EMIN_E][a] << " " 
	      << std::setw(7) << std::setprecision(1) << ss.sstats[EMIN_S][a] << std::endl;
	}	
	Out << "#-------------------------------------------------"
	    << "------------------------------------------------------------"
	    << "-----------------------------------------------" << std::endl << std::endl 
	    << std::endl;
   
	ssm_i++;
      }
      return(0);
    }
  };
};
