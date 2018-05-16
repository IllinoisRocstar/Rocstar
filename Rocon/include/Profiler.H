///
///
/// @file
/// @ingroup irad_group
/// @brief Performance Profiling interface definition
///
#ifndef _PROFILER_H_
#define _PROFILER_H_
#include <map>
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <sys/time.h>

namespace IRAD {

  ///
  /// @brief Performance profiling
  /// @ingroup irad_group
  ///
  namespace Profiler {

    ///
    /// Utility struct
    ///
    struct cumulative_stats {
      double incl;
      double excl;
      unsigned int ncalls;
      double incl_dev;
      double excl_dev;
      cumulative_stats(){
	incl = excl = incl_dev = excl_dev = 0.0;
	ncalls = 0;
      };
    };
  

    ///
    /// Utility struct
    ///
    struct parallel_stats {
      double incl_min;
      double incl_max;
      unsigned int incl_minrank;
      unsigned int incl_maxrank;
      double incl_mean;
      double incl_stdev;
      double excl_min;
      double excl_max;
      unsigned int excl_minrank;
      unsigned int excl_maxrank;
      double excl_mean;
      double excl_stdev;
      double call_mean;
      unsigned int call_max;
      unsigned int call_min;
      unsigned int call_maxrank;
      unsigned int call_minrank;
      double call_stdev;
      parallel_stats(){
	incl_mean = incl_stdev = excl_min = excl_max = 0.0;
	excl_mean = excl_stdev = call_mean = call_stdev = 0.0;
	incl_minrank = incl_maxrank = excl_minrank = 0;
	excl_maxrank = call_max = call_min = call_maxrank = 0;
	call_minrank = 0;
	incl_min = incl_max = 0.0;
      };
    };
  
    ///
    /// Utility struct
    ///
    struct scalability_stats 
    {
      std::vector<unsigned int> nprocs;
      std::vector<double> sstats[18];
    };
  
    ///
    /// \brief Marks construct entry/exit
    ///
    /// Events mark the beginning, end, and duration of instances of user 
    /// defined code constructs.
    ///
    class Event {
      friend std::ostream &operator<<(std::ostream &,const Event &);
      friend std::istream &operator>>(std::istream &,Event &);
    protected:
      /// unique identifyer
      unsigned int _id;
      /// self time
      double _exclusive;
      /// tree time
      double _inclusive;
      /// raw timestamp
      double _timestamp;
#ifdef WITH_PAPI
      /// number of hw counters monitored
      int _nhwc;
      /// tree count
      long long *_inc_hwc;
      /// self count
      long long *_exc_hwc;
      bool _hwc;
#endif
    public:
      Event()
	: _id(0),_exclusive(0.),_inclusive(0.),_timestamp(0.)
      {
#ifdef WITH_PAPI
	init_papi();
#endif
      };
      Event(unsigned int i)
	: _id(i),_exclusive(0.),_inclusive(0.),_timestamp(0.)
      {
#ifdef WITH_PAPI
	init_papi();
#endif
      };
      Event(unsigned int i,double ts)
	: _id(i),_exclusive(0.),_inclusive(0.),_timestamp(0.)
      {
#ifdef WITH_PAPI
	init_papi();
#endif
      };
      Event(unsigned int i,double e,double it)
	: _id(i),_exclusive(e),_inclusive(it),_timestamp(0.)
      {
#ifdef WITH_PAPI
	init_papi();
#endif
      };
      Event(const Event &e)
	: _id(e._id),_exclusive(e._exclusive),_inclusive(e._inclusive),
	  _timestamp(e._timestamp)
      {
#ifdef WITH_PAPI
	_nhwc = e._nhwc;
	_inc_hwc = new long long [_nhwc];
	_exc_hwc = new long long [_nhwc];
	hwc = true;
#endif
      };
#ifdef WITH_PAPI
      ~Event()
      {
	if(_hwc){
	  delete [] _inc_hwc;
	  delete [] _exc_hwc;
	}
      };
#endif
      Event &
      operator=(const Event &e)
      {
	_id = e._id;
	_exclusive = e._exclusive;
	_inclusive = e._inclusive;
	_timestamp = e._timestamp;
#ifdef WITH_PAPI
	if(!_hwc){
	  _nwhc = e._nhwc;
	  _inc_hwc = new long long [_nhwc];
	  _exc_hwc = new long long [_nhwc];
	}
	for(int i = 0;i < _nwhc;i++){
	  _inc_hwc[i] = e._inc_hwc[i];
	  _exc_hwc[i] = e._exc_hwc[i];
	}
#endif
	return(*this);
      };
      double &exclusive()
      {
	return (_exclusive);
      };
      double exclusive() const
      {
	return (_exclusive);
      };
      double &inclusive()
      {
	return (_inclusive);
      };
      double inclusive() const 
      {
	return (_inclusive);
      };
      void exclusive(double e)
      {
	_exclusive = e;
      };
      void inclusive(double i)
      {
	_inclusive = i;
      };
      double &timestamp()
      {
	return (_timestamp);
      };
      double timestamp() const
      {
	return (_timestamp);
      };
      void timestamp(double t)
      {
	_timestamp = t;
      };
      unsigned int &id()
      {
	return(_id);
      };
      unsigned int id() const
      {
	return(_id);
      };
      void id(unsigned int i)
      {
	_id = i;
      };
      bool operator<(const Event &e) const 
      {
	return(_timestamp < e._timestamp);
      };
#ifdef WITH_PAPI
      void init_papi();
      int nhwc() comst { return (_nhwc); };
      void update_hwc(long long *,long long *);
      long long &inc_hwc(int i);
      long long inc_hwc(int i) const;
      long long &exc_hwc(int i);
      long long exc_hwc(int i) const;
#endif
    };
  
    ///
    /// EVENT's stream operator
    ///
    std::ostream &
    operator<<(std::ostream &ost,const Event &e);
  
    ///
    /// EVENT's stream operator
    ///
    std::istream &
    operator>>(std::istream &ist,Event &e);

    ///
    /// \brief Simple timer
    ///
    /// Time returns the time in seconds from the system time.  The user is 
    /// responsible for adding the appropriate offset to make this number 
    /// relative to the current execution.
    ///
    inline double
    Time()
    {
      struct timeval tv;
      gettimeofday(&tv,NULL);
      double t = tv.tv_sec + tv.tv_usec/1000000.;
      return(t);
    }
  
    ///
    /// construct name to unique id.
    ///
    typedef std::map<std::string,unsigned int> FunctionMap;
    ///
    /// construct name to unique id.
    ///
    typedef std::map<unsigned int,std::string> ConfigMap;
    typedef std::map<unsigned int,cumulative_stats> StatMap;
    typedef std::list<std::pair<unsigned int,StatMap> > PStatList;
    typedef std::map<unsigned int,parallel_stats> PStatMap;
    typedef std::map<unsigned int,PStatMap> ScalaMap;
    typedef std::list<std::pair<unsigned int,std::list<Event> > > PEventList;
    typedef std::map<std::string,unsigned int> FunctionMap;
    typedef std::map<unsigned int,scalability_stats> ScalaStatMap;

    ///
    /// noop profiler
    ///
    class DummyProfilerObj {

    public:
      int Init(int id){return 0;};
      int Init(const std::string &name,int id){return 0;};
      int FunctionEntry(const std::string &name){return 0;};
      int FunctionEntry(int id){return 0;};
      int FunctionExit(const std::string &name){return 0;};
      int FunctionExit(int id){return 0;};
      int FunctionExitAll(){return 0;};
      int Finalize(){return 0;};
      int Dump(std::ostream &Ostr){return 0;}; 
      bool Ready() {return true;};
    };
    
    /// \brief Performance profiling object
    ///
    /// Profiler is a profiling object that keeps track of call durations 
    /// for user defined code constructs.  Profiler will also monitor platform 
    /// hardware counters if they are available in PAPI or HPM.  The PAPI 
    /// interface is currently under construction to use the version independent 
    /// implementation
    /// 
    class ProfilerObj {
    protected:
      /// parallel processor id
      unsigned int profiler_rank;
      /// verbosity level
      unsigned int verblevel;
      /// stream for regular output
      std::ostream *Out;
      /// stream for errors
      std::ostream *Err;
      /// creation/init time
      double time0;
      /// construct entry events
      std::list<Event> open_event_list;
      /// completed events
      std::list<Event> event_list;
      /// map from construct name to unique id
      FunctionMap function_map;
      /// map from unique id to construct name
      ConfigMap configmap;
      /// total number of constructs profiled
      unsigned int nfunc;

    public:
      ProfilerObj();

      ///
      /// \brief integer only inteface for init
      ///
      /// For compatibility with languages other than C/C++
      ///
      int Init(int id);

      ///
      /// \brief initialization 
      ///
      /// Gets the profiler ready to produce profiles of
      /// the application named "name".  The id parameter
      /// typically specifies the mpi rank of the process 
      /// upon which this instance of the ProfilerObj lives.
      ///
      int Init(const std::string &name,int id);

      ///
      /// \brief mark construct entry
      ///
      /// Programs call this upon entry into a code construct 
      /// which needs to be profiled.
      ///
      int FunctionEntry(const std::string &name);

      ///
      /// \brief mark construct entry (int only interface)
      ///
      /// Programs call this upon entry into a code construct 
      /// which needs to be profiled.
      ///
      int FunctionEntry(int id);

      ///
      /// \brief mark construct exit
      ///
      /// Programs call this upon exit from a code construct 
      /// which needs to be profiled. The resulting Event is
      /// added to the completed Events list.
      ///
      int FunctionExit(const std::string &name);

      ///
      /// \brief mark construct exit (int only)
      ///
      /// Programs call this upon exit from a code construct 
      /// which needs to be profiled. The resulting Event is
      /// added to the completed Events list.
      ///
      int FunctionExit(int id);

      ///
      /// \brief Force all open profiling Events to close (emergency)
      ///
      /// Programs call this when they want to exit unexpectedly,
      /// for example.
      ///
      int FunctionExitAll();

      ///
      /// \brief dumps closed events, clears memory
      ///
      int  Dump(std::ostream &Ostr);

      ///
      /// \brief Set outstream
      ///
      void SetOut(std::ostream *Os){Out = Os;};

      ///
      /// \brief Set errstream
      ///
      void SetErr(std::ostream *Oe){Err = Oe;};

      ///
      /// \brief Profiling output for serial application
      ///
      void SummarizeSerialExecution(std::ostream &Ostr);

      ///
      /// \brief Writes final even file
      ///
      void WriteEventFile();
    
      ///
      /// ?Not sure?
      ///
      void DumpEvents(std::ostream &Ostr);

      /// 
      /// \brief Ready to finalize?
      ///
      bool  FinalizeReady(){return (open_event_list.size() == 1); };

      /// 
      /// \brief Shut down profiler
      ///
      int  Finalize();

      /// 
      /// \brief Read configuration from file
      ///
      int  ReadConfig(const std::string &fname);

      ///
      /// \brief Read serial event file
      ///
      int  ReadEventsFromFile(const std::string &filename);

      /// 
      /// \brief Read event files from parallel run
      ///
      int  ReadParallelEventFiles(const std::vector<std::string> &infiles,
				  PEventList &par_event_list);

      ///
      /// \brief Profiling output for single parallel run
      ///
      int SummarizeParallelExecution(std::ostream &Ostr,
				     std::ostream &Ouf,
				     PEventList &parallel_event_list);
      ///
      /// \brief Read summary files from multiple parallel runs
      ///
      int ReadSummaryFiles(const std::vector<std::string> &input_files,
			   ScalaMap &scala_map);

      ///
      /// \brief Build scalability stats for multiple parallel runs
      ///
      int PopulateScalaMap(ScalaMap &scala_map,
			   ScalaStatMap &scala_statmap,
			   bool is_scaled);
      ///
      /// \brief Scalability analysis output for multiple parallel runs
      ///
      int ScalabilitySummary(ScalaStatMap &scala_statmap,std::ostream &Out);

    private:
      /// whether the profiler has been initialized
      bool _initd;
      /// whether the profiler has been finalized
      bool _finalized;

    };
  };
};

#endif
