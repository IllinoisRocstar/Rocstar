///
/// @file
/// @ingroup irad_group
/// @brief Basic utility header
///
#ifndef _BASIC_UTILITIES_
#define _BASIC_UTILITIES_

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <algorithm>
#include "PrimitiveTypes.H"


namespace IRAD {

  ///
  /// @brief General support utilities.
  /// @ingroup irad_group
  ///
  namespace Util {
    /// Error type
    typedef int Error;
    //  void Vectorize(std::vector<std::string> &retVal,const char **in);
    //  void Vectorize(std::vector<std::string> &retVal,const char **in,int);
    //  int String2Buf(const std::string &instr,void **buf);
    
    
    //  void GetContentUntil(std::istream &In,
    //		       std::string ret,
    //		       const std::string &etag);
    //  
    //  std::string GetNextContent(std::istream &In);
    
    
    // must be sorted, or this doesn't work out
    template<typename Container>
    bool LessThan(const Container &c1,const Container &c2)
    {
      typename Container::const_iterator c1i = c1.begin();
      typename Container::const_iterator c2i = c2.begin();
      while(c1i != c1.end() && c2i != c2.end()){
	if(*c1i < *c2i)
	  return(true);
	if(*c1i++ > *c2i++)
	  return(false);
      }
      if(c1i == c1.end()){
	if(c2i == c2.end()){
	  return(false);
	}
	else {
	  return(true);
	}
      }
      return(false);
    }
  
  
    // Inverts a container
    template<typename Container>
    void InvertContainer(const Container &c1,Container &c2)
    {
      unsigned int size = c1.size();
      c2.resize(size);
      typename Container::const_iterator c1i = c1.begin();
      for(unsigned int iii = 0;iii < size;iii++)
	c2[*c1i++] = iii;
    } 
  
    ///
    /// \brief Dump container contents
    ///
    /// Dumps the contents of the container to the specified std::ostream 
    /// object.
    ///
    template<typename ContainerType>
    void DumpContents(std::ostream &Ostr,const ContainerType &c,std::string del = "\n"){
      if(!c.empty()){
	typename ContainerType::const_iterator ci = c.begin();
	Ostr << *ci++;
	while(ci != c.end())
	  Ostr << del << *ci++;
      }
    }
  
    //template<typename IndexType,typename ValueType>
    //void Invert(IndexedVectorObj<IndexType,ValueType> &ico){
    //  return;
    //}
  
  
    template<typename ContainerType,typename T>
    void CopyIntoContainer(ContainerType &cont,const T* src,size_t count)
    {
      size_t datasize = sizeof(T);
      cont.resize(count);
      memcpy(&cont[0],src,count*datasize);
    }
  
    template<typename OuterCont,typename OutCont,typename InnerCont,typename MapType>
    void MapElements(OuterCont &src,OutCont &trg,MapType &m)
    {
      trg.resize(src.size());
      typename OuterCont::iterator sci = src.begin();
      typename OutCont::iterator oci = trg.begin();
      while(sci != src.end()){
	typename InnerCont::iterator ici = sci->begin();
	while(ici != sci->end())
	  oci->push_back(m[*ici++]);
	sci++;
	oci++;
      }
    }
  
    /// Return the maximum of all elements of a multicontainer
    template<typename ListContainerType, typename ListType>
    Primitive::IndexType MaxNodeId(const ListContainerType &fc)
    {
      Primitive::IndexType maxid = 0;
      typename ListContainerType::const_iterator lci = fc.begin();
      while(lci != fc.end()){
	typename ListType::const_iterator li = lci->begin();
	while(li != lci->end()){
	  if(maxid < *li){
	    maxid = *li;
	  }
	  li++;
	}
	lci++;
      }
      return(maxid);
    }
  
    /// Populate OutCont with a flat list of entries from a multicontainer
    template<typename OuterCont,typename InnerCont,typename OutCont>
    void Flatten(OuterCont &con,OutCont &ocon){
      typename OuterCont::iterator cci = con.begin();
      ocon.resize(0);
      while(cci != con.end()){
	typename InnerCont::iterator ccni = cci->begin();
	while(ccni != cci->end())
	  ocon.push_back(*ccni++);
	cci++;
      }
    }
  
    /// Return the total number of entries in a multicontainer 
    template<typename OuterContType>
    Primitive::IndexType GetTotalSize(OuterContType &con)
    {
      Primitive::IndexType total_size = 0;
      typename OuterContType::iterator ci = con.begin();
      while(ci != con.end()){
	total_size += ci->size();
	ci++;
      }
      return(total_size);
    }
  
    template<typename OuterContType,typename InnerContType,typename RetCont,typename idxtype>
    void MultiContainer2CSR(RetCont &xadj,RetCont &adj,
			    OuterContType &source)
    {
      Primitive::IndexType number_of_elements = source.size();
      xadj.resize(number_of_elements+1);
      Primitive::IndexType number_entries = GetTotalSize<OuterContType>(source);
      adj.reserve(number_entries);
      typename OuterContType::iterator si = source.begin();
      Primitive::IndexType elm = 0;
      xadj[elm++] = 0;
      while(si != source.end()){
	typename InnerContType::iterator ei = si->begin();
	xadj[elm] = static_cast<idxtype>(si->size() + xadj[elm-1]);
	elm++;
	while(ei != si->end())
	  adj.push_back(static_cast<idxtype>(*ei++-1));
	si++;
      }
    }
  
    /// Given an array of adjacent node lists (like an array
    /// of face connectivities), this function will loop thru
    /// and create a list of unique adjacent nodes for each node
    /// with the nodes that *are actually adjacent in the lists*
    template<typename ListContainerType, typename ListType>
    void CreateAdjacentNodeList(std::vector<std::list<Primitive::IndexType> > &anodelist,
				ListContainerType &fc,Primitive::IndexType nnodes=0)
    {
      if(nnodes == 0) 
	nnodes = MaxNodeId<ListContainerType,ListType>(fc);
      anodelist.resize(nnodes);
      typename ListContainerType::iterator lci = fc.begin();
      while(lci != fc.end())
	{
	  typename ListType::iterator li = lci->begin();
	  while(li != lci->end()){
	    Primitive::IndexType this_node = *li++ - 1;
	    Primitive::IndexType next_node = 0;
	    if(li == lci->end())
	      next_node = *(lci->begin());
	    else
	      next_node = *li++;
	    anodelist[this_node].push_back(next_node);
	    anodelist[next_node-1].push_back(this_node+1);
	  }
	  lci++;
	}
      for(Primitive::IndexType node = 0;node < nnodes;node++)
	{
	  anodelist[node].sort();
	  anodelist[node].unique();
	}
    }
  
    /// Given an array of adjacent node lists (like an array
    /// of face connectivities), this function will loop thru
    /// and create a list of unique adjacent nodes for each node
    /// Note that this is a little different from the AdjNodeList
    /// because in this version *every node in each list* is considered
    /// adjacent.
    template<typename ListContainerType, typename ListType>
    void AdjEList(std::vector<std::list<Primitive::IndexType> > &aelist,
		  ListContainerType &dual_con,unsigned long nel=0)
    {
      if(nel == 0) 
	nel = MaxNodeId<ListContainerType,ListType>(dual_con);
      aelist.resize(nel);
      typename ListContainerType::iterator lci = dual_con.begin();
      while(lci != dual_con.end())
	{
	  typename ListType::iterator li = lci->begin();
	  while(li != lci->end()){
	    Primitive::IndexType this_node = *li++ - 1;
	    typename ListType::iterator li2 = li;
	    while(li2 != lci->end()){
	      Primitive::IndexType nexnode = *li2++ - 1;
	      aelist[this_node].push_back(nexnode+1);
	      aelist[nexnode].push_back(this_node+1);
	    }
	  }
	  lci++;
	}
      for(Primitive::IndexType node = 0;node < nel;node++)
	{
	  aelist[node].sort();
	  aelist[node].unique();
	}
    }

    template<typename ConType,typename IConType>
    Primitive::IndexType NumberOfEdges(ConType &con)
    {
      typename ConType::iterator ci = con.begin();
      Primitive::IndexType node = 0;
      Primitive::IndexType nedges = 0;
      while(ci != con.end())
	{
	  Primitive::IndexType tnode = ++node;
	  typename IConType::iterator eni = ci->begin();
	  while(eni != ci->end())
	    {
	      Primitive::IndexType anode = *eni++;
	      if((tnode) < anode)
		nedges++;
	    }
	  ci++;
	}
      return(nedges);
    }

    template<typename ContainerType,typename Icont>
    void FormGraph(const ContainerType &adjlist)
    {
      // empty for now
    }


    ///
    /// \brief Cyclic test
    ///
    /// Typically, this is used to determine if two representations of
    /// the same face are oriented the same.
    ///
    template<typename Container>
    bool HaveSameOrientation(Container &c1,Container &c2)
    {
      if(c1.size() != c2.size())
	return(false);
      typename Container::iterator c1i = c1.begin();
      typename Container::iterator c2i = 
	std::find(c2.begin(),c2.end(),*c1i);
      if(c2i == c2.end())
	return(false);
      while(c1i != c1.end()){
	if(*c1i != *c2i)
	  return false;
	c1i++;
	c2i++;
	if(c2i == c2.end())
	  c2i = c2.begin();
    
      }
      return(true);
    }

    ///
    /// \brief Anti-cyclic test
    ///
    /// Typically, this is used to determine if two representations of
    /// the same face are oriented the same.
    ///
    template<typename Container>
    bool HaveOppositeOrientation(Container &c1,Container &c2)
    {
      if(c1.size() != c2.size())
	return(false);
      typename Container::iterator c1i = c1.begin();
      typename Container::reverse_iterator c2i = 
	std::find(c2.rbegin(),c2.rend(),*c1i);
      if(c2i == c2.rend())
	return(false);
      while(c1i != c1.end()){
	if(*c1i++ != *c2i++)
	  return false;
	if(c2i == c2.rend())
	  c2i = c2.rbegin();
      }
      return(true);
    }

  
    template<class T>
    struct CmpPairs
    {
      CmpPairs(const std::vector<T> &v): v_(v) {}
      std::vector<T> v_;
      bool operator()(int a, int b){ return v_[a] < v_[b]; }
    };
    
    template<class T>
    CmpPairs<T> CreateCmpPairs(const std::vector<T> & v) 
    { return CmpPairs<T>(v); };
    
    template<class Vals>
    void SortPermutation(const Vals& values, std::vector<unsigned int>& v)
    {
      unsigned int size = values.size(); 
      v.clear(); 
      v.reserve(size);
      for(unsigned int i=0; i < size; ++i)
        v.push_back(i);
      std::sort(v.begin(), v.end(), CreateCmpPairs(values));
    };
  
  
    inline std::string GetNextContent(std::istream &In)
    {
      std::string line;
      std::string ret;
      while(ret.empty() && std::getline(In,line)){
	line = line.substr(0,line.find("#"));
	std::istringstream Istr(line);
	std::string token;
	Istr >> token;
	if(!token.empty() && token[0] != '#')
	  ret.assign(line);
      }
      return(ret);
    }
  
    inline void GetContentUntil(std::istream &In,
				std::string ret,
				const std::string &etag)
    {
      std::ostringstream Ostr;
      std::string line;
      while(std::getline(In,line) && ret.empty()){
	line = line.substr(0,line.find("#"));
	if(!line.empty()){
	  std::istringstream Istr(line);
	  std::string token;
	  Istr >> token;
	  if(!token.empty()){
	    if(token == etag){
	      //	    Ostr << line << std::endl;
	      ret = Ostr.str();
	    }
	    else 
	      Ostr << line << std::endl;
	  }
	}
      }
    }
  
    inline int String2Buf(const std::string &instr,void **buf)
    {
      if(instr.empty())
	return 0;
      int rval = instr.size();
      *buf = new char [rval];
      void *resl = std::memcpy(*buf,instr.c_str(),rval);
      if(!resl)
	return(0);
      return(rval);
    }
  
    ///
    /// \brief Creates space delimited tokens in place
    ///
    inline void Trim(std::string &instr,bool preserve_newline = false)
    {
      std::istringstream Istr(instr);
      std::string line;
      std::ostringstream Ostr;
      while(std::getline(Istr,line)){
	std::istringstream Line(line);
	std::string token;
	Line >> token;
	Ostr << token;
	while(Line >> token)
	  if(!token.empty())
	    Ostr << " " << token;
	if(preserve_newline)
	  Ostr << "\n";
      }
      instr = Ostr.str();
    }
  
    ///
    /// \brief Returns space delimited tokens
    ///
    inline const std::string Trimmed(const std::string &instr,bool preserve_newline = false)
    {
      std::istringstream Istr(instr);
      std::string line;
      std::ostringstream Ostr;
      while(std::getline(Istr,line)){
	std::istringstream Line(line);
	std::string token;
	Line >> token;
	Ostr << token;
	while(Line >> token)
	  if(!token.empty())
	    Ostr << " " << token;
	if(preserve_newline)
	  Ostr << "\n";
      }
      return(Ostr.str());
    }
  
    ///
    /// \brief Strip absolute path
    /// 
    /// Strips the absolute path to recover the base file name or
    /// the executable file name.
    ///
    inline const std::string 
    stripdirs(const std::string &pname)
    {
      std::string retval;
      std::string::size_type x = pname.find("/");
      if(x == std::string::npos)
	return(pname);
      return(pname.substr(pname.find_last_of("/")+1));
    }
  
    ///
    /// \brief Tokenize string
    ///
    /// Breaks source string up into a vector of space delimited tokens
    ///
    inline void
    TokenizeString(std::vector<std::string> &tokens,const std::string &source)
    {
      tokens.resize(0);
      std::istringstream Istr(source);
      std::string token;
      while(Istr >> token)
	tokens.push_back(token);
    }
  
    ///
    /// \brief Tokenize string w/ delimiter
    ///
    /// Breaks source string up into a vector of tokens delimited by specified delimiter
    ///  
    inline void
    TokenizeString(std::vector<std::string> &tokens,const std::string &source,const char delim)
    {
      tokens.resize(0);
      std::string::size_type ssize = source.size();
      std::string::size_type x = 0;
      std::ostringstream Ostr;
      std::string token;
      while(x < ssize){
	if(source[x] != delim)
	  token += source[x];
	if((source[x] == delim) || (x == (ssize-1))){
	  tokens.push_back(token);
	  token.erase();
	}
	x++;
      }
    }
  
    ///
    /// \brief Process a range in the format "t1:t2".
    ///
    template<typename NumType>
    void ProcessRange(NumType &t1,NumType &t2,const std::string stinter)
    {
      if(!stinter.empty()){
        std::string::size_type x = stinter.find(":");
        std::string t1s = stinter;
        std::string t2s = stinter;
        if(x != std::string::npos){
          t1s = stinter.substr(0,x);
          t2s = stinter.substr(x+1);
          if(!t2s.empty()){
            std::istringstream Istr(t2s);
            Istr >> t2;
          }
        }
        else {
          t1s.erase();
        }
        if(!t1s.empty()){
          std::istringstream Istr(t1s);
          Istr >> t1;
        }
        if(!t2s.empty()){
          std::istringstream Istr(t2s);
          Istr >> t2;
        }

      }
    }

    ///
    /// \brief File opener
    ///
    inline int OpenFile(std::ifstream &Inf,const std::string &filename)
    {
      Inf.open(filename.c_str());
      if(!Inf)
	return(-1);
      return 0;
    }
  
    // Utilities
    inline void
    Vectorize(std::vector<std::string> &retVal,const char **in)
    {
      int i = 0;
      while(in[i] != NULL)
	retVal.push_back(in[i++]);
    }
  
    inline void
    Vectorize(std::vector<std::string> &retVal,const char **in,int n)
    {
      //  retVal.resize(0);
      if(n <= 0) return;
      int i = 0;
      while((in[i] != NULL) && i < n)
	retVal.push_back(in[i++]);
    }
  };
  /// @}
};
#endif
