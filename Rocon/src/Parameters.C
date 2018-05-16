///
/// @file
/// @ingroup irad_group
/// @brief Parameters object implementation
///
#include "Parameters.H"


namespace IRAD {
  namespace Util {

    Util::ParamType *Util::Parameters::ParamPtr(const std::string &key)
    {
      Util::Parameters::const_iterator pi = this->begin();
      while(pi != this->end()){
        unsigned int index = pi - this->begin();
        if(pi->first == key)
          return(&((*this)[index]));
        pi++;
      }
      return(NULL);
    }
    
    std::string Util::Parameters::GetValue(const std::string &key) const
    {
      std::string retval;
      std::string value(this->Param(key));
      if(!value.empty()){
	std::istringstream Istr(value);
	Istr >> retval;
      }
      return(retval);
    }

    std::vector<std::string> Util::Parameters::GetValueVector(const std::string &key) const
    {
      std::vector<std::string> retval;
      std::string value(this->Param(key));
      if(!value.empty()){
	std::istringstream Istr(value);
	std::string tmpval;
	while(Istr >> tmpval)
	  retval.push_back(tmpval);
      }
      return(retval);
    }

    std::istream &Util::Parameters::ReadFromStream(std::istream &Is)
    {
      std::string line;
      int n = 0;
      while(Is){
	std::getline(Is,line);
	n++;
	// Removes any leading whitespace
	std::string::size_type x = line.find('#');
	line = line.substr(0,x);
	if(!line.empty()){
	  x = line.find('=');
	  if(x == std::string::npos)
	    return(Is);
	  Util::ParamType param;
	  std::istringstream Instr(line.substr(0,x));
	  Instr >> param.Key();
	  std::vector<std::string> tokens;
	  line = line.substr(x+1,line.size());
	  TokenizeString(tokens,line);
	  std::ostringstream Ostr;
	  std::vector<std::string>::iterator ti = tokens.begin();
	  if(ti != tokens.end())
	    Ostr << *ti++;
	  while(ti != tokens.end())
	    Ostr << " " << *ti++;
	  param.Value() = Ostr.str();
	  this->push_back(param);
	}
      }
      return(Is);
    }

    std::string Util::Parameters::Param(const std::string &key) const
    {
      std::string empty;
      Util::Parameters::const_iterator pi = this->begin();
      while(pi != this->end()){
	if(pi->first == key)
	  return(pi->second);
	pi++;
      }
      return(empty);
    }
  
    bool Util::Parameters::IsSet(const std::string &key) const
    {
      return(!Param(key).empty());
    };
    

    std::ostream &Util::Parameters::WriteToStream(std::ostream &oSt) const
    {
      Util::Parameters::const_iterator pi = this->begin();
      while(pi != this->end()){
        oSt << *pi++;
        if(pi != this->end())
          oSt << std::endl;
      }
      return(oSt);
    };

    std::ostream &operator<<(std::ostream &oSt,
			     const Util::Parameters &pv)
    {
      return(pv.WriteToStream(oSt));
    }
  
    std::istream &operator>>(std::istream &iSt,
			     Util::Parameters &pv)
    {
      return(pv.ReadFromStream(iSt));
    }
  
  
    std::ostream &operator<<(std::ostream &Ostr,const Util::ParamType &param)
    {
      Ostr << param.Key() << " = " << param.Value();
      return(Ostr);
    }

  };
};


