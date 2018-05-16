///
/// \file
/// \ingroup irad_group
/// \brief Parameters object interface (for config files, etc)
///
#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include "primitive_utilities.H"

namespace IRAD { 

  namespace Util {

    typedef Primitive::KeyValuePairObj<std::string,std::string> ParamType;
    std::ostream &operator<<(std::ostream &Ostr,
			     const Util::ParamType &param);
    
    class Parameters : public std::vector<Util::ParamType>
    
    {
      friend std::ostream &operator<<(std::ostream &oSt,
				      const Util::Parameters &pv);
      friend std::istream &operator>>(std::istream &iSt,
				      Util::Parameters &pv);
    public:
      std::string GetValue(const std::string &key) const;
      std::vector<std::string> GetValueVector(const std::string &key) const;
      
      Util::ParamType *ParamPtr(const std::string &key);

      template<typename T>
      T GetValue(const std::string &key) const
      {
	T retval;
	std::string value(this->Param(key));
	if(!value.empty()){
	  std::istringstream Istr(value);
	  Istr >> retval;
	}
	return(retval);
      };
      
      template<typename T>
      std::vector<T> GetValueVector(const std::string &key) const
      {
	std::vector<T> retval;
	std::string value(this->Param(key));
	if(!value.empty()){
	  std::istringstream Istr(value);
	  T tmpval;
	  while(Istr >> tmpval)
	    retval.push_back(tmpval);
	}
	return(retval);
      };
      virtual int SetValue(const std::string &key,const std::string &value){
        if(ParamPtr(key)){
          ParamPtr(key)->Value(value);
          return(0);
        }
        return(1);
      };
      virtual std::istream &ReadFromStream(std::istream &Is);
      virtual std::ostream &WriteToStream(std::ostream &Os) const;
      //      int ReadFromStream(std::istream &Is);
      std::string Param(const std::string &Key) const;
      bool IsSet(const std::string &Key) const;
      virtual ~Parameters(){};
    };
    std::ostream &operator<<(std::ostream &oSt,const Util::Parameters &pv);
    std::istream &operator>>(std::istream &iSt,Util::Parameters &pv);
  };
};


#endif
