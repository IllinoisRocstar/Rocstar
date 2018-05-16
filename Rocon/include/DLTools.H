///
/// \file
/// \ingroup support
/// \brief Dynamic Loading Tools interface
///
#ifndef _DL_TOOLS_H_
#define _DL_TOOLS_H_
#include <dlfcn.h>
#include <string>
#include <sstream>

namespace Sys {

  typedef void* dlHandleType;

  template<typename ObjectType>
  class ObjectFactory {
  protected:
    ObjectType* (*_create)();
    void        (*_destroy)(ObjectType*);
  public:
    ObjectFactory() : _create(NULL), _destroy(NULL)
    {};
    virtual ObjectType*          Create(){return((*_create)());                    };
    virtual ObjectType*        (*CreatePtr())(){return(_create);                   };
    virtual void               (*DestroyPtr())(ObjectType*){return(_destroy);      };
    virtual void                 Destroy(ObjectType *objptr){(*_destroy)(objptr);  };
    virtual void                 SetCreate(ObjectType* (*increate)()){_create = increate;};
    virtual void                 SetDestroy(void (*indestroy)(ObjectType*)){_destroy = indestroy; };
    virtual ~ObjectFactory(){};
  };

  template<typename ObjectType>
  class DynamicObjectLib : public std::pair<std::string,dlHandleType>
  {
  protected:
    std::string                 _path;
    ObjectFactory<ObjectType>   _factory;
  public:
    virtual dlHandleType       Handle() const { return(this->second);      };
    virtual std::string        Name()   const { return(this->first);       };
    virtual std::string        File()   const { return(_path);             };
    virtual std::string       &Name() { return(this->first);               };
    virtual std::string       &File() { return(_path);                     };
    virtual std::string        Error(){ return(std::string(dlerror()));    };
    virtual ObjectFactory<ObjectType> &Factory() { return(_factory);       };
    virtual ObjectFactory<ObjectType>  Factory() const { return(_factory); };
    virtual ObjectType* CreateObject(const std::string &objectname)
    {
      ObjectType* retval = static_cast<ObjectType*>(NULL);
      if(this->second){
	std::ostringstream Ostr;
	Ostr << "create_" << objectname;
	_factory.SetCreate((ObjectType * (*)())dlsym(this->second,Ostr.str().c_str()));
	if(_factory.CreatePtr() != NULL)
	  retval = _factory.Create();
      }
      return(retval);
    };
    virtual void DestroyObject(const std::string &objectname,ObjectType *objptr)
    {
      if(this->second){
	std::ostringstream Ostr;
	Ostr << "destroy_" << objectname;
	_factory.SetDestroy((void (*)(ObjectType *))dlsym(this->second,Ostr.str().c_str()));
	if(_factory.DestroyPtr() != NULL)
	  _factory.Destroy(objptr);
      }
    };
    virtual int Load(const std::string &name,const std::string &path)
    {
      void *handle = dlopen(path.c_str(),RTLD_LAZY);
      if(handle){
	if(this->second) 
	  dlclose(this->second);
	this->first  = name;
	this->second = handle;
      }
      else
	return(1);
      return(0);
    };
    virtual void Unload()
    { 
      if(this->second) 
	dlclose(this->second);
      this->second = NULL;
      this->first.erase();
      _path.erase();
      _factory.SetCreate(NULL);
      _factory.SetDestroy(NULL);
    };
    virtual ~DynamicObjectLib(){ if(this->second) dlclose(this->second);};
  };
}
#endif
