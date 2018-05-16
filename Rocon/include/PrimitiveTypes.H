///
/// \file
/// \ingroup irad_group
/// \brief Simple utility classes
///
#ifndef _PRIMITIVETYPES_H_
#define _PRIMITIVETYPES_H_
#include <vector>
#include <list>
#include <map>
#include <cassert>
#include <sstream>

namespace IRAD {
  
  /// Basic primitive types
  namespace Primitive {

    ///
    /// \brief Indexed container
    ///
    template<typename IndexType,typename ValType>
    class IndexedVectorObj : public std::vector<std::pair<IndexType,ValType> >
    {};

    ///
    /// \brief Type agnostic dummy stepper
    ///
    template<typename DataType>
    class DummyStepper {
    private:
      DataType _l;
      DataType _cur;
      DataType _interval;
      size_t   _N;
    public:
      virtual ~DummyStepper(){};
      DummyStepper() : _l(static_cast<DataType>(0)), 
		       _cur(static_cast<DataType>(0)),
		       _interval(static_cast<DataType>(1)),
		       _N(0)
      {};
      DummyStepper(const DataType &lval) : _l(lval), 
					   _cur(lval), 
					   _interval(static_cast<DataType>(1)),
					   _N(0)
      {};
      DummyStepper(const DataType &lval,const DataType &inint) : 
	_l(lval), 
	_cur(lval), 
	_interval(inint),
	_N(0)
      {};
      DataType SetInterval(const DataType &inint) 
      { 
	_interval = inint;
	return(_cur);
      };
      DataType Interval() const { return(_interval); };
      size_t Reset(const DataType &lval)
      {
	_l = lval;
	_cur=lval;
	_N = 0;
	return(_N);
      };
      size_t Reset()
      { 
	_cur = _l;
	_N = 0; 
	return(_N);
      };
      size_t GetStep(){return(_N);};
      DataType Value() const { return(_cur); };
      virtual DataType operator++() 
      { 
	_N++;
	return(_cur+=_interval);
      };
      virtual DataType operator++(int) 
      { 
	DataType rval = _cur;
	_cur+=_interval;
	_N++;
	return(rval);
      };
      virtual DataType operator--() 
      { 
	_N--;
	return(_cur-=_interval);
      };
      virtual DataType operator--(int) 
      { 
	DataType rval = _cur;
	_cur-=_interval;
	_N--;
	return(rval);
      };
    };


    ///
    /// \brief Type agnostic stepper
    ///
    template<typename DataType>
    class Stepper {
    private:
      bool                            _good;
      bool                            _done;
      std::pair<DataType,DataType>    _limits;
      std::vector<DataType>           _values;
      typename std::vector<DataType>::iterator _curstep;
    public:
      virtual ~Stepper(){};
      Stepper() : _good(false), _done(false) {};
      Stepper(const DataType &l,const DataType &u,size_t N)
      {
	Init(l,u,N);
      }
      Stepper(const DataType &l,const DataType &u,const DataType &interval)
      {
	Init(l,u,interval);
      };
      virtual DataType Init(const DataType &l,const DataType &u,size_t N)
      {
	_limits.first  = l;
	_limits.second = u;
	_values.resize(N);
	DataType interval = (u - l)/static_cast<DataType>(N);
	for(size_t i = 0;i < (N-1);i++)
	  _values[i] = l + i*interval;
	_values[N-1] = u;
	_curstep = _values.begin();
	_good = true;
	_done = false;
	return(*_curstep);
      };
      virtual DataType Init(const DataType &l,const DataType &u,const DataType &interval)
      {
	_limits.first  = l;
	_limits.second = u;
	size_t N =  static_cast<size_t>((u - l)/interval);
	for(size_t i = 0;i < (N-1);i++)
	  _values[i] = l + i*interval;
	_values[N-1] = u;
	_curstep = _values.begin();
	_good = true;
	_done = false;
	return(*_curstep);
      };
      virtual size_t GetStep() const { return(_curstep-_values.begin()+1); };
      virtual DataType GetNext() { return (*++_curstep);};
      virtual DataType Value() const { return (*_curstep);};
      virtual DataType Value(size_t n) const 
      { 
	assert(n <= _values.size());
	return (_values[n-1]);
      };
      virtual DataType Restart()
      {
	_curstep = _values.begin(); 
	return(*_curstep);
      };
      virtual DataType Refine(size_t N){return Value();};
      virtual DataType Reset(size_t N) {return(Init(_limits.first,_limits.second,N));};
      virtual void Finish() 
      {
	_done = true;
	_curstep = _values.end();
      };
      virtual DataType operator++() 
      { 
	return(*++_curstep);
      };
      virtual DataType operator++(int) { 
	return(*(_curstep++));};
      virtual DataType operator--() { return(*--_curstep);};
      virtual DataType operator--(int) { return(*(_curstep--));};
      virtual bool Done() const {
	if(_curstep == _values.end())
	  return(true);
	return (_done); 
      };
      virtual bool Good() const { return (_good); };
      virtual bool Bad() const { return  (!_good); };
      virtual void SetBad() { _good = false; };
      virtual const std::vector<DataType> &Values() { return(_values); };
    };

    ///
    /// \brief Simple key-value pair
    ///
    template<typename K,typename V>
    class KeyValuePairObj : public std::pair<K,V>
    {
    public:
      KeyValuePairObj(const std::pair<K,V> &p){
	this->first = p.first;
	this->second = p.second;
      };
      KeyValuePairObj(const K &k,const V &v){
	this->first = k;
	this->second = v;
      };
      KeyValuePairObj(){};
      K Key() const{ return this->first; };
      V Value() const{ return this->second; };
      K &Key() { return this->first; };
      V &Value() { return this->second; };
      K &Key(const K &k) { this->first = k; return(this->first); };
      V &Value(const V &v) { this->second = v; return(this->second); };
      //  virtual std::ostream &Dump(std::ostream &Ostr){
      //    Ostr << Key() << ":" << Value();
      //    return(Ostr);
      //  };
    };

    ///
    /// \brief Generic vector
    ///
    template<typename DataType>
    class Vector {
    private:
      bool _mydata;
      size_t _size;
    protected:
      DataType *D;
    public:
      Vector(size_t size = 0)
	: _size(size)
      {
	_create();
      };
      Vector(const Vector &p1,const Vector &p2)
      {
	assert(p1._size == p2._size);
	_size = p1._size;
	_create();
	Copy(p2-p1);
      };
      Vector(DataType *a,size_t size=0)
	: _size(size)
      {
	_mydata = false;
	D = a;
      };
      Vector(const Vector<DataType> &v)
      {
	_size = v._size;
	_create();
	Copy(v);
      };
      Vector(Vector<DataType> &v)
      {
	_mydata = false;
	_size = v._size;
	D = v.D;
      };
      Vector(std::vector<DataType> &v)
      {
	_mydata = false;
	_size = v.size();
	D = &v[0];
      };
      Vector(const std::vector<DataType> &v)
      {
	_size = v.size();
	_create();
	size_t i = 0;
	typename std::vector<DataType>::const_iterator vi = v.begin();
	while(vi != v.end())
	  D[i++] = *vi++;
      };
      ~Vector()
      {
	_destroy();
      };
      Vector &Copy(const Vector<DataType> &v)
      {
	for(size_t i = 0;i < _size;i++)
	  D[i] = v[i];
	return(*this);
      };
      //     C3Vector &Init(const double *a)
      //     {
      //       V[0] = a[0];
      //       V[1] = a[1];
      //       V[2] = a[2];
      //       return(*this);
      //     };
      //     C3Vector &Init(const C3Vector &v)
      //     {
      //       V[0] = v[0];
      //       V[1] = v[1];
      //       V[2] = v[2];
      //       return(*this);
      //     };
      //     C3Vector &Init(double ix,double iy,double iz)
      //     {
      //       V[0] = ix;
      //       V[1] = iy;
      //       V[2] = iz;
      //       return(*this);
      //     };
      //     C3Vector &Init(const C3Vector &p1,const C3Vector &p2)
      //     {
      //       Copy(p2 - p1);
      //       return(*this);
      //     };
      //     C3Vector &Init()
      //     {
      //       V[0] = V[1] = V[2] = 0.0;
      //       return(*this);
      //     };
      //     C3Vector &Clear()
      //     {
      //       V[0] = V[1] = V[2] = 0.0;
      //       return(*this);
      //     };
      //     double Norm() const
      //     {
      //       return(sqrt((V[0] * V[0]) + (V[1] * V[1]) + (V[2] * V[2])));
      //     };
      //     double Mag() const
      //     {
      //       return (Norm());
      //     };
      //     double Mag2() const
      //     {
      //       return ((V[0]*V[0])+(V[1]*V[1])+(V[2]*V[2]));
      //     };
      //     C3Vector &Normalize()
      //     {
      //       double n = Norm();
      //       V[0] /= n;
      //       V[1] /= n;
      //       V[2] /= n;
      //       return(*this);
      //     };
      //     C3Vector Unit() const
      //     {
      //       C3Vector v(*this);
      //       return(v.Normalize());
      //     };
      //     ///
      //     /// \brief Dot Product
      //     ///
      //     double operator*(const C3Vector &v2) const
      //     {
      //       return((V[0] * v2.V[0]) +
      // 	     (V[1] * v2.V[1]) +
      // 	     (V[2] * v2.V[2]));
      //     };
      //     ///
      //     /// \brief Scaling
      //     ///
      //     C3Vector operator*(double num) const
      //     {
      //       C3Vector rv(num*V[0],num*V[1],num*V[2]);
      //       return(rv);
      //     };
      //     C3Vector &operator*=(double num) 
      //     {
      //       V[0] *= num;
      //       V[1] *= num;
      //       V[2] *= num;
      //       return(*this);
      //     };
      //     ///
      //     /// \brief Cross Product
      //     ///
      //     C3Vector operator%(const C3Vector &v2) const
      //     {
      //       C3Vector v;
      //       v.V[0] = ((v2.V[2] * V[1]) - (v2.V[1] * V[2]));
      //       v.V[1] = ((v2.V[0] * V[2]) - (v2.V[2] * V[0]));
      //       v.V[2] = ((v2.V[1] * V[0]) - (v2.V[0] * V[1]));
      //       return(v);
      //     };
      //     C3Vector &operator=(const C3Vector &v2)
      //     {
      //       return(Copy(v2));
      //     };
      //     C3Vector &operator%=(const C3Vector &v2)
      //     {
      //       double t0 = V[0];
      //       double t1 = V[1];
      //       double t2 = V[2];
      //       V[0] = ((v2.V[2] * t1) - (v2.V[1] * t2));
      //       V[1] = ((v2.V[0] * t2) - (v2.V[2] * t0));
      //       V[2] = ((v2.V[1] * t0) - (v2.V[0] * t1));
      //       return(*this);
      //     };
      //     C3Vector operator+(const C3Vector &v2) const
      //     {
      //       C3Vector v(V[0] + v2.V[0], V[1] + v2.V[1], V[2] + v2.V[2]);
      //       return(v);
      //     };
      //     C3Vector operator-(const C3Vector &v2) const
      //     {
      //       return((*this) + (-1.0 * v2));
      //     };
      //     C3Vector &operator+=(const C3Vector &v2)
      //     {
      //       V[0] += v2.V[0];
      //       V[1] += v2.V[1];
      //       V[2] += v2.V[2];
      //       return(*this);
      //     };
      //     C3Vector &operator-=(const C3Vector &v2)
      //     {
      //       V[0] -= v2.V[0];
      //       V[1] -= v2.V[1];
      //       V[2] -= v2.V[2];
      //       return(*this);
      //     };
      //     double &x() 
      //     {
      //       return(V[0]);
      //     };
      //     double x() const
      //     {
      //       return(V[0]);
      //     };
      //     C3Vector &x(double i)
      //     {
      //       V[0] = i;
      //       return(*this);
      //     };
      //     double &y()
      //     {
      //       return(V[1]);
      //     };
      //     double y() const
      //     {
      //       return(V[1]);
      //     };
      //     C3Vector &y(double i)
      //     {
      //       V[1] = i;
      //       return(*this);
      //     };
      //     double &z()
      //     {
      //       return(V[2]);
      //     };
      //     double z() const
      //     {
      //       return(V[2]);
      //     };
      //     C3Vector &z(double i)
      //     {
      //       V[2] = i;
      //       return(*this);
      //     };
      //     const double &operator[](unsigned int i) const
      //     {
      //       assert(i>0 && i<4);
      //       return(V[i-1]);
      //     };
      //     double &operator[](unsigned int i)
      //     {
      //       assert(i>0 && i<4);
      //       return(V[i-1]);
      //     };
      //     bool operator==(const C3Vector &v) const
      //     {
      //       return( (fabs(V[0] - v.V[0]) < TOL) &&
      // 	      (fabs(V[1] - v.V[1]) < TOL) &&
      // 	      (fabs(V[2] - v.V[2]) < TOL) );
      //     };
      //     bool operator<(const C3Vector &v) const
      //     {
      //       return (Mag2() < v.Mag2());
      //     };
      //     bool operator>(const C3Vector &v) const
      //     {
      //       return (Mag2() > v.Mag2());
      //     };
      //     bool operator<=(const C3Vector &v) const
      //     {
      //       return ( (*this < v) ||
      // 	       (*this == v));
      //     };
      //     bool operator>=(const C3Vector &v) const
      //     {
      //       return ( (*this > v) ||
      // 	       (*this == v) );
      //     };
      //     void DestroyData() { _destroy(); };
      //     void SetData(double *d) { 
      //       V = d; 
      //       _mydata = false;
      //     };
    private:
      void _create()
      {
	if(_size > 0){
	  D = new DataType [_size];
	  _mydata = true;
	}
	else{
	  _size = 0;
	  _mydata = false;
	  D = NULL;
	}
      };
      void _destroy()
      {
	if(_mydata){
	  delete [] D;
	  _mydata = false;
	  D = NULL;
	}
      };
    };

    typedef char byte;
    typedef unsigned char ubyte;


    typedef unsigned int IndexType;
    typedef std::vector<std::list<IndexType> > IndexVecList;
    typedef std::list<IndexType> IndexList;
    typedef std::vector<IndexType> IndexVec;
    template<typename DataType,typename IndexType = IndexType> 
    class MultiContainer {
    public:
      typedef std::vector<std::vector<DataType> > VecVec;
      typedef std::vector<std::map<DataType,IndexType> > VecMap;
      typedef std::vector<std::list<DataType> > VecList;
    };
    //
    //  No such thing as typedef templates - but there 
    //  should be!
    // 
    //  template<typename DataType>
    //  typedef std::vector<std::map<DataType> > VecMap;
    //  typedef std::vector<std::vector<DataType> > MultiVec;
  };
  
  /// Metadata constructs
  namespace Meta {
  
    class DataObject : public std::map<std::string,std::string>
    {
      //    friend std::ostream &operator<<(std::ostream &oSt,
      //				    const Meta::DataObject &md);
      //    friend std::istream &operator>>(std::istream &iSt,
      //				    Meta::DataObject &pv);
    private:
      Meta::DataObject *_mymd;
      bool _own;
    public:
      DataObject() : std::map<std::string,std::string>(), _mymd(NULL), _own(false)
      {};
      ~DataObject() { 
	if(_mymd && _own)
	  delete _mymd;
      };
      Meta::DataObject *InternalMetaData() { return (_mymd); };
      const Meta::DataObject *InternalMetaData() const { return ((const Meta::DataObject *)_mymd); };
      Meta::DataObject *SetInternalMetaData(Meta::DataObject *indata = NULL){
	if(_mymd && _own)
	  delete _mymd;
	if(indata){
	  _mymd = indata;
	  _own = false;
	}
	else {
	  _mymd = new Meta::DataObject;
	  _own = true;
	}
	return(_mymd);
      };
    
      std::string GetValue(const std::string &key) const
      {
	std::string retval;
	Meta::DataObject::const_iterator vi = this->find(key);
	if(vi != this->end())
	  retval = vi->second;
	return(retval);
      };

      template<typename T>
      T GetValue(const std::string &key) const
      {
	T retval;
	Meta::DataObject::const_iterator vi = this->find(key);
	if(vi != this->end()){
	  std::istringstream Istr(vi->second);
	  Istr >> retval;
	}
	return(retval);
      };
    
      template<typename T>
      T ReadValue(const std::string &key) const
      {
	T retval;
	Meta::DataObject::const_iterator vi = this->find(key);
	if(vi != this->end()){
	  std::istringstream Istr(vi->second);
	  Istr.read(reinterpret_cast<char *>(&retval),sizeof(T));
	}
	return(retval);
      };
    
      template<typename T>
      std::vector<T> GetValueVector(const std::string &key) const
      {
	std::vector<T> retval;
	Meta::DataObject::const_iterator vi = this->find(key);
	if(vi != this->end()){
	  std::istringstream Istr(vi->second);
	  T tmpval;
	  while(Istr >> tmpval)
	    retval.push_back(tmpval);
	}
	return(retval);
      };

      template<typename T>
      std::vector<T> ReadValueVector(const std::string &key) const
      {
	std::vector<T> retval;
	Meta::DataObject::const_iterator vi = this->find(key);
	if(vi != this->end()){
	  std::istringstream Istr(vi->second);
	  T tmpval;
	  while(Istr.read(reinterpret_cast<char *>(&tmpval),sizeof(T)))
	    retval.push_back(tmpval);
	}
	return(retval);
      };
      //    int ReadFromStream(std::istream &Is);
      //    std::string Param(const std::string &Key) const;
    };
    //  std::ostream &operator<<(std::ostream &oSt,const Util::Parameters &pv);
    //  std::istream &operator>>(std::istream &iSt,Util::Parameters &pv);
  };
};
#endif
