#ifndef TSTT_UTIL_HPP
#define TSTT_UTIL_HPP

#include <sidl_cxx.hh>
#include "TSTTB.hh"

namespace Mesquite {

static inline msq_std::string process_tstt_error( TSTTB::Error &tstt_err )
{
  msq_std::string str;
  msq_std::string result("TSTT ERROR: ");
  result += tstt_err.getNote();
  MSQ_DBGOUT(1) << "TSTT Error:" << msq_std::endl;
  MSQ_DBGOUT(1) << tstt_err.getNote() << msq_std::endl;
  tstt_err.getDescription(str);
  MSQ_DBGOUT(1) << str << msq_std::endl;
  MSQ_DBGOUT(1) << tstt_err.getTrace() << msq_std::endl;
  return result;
}

template <class T> static inline
T* convert_from_sidl_vector( sidl::array<T>& array )
{  return reinterpret_cast<T*>(array._get_ior()->d_firstElement); }

template <class T> static inline sidl::array<T> alloc_sidl_vector( size_t size )
{
  int32_t lower = 0;
  int32_t upper = size - 1;
  return sidl::array<T>::createCol( 1, &lower, &upper );
}

template <class T> static inline sidl::array<T> alloc_sidl_vector( size_t size, T init )
{
  sidl::array<T> result = alloc_sidl_vector<T>(size);
  T* ptr = convert_from_sidl_vector( result );
  for ( T* const end = ptr + size; ptr != end; ++ptr)
    *ptr = init;
  return result;
}

template <class S, class T> static inline void copy_from_sidl( sidl::array<S>& source,
                                                        T* target )
{
  typename sidl::array<S>::iterator i = source.begin();
  for (; i != source.end(); ++i, ++target)
    *target = (T)*i;
}

template <class T> static inline 
sidl::array<T> convert_to_sidl_vector( T* array, size_t size )
{
  sidl::array<T> result;
  int32_t lower = 0, upper = size - 1, stride = 1;
  result.borrow( array, 1, &lower, &upper, &stride );
  return result;
}

} // namespace Mesquite

#endif
