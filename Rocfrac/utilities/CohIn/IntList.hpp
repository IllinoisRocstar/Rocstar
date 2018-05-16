#ifndef INTLIST_HPP
#define INTLIST_HPP

#include "general.hpp"
#include <assert.h>

class IntList {
public:

  IntList();
  IntList( const IntList& olist );
  ~IntList();

  const IntList& operator=( const IntList& olist);

  void append(int ed);
  void insert(int ed);
  void insert_first(int ed);

  void reset(); 
  void next();

  int& get();
  int remove();

  boolean empty();
  int size();

  boolean move_to(int ed);
  int index() const;
  void index( int ind );

private:

  struct Elem {
   int d_val;
   Elem *d_next;

   Elem() : d_val(0), d_next(0) {}
   Elem(int val, Elem* next = 0) : 
     d_val(val), d_next(next) {}

   ~Elem() {
     if( d_next ) delete d_next;
   }
  };

  Elem* d_first;
  Elem* d_current;

  int d_size;
};


inline void IntList::insert(int val) { // at d_current
  assert( d_current || !d_first);
  d_size++;
  if( !d_current ){
    d_first = new IntList::Elem( val );
    d_current = d_first;
    return;
  }
  d_current->d_next = new Elem( val, d_current->d_next);
}

inline void IntList::insert_first(int val) { // at d_current
  d_size++;
  d_first = new IntList::Elem( val, d_first);
  d_current = d_first;
}


inline void IntList::reset() {
  d_current = d_first;
}

inline void IntList::next() {
  d_current = d_current->d_next;
  if( !d_current ){
    d_current = d_first;
  }
}

inline int& IntList::get() {
  return d_current->d_val;
}

inline boolean IntList::empty(){
  return ( d_size <= 0 ? TRUE : FALSE);
}

inline int IntList::size(){
  return d_size;
}


#endif
  

 

