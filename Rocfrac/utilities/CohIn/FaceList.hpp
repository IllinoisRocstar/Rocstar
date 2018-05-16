#ifndef FACELIST_HPP
#define FACELIST_HPP

#include "general.hpp"
#include "Face.hpp"

class FaceList {
public:

  FaceList();
  FaceList( const FaceList& olist );
  ~FaceList();

  const FaceList& operator=( const FaceList& olist);

  void append(Face* val);
  void insert(Face* val);
  void insert_first(Face* val);

  void reset(); 
  void next();

  Face* get();
  Face* remove();

  boolean empty();
  int size();

  boolean move_to(Face* val);
  int index() const;
  void index( int ind );

private:

  struct Elem {
   Face *d_val;
   Elem *d_next;

   Elem() : d_val(0), d_next(0) {}
   Elem(Face* val, Elem* next = 0) : 
     d_val(val), d_next(next) {}

   ~Elem() {
     if( d_next ) delete d_next;
   }
  };

  Elem* d_first;
  Elem* d_current;

  int d_size;
};


inline void FaceList::insert(Face* val) { // at d_current
  if( !d_current && d_first){
    cerr << " Gevalt - insert\n";
    return;
  }
  d_size++;
  if( !d_current ){
    d_first = new FaceList::Elem( val );
    d_current = d_first;
    return;
  }
  d_current->d_next = new Elem( val, d_current->d_next);
}

inline void FaceList::insert_first(Face* val) { // at d_current
  d_size++;
  d_first = new FaceList::Elem( val, d_first);
  d_current = d_first;
}


inline void FaceList::reset() {
  d_current = d_first;
}

inline void FaceList::next() {
  d_current = d_current->d_next;
  if( !d_current ){
    d_current = d_first;
  }
}

inline Face* FaceList::get() {
  return d_current->d_val;
}

inline boolean FaceList::empty(){
  return ( d_size <= 0 ? TRUE : FALSE);
}

inline int FaceList::size(){
  return d_size;
}


#endif
  

 

