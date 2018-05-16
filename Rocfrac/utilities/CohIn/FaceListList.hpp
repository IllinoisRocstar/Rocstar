#ifndef FACELISTLIST_HPP
#define FACELISTLIST_HPP

#include "general.hpp"
#include "FaceList.hpp"

class FaceListList {
public:

  FaceListList();
  FaceListList( const FaceListList& olist );
  ~FaceListList();

  const FaceListList& operator=( const FaceListList& olist);

  void append(FaceList* val);
  void insert(FaceList* val);
  void insert_first(FaceList* val);

  void reset(); 
  void next();

  FaceList* get();
  FaceList* remove();

  boolean empty();
  int size();

  boolean move_to(FaceList* val);
  int index() const;
  void index( int ind );

private:

  struct Elem {
   FaceList *d_val;
   Elem *d_next;

   Elem() : d_val(0), d_next(0) {}
   Elem(FaceList* val, Elem* next = 0) : 
     d_val(val), d_next(next) {}

   ~Elem() {
     if( d_next ) delete d_next;
   }
  };

  Elem* d_first;
  Elem* d_current;

  int d_size;
};


inline void FaceListList::insert(FaceList* val) { // at d_current
  if( !d_current && d_first){
    cerr << " Gevalt - insert\n";
    return;
  }
  d_size++;
  if( !d_current ){
    d_first = new FaceListList::Elem( val );
    d_current = d_first;
    return;
  }
  d_current->d_next = new Elem( val, d_current->d_next);
}

inline void FaceListList::insert_first(FaceList* val) { // at d_current
  d_size++;
  d_first = new FaceListList::Elem( val, d_first);
  d_current = d_first;
}


inline void FaceListList::reset() {
  d_current = d_first;
}

inline void FaceListList::next() {
  d_current = d_current->d_next;
  if( !d_current ){
    d_current = d_first;
  }
}

inline FaceList* FaceListList::get() {
  return d_current->d_val;
}

inline boolean FaceListList::empty(){
  return ( d_size <= 0 ? TRUE : FALSE);
}

inline int FaceListList::size(){
  return d_size;
}


#endif
  

 

