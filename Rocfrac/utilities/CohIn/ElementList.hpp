#ifndef ELEMENTLIST_HPP
#define ELEMENTLIST_HPP

#include "general.hpp"
#include "Element.hpp"

class ElementList {
public:

  ElementList();
  ElementList( const ElementList& olist );
  ~ElementList();

  const ElementList& operator=( const ElementList& olist);

  void append(Element* val);
  void insert(Element* val);
  void insert_first(Element* val);

  void reset(); 
  void next();

  Element* get();
  Element* remove();

  boolean empty();
  int size();

  boolean move_to(Element* val);
  int index() const;
  void index( int ind );

private:

  struct Elem {
   Element *d_val;
   Elem *d_next;

   Elem() : d_val(0), d_next(0) {}
   Elem(Element* val, Elem* next = 0) : 
     d_val(val), d_next(next) {}

   ~Elem() {
     if( d_next ) delete d_next;
   }
  };

  Elem* d_first;
  Elem* d_current;

  int d_size;
};


inline void ElementList::insert(Element* val) { // at d_current
  if( !d_current && d_first){
    cerr << " Gevalt - insert\n";
    return;
  }
  d_size++;
  if( !d_current ){
    d_first = new ElementList::Elem( val );
    d_current = d_first;
    return;
  }
  d_current->d_next = new Elem( val, d_current->d_next);
}

inline void ElementList::insert_first(Element* val) { // at d_current
  d_size++;
  d_first = new ElementList::Elem( val, d_first);
  d_current = d_first;
}


inline void ElementList::reset() {
  d_current = d_first;
}

inline void ElementList::next() {
  d_current = d_current->d_next;
  if( !d_current ){
    d_current = d_first;
  }
}

inline Element* ElementList::get() {
  return d_current->d_val;
}

inline boolean ElementList::empty(){
  return ( d_size <= 0 ? TRUE : FALSE);
}

inline int ElementList::size(){
  return d_size;
}


#endif
  

 

