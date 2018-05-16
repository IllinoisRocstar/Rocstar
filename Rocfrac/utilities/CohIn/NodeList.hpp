#ifndef NODELIST_HPP
#define NODELIST_HPP

#include "general.hpp"
#include "Node.hpp"

class NodeList {
public:

  NodeList();
  NodeList( const NodeList& olist );
  ~NodeList();

  const NodeList& operator=( const NodeList& olist);

  void append(Node* val);
  void insert(Node* val);
  void insert_first(Node* val);

  void reset(); 
  void next();

  Node* get();
  Node* remove();

  boolean empty();
  int size();

  boolean move_to(Node* val);
  int index() const;
  void index( int ind );

private:

  struct Elem {
   Node *d_val;
   Elem *d_next;

   Elem() : d_val(0), d_next(0) {}
   Elem(Node* val, Elem* next = 0) : 
     d_val(val), d_next(next) {}

   ~Elem() {
     if( d_next ) delete d_next;
   }
  };

  Elem* d_first;
  Elem* d_current;

  int d_size;
};


inline void NodeList::insert(Node* val) { // at d_current
  if( !d_current && d_first){
    cerr << " Gevalt - insert\n";
    return;
  }
  d_size++;
  if( !d_current ){
    d_first = new NodeList::Elem( val );
    d_current = d_first;
    return;
  }
  d_current->d_next = new Elem( val, d_current->d_next);
}

inline void NodeList::insert_first(Node* val) { // at d_current
  d_size++;
  d_first = new NodeList::Elem( val, d_first);
  d_current = d_first;
}


inline void NodeList::reset() {
  d_current = d_first;
}

inline void NodeList::next() {
  d_current = d_current->d_next;
  if( !d_current ){
    d_current = d_first;
  }
}

inline Node* NodeList::get() {
  return d_current->d_val;
}

inline boolean NodeList::empty(){
  return ( d_size <= 0 ? TRUE : FALSE);
}

inline int NodeList::size(){
  return d_size;
}


#endif
  

 

