
#include "FaceList.hpp"


FaceList::FaceList() :
  d_first(0), 
  d_current(0),
  d_size(0) {}

FaceList::FaceList( const FaceList& olist) :
  d_size(olist.d_size) {

  // just add all the list
  Elem* oval = olist.d_first;
  Elem* prev(0);
  while( oval ){
    Elem *nelem = new FaceList::Elem( oval->d_val );
    if( prev ){
      prev->d_next = nelem;
    }
    else {
      d_first = nelem;
    }
    prev = nelem;
    oval = oval->d_next;
  }
  d_current = d_first;
}


FaceList::~FaceList() {
  if( d_first ) delete d_first;
}
    
const FaceList& FaceList::operator=( const FaceList& olist){

  if( d_first ) delete d_first;

  d_size = olist.d_size;
  // just add all the list
  Elem* oval = olist.d_first;
  Elem* prev(0);
  while( oval ){
    Elem *nelem = new FaceList::Elem( oval->d_val );
    if( prev ){
      prev->d_next = nelem;
    }
    else {
      d_first = nelem;
    }
    prev = nelem;
    oval = oval->d_next;
  }
  d_current = d_first;
  return *this;
}

void FaceList::append(Face* val){

  Elem* bef(0);
  Elem* af( d_first );
  while( af ){
    bef = af;
    af = af->d_next;
  }
  if( bef ){
    bef->d_next = new FaceList::Elem( val );
  }
  else {
    d_first = new Elem( val );
    d_current = d_first;
  }
  d_size++;
}

Face* FaceList::remove(){
  Face* val = d_current->d_val;
  
  d_size--;

  Elem* bef(0);
  Elem* af( d_first );
  while( af != d_current){  
    bef = af;
    af = af->d_next;
  }
  if( !bef ){
    d_first = d_current->d_next;
  }
  else {
    bef->d_next = d_current->d_next;
  }
  d_current = d_current->d_next;
  af->d_next = 0;
  delete af;// the old current - remove only it
  if( !d_current ) d_current = d_first;
  return val;
}



boolean FaceList::move_to(Face* val){
  d_current = d_first;
  while( d_current && d_current->d_val != val ){
    d_current = d_current->d_next;
  }
  return( d_current ? TRUE : FALSE );
}

int  FaceList::index() const {

  int ind = 0;
  Elem* pass = d_first;
  while( pass != d_current ){
    ind++;
    pass = pass->d_next;
  }
  return ind;
}

void FaceList::index( int ind ){

  d_current = d_first;
  int i;
  for( i = 0; i < ind; i++ ){
    d_current = d_current->d_next;
  }
}


