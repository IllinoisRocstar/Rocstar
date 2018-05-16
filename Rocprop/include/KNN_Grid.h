#ifndef SPACEGRID_INCLUDED 
#define SPACEGRID_INCLUDED


#include <iostream>
#include <vector>
#include <algorithm>
#include <ext/hash_map>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <queue>
#include "NVec.h"

//FIX -- use hash_map
//FIX -- hash_map already keeps <key,data> pairs. keeping data with key in them 
//wastes space.

//map<unsigned int, void *> mymap;

using namespace std;
using namespace nvc;

const int KNNDIM = 3;

class UniformGrid
{
public:
    typedef unsigned int index_t;
    typedef float coord_t;
    typedef NVec<3,index_t> Index;
    typedef NVec<3,coord_t> Vec;

    Index size;
    Point3D bb_min, bb_max, cellsize;

    UniformGrid(index_t w, index_t h, index_t d) { size = Index(w,h,d); }

    void set_extent(const Point3D & m, const Point3D & M)
    {
	bb_min = m;  bb_max = M;
	cellsize = M - m;
	for(int k=0; k<KNNDIM; k++)
	{
	    cellsize[k] /= (float)(size[k]);
	}
    }
  
  void get_cellsize(Point3D & c){c =  cellsize;};

  UniformGrid::Index compute_index(const Point3D & v)
  {
	UniformGrid::Index i;
	for(int k=0; k<KNNDIM; k++)
	  {
	  int x = (int)floor(((v[k] - bb_min[k]) / cellsize[k]));
	  if (x <0)
		i[k] = 0;
	  else
		i[k] = x;
	  if( i[k] >= size[k])
		{ 
		  i[k] = size[k]-1;
		}
	  }
	return i;
    }
};


template<class CellData>
class SparseGrid : public UniformGrid
{
    protected:
    typedef map<unsigned int, CellData *> Map;
  
    //
    // ASSUME: Each component of the index must be AT MOST 10 bits (0..2047).
    //         We're assuming that this should be sufficient for most
    //         reasonable uses for which we'll be using this grid.
    //
    unsigned int pack_index(const UniformGrid::Index& i)
    	{ return (i[0]<<20) | (i[1]<<10) | i[2]; }

    

    Map cellmap;
    
public:
    SparseGrid(int w, int h, int d) : UniformGrid(w,h,d) { }
    ~SparseGrid()
    {
	// Delete all the entries that were allocated by locate()
	for(typename Map::iterator i=cellmap.begin(); i!=cellmap.end(); ++i)
	    delete i->second;
    }

    //
    // Iteration interface just reflects definitions from Map
    //
    typedef typename Map::iterator iterator;
    typedef typename Map::const_iterator const_iterator;

    iterator begin() { return cellmap.begin(); }
    iterator end() { return cellmap.end(); }
    CellData *iterator_data(iterator i) { return i->second; }

    int usage_count() const { return cellmap.size(); }
 

    //
    // Accessors become lookups in the internal Map
    //

  CellData *locate_index_only(const UniformGrid::Index& i)
  {
    iterator iter = cellmap.find(pack_index(i));
    return iter==cellmap.end() ? NULL : iter->second;
  }
  
  CellData *locate_packed_index_only(unsigned int i)
  {
    iterator iter = cellmap.find(i);
    return iter==cellmap.end() ? NULL : iter->second;
  }

    CellData *locate_only(const Point3D& v)
    	{ return locate_index_only(compute_index(v)); }

   CellData *locate(const Point3D& v)
    {
	UniformGrid::Index i = compute_index(v);
	CellData *d = locate_index(i);//locate_index_only(i);
	if( !d )
	  {
	    cellmap[pack_index(i)] = d = new CellData;
	  }
	return d;
    }

  CellData *locate_index(const UniformGrid::Index& i)
  {
    CellData *d = locate_index_only(i);
    if( !d )
      {
	cellmap[pack_index(i)] = d = new CellData;
      }
    return d;
  }

  CellData *locate_packed_index(const unsigned int i)
  {
    CellData *d = locate_packed_index_only(i);
    if( !d )
      {
	cellmap[i] = d = new CellData;
      }
    return d;
  }
};


template <typename T>
class KNN_Grid
{  

  //----------------------------------------------------- 
  //Helper classes
    class KNNnbr
    {
    public:
      float distance;
      T  item;

      KNNnbr(){};
      bool operator<(const KNNnbr & p)const
      {
	return (distance < p.distance);
      }

      KNNnbr & operator = (const KNNnbr & p){ 
 	distance = p.distance; item = p.item; 
 	return *this;
      }

    };

    typedef vector< T > SpaceCell;
	typedef unsigned int index_t;
	typedef NVec<3,index_t> Index;
  //-------------------------------------------------------


  //private data
  Point3D bbmin;
  Point3D bbmax;
  int     cell_dim[3];
  float   cellsz[3];

  SparseGrid< SpaceCell > * grid;
  
public:
  KNN_Grid(Point3D m, Point3D M, int long_dim)
    :bbmin(m),bbmax(M)
  {
    //Construct a grid with square elements
    Point3D box_len = bbmax-bbmin;
    int max_dim;
    // Find max length dimension of the grid
    if ( box_len[0] > box_len[1])
      max_dim = 0;
    else
      max_dim = 1;
    if (box_len[2]>box_len[max_dim])
      max_dim = 2;
    
    cellsz[max_dim] = box_len[max_dim]/((float) long_dim);
    for(int i=0; i<3;i++){
      cellsz[i] = cellsz[max_dim];
      cell_dim[i] = (int)ceil(box_len[0]/cellsz[max_dim]); // ** maybe box_len[i] **
      bbmax[i] = bbmin[i] + ((float) cell_dim[i]) * cellsz[i];  
    }
    grid = new SparseGrid< SpaceCell >(cell_dim[0],cell_dim[1],cell_dim[2]);
    grid->set_extent(bbmin,bbmax);	
	//cout << "Grid Extent \n";	
	
  }

  KNN_Grid(Point3D m, Point3D M, int x_dim, int y_dim, int z_dim)
    :bbmin(m),bbmax(M)
  {
    cell_dim[0]=x_dim;cell_dim[1]=y_dim;cell_dim[2]=z_dim;
    grid = new SparseGrid< SpaceCell >(cell_dim[0],cell_dim[1],cell_dim[2]);
    grid->set_extent(bbmin,bbmax); 
  }

  int insert(T & item, Point3D &p)
  {
    SpaceCell * cell = grid->locate(p);
    cell->push_back(item);
	//cout << "Size: "<< cell->size() << "\n" ;
    return cell->size();
  }

  int insertRange(T & item, Point3D &p, Point3D &q, Point3D &r){
	int count = 0;
	
	Index min,max;
	
	Index a = grid->compute_index(p);
	Index b = grid->compute_index(q);
	Index c = grid->compute_index(r);
	
	for(int i=0;i<3;++i)
	{
		min[i] = (a[i]>b[i])?b[i]:a[i];
		max[i] = (a[i]>b[i])?a[i]:b[i];
	}

	for(int i=0;i<3;++i)
	{
		min[i] = (min[i]>c[i])?c[i]:min[i];
		max[i] = (max[i]>c[i])?max[i]:c[i];
	}

	for(int i=min[0]; i<=max[0]; ++i)
	{
		for(int j=min[1]; j<=max[1]; ++j)
		{
			for(int k=min[2]; k<=max[2]; ++k)
			{
				Index x(i,j,k);
				SpaceCell * cell = grid->locate_index(x);
				typename vector<T>::iterator temp ;
				//temp = find(cell->begin(), cell->end(), item);
				//if(temp == cell->end() || cell->size() == 0)
				//{
					cell->push_back(item);
					count ++;//= cell->size();
				//}
			}
		}
	}
	
	return count;
  }
  
  void clear(){
    delete grid;
    grid = new SparseGrid< SpaceCell >(cell_dim[0],cell_dim[1],cell_dim[2]);
    grid->set_extent(bbmin,bbmax);
  }

  unsigned int k_nearest(Point3D &p, int k, vector<T> & nbrs){
    KNNnbr pnew;
    priority_queue<KNNnbr> pts;
    vector< T > * cell_pts = NULL;
 
    UniformGrid::Index grid_coord = grid->compute_index(p);
    for(int x = grid_coord[0] - 2; x<grid_coord[0]+2;x++)
      for(int y = grid_coord[1] - 2; y<grid_coord[1]+2;y++)
	for(int z = grid_coord[2] - 2; z<grid_coord[2]+2;z++)
	  {
	    //get cell points
	    UniformGrid::Index i(x,y,z); 
	    //cerr << "Checking cell: [" << x <<  " , " << y << " , " 
	    // << z << "]" << endl; 
	    cell_pts = get_cell(i);
	    if (cell_pts != NULL)
	      {
		for(int j=0;j<cell_pts->size();j++)
		  {
		    pnew.distance = euclid_distance((*cell_pts)[j].coord,p);
		    pnew.item = (*cell_pts)[j];
		    pts.push(pnew);
		// insert candidates into queue
		  }
	      }
	  }
    int j=0;
    int ptnum = pts.size();
    nbrs.clear();
    while( (j<k) && (!pts.empty()))
      {
	for(int j=0; j<k;j++)
	  {
	    KNNnbr n = pts.top();
	    nbrs.push_back(n.item);
	    pts.pop();
	  }
      }
    return ptnum;
  }
 
vector<T> * get_cell_range(Point3D &p, Point3D &q){
	vector<T> * ret = new vector<T>;
	vector<T> * temp = NULL;
	
	Index a = grid->compute_index(p);
	Index b = grid->compute_index(q);
	// cout << "A: "<< a[0]<< ", "<< a[1]<< ", "<< a[2]<< "\n";
	// cout << "B: "<< b[0]<< ", "<< b[1]<< ", "<< b[2]<< "\n";
	Index min,max;
	
	for(int i=0;i<3;++i)
	{
		min[i] = (a[i]>b[i])?b[i]:a[i];
		max[i] = (a[i]>b[i])?a[i]:b[i];
	}

	for(int i=min[0]; i<=max[0];++i)
	{
		for(int j=min[1]; j<=max[1]; ++j)
		{
			for(int k=min[2]; k<=max[2]; ++k)
			{
				Index p(i,j,k);
				//cout<<"Index x: "<< i <<" y:"<< j <<" z:"<< k <<"\n";
				temp = get_cell(p);
				if(temp!=NULL){
					//cout<< "temp size: " << temp->size() <<"\n";
					ret->insert(ret->begin(), temp->begin(), temp->end());
				}
			}
		}
	}
	return ret;
  }
  
  vector<T> * get_cell(Point3D &p){
    return grid->locate_only(p);
  }
	
  vector<T> * get_cell(UniformGrid::Index &p){
    for(int i=0;i<3;i++)
		if ((p[i] < 0) || (p[i] >= cell_dim[i]))
			return NULL;
	//cerr << "i get in here\n";
    return grid->locate_index_only(p);
  }
  
};

  
#endif
