/****************************************************************************************/
/*                                                                                      */
/*                                      vinci.h                                         */
/*                                                                                      */
/****************************************************************************************/
/*                                                                                      */
/* Authors: Benno Bueeler (bueeler@ifor.math.ethz.ch)                                   */
/*          and                                                                         */
/*          Andreas Enge (enge@ifor.math.ethz.ch)                                       */
/*          Institute for Operations Research                                           */
/*	    Swiss Federal Institute of Technology Zurich                                */
/*	    Switzerland                                                                 */
/*                                                                                      */
/* Last Changes: April 7, 1998                                                          */
/*                                                                                      */
/****************************************************************************************/
/*                                                                                      */
/* global constants, types, variable and function declarations                          */
/*                                                                                      */
/****************************************************************************************/

#define COPYRIGHT1 "(C) Benno Bueeler and Andreas Enge, {bueeler,enge}@ifor.math.ethz.ch"
#define COPYRIGHT2 "    Developed in the ETHZ-EPFL project on Optimisation and Computational"
#define COPYRIGHT3 "    Geometry with Komei Fukuda"
#define VERSION "1"
#define VERSION_DATE "March 1998"

#define T01 "The package computes the volume of a polytope whose vertices, defining"
#define T02 "hyperplanes and/or incidences of vertices and facettes are stored in files"
#define T03 "following Avis' and Fukuda's polytope format. The vertices are supposed to"
#define T04 "be in a file with extension '.ine', the hyperplanes in '.ext' and the"
#define T05 "incidences in '.icd' (see the sample files 'square.ext', 'square.icd' and"
#define T06 "'square.ine').\n"
#define T07 "Its basic call is 'vinci file' where 'file' stands for the polyhedron file"
#define T08 "name without extension, e. g. 'vinci square'. In this case the existing files"
#define T09 "and installed additional packages are analysed and according to the problem"
#define T10 "type an appropriate volume computation method is chosen automatically."
#define T11 ""
#define T12 "The following additional parameters are possible:"
#define T13 "-m followed by the label of a specific volume computation method; to see a"
#define T14 "   list of the implemented methods, use 'vinci -m' without any label."
#define T15 "-p1, -p2 or -p3 which are valid for some methods. They determine the space"
#define T16 "   and time complexity. Higher values lead to more space consume and a fast-"
#define T17 "   er execution. If none of the parameters is specified the value"
#define T18 "   DEFAULT_PRECOMP defined in 'vinci.h' is chosen."
#define T19 "-s directly followed by a natural integer. The value determines for how many"
#define T20 "   recursion levels intermediate results are stored. A higher value speeds"
#define T21 "   up certain methods considerably while needing more storage space."
#define T22 "-r directly followed by an integer. The value sets the random seed used for"
#define T23 "   determining the objective function for Lawrence's formula."
#define T24 "\nFor more information please consult the manual."

/****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/****************************************************************************************/
/* The following constants may be changed by a user                                     */
/****************************************************************************************/

#define LRS_EXEC      "lrs"
#define CDD_EXEC      "cddf+"
#define QHULL_EXEC    "qhull"
#define QHULL_OPTIONS "d i Q0 Qz"
   /* location of external programmes with paths and options for qhull                  */
   
#define PIVOTING 1
#define MIN_PIVOT 0.5
#define PIVOTING_LASS 0
#define MIN_PIVOT_LASS 0.1
   /* for choosing a pivoting strategy whenever this is needed, e.g. for determinant    */
   /* computation. The last two constants are valid for Lasserre's method, the first    */
   /* two for any other method. PIVOTING and PIVOTING_LASS can take the following       */
   /* values:                                                                           */
   /*  0: The first row with entry bigger than MIN_PIVOT (in its absolute value) is     */
   /*     chosen.                                                                       */
   /*     bigger than EPSILON or smaller than -EPSILON) is chosen                       */
   /*  1: partial pivoting; the row with maximal entry (absolute value) is chosen;      */
   /*     MIN_PIVOT is ignored.                                                         */
   /*  2: total pivoting; the row and column with maximal entry (absolute value) are    */
   /*     chosen, MIN_PIVOT is ignored. If total pivoting is not possible (e. g. for    */
   /*     solving linear equations) partial pivoting is performed.                      */
   
#define MAX_PRECOMP 3
#define DEFAULT_PRECOMP 1
   /* The global variable "G_Precomp" specifies the precomputation effort of the pro-   */
   /* gram; it can be set to 1, 2 or 3 by using the options "-p1", "-p2", "-p3" resp.   */
   /*  1: Each intersection is done with the hyperplane itself.                         */
   /*  2: The intersections of all possible choices of two hyperplanes are computed     */
   /*     once. Then instead of intersecting with one hyperplane, intersection is done  */
   /*     with two hyperplanes at the same time: The new one and another one which com- */
   /*     pletely contains the set. So the same result is obtained, but intersection is */
   /*     faster because the sets are smaller.                                          */
   /*  3: The same as 2, but with three hyperplanes.                                    */
   /* MAX_PRECOMP is the maximum value the variable can take, DEFAULT_PRECOMP the value */
   /* used when nothing is specified by the user.                                       */

#define DEFAULT_STORAGE 20
   /* For so many recursion levels intermediate volumes are stored, starting from the   */
   /* level with two planes fixed. The value is only active if no option -s is speci-   */
   /* fied; it may be set to 0, for instance, by using the option -s0.                  */
   
#define STATISTICS 
   /* If STATISTICS is defined, during volume computation, some statistical variables   */
   /* like the number of simplices and their volume distribution are withheld.          */

#define NO_RATIONAL
   /* If RATIONAL is defined, exact computation will be used in some places. In this    */
   /* case, some real variables will be defined as rational. This option only works if  */
   /* the g++-compiler is used, and its exact consequences are unknown...               */
   
#define NO_COG
   /* If COG is defined, the programme makes functions available which compute the      */
   /* center of gravity as well as the volume. Their names start with                   */
   /* "center_of_gravity". Otherwise basically the same code is used for computing only */
   /* the volume. The function parameters are still the same (when dropping the         */
   /* variable supposed to hold the center of gravity), but their names now start with  */
   /* "volume". If COG is undefined, some volume computation functions using Lawrence's */
   /* formula or Lasserre's method are added. These functions are not appropriate for   */
   /* computing the center of gravity.                                                  */
   
/****************************************************************************************/
/* The following constants should not be changed by a user                              */
/****************************************************************************************/

#define	EPSILON 1e-10
#define ARRAYSIZESTEP 5
                          /* the dynamic arrays are increased in steps of this size */

#define FALSE 0
#define TRUE 1
#define EOLN 10           /* end of line character */

#define NONE       0      /* constants for data types in files */
#define INTEGER_T  1
#define REAL_T     2
#define RATIONAL_T 3

#define RCH    1          /* constants for the volume computation methods */
#define GS     2
#define GSH    3
#define HOT    4
#define HOTH   5
#define CDDP   6
#define CDDM   7
#define QHULL  8
#define LAWD   9
#define LAWND 10
#define RLASS 11
#define LRS   12

#define PREPARE  1        /* constants for precomputation of hyperplane intersections   */
#define RETRIEVE 2
#define FREE     3

#define KEY_PLANES     1  /* constants for the key type actually used in the balanced   */
#define KEY_VERTICES   2  /* tree routines                                              */
#define KEY_PLANES_VAR 3

#define SCHMIDT     1     /* constants for the orthonormalisation technique used */
#define HOUSEHOLDER 2

#ifdef STATISTICS
#define STAT_SMALLEST_EXP -200
#define STAT_BIGGEST_EXP 200
#endif
   /* The simplices in the size classes from 1e(STAT_SMALLEST_EXP) to                   */
   /* 1e(STAT_BIGGEST_EXP) are counted separately for each class. Simplices too small   */
   /* or too big are summarized in one variable. */
   
#ifdef RATIONAL
#include <Rational.h>
#endif

/****************************************************************************************/
/*                                  type definitions                                    */
/****************************************************************************************/

typedef double        real;
typedef int           boolean;
typedef unsigned char T_LassInt;

#ifdef RATIONAL
   typedef Rational rational;
#else
   typedef real     rational;
#endif

struct T_Vertex
     {real *coords;
          /* the coordinates of a vertex */
      long int no;           
     }; /* type of a vertex */
typedef struct T_Vertex T_Vertex;

/* The sets of vertices are implemented as ordered lists */
struct T_VertexSet
       {int      maxel;  /* The list of vertices loe may contain elements 0 to maxel    */
        int      lastel; /* Effectively, loe contains elements from 0 to lastel;        */
                         /* an empty set is indicated by lastel == -1.                  */
        T_Vertex **loe;  /* The elements of the set are pointers to vertices which are  */
                         /* stored in loe, in ascending order following the numbering   */
                         /* of the vertices */
       };
typedef struct T_VertexSet T_VertexSet;


/* The sets of VertexSets are implemented as linked chain because only simple           */
/* operations are needed on them. */
struct T_VertexSuperset
       {struct T_VertexSuperset *next;
        struct T_VertexSet       content;
       };
typedef struct T_VertexSuperset T_VertexSuperset;


/* Types for storing face volumes in balanced trees. For storing and retrieving a key   */
/* is needed for every face; all possibilities for keys are defined by the union T_Key. */
/* Only one of the keys is active at a time. T_Tree defines the tree itself.            */
union T_Key
      {int *hyperplanes;
          /* The indices of the hyperplanes which are intersected to obtain the face,   */
          /* in ascending order. The list is terminated by G_m + 1.                     */
       struct {T_VertexSet set;
               int         d;
              } vertices;
          /* The vertices contained in the face and the dimension for which the volume  */
          /* is stored.                                                                 */
       struct {T_LassInt *hyperplanes, *variables;} hypervar;
          /* A list of fixed constraints and variables onto which the face has been     */
          /* projected, used for Lasserre's formula                                     */
      };
typedef union T_Key T_Key;

struct T_Tree
       {struct T_Tree *tree_l, *tree_r; /* the left and right subtrees */
        int           tree_b;
        T_Key         key;
        rational      vol;              /* the stored volume */
       };
typedef struct T_Tree T_Tree;
        
/****************************************************************************************/
/*                            global variable declarations                              */
/****************************************************************************************/

extern int G_d;
extern int G_m;
extern real **G_Hyperplanes;
   /* The first d components are the coordinates of a normal vector, the d+first com-   */
   /* ponent is the right hand side. */
extern T_VertexSet *G_Incidence;
   /* sets of the vertices contained in the hyperplanes. The indices are the same as in */
   /* 'G_Hyperplanes'. */
extern T_VertexSet G_Vertices;
   /* set of the vertices of the polytope                                               */
extern int G_Precomp;
   /* see the annotations for MAX_PRECOMP and DEFAULT_PRECOMP                           */
extern int G_Storage;
   /* see the annotations for DEFAULT_STORAGE                                           */
extern int G_RandomSeed;
extern boolean G_OutOfMem;
extern void *G_MemRes;
   /* Memory reserve; in the storage routine a certain part of the memory is blocked by */
   /* allocating it to this variable. If the programme runs out of memory this part is  */
   /* freed and G_OutOfMem is set to TRUE.                                              */
extern rational G_Minus1;
   /* dummy to have a pointer, namely &G_Minus1, to the value -1 */

#ifdef STATISTICS
   /* some variables for statistical purpose */
   extern unsigned int Stat_Count;
      /* for counting the number of partial volumes computed */
   extern real Stat_Smallest, Stat_Biggest; 
      /* the volumes of the smallest and the biggest simplex */
   extern unsigned int Stat_CountPos [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 3]; 
      /* Stat_CountPos [i] counts the number of positive partial volumes (usually       */
      /* simplex volumes except for Lawrence's formula) which are contained between     */
      /* 1e(i - 1 + STAT_SMALLEST_EXP) and 1e(i + STAT_SMALLEST_EXP). Smaller volumes   */
      /* are counted in Stat_CountPos [0], bigger ones in the last entry of             */
      /* Stat_CountPos. */
   extern unsigned int Stat_CountNeg [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 3];
      /* the same for the absolute values of "negative volumes" in Lawrence's formula   */
   extern unsigned int *Stat_CountStored, *Stat_CountRetrieved;
      /* counts the number of volumes stored in and retrieved from the tree             */
   extern unsigned int Stat_CountShifts;
      /* counts the number of shifts performed in Lasserre's method                     */
   extern long int Stat_ActualMem;
      /* the memory actually used on the heap                                           */
   extern long int Stat_MaxMem;
      /* the maximal heap memory used during execution                                  */
      /* Both values do not comprise the memory reserve G_MemRes.                       */
#endif

/****************************************************************************************/
/*                    functions and procedures from 'vinci_memory'                      */
/****************************************************************************************/

void *my_malloc (long int size);
void *my_realloc (void *pointer, long int new_size, long int size_diff);
void my_free (void *pointer, long int size);
T_Vertex *create_vertex ();
void free_vertex (T_Vertex *v);
void create_hyperplanes ();
void free_hyperplanes ();
void create_incidence ();
void free_incidence ();
T_VertexSet *create_faces ();
void free_faces (T_VertexSet *face);
rational ***create_basis ();
void free_basis (rational *** basis);
rational *create_fact ();
void free_fact (rational *fact);
rational *create_vector ();
void free_vector (rational *v);
int *create_int_vector (int n);
void free_int_vector (int *v, int n);
rational **create_matrix (int m, int n);
void redim_matrix (rational ***A, int m_alt, int m_neu, int n);
void free_matrix (rational **A, int m, int n);
void create_key (T_Key *key, int key_choice);
void free_key (T_Key key, int key_choice);

void tree_out (T_Tree **ppr , int *pi_balance, T_Key key, rational **volume,
   T_Key **keyfound, int key_choice);
void add_hypervar (T_LassInt hyperplane, T_LassInt variable, T_Key *key);
void delete_hypervar (T_LassInt hyperplane, T_LassInt variable, T_Key *key);

/****************************************************************************************/
/*                      functions and procedures from 'vinci_set'                       */
/****************************************************************************************/

T_VertexSet create_empty_set (void);
T_VertexSet duplicate_set (T_VertexSet s);
rational normalise_vertices ();
void print_set (FILE *f, T_VertexSet s);
boolean empty (T_VertexSet s);
boolean is_in_set (T_Vertex *e, T_VertexSet s);
boolean is_contained (T_VertexSet s1, T_VertexSet s2);
boolean are_equal_sets (T_VertexSet s1, T_VertexSet s2);
void add_element (T_VertexSet *s, T_Vertex *e);
boolean delete_element (T_VertexSet *s, T_Vertex *e);
void intersect (T_VertexSet s1, T_VertexSet s2, T_VertexSet *inter);
T_VertexSuperset *create_empty_superset (void);
void free_superset (T_VertexSuperset **S);
void print_superset (FILE *f, T_VertexSuperset *S);
boolean is_in_superset (T_VertexSet s, T_VertexSuperset *S);
void add_superelement (T_VertexSuperset **S, T_VertexSet s);

/****************************************************************************************/
/*                  functions and procedures from 'vinci_computation'                   */
/****************************************************************************************/

rational factorial (int n);
#ifdef COG
void simplex_center (T_VertexSet S, real *center);
#endif
rational det_and_invert (rational **A, int rows, int columns, boolean verbose);
void simplex_volume (T_VertexSet S, rational *volume, boolean verbose);
T_VertexSet* do_intersections (int precomp, int i1, int i2, int i3, int action);
rational add_orthonormal_schmidt (int d, T_VertexSet face, rational **B,
   T_Vertex *vertex);
rational orthonormal_schmidt (int d, T_VertexSet face, rational **B);
rational add_orthonormal_householder (int d, T_VertexSet face, rational **H,
   T_Vertex *vertex);
rational orthonormal_householder (int d, T_VertexSet face, rational **H);

/****************************************************************************************/
/*                     functions and procedures from 'vinci_volume'                     */
/****************************************************************************************/

#ifdef COG
void center_of_gravity_ch_file (rational *volume, real *center,
   char vertexfile [255], char incidencefile [255]);
void center_of_gravity_cdd_file (rational *volume, real *center, int direction,
   char vertexfile [255]);
void center_of_gravity_qhull_file (rational *volume, real *center, char vertexfile [255]);
#else
void volume_ch_file (rational *volume, char vertexfile [255], char incidencefile [255]);
void volume_ortho_file (rational *volume, char vertexfile [255], char incidencefile [255],
   int key_choice, int ortho_choice);
void volume_cdd_file (rational *volume, int direction, char vertexfile [255]); 
void volume_qhull_file (rational *volume, char vertexfile [255]);
void volume_lawrence_file (rational *volume, char vertexfile [255], char planesfile [255],
                           char incidencefile [255]);
void volume_lawrence_lrs_file (rational *volume, char planesfile [255]);
void volume_lrs_file (rational *volume, char *rational_volume, char vertexfile [255]);
#endif

/****************************************************************************************/
/*                 functions and procedures from 'vinci_lass'                           */
/****************************************************************************************/

#ifndef COG
void volume_lasserre_file (double planes2[7][4], rational *volume);
#endif

