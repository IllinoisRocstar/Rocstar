/****************************************************************************************/
/*                                                                                      */
/*                              vinci_lass.c                                            */
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
/* Last Changes: October 18, 1997                                                       */
/*                                                                                      */
/****************************************************************************************/
/*                                                                                      */
/* Lasserre's volume computation method                                                 */
/*                                                                                      */
/****************************************************************************************/

#include "vinci.h"


#define MAXIMUM 1.0e150   /* define the maximum used for testing infinity */
#define	EPSILON_LASS EPSILON   /* Numbers smaller than this are treated as zero in rhs*/
#define	EPS1    EPSILON    /* Numbers smaller than this are treated as zero in coefficient */
#define	EPS_NORM EPSILON   /* EPSILON used for constraint normalization*/
#define LaShiftLevel 0    /* Shifting is possible if d>=LaShiftLevel */
#define LaShift 1         /* shift polytope to make at least d components of
                             the rhs zero if there are less than d-LaShift zeros */
/* #define ReverseLass */ /* perform recursion from last to first constraint;
                             if undefined the recursion starts with the first constraint */
#define verboseFirstLevel /* output the intermediate volume of the first level */


/******************/
/*global variables*/
/******************/

/****************************************************************************************/
/*                                                                                      */
/* definitions of the global variables                                                  */
/*                                                                                      */
/****************************************************************************************/


/****************************************************************************************/

int G_d;
int G_m;

real **G_Hyperplanes = NULL; 
T_VertexSet *G_Incidence = NULL;
T_VertexSet G_Vertices;

int G_Precomp = 0;
int G_Storage = -1;
int G_RandomSeed = 4;

boolean G_OutOfMem = TRUE;
void *G_MemRes = NULL;
rational G_Minus1 = -1;

/****************************************************************************************/

rational *A;
rational *pivotrow;             /* copy of pivot row */
T_LassInt *All_index;  /* All eliminated and superfluous indices (sorted) */
T_LassInt *Pivot;      /* All substituted variables (sorted) */
int **p2c;        /* pivot to constraints: which variable is fixed in which constraint;
                     the variable index is given in the leading column, the constraint 
		     index in the second */

rational  * planescopy;               /* needed in shift_P in the lasserre-code */
static T_Key   key, *keyfound;    /* key for storing the actually considered face and   */
                                  /* found key when a volume could be retrieved */
static T_Tree  *tree_volumes;     /* tree for storing intermediate volumes */
static T_VertexSet *face;
   /* face considered at each recursion level */

/***************/
/*help routines*/
/***************/

/****************************************************************************************/
/*                                 memory allocation                                    */
/****************************************************************************************/

void *my_malloc (long int size)
   /* works exactly like malloc, but keeps trace of the used memory in the statistical  */
   /* variables and eventually frees the memory reserve                                 */
   
{  void *pointer;

   pointer = malloc (size);

   if (pointer == NULL)
      if (!G_OutOfMem)
      {
         fprintf (stderr, "\n***** WARNING: Out of memory; freeing memory reserve\n");
         free (G_MemRes);
         G_OutOfMem = TRUE;
         pointer = my_malloc (size);
      }
      else
      {
         fprintf (stderr, "\n***** ERROR: Out of memory\n");
         exit (0);
      }
   else
   {

   }

   return pointer;
}

/****************************************************************************************/

void my_free (void *pointer, long int size)
   /* frees the memory space used by pointer and keeps track of the freed space in the  */
   /* statistical variables                                                             */

{
   free (pointer);
}

/****************************************************************************************/

void create_key (T_Key *key, int key_choice)
   /* creates the dynamic parts of the key; G_d must be set correctly */

{
   switch (key_choice)
   {
   case KEY_PLANES:
      key -> hyperplanes = create_int_vector (G_d + 1);
      break;
   case KEY_PLANES_VAR:
      key -> hypervar.hyperplanes =
                                 (T_LassInt *) my_malloc ((G_Storage + 2) * sizeof (T_LassInt));
      key -> hypervar.variables =
                                 (T_LassInt *) my_malloc ((G_Storage + 2) * sizeof (T_LassInt));
      break;
   }
}

/****************************************************************************************/

void free_key (T_Key key, int key_choice)
   /* frees the dynamic parts of the key */

{
   switch (key_choice)
   {
   case KEY_PLANES:
      free_int_vector (key.hyperplanes, G_d + 1);
      break;
   case KEY_PLANES_VAR:
      my_free (key.hypervar.hyperplanes, (G_Storage + 2) * sizeof (T_LassInt));
      my_free (key.hypervar.variables, (G_Storage + 2) * sizeof (T_LassInt));
      break;
   }
}
/****************************************************************************************/

void add_hypervar (T_LassInt hyperplane, T_LassInt variable, T_Key *key)
   /* adds the specified hyperplane and variable index to the variable "key" maintain-  */
   /* ing the ascending orders; if one index is G_m+1 resp. G_d+1 it is omitted.        */
   /* For the working of the procedure it is necessary that the last array entry is     */
   /* G_m + 1 resp. G_d + 1 and that there is still some space left in the arrays.      */
   /* Attention: Only use this function if you work with the planes and variables as    */
   /* key!                                                                              */
   
{  int i, j;

   if (hyperplane != G_m+1)
   {  
      for (i = 0; (*key).hypervar.hyperplanes [i] < hyperplane; i++);
      if ((*key).hypervar.hyperplanes [i] != hyperplane)
      {  /* insert index */
         for (j = G_d; j > i; j--)
            (*key).hypervar.hyperplanes [j] = (*key).hypervar.hyperplanes [j-1];
         (*key).hypervar.hyperplanes [i] = hyperplane;
      }
   }

   if (variable != G_d+1)
   {  
      for (i = 0; (*key).hypervar.variables [i] < variable; i++);
      if ((*key).hypervar.variables [i] != variable)
      {  /* insert index */
         for (j = G_d; j > i; j--)
            (*key).hypervar.variables [j] = (*key).hypervar.variables [j-1];
         (*key).hypervar.variables [i] = variable;
      }
   }
}


/****************************************************************************************/

void delete_hypervar (T_LassInt hyperplane, T_LassInt variable, T_Key *key)
   /* deletes the indices from the variable key; if one index is -1 it is omitted.      */
   /* Attention: Only use this function if you work with the planes and variables as    */
   /* key!                                                                              */
   
{  int i, j;

   if (hyperplane != G_m+1)
   {  
      for (i = 0; i <= G_d && (*key).hypervar.hyperplanes [i] != hyperplane; i++);
      if (i != G_d + 1)
      {  /* index found, delete it */
         for (j = i; (*key).hypervar.hyperplanes [j] != G_m + 1; j++)
            (*key).hypervar.hyperplanes [j] = (*key).hypervar.hyperplanes [j+1];
      }
   }

   if (variable != G_d+1)
   {  
      for (i = 0; i <= G_d && (*key).hypervar.variables [i] != variable; i++);
      if (i != G_d + 1)
      {  /* index found, delete it */
         for (j = i; (*key).hypervar.variables [j] != G_d + 1; j++)
            (*key).hypervar.variables [j] = (*key).hypervar.variables [j+1];
      }
   }
}

/****************************************************************************************/

/****************************************************************************************/

int *create_int_vector (n)
   /* reserves memory space for a vector with n entries                                 */
   
{  int *v;
   
   v = (int *) my_malloc (n * sizeof (int));
   
   return v;
}

/****************************************************************************************/

void free_int_vector (int *v, int n)
   /* frees the memory space needed by the vector v of length n                         */

{  
   my_free (v, n * sizeof (int));
}

/****************************************************************************************/


static T_LassInt add_reduced_index(T_LassInt red, T_LassInt * indices, 
                                 T_LassInt * ref_indices)
/* insert new index into ref_indices maintaining sorting; if indices!=NULL this index
   is also inserted into indices. returns the original index base.
   assumption: red counted in the original systen is not contained in ref_indices */

{ register int i;
  T_LassInt xch, base;
    
  for (i=0; red>=ref_indices[i]; i++) red++;  /* reduced index -> original index */
  base=red;  
  while (ref_indices[i]<=G_m) {
      xch=ref_indices[i];
      ref_indices[i]=red;
      red=xch;
      i++;
  };
  ref_indices[i]=red;
  ref_indices[i+1]=G_m+2;
  if (indices==NULL) return base;
  red=base;
  for (i=0; base>indices[i]; i++);
  while (indices[i]<=G_m) {
      xch=indices[i];
      indices[i]=red;
      red=xch;
      i++;
  };
  indices[i]=red;
  indices[i+1]=G_m+2;
  return base;
}


static void del_original_indices(T_LassInt *indices, T_LassInt *org_indices)
/* delete original indices in org_indices maintaining sorting.
   assumption: all the indices are contained in org_indices.
   the end of indices is marked by G_m+2 */

{   register int i, cnt;

    i=cnt=0;
    while (org_indices[i]<=G_m) {
	while ((org_indices[cnt+i]==indices[cnt])&&(indices[cnt]<=G_m)) cnt++;
	org_indices[i]=org_indices[i+cnt];
	i++;
    }
}


static void del_original(T_LassInt base, T_LassInt * indices)
/* delete base in indices maintaining sorting.
   assumption: base is contained in indices. */

{ int i;

  for (i=0; ((base!=indices[i])&&(indices[i]<=G_m)); i++);  /* search original index */
  if (base!=indices[i]) {
      fprintf(stderr, "ERROR: Deletion index not found!\n");
      exit(0);
  };
  for (;indices[i]<=G_m;i++) indices[i]=indices[i+1];
}  


static void rm_original_inElAll_index(T_LassInt baserow)
/* delete baserow in All_index maintaining sorting. */
{   del_original(baserow, All_index); }


static void rm_constraint(rational* A, int *LastPlane_, int d, int rm_index)
/* removes the constraints given in rm_index and adjusts *LastPlane */

{   register rational *p1, *p2; 
    register int i;

    p1=A+rm_index*(d+1);
    p2=A+(rm_index+1)*(d+1);
    for (i=0; i<(((*LastPlane_)-rm_index)*(d+1)); i++) {
	*p1=*p2;
	p1++;
	p2++;
    };
    (*LastPlane_)--;
}


static rational * compact()
{   register int i, j;
    register rational *po,*pc;

    if (!(pc = (rational *) my_malloc (G_m*(G_d+1)*sizeof(rational)))) {
	fprintf (stderr, "\n***** ERROR: Out of memory in 'compact.*pc'");
	exit(0); 
    }
    po=pc;
    for (i=0; i<G_m; i++) {
	for (j=0; j<=G_d; j++,pc++) *pc= G_Hyperplanes [i][j];
    };
    return po;
}


/***************/
/*Core routines*/
/***************/


static int notInPivot(int * pivot, int col, int i)
{ register int h;
  for (h=0;h<col;h++)
   if (pivot[h]==i) return FALSE;
  return TRUE;
}


static void shift_P(rational *A, int LastPlane_, int d)
/*  shift one vertex of the polytope into the origin, that
    is, make at least d components of b equal zero */

{   register rational  *p1, *p2, *p3, d1, d2, d3;
    register int col, i, j;
    static int *pivot = NULL;
                 /* contains the pivot row of each column */



    if (pivot == NULL) pivot = create_int_vector (G_d + 1);
    
    p1=A;                         /* search pivot of first column */
    pivot[0]=0; 
    d3=fabs(d1=*p1);
    for (i=0; i<=LastPlane_; i++) {
        d2=fabs(*p1);
#if PIVOTING_LASS == 0
	if (d2>=MIN_PIVOT_LASS) {pivot[0]=i; d1=*p1; break;};
#endif
	if (d2>d3) { pivot[0]=i; d1=*p1; d3=d2; };
	p1+=(d+1);
    }
    /* copy pivot row into planescopy */
    p1=A+pivot[0]*(d+1)+1;   
    p2=planescopy+pivot[0]*(d+1)+1;
    for (i=1,d2=1.0/d1; i<=d; i++,p1++,p2++) *p2 = (*p1)*d2;
    /* complete first pivoting and copying */
    p1=A+1;                          
    p2=planescopy+1;
    for (i=0; i<=LastPlane_; i++, p1++, p2++) {
	if (i==pivot[0]) {
	    p1+=d;
	    p2+=d;
	    continue;   /* pivot row already done */
	}
	d1=*(p1-1); 
	p3=planescopy+pivot[0]*(d+1)+1;
	for (j=1; j<=d; j++, p1++, p2++, p3++) (*p2)=(*p1)-d1*(*p3);
    }
    
    /* subsequent elimination below */
  
    for (col=1;col<d;col++) {
	for (i=0;i<=LastPlane_;i++)       /* search first row not already used as pivot row*/
	    if (notInPivot(pivot,col,i)) {
		pivot[col]=i; 
		break; 
	    }
	p1=planescopy+i*(d+1)+col;               /* search subsequent pivot row */
	d3=fabs(d1=*p1);
	for (; i<=LastPlane_; i++, p1+=(d+1))  
	    if (notInPivot(pivot,col,i)) {
	        d2=fabs(*(p1));
#if PIVOTING_LASS == 0
		if (d2>=MIN_PIVOT_LASS) {
		    pivot[col]=i; 
		    d1=*p1;
		    break; 
		}
#endif
		if (d2>d3) { 
		    pivot[col]=i;
		    d1=*p1;
		    d3=d2;
		}
	    };
	/* update pivot row */
	p1=planescopy+pivot[col]*(d+1)+col+1;
	d2=1.0/d1;
	for (j=col+1; j<=d; j++, p1++) (*p1) *= d2;
	if (col==(d-1)) break;   /* the rest is not needed in the last case */
        /* update rest of rows */
        p1=planescopy+col+1;
        p2=planescopy+pivot[col]*(d+1)+col+1;
	for (i=0; i<=LastPlane_; i++, p1+=(col+1)) {
	    if (!notInPivot(pivot,col+1,i)) {
	        p1+=d-col;
		continue;
	    }
	    d1=*(p1-1);
	    for (j=col+1; j<=d; j++, p1++, p2++) *p1=(*p1)-d1*(*p2);
	    p2-=d-col;
	}
    };

    /* compute x* by backward substitution; result goes into rhs of planescopy */
  
    for (i=d-2; 0<=i; i--){
        p1=planescopy+pivot[i]*(d+1)+d;
	p2=p1-d+i+1;
	for (j=i+1; j<d; j++, p2++)
	    *(p1)-= (*p2)*(*(planescopy+pivot[j]*(d+1)+d));
    }
 
    /* compute shifted b  */

    for (i=0; i<=LastPlane_; i++) {
        p1=A+i*(d+1);
        p2=p1+d;
	if (notInPivot(pivot,d,i)) 
	    for (j=0; j<d; j++,p1++) {
		*p2 -= (*p1)*(*(planescopy+pivot[j]*(d+1)+d));
	    }
	else *p2=0;
    }
}

static int norm_and_clean_constraints(rational* A, int *LastPlane_, int d, 
                               T_LassInt *Del_index, int Index_needed)
/* Other (simpler) implementation of version lasserre-v15. 
   Finally (up to the sign) identical constraints in A are detected. If they are
   identical the back one is removed, otherwise the system is infeasible. LastPlane_
   is reduced accordingly to the elimination process as well as insertion of the
   corresponding original indices into Del_index if Index_needed is true. */

{   register int i, j, row = 0;
    register rational r0, *p1, *p2;

    /* find nonzero[][] and maximal elements and normalize */
  
    p1=A;                                  /* begin of first constraint */
    while (row<=(*LastPlane_)) {           /* remove zeros and normalize */
	r0=0.0;                            /* norm of vector */
        for (j=0; j<d; j++,p1++) 
	    r0+=(*p1)*(*p1);               /* compute euclidean norm */
        r0=sqrt(r0);
	if (r0<EPS_NORM) {
            if ((*p1)<-100000*EPS1){      /* if negative rhs */
		return 1;                  /* infeasible constraint */
	    }
	    rm_constraint(A, LastPlane_, d,row);
	    if (Index_needed) add_reduced_index(row, Del_index, All_index);
	    p1-=d;
	}
	else {
	    r0=1.0/r0;
	    p1-=d;
	    for (j=0; j<=d; j++,p1++)
		(*p1)*=r0;
	    row++; 
	}
    }

    /* detect identical or reverse constraints */
    
    for (row=0; row<(*LastPlane_); row++) {
	i=row+1;
	while (i<=*LastPlane_) {        /* test all subsequent rows i if equal to row */
            r0=0.0;
 	    p1=A+row*(d+1);
	    p2=A+i*(d+1);
            for (j=0;j<d;j++,p1++,p2++)
	        r0+=(*p1)*(*p2);        /* cosinus of arc among those two vectors */
	    if (r0>0) {
	        /* NEW VERSION of removing constraints */ 
	        if (fabs(r0-1.0)<EPS_NORM) {
		    if ((*p1)>(*p2)){
		    	if (Index_needed) add_reduced_index(row, Del_index, All_index);
			rm_constraint(A, LastPlane_, d,row);
			i=row+1;
                    }
		    else {
			if (Index_needed) add_reduced_index(i, Del_index, All_index);
			if (i<(*LastPlane_)) 
			    rm_constraint(A, LastPlane_, d,i);
			else (*LastPlane_)--;
                    }
		}
                else i++;

                /* OLD VERSION :
	        if ((fabs(r0-1.0)<EPS_NORM) && (fabs((*p1)-(*p2))<EPS1)){
		    if (Index_needed) add_reduced_index(i, Del_index, All_index);
		    if (i<(*LastPlane_)) 
			rm_constraint(A, LastPlane_, d,i);
		    else (*LastPlane_)--;
                }
                else i++;
		*/
	    }
	    else {
	        if (fabs(r0+1.0)<EPS_NORM){
		    if ((*p1)>0){
		        if ((*p2)<(EPS1-(*p1))) return 1; 
		     }
		     else {
		         if ((*p1)<(EPS1-(*p2))) return 1; 
		     }
		}
	        i++;
	    }
	}
    }
    return 0;  /* elimination succesful */
}

static rational lass(rational *A, int LastPlane_, int d)
/* A has exact dimension (LastPlane_+1)*(d+1). The function returns
   the volume; an underscore is appended to LastPlane_ and d */
   
{   rational * redA;            /* A reduced by one dimension and constraint */
    int i, j;
    T_LassInt baserow, basecol, col;
    int dimdiff, row;         /* dimension difference and boolean */
    int i_balance = FALSE;
    rational ma, mi, *volume, *realp1, *realp2;
    int Index_needed;         /* Boolean, if index operations are needed */
    T_LassInt * Del_index; /* contains the indices of the deleted planes */

    /* test if volume is already known and return it if so */

    dimdiff = G_d-d;
    if ((G_Storage > (dimdiff-2)) && (dimdiff >= 2)) {
        tree_out (&tree_volumes, &i_balance, key, &volume, &keyfound, KEY_PLANES_VAR);
        if ((*volume)>=0)  {  /* this volume has already been computed */	
	  printf("\nNeed Scale routine!\n");
	}
	if (!G_OutOfMem) {
	    (*volume)=0;      /* initialize */
	    dimdiff=TRUE;
	}
	else dimdiff=FALSE;
    }
    else dimdiff=FALSE;

    /* if d==1 compute the volume and give it back */

    if (d == 1) {
	ma=-MAXIMUM;
	mi= MAXIMUM;
	for (i=0; i<=LastPlane_; i++,A+=2) { 
	    if (*A>EPSILON_LASS) { if ((*(A+1)/ *A)<mi) mi=(*(A+1)/ *A); }
	    else if (*A<-EPSILON_LASS) { if ((*(A+1)/ *A)>ma) ma=*(A+1)/ *A; } 
            else if ((*(A+1))<-(100000*EPSILON_LASS)) return 0; 
	}
	if ((ma<-.5*MAXIMUM)||(mi>.5*MAXIMUM)) {
	    printf("\nVolume is unbounded!\n");
	    exit(0);
	}
	if ((mi-ma)>EPSILON_LASS) {
	    if (dimdiff) (*volume)=mi-ma;
	    return mi-ma;
	}
	return 0;
    }

    /* if d>1 apply the recursive scheme by fixing constraints. */

    Index_needed = (G_Storage>(G_d-d-1));
    if (Index_needed){
	if (!(Del_index = (T_LassInt *) my_malloc ((LastPlane_ + 2) * sizeof (T_LassInt)))){
	    fprintf (stderr, "\n***** ERROR/WARNING: Out of memory in 'lass'\n");
	    if (G_OutOfMem) exit(0);
	    G_OutOfMem = TRUE;
	    free(G_MemRes);
            Del_index = (T_LassInt *) malloc ((LastPlane_ + 2) * sizeof (T_LassInt));
	};
        Del_index[0]=G_m+2;   /* initialize: mark end */
    }
    ma=0;                                         /* used to sum up the summands */
    if (norm_and_clean_constraints(A, &LastPlane_, d, Del_index, Index_needed)!=0)
        goto label2;

    /* if appropriate shift polytope */

    if (d>=LaShiftLevel) {
	realp1=A+d;
	realp2=realp1+LastPlane_*(d+1);
	j=0;
	while (realp1<=realp2) {
	    if (fabs(*realp1)<EPSILON_LASS) j++;
	    realp1+=d+1;
	}
	if (d-j>=LaShift) shift_P(A, LastPlane_, d);
    }


    redA = (rational *) my_malloc (LastPlane_* d*sizeof(rational));
    if (redA == NULL) {
	fprintf (stderr, "\n***** ERROR/WARNING: Out of memory in 'lass.*redA'\n");
	if (G_OutOfMem) exit(0);
	G_OutOfMem = TRUE;
	free(G_MemRes);
	redA = (rational *) malloc (LastPlane_* d*sizeof(rational));
    }
#ifdef ReverseLass
    for (row=LastPlane_; row>=0; row--) {
#else
    for (row=0; row<=LastPlane_; row++) {
#endif
	if (fabs(*(A+row*(d+1)+d))<EPSILON_LASS) 
            continue;                        /* skip this constraint if b_row == 0 */
	if (Index_needed)
	{  baserow=add_reduced_index(row, NULL, All_index);
           p2c[G_d-d][1] = baserow;
	   add_hypervar (baserow, G_d+1, &key);
	}	
	memcpy(&pivotrow[0], A+row*(d+1), sizeof(rational)*(d+1));
	col=0;                               /* search for pivot column */
	for (i=0; i<d; i++) {        
#if PIVOTING_LASS == 0
	    if (fabs(pivotrow[i])>=MIN_PIVOT_LASS) {col=i; break;};
#endif
	    if (fabs(pivotrow[i])>fabs(pivotrow[col])) col=i;
	};
	if (G_Storage>(G_d-d-1))
	{  basecol=add_reduced_index(col, NULL, Pivot);
           p2c[G_d-d][0] = basecol;
	   add_hypervar (G_m+1, basecol, &key);
	}

        /* copy A onto redA and at the same time perform pivoting */
	 
	mi=1.0/pivotrow[col];
	for (i=0; i<=d; i++) pivotrow[i]*=mi;
	realp1=A;
	realp2=redA;
	for (i=0; i<=LastPlane_; i++) {
	    if (i==row) {
		realp1+=d+1;
		continue;
	    };
	    mi=*(A+(i*(d+1))+col);
	    for (j=0; j<=d; j++) {
		if (j==col) {
		    realp1++;
		    continue;
		};
		*realp2=(*realp1)-pivotrow[j]*mi;
		realp1++;
		realp2++;
	    };
	};
	ma+= *(A+row*(d+1)+d)/(d*fabs(*(A+row*(d+1)+col)))
	     *lass(redA, LastPlane_-1, d-1);
        if (Index_needed)
        {  rm_original_inElAll_index(baserow);
           delete_hypervar (baserow, G_d+1, &key);
        }
	if (G_Storage>(G_d-d-1))
	{  del_original(basecol, Pivot);
	   delete_hypervar (G_m+1, basecol, &key);
	}
        /*$#ifdef verboseFirstLevel
            if (d==G_d) 
	        printf("\nVolume accumulated to iteration %i is %20.12f",row,ma );
        #endif$*/
    };
    my_free (redA, LastPlane_* d * sizeof (rational));
    label2: 
    if (Index_needed) {
	del_original_indices(Del_index, All_index);
        my_free (Del_index, (LastPlane_ + 2) * sizeof (T_LassInt));
    };
    if (dimdiff)(*volume)=ma;
    return ma;
}

/****************************************************************************************/

#ifndef COG
void VOLUME_LASSERRE_FILE ( double planes2[7][4], rational *volume )

  {  int i, j;

  G_d = 3;
  G_m = 7;
  G_Storage = 0;

   /* reserves memory space for the global variable G_Hyperplanes; G_m and G_d must be  */
   /* set correctly                                                                     */
   
  G_Hyperplanes = (real **) my_malloc (G_m * sizeof (real *));
   
  for (i = 0; i < G_m; i++)
    G_Hyperplanes [i] = (real *) my_malloc ((G_d + 1) * sizeof (real));
  /* The last entry is needed for the right hand side. */

  /* Transfer input from Fortran to C array */

  for( i = 0; i <= G_d; i++)
    for( j= 0; j < G_m; j++){
      G_Hyperplanes [j][i] = planes2[j][i];
   /*     printf("\n%10.3e",planes2[j][i]); */
    }
  
/*    fprintf (stderr, "\n****** Hyperplanes are: ******\n"); */
/*    for (j = 0; j < G_m; j++) */
/*      {  fprintf (stderr, "  Hyperplane [%i]:  ", j); */
/*      for (i=0; i < G_d; i++) fprintf (stderr, "%10.3lf", G_Hyperplanes [j] [i]); */
/*      fprintf (stderr,"\n"); */
/*      fprintf (stderr, " :%10.3lf", G_Hyperplanes [j] [G_d]); */
/*      fprintf (stderr,"\n"); */
/*      } */

/*     if (G_Storage > G_d - 3) */
/*        G_Storage = G_d - 3; */
      /* necessary to prevent memory waste because in the tree arrays of length         */
      /* G_Storage + 2 are allocated                                                    */
      
   pivotrow = (rational *) my_malloc ((G_d + 1) * sizeof (rational));
   All_index = (T_LassInt *) my_malloc ((G_m + 1) * sizeof (T_LassInt));
   Pivot = (T_LassInt *) my_malloc ((G_d + 1) * sizeof (T_LassInt));
   p2c = (int **) my_malloc (G_d * sizeof (int *));
   for (i=0; i<G_d; i++){
       p2c[i] = (int *) my_malloc (2 * sizeof (int));
   }
   G_OutOfMem=FALSE;
   G_MemRes=malloc(G_d*G_d*G_m*sizeof(rational)); /* memory reserve; */

   A=compact();
   planescopy=compact();
   tree_volumes = NULL;
   create_key (&key, KEY_PLANES_VAR);
   key.hypervar.hyperplanes [0] = G_m + 1;
   key.hypervar.variables [0] = G_d + 1;
   All_index[0]=G_m+2;  /* initialization (end mark) */
   Pivot[0]=G_m+2;	/* initialization (end mark) */
   *volume = lass (A, G_m-1, G_d);
    
   if (!G_OutOfMem) free (G_MemRes);
   free_key (key, KEY_PLANES_VAR);
}

void volume_lasserre_file ( double planes2[7][4], rational *volume ) {
  VOLUME_LASSERRE_FILE( planes2, volume);
}

void VOLUME_LASSERRE_FILE_( double planes2[7][4], rational *volume ) {
  VOLUME_LASSERRE_FILE( planes2, volume);
}

void volume_lasserre_file_( double planes2[7][4], rational *volume ) {
  VOLUME_LASSERRE_FILE( planes2, volume);
}

#endif

/****************************************************************************************/

/****************************************************************************************/

void tree_out (T_Tree **ppr , int *pi_balance, T_Key key, rational **volume,
   T_Key **keyfound, int key_choice)
   /* looks up the node corresponding to the variable "key" in the specified tree. If   */
   /* it is found the volume is returned via the variable of the same name.             */
   /* At the same time the found key is returned in "foundkey"; this is important for   */
   /* Lasserre, where only a part of the key is compared, but the whole key is needed   */
   /* later.                                                                            */
   /* Otherwise the action taken depends on the memory reserves: If G_OutOfMem is FALSE */
   /* a new node is created and a pointer to its volume part is returned via the varia- */
   /* ble "volume" so that the computed volume can be inserted by the calling routine.  */
   /* If there is no more memory the volume -1 is returned.                             */
   /* As in the previous routine "key_choice" determines the active part of the keys.   */
   
{  T_Tree *p1, *p2;
   int	  cmp;

   /* Are we grounded? If so, add the node here and set the rebalance flag, then exit.  */
   if (!*ppr)
   {
      if (G_OutOfMem)
      {  /* don't allocate if out of memory */
	 *volume = &G_Minus1;
      }
      else
      {  (*ppr) = (T_Tree *) my_malloc (sizeof (T_Tree));
         (*ppr) -> tree_l = NULL;
	 (*ppr) -> tree_r = NULL;
	 (*ppr) -> tree_b = 0;
	 /* copy the key into the new node */
	 create_key (&((*ppr) -> key), key_choice);
	 switch (key_choice)
	 {
	 case KEY_PLANES:
	    memcpy ((*ppr) -> key.hyperplanes, key.hyperplanes,
	            (G_d + 1) * sizeof (int));
	    break;
	 case KEY_VERTICES:
	    (*ppr) -> key.vertices.set = duplicate_set (key.vertices.set);
	    (*ppr) -> key.vertices.d = key.vertices.d;
	    break;
	 case KEY_PLANES_VAR:
	    memcpy ((*ppr) -> key.hypervar.hyperplanes, key.hypervar.hyperplanes,
	            (G_Storage + 2) * sizeof (T_LassInt));
	    memcpy ((*ppr) -> key.hypervar.variables, key.hypervar.variables,
	            (G_Storage + 2) * sizeof (T_LassInt));
	    break;
	 }
	 (*ppr) -> vol = -1;       /* to recognise that element is newly created */
	 *volume = &((*ppr) -> vol);
	 *pi_balance = TRUE;
      }
      return;
   }

   cmp = compare_key ((*ppr) -> key, key, key_choice);

   /* if LESS, prepare to move to the left. */
   if (cmp < 0)
   {
      tree_out (&((*ppr) -> tree_l), pi_balance, key, volume, keyfound, key_choice);
      if (*pi_balance)
      {  /* left branch has grown longer */
         switch ((*ppr) -> tree_b)
         {
         case 1:  /* right branch WAS longer; balance is ok now */
                  /* LESS: case 1.. balance restored implicitly */
            (*ppr) -> tree_b = 0;
            *pi_balance = FALSE;
            break;
         case 0:  /* balance WAS okay; now left branch longer */
                  /* LESS: case 0.. balance bad but still ok */
            (*ppr) -> tree_b = -1;
            break;
         case -1: /* left branch was already too long. rebalance */
            p1 = (*ppr) -> tree_l;
            if (p1 -> tree_b == -1)
            {  /* LESS: single LL */
               (*ppr) -> tree_l = p1->tree_r;
               p1 -> tree_r = (*ppr);
               (*ppr) -> tree_b = 0;
               (*ppr) = p1;
            }
            else
            {  /* LESS: real LR */
               p2 = p1 -> tree_r;
               p1 -> tree_r = p2 -> tree_l;
               p2 -> tree_l = p1;
               (*ppr) -> tree_l = p2 -> tree_r;
               p2 -> tree_r = (*ppr);
               if (p2 -> tree_b == -1)
                  (*ppr) -> tree_b = 1;
               else (*ppr) -> tree_b = 0;
               if (p2->tree_b == 1)
                  p1 -> tree_b = -1;
               else p1 -> tree_b = 0;
               (*ppr) = p2;
            }
            (*ppr) -> tree_b = 0;
            *pi_balance = FALSE;
         } /* switch */
      } /* if */
   } /* cmp < 0 */

   /* if MORE, prepare to move to the right. */
   else if (cmp > 0)
   {
      tree_out (&((*ppr) -> tree_r), pi_balance, key, volume, keyfound, key_choice);
      if (*pi_balance)
      {  /* right branch has grown longer */
         switch ((*ppr) -> tree_b)
         {
         case -1: /* MORE: balance was off, fixed implicitly */
            (*ppr) -> tree_b = 0;
            *pi_balance = FALSE;
            break;
         case 0:  /* MORE: balance was okay, now off but ok */
            (*ppr)->tree_b = 1;
            break;
         case 1:  /* MORE: balance was off, need to rebalance */
            p1 = (*ppr) -> tree_r;
            if (p1 -> tree_b == 1)
            {  /* MORE: single RR */
               (*ppr) -> tree_r = p1 -> tree_l;
               p1 -> tree_l = (*ppr);
               (*ppr) -> tree_b = 0;
               (*ppr) = p1;
            }
            else
            {  /* MORE: real RL */
               p2 = p1 -> tree_l;
               p1 -> tree_l = p2 -> tree_r;
               p2 -> tree_r = p1;
               (*ppr) -> tree_r = p2 -> tree_l;
               p2 -> tree_l = (*ppr);
               if (p2 -> tree_b == 1)
                  (*ppr) -> tree_b = -1;
               else (*ppr) -> tree_b = 0;
               if (p2 -> tree_b == -1)
                  p1 -> tree_b = 1;
               else p1 -> tree_b = 0;
               (*ppr) = p2;
            }
            (*ppr) -> tree_b = 0;
            *pi_balance = FALSE;
         } /* switch */
      } /* if */
   } /* cmp > 0 */

   /* not less, not more: this is the same key! give volume back! */
   else
   {
      *pi_balance = FALSE;
      *volume = &((*ppr) -> vol);
      *keyfound = &((*ppr) -> key);
   }
}
/****************************************************************************************/

T_VertexSet duplicate_set (T_VertexSet s)
   /* creates a new set with the same elements as s */
   
{  T_VertexSet newset;

   newset.maxel = s.lastel;
   newset.lastel = s.lastel;
   newset.loe = (T_Vertex **) my_malloc ((s.lastel + 1) * sizeof (T_Vertex *));
   memcpy (newset.loe, s.loe, (s.lastel + 1) * sizeof (T_Vertex *));
   return newset;
}
/****************************************************************************************/
/*       routines for storing intermediate volumes in balanced trees (avl-trees)        */
/****************************************************************************************/

static int compare_key (T_Key key1, T_Key key2, int key_choice)
   /* compares key1 with key2; if key1 is smaller, -1 is returned, if key1 is bigger +1 */
   /* and if both are equal 0. key_choice determines which component of the key is      */
   /* relevant for comparing.                                                           */
   
{  int       i;
   T_LassInt *p1, *p2;
   
   switch (key_choice)
   {
   case KEY_PLANES:
      for (i=0; ; i++)
         if      (key1.hyperplanes [i] < key2.hyperplanes [i]) return -1;
	 else if (key1.hyperplanes [i] > key2.hyperplanes [i]) return 1;
         else if (key1.hyperplanes [i] > G_m)                  return 0;
      break;
   case KEY_VERTICES:
      if      (key1.vertices.d < key2.vertices.d) return -1;
      else if (key1.vertices.d > key2.vertices.d) return 1;
      else /* both volumes are for the same dimension */
      if      (key1.vertices.set.lastel < key2.vertices.set.lastel) return -1;
      else if (key1.vertices.set.lastel > key2.vertices.set.lastel) return 1;
      else
      {  /* both sets have the same cardinality */
         for (i=0; i <= key1.vertices.set.lastel; i++)
            if      (key1.vertices.set.loe [i] -> no < key2.vertices.set.loe [i] -> no)
               return -1;
            else if (key1.vertices.set.loe [i] -> no > key2.vertices.set.loe [i] -> no)
               return 1;
      }
      return 0;
      break;
   case KEY_PLANES_VAR:
/*
      for (i=0; ; i++)
         if      (key1.hypervar.hyperplanes [i] < key2.hypervar.hyperplanes [i])
                                                                        return -1;
	 else if (key1.hypervar.hyperplanes [i] > key2.hypervar.hyperplanes [i])
	                                                                return 1;
         else if (key1.hypervar.hyperplanes [i] > G_m)                  return 0;
*/
      for (p1 = key1.hypervar.hyperplanes, p2 = key2.hypervar.hyperplanes;;
           p1++, p2++)
         if      ((*p1) < (*p2)) return -1;
         else if ((*p1) > (*p2)) return  1;
         else if ((*p1) > G_m)   return  0;
      break;
   }
}
