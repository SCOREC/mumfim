//#include <stdlib.h> 
//gcc2.96 
#include <iostream> 
#include <memory> 
#include <math.h>
#include <assert.h>
#include "CSR_Matrix.h"
#include "listman.h"
#include <cstring>
#include <stdio.h>
#include <cstdlib>


#ifndef sgi
#ifndef MSC
 using std::memcpy;
#endif
#endif

#define SWAP(a,b)  temp=(a);(a)=(b);(b)=temp;
#define SWAPI(a,b) tempi=(a);(a)=(b);(b)=tempi;

const int M_sort2  = 7;
const int NSTACK   = 50;


CSR_Matrix::CSR_Matrix() : 
  a_(NULL), ai_(NULL), ptr_(NULL), jptr_(NULL), 
  NbLines(0), AllocationDone(false), storage(UNKNOWN_STORAGE),
  MatrixWasReordered(false), 
  permr_(NULL), permp_(NULL), rpermr_(NULL)  { }


void CSR_Matrix::Allocate(const int n) 
{
  NbLines = n;  
  storage = HARWELL_BOEING_STORAGE ;
  a_    = List_Create (NbLines, NbLines, sizeof(double));
  ai_   = List_Create (NbLines, NbLines, sizeof(int));
  ptr_  = List_Create (NbLines, NbLines, sizeof(int));
  jptr_ = List_Create (NbLines+1, NbLines, sizeof(int));
  int i, j=0;
  for(i=0; i<NbLines; i++) List_Add (jptr_, &j);
  ZeroMatrix ( );

  // Add a zero to each entry on the diagonal -- this gives better structure
  // and avoids a bug in the matrix vector multiply
  //printf("Before addmatrix in Allocate\n");

  for (i = 1; i <= NbLines; i++) {
    AddMatrix(i,i,0.0);
  }
  AllocationDone = true;

  delete []rpermr_; rpermr_ = NULL;
  delete []permr_;  permr_ = NULL;
  delete []permp_;  permp_ = NULL;
  MatrixWasReordered = false;
}

void CSR_Matrix::ReAllocate(const int n, const int nbz) 
{
  List_Delete(a_);
  List_Delete(ai_);
  List_Delete(ptr_);
  List_Delete(jptr_);
  delete []rpermr_; rpermr_ = NULL;
  delete []permr_;  permr_ = NULL;
  delete []permp_;  permp_ = NULL;

  NbLines = n;  
  storage = HARWELL_BOEING_STORAGE ;
  a_    = List_Create (nbz, NbLines, sizeof(double));
  ai_   = List_Create (nbz, NbLines, sizeof(int));
  ptr_  = List_Create (nbz, NbLines, sizeof(int));
  jptr_ = List_Create (NbLines+1, NbLines, sizeof(int));
  int i, j=0;
  for(i=0; i<NbLines; i++) List_Add (jptr_, &j);
  ZeroMatrix ( );

  for (i = 1; i <= NbLines; i++) 
    {
      AddMatrix(i,i,0.0);
    }
  AllocationDone = true;
  MatrixWasReordered = false;
}


CSR_Matrix::CSR_Matrix(const int n) : 
a_(NULL), ai_(NULL), ptr_(NULL), jptr_(NULL), 
NbLines(0), AllocationDone(false), 
storage(UNKNOWN_STORAGE),
MatrixWasReordered(false), 
permr_(NULL), permp_(NULL), rpermr_(NULL)
{
  Allocate(n);
  return ;
}

void CSR_Matrix::SetSize(const int& n) {

  if (AllocationDone) { 
    fprintf(stderr, "Error: Matrix already allocated\n");
    assert(0); 
  }
  Allocate(n);
  return ;  
}

void CSR_Matrix::Clear(){

    NbLines = 0;
    storage = UNKNOWN_STORAGE;
    AllocationDone = false;
    MatrixWasReordered = false;

//printf("Hello1\n");
    List_Delete(a_);  a_ = NULL;
//printf("Hello2\n");
    List_Delete(ai_); ai_ = NULL;
//printf("Hello3\n");
    List_Delete(ptr_); ptr_ = NULL;
//  free(ptr_->array);
////printf("Hello3.1\n");
//  assert(ptr_ != 0);
////printf("Hello3.2\n");
//  free(ptr_);
//printf("Hello4\n");
//if (jptr_) printf("jptr_ n'est pas nul\n");
//else       printf("jptr_   est     nul\n");
    List_Delete(jptr_); jptr_ = NULL;
//printf("Hello5\n");
    delete []rpermr_; rpermr_ = NULL;
//printf("Hello6\n");
    delete []permr_; permr_ = NULL;
//printf("Hello7\n");
    delete []permp_; permp_ = NULL;
//printf("Hello8\n");
}


CSR_Matrix::~CSR_Matrix(){
//printf("Hellui1\n");
    List_Delete(a_);
//printf("Hellui2\n");
    List_Delete(ai_); 
//printf("Hellui3\n");
    List_Delete(ptr_);
//  free(ptr_->array);
////printf("Hellui3.1\n");
//  assert(ptr_ != 0);
////printf("Hellui3.2\n");
//  free(ptr_);
//printf("Hellui4\n");
    List_Delete(jptr_); 
//printf("Hellui5\n");
if ( MatrixWasReordered == true ) {
    delete []rpermr_;
//printf("Hellui6\n");
    delete []permr_;
//printf("Hellui7\n");
    delete []permp_;
//printf("Hellui8\n");
}
}



 int CSR_Matrix :: GetNbUnknown() {

 return NbLines ;

 }

void  CSR_Matrix :: ZeroMatrix ( ){


  List_Reset (a_);
  List_Reset (ai_);
  List_Reset (ptr_);
  List_Reset (jptr_);
  int i,j=0;
  for (i=0; i< NbLines; i++) List_Add ( jptr_, &j);
}

void  CSR_Matrix::AddMatrix ( int ic, int il, double val) {

   assert( MatrixWasReordered == false );

  // attention a la transposition 

  int     *ai, *pp, n, iptr, iptr2, jptr, *ptr, zero = 0;
  double   *a;

  if ( storage == COMPRESSED_ROW_STORAGE ) PutInHarwellBoeing();

  //printf( "A(%d,%d) = %12.8e Number of Non Zero %d \n",ic,il,val,List_Nbr(a_));

   il--;
   pp  = (int*) jptr_->array;
   ptr = (int*) ptr_->array;
   ai  = (int*) ai_->array;
   a   = (double*) a_->array;
   
   iptr = pp[il];
   iptr2 = iptr-1;
   
   while(iptr>0){
     iptr2 = iptr-1;
     jptr = ai[iptr2];
     if(jptr == ic){
       a[iptr2] += val;
       return;
     }
     iptr = ptr[iptr2];
   }
   
   List_Add (a_, &val);
   List_Add (ai_, &ic);
   List_Add (ptr_, &zero);
   
   // The pointers amy have been modified
   //   if there has been a reallocation in List_Add  
   
   ptr = (int*) ptr_->array;
   ai  = (int*) ai_->array;
   a   = (double*) a_->array;
   
   n = List_Nbr(a_);

   if(!pp[il]) pp[il] = n;
   else ptr[iptr2] = n;
   
}

double CSR_Matrix::GetMatrix ( int ic, int il) {

  if ( storage == HARWELL_BOEING_STORAGE && MatrixWasReordered == false) { 

   //printf("HellllllllloBow\n");

    double val;
    int     *ai, *pp, iptr, iptr2, jptr, *ptr; 
    double   *a;

     int iptr_last,iptr_last2;

      il--;
      //printf("Hela1\n");
      pp  = (int*) jptr_->array;
      //printf("Hela2\n");
      ptr = (int*) ptr_->array;
      //printf("Hela3\n");
      ai  = (int*) ai_->array;
      //printf("Hela4\n");
      a   = (double*) a_->array;
      //printf("Hela5\n");
    
      iptr = pp[il];
      //printf("Hela6\n");
      iptr2 = iptr-1;
    
     iptr_last = -1;  
     iptr_last2 = -1;  

      while(iptr>0){
        iptr2 = iptr-1;
      //printf("Hela7\n");
        jptr = ai[iptr2];
      //printf("Hela8\n");


        if(jptr == ic){

      //printf("Hela9\n");
      assert(a != NULL);
      //printf("iptr2 = %d\n", iptr2);
	val = a[iptr2];
      //printf("Hela10\n");
      //printf("val = %12.5e\n", val);
        return val;
      }
      //printf("Hela11\n");
      iptr = ptr[iptr2];
      //printf("Hela12\n");
    }
    return 0.0;

  } else if ( storage == COMPRESSED_ROW_STORAGE && MatrixWasReordered == false) { 


//  double val;
  int     *row_ptr, *col_ptr;
  double   *val_ptr;
//  int iptr_last,iptr_last2;

    val_ptr    = (double*) a_->array;
    row_ptr    = (int*)    jptr_->array;
    col_ptr    = (int*)    ptr_->array;

        for (int t=row_ptr[ic-1]-1; t<row_ptr[ic]-1; t++)
           if (col_ptr[t] == il ) return val_ptr[t];

        if ( 1 <=ic && ic <= NbLines && 1 <=il && il <= NbLines ) return 0.0;

        else
        {
            printf( "Array accessing exception -- out of bounds.");
            assert( 1==0 );
        }

  } else if ( storage == COMPRESSED_ROW_STORAGE && MatrixWasReordered == true ) { 

//  double val;
  int     *row_ptr, *col_ptr;
  double   *val_ptr;
//  int iptr_last,iptr_last2;


    val_ptr    = (double*) a_->array;
    row_ptr    = (int*)    jptr_->array;
    col_ptr    = (int*)    ptr_->array;

  for (int t=row_ptr[ic-1]-1; t<row_ptr[ic]-1; t++)
    if (col_ptr[t] == il ) return val_ptr[t];
  
  if ( 1 <=ic && ic <= NbLines && 1 <=il && il <= NbLines ) return 0.0;
  
  else
    {
      printf( "Array accessing exception -- out of bounds.");
      assert( 1==0 );
    }
  
} else { 
  
  printf( "CSR_Matrix : Matrix storage is corrupted. \n"); 
  
  assert (1==0); };
  
}


int CSR_Matrix::ReturnNumberOfNonZero ( void ) { 


  //printf(" Number of non zero %d\n",List_Nbr(a_) );

   return  List_Nbr(a_);

}

double * CSR_Matrix::ReturnValPointer ( void ) 
{   
  return  (double*) a_->array;
}

int * CSR_Matrix::ReturnRowPointer  ( void ) 
{   
  return  (int*) jptr_->array;  
}

int * CSR_Matrix::ReturnColumIndexPointer( void )
{ 
  return (int*)  ptr_->array;  
}


void  CSR_Matrix :: EndOfAssembly() {

 PutInCompressedRow();

 }


void CSR_Matrix::PutInCompressedRow() {
//void csr_format in sort.c

 if (storage == HARWELL_BOEING_STORAGE) {

   storage = COMPRESSED_ROW_STORAGE;

   int    i,*ptr,*jptr,*ai,iptr,iptr2;
   double *a;

   ptr  = (int*)ptr_->array;
   jptr = (int*)jptr_->array;
   ai   = (int*)ai_->array;
   a    = (double*)a_->array;
   

   for(i=0;i<NbLines;i++){
     iptr = jptr[i];
     while(iptr){
       iptr2 = iptr - 1;
       iptr = ptr[iptr2];
       ptr[iptr2] = i+1;
     }
   }
   
   sort2(List_Nbr(a_),a,ai,ptr);
   
   deblign(List_Nbr(a_),ptr,jptr,ai);
   
   jptr[NbLines]=List_Nbr(a_)+1;

 }

 return;
}


void CSR_Matrix::PutInHarwellBoeing() {
                       
 if (storage == COMPRESSED_ROW_STORAGE) {
  storage = HARWELL_BOEING_STORAGE;
  char* temp  = ptr_->array;
  ptr_->array = ai_->array;
  ai_->array = temp;
}
  return;
} 


void  CSR_Matrix::sort2(unsigned long n, double arr[], int ai[] , int aj [] ) {

  unsigned long i,ir=n,j,k,l=1;
  int *istack,jstack=0,tempi;
  double a,temp;
  int    b,c;
    
  istack=ivector(1,NSTACK);
  for (;;) {
    if (ir-l < M_sort2) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j -1];
	b=ai[j -1];
	c=aj[j -1];
	for (i=j-1;i>=1;i--) {
	  if (cmpij(ai[i -1],aj[i -1],b,c) <= 0) break;
	  arr[i+1 -1]=arr[i -1];
	  ai[i+1 -1]=ai[i -1];
	  aj[i+1 -1]=aj[i -1];
	}
	arr[i+1 -1]=a;
	ai[i+1 -1]=b;
	aj[i+1 -1]=c;
      }
      if (!jstack) {
	free_ivector(istack,1,NSTACK);
	return;
      }
      ir=istack[jstack];
      l=istack[jstack-1];
      jstack -= 2;
    } 
    else {
      k=(l+ir) >> 1;
      SWAP(arr[k -1],arr[l+1 -1])
      SWAPI(ai[k -1],ai[l+1 -1])
      SWAPI(aj[k -1],aj[l+1 -1])
      if (cmpij(ai[l+1 -1],aj[l+1 -1],ai[ir -1],aj[ir -1])>0){
	SWAP(arr[l+1 -1],arr[ir -1])
	SWAPI(ai[l+1 -1],ai[ir -1])
	SWAPI(aj[l+1 -1],aj[ir -1])
      }
      if (cmpij(ai[l -1],aj[l -1],ai[ir -1],aj[ir -1])>0){
	SWAP(arr[l -1],arr[ir -1])
	SWAPI(ai[l -1],ai[ir -1])
	SWAPI(aj[l -1],aj[ir -1])
      }
      if (cmpij(ai[l+1 -1],aj[l+1 -1],ai[l -1],aj[l -1])>0){
	SWAP(arr[l+1 -1],arr[l -1])
	SWAPI(ai[l+1 -1],ai[l -1])
	SWAPI(aj[l+1 -1],aj[l -1])
      }
      i=l+1;
      j=ir;
      a=arr[l -1];
      b=ai[l -1];
      c=aj[l -1];
      for (;;) {
	do i++; while (cmpij(ai[i -1],aj[i -1],b,c) < 0);
	do j--; while (cmpij(ai[j -1],aj[j -1],b,c) > 0);
	if (j < i) break;
	SWAP(arr[i -1],arr[j -1])
	SWAPI(ai[i -1],ai[j -1])
	SWAPI(aj[i -1],aj[j -1])
	}
      arr[l -1]=arr[j -1];
      arr[j -1]=a;
      ai[l -1]=ai[j -1];
      ai[j -1]=b;
      aj[l -1]=aj[j -1];
      aj[j -1]=c;
      jstack += 2;
      if (jstack > NSTACK) {
	fprintf(stderr,"NSTACK too small in sort2.\n");
	exit(1);
      }
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } 
      else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}


int  CSR_Matrix::cmpij(int ai,int aj,int bi,int bj){
  if(ai<bi)return -1;
  if(ai>bi)return 1;
  if(aj<bj)return -1;
  if(aj>bj)return 1;
  return 0;
}

int * CSR_Matrix::ivector(long nl, long nh) {
  // allocate an int vector with subscript range v[nl..nh] 
  int *v;
  
  v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
  if (!v) fprintf(stderr, "allocation failure in ivector()\n");
  return v-nl+1;
}

 void CSR_Matrix :: free_ivector(int *v, long nl, long nh){
  // free an int vector allocated with ivector() 
  free((char*)(v+nl-1));
}


void CSR_Matrix::ExecuteReordering() {

  if ( MatrixWasReordered == false ){

  PutInCompressedRow();

  double *a;
  int    *jptr,*ai;
  int    *mask, *levels;

  int i,j,k;
  a    = (double*) a_->array;
  jptr = (int*)    jptr_->array;
  ai   = (int*)    ptr_->array; 

//  permr_  = (int*) Malloc (NbLines * sizeof(int));
//  rpermr_ = (int*) Malloc (NbLines * sizeof(int));
//  permp_  = (int*) Malloc (2 * NbLines * sizeof(int));
  permr_  = new int[NbLines];
  rpermr_ = new int[NbLines];
  permp_  = new int[2 * NbLines];

  for (i=0 ; i<NbLines ; i++) {
    permr_[i] = rpermr_[i] = permp_[i+NbLines] = i+1;
  }

  //  fprintf(stdout, "# --> RCMK renumbering\n");

  double *a_rcmk;
  int    *jptr_rcmk, *ai_rcmk;

  a_rcmk    = (double*) Malloc (List_Nbr(a_) * sizeof(double));
  jptr_rcmk = (int*)    Malloc ((NbLines + 1) * sizeof(int));
  ai_rcmk   = (int*)    Malloc (List_Nbr(a_) * sizeof(int));

  mask   = (int*) Malloc (List_Nbr(a_) * sizeof(int));
  levels = (int*) Malloc ((NbLines + 1) * sizeof(int));
  
  i = j = k = 1;
  
  cmkreord  (NbLines, a, ai, jptr, a_rcmk, ai_rcmk, jptr_rcmk, &i,
	     permr_, mask, &j, &k, rpermr_, levels);

  double *vv;
  vv = (double*) Malloc (List_Nbr(a_) * sizeof(double));

  sort_col  ( NbLines, a_rcmk, ai_rcmk, jptr_rcmk, mask, vv);
  
  Free(vv);
  Free(mask);
  Free(levels);

  for ( i=0 ; i < List_Nbr(a_) ; i++ ){  
   a[i]    = a_rcmk[i];
   ai[i]   = ai_rcmk[i];   
 }

  for ( i=0 ; i < NbLines+1 ; i++ ) jptr[i] = jptr_rcmk[i];

    Free (a_rcmk);
    Free (jptr_rcmk);
    Free (ai_rcmk);

    //  printf("The number of entries in the matrix is %d\n", List_Nbr(a_) );


#if DEBUGNIC>=1
   printf("CSR::The row pointer is\n");
   for (i = 0; i < NbLines + 1; i++) {
     printf("%d ", jptr[i]);
   }
   printf("\n");
   printf("The col pointer and value are\n");
   for (i = 0; i < List_Nbr(a_); i++) {
     printf("%d %12.5e\n", ai[i], a[i]);
   }
#endif

 MatrixWasReordered = true ;
}

}

void CSR_Matrix::ReorderArray(double * array) {

 if ( MatrixWasReordered  == true ){
 
 int i;
 double * temp;
 temp = new double [ NbLines ];

  for ( i=0;i<NbLines;i++) temp[i] = array[rpermr_[i] - 1];
 
  for ( i=0;i<NbLines;i++) array[i] = temp[i];
  delete [] temp;
 }

}

void CSR_Matrix::InverseReorderArray(double * array) {


 if ( MatrixWasReordered == true ){

 double * temp;
 temp = new double [ NbLines ];
 
 int i,j,k;

   for ( i=0;i<NbLines;i++) {
    j = permr_[i] - 1;
    k = permp_[j+NbLines] - 1;        
    temp[i] = array[k];
   }  

    for ( i=0;i<NbLines;i++) array[i] = temp[i];
 
 delete [] temp;

 }
}

void CSR_Matrix::cmkreord(int n, double *a, int *ja, int *ia, double *a0,
			   int *ja0, int *ia0, int *init, int * iperm, int * mask,
			   int * maskval, int * nlev, int * riord, int * levels) {

    int c__1 = 1;
    // System generated locals 
    int i__1;

    // Local variables 
//    extern // Subroutine  int exchange_();
    int i;
//    extern // Subroutine  int dperm_(), perphn_(), rversp_();

// -----------------------------------------------------------------------
 
//      Cuthill-McKee Reordering : call to perphn 
//                                 renumber the nodes 
// -------------------------parameters------------------------------------
 
// on entry: 
// ---------- 
// n      = number of nodes in the graph 
// ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data 

//          structure) 
// init   = initial node 
// iperm  = integer array indicating in which order to  traverse the graph
 
//          in order to generate all connected components. 
//          The nodes will be traversed in order iperm(1),....,iperm(n) 
//          Convention: 
//          if iperm(1) .eq. 0 on entry then BFS will traverse the 
//          nodes in the  order 1,2,...,n. 

// riord  = (also an ouput argument). on entry riord contains the labels 

//          of the nfirst nodes that constitute the first level. 

// mask   = array used to indicate whether or not a node should be 
//          condidered in the graph. see maskval. 
//          mask is also used as a marker of  visited nodes. 

// maskval= consider node i only when:  mask(i) .eq. maskval 
//          maskval must be .gt. 0. 
//          thus, to consider all nodes, take mask(1:n) = 1. 
//          maskval=1 (for example) 

// on return 
// --------- 
// mask   = on return mask is restored to its initial state. 
// riord  = `reverse permutation array'. Contains the labels of the nodes 

//          constituting all the levels found, from the first level to 
//          the last. 
// levels = pointer array for the level structure. If lev is a level 
//          number, and k1=levels(lev),k2=levels(lev+1)-1, then 
//          all the nodes of level number lev are: 
//          riord(k1),riord(k1+1),...,riord(k2) 
// nlev   = number of levels found 
// -----------------------------------------------------------------------
 
//    local variables 
//     print*,'   ... Searching levels' 
    // Parameter adjustments 
    --levels;
    --riord;
    --mask;
    --iperm;
    --ia0;
    --ja0;
    --a0;
    --ia;
    --ja;
    --a;

    // Function Body 
    *maskval = 1;
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
	mask[i] = *maskval;
    }
    iperm[1] = 0;
    perphn(n, &ja[1], &ia[1], init, &iperm[1], &mask[1], maskval, nlev, &
	    riord[1], &levels[1]);
    rversp(n, &riord[1]);
//     print*,'   ... Renumbering' 
    exchange(n, &riord[1], &iperm[1]);
    dperm(n, &a[1], &ja[1], &ia[1], &a0[1], &ja0[1], &ia0[1], &iperm[1], &
	    iperm[1], &c__1);

} 

int CSR_Matrix::perphn(int n, int *ja, int *ia, int *init, int *iperm,
			  int * mask, int * maskval, int * nlev, 
			  int * riord, int * levels){
    // System generated locals 
    int i__1;

    // Local variables 
    int j, nlevp, mindeg, nfirst, deg;
//    extern // Subroutine  int bfs_();
    int nod;
//    extern integer maskdeg_();

// -----------------------------------------------------------------------
 
//     finds a pseudo-peripheral node and does a BFS search from it. 
// -----------------------------------------------------------------------
 
// see routine  dblstr for description of parameters 
// input: 
// ------- 
// ja, ia  = list pointer array for the adjacency graph 
// mask    = array used for masking nodes -- see maskval 
// maskval = value to be checked against for determing whether or 
//           not a node is masked. If mask(k) .ne. maskval then 
//           node k is not considered. 
// init    = init node in the pseudo-peripheral node algorithm. 

// output: 
// ------- 
// init    = actual pseudo-peripherial node found. 
// nlev    = number of levels in the final BFS traversal. 
// riord   = 
// levels  = 
// -----------------------------------------------------------------------
 
    // Parameter adjustments 
    --levels;
    --riord;
    --mask;
    --iperm;
    --ia;
    --ja;

    // Function Body 
    nlevp = 0;
L1:
    riord[1] = *init;
    nfirst = 1;
    bfs(n, &ja[1], &ia[1], &nfirst, &iperm[1], &mask[1], maskval, &riord[1], 
	    &levels[1], &(*nlev));

    if (*nlev > nlevp) {
	mindeg = levels[*nlev + 1] - 1;
	i__1 = levels[*nlev + 1] - 1;
	for (j = levels[*nlev]; j <= i__1; ++j) {
	    nod = riord[j];
	    deg = maskdeg(&ja[1], &ia[1], &nod, &mask[1], maskval);
	    if (deg < mindeg) {
		*init = nod;
		mindeg = deg;
	    }
	}
	nlevp = *nlev;
	goto L1;
    }
    return 0;
} 

int  CSR_Matrix::rversp(int n, int *riord) {

    // System generated locals 
    int i__1;

    // Local variables 
    int j, k;

// -----------------------------------------------------------------------
//
//     this routine does an in-place reversing of the permutation array 
//     riord -- 
// -----------------------------------------------------------------------

    // Parameter adjustments //
    --riord;

    // Function Body
    i__1 = n / 2;
    for (j = 1; j <= i__1; ++j) {
	k = riord[j];
	riord[j] = riord[n - j + 1];
	riord[n - j + 1] = k;
// L26: 
    }
    return 0;
}


int CSR_Matrix::maskdeg(int *ja, int *ia, int *nod,
			  int *mask, int *maskval) {

    // System generated locals 
    int ret_val, i__1;

    // Local variables 
    int k, deg;

// ----------------------------------------------------------------------- 
    // Parameter adjustments 
    --mask;
    --ia;
    --ja;

    // Function Body 
    deg = 0;
    i__1 = ia[*nod + 1] - 1;
    for (k = ia[*nod]; k <= i__1; ++k) {
	if (mask[ja[k]] == *maskval) {
	    ++deg;
	}
    }
    ret_val = deg;
    return ret_val;
}

void CSR_Matrix::exchange(int n, int *iriord, int *iperm) {
    // System generated locals 
    int i__1;

    // Local variables 
    int i;

//-----------------------------------------------------------------------
//  Reverse a permutation vector 
//  On entry : 
//  ---------- 
//            n      : dimension 
//            iriord : initial reordering vector 
//  On return : 
//  ----------- 
//           iperm  : permutation vector to be used with SPARSKIT (dperm .
//-----------------------------------------------------------------------
//     local variable 
    // Parameter adjustments 
    --iperm;
    --iriord;

    // Function Body 
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
	iperm[iriord[i]] = i;
    }
} 

void CSR_Matrix::sort_col(int n,double * a, int * ja, int *ia,
			    int * iw, double * rw){
    // System generated locals 
    int i__1, i__2;

    // Local variables 
    int ideb, ifin;
//    extern  int sort_irv__();
    int i, j, k;

    // Parameter adjustments 
    --rw;
    --iw;
    --ia;
    --ja;
    --a;

    // Function Body 
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
	ideb = ia[i];
	ifin = ia[i + 1] - 1;
	k = 0;
	i__2 = ifin;
	for (j = ideb; j <= i__2; ++j) {
	    ++k;
	    iw[k] = ja[j];
	    rw[k] = a[j];
	}
	i__2 = ifin - ideb + 1;
	sort_irv(&iw[1], &rw[1], &i__2);
	k = 0;
	i__2 = ifin;
	for (j = ideb; j <= i__2; ++j) {
	    ++k;
	    ja[j] = iw[k];
	    a[j] = rw[k];
	}
    }
    return;
}


int CSR_Matrix::dperm(int nrow, double *a, int *ja, int * ia, double * ao,
			int * jao, int * iao, int * perm, int *qperm, int * job) {

    int locjob;

// -----------------------------------------------------------------------
 
// This routine permutes the rows and columns of a matrix stored in CSR 
// format. i.e., it computes P A Q, where P, Q are permutation matrices. 

// P maps row i into row perm(i) and Q maps column j into column qperm(j):
 
//      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix 
// In the particular case where Q is the transpose of P (symmetric 
// permutation of A) then qperm is not needed. 
// note that qperm should be of length ncol (number of columns) but this 

// is not checked. 
// -----------------------------------------------------------------------
 
// Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991. 
// -----------------------------------------------------------------------
 
// on entry: 
// ---------- 
// n 	= dimension of the matrix 
// a, ja, 
//    ia = input matrix in a, ja, ia format 
// perm 	= integer array of length n containing the permutation arrays 
// 	  for the rows: perm(i) is the destination of row i in the 
//         permuted matrix -- also the destination of column i in case 
//         permutation is symmetric (job .le. 2) 

// qperm	= same thing for the columns. This should be provided only 
//         if job=3 or job=4, i.e., only in the case of a nonsymmetric 
// 	  permutation of rows and columns. Otherwise qperm is a dummy 

// job	= integer indicating the work to be done: 
// * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P) 
// 		job = 1	permute a, ja, ia into ao, jao, iao 
// 		job = 2 permute matrix ignoring real values. 
// * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q 
// 		job = 3	permute a, ja, ia into ao, jao, iao 
// 		job = 4 permute matrix ignoring real values. 

// on return: 
// ----------- 
// ao, jao, iao = input matrix in a, ja, ia format 

// in case job .eq. 2 or job .eq. 4, a and ao are never referred to 
// and can be dummy arguments. 
// Notes: 
// ------- 
//  1) algorithm is in place 
//  2) column indices may not be sorted on return even  though they may be
 
//     on entry. 
// ----------------------------------------------------------------------c
 
// local variables 

//     locjob indicates whether or not real values must be copied. 

    // Parameter adjustments 
    --qperm;
    --perm;
    --iao;
    --jao;
    --ao;
    --ia;
    --ja;
    --a;

    // Function Body 
    locjob = *job % 2;

// permute rows first 

    rperm(nrow, &a[1], &ja[1], &ia[1], &ao[1], &jao[1], &iao[1], &perm[1], &
	    locjob);

// then permute columns 

    locjob = 0;

    if (*job <= 2) {
	cperm(nrow, &ao[1], &jao[1], &iao[1], &ao[1], &jao[1], &iao[1], &
		perm[1], &locjob);
    } else {
	cperm(nrow, &ao[1], &jao[1], &iao[1], &ao[1], &jao[1], &iao[1], &
		qperm[1], &locjob);
    }

    return 0;
}

int CSR_Matrix::cperm(int nrow, double * a, int *ja, int * ia, double * ao,
			int * jao, int * iao, int * perm, int * job) {

    // System generated locals 
    int i__1;

    // Local variables 
    int i, k, nnz;

// -----------------------------------------------------------------------
 
// this subroutine permutes the columns of a matrix a, ja, ia. 
// the result id written in the output matrix  ao, jao, iao. 
// cperm computes B = A P, where  P is a permutation matrix 
// that maps column j into column perm(j), i.e., on return 
//      a(i,j) becomes a(i,perm(j)) in new matrix 
// Y. Saad, May 2, 1990 / modified Jan. 28, 1991. 
// -----------------------------------------------------------------------
 
// on entry: 
// ---------- 
// nrow 	= row dimension of the matrix 

// a, ja, ia = input matrix in csr format. 

// perm	= integer array of length ncol (number of columns of A 
//         containing the permutation array  the columns: 
//         a(i,j) in the original matrix becomes a(i,perm(j)) 
//         in the output matrix. 

// job	= integer indicating the work to be done: 
// 		job = 1	permute a, ja, ia into ao, jao, iao 
//                       (including the copying of real values ao and 
//                       the array iao). 
// 		job .ne. 1 :  ignore real values ao and ignore iao. 

// ------------ 
// on return: 
// ------------ 
// ao, jao, iao = input matrix in a, ja, ia format (array ao not needed) 


// Notes: 
// ------- 
// 1. if job=1 then ao, iao are not used. 
// 2. This routine is in place: ja, jao can be the same. 
// 3. If the matrix is initially sorted (by increasing column number) 
//    then ao,jao,iao  may not be on return. 

// ----------------------------------------------------------------------c
 
// local parameters: 

    // Parameter adjustments 
    --perm;
    --iao;
    --jao;
    --ao;
    --ia;
    --ja;
    --a;

    // Function Body 
    nnz = ia[nrow + 1] - 1;
    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	jao[k] = perm[ja[k]];
// L100: 
    }

//     done with ja array. return if no need to touch values. 

    if (*job != 1) {
	return 0;
    }

// else get new pointers -- and copy values too. 

    i__1 = nrow + 1;
    for (i = 1; i <= i__1; ++i) {
	iao[i] = ia[i];
// L1: 
    }

    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	ao[k] = a[k];
//  L2:      
    }

    return 0;

}

int CSR_Matrix::bfs(int n, int * ja, int * ia, int * nfirst, int * iperm,
		      int * mask, int * maskval, int * riord, int * levels, int * nlev) {

    // System generated locals 
    int i__1;

    // Local variables 
    int iend;
//    extern // Subroutine  int add_lvst__();
    int j, ii, istart;

//    static logical permut;
    int permut;

    int nod;

// -----------------------------------------------------------------------
 
// finds the level-structure (breadth-first-search or CMK) ordering for a 

// given sparse matrix. Uses add_lvst. Allows a  set of nodes to be 
// the initial level (instead of just one node). Allows masked nodes. 
// -------------------------parameters------------------------------------
 
// on entry: 
// ---------- 
// n      = number of nodes in the graph 
// ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data 

//          structure) 
// nfirst = number of nodes in the first level that is input in riord 
// iperm  = integer array indicating in which order to  traverse the graph
 
//          in order to generate all connected components. 
//          The nodes will be traversed in order iperm(1),....,iperm(n) 
//          Convention: 
//          if iperm(1) .eq. 0 on entry then BFS will traverse the 
//          nodes in the  order 1,2,...,n. 

// riord  = (also an ouput argument). on entry riord contains the labels 

//          of the nfirst nodes that constitute the first level. 

// mask   = array used to indicate whether or not a node should be 
//          condidered in the graph. see maskval. 
//          mask is also used as a marker of  visited nodes. 

// maskval= consider node i only when:  mask(i) .eq. maskval 
//          maskval must be .gt. 0. 
//          thus, to consider all nodes, take mask(1:n) = 1. 
//          maskval=1 (for example) 

// on return 
// --------- 
// mask   = on return mask is restored to its initial state. 
// riord  = `reverse permutation array'. Contains the labels of the nodes 

//          constituting all the levels found, from the first level to 
//          the last. 
// levels = pointer array for the level structure. If lev is a level 
//          number, and k1=levels(lev),k2=levels(lev+1)-1, then 
//          all the nodes of level number lev are: 
//          riord(k1),riord(k1+1),...,riord(k2) 
// nlev   = number of levels found 
// -----------------------------------------------------------------------
 
// Notes on possible usage 
// ------------------------- 
// 1. if you want a CMK ordering from a known node, say node init then 
//    call BFS with nfirst=1,iperm(1) =0, mask(1:n) =1, maskval =1, 
//    riord(1) = init. 
// 2. if you want the RCMK ordering and you have a preferred initial node 

//     then use above call followed by reversp(n,riord) 
// 3. Similarly to 1, and 2, but you know a good LEVEL SET to start from 

//    (nfirst = number if nodes in the level, riord(1:nfirst) contains 
//    the nodes. 
// 4. If you do not know how to select a good initial node in 1 and 2, 
//    then you should use perphn instead. 

// -----------------------------------------------------------------------
 
//     local variables -- 
    // Parameter adjustments 
    --levels;
    --riord;
    --mask;
    --iperm;
    --ia;
    --ja;

    // Function Body 
    permut = iperm[1] != 0;

//     start pointer structure to levels 

    *nlev = 0;

//     previous end 

    istart = 0;
    ii = 0;

//     current end 

    iend = *nfirst;

//     intialize masks to zero -- except nodes of first level -- 

    i__1 = *nfirst;
    for (j = 1; j <= i__1; ++j) {
	mask[riord[j]] = 0;
// L12: 
    }
// -----------------------------------------------------------------------
 
// L13: 

L1:
    ++(*nlev);
    levels[*nlev] = istart + 1;
    add_lvst(&istart, &iend, nlev, &riord[1], &ja[1], &ia[1], &mask[1], 
	    maskval);
    if (istart < iend) {
	goto L1;
    }
L2:
    ++ii;
    if (ii <= n) {
	nod = ii;

//	if (permut) {
        if (permut == 1) {
	    nod = iperm[nod];
	}
	if (mask[nod] == *maskval) {

//     start a new level 

	    istart = iend;
	    ++iend;
	    riord[iend] = nod;
	    mask[nod] = 0;
	    goto L1;
	} else {
	    goto L2;
	}
    }
// -----------------------------------------------------------------------
 
// L3: 
    levels[*nlev + 1] = iend + 1;
    i__1 = iend;
    for (j = 1; j <= i__1; ++j) {
	mask[riord[j]] = *maskval;
    }
    return 0;
} 

int CSR_Matrix::sort_irv(int * itmp, double *rtmp, int * n) {

    // System generated locals 
    int i__1, i__2;

    // Local variables 
    int jmin, i, j, it;
    double rt;
    int itmpmin;

//-----------------------------------------------------------------------
//   Tri selon une methode on ne peut plus simple */
//-----------------------------------------------------------------------
    // Parameter adjustments */
    --rtmp;
    --itmp;

    //      Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	itmpmin = itmp[i];
	jmin = i;
	i__2 = *n;
	for (j = i + 1; j <= i__2; ++j) {
	    if (itmp[j] < itmpmin) {
		jmin = j;
		itmpmin = itmp[jmin];
	    }
	}
	it = itmp[i];
	itmp[i] = itmp[jmin];
	itmp[jmin] = it;
	rt = rtmp[i];
	rtmp[i] = rtmp[jmin];
	rtmp[jmin] = rt;
    }
    return 0;
} 

void CSR_Matrix::rperm(int nrow,double * a,int * ja,int * ia, double * ao, 
			int * jao, int * iao, int * perm, int * job) {
 
   // System generated locals 
    int i__1, i__2;

    // Local variables 
    int i, j, k, ii, ko;

//    static logical values;
   
    int  values;

// -----------------------------------------------------------------------
 
// this subroutine permutes the rows of a matrix in CSR format. 
// rperm  computes B = P A  where P is a permutation matrix. 
// the permutation P is defined through the array perm: for each j, 
// perm(j) represents the destination row number of row number j. 
// Youcef Saad -- recoded Jan 28, 1991. 
// -----------------------------------------------------------------------
 
// on entry: 
// ---------- 
// n 	= dimension of the matrix 
// a, ja, ia = input matrix in csr format 
// perm 	= integer array of length nrow containing the permutation arrays 

// 	  for the rows: perm(i) is the destination of row i in the 
//         permuted matrix. 
//         ---> a(i,j) in the original matrix becomes a(perm(i),j) 
//         in the output  matrix. 

// job	= integer indicating the work to be done: 
// 		job = 1	permute a, ja, ia into ao, jao, iao 
//                       (including the copying of real values ao and 
//                       the array iao). 
// 		job .ne. 1 :  ignore real values. 
//                     (in which case arrays a and ao are not needed nor 
//                      used). 

// ------------ 
// on return: 
// ------------ 
// ao, jao, iao = input matrix in a, ja, ia format 
// note : 
//        if (job.ne.1)  then the arrays a and ao are not used. 
// ----------------------------------------------------------------------c
//           Y. Saad, May  2, 1990                                      c 
// ----------------------------------------------------------------------c
    // Parameter adjustments 
    --perm;
    --iao;
    --jao;
    --ao;
    --ia;
    --ja;
    --a;

    // Function Body 
    values = *job == 1;

//     determine pointers for output matix. 

    i__1 = nrow;
    for (j = 1; j <= i__1; ++j) {
	i = perm[j];
	iao[i + 1] = ia[j + 1] - ia[j];
// L50: 
    }

// get pointers from lengths 

    iao[1] = 1;
    i__1 = nrow;
    for (j = 1; j <= i__1; ++j) {
	iao[j + 1] += iao[j];
// L51: 
    }

// copying 

    i__1 = nrow;
    for (ii = 1; ii <= i__1; ++ii) {

// old row = ii  -- new row = iperm(ii) -- ko = new pointer 

	ko = iao[perm[ii]];
	i__2 = ia[ii + 1] - 1;
	for (k = ia[ii]; k <= i__2; ++k) {
	    jao[ko] = ja[k];
	    if (values==1) {
		ao[ko] = a[k];
	    }
	    ++ko;
// L60: 
	}
// L100: 
    }

    return;

}

int CSR_Matrix::add_lvst(int * istart, int * iend, int * nlev, int * riord,
			   int * ja, int * ia, int * mask, int * maskval) {

    // System generated locals
    int i__1, i__2;

    // Local variables
    int i, j, k, ir, nod;

// ---------------------------------------------------------------------- 
// adds one level set to the previous sets. span all nodes of previous
// set. Uses Mask to mark those already visited.
// -----------------------------------------------------------------------

    // Parameter adjustments 
    --mask;
    --ia;
    --ja;
    --riord;

    // Function Body
    nod = *iend;
    i__1 = *iend;
    for (ir = *istart + 1; ir <= i__1; ++ir) {
	i = riord[ir];
	i__2 = ia[i + 1] - 1;
	for (k = ia[i]; k <= i__2; ++k) {
	    j = ja[k];
	    if (mask[j] == *maskval) {
		++nod;
		mask[j] = 0;
		riord[nod] = j;
	    }
// L24: 
	}
// L25: 
    }
    *istart = *iend;
    *iend = nod;
    return 0;
//  -----------------------------------------------------------------------
 
}

void  CSR_Matrix::deblign(int nz , int *ptr , int *jptr , int *ai) {

  int i,ilign,indptr;

  ilign = 1;
  indptr = -1;

  jptr[0] = 1;
  for(i=1; i<nz; i++) {
    if (ai[i-1] < ai[i]) {
      jptr[ilign++]=i+1;
      ai[i-1] = 0;
    }    
    else{
      ai[i-1] = i+1;
    }
  }
  ai[nz-1] = 0;
}
/*
CompRow_Mat_double * CSR_Matrix::ConvertToIMLFormat(void) {

  printf("NbLines is %d\n", NbLines);
  if (NbLines == 0) {fprintf(stderr, "Matrix was not converted because
                                      its number of lines is zero\n");
		     assert(0);}

  
  PutInCompressedRow();

  int       nb_val = ReturnNumberOfNonZero();
  double * val_ptr = ReturnValPointer();
  int    * row_ptr = ReturnRowPointer(); 
  int    * col_ptr = ReturnColumIndexPointer();


//old
//  int i;
//  for (i = 0; i < NbLines + 1; i++)   row_ptr[i]--; 
//  for (i = 0; i < nb_val; i++)        col_ptr[i]--;

//new NICOLAS
    int* row_ptr_iml = new int[NbLines + 1];
    int* col_ptr_iml = new int[nb_val];
    int i;
    for (i = 0; i < NbLines + 1; i++)   row_ptr_iml[i] = row_ptr[i]-1; 
    for (i = 0; i < nb_val; i++)        col_ptr_iml[i] = col_ptr[i]-1;

#if DEBUGNIC>=1
   printf("Number of Non Zero %d\n",nb_val);
   printf("The row pointer is\n");
   for (i = 0; i < NbLines + 1; i++) {
     printf("%d ", row_ptr[i]);
   }
   printf("\n");
   printf("The col pointer and value are\n");
   for (i = 0; i < nb_val; i++) {
     printf("%d %12.5e\n", col_ptr[i], val_ptr[i]);
   }
#endif
  CompRow_Mat_double * mat;
//NICOLAS
  mat = new CompRow_Mat_double(NbLines, NbLines, nb_val, val_ptr, row_ptr_iml, col_ptr_iml);
  delete []row_ptr_iml;
  delete []col_ptr_iml;


   if (mat == NULL) {fprintf(stderr, "The matrix CSR was not converted properly to iml format\n"
                                     "the pointer is null\n");
                     assert(0);}
   return mat;
}
*/

CSR_Matrix& CSR_Matrix::operator=(const CSR_Matrix& in) {

  if (this != &in) {
//printf("Hellu1\n");
    Clear();

//printf("Hellu2\n");
//printf("Hellu3\n");
    MatrixWasReordered = in.MatrixWasReordered;
//printf("Hellu4\n");
    AllocationDone     = in.AllocationDone;
//printf("Hellu5\n");
    storage            = in.storage;
    NbLines            = in.NbLines;


    if (in.AllocationDone) {
//printf("Hellu6\n");
      //A.M.S.ILU_Exists = in.A.M.S.ILU_Exists;
//printf("Hellu10\n");
      a_    = List_Copy (in.a_);
//printf("Hellu11\n");
      ai_   = List_Copy (in.ai_);
//printf("Hellu12\n");
      ptr_  = List_Copy (in.ptr_);
//printf("Hellu13\n");
      jptr_ = List_Copy (in.jptr_);
//printf("Hellu14\n");
      
    }
    if (in.MatrixWasReordered) {
      
//  //printf("Hellu15\n");
//        permr_  = (int*) Malloc (NbLines * sizeof(int));
//  //printf("Hellu16\n");
//        rpermr_ = (int*) Malloc (NbLines * sizeof(int));
//  //printf("Hellu17\n");
//        permp_  = (int*) Malloc (2 * NbLines * sizeof(int));
//  //printf("Hellu18\n");
      
//printf("Hellu15\n");
  permr_  = new int[NbLines];
//printf("Hellu16\n");
  rpermr_ = new int[NbLines];
//printf("Hellu17\n");
  permp_  = new int[2 * NbLines];
//printf("Hellu18\n");

//gcc2.96
      memcpy(permr_, in.permr_, NbLines * sizeof(int));
      memcpy(rpermr_, in.rpermr_, NbLines * sizeof(int));
      memcpy(permp_, in.permp_, 2 * NbLines * sizeof(int));
      
    }
    
//    a_lu = in.a_lu;
//    A.M.S.jlu = in.A.M.S.jlu;
//    A.M.S.ju  = in.A.M.S.ju;
  }
  return *this;
}
  

CSR_Matrix::CSR_Matrix(const CSR_Matrix& in) {

  NbLines = in.NbLines;
  MatrixWasReordered = in.MatrixWasReordered;
  storage = in.storage;
  AllocationDone = in.AllocationDone;

  if (in.AllocationDone) {
    a_    = List_Copy (in.a_);
    ai_   = List_Copy (in.ai_);
    ptr_  = List_Copy (in.ptr_);
    jptr_ = List_Copy (in.jptr_);
  }
  
  if (in.MatrixWasReordered) {

//    permr_  = (int*) Malloc (NbLines * sizeof(int));
//    rpermr_ = (int*) Malloc (NbLines * sizeof(int));
//    permp_  = (int*) Malloc (2 * NbLines * sizeof(int));
  permr_  = new int[NbLines];
  rpermr_ = new int[NbLines];
  permp_  = new int[2 * NbLines];

    memcpy(permr_, in.permr_, NbLines * sizeof(int));
    memcpy(rpermr_, in.rpermr_, NbLines * sizeof(int));
    memcpy(permp_, in.permp_, 2 * NbLines * sizeof(int));

  }

//   a_lu = in.a_lu;
//   A.M.S.jlu = in.A.M.S.jlu;
//   A.M.S.ju  = in.A.M.S.ju;
    
}

void CSR_Matrix::OutputMatrixMatlabFormat(std::ostream &out) 
{
    double tol = 1E-15;
    double val;
    
    int i,j;
    const int Asize = NbLines;
    for(i=1;i<=Asize;i++){
      for(j=1;j<=Asize;j++){
	val = GetMatrix(i,j) ;
	if ( fabs(val) > tol )
	  out << "A( " << i << "," << j << ")=" << val <<"\n"; 
      }
    }
    return;
}



void CSR_Matrix::Load (const CSR_Matrix& in) {
  assert(0); //it is bugged I do not know why, God help me
  assert(AllocationDone && in.AllocationDone && NbLines == in.NbLines 
         && in.MatrixWasReordered == false);
  //we copy thr matrix
  printf("a_->nmax %d, in.a_->nmax %d\n", a_->nmax, in.a_->nmax);
  printf("ai_->nmax %d, in.ai_->nmax %d\n", ai_->nmax, in.ai_->nmax);
  printf("ptr_->nmax %d, in.ptr_->nmax %d\n", ptr_->nmax, in.ptr_->nmax);
  printf("jptr_->nmax %d, in.ajptr_->nmax %d\n", jptr_->nmax, in.jptr_->nmax);
  printf("a_->size %d, in.a_->size %d\n", a_->size, in.a_->size);
  printf("ai_->size %d, in.ai_->size %d\n", ai_->size, in.ai_->size);
  printf("ptr_->size %d, in.ptr_->size %d\n", ptr_->size, in.ptr_->size);
  printf("jptr_->size %d, in.jptr_->size %d\n", jptr_->size, in.jptr_->size);
  assert(a_->nmax == in.a_->nmax);
  assert(ai_->nmax == in.ai_->nmax);
  assert(ptr_->nmax == in.ptr_->nmax);
  assert(jptr_->nmax == in.jptr_->nmax);
  assert(a_->size == in.a_->size);
  assert(ai_->size == in.ai_->size);
  assert(ptr_->size == in.ptr_->size);
  assert(jptr_->size == in.jptr_->size);

  memcpy(a_->array,    in.a_->array,    a_->nmax * a_->size);
  memcpy(ai_->array,   in.ai_->array,   ai_->nmax * ai_->size);
  memcpy(ptr_->array,  in.ptr_->array,  ptr_->nmax * ptr_->size);
  memcpy(jptr_->array, in.jptr_->array, jptr_->nmax * jptr_->size);
  //a_    = List_Paste (in.a_); 
  //ai_   = List_Paste (in.ai_);
  //ptr_  = List_Paste (in.ptr_);
  //jptr_ = List_Paste (in.jptr_);

  MatrixWasReordered = false;
  delete []rpermr_; rpermr_ = NULL;
  //printf("Helly6\n");
  delete []permr_; permr_ = NULL;
  //printf("Helly7\n");
  delete []permp_; permp_ = NULL;
  //printf("Helly8\n");
}






