/* Auteur: Marc UME  **  Date: 01/01/94  **  Commentaire: Rien */

#include <string>
#include <stdlib.h>
#include <search.h>
#include <sys/types.h>
#include "listman.h"
#include <cstring>

void *Malloc(size_t size)
{
  void *ptr;

  if (!size) return(NULL);
  ptr = malloc(size);
  return(ptr);
}

void *Realloc(void *ptr, size_t size)
{
  if (!size) return(NULL);
  ptr = realloc(ptr,size);
  return(ptr);
}

void Free(void *ptr)
{
//printf("Hellt1\n");
  if (ptr == NULL) return;
//printf("Hellt2\n");
  free(ptr);
//printf("Hellt3\n");
  ptr = NULL;
//printf("Hellt4\n");
}

/* - List ---------------------------------------------------------- */

List_T *List_Create(int n, int incr, int size)
{
  List_T *liste;

  if (n <= 0)  n = 1 ;
  if (incr <= 0) incr = 1;

  liste = (List_T *)Malloc(sizeof(List_T));

  liste->nmax    = 0;
  liste->incr    = incr;
  liste->size    = size;
  liste->n       = 0;
  liste->isorder = 0;
  liste->array   = NULL;

  List_Realloc(liste,n);
  return(liste);
}

void List_Delete(List_T *liste)
{
//printf("Hella0\n");
  if (liste != 0) {
//printf("Hella1\n");
//if (liste->array) printf("liste->array pas nul\n");
//else              printf("liste->array nul\n");
    Free(liste->array);
//printf("Hella2\n");
    Free(liste);
//printf("Hella3\n");
  }
}

void List_Realloc(List_T *liste,int n)
{
  char* temp;
  if (n <= 0) return;
  if (liste->array == NULL) {
    liste->nmax = ((n - 1) / liste->incr + 1) * liste->incr;
/*if (liste->nmax * liste->size == 0) printf("dans realloc probleme avec array null\n");*/
    liste->array = (char *)Malloc(liste->nmax * liste->size);
  }
  else {
    if (n > liste->nmax) {
      liste->nmax = ((n - 1) / liste->incr + 1) * liste->incr;
/*      liste->array = (char *)Realloc(liste->array, liste->nmax * liste->size); */
/*if (liste->nmax * liste->size == 0) printf("dans realloc probleme second case\n");*/
      temp = (char *)Realloc(liste->array, liste->nmax * liste->size);
      liste->array = temp;
    }
  }
}

void List_Add(List_T *liste, void *data)
{
  liste->n++;

  List_Realloc(liste,liste->n);
  liste->isorder = 0;
  memcpy(&liste->array[(liste->n - 1) * liste->size],data,liste->size);
}

int List_Nbr(List_T *liste)
{
  return(liste->n);
}

int List_Nbr0(List_T *liste)
{
  return (liste)? liste->n : 0 ;
}

void List_Insert(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  if (List_Search(liste,data,fcmp) == 0)
    List_Add(liste,data);
}

int List_Replace(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) List_Tri(liste,fcmp);
  liste->isorder = 1;
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  if (ptr == NULL) {
    List_Add(liste,data);
    return(0);
  }
  else {
    memcpy(ptr,data,liste->size);
    return (1);
  }
}

void List_Read(List_T *liste, int index, void *data)
{
  memcpy(data,&liste->array[index * liste->size],liste->size);
}

void List_Write(List_T *liste, int index, void *data)
{
  liste->isorder = 0;
  memcpy(&liste->array[index * liste->size],data,liste->size);
}

void List_Put(List_T *liste, int index, void *data)
{
  if (index >= liste->n) {
    liste->n = index + 1;
    List_Realloc(liste,liste->n);
    List_Write(liste,index,data);
  } else {
    List_Write(liste,index,data);
  }
}

void List_Pop(List_T *liste)
{
  liste->n -- ;
}

void *List_Pointer(List_T *liste, int index)
{
/*  liste->isorder = 0; */
  return(&liste->array[index * liste->size]);
}

void *List_Pointer_NoChange(List_T *liste, int index)
{
  return(&liste->array[index * liste->size]);
}

void List_Tri(List_T *liste, int (*fcmp)(const void *a, const void *b))
{
  qsort(liste->array,liste->n,liste->size,fcmp);
}

int List_Search(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) { List_Tri(liste,fcmp) ; liste->isorder = 1 ; }
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  if (ptr == NULL) return(0);
  return (1);
}

int List_ISearch(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) List_Tri(liste,fcmp);
  liste->isorder = 1;
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  if (ptr == NULL) return(-1);
  return (((long)ptr - (long)liste->array) / liste->size);
}

int List_Query(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) List_Tri(liste,fcmp);
  liste->isorder = 1;
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  if (ptr == NULL) return(0);

  memcpy(data,ptr,liste->size);
  return (1);
}


int List_LQuery(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b), int first)
{
  char *ptr;  
  char *startptr;
  startptr =  ptr + liste->size;
  if ( startptr >= ( liste->array + liste->n * liste->size))
    startptr = NULL;
  memcpy(data,ptr,liste->size);
  return (1);
}

void *List_PQuery(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) List_Tri(liste,fcmp);
  liste->isorder = 1;
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  return(ptr);
}

int List_Suppress(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  char *ptr;
  int len;
  
  ptr = (char*)List_PQuery(liste,data,fcmp) ;
  if (ptr == NULL) return(0);
  
  liste->n--;
  len = liste->n - (((long)ptr - (long)liste->array) / liste->size);
  if (len > 0) memmove(ptr, ptr + liste->size, len * liste->size);
  return(1);
}

void List_Reset(List_T *liste)
{
  liste->n = 0;
}





//Added by Nicolas to define a copy constructeur of CSR_Matrix

List_T *List_Copy(List_T *in) {
  List_T * liste;
  if (in == 0) liste = 0;
  else {
    liste = (List_T *)Malloc(sizeof(List_T));
    
    liste->nmax    = in->nmax;
    liste->incr    = in->incr;
    liste->size    = in->size;
    liste->n       = in->n;
    liste->isorder = in->isorder;
    liste->array   = NULL;
    
    List_Realloc(liste,liste->nmax);
    
    memcpy(liste->array, in->array, liste->nmax * liste->size);
  }

  return liste;
}


