#ifndef _LISTMAN_H_
#define _LISTMAN_H_


typedef struct {
  int nmax;
  int size;
  int incr;
  int n;
  int isorder;
  char *array;
} List_T;

void *Malloc(size_t size);
void Free(void *ptr);

List_T *List_Create(int n, int incr, int size);
void    List_Delete(List_T *liste);
void    List_Realloc(List_T *liste,int n);
void    List_Add(List_T *liste, void *data);
int     List_Nbr(List_T *liste);
int     List_Nbr0(List_T *liste);
void    List_Insert(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
int     List_Replace(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
void    List_Read(List_T *liste, int index, void *data);
void    List_Write(List_T *liste, int index, void *data);
void    List_Put(List_T *liste, int index, void *data);
void    List_Pop(List_T *liste);
void   *List_Pointer(List_T *liste, int index);
void    List_Tri(List_T *liste, int (*fcmp)(const void *a, const void *b));
int     List_Search(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
int     List_ISearch(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
int     List_Query(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
int     List_LQuery(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b), int first);
void   *List_PQuery(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
int     List_Suppress(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
void    List_Reset(List_T *liste);
//Added by Nicolas to define a copy constructeur of CSR_Matrix
//Copy allocate and paste
List_T *List_Copy(List_T *in);
//this obe just paste
//List_T *List_Paste(List_T *in);
#endif




