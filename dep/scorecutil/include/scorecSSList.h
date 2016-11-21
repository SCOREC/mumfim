/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/
#ifndef H_scorecSSList
#define H_scorecSSList

#include "scorecSListBase.h"

  {
namespace Util {


template<class T>
struct scorecSSLink : public scorecSLinkBase{
  T info;
  scorecSSLink(T a) : info(a) {}
};

template<class T> class scorecSSListCIter;
template<class T> class scorecSSListIter;

template<class T>
class scorecSSList : private scorecSListBase{
friend class scorecSSListCIter<T>;
friend class scorecSSListIter<T>;
public:
  scorecSSList() : scorecSListBase() {}
  scorecSSList(scorecSSLink<T> *a) : scorecSListBase(a) {}

  scorecSSList(const scorecSSList<T> &c);

  int size(void) const
    { return scorecSListBase::size(); }
  void insert(const scorecSSList<T> &al);
  void insert(const T& a)
    { scorecSListBase::insert( new scorecSSLink<T>(a)); }
  void append(const T& a)
    { scorecSListBase::append( new scorecSSLink<T>(a)); }
  void append(const scorecSSList<T> &al);
  void remove(const T &a);
  int inList(const T& a) const;
  void appendUnique(const T& a);
  void appendUnique(const scorecSSList<T> &al);
  T& nth(int i) const
    { return ((scorecSSLink<T>*)(scorecSListBase::nth(i)))->info;}
  void clear()
    { scorecSListBase::clear(); }
  void reverse()
    { scorecSListBase::reverse(); }
  
  void operator=(const scorecSSList<T> &l);
};

template<class T>
class scorecSSListCIter : private scorecSListBaseIter {
public:
  scorecSSListCIter(const scorecSSList<T>& s) : scorecSListBaseIter( (scorecSSList<T>&)s) {}

  inline int operator() (T &el)
  { 
    scorecSSLink<T> *lnk = (scorecSSLink<T>*) scorecSListBaseIter::operator() ();
    if(lnk)
      el = lnk->info;
    return (lnk!=0);
  }
  int next(T &el);
  inline void reset(void)
    { scorecSListBaseIter::reset(); }
  scorecSSList<T> const * list() const;
};


template<class T>
class scorecSSListIter : private scorecSListBaseIter {
public:
  scorecSSListIter(scorecSSList<T>& s) : scorecSListBaseIter(s) {}

  inline T* operator() ()
    { scorecSSLink<T> *lnk = (scorecSSLink<T>*) scorecSListBaseIter::operator() ();
      return lnk ? &lnk->info : 0; }
  inline int operator() (T &el)
    { scorecSSLink<T> *lnk = (scorecSSLink<T>*) scorecSListBaseIter::operator() ();
      if(lnk){
	el = lnk->info;
      }
      return (lnk!=0);
    }
  //return lnk ? el == (el = lnk->info) : 0; }
  inline void reset(void)
    { scorecSListBaseIter::reset(); }
  inline void insert(T ne)
    { scorecSListBaseIter::insert( new scorecSSLink<T>(ne) ); }
  inline void insertBefore(T ne)
    { scorecSListBaseIter::insertBefore( new scorecSSLink<T>(ne) ); }
  inline void remove()
    { scorecSListBaseIter::remove(); }
  void replaceCurrent( T ne );
};

template<class T>
scorecSSList<T>::scorecSSList(const scorecSSList<T> &c) : scorecSListBase()
{
  scorecSSLink<T> *ce = (scorecSSLink<T> *)(c.last);
  
  while(ce){
    ce = ce ? (scorecSSLink<T> *)(ce->next) : 0;
    append( ce->info );
    if(ce == (scorecSSLink<T> *)(c.last))
      ce = 0;
  }
}


template<class T>
void scorecSSList<T>::operator=(const scorecSSList<T> &c)
{
  clear();
  scorecSSLink<T> *ce = (scorecSSLink<T> *)(c.last);
  
  while(ce){
    ce = ce ? (scorecSSLink<T> *)(ce->next) : 0;
    append( ce->info );
    if(ce == (scorecSSLink<T> *)(c.last))
      ce = 0;
  }
}

template<class T>
int scorecSSList<T>::inList(const T& a) const
{
  scorecSSLink<T> *ce = (scorecSSLink<T> *)last;
  
  while(ce){
    ce = ce ? (scorecSSLink<T> *)(ce->next) : 0;
    if (ce->info == a)
      return 1;
    if(ce == (scorecSSLink<T> *)last)
      ce = 0;
  }
  return 0;
}

template<class T>
void scorecSSList<T>::remove(const T& remItem)
{
  scorecSSListIter<T> iter(*this);
  T a;
  while(iter(a)){
    if(a == remItem){
      iter.remove();
      break;
    }
  }
}

template<class T>
void scorecSSList<T>::insert(const scorecSSList<T> &al)
{
  scorecSSListCIter<T> aIter(al);
  scorecSSListIter<T> thisIter(*this);
  T item;
  while(aIter(item)){
    thisIter.insert(item);
    thisIter();
  }
}

template<class T>
void scorecSSList<T>::append(const scorecSSList<T> &al)
{
  scorecSSLink<T> *ce = (scorecSSLink<T> *)(al.last);
  
  while(ce){
    ce = ce ? (scorecSSLink<T> *)(ce->next) : 0;
    append(ce->info);
    if(ce == (scorecSSLink<T> *)(al.last))
      ce = 0;
  }
}

template<class T>
void scorecSSList<T>::appendUnique(const T& a)
{
  if( ! inList(a) )
    append(a);
}

template<class T>
void scorecSSList<T>::appendUnique(const scorecSSList<T> &al)
{
  scorecSSLink<T> *ce = (scorecSSLink<T> *)(al.last);
  
  while(ce){
    ce = ce ? (scorecSSLink<T> *)(ce->next) : 0;
    appendUnique(ce->info);
    if(ce == (scorecSSLink<T> *)(al.last))
      ce = 0;
  }
}

template<class T>
int scorecSSListCIter<T>::next(T &el)
{ 
  scorecSSLink<T> *lnk = (scorecSSLink<T>*)scorecSListBaseIter::next();
  if(lnk){
    el = lnk->info;
    return 1;
  } else
    return 0;
}

 template<class T>
 void scorecSSListIter<T>::replaceCurrent( T ne )
 { ((scorecSSLink<T>*)ce)->info = ne; }

template<class T>
scorecSSList<T> const * scorecSSListCIter<T>::list() const
{ return (scorecSSList<T> *)scorecSListBaseIter::list(); }


}} // End namespaces.
#endif
