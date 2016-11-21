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
#ifndef H_scorecSListBase
#define H_scorecSListBase

struct scorecSLinkBase {
  scorecSLinkBase * next;
  scorecSLinkBase() { next = 0; }
  scorecSLinkBase(scorecSLinkBase *p) { next = p; }
};

class scorecSListBase{
protected:
  scorecSLinkBase *last; //last->next is head
  int Size;
public:
  scorecSListBase(void)
    : last(0), Size(0) {} 
  scorecSListBase(scorecSLinkBase *a)
    : last(a->next=a), Size(1) {}
  
  int size(void) const { return Size; }
  void insert(scorecSLinkBase *a);
  void append(scorecSLinkBase *a);
  
  scorecSLinkBase* get(void);
  
  void clear(void);
  scorecSLinkBase* nth(int i) const;
  void reverse();
  
  ~scorecSListBase();
  
  friend class scorecSListBaseIter;
};

class scorecSListBaseIter{
protected:
  scorecSListBase *cs;
  scorecSLinkBase *ce,*cl;
public:
  inline scorecSListBaseIter(scorecSListBase &s)
  { cs = &s; ce = cs->last; cl = 0;}
  inline void reset(void) { ce = cs->last; cl = 0; }
  void insert(scorecSLinkBase *);
  inline void insertBefore(scorecSLinkBase *);
  inline scorecSLinkBase* operator() ();
  scorecSLinkBase * next();
  void remove();
  scorecSListBase const * list() const;
};

inline scorecSLinkBase * scorecSListBaseIter::operator() ()
{
  scorecSLinkBase *ret = 0;
  if(ce)
    ret = (!(ce==cs->last && cl)) ? (cl=ce,ce=ce->next) : 0;
  return ret;
}

inline void scorecSListBaseIter::insertBefore(scorecSLinkBase *ne)
{
  if(cl){  // if there is a last element
    ne->next = ce;
    cl->next = ne;
    cs->Size++;
  } else // otherwise insert at start of list
    cs->insert(ne);
}

#endif
