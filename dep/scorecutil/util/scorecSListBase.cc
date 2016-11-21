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
#include "scorecSListBase.h"

scorecSListBase::~scorecSListBase()
{
  scorecSLinkBase* cur = 0;
  while(cur = scorecSListBase::get())
    delete cur;
}


void scorecSListBase::insert(scorecSLinkBase *a)
{
  if(last)
    a->next = last->next;
  else
    last = a;
  last->next = a;
  Size ++;
}

void scorecSListBase::append(scorecSLinkBase *a)
{
  if(last){
    a->next = last->next;
    last = last->next = a;
  } else
    last = a->next = a;
  Size ++;
}

scorecSLinkBase * scorecSListBase::get(void)
{
  if(last == 0)
    return 0;
  scorecSLinkBase *f = last->next;
  if(f == last)
    last = 0;
  else
    last->next = f->next;
  Size--;
  return f;
}

scorecSLinkBase * scorecSListBase::nth(int i) const
{
  scorecSLinkBase *ce = last->next;
  int n = 0;
  while(ce){
    if(n==i)
      return ce;
    ce = ce ? ce->next : 0;
    if(ce == last->next)
      ce = 0;
    n++;
  }
  return 0;
}

void scorecSListBase::reverse()
{
  scorecSListBase temp;
  scorecSLinkBase *ce;
  while(ce = get()){
    temp.insert(ce);
  }
  last = temp.last;
  Size = temp.Size;
  temp.last = 0;
  temp.Size = 0;
}


void scorecSListBaseIter::remove()
{
  if(cl){
    if(cs->Size==1){ // deleting only item in list
      delete ce;
      cs->last = ce = cl = 0;
    } else {
      if(ce == cs->last)
	cs->last = cl;
      cl->next = ce->next;
      delete ce;
      ce = cl;
      cl = 0;
    }
    cs->Size--;
  }
}

void scorecSListBase::clear()
{
  scorecSLinkBase *cur = 0;
  while(cur = scorecSListBase::get())
    delete cur;
}

void scorecSListBaseIter::insert(scorecSLinkBase *ne)
{
  if(!cs->last){ // empty list
    cs->last = ne;
    ne->next = ne;
    ce = ne;
    cl = ne;
  } else {
    if((!(ce==cs->last && cl))){
      ne->next = ce->next;
      ce->next = ne;
    } else {
      ne->next = cs->last->next;
      cs->last->next = ne;
      cs->last = ne;
    }
  }
  cs->Size++;
}

scorecSLinkBase * scorecSListBaseIter::next() // returns next without advancing iterator
{
  if(ce){
    scorecSLinkBase * ret = (!(ce==cs->last && cl)) ? ce->next : 0;
    return ret;
  }
  return 0;
}

scorecSListBase const * scorecSListBaseIter::list() const
{ return cs; }

