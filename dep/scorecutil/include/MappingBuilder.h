#ifndef H_MappingBuilder
#define H_MappingBuilder

#include <MeshSim.h>

#include <string>
#include <vector>

#include "mPoint.h"

  class Mapping;
  class MappingBuilder 
    {
/*     protected: */
/*       bool isCubeRegular (pEntity ent); */
/*       { */
/* 	// WILL DO BETTER !!! */
/* 	//    for(int i=0;i<ent->size(0);i++) */
/* 	//  { */
/* 	//	mPoint p = ((mVertex*)(ent->get(0,i)))->point(); */
/* 	//	if(p(2) != 0.0)return false; */
/* 	// } */
/* 	return true; */
/*       } */
    public:
      virtual Mapping *BuildMapping (pEntity, int SystemOfCoordinates = 0) = 0;
      virtual std::string getMappingBuilderName() const = 0;
      virtual Mapping *BuildMappingByKnots(pEntity, int SystemOfCoordinates = 0, std::vector<mPoint> *knots = 0) {};
    };

#endif
