.PHONY:all

CXXFLAGS+=-ggdb
RUN:=valgrind

LD_LIBRARY_PATH+=:../../lib/x86_linux/
LATTICE_INCLUDES= -I/users/SCOREC/public/meshSim/latest/include/ -I/users/bfrantz/develop/scorecutil/include -I/users/bfrantz/develop/lattice-6.0-061119/include
LATTICE_LIBS= -L ./lib/x86_rhel3_gcc32 -L/users/SCOREC/public/meshSim/latest/lib/x86_rhel3_gcc32  -L/users/bfrantz/develop/lattice-6.0-061119/lib/x86_rhel3_gcc32 -lSimLatticeModel -lSimModel -lpthread

all: simmetrix_iterator_adaptor_unittest VolumeBuckets_unittest-run mVector_unittest-run mTensor2_unittest-run 

mVector_unittest-run: mVector_unittest
	./mVector_unittest

mVector_unittest: mVector_unittest.cpp ../include/mVector.h
	g++ -ggdb mVector_unittest.cpp -I ../include -o mVector_unittest

mTensor2_unittest: mTensor2_unittest.cpp
	g++ mTensor2_unittest.cpp -I ../include -o mTensor2_unittest  -L../../lib/x86_linux/ -lscorecutil_MeshSim  -L/users/SCOREC/public/meshSim/5.4-060520/lib/x86_rhel3_gcc32 -lSimMeshing -lSimMeshTools -lSimModel -lSimParasolid160 -lSimMeshing -lSimMeshTools -lSimModel -lpthread
# mTensor2_unittest-run: mTensor2_unittest
# 	./mTensor2_unittest

# VolumeBuckets_unittest-run: VolumeBuckets_unittest
# 	./VolumeBuckets_unittest

VolumeBuckets_unittest: VolumeBuckets_unittest.cpp ../include/VolumeBuckets.h
	g++ VolumeBuckets_unittest.cpp -I ../include -o VolumeBuckets_unittest


nested_iterator_unittest: nested_iterator_unittest.cpp ../include/beta/nested_iterator.h
	g++ $(CXXFLAGS) nested_iterator_unittest.cpp -I ../include -o nested_iterator_unittest


volume_brick_map_unittest: volume_brick_map_unittest.cpp ../include/beta/volume_brick_map.h   ../include/beta/nested_iterator.h 
	g++ $(CXXFLAGS) volume_brick_map_unittest.cpp -I ../include -o volume_brick_map_unittest
volume_brick_map_unittest-run: volume_brick_map_unittest
	$(RUN) ./volume_brick_map_unittest


simmetrix_iterator_adaptor_unittest: simmetrix_iterator_adaptor_unittest.cpp ../include/beta/simmetrix_iterator_adaptor.hpp    ../include/beta/nested_iterator.h 
	g++ $(CXXFLAGS) -finline-limit=0 $@.cpp  -I ../include $(LATTICE_INCLUDES) $(LATTICE_LIBS) -o $@

affine_space_unittest: affine_space_unittest.cpp ../include/affine_space.h
	g++ $(CXXFLAGS) $@.cpp  -I ../include -o $@

log_unittest: log_unittest.cpp ../include/Log.h 
	g++ $(CXXFLAGS) $@.cpp -I ../include -o $@ -L ../../scorecutil -lSCORECUtil



%-run: %
	$(RUN) ./$<
