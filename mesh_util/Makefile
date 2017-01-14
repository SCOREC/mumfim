
BUILD := 10.0-160902
PLATFORM := x64_rhel7_gcc48
FNAME := partitionMesh
PARALLEL := openmpi110
CC := /fasttmp/wtobin/develop/install/openmpi/1.10.0/bin/mpicxx
LINK := $(CC)
#  CCFLAGS := -O2
CCFLAGS := -g
LINKFLAGS :=
LIBS := $(LIBS) -lpthread -lm


EXECNAME := bin/$(PLATFORM)/$(FNAME)

PROJECTROOT := /net/common/meshSim
INCLUDES := -I$(PROJECTROOT)/$(BUILD)/include
LIBPATH := -L$(PROJECTROOT)/$(BUILD)/lib/$(PLATFORM) -L$(PROJECTROOT)/$(BUILD)/lib/$(PLATFORM)/psKrnl
LIBS = $(LIBPATH) -lSimPartitionedMesh-mpi -lSimMeshing -lSimMeshTools -lSimModel \
  -lSimPartitionWrapper-$(PARALLEL) -lSimParasolid290 -lpskernel

COMPFLAG := -c
OBJNAME := obj/$(PLATFORM)/$(FNAME).o
OBJTARGET := -o $(OBJNAME)
LINKTARGET := -o $(EXECNAME)

$(EXECNAME): $(OBJNAME)
	-mkdir -p bin/$(PLATFORM)
	$(LINK) $(LINKFLAGS) $(OBJNAME) $(LINKTARGET) $(LIBS)

$(OBJNAME): $(FNAME).cc
	-mkdir -p obj/$(PLATFORM)
	$(CC) $(CCFLAGS) $(INCLUDES) $(COMPFLAG) $(FNAME).cc $(OBJTARGET)

clean:
	rm -f $(EXECNAME) $(OBJNAME)
