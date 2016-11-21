
BUILD := latest
PLATFORM := x64_rhel6_gcc44
FNAME := partitionMesh
PARALLEL := openmpi16

EXECNAME := bin/$(PLATFORM)/$(FNAME)

PROJECTROOT := /net/common/meshSim
INCLUDES := -I$(PROJECTROOT)/$(BUILD)/include
LIBPATH := -L$(PROJECTROOT)/$(BUILD)/lib/$(PLATFORM) -L$(PROJECTROOT)/$(BUILD)/lib/$(PLATFORM)/psKrnl
LIBS = $(LIBPATH) -lSimPartitionedMesh-mpi -lSimMeshing -lSimMeshTools -lSimModel \
  -lSimPartitionWrapper-$(PARALLEL) -lSimParasolid270 -lpskernel

COMPFLAG := -c
OBJNAME := obj/$(PLATFORM)/$(FNAME).o
OBJTARGET := -o $(OBJNAME)
LINKTARGET := -o $(EXECNAME)

ifeq ($(PLATFORM),x64_rhel6_gcc44)
  # May need to change these to get correct MPI implementation, paths etc
  PARALLEL := openmpi16
  CC := /usr/local/openmpi/1.6.5-ib/bin/mpicxx
  LINK := $(CC)
  CCFLAGS := -O2
  LINKFLAGS :=
  LIBS := $(LIBS) -lpthread -lm
endif

$(EXECNAME): $(OBJNAME)
	-mkdir -p bin/$(PLATFORM)
	$(LINK) $(LINKFLAGS) $(OBJNAME) $(LINKTARGET) $(LIBS)

$(OBJNAME): $(FNAME).cc
	-mkdir -p obj/$(PLATFORM)
	$(CC) $(CCFLAGS) $(INCLUDES) $(COMPFLAG) $(FNAME).cc $(OBJTARGET)

clean:
	rm -f $(EXECNAME) $(OBJNAME)
