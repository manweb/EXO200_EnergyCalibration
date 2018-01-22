CXX = g++
CXXFLAGS = -g -O2 -Wall -fPIC

# --- ROOT --------------------------------------------------------------
CXXFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs) -lXMLParser -lThread -lSpectrum

# --- EXOSoftware -------------------------------------------------------
CXXFLAGS +=-I/nfs/slac/g/exo/software/hudson/builds-rhel5/current/include
EXOLIBS = -L/nfs/slac/g/exo/software/hudson/builds-rhel5/current/lib -lEXOUtilities

all: GetData

GetData: GetData.o
	@ echo "Linking $@..."
	@ ${CXX} ${CXXFLAGS} -o $@ $^ ${ROOTLIBS} ${EXOLIBS}

GetData.o: GetData.cc
	@ ${CXX} ${CXXFLAGS} -c $^ -o $@

.PHONY : clean

clean:
	@echo "Cleaning up..."
	@ rm -f GetData GetData.o
