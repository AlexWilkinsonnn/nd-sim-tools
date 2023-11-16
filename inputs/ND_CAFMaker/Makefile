CXX = g++
CXXFLAGS = -g -Wall -fPIC -DNO_ART
ROOTFLAGS = `root-config --cflags --glibs`
INCLUDE = -I$(GENIE_INC)/GENIE

LDLIBS += -L$(LOG4CPP_LIB) -llog4cpp
LDLIBS += -L/usr/lib64 -lxml2
LDLIBS += -L$(PYTHIA6) -lPythia6
LDLIBS += -L$(ROOTSYS)/lib -lGeom -lEGPythia6
LDLIBS += -L$(GENIE)/lib \
                                        -lGAlgorithm \
                                        -lGBaryonResonance \
                                        -lGBase \
                                        -lGBodekYang \
                                        -lGCharm \
                                        -lGCoh \
                                        -lGCrossSections \
                                        -lGDecay \
                                        -lGDfrc \
                                        -lGDIS \
                                        -lGElas \
                                        -lGElFF \
                                        -lGEVGCore \
                                        -lGEVGDrivers \
                                        -lGEVGModules \
                                        -lGFluxDrivers \
                                        -lGFragmentation \
                                        -lGGeo \
                                        -lGGiBUU \
                                        -lGHadronTransp \
                                        -lGHEP \
                                        -lGInteraction \
                                        -lGLlewellynSmith \
                                        -lGMEC \
                                        -lGReinSehgal \
                                        -lGSingleKaon \
                                        -lGMessenger \
                                        -lGMuELoss \
                                        -lGNtuple \
                                        -lGNuclear \
                                        -lGNuE \
                                        -lGNuGamma \
                                        -lGNumerical \
                                        -lGPDF \
                                        -lGPDG \
                                        -lGQEL \
                                        -lGQPM \
                                        -lGRegistry \
                                        -lGRES \
                                        -lGUtils \
                                        -lGReWeight

# make a binary for every .cxx file
all : $(patsubst %.cxx, %.o, $(wildcard *.cxx))

# rule for each target
%.o : %.cxx
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -o $*.o $(LDLIBS) -c $*.cxx #compile
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -o $* $*.o        #link

clean:
	rm -f $(wildcard *.o) $(patsubst %.cxx, %, $(wildcard *.cxx))
	rm -f $(wildcard AutoDict_*)
