CPPFLAGS=-I $(BOOST_INC) \
         -I $(CANVAS_INC) \
         -I $(CANVAS_ROOT_IO_INC) \
         -I $(CETLIB_INC) \
         -I $(CETLIB_EXCEPT_INC) \
         -I $(FHICLCPP_INC) \
         -I $(GALLERY_INC) \
         -I $(LARCOREOBJ_INC) \
         -I $(LARDATAOBJ_INC) \
         -I $(NUSIMDATA_INC) \
	 -I $(CLHEP_INC) \
         -I $(ROOT_INC) \
         -I $(LARCOREALG_INC) \

# Get root-config output and filter out -stdlib=libc++
#ROOT_CFLAGS := $(shell root-config --cflags | sed 's/-stdlib=libc++//g')
#CXXFLAGS += $(ROOT_CFLAGS)

LIBCORE_DIR=$(patsubst %/include,%/lib,$(ROOT_INC))
LIBCLHEP_DIR=$(patsubst %/include,%/lib,$(CLHEP_INC))

CXXFLAGS=$$(root-config --cflags --libs) -Wall -Werror -pedantic -Wno-unused-variable -Wno-unused-but-set-variable -I./inc/ 
CXX=g++
LDFLAGS=-L $(CANVAS_LIB) -l canvas \
        -L $(CANVAS_ROOT_IO_LIB) -l canvas_root_io \
	-L $(CETLIB_LIB) -l cetlib \
        -L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
        -L $(LARCOREALG_LIB) -l larcorealg_Geometry \
        -L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase \
	-L $(GALLERY_LIB) -l gallery \
        -L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
        -L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
	-L $(LARDATAOBJ_LIB) -l lardataobj_MCBase -l lardataobj_RawData -l lardataobj_AnalysisBase -lm -L $(LIBCORE_DIR) -l Core \
	-L $(LIBCORE_DIR) -lTree \
	-L $(LIBCORE_DIR) -lRIO \
	-L $(LIBCORE_DIR) -lHist
	#-L $(LIBCLHEP_DIR) -l CLHEP-Geometry-2.4.7.1 \
	#-L $(LIBCLHEP_DIR) -l CLHEP-Vector-2.4.7.1
	#-L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase -l lardataobj_MCBase -l lardataobj_RawData -l lardataobj_OpticalDetectorData -l lardataobj_AnalysisBase -lm
SOURCES=src/main.cpp src/Waveforms.cpp src/PointSignalFunctions.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
DEPENDENCIES=$(SOURCES:.cpp=.d)
EXECUTABLE=Ar39Study

$(EXECUTABLE): $(OBJECTS)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)
	@rm $(OBJECTS)

.cpp.o:
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@

clean:
	rm $(EXECUTABLE)
