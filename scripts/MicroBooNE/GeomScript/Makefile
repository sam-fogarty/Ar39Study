CPPFLAGS=-I $(BOOST_INC) \
         -I $(CETLIB_INC) \
         -I $(CETLIB_EXCEPT_INC) \

CXXFLAGS=-std=c++14 -Wall -Werror -pedantic -Wno-unused-variable -Wno-unused-but-set-variable
CXX=g++

GeomScript: GeomScript.cpp
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $<
	@printf "Successfully compiled GeomScript!\n\n"
	@echo GeomScript requires that the ChannelWireGeometry_V2.txt file has the top text chopped off.
	@echo The file should only be a table of numbers -- no text.
	@echo Record intersections with ./GeomScript \>\> CollectionIntersections.txt

all: GeomScript

clean:
	rm *.o GeomScript
