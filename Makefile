
# C++ compiler 
CXX = g++

# Standard flags for C++ 
CXXFLAGS ?= -I/usr/local/include/openbabel-2.0 -g

# Standard preprocessor flags (common for CC and CXX) 
CFLAGS ?= -I/usr/local/include/openbabel-2.0 -g -std=c++11

# Standard linker flags 
LDFLAGS ?= -L/usr/local/lib -lopenbabel -pthread -g

# -------------------------------------------------------------------------
# Do not modify the rest of this file!
# -------------------------------------------------------------------------

### Variables: ###

CPPDEPS ?= -MT$@ -MF`echo $@ | sed -e 's,\.o$$,.d,'` -MD -MP
MULTIPLE_CXXFLAGS ?= -I. $(CPPFLAGS) $(CXXFLAGS)


MULTIPLE_OBJECTS =  \
	     dzme_mathlib.o \
	     dzme_result.o \
             dzme.o 

### Conditionally set variables: ###



### Targets: ###

all: dzme

install: 

uninstall: 

clean: 
	rm -f ./*.o
	rm -f g/*.d
	rm -f dzme

dzme: $(MULTIPLE_OBJECTS)
	$(CXX) -o $@ $(MULTIPLE_OBJECTS)  $(MULTIPLE_CXXFLAGS) $(LDFLAGS)

dzme.o: ./dzme.cpp
	$(CXX) -c -o $@ $(MULTIPLE_CXXFLAGS) $(CPPDEPS) $<

dzme_result.o: ./result.cpp
	$(CXX) -c -o $@ $(MULTIPLE_CXXFLAGS) $(CPPDEPS) $<
	
dzme_mathlib.o: ./mathlib.cpp
	$(CXX) -c -o $@ $(MULTIPLE_CXXFLAGS) $(CPPDEPS) $<


#.PHONY: all install uninstall clean


# Dependencies tracking:
-include ./*.d
