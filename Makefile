src = $(wildcard src/*.c) 
ccsrc = $(wildcard src/*.cpp)
obj = $(src:.c=.o) $(ccsrc:.cpp=.o)

#CXXFLAGS =-D USE_EIGEN -D NDEBUG -O3 -fopenmp -std=gnu++11 
CXXFLAGS =-D USE_EIGEN -I/home/strebdom/git/spectra/include/Spectra -I/usr/include/eigen3 -D NDEBUG -O3 -mtune=native -fopenmp -std=gnu++11

LDFLAGS =-g -O3 -fopenmp -std=gnu++11 -lgsl -lgslcblas -lnlopt 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

a.out: $(obj)
	$(CXX) -o $@ $^ $(LDFLAGS)

all: a.out

.PHONY: clean
clean:
	rm -f $(obj) a.out

