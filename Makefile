#
CXX = g++
CXXFLAGS = -Wextra -g3 -O3 -std=c++14
#
Plasma: Plasma.o Plasma_funcs.o FFT.o
	$(CXX) $(CXXFLAGS) -o Plasma Plasma.o Plasma_funcs.o FFT.o -lfftw3
clean:
	rm *.o
