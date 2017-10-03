CXX = g++
CXXFLAGS = -std=c++11 -fopenmp

all: reduction

reduction: main.o reduction.o input_output.o vector_operations.o
	$(CXX) $(CXXFLAGS) main.o reduction.o input_output.o vector_operations.o -o reduction

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

reduction.o: reduction.cpp
	$(CXX) $(CXXFLAGS) -c reduction.cpp

input_output.o: input_output.cpp
	$(CXX) $(CXXFLAGS) -c input_output.cpp

vector_operations.o: vector_operations.cpp
	$(CXX) $(CXXFLAGS) -c vector_operations.cpp

clean:
	rm -rf *.o reduction