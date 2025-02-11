CXXFLAGS=-O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
CXX=g++

all: main.o tridiagonal_matrix.o systems_builder.o solver.o thomas_method.o functions.o arguments.o
	$(CXX) main.o tridiagonal_matrix.o systems_builder.o solver.o functions.o thomas_method.o arguments.o -o a
	
main.o: main.cpp arguments.h enums.h
	$(CXX) -c $(CXXFLAGS) -o main.o main.cpp

tridiagonal_matrix.o: tridiagonal_matrix.cpp tridiagonal_matrix.h arguments.h functions.h enums.h
	$(CXX) -c $(CXXFLAGS) -o tridiagonal_matrix.o tridiagonal_matrix.cpp

systems_builder.o: systems_builder.cpp systems_builder.h enums.h
	$(CXX) -c $(CXXFLAGS) -o systems_builder.o systems_builder.cpp
	
solver.o: solver.cpp solver.h enums.h tridiagonal_matrix.h systems_builder.h thomas_method.h
	$(CXX) -c $(CXXFLAGS) -o solver.o solver.cpp
	
thomas_method.o: thomas_method.cpp thomas_method.h enums.h tridiagonal_matrix.h
	$(CXX) -c $(CXXFLAGS) -o thomas_method.o thomas_method.cpp
	
functions.o: functions.cpp functions.h enums.h
	$(CXX) -c $(CXXFLAGS) -o functions.o functions.cpp
	
arguments.o: arguments.cpp arguments.h enums.h
	$(CXX) -c $(CXXFLAGS) -o arguments.o arguments.cpp
	
	
clean:
	rm -f *.o *.out
