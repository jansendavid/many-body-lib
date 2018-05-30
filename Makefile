all:
CXX= g++
FLAGS= -std=c++11
FLAGS+=-g -Wall 
CFLAGS+= -g -Wall 
CFLAGS+= std=c++11 
OBJECTS= main.cpp
LIBS+= -larmadillo
TESTLIBS+= -lboost_unit_test_framework

main: $(OBJECTS) basis.hpp operators.hpp
	$(CXX) $(FLAGS) -o main $(OBJECTS)  $(LIBS) $(TESTLIBS)
#-static

test: test.cpp basis.hpp operators.hpp
	$(CXX) $(FLAGS) -o test test.cpp  $(LIBS) $(TESTLIBS)
#-static

clean:
	rm *.o main test *.bin
