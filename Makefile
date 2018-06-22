all:
CXX= g++
FLAGS= -std=c++14 
FLAGS+=-g -Wall 
CFLAGS+= -g -Wall 
CFLAGS+= std=c++11 
OBJECTS= main.cpp
LIBS+= -larmadillo
TESTLIBS+= -lboost_unit_test_framework

main: $(OBJECTS) basis.hpp operators.hpp
	$(CXX) $(FLAGS) -o main $(OBJECTS)  $(LIBS) 
#-static

test: test.cpp basis.hpp operators.hpp accesfunctions.hpp
	$(CXX) $(FLAGS) -o test test.cpp  $(LIBS) $(TESTLIBS)
#-static

clean:
	rm *.o main test *.bin
