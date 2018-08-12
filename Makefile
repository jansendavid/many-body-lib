#all:
#CXX= icc
CXX= icc
# PATHS

MKLPATH=${MKLROOT}/lib/intel64_lin
MKLINCLUDE=${MKLROOT}/include

INCS+=-I${MANYBODY}/src
# FLAGS

#FLAGS= -std=c++17
FLAGS+= -Wall 
#MKLLINGING
MKLLINK+= -I${MKLROOT}/include -DMKL_Complex16="std::complex<double>" -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
OBJECTS= main.cpp
# OTHER LIBS

TESTLIBS+= -lboost_unit_test_framework
#  DEBUG
ND+= -DNDEBUG -O3
D+= -g
main: $(OBJECTS) src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	g++ $(FLAGS) $(INCS) -o main $(OBJECTS) $(MKLLINK) $(ND)
#-static
test: testdir/test.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp
	$(CXX) $(D)  $(FLAGS) $(INCS) testdir/test.cpp -o test  $(MKLLINK) 
basistest: testdir/basistest.cpp src/basis.hpp src/accesfunctions.hpp src/numerics.hpp
	$(CXX) $(FLAGS) $(INCS) testdir/basistest.cpp -o basistest  $(MKLLINK) $(LIBS) $(TESTLIBS)

operatortest: testdir/basistest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	$(CXX) $(FLAGS) $(INCS) testdir/operatortest.cpp -o operatortest  $(MKLLINK) $(LIBS) $(TESTLIBS)

diagtest: testdir/diagtest.cpp src/basis.hpp  src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.hpp
	$(CXX) $(FLAGS) $(INCS) testdir/diagtest.cpp -o diagtest  $(MKLLINK) $(LIBS) $(TESTLIBS)

timeevtest: testdir/timeevtest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.hpp
	$(CXX) $(FLAGS) $(INCS) testdir/timeevtest.cpp -o timeevtest  $(MKLLINK) $(LIBS) $(TESTLIBS)
#-static

clean:
	rm *.o main test timeevtest operatortest diagtest basistest *.bin
