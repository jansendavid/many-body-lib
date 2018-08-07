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
MKLLINK+= -DMKL_Complex16="std::complex<double>" -DMKL_ILP64 -m64 -I${MKLROOT}/include -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
OBJECTS= main.cpp
# OTHER LIBS

TESTLIBS+= -lboost_unit_test_framework
#  DEBUG
ND=+ -DNDEBUG -O3
D+= -g
main: $(OBJECTS) src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	g++ $(FLAGS) $(INCS) -o main $(OBJECTS) $(MKLLINK)
#-static

basistest: testdir/basistest.cpp src/basis.hpp src/accesfunctions.hpp src/numerics.hpp
	$(CXX) $(FLAGS) $(INCS) testdir/basistest.cpp -o basistest  $(MKLLINK) $(LIBS) $(TESTLIBS)

operatortest: testdir/basistest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	$(CXX) $(FLAGS) $(INCS) testdir/operatortest.cpp -o operatortest  $(MKLLINK) $(LIBS) $(TESTLIBS)

timeevtest: testdir/timeevtest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp
	$(CXX) $(FLAGS) $(INCS) testdir/timeevtest.cpp -o timeevtest  $(MKLLINK) $(LIBS) $(TESTLIBS)
#-static

clean:
	rm *.o main test *.bin
