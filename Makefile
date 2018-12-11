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
MKLLINKLIB32+= -DMKL_Complex16="std::complex<double>" -I${MKLROOT}/include  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
MKLLINKLIB64+= -DMKL_Complex16="std::complex<double>" -DMKL_ILP64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
MKLLINKPROG32+=  -I${MKLROOT}/include  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
MKLLINKPROG64+= -DMKL_ILP64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
OBJECTS= main.cpp
# OTHER LIBS
LIBSLINK+=-L${MANYBODY}/libs/
LIBS32+=-leigenmkl32
LIBS64+=-leigenmkl64
TESTLIBS+= -lboost_unit_test_framework
#  DEBUG
ND+= -DNDEBUG -O3
D+= -g
# executables
EXEC+=test timeevtest exvalcompare operatortest diagtest32 diagtest64 basistest basisbenchmark
diag32.o: src/diag.h src/diag.cpp
	$(CXX) $(ND)  $(FLAGS) $(INCS) src/diag.cpp -c -o diag32.o  $(MKLLINKLIB32)
diag64.o: src/diag.h src/diag.cpp
	$(CXX) $(ND)  $(FLAGS) $(INCS) src/diag.cpp -c -o diag64.o  $(MKLLINKLIB64)
# LIBS
libs/libeigenmkl64.a: diag64.o
	ar cr libs/libeigenmkl64.a diag64.o
libs/libeigenmkl32.a: diag32.o
	ar cr libs/libeigenmkl32.a diag32.o

main: $(OBJECTS) src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	$(CXX) $(FLAGS) $(INCS) main.cpp -o main  $(MKLLINKPROG32) $(LIBSLINK) $(LIBS32) 
#-static
test: testdir/test.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp
	$(CXX) $(D)  $(FLAGS) $(INCS) testdir/test.cpp -o test  $(MKLLINKPROG32) 
basistest: testdir/basistest.cpp src/basis.hpp src/accesfunctions.hpp src/numerics.hpp
	$(CXX) $(FLAGS) $(INCS) testdir/basistest.cpp -o basistest  $(MKLLINKPROG32) $(LIBSLINK) $(TESTLIBS)

operatortest: testdir/basistest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	$(CXX) $(FLAGS) $(INCS) testdir/operatortest.cpp -o operatortest  $(MKLLINKPROG32+) $(LIBSLINK) $(TESTLIBS)

diagtest32: testdir/diagtest.cpp src/basis.hpp  src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h libs/libeigenmkl32.a
	$(CXX) $(FLAGS) $(INCS) testdir/diagtest.cpp -o diagtest32  $(MKLLINKPROG32) $(LIBSLINK) $(LIBS32) $(TESTLIBS)

diagtest64: testdir/diagtest.cpp src/basis.hpp  src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h libs/libeigenmkl64.a
	$(CXX) $(FLAGS) $(INCS) testdir/diagtest.cpp -o diagtest64  $(MKLLINKPROG64) $(LIBSLINK) $(LIBS64) $(TESTLIBS)
timeevtest: testdir/timeevtest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) testdir/timeevtest.cpp -o timeevtest  $(MKLLINKPROG32) $(LIBSLINK) $(LIBS32) $(TESTLIBS)

exvalcompare: testdir/exvalcompare.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) testdir/exvalcompare.cpp -o exvalcompare  $(MKLLINK) $(LIBSLINK) $(LIBS32) $(TESTLIBS) -DMOM
#-static

basisbenchmark: benchmarkdir/basisbenchmark.cpp src/basis.hpp src/accesfunctions.hpp src/numerics.hpp
	$(CXX) $(FLAGS) $(INCS) benchmarkdir/basisbenchmark.cpp -o basisbenchmark  $(MKLLINK) $(LIBS)  $(ND) -g
clean:
	rm *.o main $(EXEC)
