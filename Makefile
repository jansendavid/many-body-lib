
#all:
#CXX= icc
IXX=g++
# PATHS

MKLPATH=${MKLROOT}/lib/intel64_lin
MKLINCLUDE=${MKLROOT}/include

INCS+=-I${MANYBODY}/src
INCS+=-I$(EIGEN)
# FLAGS

#FLAGS= -std=c++17
FLAGS+= -Wall -std=c++17
#MKLLINGING
MKLLINK32+=-DMKL_Complex16="std::complex<double>" -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

MKLLINK64+= -DMKL_Complex16="std::complex<double>"  -DMKL_ILP64 -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
OBJECTS= main.cpp
# OTHER LIBS
LIBSLINK+=-L${MANYBODY}/libs/
LIBS32+=-leigenmkl32
LIBS64+=-leigenmkl64
TESTLIBS+= -lboost_unit_test_framework
#  DEBUG
ND+= -DNDEBUG -O2
D+= -g
# executables
EXEC+=test timeevtest exvalcompare operatortest diagtest32 diagtest64 basistest basisbenchmark
diag32.o: src/diag.h src/diag.cpp
	$(IXX) $(ND)  $(FLAGS) $(INCS) src/diag.cpp -c -o diag32.o  $(MKLLINK32)
diag64.o: src/diag.h src/diag.cpp
	$(IXX) $(ND)  $(FLAGS) $(INCS) src/diag.cpp -c -o diag64.o  $(MKLLINK64)
# LIBS
libs/libeigenmkl64.a: diag64.o
	ar cr libs/libeigenmkl64.a diag64.o
libs/libeigenmkl32.a: diag32.o
	ar cr libs/libeigenmkl32.a diag32.o

main: $(OBJECTS) src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	$(IXX) $(FLAGS) $(INCS) main.cpp -o main  $(MKLLINK32) $(LIBSLINK) $(LIBS32) 
#-static
test: testdir/test.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp
	$(IXX) $(D)  $(FLAGS) $(INCS) testdir/test.cpp -o test  $(MKLLINK32) 
basistest: testdir/basistest.cpp src/basis.hpp src/accesfunctions.hpp src/numerics.hpp
	$(IXX) $(FLAGS) $(INCS) testdir/basistest.cpp -o basistest  $(MKLLINK32) $(LIBSLINK) $(TESTLIBS)

operatortest: testdir/basistest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	$(IXX) $(FLAGS) $(INCS) testdir/operatortest.cpp -o operatortest  $(MKLLINK32) $(LIBSLINK) $(TESTLIBS) $(LIBS32)

diagtest32: testdir/diagtest.cpp src/basis.hpp  src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h libs/libeigenmkl32.a
	$(IXX) $(FLAGS) $(INCS) testdir/diagtest.cpp -o diagtest32  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(TESTLIBS)

diagtest64: testdir/diagtest.cpp src/basis.hpp  src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h libs/libeigenmkl64.a
	$(IXX) $(FLAGS) $(INCS) testdir/diagtest.cpp -o diagtest64  $(MKLLINK64) $(LIBSLINK) $(LIBS64) $(TESTLIBS)
timeevtest: testdir/timeevtest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h
	$(IXX) $(FLAGS) $(INCS) testdir/timeevtest.cpp -o timeevtest  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(TESTLIBS)
reddmtest: testdir/reddmtest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(IXX) $(FLAGS) $(INCS) testdir/reddmtest.cpp -o reddmtest  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(TESTLIBS)
reddm: examples/reddm.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(IXX) $(FLAGS) $(INCS) examples/reddm.cpp -o reddm  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND)

exvalcompare: testdir/exvalcompare.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h
	$(IXX) $(FLAGS) $(INCS) testdir/exvalcompare.cpp -o exvalcompare  $(MKLLINK) $(LIBSLINK) $(LIBS32) $(TESTLIBS) -DMOM
#-static

basisbenchmark: benchmarkdir/basisbenchmark.cpp src/basis.hpp src/accesfunctions.hpp src/numerics.hpp
	$(IXX) $(FLAGS) $(INCS) benchmarkdir/basisbenchmark.cpp -o basisbenchmark  $(MKLLINK) $(LIBS)  $(ND) -g
clean:
	rm *.o main $(EXEC)
