
#all:
#CXX= icc

IMPI=mpic++ -std=c++17

IXX=g++

# PATHS
MPILINK32= -DMKL_Complex16="std::complex<double>" -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -lboost_serialization -lboost_mpi  -I$MPI_INCLUDE -L$MPI_LIB -lmpi_cxx
#BOOST_INCLUDE_LINE=-I/usr/include/boost/mpi



MKLPATH=${MKLROOT}/lib/intel64_lin
MKLINCLUDE=${MKLROOT}/include

INCS+=-I${MANYBODY}/src
INCS+=-I${EIGEN}
INCS+=-I${EINC}
# FLAGS

#FLAGS= -std=c++17
FLAGS+= -Wall -std=c++17
#MKLLINGING
MKLLINK32+=-DMKL_Complex16="std::complex<double>" -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

MKLLINK64+= -DMKL_Complex16="std::complex<double>"  -DMKL_ILP64 -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
OBJECTS= main.cpp
# OTHER LIBS
LIBSLINK+=-L${MANYBODY}/libs
LIBSLINK+=-L${ELIBS}
LIBS32+=-leigenmkl32 -lboost_program_options
LIBS64+=-leigenmkl64 -lboost_program_options
TESTLIBS+= -lboost_unit_test_framework
#  DEBUG
ND+= -DNDEBUG
OP+=-O2
D+= -g
# executables
EXEC+=test timeevtest exvalcompare operatortest diagtest32 diagtest64 basistest basisbenchmark
diag32.o: src/diag.h src/diag.cpp
	$(IXX) $(ND)  $(FLAGS) $(INCS) src/diag.cpp -c -o diag32.o  $(MKLLINK32) $(OP)
diag64.o: src/diag.h src/diag.cpp
	$(IXX) $(ND)  $(FLAGS) $(INCS) src/diag.cpp -c -o diag64.o  $(MKLLINK64) $(OP)
# LIBS
libs/libeigenmkl64.a: diag64.o
	ar cr libs/libeigenmkl64.a diag64.o
libs/libeigenmkl32.a: diag32.o
	ar cr libs/libeigenmkl32.a diag32.o

main: $(OBJECTS) src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	$(IXX) $(FLAGS) $(INCS) main.cpp -o main  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(OP)
#-static
test: testdir/test.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp
	$(IXX) $(D)  $(FLAGS) $(INCS) testdir/test.cpp -o test  $(MKLLINK32) $(OP)
basistest: testdir/basistest.cpp src/basis.hpp src/accesfunctions.hpp src/numerics.hpp
	$(IXX) $(FLAGS) $(INCS) testdir/basistest.cpp -o basistest  $(MKLLINK32) $(LIBSLINK) $(TESTLIBS) $(OP)

operatortest: testdir/basistest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp
	$(IXX) $(FLAGS) $(INCS) testdir/operatortest.cpp -o operatortest  $(MKLLINK32) $(LIBSLINK) $(TESTLIBS) $(LIBS32)  $(OP)

diagtest32: testdir/diagtest.cpp src/basis.hpp  src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h libs/libeigenmkl32.a
	$(IXX) $(FLAGS) $(INCS) testdir/diagtest.cpp -o diagtest32  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(TESTLIBS) $(OP)

diagtest64: testdir/diagtest.cpp src/basis.hpp  src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h libs/libeigenmkl64.a
	$(IXX) $(FLAGS) $(INCS) testdir/diagtest.cpp -o diagtest64  $(MKLLINK64) $(LIBSLINK) $(LIBS64) $(TESTLIBS) $(OP)
timeevtest: testdir/timeevtest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) testdir/timeevtest.cpp -o timeevtest  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(TESTLIBS) $(OP)
FTLMtest: testdir/FTLMtest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h
	$(CXX)  $(FLAGS) $(INCS)  testdir/FTLMtest.cpp -o FTLMtest  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(TESTLIBS) $(ND) $(OP)

reddmtest: testdir/reddmtest.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(IXX) $(FLAGS) $(INCS) testdir/reddmtest.cpp -o reddmtest  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(TESTLIBS) $(OP)
reddm: examples/reddm.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h

	$(CXX) $(FLAGS) $(INCS) examples/reddm.cpp -o reddm  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)

holstexDi: examples/holstexDi.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp  src/diag.h

	$(CXX) $(FLAGS) $(INCS) examples/holstexDi.cpp -o holstexDi  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)

holstFTexact: examples/holstFTexact.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) examples/holstFTexact.cpp -o holstFTexact  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)
holstSpekfTex: examples/holstSpekfTex.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) examples/holstSpekfTex.cpp -o holstSpekfTex  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)
hetFerex: examples/hetFerex.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) examples/hetFerex.cpp -o hetFerex  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)
hetFerexSP: examples/hetFerexSP.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) examples/hetFerexSP.cpp -o hetFerexSP  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)

holsttimeexact: examples/holsttimeexact.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) examples/holsttimeexact.cpp -o holsttimeexact  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)

holstLTLM: examples/holstLTLM.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) examples/holstLTLM.cpp -o holstLTLM  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)
holstFTLM: examples/holstFTLM.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) examples/holstFTLM.cpp -o holstFTLM  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)

holstexDi: examples/holstexDi.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(CXX) $(FLAGS) $(INCS) examples/holstexDi.cpp -o holstexDi  $(MKLLINK32) $(LIBSLINK) $(LIBS32) $(ND) $(OP)


holstFTLMpara: examples/holstFTLMpara.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(IMPI) $(FLAGS) $(INCS) examples/holstFTLMpara.cpp -o holstFTLMpara  $(LIBSLINK) $(MPILINK32)  $(LIBS32) $(ND)  $(OP)

holstLTLMpara: examples/holstLTLMpara.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/reddm.hpp src/diag.h
	$(IMPI) $(FLAGS) $(INCS) examples/holstLTLMpara.cpp -o holstLTLMpara  $(LIBSLINK) $(MPILINK32)  $(LIBS32) $(ND)  $(OP)


exvalcompare: testdir/exvalcompare.cpp src/basis.hpp src/operators.hpp src/accesfunctions.hpp src/numerics.hpp src/timeev.hpp src/diag.h
	$(IXX) $(FLAGS) $(INCS) testdir/exvalcompare.cpp -o exvalcompare  $(MKLLINK) $(LIBSLINK) $(LIBS32) $(TESTLIBS) -DMOM $(OP)
#-static

basisbenchmark: benchmarkdir/basisbenchmark.cpp src/basis.hpp src/accesfunctions.hpp src/numerics.hpp
	$(IXX) $(FLAGS) $(INCS) benchmarkdir/basisbenchmark.cpp -o basisbenchmark  $(MKLLINK) $(LIBS)  $(ND) -g $(OP)
clean:
	rm *.o main $(EXEC)
