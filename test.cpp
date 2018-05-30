#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>
#include"basis.hpp"
#include"operators.hpp"
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;

BOOST_AUTO_TEST_SUITE(basistesting)
BOOST_AUTO_TEST_CASE(phonondimension)
{
  const size_t L1=3;
  const size_t L2=4;
  const size_t L3=5;
  const size_t L4=6;

  for(size_t i=2; i<3; i++)
    {
      PhononBasis<L1> g1{i};
      PhononBasis<L2> g2{i};
      PhononBasis<L3> g3{i};
      PhononBasis<L4> g4{i};

      BOOST_CHECK(g1.dim==static_cast<size_t>(std::pow((i+1), L1)));
      BOOST_CHECK(g2.dim==static_cast<size_t>(std::pow((i+1), L2)));
      BOOST_CHECK(g3.dim==static_cast<size_t>(std::pow((i+1), L3)));
      BOOST_CHECK(g4.dim==static_cast<size_t>(std::pow((i+1), L4)));
      
  }
}
BOOST_AUTO_TEST_CASE(electrondimension)
{
   const size_t L1=3;
   const size_t L2=4;
   const size_t L3=5;

   ElectronBasis<L1> g1;
   const ElectronBasis<L2> g2{2};
   ElectronBasis<L3> g3{1};


      
   BOOST_CHECK(g1.dim==static_cast<size_t>(std::pow((2), L1)));
   BOOST_CHECK(g2.dim==Factorial(L2)/(Factorial(L2-2)*Factorial(2)));
   BOOST_CHECK(g3.dim== L3);



 }

BOOST_AUTO_TEST_CASE(bosondimension)
{

   const size_t L2=4;
   const size_t L3=5;


   const BosonBasis<L2> g2{2};
   BosonBasis<L3> g3{1};


      
   //BOOST_CHECK(g1.dim==Factorial(L2)/(Factorial(L2-2)*Factorial(2)));
   BOOST_CHECK(g2.dim==Factorial(L2+2-1)/(Factorial(L2-1)*Factorial(2)));
   BOOST_CHECK(g3.dim==Factorial(L3+1-1)/(Factorial(L3-1)*Factorial(1)));



 }

BOOST_AUTO_TEST_CASE(tensorproductdimension)
{
  const size_t L1=3;
  const size_t L2=4;
  const size_t L3=5;
 //  const size_t L4=6;


  PhononBasis<L1> g{1};
  ElectronBasis<L1> e{1};
  TensorProduct<  ElectronBasis<L1>,  PhononBasis<L1>> TP(e, g);

  
   for(size_t i=1; i<2; i++)
       {
	 PhononBasis<L1> g1{i};
         PhononBasis<L2> g2{i};
	 PhononBasis<L3> g3{i};


	 ElectronBasis<L1> e1;
	 ElectronBasis<L2> e2{2};
	 ElectronBasis<L3> e3{1};

	 TensorProduct<  ElectronBasis<L1>,  PhononBasis<L1>> TP1(e1, g1);
	 TensorProduct<  ElectronBasis<L2>,  PhononBasis<L2>> TP2(e2, g2);
         TensorProduct<  ElectronBasis<L3>,  PhononBasis<L3>> TP3(e3, g3);

	 BOOST_CHECK(g1.dim*e1.dim==TP1.dim);
	 BOOST_CHECK(g2.dim*e2.dim==TP2.dim);
	 BOOST_CHECK(g3.dim*e3.dim==TP3.dim);


  }
}
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE(operatortesting)
BOOST_AUTO_TEST_CASE(operatoroperations)
{
   // const size_t L=3;
   // ElectronBasis<L> e;
   // operators::numberoperator<ElectronBasis<L>> n1(e);
  //   numberoperator<ElectronBasis> n2;
  //       kineticoperator<ElectronBasis> cdagc1(e);
  // 	kineticoperator<ElectronBasis> cdagc2;
  
  // 	double t=3;
  // 	souble l=1;
  // 	hamiltonian<ElectronBasis> H1{t*n1+ l*cdagc1};
  // 		hamiltonian<ElectronBasis > H2{t*n2+ l*cdagc2};
  // 		H2(e);
  // 		cout<< H2.energies();
  // 				cout<< H1.energies(); 

       
}
// BOOST_AUTO_TEST_CASE(electrondimension)
// {
// }
BOOST_AUTO_TEST_SUITE_END()
// EOF


