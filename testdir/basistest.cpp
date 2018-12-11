#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>
#include"basis.hpp"

using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;
using namespace Many_Body;
BOOST_AUTO_TEST_SUITE(basistesting)
BOOST_AUTO_TEST_CASE(phonondimension)
{
  // const size_t L1=3;
  // const size_t L2=4;
  // const size_t L3=5;
  // const size_t L4=6;
  BosonState b1(5);
   std::vector<int> state {1, 4, 1, 0};
     BosonState b2(state);
     std::cout<< b2[1]<< std::endl;
     BOOST_CHECK(std::abs(b1.GetId()-0)<0.0001);
          BOOST_CHECK(std::abs(b2[1]-4)<0.0001);
	  PhononBasis g1{1, 4};

   for(size_t i=2; i<3; i++)
     {
       PhononBasis g1{i, 4};
       PhononBasis g2{i, 3};
  //     PhononBasis<L3> g3{i};
  //     PhononBasis<L4> g4{i};

       BOOST_CHECK(g1.dim==static_cast<size_t>(std::pow((i+1), 4)));
       BOOST_CHECK(g2.dim==static_cast<size_t>(std::pow((i+1), 3)));
  //     BOOST_CHECK(g3.dim==static_cast<size_t>(std::pow((i+1), L3)));
  //     BOOST_CHECK(g4.dim==static_cast<size_t>(std::pow((i+1), L4)));
      
  }
}
BOOST_AUTO_TEST_CASE(electrondimension)
{
  ElectronState df;
  std::vector<int> state {1, 1, 1, 0};
  ElectronState constru1(state);
    BOOST_CHECK(std::abs(constru1.GetId()-14.)<0.0001);
    constru1.flip(1);
      BOOST_CHECK(std::abs(constru1.GetId()-10.)<0.0001);
      ElectronState constru2(3, 4);
      BOOST_CHECK(std::abs(constru2.GetId()-3.)<0.0001);
ElectronBasis g1( 3);
      const ElectronBasis g2(2,5);
      ElectronBasis g3(1, 3);
      
      BOOST_CHECK(g1.dim==std::pow(2, 3));
    BOOST_CHECK(g2.dim==Factorial(5)/(Factorial(5-2)*Factorial(2)));
    BOOST_CHECK(g3.dim== 3);
    



 }

BOOST_AUTO_TEST_CASE(bosondimension)
{

   // const size_t L2=4;
   // const size_t L3=5;
   // std::array<size_t, L2> stateArray{0, 1, 2};
   
   //   BosonState<L2> a{stateArray};

   // const BosonBasis<L2> g2{2};
   // BosonBasis<L3> g3{1};


      
   // //BOOST_CHECK(g1.dim==Factorial(L2)/(Factorial(L2-2)*Factorial(2)));
   // BOOST_CHECK(g2.dim==Factorial(L2+2-1)/(Factorial(L2-1)*Factorial(2)));
   // BOOST_CHECK(g3.dim==Factorial(L3+1-1)/(Factorial(L3-1)*Factorial(1)));



 }

BOOST_AUTO_TEST_CASE(tensorproductdimension)
{
  // const size_t L1=3;
 //  const size_t L2=4;
 //  const size_t L3=5;
 // //  const size_t L4=6;


  PhononBasis g{1, 4};
  ElectronBasis e{1, 4};
 TensorProduct<  ElectronBasis,  PhononBasis> TP(e, g);


    for(size_t i=1; i<2; i++)
        {
	  PhononBasis g1{i, 3};
	  PhononBasis g2{i, 4};
	  PhononBasis g3{i, 5};


	  ElectronBasis e1{1, 3};
 	 ElectronBasis e2{2, 4};
 	 ElectronBasis e3{3, 5};

 TensorProduct<  ElectronBasis,  PhononBasis> TP1(e1, g1);
  	 TensorProduct<  ElectronBasis,  PhononBasis> TP2(e2, g2);
          TensorProduct<  ElectronBasis,  PhononBasis> TP3(e3, g3);

  	 BOOST_CHECK(g1.dim*e1.dim==TP1.dim);
 	 BOOST_CHECK(g2.dim*e2.dim==TP2.dim);
  	 BOOST_CHECK(g3.dim*e3.dim==TP3.dim);


  }
}
BOOST_AUTO_TEST_SUITE_END()
