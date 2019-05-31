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
  const size_t L1=3;
  const size_t L2=4;
  const size_t L3=5;
  const size_t L4=6;

  BosonState b1( 4, 5);
    
  std::vector<int> state {1, 4, 1, 0};
    BosonState b2(state);
    
    
 
      BOOST_CHECK(b1.GetId()==0);
      BOOST_CHECK(b2[1]==4);
      PhononBasis g1{ 4, 1};

    for(size_t i=2; i<3; i++)
      {
	PhononBasis g1{ 4, i};
        PhononBasis g2{ 3, i};
   

        BOOST_CHECK(g1.dim==static_cast<size_t>(std::pow((i+1), 4)));
        BOOST_CHECK(g2.dim==static_cast<size_t>(std::pow((i+1), 3)));
   
      
   }
}
BOOST_AUTO_TEST_CASE(electrondimension)
{
  ElectronState df;
  std::vector<int> state {1, 1,      1, 0};
  ElectronState constru1(state);

     BOOST_CHECK(constru1.GetId()==14);
      constru1.flip(1);

      BOOST_CHECK(constru1.GetId()==10);
        ElectronState constru2(3, 4);
	
//        BOOST_CHECK(constru2.GetId()==3);
 ElectronBasis g1( 3);

 const ElectronBasis g2(5, 2);
       ElectronBasis g3( 3, 1);
      
       BOOST_CHECK(g1.dim==std::pow(2, 3));
     BOOST_CHECK(g2.dim==Factorial(5)/(Factorial(5-2)*Factorial(2)));
     BOOST_CHECK(g3.dim== 3);
    



 }

// BOOST_AUTO_TEST_CASE(bosondimension)
// {

   // const size_t L2=4;
   // const size_t L3=5;
   // std::array<size_t, L2> stateArray{0, 1, 2};
   
   //   BosonState<L2> a{stateArray};

   // const BosonBasis<L2> g2{2};
   // BosonBasis<L3> g3{1};


      
   // //BOOST_CHECK(g1.dim==Factorial(L2)/(Factorial(L2-2)*Factorial(2)));
   // BOOST_CHECK(g2.dim==Factorial(L2+2-1)/(Factorial(L2-1)*Factorial(2)));
   // BOOST_CHECK(g3.dim==Factorial(L3+1-1)/(Factorial(L3-1)*Factorial(1)));



// }

 BOOST_AUTO_TEST_CASE(tensorproductdimension)
 {
//   // const size_t L1=3;
//  //  const size_t L2=4;
//  //  const size_t L3=5;
//  // //  const size_t L4=6;


   PhononBasis g{ 4, 1};

   

 ElectronBasis e{ 4, 1};

 TensorProduct<  ElectronBasis,  PhononBasis> TP(e, g);


     for(size_t i=1; i<2; i++)
         {
	   PhononBasis g1{3, i};
	   PhononBasis g2{ 4, i};
	   PhononBasis g3{ 5, i};


	   ElectronBasis e1{ 3, 1};
	   ElectronBasis e2{ 4, 2};
	   ElectronBasis e3{ 5, 3};

 TensorProduct<  ElectronBasis,  PhononBasis> TP1(e1, g1);
  	 TensorProduct<  ElectronBasis,  PhononBasis> TP2(e2, g2);
          TensorProduct<  ElectronBasis,  PhononBasis> TP3(e3, g3);

  	 BOOST_CHECK(g1.dim*e1.dim==TP1.dim);
 	 BOOST_CHECK(g2.dim*e2.dim==TP2.dim);
  	 BOOST_CHECK(g3.dim*e3.dim==TP3.dim);


   }
}
BOOST_AUTO_TEST_SUITE_END()
