#include<iostream>
//#include"basis.hpp"

int main()
{
// const size_t L1=3;
//   const size_t L2=4;
//   const size_t L3=5;
 //  const size_t L4=6;


  //  many_body::PhononBasis<L1> g{1};
  
  std::cout << Many_Body::PrimeNumber(3) << '\n';
//     BosonBasis<L1> gg{1};
//   std::cout << gg << '\n';
// ElectronBasis<L1> e{1};
//      std::cout << e << '\n';
  




//    TensorProduct<  ElectronBasis<L1>,  PhononBasis<L1>> TP(e, g);

//    std::cout << TP;
  //  for(size_t i=1; i<2; i++)
  //      {
  // 	  PhononBasis<L1> g1{i};
  //         PhononBasis<L2> g2{i};
  // 	 PhononBasis<L3> g3{i};


  // 	 ElectronBasis<L1> e1;
  // 	 ElectronBasis<L2> e2{2};
  // 	 ElectronBasis<L3> e3{1};

  // 	 TensorProduct<  ElectronBasis<L1>,  PhononBasis<L1>> TP1(e1, g1);
  // 	 TensorProduct<  ElectronBasis<L2>,  PhononBasis<L2>> TP2(e2, g2);
  //        TensorProduct<  ElectronBasis<L3>,  PhononBasis<L3>> TP3(e3, g3);

  // 	 (g1.dim*e1.dim==TP1.dim);
  // 	 // BOOST_CHECK(g2.dim*e2.dim==TP2.dim);
  // 	 // BOOST_CHECK(g3.dim*e3.dim==TP3.dim);


  // }


  return 0;
}
