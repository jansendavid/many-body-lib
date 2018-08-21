#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>
#include"basis.hpp"
#include"operators.hpp"
//#include"diag.hpp"

using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;
using namespace Many_Body;
BOOST_AUTO_TEST_SUITE(basistesting)

BOOST_AUTO_TEST_SUITE(operatortesting)
BOOST_AUTO_TEST_CASE(operatoroperations)
{
  using namespace Eigen;
  using Mat=Operators::Mat;
   const size_t L=3;
   ElectronBasis<L> e(2);
   Mat NR=Operators::NumberOperator(e);
   std::cout << NR;
   Mat EK=Operators::EKinOperator(e);
   
      std::cout << EK;
   // VectorXd v(5);
  //  MatrixXcd mat = MatrixXcd::Random(5, 5);
  // MatrixXcd mat4= mat+mat.adjoint();
  //  MatrixXcd mat2=mat4;

  //  Many_Body::diagherm(mat4, v);
  //     MatrixXcd mat3=mat4.adjoint();
  //  std::cout<< mat3*mat2*mat4<< std::endl;
  //   std::cout<< v<< std::endl;
  //   numberoperator<ElectronBasis> n2;
  //      kineticoperator<ElectronBasis> cdagc1(e);
  // 	kineticoperator<ElectronBasis> cdagc2;
  
  // 	double t=3;
  // 	souble l=1;
  // 	hamiltonian<ElectronBasis> H1{t*n1+ l*cdagc1};
  // 		hamiltonian<ElectronBasis > H2{t*n2+ l*cdagc2};
  // 		H2(e);
  // 		cout<< H2.energies();
 // 				cout<< H1.energies(); 

}      
}
// BOOST_AUTO_TEST_CASE(electrondimension)
// {
// }
BOOST_AUTO_TEST_SUITE_END()
// EOF


