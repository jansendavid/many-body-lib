#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#define EIGEN_USE_MKL_ALL 
#include <boost/test/unit_test.hpp>
#include"basis.hpp"
#include"operators.hpp"
#include"tpoperators.hpp"

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

   ElectronBasis<L> e1(1);
   //std::cout << int(L/2) << std::endl;
   ElectronBasis<L> e2(2);
      double t1=1;
   double t2=1;
   double u=1;
   // std::cout<< e1;
   //   std::cout<< e2;
   TensorProduct<ElectronBasis<L>, ElectronBasis<L>> TP(e1, e2);
   //        std::cout<< TP;
   //Mat E11=Operators::EKinOperator(TP,  t1);
   //   std::cout<< E11;
   Mat E1=Operators::EKinOperatorL(TP, e1, t1);
   Mat E2=Operators::EKinOperatorR(TP, e2, t2);
    Mat C=Operators::CalculateCouplungOperator(TP, e2, u);
    Mat H=E1+E2+ C;

    Eigen::MatrixXd HH =Eigen::MatrixXd(H);
    std::cout << HH << std::endl;
    Eigen::VectorXd ev(TP.dim);
    // diag(HH, ev);

    // std::cout << ev << std::endl;
   // std::cout << NR;
   // Mat EK=Operators::EKinOperator(e);
   
   //    std::cout << EK;
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


