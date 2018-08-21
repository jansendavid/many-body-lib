// #define BOOST_TEST_DYN_LINK
// #define BOOST_TEST_MODULE Suites
// #include <boost/test/unit_test.hpp>

#include<iostream>
#include "basis.hpp"
#include "operators.hpp"
#include"numerics.hpp"
#include"diag.hpp"
#include "files.hpp"
// using namespace boost::unit_test;
// using boost::unit_test_framework::test_suite;
// using namespace Many_Body;
// BOOST_AUTO_TEST_SUITE(timeevesting)
// BOOST_AUTO_TEST_CASE(timeev)
using namespace Many_Body;
using namespace Eigen;
  template <typename Matrix>
    double diag_with_lancz(Eigen::VectorXcd& initialState, const Matrix& ham, const Matrix& obs, size_t lanczosDim)
  { double sum=0;

    Eigen::VectorXcd initialState2=initialState;
    Eigen::MatrixXcd Q(ham.rows(), lanczosDim);
    Many_Body::TriDiagMat tri=Many_Body::Lanczos(ham, initialState2, lanczosDim, Q);
    Eigen::MatrixXd S(Q.cols(), Q.cols());
              Eigen::VectorXd eigenVals(Q.cols());
    Many_Body::diag(tri, S, eigenVals);
   Eigen::MatrixXcd S2=S.cast<std::complex<double>>();  
   VectorXcd initialStateTemp(Q.cols());
    initialStateTemp= (S2.row(0).transpose());
   //  // initialStateTemp= S2*initialStateTemp;
   //  // initialState= Q*initialStateTemp;
   //  assert(std::abs(initialState.norm() -1.) < Many_Body::err);
     Eigen::MatrixXcd obs2=S.adjoint()*Q.adjoint()*obs*Q*S;
    for (size_t i = 0; i < lanczosDim; ++i)
      {
	
      std::complex<double> c=std::abs(initialStateTemp(i))*std::abs(initialStateTemp(i))*obs2(i, i);
    	 sum+=real(c);
    }
    return sum;
    
  }
int main()
{

    size_t numberOfSteps=5;
   const size_t L=8;
   ElectronBasis<L> e(4);
   size_t D=700;
   //std::cout << e << std::endl;
   Eigen::MatrixXcd AA = Eigen::MatrixXcd::Random(D, D);
   	 Eigen::MatrixXcd H = AA + AA.adjoint();
   //  Operators::Mat H= Operators::NumberOperator(e)+Operators::EKinOperator(e);
   // Eigen::VectorXcd inistate(e.D);
  
   // inistate.setZero();
   // inistate[0]=1;
	 Eigen::VectorXcd inistate=Eigen::VectorXcd::Random(D);
	 inistate=inistate/inistate.norm();
   Eigen::VectorXd eigenVals(D);
   Eigen::MatrixXcd BB = Eigen::MatrixXcd::Random(D, D);
   	 Eigen::MatrixXcd O = BB + BB.adjoint();
   // Operators::Mat O(e.dim, e.dim);
   //  O.setZero();
   //  O.coeffRef(1, 1)=0.5;
   //  O.coeffRef(2, 2)=0.5;
    
    // HAmiltonian =H
   Eigen::MatrixXcd HH=Eigen::MatrixXcd(H);
   // Eigen::MatrixXd HH2=Eigen::MatrixXd(H);

   
   
   
   Many_Body::diag(HH, eigenVals);
Eigen::VectorXcd newIn=inistate;
Eigen::VectorXcd newIn2=inistate;
  double evalL=diag_with_lancz(newIn2, H, O, 20);
 double evalN=0;
  MatrixXcd M=HH.adjoint()*O*HH;
  newIn=HH.adjoint()*newIn;
  for(size_t i=0; i<D; i++)
   {

        std::complex<double> c=M(i, i)*std::abs(newIn(i))*std::abs(newIn2(i));
        evalN+=real(c);
   }
   std::cout<< evalN << "  " << evalL << "  " << D <<std::endl;
  return 0;
}


//BOOST_AUTO_TEST_SUITE_END()
// EOF
