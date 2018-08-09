#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>

#include<iostream>

#include"numerics.hpp"
#include"diag.hpp"
#include "files.hpp"
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;
using namespace Many_Body;
BOOST_AUTO_TEST_SUITE(timeevesting)
BOOST_AUTO_TEST_CASE(timeev)
{
  size_t size=10;
  Eigen::VectorXd x = Eigen::VectorXd::Random(size);
  x=x/x.norm();

    Eigen::MatrixXd AA = Eigen::MatrixXd::Random(size, size);
     Eigen::MatrixXd A = AA + AA.transpose();
  std::cout<< A<<std::endl;
  Eigen::MatrixXd B=A;
  Eigen::VectorXd evA(size);
  Eigen::VectorXd evB(size);
  Many_Body::TriDiagMat tri=Many_Body::Lanczos(B, x, size);
  // std::cout<< "en"<<std::endl;
  //   std::cout<< tri.diagel<<std::endl;
  //     std::cout<< "to"<<std::endl;
  //       std::cout<< tri.offDiag<<std::endl;
  // 	 std::cout<< "to"<<std::endl;
Eigen::MatrixXd S(size, size);
 S.setZero();
 diag(tri, S, evB);
  diag(A, evA);
  //Many_Body::diag(A, evA);
  for (int i = 0; i < size; ++i)
  {
    std::cout << std::abs(evA(i)-evB(i)) << '\n';
  }
  
  //
}


BOOST_AUTO_TEST_SUITE_END()
// EOF
