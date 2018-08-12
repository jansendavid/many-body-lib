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
  {  size_t size=5;
  Eigen::VectorXcd x = Eigen::VectorXcd::Random(size);
  x=x/x.norm();

    Eigen::MatrixXd AA = Eigen::MatrixXd::Random(size, size);
     Eigen::MatrixXd A = AA + AA.transpose();

  Eigen::MatrixXd B=A;
  Eigen::VectorXd evA(size);
  Eigen::VectorXd evB(size);
  Eigen::MatrixXcd Q(size, size);
  Many_Body::TriDiagMat tri=Many_Body::Lanczos(B, x, size, Q);

Eigen::MatrixXd S(size, size);
 S.setZero();
 diag(tri, S, evB);
  diag(A, evA);

  for (int i = 0; i < size; ++i)
  {
       BOOST_CHECK(std::abs(evA(i)-evB(i))<Many_Body::err);
  }

  }{
    size_t size=5;
  Eigen::VectorXcd x = Eigen::VectorXcd::Random(size);
  x=x/x.norm();

    Eigen::MatrixXcd AA = Eigen::MatrixXcd::Random(size, size);
    Eigen::MatrixXcd A = AA.adjoint()*AA;


  Eigen::MatrixXcd B=A;
  Eigen::VectorXd evA(size);
  Eigen::VectorXd evB(size);
  Eigen::MatrixXcd Q(size, size);
  Many_Body::TriDiagMat tri=Many_Body::Lanczos(B, x, size, Q);

Eigen::MatrixXd S(size, size);
  S.setZero();
  diag(tri, S, evB);
   diag(A, evA);
        std::cout << '\n';
  for (int i = 0; i < size; ++i)
  {
    BOOST_CHECK(std::abs(evA(i)-evB(i))<Many_Body::err);


  }
  }
  //
}


BOOST_AUTO_TEST_SUITE_END()
// EOF
