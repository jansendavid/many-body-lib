#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>
#include"basis.hpp"
#include"operators.hpp"
#include"timeev.hpp"
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;
using namespace Many_Body;
BOOST_AUTO_TEST_SUITE(timeevesting)
BOOST_AUTO_TEST_CASE(timeev)
{
  size_t numberOfSteps=100;
   const size_t L=10;
   ElectronBasis<L> e(5);
   using NO=Operators::NumberOperator<   ElectronBasis<L>>;
   using EKO= Operators::EKinOperator<   ElectronBasis<L>>;
   Operators::Mat H( EKO(e).mat+NO(e).mat);
  Eigen::VectorXd inistate;
  Eigen::VectorXd eigenVals;
  
  Operators::mat Observable = NO(e).mat;
    // HAmiltonian =H
   Many_Body::diagsym(H, eigenVals);
 
  double dt= 0.01;
  
  Eigen::VectorXd output(numberOfSteps);
  TimeEv::timeev_exact(inistate, H, O, output, dt, numberOfSteps);
}


BOOST_AUTO_TEST_SUITE_END()
// EOF


