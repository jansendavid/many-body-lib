#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>
#include"basis.hpp"
#include"operators.hpp"
#include"timeev.hpp"
#include"diag.hpp"
#include "files.hpp"
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;
using namespace Many_Body;
BOOST_AUTO_TEST_SUITE(timeevesting)
BOOST_AUTO_TEST_CASE(timeev)
{
  size_t numberOfSteps=400;
   const size_t L=3;
   ElectronBasis<L> e(2);
   std::cout << e << std::endl;
   Operators::Mat H= Operators::NumberOperator(e)+Operators::EKinOperator(e);
   Eigen::VectorXd inistate(e.dim);
   inistate[0]=1;
   Eigen::VectorXd eigenVals(e.dim);
    std::cout <<H<< std::endl;
    //Operators::Mat O = Operators::EKinOperator(e);
    Eigen::MatrixXd O(e.dim, e.dim);
    O.setZero();
    O(1, 1)=0.5;
    O(2, 2)=0.5;
    
    // HAmiltonian =H
   Eigen::MatrixXd HH=Eigen::MatrixXd(H);
   Eigen::MatrixXd OO=Eigen::MatrixXd(O);
   Many_Body::diag(HH, eigenVals);

   double dt= 0.1;
   //std::cout <<HH.adjoint()*HH<< std::endl;
   std::cout <<eigenVals<< std::endl;
   Eigen::VectorXd outputTime(numberOfSteps);
     Eigen::VectorXd outputVals(numberOfSteps);
     // TimeEv::timeev_exact(inistate, HH, OO, eigenVals, outputVals, outputTime, dt, numberOfSteps);
     // Many_Body::ToFile(outputTime, outputVals, "timetestexact.dat", numberOfSteps);
     Eigen::VectorXd q(numberOfSteps);
     Many_Body::Lanczos(HH, inistate, 12);
}


BOOST_AUTO_TEST_SUITE_END()
// EOF


