#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>

#include<iostream>
#include <eigen3/Eigen/Eigenvalues> 
#include"numerics.hpp"
#include"reddm.hpp"
#include"tpoperators.hpp"
#include "files.hpp"
#include"FTLanczos.hpp"
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;
using namespace Many_Body;
BOOST_AUTO_TEST_SUITE(timeevesting)
BOOST_AUTO_TEST_CASE(timeev)
{
    using Mat= Operators::Mat;
  int L=4;
  double omega=1;
  double gamma=1;
  double t0=1;
  double T=2;
   using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
        PhononBasis g2{ 2, 1};
  ElectronBasis e( L, 1);
  //  std::cout<< e<<std::endl;
  
  PhononBasis ph(L, 3);
  //std::cout<< ph<<std::endl;
  HolsteinBasis TP(e, ph);
        Mat E1=Operators::EKinOperatorL(TP, e, t0, true);
      Mat Ebdag=Operators::BosonCOperator(TP, ph, gamma, true);
      Mat Eb=Operators::BosonDOperator(TP, ph, gamma, true);
      Mat Eph=Operators::NumberOperator(TP, ph, omega,  true);
      
      //Mat E=Operators::NumberOperatore(TP, e, 1, false);
      //    std::cout<< HH << std::endl;
      Eigen::VectorXd eigenVals(TP.dim);
      Mat H=E1+Eph  +Ebdag + Eb;
      Mat O=E1;

	LTLM(H, O);
	
 }


BOOST_AUTO_TEST_SUITE_END()
