#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>
#define EIGEN_USE_MKL_ALL
#include<iostream>
#include <Eigen/Eigenvalues> 
#include"numerics.hpp"
#include"reddm.hpp"
#include"tpoperators.hpp"
#include "files.hpp"
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
  std::cout<< TP.dim << std::endl;
        Mat E1=Operators::EKinOperatorL(TP, e, t0, true);
      Mat Ebdag=Operators::BosonCOperator(TP, ph, gamma, true);
      Mat Eb=Operators::BosonDOperator(TP, ph, gamma, true);
      Mat Eph=Operators::NumberOperator(TP, ph, omega,  true);
      
      //Mat E=Operators::NumberOperatore(TP, e, 1, false);
      //    std::cout<< HH << std::endl;
      Eigen::VectorXd eigenVals(TP.dim);
            Mat H=E1+Eph  +Ebdag + Eb;
          Eigen::MatrixXd HH=Eigen::MatrixXd(H);
	  auto optModes=makeThermalRDMTP(HH, TP, T,2);
	  double sum=0;
	  for(auto& l : optModes)
	    { 
	      std::cout<< l<<std::endl;
	      sum+=l.sum();
	    }
	  std::cout << "and the sum of all was "<< sum<<std::endl;
  //   std::cout<< TP<< std::endl;
  // Eigen::VectorXd v1=Eigen::VectorXd::Zero(TP.dim);
  // v1[1]=1/std::sqrt(2);
  // Eigen::VectorXd v2=Eigen::VectorXd::Zero(TP.dim);
  //   v2[5]=1/std::sqrt(2);
  // 	Eigen::VectorXd V1=v1 +v2;
  // 	Eigen::MatrixXd V2=V1.transpose();
  // 	// std::cout<< V1<<std::endl;
  // 	// std::cout<< V2<<std::endl;	
  // 	MatrixXd M=V1*(V2);
  // 	// std::cout<< M<<std::endl;
  // 	// 	std::cout<< g2<<std::endl;
  // 	//	makeRedDM(g2, 0, M);
  // 	makeRedDMTP(TP, ph,  0, M);
	
 }


BOOST_AUTO_TEST_SUITE_END()
