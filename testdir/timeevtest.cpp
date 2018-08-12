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
  size_t numberOfSteps=10;
   const size_t L=3;
   ElectronBasis<L> e(2);
   //std::cout << e << std::endl;
   Operators::Mat H= Operators::NumberOperator(e)+Operators::EKinOperator(e);
   Eigen::VectorXd inistate=Eigen::VectorXd::Random(e.dim);
   inistate=inistate/inistate.norm();
   // inistate.setZero();
   // inistate[0]=1;
   Eigen::VectorXd eigenVals(e.dim);
   //std::cout <<H<< std::endl;
    //Operators::Mat O = Operators::EKinOperator(e);
    Eigen::MatrixXd O(e.dim, e.dim);
    O.setZero();
    O(1, 1)=0.5;
    O(2, 2)=0.5;
    
    // HAmiltonian =H
   // Eigen::MatrixXd HH=Eigen::MatrixXd(H);
   // Eigen::MatrixXd HH2=Eigen::MatrixXd(H);
    Eigen::MatrixXcd AA = Eigen::MatrixXcd::Random(e.dim, e.dim);
    Eigen::MatrixXcd HH = AA.adjoint()*AA;
        Eigen::MatrixXcd HH2=HH;
   Eigen::MatrixXd OO=O;
   // .cast<std::complex<double>>();
   Many_Body::diag(HH, eigenVals);

   double dt= 5;
   Eigen::VectorXcd newIn=inistate.cast<std::complex<double>>();
      Eigen::VectorXcd newIn2=inistate.cast<std::complex<double>>();
   std::cout <<e.dim<< std::endl;
   Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals, dt);
   Eigen::MatrixXcd cEVec=HH.cast<std::complex<double>>();
   Eigen::VectorXd outputTime(numberOfSteps);
   Eigen::VectorXd outputVals(numberOfSteps);
       Eigen::VectorXd outputVals2(numberOfSteps);
           for (size_t i = 0; i < numberOfSteps; ++i)
       {
	    // Eigen::VectorXd inistate2=Eigen::VectorXd::Random(e.dim);
	    //     Eigen::VectorXcd newIn2=inistate2.cast<std::complex<double>>();
	 // Eigen::MatrixXd HH2=AA + AA.transpose();
	// 		newIn2=newIn2/newIn2.norm();
       	 TimeEv::timeev_exact(newIn, cEVec, evExp);
       	TimeEv::timeev_lanzcos(newIn2, HH2, 3, dt); 
     		std::complex<double> c=(newIn.adjoint()*(OO*newIn))(0);
       		std::complex<double> c2=(newIn2.adjoint()*(OO*newIn2))(0);
       		// assert(imag(c)<Many_Body::err);
   //    // 				outputVals(i)=real(c);
      							outputVals2(i)=real(c2);
        		outputTime(i)=i*dt;
			std::cout << real(c) << "  " << real(c2) << '\n';
        }
     
      // Many_Body::ToFile(outputTime, outputVals, "timetestexact.dat", numberOfSteps);
      // Many_Body::ToFile(outputTime, outputVals2, "timetestexact2.dat", numberOfSteps);
     // Eigen::VectorXd q(numberOfSteps);
     
}


BOOST_AUTO_TEST_SUITE_END()
// EOF


