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
  {
  size_t numberOfSteps=10;
   const size_t L=4;
   ElectronBasis<L> e(2);
   //std::cout << e << std::endl;
   // Eigen::MatrixXd AA = Eigen::MatrixXd::Random(e.dim, e.dim);
   // 	 Eigen::MatrixXd H = AA + AA.transpose();
        Operators::Mat H= Operators::NumberOperator(e)+Operators::EKinOperator(e);
   Eigen::VectorXcd inistate(e.dim);
  
   inistate.setZero();
   inistate[0]=1;
   Eigen::VectorXd eigenVals(e.dim);

   Operators::Mat O(e.dim, e.dim);
    O.setZero();
    O.coeffRef(1, 1)=0.5;
    O.coeffRef(2, 2)=0.5;
    
    // HAmiltonian =H
   Eigen::MatrixXd HH=Eigen::MatrixXd(H);
   // Eigen::MatrixXd HH2=Eigen::MatrixXd(H);

   
   
   
   Many_Body::diag(HH, eigenVals);

   double dt= 5;

	

   Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals, dt);
   Eigen::MatrixXcd cEVec=HH.cast<std::complex<double>>();
   Eigen::VectorXd outputTime(numberOfSteps);
   Eigen::VectorXd outputVals(numberOfSteps);
       Eigen::VectorXd outputVals2(numberOfSteps);
       Eigen::VectorXcd newIn=inistate;
              Eigen::VectorXcd newIn2=inistate;
           for (size_t i = 0; i < numberOfSteps; ++i)
       {

       	 TimeEv::timeev_exact(newIn, cEVec, evExp);
  	 TimeEv::timeev_lanzcos(newIn2, H, 3, dt); 
  	 std::complex<double> c=(newIn.adjoint()*(O*newIn))(0);
  	 std::complex<double> c2=(newIn2.adjoint()*(O*newIn2))(0);

  			// 	outputVals(i)=real(c);
  			// 	outputVals2(i)=real(c2);
        		// outputTime(i)=i*dt;
	 
  	 BOOST_CHECK(std::abs(real(c2)-real(c))<Many_Body::err);
        }
     
      // Many_Body::ToFile(outputTime, outputVals, "timetestexact.dat", numberOfSteps);
      // Many_Body::ToFile(outputTime, outputVals2, "timetestexact2.dat", numberOfSteps);
     // Eigen::VectorXd q(numberOfSteps);
       }

	   {
  size_t numberOfSteps=2;
   const size_t L=4;
   ElectronBasis<L> e(2);
   // Operators::Mat H= Operators::NumberOperator(e)+Operators::EKinOperator(e);
   //std::cout << e << std::endl;
   //   Eigen::MatrixXd HH=Eigen::MatrixXd(H);
    Eigen::MatrixXcd AA = Eigen::MatrixXcd::Random(e.dim, e.dim);
   Eigen::MatrixXcd H = AA.adjoint()*AA;
Eigen::MatrixXcd HH=Eigen::MatrixXcd(H);
   Eigen::VectorXcd inistate=Eigen::VectorXcd::Random(e.dim);
   inistate=inistate/inistate.norm();
  
   Eigen::VectorXd eigenVals(e.dim);
   // Eigen::MatrixXcd BB = Eigen::MatrixXcd::Random(e.dim, e.dim);
   // Eigen::MatrixXcd B = BB.adjoint()*BB;
  Operators::Mat B(e.dim, e.dim);
    B.setZero();
    B.coeffRef(1, 1)=0.5;
    B.coeffRef(2, 2)=0.5;
    

    //Eigen::MatrixXcd HH=Eigen::MatrixXcd(H);


   
   
   
   Many_Body::diag(HH, eigenVals);

   double dt= 5;

	   Eigen::MatrixXcd cEVec=HH.cast<std::complex<double>>();

   Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals, dt);

   Eigen::VectorXd outputTime(numberOfSteps);
   Eigen::VectorXd outputVals(numberOfSteps);
   Eigen::VectorXd outputVals2(numberOfSteps);
   Eigen::VectorXcd newIn=inistate;
   Eigen::VectorXcd newIn2=inistate;
           for (size_t i = 0; i < numberOfSteps; ++i)
       {

       	 TimeEv::timeev_exact(newIn, cEVec, evExp);
	 TimeEv::timeev_lanzcos(newIn2, H, e.dim, dt); 
	 std::complex<double> c=(newIn.adjoint()*(B*newIn))(0);
	 std::complex<double> c2=(newIn2.adjoint()*(B*newIn2))(0);

	 	 BOOST_CHECK(std::abs(real(c2)-real(c))<0.001);
        }

       }
}


BOOST_AUTO_TEST_SUITE_END()
// EOF


