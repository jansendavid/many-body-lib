#include<iostream>
#include"basis.hpp"
#include"operators.hpp"
#include"diag.h"
#include"tpoperators.hpp"
#include"files.hpp"
#include"timeev.hpp"
#include"ETH.hpp"
int main(int argc, char *argv[])
{
  using namespace Eigen;
using namespace std;
using namespace Many_Body;
  using Fermi_HubbardBasis= TensorProduct<ElectronBasis, ElectronBasis>;
   using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
  using Mat= Operators::Mat;
 //  const size_t L=10;

   
 //  const double t1=std::stod(argv[1]);
 //  // const double t11=std::stoi(argv[2]);
 //  const double t2= std::stod(argv[2]);
 //  // const double t22= std::stod(argv[4]);
 //  const double u= std::stod(argv[3]);
 //  size_t electrons=int(L/2);
 // std::vector<int> ee(L, 0);
 //      ee[L-1]=1;
 //      ElectronState estate1(ee);
 //      ElectronBasis e(0, L);


 // ElectronBasis e2(electrons, L);
 // Fermi_HubbardBasis TP(e, e2);
 // // std::cout<< e<< std::endl;
 // // std::cout<< e2 << std::endl;
 
 //   Mat E1=Operators::EKinOperatorL(TP, e, t1);
 //     Mat E2=Operators::EKinOperatorR(TP, e2, t2);
 // // //   //   Mat E11=Operators::EKinOperatorLNNN(TP, e1, t11);
 // // //   // Mat E22=Operators::EKinOperatorRNNN(TP, e2, t22);
 //     Mat C=Operators::CalculateCouplungOperator(TP, e2, u);
 //    Mat H=E1+E2+C;
   const size_t L=4;

  
  const double t1=std::stod(argv[1]);
  const double omg=std::stod(argv[2]);
  const double gamma=std::stod(argv[3]);
  //std::vector<int> ee(L, 0);
  // ee[L-1]=1;
  double dt=0.0001;
  double ttot=0.1;
  ElectronState estate1(L, 1);
  ElectronBasis e( L, 1);
  size_t M=3;
  //std::cout<< e<<std::endl;
std::vector<size_t> es(L, 0);
      es[0]=1;
      ElectronState aa(es);
  PhononBasis ph(L, M);
  //               std::cout<< ph<<std::endl;
      HolsteinBasis TP(e, ph);
      std::vector<size_t> state(L, 0);
  std::fill(state.begin(), state.end(), 1);
  BosonState b2(state, M);
 
     auto it3=TP.lbasis.find(aa.GetId());
     auto it4=TP.rbasis.find(b2.GetId());
     size_t StateNr= Position(*it4)*TP.lbasis.dim +Position(*it3);
     //        std::cout<< TP << std::endl;
     //    std::cout<< b2.GetId()<<std::endl;
     // 	std::cout<< b2.GetId()<<std::endl;
     // std::cout<<"state nr "<< StateNr<<std::endl;               
         Eigen::VectorXcd inistate(TP.dim);

 	 Mat O=Operators::NumberOperator(TP, ph, 1,  false);
   inistate.setZero();

   inistate[StateNr]=1;
   	          Eigen::VectorXcd i0=inistate;
      //   std::cout<< e << std::endl;
		  
           std::cout<< TP.dim << std::endl;     
      Mat E1=Operators::EKinOperatorL(TP, e, t1, false);
      Mat Ebdag=Operators::BosonCOperator(TP, ph, gamma, false);
      Mat Eb=Operators::BosonDOperator(TP, ph, gamma, false);
      Mat Eph=Operators::NumberOperator(TP, ph, omg,  false);
      
      //Mat E=Operators::NumberOperatore(TP, e, 1, false);
      //    std::cout<< HH << std::xbendl;
      Eigen::VectorXd eigenVals(TP.dim);
      Mat H=E1+Eph  +Ebdag + Eb;
      //    
 
    Eigen::MatrixXd HH=Eigen::MatrixXd(H);
    //        std::cout<< HH<<std::endl;
        Eigen::MatrixXd N=Eigen::MatrixXd(Eph);
     Many_Body::diagMat(HH, eigenVals);
         Eigen::VectorXd energy(TP.dim);
    	   std::vector<double> obs(TP.dim);
	   std::vector<double> obs2(TP.dim);
	   std::vector<double> ensvec(TP.dim);
	   Eigen::MatrixXd PHD2=HH.adjoint()*N.selfadjointView<Lower>()*HH;

     Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals, dt);
   Eigen::MatrixXcd cEVec=HH.cast<std::complex<double>>();
    Eigen::VectorXcd newIn=inistate;
   int i=0;
   // O=H;
        while(i*dt<ttot)
       {

       	 TimeEv::timeev_exact(newIn, cEVec, evExp);
  	 
   	  	 std::complex<double> c=(newIn.adjoint()*(O*newIn))(0);
		 //	 std::complex<double> c=(inistate.adjoint()*(newIn))(0);
  	 
   	 i++;

   			// 	outputVals(i)=real(c);
   			// 	outputVals2(i)=real(c2);
        		// outputTime(i)=i*dt;
   	 std::cout<< real(c)<< " dt "<< i*dt<< std::endl;
	 //	  	 BOOST_CHECK(std::abs(real(c2)-real(c))<Many_Body::err);
	     	             }

 //std::cout<< MatrixXd(OBS2) << std::endl;
 // Eigen::MatrixXd OBS22=HH.adjoint()*(OBS2.selfadjointView<Lower>())*HH;
 // Eigen::VectorXd diagobs2=OBS22.diagonal();
	// std::cout<< en(0) << std::endl;
 // std::cout<<diagobs2(0)<<std::endl; 
    std::string filename=".bin";
    //     bin_write("E"+filename, en);
  return 0;
}
 
