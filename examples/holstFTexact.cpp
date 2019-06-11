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

  const size_t L=4;

  
  const double t1=std::stod(argv[1]);
  const double omg=std::stod(argv[2]);
  const double gamma=std::stod(argv[3]);
    std::string somg=std::string(argv[3]).substr(0,3);
  std::string sgam=std::string(argv[2]).substr(0,3);
  int M=3;
  
  double mean=1;
  bool PB=false;
  //*M*omg*0.5; 
  //double beta=0.1;
    ///mean;
  //    	      double T=1./beta;

  //std::vector<int> ee(L, 0);
  // ee[L-1]=1;
  double dt=0.1;
  ElectronState estate1(L, 1);
  ElectronBasis e( L, 1);
  std::cout<< e<<std::endl;
  
  PhononBasis ph(L, M);
  //                 std::cout<< ph<<std::endl;
      HolsteinBasis TP(e, ph);
      //      std::cout<< TP << std::endl;             
         Eigen::VectorXcd inistate(TP.dim);

 	 Mat O=Operators::NumberOperator(TP, ph, 1,  PB);
   inistate.setZero();

   inistate[0]=1;
   	          Eigen::VectorXcd i0=inistate;
      //   std::cout<< e << std::endl;
      //std::cout<< ph << std::endl;    
		             std::cout<< TP.dim << std::endl;     
      Mat E1=Operators::EKinOperatorL(TP, e, t1,PB);
       Mat Ebdag=Operators::BosonCOperator(TP, ph, gamma, PB);
       Mat Eb=Operators::BosonDOperator(TP, ph, gamma, PB);
      Mat Eph=Operators::NumberOperator(TP, ph, omg,  PB);
      
      //Mat E=Operators::NumberOperatore(TP, e, 1, false);
      //    std::cout<< HH << std::endl;
      Eigen::VectorXd eigenVals(TP.dim);
      Mat H=E1+Eph+Ebdag + Eb;
      

    Eigen::MatrixXd HH=Eigen::MatrixXd(H);
    //         std::cout<< HH<<std::endl;
        Eigen::MatrixXd N=Eigen::MatrixXd(Eph);
     Many_Body::diagMat(HH, eigenVals);
     //         	      std::cout<< eigenVals<< std::endl;
     //std::cout<< std::endl;
	      //	      std::cout<< " mean "<< eigenVals.mean() << "  ana min " << 0.5*L*M<< std::endl;
         Eigen::VectorXd energy(TP.dim);
    	   std::vector<double> obs(TP.dim);
    	   	   std::vector<double> obs2(TP.dim);
    	   	   std::vector<double> ensvec(TP.dim);
    		       	      Eigen::MatrixXd PHD2=HH.adjoint()*N.selfadjointView<Lower>()*HH;
    	      Eigen::VectorXd v=PHD2.diagonal(); 
    	   for(int i=0; i<energy.size(); i++)
    {
      //obs[i]=real(obstot[i]);
      
      ensvec[i]=eigenVals(i);
      obs[i]=v(i);

      //std::cout<< ensvec[i]<<'\n';
    }
	   	    std::cout<< "gS "<<   eigenVals[0] <<std::endl;
		    std::cout<< "mean "<<   eigenVals.mean() <<std::endl;
	   std::vector<double> Ovec;

	   	   std::vector<double> Evec;
		   std::vector<double> Tr={0.01, 0.05, 0.1, 0.15, 0.2, 0.25};
		  	  for(auto t: Tr){
	     
		 // //	      std::cout<< l<<std::endl;
	         // bin_write(filename, l);
		 // n++;
		  

		      Ovec.push_back(expvalCan(ensvec, obs,  t));
		      	 Evec.push_back(expvalCan(ensvec, ensvec,  t)); 

	    }
 
			  int pb=PB;
			  std::string Hs="H";
			  	  std::string Ts="T";
				  std::string phds="PHD";
		  	  std::string filename=Hs+"exactL"+std::to_string(L)+"M"+std::to_string(M)+"t0_"+"1.0"+"gam"+sgam+"omg"+ somg+ "PB" +std::to_string(pb)  + ".bin";
		       std::string filename2=phds+"exactL"+std::to_string(L)+"M"+std::to_string(M)+"t0_"+"1.0"+"gam"+sgam+"omg"+ somg+ "PB" +std::to_string(pb)  + ".bin";
		      std::string filename3=Ts+"exactL"+std::to_string(L)+"M"+std::to_string(M)+"t0_"+"1.0"+"gam"+sgam+"omg"+ somg+ "PB" +std::to_string(pb)  + ".bin";
	  
	          bin_write(filename, Evec);
		  	         bin_write(filename2, Ovec);
		  		 	         bin_write(filename3, Tr);

    	    

//	      Std::cout<<"for beta = "<< beta*mean <<std::endl;
//    	      std::cout<<"for beta/L = "<< beta << " we get "<<std::endl;


//    	     std::cout<<   expvalCan(ensvec, obs,  T)/L <<std::endl;
    // 	     std::cout<< "obs 2" << std::endl; 
 	      // std::cout<<   expvalCan(ensvec, obs2,  T)/L <<std::endl;
	      // 	     		     std::cout<<  v.mean()/L <<std::endl;
   //   Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals, dt);
   // Eigen::MatrixXcd cEVec=HH.cast<std::complex<double>>();
   //  Eigen::VectorXcd newIn=inistate;
   // int i=0;
   //      while(i*dt<1.5)
   //     {

   //     	 TimeEv::timeev_exact(newIn, cEVec, evExp);
  	 
   // 	 //   	 std::complex<double> c=(newIn.adjoint()*(O*newIn))(0);
   // 	 std::complex<double> c=(inistate.adjoint()*(newIn))(0);
  	 
   // 	 i++;

   // 			// 	outputVals(i)=real(c);
   // 			// 	outputVals2(i)=real(c2);
   //      		// outputTime(i)=i*dt;
   // 	 std::cout<< real(c)<< " "<< i*dt<< std::endl;
	 //  	 BOOST_CHECK(std::abs(real(c2)-real(c))<Many_Body::err);
	     //	             }

 //std::cout<< MatrixXd(OBS2) << std::endl;
 // Eigen::MatrixXd OBS22=HH.adjoint()*(OBS2.selfadjointView<Lower>())*HH;
 // Eigen::VectorXd diagobs2=OBS22.diagonal();
	// std::cout<< en(0) << std::endl;
 // std::cout<<diagobs2(0)<<std::endl; 

    //     bin_write("E"+filename, en);
  return 0;
}
 
