#include<iostream>
#include"basis.hpp"
#include"operators.hpp"
#include"diag.h"
#include"files.hpp"
#include"timeev.hpp"
#include"ETH.hpp"
#include <boost/program_options.hpp>
int main(int argc, char *argv[])
{
  using namespace Eigen;
using namespace std;
using namespace Many_Body;
  using Fermi_HubbardBasis= TensorProduct<ElectronBasis, ElectronBasis>;
   using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
  using Mat= Operators::Mat;
  using boost::program_options::value;
   size_t Llead{};
  size_t Lchain{};
  double tint{};
  double t0{};
  double tl{};
    double V{};
  double dt{};
  double tot{};
  bool PB{};
  try
  {
    boost::program_options::options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("Ll", value(&Llead)->default_value(4), "Ll")
      ("Lc", value(&Lchain)->default_value(2), "Lc")
      ("t0", value(&t0)->default_value(1.), "t0")
      ("tl", value(&tl)->default_value(1.), "tl")
      ("V", value(&V)->default_value(1.), "V")
      ("tint", value(&tint)->default_value(1.), "tint")
      ("dt", value(&dt)->default_value(0.1), "dt")
      ("tot", value(&tot)->default_value(1.), "tot")
      ("pb", value(&PB)->default_value(false), "PB");
  


    boost::program_options::variables_map vm;
    boost::program_options::store(parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    if (vm.count("help"))
      {std::cout << desc << '\n'; return 0;}
    else{
      if (vm.count("L"))
      {      std::cout << "L: " << vm["L"].as<size_t>() << '\n';
	
      }
     if (vm.count("M"))
      {
	std::cout << "M: " << vm["M"].as<size_t>() << '\n';
	
      }
      if (vm.count("t0"))
      {
	std::cout << "t0: " << vm["t0"].as<double>() << '\n';	
      }
       if (vm.count("omg"))
      {
	std::cout << "omega: " << vm["omg"].as<double>() << '\n';
      }
       if (vm.count("V"))
      {
	std::cout << "V: " << vm["V"].as<double>() << '\n';
      }
       if (vm.count("dt"))
      {
	std::cout << "dt: " << vm["dt"].as<double>() << '\n';
      }
       if (vm.count("tot"))
      {
	std::cout << "total time: " << vm["tot"].as<double>() << '\n';
      }
       if (vm.count("pb"))
      {
	std::cout << "PB: " << vm["pb"].as<bool>() << '\n';
      }
    }
  }
  catch (const boost::program_options::error &ex)
  {
    std::cerr << ex.what() << '\n';
    return 0;
  }

  //std::vector<int> ee(L, 0);
  // ee[L-1]=1;

  int Ltot=2*Llead+Lchain;
  ElectronBasis e( Ltot, int(Ltot/2));
  
  //std::cout<< e<<std::endl;



   // inistate.setZero();

   // inistate[StateNr]=1;
   // 	          Eigen::VectorXcd i0=inistate;
   //    //   std::cout<< e << std::endl;
		  
            std::cout<< e.dim << std::endl;  
	   
  Mat E1=Operators::EKinOperator(e, tl, PB, 0, Llead-1);
  	   
   Mat E2=Operators::EKinOperator(e, tl, PB, Llead+Lchain, Ltot-1);
  Mat EI1=Operators::EKinOperator(e, tint, PB, Llead-1, Llead);
  Mat EI2=Operators::EKinOperator(e, tint, PB, Llead+Lchain-1, Llead+Lchain);
    Mat NI1=Operators::NumberOperator(e, -V/2, PB,  0, Llead-1);
  Mat NI2=Operators::NumberOperator(e, V/2, PB, Llead+Lchain, Ltot-1);
  // // Mat E4=Operators::EKinOperator(e, tl, PB, Llead+Lchain, Ltot-1);
   Mat E3=Operators::EKinOperator(e, t0, PB, Llead, Llead+Lchain-1);
   Mat C1=Operators::CurrOperator(e, tint, PB, Llead-1, Llead);
   Mat C2=Operators::CurrOperator(e, tint, PB, Llead+Lchain-1, Llead+Lchain);

       Eigen::VectorXd eigenVals1(e.dim);
       Eigen::VectorXd eigenVals2(e.dim);
       Mat H1=E1+E2+E3+ EI1+ EI2;
       Mat H2=E1+E2+E3+ EI1+ EI2 +NI2 +NI1;
       //;
      //    
 
    Eigen::MatrixXd HH2=Eigen::MatrixXd(H2);
       Eigen::MatrixXd HH1=Eigen::MatrixXd(H1);
    //std::cout<< HH<<std::endl;
   //      Eigen::MatrixXd N=Eigen::MatrixXd(Eph);
      Many_Body::diagMat(HH1, eigenVals1);
      Many_Body::diagMat(HH2, eigenVals2);
      //      std::cout<< eigenVals[0]<< std::endl;
   //       Eigen::VectorXd energy(TP.dim);
   //  	   std::vector<double> obs(TP.dim);
   // 	   std::vector<double> obs2(TP.dim);
   // 	   std::vector<double> ensvec(TP.dim);
   // 	   Eigen::MatrixXd PHD2=HH.adjoint()*N.selfadjointView<Lower>()*HH;

      Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals2, dt);
    Eigen::MatrixXcd cEVec=HH2.cast<std::complex<double>>();
    Eigen::VectorXcd newIn=HH1.col(0);
 int i=0;
 auto O=C1+C2;
        while(i*dt<tot)
       {

       	 TimeEv::timeev_exact(newIn, cEVec, evExp);
  	 
       	  	 std::complex<double> c=im*(newIn.adjoint()*(O*newIn))(0);
		 std::complex<double> c2=(newIn.adjoint()*(H1*newIn))(0);
       		 //	 std::complex<double> c=(inistate.adjoint()*(newIn))(0);
  	 
       	 i++;

       			// 	outputVals(i)=real(c);
       			// 	outputVals2(i)=real(c2);
        		// outputTime(i)=i*dt;
       	 std::cout<< std::setprecision(8)<<c<< " dt "<< i*dt<< "  E  "<< c2<< std::endl;
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
 
