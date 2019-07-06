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
   int Llead{};
  int Lchain{};
  double tint{};
  double t0{};
  double tl{};
    double V{};
  double dt{};
  double tot{};
  bool PB{};
  std::string sLlead{};
  std::string sLchain{};
  std::string stint{};
  std::string st0{};
  std::string stl{};
  std::string sV{};
  std::string sdt{};
  std::string stot{};
  std::string filename;
  std::string sPB{};
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
      if (vm.count("Ll"))
      {      std::cout << "Llead: " << vm["Ll"].as<int>() << '\n';
      	sLlead="Ll"+std::to_string(vm["Ll"].as<int>());
  	filename+=sLlead;
      }
         if (vm.count("Lc"))
      {      std::cout << "Lchain: " << vm["Lc"].as<int>() << '\n';
      	sLchain="Lc"+std::to_string(vm["Lc"].as<int>());
  	filename+=sLchain;
      }
      	 if (vm.count("t0"))
      {      std::cout << "t0: " << vm["t0"].as<double>() << '\n';
      	st0="t0"+std::to_string(vm["t0"].as<double>()).substr(0, 3);
      	filename+=st0;
      }
      	 	 if (vm.count("tl"))
      {      std::cout << "tl: " << vm["tl"].as<double>() << '\n';
      	stl="tl"+std::to_string(vm["tl"].as<double>()).substr(0, 3);
      		filename+=stl;
      }
      		 if (vm.count("tint"))
      {      std::cout << "tint: " << vm["tint"].as<double>() << '\n';
      	stint="tint"+std::to_string(vm["tint"].as<double>()).substr(0, 3);
      		filename+=stint;
      }
      		 if (vm.count("V"))
      {      std::cout << "V: " << vm["V"].as<double>() << '\n';
      	sV="V"+std::to_string(vm["V"].as<double>()).substr(0, 3);
      		filename+=sV;
      }
      		 if (vm.count("tot"))
      {      std::cout << "tot: " << vm["tot"].as<double>() << '\n';
      	stot="tot"+std::to_string(vm["tot"].as<double>()).substr(0, 3);
      		filename+=stot;
      }
		 if (vm.count("dt"))
      {      std::cout << "dt: " << vm["dt"].as<double>() << '\n';
      	sdt="dt"+std::to_string(vm["dt"].as<double>()).substr(0, 6);
      		filename+=sdt;
      }
       if (vm.count("pb"))
      {
  	std::cout << "PB: " << vm["pb"].as<bool>() << '\n';
  		sPB="PB"+std::to_string(vm["pb"].as<bool>());
  	filename+=sPB;
      }
     }
  }
  catch (const boost::program_options::error &ex)
  {
    std::cerr << ex.what() << '\n';
    return 0;
  }

  int Ltot=2*Llead+Lchain;
  ElectronBasis e( Ltot, 1);
 
		filename+=".bin";
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
      std::cout << "done diag"<<std::endl;

      Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals2, dt);
    Eigen::MatrixXcd cEVec=HH2.cast<std::complex<double>>();
    Eigen::VectorXcd newIn=HH1.col(0);
    std::vector<Eigen::VectorXcd> statesVec(int(Ltot/2), Eigen::VectorXcd::Zero(e.dim));
    for(int i=0; i<int(Ltot/2); i++)
      {

	statesVec[i]=HH1.col(i);

      }
 int i=0;
 auto O=C1+C2;
    std::vector<double> obstebd;
  std::vector<double> time;
        while(i*dt<tot)
       {
	 double sum1{0};
	 double sum2{0};
	   	 i++;
	 for(auto& state : statesVec)
	   {
       	 TimeEv::timeev_exact(state, cEVec, evExp);
  	 
       	  	 std::complex<double> c=im*(state.adjoint()*(O*state))(0);
		 std::complex<double> c2=(state.adjoint()*(H1*state))(0);
		 sum1+=real(c);
		 sum2+=real(c2);
	   }
  		  time.push_back(i*dt);
  		  obstebd.push_back(sum1);
     

       			// 	outputVals(i)=real(c);
       			// 	outputVals2(i)=real(c2);
        		// outputTime(i)=i*dt;
       	 std::cout<< std::setprecision(8)<<sum1<< " dt "<< i*dt<< "  E  "<< sum2<< std::endl;
       	 //	  	 BOOST_CHECK(std::abs(real(c2)-real(c))<Many_Body::err);
       	     	             }

       Many_Body::bin_write("SPextimeFermi"+filename, time);
    Many_Body::bin_write("SPexJFermi"+filename, obstebd);
  return 0;
}
 
