
#include<iostream>
#include <Eigen/Eigenvalues> 
#include"numerics.hpp"
#include"reddm.hpp"
#include"tpoperators.hpp"
#include "files.hpp"


using namespace Many_Body;
int main(int argc, char *argv[])
{
  
  size_t M{};
  size_t L{};
  double t0{};
  double omega{};
  double gamma{};
  try
  {
    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("L", value(&L)->default_value(4), "L")
      ("M", value(&M)->default_value(2), "M")
      ("t", value(&t0)->default_value(1.), "t0")
      ("gam", value(&gamma)->default_value(1.), "gamma")
      ("omg", value(&omega)->default_value(1.), "omega");
  


    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help"))
      {std::cout << desc << '\n'; return 0;}
    else{
      if (vm.count("L"))
      {      std::cout << "L: " << vm["L"].as<size_t>() << '\n';
	
      }
     if (vm.count("M,m"))
      {
	std::cout << "M: " << vm["M"].as<size_t>() << '\n';
	
      }
      if (vm.count("t"))
      {
	std::cout << "t0: " << vm["t"].as<double>() << '\n';	
      }
       if (vm.count("omg"))
      {
	std::cout << "omega: " << vm["omg"].as<double>() << '\n';
      }
       if (vm.count("gam"))
      {
	std::cout << "gamma: " << vm["gam"].as<double>() << '\n';
      }
    }
  }
  catch (const error &ex)
  {
    std::cerr << ex.what() << '\n';
    return 0;
  }

    using Mat= Operators::Mat;

  
   using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;

  ElectronBasis e( L, 1);
    std::cout<< e<<std::endl;
  
  PhononBasis ph(L, M);
  std::cout<< ph<<std::endl;
  HolsteinBasis TP(e, ph);
    std::cout<< TP << std::endl;
    std::cout<< TP.dim << std::endl;
        Mat E1=Operators::EKinOperatorL(TP, e, t0, true);
      Mat Ebdag=Operators::NBosonCOperator(TP, ph, gamma, true);
      Mat Eb=Operators::NBosonDOperator(TP, ph, gamma, true);
      Mat Eph=Operators::NumberOperator(TP, ph, omega,  true);
      
      //Mat E=Operators::NumberOperatore(TP, e, 1, false);
      //    std::cout<< HH << std::endl;
      Eigen::VectorXd eigenVals(TP.dim);
      Mat H=E1+Eph  +Ebdag + Eb;

      // state
//         Eigen::VectorXcd evec=Eigen::VectorXcd::Zero(H.rows());
// 	evec(3)=1./std::sqrt(3*2);
// 	evec(7)=std::exp(-im*static_cast<std::complex<double>>(2*pi/3))*1./std::sqrt(3*2);
// 	   evec(11)=std::exp(-im*static_cast<std::complex<double>>(4*pi/3))*1./std::sqrt(3*2);
// 	   evec(21)=1./std::sqrt(3*2);
// 	evec(22)=std::exp(-im*static_cast<std::complex<double>>(2*pi/3))*1./std::sqrt(3*2);
// 	evec(23)=std::exp(-im*static_cast<std::complex<double>>(4*pi/3))*1./std::sqrt(3*2);
// 	std::cout<< " E in "<<(evec.adjoint()*H*evec)(0)<<std::endl;
// 	std::cout<< " norm "<<(evec.adjoint()*evec)(0)<<std::endl;
// auto ov=makeRedDMTP(TP, ph,  0, evec);
//  for(auto g : ov)
//    {
//      std::cout<< g << std::endl;
//      std::cout<< std::endl<<std::endl;
//    }

      Eigen::MatrixXd HH=Eigen::MatrixXd(H);
	  Eigen::MatrixXd HH2=Eigen::MatrixXd(H);
	  Eigen::VectorXd ev=Eigen::VectorXd(TP.dim);
	  bool isDiag=false;
	 //  std::vector<double> Tr={0.0100, 0.1000, 1.000, 2.0000, 5.0000, 10.0000};
	 //  for(auto t: Tr){
	 //      std::string sT=std::string(std::to_string(t)).substr(0,6);
 auto t=T;
	       auto optModes=makeThermalRDMTP(HH,ev,  TP, t, isDiag, 0);
	 //  //	  std::cout << "sum of all eigenvalues "<< optModes.sum()<< std::endl;
	 //      isDiag=true;
	  int n=0;
	  double ent{0};	 
	  for(auto& l : optModes)
	    {
	      	      std::string sn=std::string(std::to_string(n)).substr(0,1);
		      //std::string filename="OML"+std::to_string(L)+"M"+std::to_string(M)+"t0_"+"1.0"+"gam"+sgam+"omg"+ somg+"T"+ sT+ "esec" +sn  + ".bin";
		 std::cout << n << std::endl;
		 for(int i=0; i<l.rows(); i++)
		   {
		     ent-=l(i)*std::log(l(i));
		   }
		 std::cout<< "value "<< " fot T = "<< t <<std::endl;
		 	      std::cout<< l<<std::endl;
			      std::cout << " and sum "<< l.sum()<< std::endl;
			      //bin_write(filename, l);
		 n++;
	    }
	
	  std::cout<< "ENTROPY "<< ent << std::endl;
  return 0;
}
