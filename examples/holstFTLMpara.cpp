#include<iostream>
#include <Eigen/Eigenvalues> 
#include"numerics.hpp"
#include"reddm.hpp"
#include"tpoperators.hpp"
#include "files.hpp"
#include"FTLanczos.hpp"
#include<boost/mpi.hpp>
#include<boost/mpi/communicator.hpp>
#include<boost/mpi/environment.hpp>
#include <boost/program_options.hpp>
using namespace Many_Body;
using Mat= Operators::Mat;
namespace mpi = boost::mpi;
using namespace boost::program_options;
int main(int argc, char *argv[])
{


  int M{};
  int L{};
  double t0{};
  double omega{};
  double gamma{};
  bool PB{0};
   int runs={};
   int Ldim={};
    double T{};
    double err{};
  std::string sLdim={};
  std::string sruns={};
std::string sT{};
std::string sM{};
std::string sL{};
std::string st0{};
std::string somega{};
std::string sgamma{};
std::string sPB{};
 std::string serr{};
 std::string filename="FTLM";
   

  try
  {
    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("L", value(&L)->default_value(4), "L")
      ("M", value(&M)->default_value(2), "M")
      ("r", value(&runs)->default_value(20), "r")
      ("Ld", value(&Ldim)->default_value(20), "Ld")
      ("t0", value(&t0)->default_value(1.), "t0")
      ("gam", value(&gamma)->default_value(1.), "gamma")
      ("omg", value(&omega)->default_value(1.), "omega")
      ("T", value(&T)->default_value(1.), "T")
      ("pb", value(&PB)->default_value(false), "PB")
      ("err", boost::program_options::value(&err)->default_value(1E-9), "err");
    boost::program_options::variables_map vm;
    boost::program_options::store(parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
 if (vm.count("help"))
      {std::cout << desc << '\n'; return 0;}
    else{
      if (vm.count("L"))
      {
  	std::cout << "L: " << vm["L"].as<int>() << '\n';
	sL="L"+std::to_string(vm["L"].as<int>());
	filename+=sL;
	
      }
     if (vm.count("M"))
      {
  	std::cout << "M: " << vm["M"].as<int>() << '\n';
	sM="M"+std::to_string(vm["M"].as<int>());
	filename+=sM;
	
      }
     if (vm.count("r"))
      {
  	std::cout << "runs: " << vm["r"].as<int>() << '\n';
	sruns="r"+std::to_string(vm["r"].as<int>());
	filename+=sruns;
	
      }
if (vm.count("Ld"))
      {
  	std::cout << "lanczos dim: " << vm["Ld"].as<int>() << '\n';

	
      }
      if (vm.count("t"))
      {
      	std::cout << "t0: " << vm["t"].as<double>() << '\n';
      	      	st0="t0"+std::to_string(vm["t0"].as<double>()).substr(0, 3);
      	filename+=st0;
      }
       if (vm.count("omg"))
      {
      	std::cout << "omega: " << vm["omg"].as<double>() << '\n';
      	      	somega="omg"+std::to_string(vm["omg"].as<double>()).substr(0, 3);
      	filename+=somega;
      }
       if (vm.count("gam"))
      {
      	std::cout << "gamma: " << vm["gam"].as<double>() << '\n';
      	sgamma="gam"+std::to_string(vm["gam"].as<double>()).substr(0, 3);
      	filename+=sgamma;
      }
       if (vm.count("T"))
      {
      	std::cout << "T: " << vm["T"].as<double>() << '\n';
      }if (vm.count("pb"))
      {
      	std::cout << "PB: " << vm["pb"].as<bool>() << '\n';
      	sPB="PB"+std::to_string(vm["pb"].as<bool>());
      	filename+=sPB;
      }
      		     if (vm.count("err"))
      {      std::cout << "error: " << vm["err"].as<double>() << '\n';
      	 std::stringstream ss;
      	 ss<<vm["err"].as<double>();
      	 serr="err"+ss.str();
      	filename+=serr;
      }
    }
  }
  catch (const error &ex)
  {
    std::cerr << ex.what() << '\n';
    return 0;
  }

     using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
 
  filename+=".bin";
     std::vector<double> Tem;
     std::vector<double> beta;
     for(int i=1; i<11; i++)
       {
	 	  Tem.push_back((0.1*i));
	  beta.push_back(1./(0.1*i));
       }
     
     
   
  std::vector<Mat> obs;
  {
    
    ElectronBasis e( L, 1);
  PhononBasis ph(L, M);
   HolsteinBasis TP(e, ph);
   std::cout<< TP.dim<<std::endl;
  
  
    obs.push_back(Operators::EKinOperatorL(TP, e, t0,PB));
    
         obs.push_back(Operators::NumberOperator(TP, ph, 1,  PB));
         obs.push_back(Operators::EKinOperatorL(TP, e, t0,PB));
	 obs.push_back(Operators::NBosonDOperator(TP, ph, 1, PB));
	 obs[obs.size()-1]+=Operators::NBosonCOperator(TP, ph, 1, PB);
	  obs[0]+=gamma*obs[obs.size()-1];
	  obs[0]+=omega*obs[1];
	
  }


  // mpi::environment env;
  // mpi::communicator world;
  // //  std::vector<double> As(obs.size(), 0);
  // Eigen::MatrixXd As=Eigen::MatrixXd::Zero(beta.size(), obs.size());
  //   Eigen::VectorXd Zs=Eigen::VectorXd::Zero(beta.size());

  //   //  double Z{0.};
  // if(world.rank()==0)
  //   {

  //     Eigen::MatrixXd Astot=Eigen::MatrixXd::Zero(beta.size(), obs.size());
  //     Eigen::VectorXd Zstot=Eigen::VectorXd::Zero(beta.size());
  //         for(int i=0; i<runs/world.size(); i++)
  //      {
  //     auto [Observables, SUMs]=calculate_lanczFT_fast(obs[0], obs, beta, Ldim, err);

  //   	As+=Observables;
  //   	Zs+=SUMs;
  //     }
  //      std::cout<< " process # " << world.rank() << " got meanZ "<< Zs.mean() << std::endl;
  //      for(size_t i=0; i<beta.size(); i++)
  //   	 {
  //   	   reduce(world, Zs(i), Zstot(i), std::plus<double>(), 0);

  //   for(size_t k=0; k<obs.size(); k++)
  //   	   {
  //   	     reduce(world, As(i, k), Astot(i, k), std::plus<double>(), 0);

  //   	   }
  //   	 }

  //     for(size_t i=0; i<beta.size(); i++)
  //   	 {
  //   	   Astot.row(i)/=Zstot(i);
  //   	 }
  //     std::cout<< "Astot  "<<std::endl<< Astot<< std::endl;
  //           for(size_t i=0; i<beta.size(); i++)
  //   	 {
  //   	   std::cout<<" T "<< 1./beta[i] << "  "<<Astot(i, 0)<<" SUM "<< Astot(i, 1)+Astot(i, 2)+Astot(i, 3)*gamma<<std::endl;
  //   	 }
  // 	    	    bin_write("E"+filename, Eigen::VectorXd(Astot.col(0)));
  // 	    bin_write("Nph"+filename,  Eigen::VectorXd(Astot.col(1)));
  // 	    bin_write("EK"+filename, Eigen::VectorXd(Astot.col(2)));
  // 	    bin_write("nX"+filename, Eigen::VectorXd(Astot.col(3)));
  // 	    bin_write("temp"+filename, Tem);
      
  //   }
  // else{

   
  //    for(int i=0; i<runs/world.size(); i++)
  //      {
  //        		auto [Observables, SUMs]=calculate_lanczFT_fast(obs[0], obs, beta, Ldim, err);
  //   	As+=Observables;
  //   	Zs+=SUMs;
	
  //     }
  //   	std::cout<< " process # " << world.rank() << " got meanZ "<< Zs.mean() << std::endl;
  //      for(size_t i=0; i<beta.size(); i++)
  //   	 {
  //   	   reduce(world, Zs[i], std::plus<double>(), 0);

  //   for(size_t k=0; k<obs.size(); k++)
  //   	   {
  //   	     reduce(world, As(i, k), std::plus<double>(), 0);

  //   	   }

    


  //   	}
  // }


  return 0;
}
