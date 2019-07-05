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
using namespace Many_Body;
    using Mat= Operators::Mat;
namespace mpi = boost::mpi;
int main(int argc, char *argv[])
{
    size_t M{};
  size_t L{};
  double t0{};
  double omega{};
  double gamma{};
    bool PB{};
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
    ("pb", value(&PB)->default_value(true), "PB");
  


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
                     if (vm.count("pb"))
      {
	std::cout << "PB: " << vm["pb"].as<bool>() << '\n';
      }
    }
  }
  catch (const error &ex)
  {
    std::cerr << ex.what() << '\n';
    return 0;
  }

     using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;

   bool PB=true;
 
  double beta=1;
  int Ldim=20;
 
    Mat H;
    Mat N;
  std::vector<Mat> obs;
  {
    
    ElectronBasis e( L, 1);
  PhononBasis ph(L, M);
  HolsteinBasis TP(e, ph);
  std::cout<< TP.dim<<std::endl;
        Mat E1=Operators::EKinOperatorL(TP, e, t1,PB);
       Mat Ebdag=Operators::BosonCOperator(TP, ph, gamma, PB);
       Mat Eb=Operators::BosonDOperator(TP, ph, gamma, PB);
      Mat Eph=Operators::NumberOperator(TP, ph, omg,  PB);

      Eigen::VectorXd eigenVals(TP.dim);
       H=E1+Eb+Eph+Ebdag;
        N=Eph/omg;
       obs.push_back(H);
       obs.push_back(N);
  }
  int runs=100;

  mpi::environment env;
  mpi::communicator world;
  std::vector<double> As(obs.size(), 0);

  double Z{1.};
  if(world.rank()==0)
    {
      std::vector<double> Astot(obs.size(), 0);
      double Ztot{0.};
       for(int i=0; i<runs/world.size(); i++)
      {
  	    auto tup=calculate_lanczFT(obs[0], obs, beta, Ldim);
  auto obs=std::get<0>(tup);
         auto Zt=std::get<1>(tup);
  	
  	 for(int k=0; k<obs.size(); k++)
  	   {
  	     As[k]+=obs[k];
  	   }

        Z+=Zt;  
      }
             reduce(world, Z, Ztot, std::plus<double>(), 0);

    for(int k=0; k<obs.size(); k++)
    	   {
      reduce(world, As[k], Astot[k], std::plus<double>(), 0);

    	   }
    std::cout<< " process # " << world.rank() << " got A " << As[0] << std::endl;
  

      for(auto l: Astot)
	{	std::cout<< " total is "<< l/Ztot<< std::endl;}
	
    }
  else{

   
    for(int i=0; i<runs/world.size(); i++)
      {
  	    auto tup=calculate_lanczFT(obs[0], obs, beta, Ldim);
  auto obs=std::get<0>(tup);
         auto Zt=std::get<1>(tup);
  	 //	std::transform(obstot.begin(), obstot.end(), obs.begin(), obs.end(), std::plus<double>());
  	 for(int k=0; k<obs.size(); k++)
  	   {
  	     As[k]+=obs[k];
  	   }

        Z+=Zt;  
      }
        std::cout<< " process # " << world.rank() << " got A " << As[0] << std::endl;
    reduce(world, Z, std::plus<double>(), 0);
    for(int k=0; k<obs.size(); k++)
  	   {
  	     reduce(world, As[k], std::plus<double>(), 0);
  	   }
  }


  return 0;
}
