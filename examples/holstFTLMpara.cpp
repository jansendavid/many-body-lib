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
     using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
  const double t1=std::stod(argv[1]);
  const double omg=std::stod(argv[2]);
  const double gamma=std::stod(argv[3]);
   bool PB=true;
  const size_t L=5;
  double beta=1;
  int Ldim=20;
    int M=2;
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
