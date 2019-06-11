#include<iostream>
#include <eigen3/Eigen/Eigenvalues> 
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
int main(int argc, char *argv[])
{
  mpi::environment env;
  mpi::communicator world;
  double A{0};
  double Z{0};
  if(world.rank()==0)
    {
      double Atot{0};
      double Ztot{0};
      reduce(world, A, Atot, std::plus<double>(), 0);
      reduce(world, Z, Ztot, std::plus<double>(), 0);
      std::cout<< " total is "<< Atot/Ztot<< std::endl;
    }
  else{
    std::cout<< " process # " << world.rank() << " got A " << A << std::endl;
    reduce(world, A, std::plus<double>(), 0);
    reduce(world, Z, std::plus<double>(), 0);
  }

  //  using Mat= Operators::Mat;
  //   int L=5;
  // double omega=1;
  // double gamma=1;
  // double t0=1;
  // int M=2;
  // double mean= 0.5*omega*L*M;
  // double beta=1./mean;
  
  // double T=1./beta;
  //  using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
  //       PhononBasis g2{ 2, 1};
  // ElectronBasis e( L, 1);

  
  // PhononBasis ph(L, M);

  // HolsteinBasis TP(e, ph);


  //       Mat E1=Operators::EKinOperatorL(TP, e, t0, false);
  //      Mat Ebdag=Operators::BosonCOperator(TP, ph, gamma, false);
  //      std::cout<< "dim "<< TP.dim<< std::endl;
  //      Mat Eb=Operators::BosonDOperator(TP, ph, gamma, false);
  //      Mat Eph=Operators::NumberOperator(TP, ph, omega,  false);
      
  //     Mat N=Operators::NumberOperator(TP, ph, 1,  true);
  //     //    std::cout<< HH << std::endl;
  //     Eigen::VectorXd eigenVals(TP.dim);
  //      	Mat H=E1+Eph +Ebdag + Eb;
  //      Mat O=N;
  //      std::vector<Mat> v{H, O};
  //      //       auto HH=Eigen::MatrixXd(H);

  //     auto ev=Eigen::VectorXd(H.rows());
  //     //    diagMat(HH, ev);

  //     auto  o=FTLM(H, v, T, 800);

  //      std::cout<< "for beta/mean = " << beta << std::endl;
  //      for(auto& l: o)
  //  	{std::cout<< l/L <<std::endl; }
  
  return 0;
}
