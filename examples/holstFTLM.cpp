#include<iostream>
#define EIGEN_USE_MKL_ALL
#include <Eigen/Eigenvalues> 
#include"numerics.hpp"
#include"reddm.hpp"
#include"tpoperators.hpp"
#include "files.hpp"
#include"FTLanczos.hpp"
using namespace Many_Body;
    using Mat= Operators::Mat;
int main(int argc, char *argv[])
{
   using Mat= Operators::Mat;
    int L=3;
  double omega=1;
  double gamma=1;
  double t0=1;
  int M=2;
  double mean=0.5*omega*L*M;;
  
  double beta=1;
  
  double T=1./beta;
   using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
      
  ElectronBasis e( L, 1);

  
  PhononBasis ph(L, M);

  HolsteinBasis TP(e, ph);


        Mat E1=Operators::EKinOperatorL(TP, e, t0, true);
       Mat Ebdag=Operators::NBosonCOperator(TP, ph, gamma, true);
       std::cout<< "dim "<< TP.dim<< std::endl;
       Mat Eb=Operators::NBosonDOperator(TP, ph, gamma, true);
       Mat Eph=Operators::NumberOperator(TP, ph, omega,  true);
      
      Mat N=Operators::NumberOperator(TP, ph, 1,  true);
      //    std::cout<< HH << std::endl;
      Eigen::VectorXd eigenVals(TP.dim);
       	Mat H=E1+Eph +Ebdag + Eb;
       Mat O=N;
       std::vector<Mat> v{H, O};
       //       auto HH=Eigen::MatrixXd(H);

      auto ev=Eigen::VectorXd(H.rows());
      //    diagMat(HH, ev);

      auto  o=FTLM(H, v, T, 800);

       std::cout<< "for beta/mean = " << beta << std::endl;
       for(auto& l: o)
   	{std::cout<< l<<std::endl; }
  
  return 0;
}
