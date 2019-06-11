#include<iostream>
#include <eigen3/Eigen/Eigenvalues> 
#include"numerics.hpp"
#include"reddm.hpp"
#include"tpoperators.hpp"
#include "files.hpp"
#include"FTLanczos.hpp"
using namespace Many_Body;
    using Mat= Operators::Mat;
int main(int argc, char *argv[])
{
      int L=5;
  double omega=1;
  double gamma=1;
  double t0=1;
  int M=2;
  double mean= 0.5*omega*L*M;
  double beta=1./mean;
  
  double T=1./beta;
   using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
        PhononBasis g2{ 2, 1};
  ElectronBasis e( L, 1);

  
  PhononBasis ph(L, M);

  HolsteinBasis TP(e, ph);


        Mat E1=Operators::EKinOperatorL(TP, e, t0, false);
       Mat Ebdag=Operators::BosonCOperator(TP, ph, gamma, false);
       std::cout<< "dim "<< TP.dim<< std::endl;
       Mat Eb=Operators::BosonDOperator(TP, ph, gamma, false);
       Mat Eph=Operators::NumberOperator(TP, ph, omega,  false);
      
      //Mat E=Operators::NumberOperatore(TP, e, 1, false);
      //    std::cout<< HH << std::endl;
      Eigen::VectorXd eigenVals(TP.dim);
       	Mat H=E1+Eph +Ebdag + Eb;
	auto O=Operators::NumberOperator(TP, ph, omega,  false);
	std::vector<Mat> v{H, O};
       //       auto HH=Eigen::MatrixXd(H);

      auto ev=Eigen::VectorXd(H.rows());
      //    diagMat(HH, ev);


	auto o=LTLM(H, v, T, 800);


  
	           std::cout<< "for beta/mean = " << beta << std::endl;
       for(auto& l: o)
   	{std::cout<< l/L <<std::endl; }
  
  return 0;
}

