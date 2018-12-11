#include<iostream>
#include"basis.hpp"
#include"operators.hpp"
#include"diag.h"
#include"tpoperators.hpp"
int main(int argc, char *argv[])
{

using namespace Many_Body;
 using namespace Operators;
  const size_t L=4;

   ElectronBasis<L> e1(1);
   
   ElectronBasis<L> e2(L/2);
//    std::cout<< e1<< std::endl;
// std::cout<< e2<< std::endl;
  const double t1=std::stoi(argv[1]);
  const double t11=std::stoi(argv[2]);
  const double t2= std::stod(argv[3]);
  const double t22= std::stod(argv[4]);
  const double u= std::stod(argv[5]);
     std::array<size_t, L> estate;
  estate.fill(0);
  estate[0]=1;

  TensorProduct<ElectronBasis<L>, ElectronBasis<L>> TP(e1, e2);
  //std::cout<< TP<< std::endl;
Mat E1=Operators::EKinOperatorL(TP, e1, t1);
   Mat E2=Operators::EKinOperatorR(TP, e2, t2);
   Mat E11=Operators::EKinOperatorLNNN(TP, e1, t11);
   Mat E22=Operators::EKinOperatorRNNN(TP, e2, t22);
   Mat C=Operators::CalculateCouplungOperator(TP, e2, u);
    Mat H=E1 +E2+ C + E22+E11;
      Eigen::MatrixXd HH=Eigen::MatrixXd(H);
     Eigen::VectorXd en(TP.dim);
     //     std::cout<< HH<<std::endl;
  Many_Body::diagMat(HH, en);
  std::cout<< en << std::endl;
  return 0;
}
