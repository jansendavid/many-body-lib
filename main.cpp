#include<iostream>
#include"basis.hpp"
#include"operators.hpp"
#include"tpoperators.hpp"
int main()
{

using namespace Many_Body;
 using namespace Operators;
  const size_t L=5;

   ElectronBasis<L> e1(1);
   
   ElectronBasis<L> e2(1);
      double t1=1;
   double t2=1;
   double u=1;
     std::array<size_t, L> estate;
  estate.fill(0);
  estate[0]=1;
  BosonState<L>  inState(phState);
   TensorProduct<ElectronBasis<L>, ElectronBasis<L>> TP(e1, e2);
Mat E1=Operators::EKinOperatorL(TP, e1, t1);
   Mat E2=Operators::EKinOperatorR(TP, e2, t2);
   Mat C=Operators::CalculateCouplungOperator(TP, e2, u);
    Mat H=E1+E2+ C;
     
  return 0;
}
