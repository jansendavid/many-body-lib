
#include<iostream>
#include"basis.hpp"
#include"operators.hpp"
int main()
{
  using namespace Many_Body;
    using namespace Eigen;
  using Mat=Operators::Mat;
   const size_t L=3;
   ElectronBasis<L> e(2);
   std::cout << e << std::endl;
   Mat NR=Operators::NumberOperator(e);
   std::cout << NR;
   Mat EK=Operators::EKinOperator(e);
   // Operators::EKinOperator<ElectronBasis<L>> n2(e);
    std::cout << EK;
  return 0;
};
