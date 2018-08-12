
#include"basis.hpp"
#include"operators.hpp"
#include"timeev.hpp"
#include"diag.hpp"
#include "files.hpp"
using namespace Many_Body;
int main()
{
  
    size_t size=2;
  Eigen::VectorXcd x = Eigen::VectorXcd::Random(size);
  x=x/x.norm();
  std::cout<< x.norm()<<std::endl;
std::cout<< x<<std::endl;
    Eigen::MatrixXcd AA = Eigen::MatrixXcd::Random(size, size);
    Eigen::MatrixXcd A = AA.adjoint()*AA;

std::cout<< A<<std::endl;
  Eigen::MatrixXcd B=A;
  Eigen::VectorXd evA(size);
  Eigen::VectorXd evB(size);
  Eigen::MatrixXcd Q(size, size);
s  Many_Body::TriDiagMat tri=Many_Body::Lanczos(B, x, size, Q);
  std::cout<< "D " << tri.diagel<<std::endl;
  std::cout<< "O " << tri.offDiag<<std::endl;  
Eigen::MatrixXd S(size, size);
  S.setZero();
  diag(tri, S, evB);
   diag(A, evA);
        std::cout << '\n';
  for (int i = 0; i < size; ++i)
  {
    //    std::cout << std::abs(evA(i)-evB(i)) << '\n';
        std::cout << evB(i)<< '\n';
  }
  
  
  return 0;
};
