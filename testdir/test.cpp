
#include"basis.hpp"
#include"operators.hpp"
#include"timeev.hpp"
#include"diag.hpp"
#include "files.hpp"
using namespace Many_Body;
int main()
{
  
 size_t size=3;
  Eigen::VectorXd x = Eigen::VectorXd::Random(size);
  x=x/x.norm();
std::cout<< x<<std::endl;
    Eigen::MatrixXd AA = Eigen::MatrixXd::Random(size, size);
     Eigen::MatrixXd A = AA + AA.transpose();
  std::cout<< A<<std::endl;
  Eigen::MatrixXd B=A;
  Eigen::VectorXd evA(size);
  Eigen::VectorXd evB(size);
  Many_Body::TriDiagMat tri=Many_Body::Lanczos(B, x, size);
  std::cout<< "en"<<std::endl;
    std::cout<< tri.diagel<<std::endl;
      std::cout<< "to"<<std::endl;
        std::cout<< tri.offDiag<<std::endl;
	 std::cout<< "to"<<std::endl;
Eigen::MatrixXd S(size, size);
 S.setZero();
 diag(tri, S, evB);
  diag(A, evA);
  //Many_Body::diag(A, evA);
  for (int i = 0; i < size; ++i)
  {
    std::cout << evA(i) << "  " << evB(i) << '\n';
  }
  
  return 0;
};
