 #ifndef LANCZOSHELPER_H
#define LANCZOSHELPER_H
#include"basis.hpp"
#include"lanzcos.hpp"
#include"diag.h"
#include<random>
template<typename Vector, typename Ham>
auto estimateLD(int& lanczosDim, Ham& hamiltonian, Vector& iniState, double err)
  ->std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
{
Eigen::MatrixXd S(lanczosDim, lanczosDim);
 Eigen::VectorXd eigenVals(lanczosDim);
    double diss=1000;
   double Eold=100;
  while(diss>err)
     {
       
       lanczosDim+=1;
 S= Eigen::MatrixXd(lanczosDim, lanczosDim);
     eigenVals=Eigen::VectorXd(lanczosDim);
    
    Many_Body::TriDiagMat tri=Many_Body::Lanczos3Vec(hamiltonian, iniState, lanczosDim);
   
     Many_Body::diag(tri, S, eigenVals);
     diss=std::abs(Eold-eigenVals(0));
     Eold=eigenVals(0);
     }
     return {S, eigenVals};
}
template<typename Vector, typename Ham>
auto estimateLDwithQ(int& lanczosDim, Ham& hamiltonian, Vector& iniState, double err)
  ->std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd, Eigen::VectorXd>
{
Eigen::MatrixXd S(lanczosDim, lanczosDim);
 Eigen::MatrixXcd Q(hamiltonian.rows(), lanczosDim);
 Eigen::VectorXd eigenVals(lanczosDim);
    double diss=1000;
   double Eold=100;
  while(diss>err)
     {
       lanczosDim+=1;
 S= Eigen::MatrixXd(lanczosDim, lanczosDim);
 Q=Eigen::MatrixXcd(hamiltonian.rows(), lanczosDim);
 eigenVals=Eigen::VectorXd(lanczosDim);
    
     Many_Body::TriDiagMat tri=Many_Body::Lanczos(hamiltonian, iniState, lanczosDim, Q);
   
     Many_Body::diag(tri, S, eigenVals);
     diss=std::abs(Eold-eigenVals(0));
     Eold=eigenVals(0);
     }
  return {Q, S, eigenVals};
}


#endif /* LANCZOSHELPER_H */
