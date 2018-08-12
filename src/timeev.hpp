#pragma once
#include "timeev.hpp"
#include "diag.hpp"
#include <complex>
namespace TimeEv{
    using Many_Body::im;
  using namespace Eigen;
  MatrixXcd EigenvalExponent(Eigen::VectorXd& eigenVals, double dt)
  {
    using std::exp;
    Eigen::VectorXcd expEigenVals=(-eigenVals*im*dt);
    MatrixXcd expEigenValsMatrix=MatrixXcd::Zero(eigenVals.rows(), eigenVals.rows());
        for (int i = 0; i < eigenVals.rows(); ++i)
       {
	 expEigenValsMatrix(i, i)=std::exp( expEigenVals(i) );
	
}
	//    expEigenValsMatrix.diagonal()=(expEigenVals);
    return expEigenValsMatrix;
  }
  // template <typename State, typename Matrix >
  void timeev_exact(Eigen::VectorXcd& initialState, Eigen::MatrixXcd& eigenVectors, Eigen::MatrixXcd& expEigenVals)
  {
    	

    //	Eigen::VectorXcd in=initialState.cast<std::complex<double>>();
        

	Eigen::MatrixXcd hamiltonianAdj=Eigen::MatrixXcd(eigenVectors.adjoint());
	 
 	initialState= (hamiltonianAdj*initialState);
	initialState=expEigenVals*initialState;
	initialState=eigenVectors*initialState;
		
			std::cout << "norm E " <<initialState.norm() << std::endl;
    return ; 
  }

  void timeev_lanzcos( Eigen::VectorXcd& initialState,  Eigen::MatrixXcd& hamiltonian, size_t lanczosDim, double dt)
  {
    Eigen::VectorXcd initialState2=initialState;
    Eigen::MatrixXcd Q(hamiltonian.rows(), lanczosDim);
     Many_Body::TriDiagMat tri=Many_Body::Lanczos(hamiltonian, initialState2, lanczosDim, Q);
     	 Eigen::MatrixXd S(lanczosDim, lanczosDim);
     	 Eigen::VectorXd eigenVals(lanczosDim);
     	 Many_Body::diag(tri, S, eigenVals);
	 //    std::cout <<    tri.diagel.rows() << "  " <<    tri.diagel.cols() <<std::endl;


    Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals, dt);
    Eigen::MatrixXcd S2=S.cast<std::complex<double>>();
       Eigen::MatrixXcd Q2=Q.cast<std::complex<double>>();
VectorXcd initialStateTemp(lanczosDim);
// std::cout << "UN " <<(evExp.adjoint()*evExp ).sum()-lanczosDim  << std::endl;
   initialStateTemp=Q2.adjoint()*initialState;
   initialStateTemp= S2.adjoint()*initialStateTemp;
         initialStateTemp=evExp*initialStateTemp;
   //   initialStateTemp=std::exp(std::complex<double> (0,1)*dt)*initialStateTemp;
   // std::cout <<    initialStateTemp.rows() << "  " <<    initialStateTemp.cols() <<std::endl;
   //     std::cout <<    evExp.rows() << "  " <<    evExp.cols() <<std::endl;
//  //	initialStateTemp= evExp*S.col(0);
   	initialStateTemp= S2*initialStateTemp;
 	  	initialState= Q2*initialStateTemp;
		//		initialState=initialState/initialState.norm();
		std::cout << "norm " <<initialState.norm() << std::endl;
    return ; 
  }
}
