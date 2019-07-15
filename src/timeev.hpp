#pragma once
#include "lanzcos.hpp"
#include "diag.h"
#include <complex>
namespace TimeEv{
    using Many_Body::im;
  using namespace Eigen;
  Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> EigenvalExponent(Eigen::VectorXd& eigenVals, double dt)
  {
    using std::exp;
    Eigen::VectorXcd expEigenVals=(-eigenVals*im*dt);
    Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> expEigenValsMatrix(eigenVals.rows(), eigenVals.rows());
        for (int i = 0; i < eigenVals.rows(); ++i)
       {
	 expEigenValsMatrix.coeffRef(i, i)=std::exp( expEigenVals(i) );
	
}
	//    expEigenValsMatrix.diagonal()=(expEigenVals);
    return expEigenValsMatrix;
  }
  // template <typename State, typename Matrix >
  void timeev_exact(Eigen::VectorXcd& initialState, Eigen::MatrixXcd& eigenVectors, Eigen::MatrixXcd& expEigenVals)
  {
    	

    //	Eigen::VectorXcd in=initialState.cast<std::std::complex<double>>();
        

	Eigen::MatrixXcd hamiltonianAdj=Eigen::MatrixXcd(eigenVectors.adjoint());
	 
 	initialState= (hamiltonianAdj*initialState);
	initialState=expEigenVals*initialState;
	initialState=eigenVectors*initialState;
		assert(std::abs(initialState.norm() -1.) < Many_Body::err);		

    return ; 
  }

  template<typename Matrix>
  void timeev_lanzcos( Eigen::VectorXcd& initialState,  Matrix& hamiltonian, size_t lanczosDim, double dt)
  {
    
    Eigen::VectorXcd initialState2=initialState;


      Many_Body::TriDiagMat tri=Many_Body::Lanczos(hamiltonian, initialState2, lanczosDim);

      Eigen::MatrixXd S(lanczosDim, lanczosDim);
      Eigen::VectorXd eigenVals(lanczosDim);
     	 Many_Body::diag(tri, S, eigenVals);
	 Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> evExp=TimeEv::EigenvalExponent(eigenVals, dt);
	 Eigen::MatrixXcd S2=S.cast<std::complex<double>>();

	 VectorXcd initialStateTemp(lanczosDim);
       initialStateTemp= evExp*(S2.row(0).transpose());
    	initialStateTemp= S2*initialStateTemp;
	
	initialState=Many_Body::lanczTrafo(initialStateTemp, initialState, lanczosDim, hamiltonian);
	assert(std::abs(initialState.norm() -1.) < Many_Body::err);
	
    return ; 
  }
}
