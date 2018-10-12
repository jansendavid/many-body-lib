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
    Eigen::MatrixXcd Q(hamiltonian.rows(), lanczosDim);

      Many_Body::TriDiagMat tri=Many_Body::Lanczos(hamiltonian, initialState2, lanczosDim, Q);

      Eigen::MatrixXd S(Q.cols(), Q.cols());
      Eigen::VectorXd eigenVals(Q.cols());

     	 Many_Body::diag(tri, S, eigenVals);




	 Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> evExp=TimeEv::EigenvalExponent(eigenVals, dt);
    Eigen::MatrixXcd S2=S.cast<std::complex<double>>();
       Eigen::MatrixXcd Q2=Q.cast<std::complex<double>>();
       VectorXcd initialStateTemp(Q.cols());

	 initialStateTemp= evExp*(S2.row(0).transpose());

    	initialStateTemp= S2*initialStateTemp;
	
 	  	initialState= Q2*initialStateTemp;
		
		  assert(std::abs(initialState.norm() -1.) < Many_Body::err);
	
    return ; 
  }
}
