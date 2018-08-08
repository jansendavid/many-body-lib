#pragma once
#include "timeev.hpp"
#include <complex>
namespace TimeEv{
    using Many_Body::im;
  // template <typename State, typename Matrix >
  void timeev_exact(Eigen::VectorXd& initialState, Eigen::MatrixXd& hamiltonian, Eigen::MatrixXd&  observable, Eigen::VectorXd& eigenVals, Eigen::VectorXd& outputVals, Eigen::VectorXd& outputTime, double dt, size_t numberOfSteps)
  {
    	using std::exp;
	//	mat_s.template cast<double>();
	Eigen::VectorXcd in=initialState.cast<std::complex<double>>();
	Eigen::VectorXcd inad=(in.transpose());
	 Eigen::MatrixXcd ham= hamiltonian.cast<std::complex<double>>();
	 Eigen::MatrixXcd hamiltonianAdj=Eigen::MatrixXcd(ham.adjoint());
	 Eigen::MatrixXcd ex(hamiltonianAdj.rows(), hamiltonianAdj.rows());
	 ex.setZero();
    Eigen::VectorXcd expEigenVal=(-eigenVals*im*dt);
 
     // for (int i = 0; i < expEigenVal.rows(); ++i)
     //   {
     // 	 expEigenVal[i]=std::exp( expEigenVal[i] );
	
     //   }
     // ex.diagonal()=expEigenVal;
    //	Eigen::VectorXcd expEigenVal=(-eigenVals*im*dt);
 
     for (int i = 0; i < expEigenVal.rows(); ++i)
       {
	 expEigenVal[i]=std::exp( expEigenVal[i] );
	
       }
     ex.diagonal()=expEigenVal;
    for (size_t i = 0; i < numberOfSteps; ++i)
      {

	 	in= ham*ex*(hamiltonianAdj*in);
		std::complex<double> c=(in.adjoint()*(observable*in))(0);
		assert(imag(c)<Many_Body::err);
				outputVals(i)=real(c);
		outputTime(i)=i*dt;
      }
    return ; 
  }
}
