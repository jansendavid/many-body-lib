#pragma once
#include<vector>
#include"mkl_lapacke.h"
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

namespace Many_Body{
  struct TriDiagMat{
    TriDiagMat(const Eigen::VectorXd& diagel, const Eigen::VectorXd& offDiag): diagel(diagel), offDiag(offDiag){}
    Eigen::VectorXd diagel;
    Eigen::VectorXd offDiag;
    auto dData(){return diagel.data();}
      auto offData(){return offDiag.data();}
  };


  void diag(    Many_Body::TriDiagMat& tri, Eigen::MatrixXd& z, Eigen::VectorXd& ev)
  {
    MKL_INT N= ev.size();
        MKL_INT LDA= N;
        MKL_INT info= LAPACKE_dstedc(LAPACK_COL_MAJOR, 'I', N, tri.dData(), tri.offData(), z.data(), N);
    if(info!=0)
      {
  	std::cout << " diagonalization failed" << '\n';
      }
    ev=tri.diagel;



   }
  template<typename Matrix>
  TriDiagMat Lanczos(Matrix& A, Eigen::VectorXcd& state, const size_t iterations, Eigen::MatrixXcd& Q)
  {
    using namespace Eigen::internal;
    using namespace Eigen;
    	assert(std::abs(state.norm() -1.) < Many_Body::err);




     
     Q.setZero();
     Q.col(0)=state;

        long double beta=1;

    double alpha(0);
    Eigen::VectorXd bandTdiag(iterations);
        Eigen::VectorXd bandTOff(iterations-1);
    bandTdiag.setZero();
        bandTOff.setZero();

    Eigen::VectorXcd qk=state;
    Eigen::VectorXcd qkmin(A.rows());
    qkmin.setZero();
    size_t dim=0;

     for (size_t k = 1; k < iterations; ++k)
       {
     	Eigen::VectorXcd qMiddle=A*qk;
		
     	std::complex<double> c=qk.adjoint()*qMiddle;

      	alpha=real(c);

      	Eigen::VectorXcd  rk=qMiddle - alpha*qk -beta*qkmin;
     	beta=rk.norm();

	bandTdiag(k-1)=alpha;
     	bandTOff(k-1)=beta;

      	qkmin=qk;
	qk=rk/rk.norm();
	     Q.col(k)=qk;
     	 
	  	 if( std::abs(beta)<0.00001)
     	   {
	     std::cout<< "happend "<<std::endl;
	     Eigen::MatrixXcd W=Q;
	           Q.resize(A.rows(), k);
		 for (int i = 0; i < k; ++i)
		   {
		     Q.col(i)=W.col(i);
		     
		   }

		 
    		  bandTdiag.resize( k);
    		  bandTOff.resize( k-1);

		 
     	    break;
     	   }


    	 
	dim=k;   
    	 

             }
    

    
     {
       Eigen::VectorXcd qMiddle=A*qk;
       	std::complex<double> c=qk.adjoint()*qMiddle;
     	alpha=real(c);

       		        bandTdiag(dim)=(alpha);

    }


     //          std::cout<< (Q.adjoint()*Q).sum() -(dim+1)<<std::endl;
           assert(std::abs( (Q.adjoint()*Q).sum() -(dim+1)) < Many_Body::err);     

	  
	  
    TriDiagMat T(bandTdiag, bandTOff);
    return T;
  }

  
  

  
       


};

