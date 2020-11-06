#pragma once
#include<complex>
#include<vector>
#include"mkl_lapacke.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

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
        MKL_INT info= LAPACKE_dstedc(LAPACK_COL_MAJOR, 'I', N, tri.dData(), tri.offData(), z.data(), N);
    if(info!=0)
      {
  	std::cout << " diagonalization failed" << '\n';
      }
    ev=tri.diagel;



   }
  template<typename Vector, typename Matrix>
  struct threeLVec
  {
    threeLVec(size_t dim, const Matrix& A): A(A), beta(1), alpha(0) {
      Vector vecOne=Eigen::VectorXcd::Zero(dim);
    Vector vecTwo=Eigen::VectorXcd::Zero(dim);
    Vector vecThree=Eigen::VectorXcd::Zero(dim);
    }
    const Matrix& A;
    double beta{1};
    double alpha{};
    Vector vecOne;
    Vector vecTwo;
    Vector vecThree;

void iterate()
    {
      	vecThree=A*vecTwo;	
	//	std::complex<double> c=qk.adjoint()*qMiddle;
	std::complex<double> c=(vecTwo).adjoint()*vecThree;

      	alpha=real(c);

	//	Eigen::VectorXcd  rk=qMiddle - alpha*qk -beta*qkmin;

	Eigen::VectorXcd  rk=vecThree - alpha*vecTwo -beta*vecOne;
	
     	beta=rk.norm();		
	vecOne=vecTwo;

	//qk=rk/rk.norm();
	vecTwo=rk/rk.norm();
    }
    void lastIterate()
    {
      	vecThree=A*vecTwo;	

	std::complex<double> c=(vecTwo).adjoint()*vecThree;

      	alpha=real(c);


    }
  
    
 
  };
  template<typename Matrix>
  TriDiagMat Lanczos3Vec(Matrix& A, Eigen::VectorXcd& state, const size_t iterations)
  {
    using namespace Eigen::internal;
    using namespace Eigen;
    	assert(std::abs(state.norm() -1.) < Many_Body::err);



	threeLVec< Eigen::VectorXcd, Matrix>  threeLanczVec(A.rows(), A);
	threeLanczVec.vecTwo=state;
        long double beta=1;

    double alpha(0);
    Eigen::VectorXd bandTdiag(iterations);
        Eigen::VectorXd bandTOff(iterations-1);
    bandTdiag.setZero();
        bandTOff.setZero();

    Eigen::VectorXcd qk=state;
    Eigen::VectorXcd qkmin(A.rows());
    qkmin.setZero();
    threeLanczVec.vecOne=qkmin;
    size_t dim=0;

     for (size_t k = 1; k < iterations; ++k)
       {
	 threeLanczVec.iterate();

      	alpha=threeLanczVec.alpha;
     	beta=threeLanczVec.beta;

	bandTdiag(k-1)=alpha;
     	bandTOff(k-1)=beta;
     	 
	  	 if( std::abs(beta)<0.00001)
     	   {
		 
    		  bandTdiag.resize( k);
    		  bandTOff.resize( k-1);

		 
     	    break;
     	   }
    	 dim=k;   
    	 }
     
     
    

    
     {

       threeLanczVec.lastIterate();
       bandTdiag(dim)=threeLanczVec.alpha;

    }




	  
	  
    TriDiagMat T(bandTdiag, bandTOff);
    return T;
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
     	 
	  	 if( std::abs(beta)<0.0001)
     	   {
	     Eigen::MatrixXcd W=Q;
	           Q.resize(A.rows(), k);
		 for (size_t i = 0; i < k; ++i)
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
     //             assert(std::abs( (Q.adjoint()*Q).sum() -(dim+1)) < Many_Body::err);     
     std::cout<<"Qdag Q sum "<< (Q.adjoint()*Q).sum()<<std::endl;
	  
	  
     TriDiagMat T(bandTdiag, bandTOff);
    return T;
  }

  
  template<typename Vector, typename Vector2, typename Matrix>
  auto lanczTrafo(const  Vector & state, const Vector2 & iniState, size_t iteration, const Matrix& A)
    ->Eigen::VectorXcd
  {
    Many_Body::threeLVec< Eigen::VectorXcd, Matrix>  threeLanczVec(A.rows(), A);
    Eigen::VectorXcd newVec=Eigen::VectorXcd::Zero(A.rows());
    threeLanczVec.vecTwo=iniState;
    Eigen::VectorXcd qkmin(A.rows());
    qkmin.setZero();
    threeLanczVec.vecOne=qkmin;
    for(size_t i=0; i<iteration; i++)
      {
	newVec+=threeLanczVec.vecTwo*state(i);
	threeLanczVec.iterate();
	
      }
    
    return newVec;
  }
  
       


};

