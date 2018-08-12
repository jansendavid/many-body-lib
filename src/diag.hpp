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
   template <typename T>
   std::vector<MKL_INT> FeastRow(  const Eigen::SparseMatrix<T, Eigen::RowMajor> A )
  {
        MKL_INT N= A.outerSize()+1;
    std::vector<MKL_INT> row(N); 

std::copy (A.outerIndexPtr(), A.outerIndexPtr()+N, row.begin() );
 for(auto& l : row)
   {l+=1;}
//std::for_each(row.begin(), row.end(), [](MKL_INT &n){ n++; });
 return row;
  }

    template <typename T>
    std::vector<MKL_INT> FeastCol(  const Eigen::SparseMatrix<T, Eigen::RowMajor> A )
  {
    MKL_INT N= A.nonZeros();
    std::vector<MKL_INT> col(N); 

    std::copy(A.innerIndexPtr(), A.innerIndexPtr()+N, col.begin() );
    // std::for_each(col.begin(), col.end(), [](int &n){ n++; });
    for(auto& l : col)
   {l+=1;}
 return col;
  }


//   // AFTER APPLYING THIS FUCTION MAT= AH.ADJOINT*H*AH+DIAG

  void diag(Eigen::MatrixXcd& aH, Eigen::VectorXd& ev)
  {
    MKL_INT N= ev.size();
        MKL_INT LDA= N;
        MKL_INT info= LAPACKE_zheevd(LAPACK_COL_MAJOR, 'V', 'U', N, aH.data(), LDA, ev.data());
    if(info!=0)
      {
  	std::cout << " diagonalization failed" << '\n';
      }



  }


  void diag(Eigen::MatrixXd& aH, Eigen::VectorXd& ev)
  {
    MKL_INT N= ev.size();
        MKL_INT LDA= N;
	    MKL_INT info= LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', N, aH.data(), LDA, ev.data());
    if(info!=0)
      {
    	std::cout << " diagonalization failed" << '\n';
      }



  }

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
     	 
	  	 if( std::abs(beta)<0.0001)
     	   {

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
          assert(std::abs( (Q.adjoint()*Q).sum() -(dim+1)) < Many_Body::err);     


     TriDiagMat T(bandTdiag, bandTOff);
    return T;
  }

  
  

  
       


};
