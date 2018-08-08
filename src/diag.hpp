#pragma once
#include<vector>
#include"mkl_lapacke.h"
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
namespace Many_Body{
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


  // AFTER APPLYING THIS FUCTION MAT= AH.ADJOINT*H*AH+DIAG

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
  void Lanczos(Eigen::MatrixXd A, size_t iteration)
  {
    using namespace Eigen;
    BandMatrix<Scalar,Size,Supers,Subs,Options> mat;
  }
  
       


};
