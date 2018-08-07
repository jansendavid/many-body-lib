#pragma once
#include<vector>
#include"mkl_lapacke.h"
#include <eigen3/Eigen/Sparse>
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

     template< class T, class Y>  
   void diagherm(T& aH, Y& ev);
  // AFTER APPLYING THIS FUCTION MAT= AH.ADJOINT*H*AH+DIAG
    template< class T, class Y>  
  void diagherm(T& aH, Y& ev)
  {
    MKL_INT N= ev.size();
        MKL_INT LDA= N;
        MKL_INT info= LAPACKE_zheevd(LAPACK_COL_MAJOR, 'V', 'U', N, aH.data(), LDA, ev.data());
    if(info!=0)
      {
  	std::cout << " diagonalization failed" << '\n';
      }



  }


};
