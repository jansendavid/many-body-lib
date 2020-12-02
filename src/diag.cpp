#include"diag.h"
namespace Many_Body{


  //  template <typename T>
//    std::vector<MKL_INT> FeastRow(  const Eigen::SparseMatrix<T, Eigen::RowMajor> A )
//   {
//         MKL_INT N= A.outerSize()+1;
//     std::vector<MKL_INT> row(N); 

// std::copy (A.outerIndexPtr(), A.outerIndexPtr()+N, row.begin() );
//  for(auto& l : row)
//    {l+=1;}
// //std::for_each(row.begin(), row.end(), [](MKL_INT &n){ n++; });
//  return row;
//   }

//     template <typename T>
//     std::vector<MKL_INT> FeastCol(  const Eigen::SparseMatrix<T, Eigen::RowMajor> A )
//   {
//     MKL_INT N= A.nonZeros();
//     std::vector<MKL_INT> col(N); 

//     std::copy(A.innerIndexPtr(), A.innerIndexPtr()+N, col.begin() );
//     // std::for_each(col.begin(), col.end(), [](int &n){ n++; });
//     for(auto& l : col)
//    {l+=1;}
//  return col;
//   }

void diag(std::complex<double>* aH, double* ev, size_t M)
  {
    MKL_INT N= M;
    MKL_INT LDA= N;
    MKL_INT info= LAPACKE_zheevd(LAPACK_COL_MAJOR, 'V', 'U', N, aH, LDA, ev);
    if(info!=0)
      {
  	std::cout << " diagonalization failed" << '\n';
      }



  }


  void diag(double* aH, double* ev, size_t M)
  {
    MKL_INT N= M;
        MKL_INT LDA= N;
	    MKL_INT info= LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', N, aH, LDA, ev);
    if(info!=0)
      {
    	std::cout << " diagonalization failed" << '\n';
      }



  }



    void diagzheev(std::complex<double>* aH, double* ev, size_t M)
  {
    MKL_INT N= M;
        MKL_INT LDA= N;
	    MKL_INT info= LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', N, aH, LDA, ev);
    if(info!=0)
      {
    	std::cout << " diagonalization failed" << '\n';
      }



  }
  void diagzheevr(std::complex<double>* aH,std::complex<double>* z, double* ev, size_t M)
  {
    MKL_INT N= M;
        MKL_INT LDA= N;
	 MKL_INT m= N;
	 char V='V';
	 char A='A';
	 char U='U';
	 // vl and vu are not referenced
	 double vl=0;
	 double vu=0;
	 // vl and vu are not referenced
    	 MKL_INT il= 0;
    	 MKL_INT iu= 0;
    	 // tolerance
    	 double abstol=1e-17;
    // 	 // ldz
	 MKL_INT ldz= M;

	 std::vector<MKL_INT> izups(2*m, 0);
     	 MKL_INT info= LAPACKE_zheevr(LAPACK_COL_MAJOR, V, A, U,  N, aH, LDA,vl, vu, il, iu,abstol, &m,  ev, z, ldz, izups.data());
    if(info!=0)
      {
    	std::cout << " diagonalization failed" << '\n';
      }



  }
      // void diagzheevr(std::complex<double>* aH, double* ev, size_t M)
  // {
  //   MKL_INT N= M;
  //       MKL_INT LDA= N;
  // 	    MKL_INT info= LAPACKE_zheevr(LAPACK_COL_MAJOR, 'V', 'U', N, aH, LDA, ev);
  //   if(info!=0)
  //     {
  //   	std::cout << " diagonalization failed" << '\n';
  //     }



  // }
void diagOnlyEv(std::complex<double>* aH, double* ev, size_t M)
  {
    MKL_INT N= M;
    MKL_INT LDA= N;
    MKL_INT info= LAPACKE_zheev(LAPACK_COL_MAJOR, 'N', 'U', N, aH, LDA, ev);
    if(info!=0)
      {
  	std::cout << " diagonalization failed" << '\n';
      }



  }


  void diagOnlyEv(double* aH, double* ev, size_t M)
  {
    MKL_INT N= M;
        MKL_INT LDA= N;
	    MKL_INT info= LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'U', N, aH, LDA, ev);
    if(info!=0)
      {
    	std::cout << " diagonalization failed" << '\n';
      }



  }
  // void dTri( double* d, double* od, double* z,  size_t M)
  // {
  //   MKL_INT N= M;
  //       MKL_INT LDA= N;
  //       MKL_INT info= LAPACKE_dstedc(LAPACK_COL_MAJOR, 'I', N, d, od, z, N);
  //   if(info!=0)
  //     {
  // 	std::cout << " diagonalization failed" << '\n';
  //     }
  //   return;

  //  }


}
