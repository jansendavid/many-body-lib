# pragma once
#include<armadillo>
#include"basis.hpp"
namespace operators{
  typedef  arma::mat oper;
  template<class Basis>
  class numberoperator{
    
  public:
    numberoperator(const Basis& aBasis ): mat(aBasis.dim, aBasis.dim)
    {
      mat.zeros();
    }
    numberoperator(const Basis& aBasis, size_t i ): mat(aBasis.dim, aBasis.dim)
    {
       mat.zeros();
       for( typename Basis::basisit it1= aBasis.basis.begin(); it1!= aBasis.basis.end(); ++it1)
	 {
	 }
      
    }
  private:
    oper mat;
   
  };
  
}
