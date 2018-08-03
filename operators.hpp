# pragma once
#include<armadillo>
#include"basis.hpp"
#include<Eigen/Sparse>
namespace Operators{
  using Many_Body::pi;
  using  Many_Body::GetLattice;
  using  Many_Body::Translate;
  using  Many_Body::Position;
  //   using  Many_Body::State();
  using Many_Body::RightId;
  using Many_Body::LeftId;
  using Mat= Eigen::SparseMatrix<std::complex<double>>;
    template<typename T> // Help class for types
  class TD;
  template<class OtherOp=Mat>
  class Operator{
    // using BasisIt= typename TotalBasis::BasisIt;
    // using LeftBasis= typename TotalBasis::LeftBasisType;
  public:

    Operator(const OtherOp& otherOperator): mat(otherOperator) { };


    friend std::ostream& operator<<(std::ostream& os,  const Operator& otherOperator)
    {
   
      
      os << otherOperator.mat;
       	  
      return os;
    }
    Mat mat;
 
 
  private:
    size_t L_;
  

    size_t dim_;
  };
 template<class TotalBasis>
  class NumberOperator{
    using BasisIt= typename TotalBasis::BasisIt;
    using Const_BasisIt= typename TotalBasis::Const_BasisIt;
  public:
   NumberOperator( const TotalBasis& totalBasis, const double omega=1. ): mat(totalBasis.dim, totalBasis.dim), L_(totalBasis.sites), dim_(totalBasis.dim)
    {
      mat.setZero();
  


      for(const auto& tpState : totalBasis)
	{

	 
	  for(size_t i=0; i<L_; i++){

	    mat.coeffRef(Position(tpState), Position(tpState))+=static_cast<std::complex<double>> (omega*totalBasis.particlesAt(Id(tpState), i));
       
	   
	  }
	}
       	 
      
    
    }
    template<typename T, typename Y>
    
    friend  Operator<Mat> operator+( const T& m, const Y& n);

    friend std::ostream& operator<<(std::ostream& os,  NumberOperator& numberOperator)
    {
   
      
      os << numberOperator.mat;
       	  
      return os;
    }

    Mat mat;
  private:
    size_t L_;


    size_t dim_;
  };
     template<class TotalBasis >
  class EKinOperator{
    using TpBasisIt= typename TotalBasis::BasisIt;
using Lattice=typename TotalBasis::Lattice;     
    
  public:
     EKinOperator( const TotalBasis& totalBasis, double var=1. , size_t m=0):  mat(totalBasis.dim, totalBasis.dim), L_(totalBasis.sites), dim_(totalBasis.dim)
    {

      mat.setZero();
      
      
	  for( auto& tpState : totalBasis)
	    {
	      double currentID=Id(tpState);
	      
	     TpBasisIt it2=totalBasis.find(currentID);
	       for(size_t i=0; i<L_-1; i++)
                {
                    size_t j=(i+1);
		   
		    if(totalBasis.particlesAt(currentID, i)>0 && totalBasis.particlesAt(currentID, i)<totalBasis.maxParticles)
	      	     {
	     	       Lattice state=GetLattice(*it2);

		       //		       TD<state> a(state);
	    	       state[j]=totalBasis.particlesAt(currentID, j) +1;
	     	       state[i]=totalBasis.particlesAt(currentID, i) -1;
		       //		       newStateNr=Position(totalBasis.find)
	     // mat.coeffRef(newStateNr, Position(tpState))-= var*value;
	     	     }
	   

  	  }
    }
    }
    template<typename T, typename Y>
    friend  Operator<Mat> operator+( const T& m, const Y& n);	

        friend std::ostream& operator<<(std::ostream& os,  EKinOperator& ekinOperator)
    {
   
      
       os << ekinOperator.mat;
       	  
       return os;
     }

       Mat mat;
   private:
     
    size_t L_;
    

     size_t dim_;
 };
  
}
