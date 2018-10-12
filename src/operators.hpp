# pragma once

#include"basis.hpp"
#include<cmath>
#include <eigen3/Eigen/Sparse>
#ifdef MOM
using ValType= std::complex<double>; 
#else
using ValType= double;
#endif
namespace Operators{
  using Many_Body::pi;
  using  Many_Body::GetLattice;
  using  Many_Body::Translate;
  using  Many_Body::Position;
  //   using  Many_Body::State();
  using Many_Body::RightId;
  using Many_Body::LeftId;
  
  using Mat= Eigen::SparseMatrix<ValType,Eigen::RowMajor>;

  size_t NextWithBC(size_t i, size_t sites, bool PB=true)
  {
    if(PB)
      {return (i+1)%sites;}
    else{return (i+1);}
  }
  size_t Length( size_t sites, bool PB=true)
  {
    if(PB)
      {return sites;}
    else{return sites-1;}
  }
   
  template<class TotalBasis, bool PB=true>
 Mat NumberOperator(const TotalBasis& totalBasis, const double omega=1.)
 {    using BasisIt= typename TotalBasis::BasisIt;
    using Const_BasisIt= typename TotalBasis::Const_BasisIt;
    size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
  


      for(const auto& tpState : totalBasis)
	{

	 
	  for(size_t i=0; i<sites; i++)
	    {

	    op.coeffRef(Position(tpState), Position(tpState))+=ValType(omega*totalBasis.particlesAt(Id(tpState), i));
       
	   
	  }
	}
      return op;
      
       }
  template<class TotalBasis, bool PB=true >
     Mat EKinOperator(const TotalBasis& totalBasis, double var=1. , size_t m=0)
     {
	 using TpBasisIt= typename TotalBasis::BasisIt;
using Lattice=typename TotalBasis::Lattice;     
    
 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      
      
      // NOT DONE
	  for( auto& tpState : totalBasis)
	    {
	      double currentID=Id(tpState);
	      
	     TpBasisIt it2=totalBasis.find(currentID);
	     for(size_t i=0; i<Length(sites, PB); i++)
                {
	  	  size_t j=NextWithBC(i, sites, PB);

	  	    size_t P1=totalBasis.particlesAt(currentID, i);
	  	    size_t P2=totalBasis.particlesAt(currentID, i);
	  	    if(totalBasis.particlesAt(currentID, i)>0 && totalBasis.particlesAt(currentID, j)<totalBasis.maxParticles)
	      	     {
	     	       Lattice state=GetLattice(*it2);

	  	       //		       TD<state> a(state);
	    	       state[j]=totalBasis.particlesAt(currentID, j) +1;
	     	       state[i]=totalBasis.particlesAt(currentID, i) -1;

	  	       auto it3= totalBasis.find(state.GetId());
	  	       size_t  newStateNr=Position(*(totalBasis.find(state.GetId())));
      op.coeffRef(newStateNr, Position(tpState))-= 0.5*var*(std::sqrt(static_cast<double>(state[i]+1.)))*(std::sqrt(static_cast<double>(state[j])));
	  	       op.coeffRef(Position(tpState), newStateNr )-= 0.5*var*(std::sqrt(static_cast<double>(state[i]+1.)))*(std::sqrt(static_cast<double>(state[j])));

	  	     }
 	     	     }
	   

   	  }
    
  return op;
  }
  
}
