# pragma once

#include"basis.hpp"
#include<cmath>
#include <eigen3/Eigen/Sparse>
#include <type_traits>
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

  template<class TotalBasis, class SubBasis>
  Mat NumberOperator(const TotalBasis& totalBasis, const SubBasis& subBasis, const double omega=1.)
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
  template<typename T>
    size_t CheckSign(const T& state, size_t i, size_t j)
  {
    size_t m=0;
    if(i>j)
      {
  	for(size_t l=j+1; l<i; l++)
  	  {
  	    m+=size_t(state[l]);
  	  }

      }
    else
      	for(size_t l=i+1; l<j; l++)
  	  {
  	    m+=size_t(state[l]);
  	  }
    return m;
  }
  template<class TotalBasis, class SubBasis >
  Mat EKinOperatorL(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , size_t m=0)
     {
       using SubBasisIt= typename SubBasis::BasisIt;
	 using TpBasisIt= typename TotalBasis::BasisIt;
	      using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
     using RightBasisIt= typename TotalBasis::RightBasisIt;
     	      using LeftBasis= typename TotalBasis::LB;     
     using RightBasis= typename TotalBasis::RB;  
using Lattice=typename SubBasis::Lattice;     
    
 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      
      
      
	   for( auto& tpState : totalBasis)
	     {
	       	    for(size_t i=0; i<sites; i++)
                 {
	       LeftBasisIt it2=subBasis.find(LeftId(tpState));	     

	       Lattice state=GetLattice(*it2);
                     size_t j=(i+1)%sites;
	   	    		    if(state[i]==state[j])
   	       	     {
	  // 	     
    		  
   	   	     }
   	   	    else{
		        	       Lattice temp=state;
	  		       state[j]=temp[i];
	  		       state[i]=temp[j];
		 
	  		       size_t signControl=CheckSign(temp, i, j);
		       
	  		       it2= subBasis.find(state.GetId());
	  	       RightBasisIt it3=totalBasis.rbasis.find(RightId(tpState));
	  	       size_t newStateNr= Position(*it3)*totalBasis.lbasis.dim +Position(*it2);
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= var;}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= var;}
		    }

	   	     }
 	      	     }
	   

   	  
    
  return op;
  }
 template<class TotalBasis, class SubBasis >
  Mat EKinOperatorLNNN(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , size_t m=0)
     {
       using SubBasisIt= typename SubBasis::BasisIt;
	 using TpBasisIt= typename TotalBasis::BasisIt;
	      using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
     using RightBasisIt= typename TotalBasis::RightBasisIt;
     	      using LeftBasis= typename TotalBasis::LB;     
     using RightBasis= typename TotalBasis::RB;  
using Lattice=typename SubBasis::Lattice;     
    
 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      
      
      
	   for( auto& tpState : totalBasis)
	     {
	       	    for(size_t i=0; i<sites; i++)
                 {
	       LeftBasisIt it2=subBasis.find(LeftId(tpState));	     

	       Lattice state=GetLattice(*it2);
                     size_t j=(i+2)%sites;
	   	    		    if(state[i]==state[j])
   	       	     {
	  // 	     
    		  
   	   	     }
   	   	    else{
		        	       Lattice temp=state;
	  		       state[j]=temp[i];
	  		       state[i]=temp[j];
		 
	  		       size_t signControl=CheckSign(temp, i, j);
		       
	  		       it2= subBasis.find(state.GetId());
	  	       RightBasisIt it3=totalBasis.rbasis.find(RightId(tpState));
	  	       size_t newStateNr= Position(*it3)*totalBasis.lbasis.dim +Position(*it2);
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= var;}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= var;}
		    }

	   	     }
 	      	     }
	   

   	  
    
  return op;
  }


  template<class TotalBasis, class SubBasis >
  Mat EKinOperatorRNNN(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , size_t m=0)
     {
       using SubBasisIt= typename SubBasis::BasisIt;
       using TpBasisIt= typename TotalBasis::BasisIt;
       using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
       using RightBasisIt= typename TotalBasis::RightBasisIt;
       using LeftBasis= typename TotalBasis::LB;     
       using RightBasis= typename TotalBasis::RB;  
       using Lattice=typename SubBasis::Lattice;     
    
 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      
      
      
	  for( auto& tpState : totalBasis)
	    {
	      	    for(size_t i=0; i<sites; i++)
                {
	       RightBasisIt it2=subBasis.find(RightId(tpState));	     

	       Lattice state=GetLattice(*it2);
                     size_t j=(i+2)%sites;
	   	    		    if(state[i]==state[j])
   	       	     {


   	   	     }
   	   	    else{
		      Lattice temp=state;
	  		       state[j]=temp[i];
	  		       state[i]=temp[j];
		 
	  		       size_t signControl=CheckSign(temp, i, j);
		       
	  		       it2= subBasis.find(state.GetId());
	  	       LeftBasisIt it3=totalBasis.lbasis.find(LeftId(tpState));
	  	       size_t newStateNr= Position(*it2)*totalBasis.lbasis.dim +Position(*it3);
	  	              if(signControl%2==0)
   	  	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= var;}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= var;}
    	  


	   	     }
 	     	     }
	   

   	  }
    
  return op;
  }

  template<class TotalBasis, class SubBasis >
  Mat EKinOperatorR(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , size_t m=0)
     {
       using SubBasisIt= typename SubBasis::BasisIt;
       using TpBasisIt= typename TotalBasis::BasisIt;
       using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
       using RightBasisIt= typename TotalBasis::RightBasisIt;
       using LeftBasis= typename TotalBasis::LB;     
       using RightBasis= typename TotalBasis::RB;  
       using Lattice=typename SubBasis::Lattice;     
    
 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      
      
      
	  for( auto& tpState : totalBasis)
	    {
	      	    for(size_t i=0; i<sites; i++)
                {
	       RightBasisIt it2=subBasis.find(RightId(tpState));	     

	       Lattice state=GetLattice(*it2);
                     size_t j=(i+1)%sites;
	   	    		    if(state[i]==state[j])
   	       	     {


   	   	     }
   	   	    else{
		      Lattice temp=state;
	  		       state[j]=temp[i];
	  		       state[i]=temp[j];
		 
	  		       size_t signControl=CheckSign(temp, i, j);
		       
	  		       it2= subBasis.find(state.GetId());
	  	       LeftBasisIt it3=totalBasis.lbasis.find(LeftId(tpState));
	  	       size_t newStateNr= Position(*it2)*totalBasis.lbasis.dim +Position(*it3);
	  	              if(signControl%2==0)
   	  	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= var;}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= var;}
    	  


	   	     }
 	     	     }
	   

   	  }
    
  return op;
  }
  template<typename TotalBasis, typename Basis>
    Mat CalculateCouplungOperator( const TotalBasis& totalBasis, const Basis& subBasis,   const double u=1.)
    {

      using BasisIt= typename TotalBasis::BasisIt;
    using Const_BasisIt= typename TotalBasis::Const_BasisIt;
    using LeftLattice=typename TotalBasis::LeftLattice;
    using RightLattice=typename TotalBasis::RightLattice; 
    size_t dim=totalBasis.dim;
    Mat op(dim, dim);
    op.setZero(); 
    size_t sites=totalBasis.sites;
	
      for(const auto& tpState : totalBasis)
  	{
auto it21=totalBasis.lbasis.find(LeftId(tpState));
	  LeftLattice lState=GetLattice(*it21);
	  auto it22=totalBasis.rbasis.find(RightId(tpState));
	  RightLattice rState=GetLattice(*it22);
	  // RightLattice rState=GetLattice(RightId(tpState));
  	  for(int i=0; i<sites; i++)
	    {
	      op.coeffRef(Position(tpState), Position(tpState))+= (u*static_cast<double>(lState[i])*static_cast<double>(rState[i]));
	  }

       
	   

  	}
      return op;
      
       }

  
}
