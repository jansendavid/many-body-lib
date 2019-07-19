# pragma once
//#define EIGEN_USE_MKL_ALL
#include"basis.hpp"
#include<cmath>

#include <Eigen/Sparse>
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
   
  template<class TotalBasis>
  Mat NumberOperator(const TotalBasis& totalBasis, const double omega=1., bool PB=0, int start=0, int stop=0)
 {
   
 
    size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
 stop = (stop!=0) ? stop : Operators::Length( sites, PB);    


      for(const auto& tpState : totalBasis)
	{

	 
	  for(size_t i=start; i<=stop; i++)
	    {

	    op.coeffRef(Position(tpState), Position(tpState))+=ValType(omega*totalBasis.particlesAt(Id(tpState), i));
       
	   
	  }
	}
      return op;
      
       }
  template<class TotalBasis>
  Mat EKinOperator(const TotalBasis& totalBasis, double var=1. , bool PB=0, int start=0, int stop=0 )
     {
	 using TpBasisIt= typename TotalBasis::BasisIt;
using Lattice=typename TotalBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      stop = (stop!=0) ? stop : Operators::Length( sites, PB);       
      std::cout<< "stop "<< stop<<std::endl;
   for( auto& tpState : totalBasis)
	     {

	       	    for(size_t i=start; i<stop; i++)
                 {

	       

		   Lattice state=GetLattice(tpState);
                     size_t j=Operators::NextWithBC(i, sites, PB);
	   	    		    if(state[i]==state[j])
   	       	     {
	  // 	     
    		  
   	   	     }
   	   	    else{
		      

		      Lattice temp=state;
		      
		       state.setPartNr(j, temp[i]);
		      state.setPartNr(i, temp[j]);
		     		 
	  	     size_t signControl=CheckSign(temp, i, j);
		       
		    auto it2 = totalBasis.find(state.GetId());
		     
	  	      
		     
		     size_t newStateNr= Position(*it2);
		     
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= ValType{var};}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= ValType{var};}
		    }

	   	     }
	   

   	  }
    
  return op;
  }
   template<class TotalBasis>
  Mat CurrOperator(const TotalBasis& totalBasis, double var=1. , bool PB=0, int start=0, int stop=0)
     {
	 using TpBasisIt= typename TotalBasis::BasisIt;
using Lattice=typename TotalBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      stop = (stop!=0) ? stop : Operators::Length( sites, PB);       
      std::cout<< "stop "<< stop<<std::endl;
   for( auto& tpState : totalBasis)
	     {

	       	    for(size_t i=start; i<stop; i++)
                 {

	       

		   Lattice state=GetLattice(tpState);
                     size_t j=Operators::NextWithBC(i, sites, PB);
	   	    		    if(state[i]==state[j])
   	       	     {
	  // 	     
    		  
   	   	     }
   	   	    else{
		      

		      Lattice temp=state;
		      
		       state.setPartNr(j, temp[i]);
		      state.setPartNr(i, temp[j]);
		     		 
	  	     size_t signControl=CheckSign(temp, i, j);
		       
		    auto it2 = totalBasis.find(state.GetId());
		     
		    double otherSign=(state[j]==1)?+1:-1;

		      
		      
		     
		     size_t newStateNr= Position(*it2);
		     
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= otherSign*ValType{var};}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= otherSign*ValType{var};}
		    }

	   	     }
	   

   	  }
    
  return op;
  }
    template<typename TotalBasis, typename State, typename Lattice>
  void Act(int i, int j, TotalBasis& totalBasis, State& tpState, Lattice state, Mat& op, double var)
  {
    if(state[i]==state[j])
      {	
    		  
      }
    else{
      Lattice temp=state;
      
      state.setPartNr(j, temp[i]);
      state.setPartNr(i, temp[j]);
		     		 
      size_t signControl=CheckSign(temp, i, j);
      
		    auto it2 = totalBasis.find(state.GetId());
		    
	  	    
		     
		     size_t newStateNr= Position(*it2);
		     
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= ValType{var};}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= ValType{var};}
		    }
  }
    template<class TotalBasis>
    Mat totalHetOperator(const TotalBasis& totalBasis, double tint, double t0, double tl, double V, int Llead1,  int Llead2,  int Lchain )
     {
       bool PB=0;
	 using TpBasisIt= typename TotalBasis::BasisIt;
using Lattice=typename TotalBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
 
   for( auto& tpState : totalBasis)
	     {
const Lattice state=GetLattice(tpState);
	       	    for(size_t i=0; i<Llead1-1; i++)
                 {
	
		   if(state[i]==1)
		     {
		       op.coeffRef(Position(tpState), Position(tpState))-= ValType{V};
		     }
		   size_t j=Operators::NextWithBC(i, sites, PB);
		     // first lead

		     Act(i, j, totalBasis, tpState, state, op, tl);
		 }
		     for(size_t i=Llead1+Lchain; i<Llead1+Lchain+Llead2-1; i++)
                 {
		   if(state[i]==1)
		     {
		       op.coeffRef(Position(tpState), Position(tpState))+= ValType{V};
		     }
		   		  
				   size_t j=Operators::NextWithBC(i, sites, PB);
		     // second lead

		     Act(i, j, totalBasis, tpState, state, op, tl);

	   	     }
		    	   if(state[Llead1+Llead2+Lchain-1]==1)
		     {
		       op.coeffRef(Position(tpState), Position(tpState))+= ValType{V};
		     }
			   if(state[Llead1-1]==1)
		     {
		       op.coeffRef(Position(tpState), Position(tpState))-= ValType{V};
		     }
		    // chain
		    	       	    for(size_t i=Llead1; i<Llead2+Llead1+Lchain-1; i++)
                 {
		   Lattice state=GetLattice(tpState);
		  
                     size_t j=Operators::NextWithBC(i, sites, PB);
		     // first lead
		     Act(i, j, totalBasis, tpState, state, op, t0);
	   	     }
				    Act(Llead1-1, Llead1, totalBasis, tpState, state, op, tint);
				    Act(Llead1+Lchain-1, Llead1+Lchain, totalBasis, tpState, state, op, tint);	  }
    
  return op;
  }
template<class TotalBasis>
Mat totCurrOperator(const TotalBasis& totalBasis, double var ,  int Llead1, int Llead2, int Lchain)
     {
       using TpBasisIt= typename TotalBasis::BasisIt;
       using Lattice=typename TotalBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
   size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
   for( auto& tpState : totalBasis)
	     {

	       int i=Llead1-1;
	       int j=i+1;
	       

		   const Lattice state=GetLattice(tpState);
          
	   	    		    if(state[i]==state[j])
   	       	     {

    		  
   	   	     }
   	   	    else{
		      

		      Lattice temp=state;

		       temp.setPartNr(j, state[i]);
		      temp.setPartNr(i, state[j]);

	  	     size_t signControl=CheckSign(state, i, j);
		       
		    auto it2 = totalBasis.find(temp.GetId());
		     
		    double otherSign=(temp[j]==1)?+1:-1;

		      
		      
		     
		     size_t newStateNr= Position(*it2);
		     
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= otherSign*ValType{var};}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= otherSign*ValType{var};}
		    }
	 i=Lchain+Llead1-1;
	        j=i+1;

		if(state[i]==state[j])
   	       	     {
	  // 	     
    		  
   	   	     }
   	   	    else{
		      

		      Lattice temp=state;
		      
		       temp.setPartNr(j, state[i]);
		      temp.setPartNr(i, state[j]);
		     		 
	  	     size_t signControl=CheckSign(state, i, j);
		       
		    auto it2 = totalBasis.find(temp.GetId());
		     
		    double otherSign=(temp[j]==1)?+1:-1;

		      
		      
		     
		     size_t newStateNr= Position(*it2);
		     
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= otherSign*ValType{var};}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= otherSign*ValType{var};}
		    }

	   	     
	   

   	  }
    
  return op;
  }
}
