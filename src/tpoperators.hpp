# pragma once
#define EIGEN_USE_MKL_ALL
#include"basis.hpp"
#include<cmath>
#include <Eigen/Sparse>
#include <type_traits>
#include"operators.hpp"
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
  Mat NumberOperator(const TotalBasis& totalBasis, const SubBasis& subBasis, const double omega=1., const bool& PB=true)

  {


    size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
   op.setZero();
  


      for(const auto& tpState : totalBasis)
	{

	 
	  for(size_t i=0; i<sites; i++)
	    {

 op.coeffRef(Position(tpState), Position(tpState))+=(omega*subBasis.particlesAt(RightId(tpState), i));
       
	   
	  }
	}
      return op;
      
       }

  template<class TotalBasis, class SubBasis >
  Mat EKinOperatorL(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , const bool& PB=true)
     {


	      using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
     using RightBasisIt= typename TotalBasis::RightBasisIt;


using Lattice=typename SubBasis::Lattice;     
    
 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      
      
      
	   for( auto& tpState : totalBasis)
	     {

	       	    for(size_t i=0; i<Operators::Length( sites, PB); i++)
                 {

	       LeftBasisIt it2=subBasis.find(LeftId(tpState));	     

	       Lattice state=GetLattice(*it2);
                     size_t j=Operators::NextWithBC(i, sites, PB);
		     //		   std::cout<< i << "  " << j <<std::endl;
	   	    		    if(state[i]==state[j])
   	       	     {
	  // 	     
    		  
   	   	     }
   	   	    else{
		      

		      Lattice temp=state;
		      
		       state.setPartNr(j, temp[i]);
		      state.setPartNr(i, temp[j]);
		     		 
	  	     size_t signControl=CheckSign(temp, i, j);
		       
		     it2= subBasis.find(state.GetId());
		     
	  	       RightBasisIt it3=totalBasis.rbasis.find(RightId(tpState));
		     
	  	       size_t newStateNr= Position(*it3)*totalBasis.lbasis.dim +Position(*it2);
		     
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
 template<class TotalBasis, class SubBasis >
  Mat EKinOperatorLNNN(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , size_t m=0)
     {


	      using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
     using RightBasisIt= typename TotalBasis::RightBasisIt;


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
				       state.setPartNr(j, temp[i]);
				       state.setPartNr(i, temp[j]);
		 
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


       using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
       using RightBasisIt= typename TotalBasis::RightBasisIt;


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
			     state.setPartNr(j, temp[i]);
 state.setPartNr(i, temp[j]);
		 
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


       using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
       using RightBasisIt= typename TotalBasis::RightBasisIt;


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
		      state.setPartNr(j, temp[i]);
 state.setPartNr(i, temp[j]);
		 
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
  // Holstein model obc, for now this is n_0(n)
   template<typename TotalBasis, typename Basis>
  Mat NBosonCOperator(const TotalBasis& totalBasis, Basis& subBasis, double var=1. , const bool& PB=true)
  {


      using LeftBasisIt= typename TotalBasis::LeftBasisIt;
      using LeftLattice=  typename    TotalBasis::LeftLattice;     
      using RightBasisIt= typename TotalBasis::RightBasisIt;     
      using RightLattice=typename TotalBasis::RightLattice;

   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
        	for(const auto& tpState : totalBasis)   
  	  {	    

  	 
  	 
	    for(size_t i=0; i<sites; i++)
	      {

		  	    RightBasisIt it2=subBasis.find(RightId(tpState));

	    	    LeftBasisIt it33=totalBasis.lbasis.find(LeftId(tpState));
	     LeftLattice stateEl(GetLattice(*it33));

	     if(stateEl[i]==0){

	       continue;}
  	    size_t particleNumber=subBasis.particlesAt(RightId(tpState), i);

  	     if( particleNumber==subBasis.maxParticles)
  	       {
  	       }
  	     else{
	       
  	       LeftBasisIt it3=totalBasis.lbasis.find(LeftId(tpState));
  	       RightLattice state(GetLattice(*it2));


	       //	       std::cout<<" at site " << i << " we  started with e "<< stateEl<< "and ph  started with satte " << state << std::endl;
			     state.setPartNr(i, particleNumber+1);

			     //			      std::cout<< "and got " << state << std::endl;

  	       it2= subBasis.find(state.GetId());
  	      
  	       size_t newStateNr= Position(*it2) *totalBasis.lbasis.dim +Position(*it3) ;



  	       op.coeffRef(newStateNr, Position(tpState))+= ValType{var*(std::sqrt(state[i]))};
  	   

  		    }
	      }
  	  }
  		return op;
  }

 template<typename TotalBasis, typename Basis>
 Mat NBosonDOperator(const TotalBasis& totalBasis, Basis& subBasis,   double var=1. ,  const bool& PB=true)
  {
     



      using LeftBasisIt= typename TotalBasis::LeftBasisIt;
      using LeftLattice=  typename    TotalBasis::LeftLattice;     
      using RightBasisIt= typename TotalBasis::RightBasisIt;     
      using RightLattice=typename TotalBasis::RightLattice;
       
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
  	for(const auto& tpState : totalBasis)   
  	  {

  	    
	    for(size_t i=0; i<sites; i++)
	      {
		RightBasisIt it2=subBasis.find(RightId(tpState));
		//		std::cout<< "at "<<Position(tpState) <<" havinfg state BD " << i <<std::endl;

	    	    LeftBasisIt it33=totalBasis.lbasis.find(LeftId(tpState));
	     LeftLattice stateEl(GetLattice(*it33));
	     //std::cout<< "acted on e state " << stateEl<< std::endl; 
	     if(stateEl[i]==0){
	       //	       std::cout<< " evaluateing "<<stateEl << "at " << i << " gave " <<stateEl[i] <<std::endl;
										 continue;}
	     //std::cout<< "and got " << satte << std::endl;
	    size_t particleNumber=subBasis.particlesAt(RightId(tpState), i);

  	     if( particleNumber==0)
  	       {
  	       }
  	     else{
	       
  	       LeftBasisIt it3=totalBasis.lbasis.find(LeftId(tpState));
  	       RightLattice state(GetLattice(*it2));
	       //std::cout<<" at site " << i << " we  started with e "<< stateEl<< "and ph  started with satte " << state << std::endl;

	       state.setPartNr(i, particleNumber-1);
	       //      std::cout<< "and got " << state << std::endl;
  	       it2= subBasis.find(state.GetId());

  	       size_t newStateNr= Position(*it2)*totalBasis.lbasis.dim +Position(*it3);
	       

  	       op.coeffRef(newStateNr, Position(tpState))+= ValType{var*(std::sqrt(state[i]+1))};
  	   
	     
  		    }
  	  }
	  }
  		return op;
   }
  
   template<class TotalBasis, class SubBasis >
  Mat EKinOperatorR(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , const bool& PB=true)
     {


	      using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
     using RightBasisIt= typename TotalBasis::RightBasisIt;


using Lattice=typename SubBasis::Lattice;     
    
 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      
      
      
	   for( auto& tpState : totalBasis)
	     {

	       	    for(size_t i=0; i<Operators::Length( sites, PB); i++)
                 {

	       RightBasisIt it2=subBasis.find(RightId(tpState));	     
	       LeftBasisIt iL=totalBasis.lbasis.find(LeftId(tpState));
	       Lattice state=GetLattice(*it2);
                     size_t j=Operators::NextWithBC(i, sites, PB);
		     //		   std::cout<< i << "  " << j <<std::endl;
	   	    		    if(state[i]== 0 or state[j]== totalBasis.rbasis.maxParticles)
   	       	     {
	  // 	     
    		  
   	   	     }
    	   	    else{
		      

		      Lattice temp=state;
		      
		       state.setPartNr(j, temp[j]+1);
		      state.setPartNr(i, temp[i]-1);
		     		 
// 		      //	  	     size_t signControl=CheckSign(temp, i, j);
		       
		     it2= subBasis.find(state.GetId());
		     
	  	       RightBasisIt it3=totalBasis.rbasis.find(state.id);
		     
	  	       size_t newStateNr= Position(*it3)*totalBasis.lbasis.dim +Position(*iL);
		     

		       op.coeffRef(newStateNr, Position(tpState))-= ValType{var}*std::sqrt(temp[j]+1)*std::sqrt(temp[i]);
// 		
		    }
				    	     it2=subBasis.find(RightId(tpState));	     
 iL=totalBasis.lbasis.find(LeftId(tpState));
	 state=GetLattice(*it2);
				    	   	    		    if(state[j]== 0 or state[i]== totalBasis.rbasis.maxParticles)
   	       	     {
	  // 	     
    		  
   	   	     }
    	   	    else{
		      

		      Lattice temp=state;
		      
		       state.setPartNr(j, temp[j]-1);
		      state.setPartNr(i, temp[i]+1);
		     		 
// 		      //	  	     size_t signControl=CheckSign(temp, i, j);
		       
		     it2= subBasis.find(state.GetId());
		     
	  	       RightBasisIt it3=totalBasis.rbasis.find(state.id);
		     
	  	       size_t newStateNr= Position(*it3)*totalBasis.lbasis.dim +Position(*iL);
		     

		       op.coeffRef(newStateNr, Position(tpState))-= ValType{var}*std::sqrt(temp[j])*std::sqrt(temp[i]+1);
// 		
		    }
	   	     }
 	      	     }
	   

   	  
    
  return op;
  }
     template<class TotalBasis, class SubBasis >
  Mat EKinOperatorRNN(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , const bool& PB=true)
     {


	      using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
     using RightBasisIt= typename TotalBasis::RightBasisIt;


using Lattice=typename SubBasis::Lattice;     
    
 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      
      
      
	   for( auto& tpState : totalBasis)
	     {

	       	    for(size_t i=0; i<Operators::Length( sites, PB); i++)
                 {

	       RightBasisIt it2=subBasis.find(RightId(tpState));	     
	       LeftBasisIt iL=totalBasis.lbasis.find(LeftId(tpState));
	       Lattice state=GetLattice(*it2);
	       size_t j=(i+2)%sites;
	       //Operators::NextWithBC(i, sites, PB);
		     //		   std::cout<< i << "  " << j <<std::endl;
	   	    		    if(state[i]== 0 or state[j]== totalBasis.rbasis.maxParticles)
   	       	     {
	  // 	     
    		  
   	   	     }
    	   	    else{
		      

		      Lattice temp=state;
		      
		       state.setPartNr(j, temp[j]+1);
		      state.setPartNr(i, temp[i]-1);
		     		 
// 		      //	  	     size_t signControl=CheckSign(temp, i, j);
		       
		     it2= subBasis.find(state.GetId());
		     
	  	       RightBasisIt it3=totalBasis.rbasis.find(state.id);
		     
	  	       size_t newStateNr= Position(*it3)*totalBasis.lbasis.dim +Position(*iL);
		     

		       op.coeffRef(newStateNr, Position(tpState))-= ValType{var}*std::sqrt(temp[j]+1)*std::sqrt(temp[i]);
// 		
		    }
				    	     it2=subBasis.find(RightId(tpState));	     
 iL=totalBasis.lbasis.find(LeftId(tpState));
	 state=GetLattice(*it2);
				    	   	    		    if(state[j]== 0 or state[i]== totalBasis.rbasis.maxParticles)
   	       	     {
	  // 	     
    		  
   	   	     }
    	   	    else{
		      

		      Lattice temp=state;
		      
		       state.setPartNr(j, temp[j]-1);
		      state.setPartNr(i, temp[i]+1);
		     		 
// 		      //	  	     size_t signControl=CheckSign(temp, i, j);
		       
		     it2= subBasis.find(state.GetId());
		     
	  	       RightBasisIt it3=totalBasis.rbasis.find(state.id);
		     
	  	       size_t newStateNr= Position(*it3)*totalBasis.lbasis.dim +Position(*iL);
		     

		       op.coeffRef(newStateNr, Position(tpState))-= ValType{var}*std::sqrt(temp[j])*std::sqrt(temp[i]+1);
// 		
		    }
	   	     }
 	      	     }
	   

   	  
    
  return op;
  }
   template<class TotalBasis, class SubBasis >
   Mat eMOMOperator(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , int k=0, const bool& PB=true)
     {
  using Many_Body::pi;
  using Many_Body::im;

	      using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
     using RightBasisIt= typename TotalBasis::RightBasisIt;


using Lattice=typename SubBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
 auto K=2*pi/sites;      
      
      
	   for( auto& tpState : totalBasis)
	     {

	       	    for(size_t i=0; i<sites; i++)
                 {

	   
                       for(size_t j=0; j<sites; j++)
                 {
		       RightBasisIt it2=totalBasis.rbasis.find(RightId(tpState));	     
	       LeftBasisIt iL=totalBasis.lbasis.find(LeftId(tpState));
	       Lattice state=GetLattice(*iL);
	       if(j==i)
		 {
		   	      Lattice temp=state;
		     if(state[i]!= 0)
		       {
		       op.coeffRef(Position(tpState), Position(tpState))+= ValType{var}/sites;
		       }
		 }
	       else{
		     //		   std::cout<< i << "  " << j <<std::endl;
	   	    		    if(state[i]== 0)
   	       	     {
	  // 	     
    		  
   	   	     }
     	   	    else{
		      

		      Lattice temp=state;
		      
		       state.setPartNr(j, 1);
		      state.setPartNr(i, 0);
		     		 
// 		      //	  	     size_t signControl=CheckSign(temp, i, j);
		       
		      //iL= subBasis.find(state.GetId());
		     
	  	       LeftBasisIt it3=totalBasis.lbasis.find(state.id);
		     
	  	       size_t newStateNr= Position(*it2)*totalBasis.lbasis.dim +Position(*it3);
		     
		       std::complex<double> value=exp(-im*static_cast<std::complex<double>>(K*(j-i)));
		       op.coeffRef(newStateNr, Position(tpState))+= ValType{var}/sites;
		       //*std::sqrt(temp[j]+1)*std::sqrt(temp[i]);
// 		
 		    }
	       }	    	    
	 

 	   	     }
  	      	     }
	   
	     }
   	  
    
  return op;
  }

     template<class TotalBasis, class SubBasis >
   Mat phMOMOperator(const TotalBasis& totalBasis, const SubBasis& subBasis, double var=1. , int k=0, const bool& PB=true)
     {


	      using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
     using RightBasisIt= typename TotalBasis::RightBasisIt;


using Lattice=typename SubBasis::Lattice;     
    
 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      
      
      
	   for( auto& tpState : totalBasis)
	     {

	       	    for(size_t i=0; i<sites; i++)
                 {

	   
                       for(size_t j=0; j<sites; j++)
                 {
		       RightBasisIt it2=totalBasis.rbasis.find(RightId(tpState));	     
	       LeftBasisIt iL=totalBasis.lbasis.find(LeftId(tpState));
	       Lattice state=GetLattice(*it2);
	       if(j==i)
		 {
		   	      Lattice temp=state;
		     if(state[i]!= 0)
		       {
		       op.coeffRef(Position(tpState), Position(tpState))+=state[i]*ValType{var}/sites;
		       }
		 }
	       else{
		     //		   std::cout<< i << "  " << j <<std::endl;
	   	    		    if(state[i]== 0 or state[j]== totalBasis.rbasis.maxParticles)
   	       	     {
	  // 	     
    		  
   	   	     }
     	   	    else{
		      

		      Lattice temp=state;
		      
		       state.setPartNr(j, temp[j]+1);
		       state.setPartNr(i, temp[i]-1);
		     		 
// 		      //	  	     size_t signControl=CheckSign(temp, i, j);
		       
		      //iL= subBasis.find(state.GetId());
		     
	  	       RightBasisIt it3=totalBasis.rbasis.find(state.id);
		     
	  	       size_t newStateNr= Position(*it3)*totalBasis.lbasis.dim +Position(*iL);
		     

		       op.coeffRef(newStateNr, Position(tpState))+= std::sqrt(temp[j]+1)*std::sqrt(temp[i])*ValType{var}/sites;
		       //*std::sqrt(temp[j]+1)*std::sqrt(temp[i]);
// 		
 		    }
	       }	    	    
	 

 	   	     }
  	      	     }
	   
	     }
   	  
    
  return op;
  }

  // template<class TotalBasis, class SubBasis>
 //  Mat NumberOperatore(const TotalBasis& totalBasis, const SubBasis& subBasis, const double omega=1., const bool& PB=true)
 // {    using BasisIt= typename TotalBasis::BasisIt;
 //    using Const_BasisIt= typename TotalBasis::Const_BasisIt;
 //    size_t dim=totalBasis.dim;
 //    size_t sites=totalBasis.sites;
 //    Mat op(dim, dim);
 //   op.setZero();
  


 //      for(const auto& tpState : totalBasis)
 // 	{

	 
 // 	  for(size_t i=0; i<Operators::Length( sites, PB); i++)
 // 	    {

 // 	      op.coeffRef(Position(tpState), Position(tpState))+=ValType(omega*subBasis.particlesAt(LeftId(tpState), i));
       
	   
 // 	  }
 // 	}
      // return op;
      
      //  }
}
