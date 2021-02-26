
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

   template<typename TotalBasis, typename State, typename Lattice>
  void Act(int i, int j, TotalBasis& totalBasis, State& tpState, Lattice state, Mat& op, double var)
  {
    if(state[i]==state[j])
      {	
    		  
      }
    else{
      Lattice temp=state;
      //  std::cout<< "i,jx "<< i << ", "<< j << std::endl;
       state.flip(i);
       //setPartNr(j, temp[i]);
      state.flip(j);
		     		 
      size_t signControl=CheckSign(temp, i, j);
    
		    auto it2 = totalBasis.find(state.GetId());
		    
	  	    
		     
		     size_t newStateNr= Position(*it2);
		     //	         std::cout<< "old "<< temp << " at "<< Position(tpState)<< std::endl;
      
		     // std::cout<< "new "<< state <<" at  "<<newStateNr<<  std::endl;
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= ValType{var};}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= ValType{var};}
		    }
  }


  // carefull with changes, making one for cdag and one for ccdag
    template<typename T>
    size_t CheckSign2(const T& state, size_t i, size_t j)
  {
    size_t m=0;
    if(i>j)
      {
  	for(size_t l=j; l<i; l++)
  	  {
  	    m+=size_t(state[l]);
  	  }

      }
    else
      	for(size_t l=i; l<j; l++)
  	  {
  	    m+=size_t(state[l]);
  	  }
    return m;
  }
  
  size_t NextWithBC(size_t i, size_t sites, bool PB=true, size_t steps=1)
  {
    if(PB)
      {return (i+steps)%sites;}
    else{return (i+steps);}
  }
  size_t Length( size_t sites, bool PB=true)
  {
    if(PB)
      {return sites;}
    else{return sites-1;}
  }
    template<class TotalBasis>
  Mat CdagOperator(const TotalBasis& totalBasis, int i,   const bool& PB=true)
{
  
   using BasisIt= typename TotalBasis::BasisIt;     



 using Lattice=typename TotalBasis::Lattice;   
    size_t dim=totalBasis.dim;
    
    Mat op(dim, dim);
       op.setZero();
 for( auto& tpState :totalBasis)
	     {
	  	 
	    
	       BasisIt it2=totalBasis.find(Id(tpState));	     

	       Lattice state=GetLattice(*it2);
	       //	       std::cout << " stat "<< totalBasis.totalmaxPar << "\n";
	       //std::cout << " i "<< i<< "  has "<< state[i]<< std::endl;
	       
	       if(state[i]==0 )
	         	 {
			   
			   auto temp=state;
	         	   state.setPartNr(i, 1);
			   size_t signControl=CheckSign2(temp, i, 0);
			   // std::cout << "at i="<< i <<" changed to  "<< state << "\n"<< "with aign cont "<< signControl<< '\n';
	         	   it2= totalBasis.find(state.GetId());
	        	   
	     		   size_t newStateNr= Position(*it2);
	      		              if(signControl%2==0)
   	      	    	 {
	      		   op.coeffRef(newStateNr, Position(tpState))+= ValType{1};}
   	      	       else
   	      	    	 { op.coeffRef(newStateNr, Position(tpState))-= ValType{1};}
	      	    }
	         
   	         	     
		 
	      }
      return op;
      
       }
  template<class TotalBasis>
  Mat COperator(const TotalBasis& totalBasis, int i,   const bool& PB=true)
{
  
   using BasisIt= typename TotalBasis::BasisIt;     



 using Lattice=typename TotalBasis::Lattice;   
    size_t dim=totalBasis.dim;
    
    Mat op(dim, dim);
       op.setZero();
 for( auto& tpState :totalBasis)
	     {
	  	 
	    
	        BasisIt it2=totalBasis.find(Id(tpState));	     

	       Lattice state=GetLattice(*it2);
	       // std::cout << " stat "<< state << "\n";
	       //std::cout << " i "<< i<< "  has "<< state[i]<< std::endl;
	       
	       if(state[i]==1)
	         	 {
			   
			   auto temp=state;
	         	   state.setPartNr(i, 0);
			   size_t signControl=CheckSign2(temp, i, 0);
			   // std::cout << "at i="<< i <<" changed to  "<< state << "\n"<< "with aign cont "<< signControl<< '\n';
	         	   it2= totalBasis.find(state.GetId());

	     		   size_t newStateNr= Position(*it2);
	      		              if(signControl%2==0)
   	      	    	 {
	      		   op.coeffRef(newStateNr, Position(tpState))+= ValType{1};}
   	      	       else
   	      	    	 { op.coeffRef(newStateNr, Position(tpState))-= ValType{1};}
	      	    }
	         
   	         	     
		 
	      }
      return op;
      
       }
  template<class TotalBasis>
  Mat NumberOperatorE(const TotalBasis& totalBasis, const double omega=1., bool PB=0, int start=0, int stop=0)
 {
   
 
    size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
 stop = (stop!=0) ? stop : sites;    


      for(const auto& tpState : totalBasis)
	{

	 
	  for(size_t i=start; i<stop; i++)
	    {
	      
	      
  op.coeffRef(Position(tpState), Position(tpState))+=ValType(omega*totalBasis.particlesAt(Id(tpState), i));
       
	   
	  }
	}
      return op;
      
       }
  template<class TotalBasis, class Functor>
  Mat EKinLongRangeOperator(const TotalBasis& totalBasis,Functor f, int lmx=1, double var=1. , bool PB=0 )
     {

using Lattice=typename TotalBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
    int  stop = Operators::Length( sites, PB);       

   for( auto& tpState : totalBasis)
	     {
	       	    for(int i=0; i<stop; i++)
                 {                     
		   for(int j=1; j<=lmx; j++)
                 {
Lattice state_1=GetLattice(tpState);
		   Lattice state_2=GetLattice(tpState);
		   int j_1=(i+j)%sites;


	   	    		    if(state_1[i]==state_1[j_1])
   	       	     {
	  // 	     
    		  
   	   	     }
   	   	    else{

		      Lattice temp=state_1;
		      state_1.switchPartNr(i, temp[j_1], j_1, temp[i]);

  size_t signControl=CheckSign(temp, i, j_1);


		    auto it2_1 = totalBasis.find(state_1.GetId());
	     		   size_t newStateNr_1= Position(*it2_1);
		
			   		
	   	               if(signControl%2==0)
   	   	    	  {
			    op.coeffRef(newStateNr_1, Position(tpState))-= ValType{var}*f(i,j_1);}

		   	       	   	       else
		   	    {
			      op.coeffRef(newStateNr_1, Position(tpState))+= ValType{var}*f(i,j_1);
}
		   

	   	     
		    }    		  


		 }	   

   	  }
	     }    
  return op;
  }
  template<class TotalBasis>
  Mat EKinOperator(const TotalBasis& totalBasis, double var=1. , bool PB=0, int start=0, int stop=0 )
     {

using Lattice=typename TotalBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      stop = (stop!=0) ? stop : Operators::Length( sites, PB);       

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



 
Mat EKinOperator(const Many_Body::OneElectronBasis& totalBasis, double var=1. , bool PB=0, int start=0, int stop=0 )
     {

       using Lattice=typename Many_Body::OneElectronBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      stop = (stop!=0) ? stop : Operators::Length( sites, PB);       

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

		      

   		    Act(i, j, totalBasis, tpState, state, op, var);
   		    }

   	   	     }
	   

   	  }
    
  return op;
  }

   template<class TotalBasis>
  Mat CurrOperator(const TotalBasis& totalBasis, double var=1. , bool PB=0, int start=0, int stop=0)
     {

using Lattice=typename TotalBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      stop = (stop!=0) ? stop : Operators::Length( sites, PB);       

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




    
  Mat CurrOperator(const  Many_Body::OneElectronBasis& totalBasis, double var=1. , bool PB=0, int start=0, int stop=0)
     {

using Lattice=typename  Many_Body::OneElectronBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
      stop = (stop!=0) ? stop : Operators::Length( sites, PB);       

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

		          state.flip(i);
      
			  state.flip(j);
	  	     size_t signControl=CheckSign(temp, i, j);
		       
		    auto it2 = totalBasis.find(state.GetId());
		     
		    double otherSign=(temp[j]==1)?+1:-1;

		      
		      
		     
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

    template<typename TotalBasis, typename State, typename Lattice>
  void Act(int i, int j, TotalBasis& totalBasis, State& tpState, Lattice state, Mat& op, double var)
  {
    if(state[i]==state[j])
      {	
    		  
      }
    else{
      Lattice temp=state;
      //  std::cout<< "i,jx "<< i << ", "<< j << std::endl;
       state.flip(i);
       //setPartNr(j, temp[i]);
      state.flip(j);
		     		 
      size_t signControl=CheckSign(temp, i, j);
    
		    auto it2 = totalBasis.find(state.GetId());
		    
	  	    
		     
		     size_t newStateNr= Position(*it2);
		     //	         std::cout<< "old "<< temp << " at "<< Position(tpState)<< std::endl;
      
		     // std::cout<< "new "<< state <<" at  "<<newStateNr<<  std::endl;
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

using Lattice=typename TotalBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
    size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
 int L_x=Llead1+(Lchain+1)/2;
 double E=0;
 if(Lchain>1)
   {

   E=V/(Lchain+1);

   }
 std::cout<< "L x "<<L_x<< " E "<<E<<std::endl;  
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
		   //	   std::cout<< "i,j "<< i << ", "<< j << std::endl;
		     // first lead


		     Act(i, j, totalBasis, tpState, state, op, tl);

		 }
		     for(size_t i=Llead1+Lchain; i<Llead1+Lchain+Llead2-1; i++)
                 {
		   if(state[i]==1)
		     {
		       // std::cout<<" i "<<i+1<< " V "<< V <<std::endl;

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
		    	       	    for(size_t i=Llead1; i<Llead1+Lchain-1; i++)
                 {
		   Lattice state=GetLattice(tpState);
		  
                     size_t j=Operators::NextWithBC(i, sites, PB);
		     // first lead
		     Act(i, j, totalBasis, tpState, state, op, t0);
		     	   if(state[i]==1)
		     {

		       		         op.coeffRef(Position(tpState), Position(tpState))+= (static_cast<double>((i+1))-L_x)*E;

		     }
	   	     }
				      	   if(state[Llead1+Lchain-1]==1)
		     {

		       	       op.coeffRef(Position(tpState), Position(tpState))+=(static_cast<double>(Llead1+Lchain)-L_x)*E;
				       //  std::cout<<" XXi "<<Llead1+Lchain << " V "<< (static_cast<ValType>(Llead1+Lchain)-L_x)*E<<std::endl;

		     }
				    Act(Llead1-1, Llead1, totalBasis, tpState, state, op, tint);
				    Act(Llead1+Lchain-1, Llead1+Lchain, totalBasis, tpState, state, op, tint);	  }
    
  return op;
  }
template<class TotalBasis>
Mat totCurrOperator(const TotalBasis& totalBasis, double var ,  int Llead1, int Llead2, int Lchain)
     {

       using Lattice=typename TotalBasis::Lattice;     

 
   size_t dim=totalBasis.dim;
   size_t sites=totalBasis.sites;
    Mat op(dim, dim);
      op.setZero();
   for( auto& tpState : totalBasis)
	     {


	       int i=Llead1-1;

	       int j=i+1;
	       

		   Lattice state=GetLattice(tpState);
          
	   	    		    if(state[i]==state[j])
   	       	     {

    		  
   	   	     }
   	   	    else{
		      

		      Lattice temp=state;

		          state.flip(i);
      
			  state.flip(j);
	  	     size_t signControl=CheckSign(temp, i, j);
		       
		    auto it2 = totalBasis.find(state.GetId());
		     
		    double otherSign=(temp[j]==1)?+1:-1;

		      
		      
		     
		     size_t newStateNr= Position(*it2);
		     
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= otherSign*ValType{var};}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= otherSign*ValType{var};}
		    }
				    state=GetLattice(tpState);

	 i=Lchain+Llead1-1;


	        j=i+1;

		if(state[i]==state[j])
   	       	     {
	  // 	     
    		  
   	   	     }
   	   	    else{
		      

		      Lattice temp=state;
		    state.flip(i);
       //setPartNr(j, temp[i]);
      state.flip(j);
		     		 
	  	     size_t signControl=CheckSign(temp, i, j);
		       
		    auto it2 = totalBasis.find(state.GetId());
		     
		    double otherSign=(temp[j]==1)?+1:-1;

		      
		      
		     
		     size_t newStateNr= Position(*it2);
		     
	   	              if(signControl%2==0)
   	   	    	 {
	   		   op.coeffRef(newStateNr, Position(tpState))-= otherSign*ValType{var};}
   	   	       else
   	   	    	 { op.coeffRef(newStateNr, Position(tpState))+= otherSign*ValType{var};}
		    }

	   	     
	   

   	  }

   op*=0.5;   

  return op;
  }
}
