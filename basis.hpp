#pragma once
/*
indexing: the site number is counted from left to right,
such that electronstat[0] returns 1 for the state 100                                                            
 
 */
#include<iostream>
#include<vector>
#include<bitset>
#include<array>
#include<map>
#include<fstream>
#include<iomanip>
#include<cassert>
#include<cmath>
#include <algorithm> 
#include<iterator>
#include<numeric>
#include<tuple>
// computing the factorial
size_t Factorial(const size_t n)
{
  if(n>20)
    {
      std::cout << "to large number " << std::endl;
}
   size_t i, x = 1;
   for (i = 1; i <= n ; i++)
   {
      x *= i;
   }
  return x;
  }
double PrimeNumber(const size_t n)
{
      size_t check,c=0;
    for(size_t i=2;i<=1000;i++)
      {
        check=0;

        for(size_t j=2;j<=i/2;j++)
        {
            if(i%j==0)
            {
              check=1;
               break;
            }
        }

      if(check==0)
        c++;

          if(c==n)
         {
          return i;
          break;
         }
     }
    return 0;
}

 template<size_t L>
 struct ElectronState
 {
   ElectronState(): state(0), sites(L)  {};
   ElectronState( std::bitset<L> state): state(state), sites(L)  {};
   ElectronState(size_t i): state(i), sites(L)  {
 	 

	 
	  assert(i<static_cast<size_t>(std::pow(2, L)));
   };
   std::bitset<L> state;

    size_t operator[] ( size_t i) {
      return static_cast<size_t>(state[L-1-i]);
    }
 inline double GetId()
   {
     return static_cast<double> (state.to_ulong());
       
       }
   friend std::ostream& operator<<(std::ostream& os,  ElectronState& state)
   {

     for(size_t i=0; i<L; i++)
       { os<< state[i] << std::setw(3);
   	   }
     os << '\n';
     return os;
   }
   size_t sites;
 };

  template<size_t L>
struct BosonState
  {
    BosonState(): sites(L)   {state.fill(0);};
   
    BosonState(const std::array<size_t, L>& state): state{state}, sites(L) {};


    std::array<size_t, L> state;

    inline double GetId()
    {
      double id=0;
      if(std::accumulate(state.begin(), state.end(), 0)==0)
	{}
      else{
       for(size_t k=0; k<L; k++)
        {
        // defining bijective hashkey to identify the states
        id+=std::sqrt(PrimeNumber(k+1))*state[k];	   
        }
      }
      return id;
    }
    
    size_t& operator[] ( size_t i) {
        return state[i];
    }


       friend std::ostream& operator<<(std::ostream& os,  BosonState& state)
   {

     for(size_t i=0; i<L; i++)
       { os<< state[i]  << std::setw(3);
   	   }
     os << '\n';
     return os;
   }
 const  size_t sites;
  };

template<size_t L>
struct ElectronBasis
{  
   using  StateIt=typename std::array<size_t, L>::iterator;
  using BasisPair= typename std::pair<size_t, ElectronState<L>>;
  using  BasisIt= typename  std::map<double, BasisPair>::iterator;
  ElectronBasis(size_t numberOfParticles): sites(L), dim(0) {
    for(size_t i=0; i<static_cast<size_t>(std::pow(2, L)); i++)
        {

	   
   	  ElectronState<L> newState(i);
            if(newState.state.count()==numberOfParticles)
             {
	      
   	      basis.insert({newState.GetId(), {dim, newState}});
	      dim++;    
             }
        
   }
  }
  ElectronBasis():   sites(L), dim(0){
     for(size_t i=0; i<static_cast<size_t>(std::pow(2, L)); i++)
        {
   	  ElectronState<L> newState(i);
            
             
   	     
   	      basis.insert({newState.GetId(), {i, newState}});
	       dim++;    
             }
        
   }
  ElectronBasis(ElectronState<L> state): sites(L), dim(0){
     
       
   	  
            
             
   	     
   	      basis.insert({state.GetId(), {0, state}});
	       dim++;    
     
        
   }
  ElectronBasis(const ElectronBasis& e)= default;
  friend std::ostream& operator<<(std::ostream& os,  ElectronBasis& electronBasis)
   {
     assert(electronBasis.dim==electronBasis.basis.size());
     BasisIt it1= electronBasis.basis.begin();
   
     while( it1!=electronBasis.basis.end())
       {
	 os<< "id"<< '\t' << it1->first << '\t' <<it1->second.first <<  '\t' <<it1->second.second ;
     it1++;
	  }

     return os;
   }
   
  std::map<double, BasisPair> basis;
  
  size_t sites;
   size_t dim;
  //  const size_t numberofparticles;
  
};

template<size_t L>
  struct PhononBasis
{
  using StateIt= typename std::array<size_t, L>::iterator;
  using BasisPair= typename std::pair<size_t, BosonState<L>>;
  using BasisIt= typename  std::map<double, BasisPair>::iterator;
  using Type= BosonState<L>;
  PhononBasis(size_t maxPhonons);
  PhononBasis()=default;
  std::map<double, BasisPair> basis;
  size_t dim;
  const size_t sites;
  const size_t maxPhonons; // maximum number of phonons on one site
  // bosonstate<L>& operator[](double i)
  //  {
  //    return basis[i];
  //  }
  friend std::ostream& operator<<(std::ostream& os,  PhononBasis& phononBasis)
   {
     assert(phononBasis.dim==phononBasis.basis.size());
     BasisIt it1= phononBasis.basis.begin();
   
     while( it1!=phononBasis.basis.end())
       {
	 os<< "id"<< '\t' << it1->first << '\t' <<it1->second.first<< '\t' <<it1->second.second;
     it1++;
	  }

     return os;
   }

   
  };

template<size_t L>
  struct BosonBasis
{
  using StateIt=typename std::array<size_t, L>::iterator;
  using Type= BosonState<L>;
  using BasisPair= typename std::pair<size_t, BosonState<L>>;
  using BasisIt= typename  std::map<double, BasisPair>::iterator;
  
  BosonBasis(size_t maxBosons);
  std::map<double, BasisPair> basis;
  const size_t numberOfParticles; // maximum number of bosons on one site
  size_t dim;
  
  
  BosonState<L>& operator[](double i)
   {
     return basis[i];
   }
  friend std::ostream& operator<<(std::ostream& os,  BosonBasis& bosonBasis)
   {
     assert(bosonBasis.dim==bosonBasis.basis.size());
     BasisIt it1= bosonBasis.basis.begin();
   
     while( it1!=bosonBasis.basis.end())
       {
	 os<< "id"<< '\t' << it1->first << '\t' <<it1->second.second<< '\t' <<it1->second.second;
     it1++;
	  }

     return os;
   }

   
  };


 template<class LeftBasis, class RightBasis>
 struct TensorProduct{

   struct BasisIds
   {
     BasisIds(double leftId, double rightId): leftId(leftId), rightId(rightId) {}
     double leftId;
     double rightId;
   };
        using  BasisIt= typename std::map<size_t, BasisIds>::iterator;
        using  LeftBasisType=  LeftBasis;
        using RightBasisType= RightBasis;
   
 
  
   TensorProduct( LeftBasis&  lbasis,   RightBasis& rbasis ) : sites(lbasis.sites), dim(lbasis.dim*rbasis.dim), dimfirst(lbasis.dim), dimsecond(rbasis.dim), lbasis(lbasis), rbasis(rbasis)  {
          auto it1= lbasis.basis.begin();

       		size_t i=0;
       		  size_t j=0;
		 
       		while(it1!=lbasis.basis.end())
       		  {
       		    auto  it2= rbasis.basis.begin();
      		    
 		    while(it2!=rbasis.basis.end()){
		      //	 BasisIds 
 		
 	     basis.insert({it2->second.first*lbasis.dim+ it1->second.first, {it1->second.second.GetId(), it2->second.second.GetId()}});
 				   ++it2;
				
 				 	j++;
       		  }
 		 ++it1;
 		 i++;
 		  }
 		// 		std::cout << basis.size() << '\n';
 		assert(dim==basis.size());
		
       }
   const size_t sites;    
   size_t dim;
   size_t dimfirst;
   size_t dimsecond;
   LeftBasis lbasis;
   RightBasis rbasis;
   std::map<size_t, BasisIds > basis;
     friend std::ostream& operator<<(std::ostream& os,  TensorProduct& tensorProduct)
    {
      assert(tensorProduct.dim==tensorProduct.basis.size());
      BasisIt it1= tensorProduct.basis.begin();
   
      while( it1!=tensorProduct.basis.end())
        {
      	 os<< "id"<< '\t' << it1->first << '\t' <<it1->second.leftId<< '\t' <<it1->second.rightId << '\n';
      it1++;
      	  }

      return os;
    }


 };

  // MEMBER FUNCTION:
template<size_t L>
PhononBasis<L>::PhononBasis(size_t maxPhonons): dim{0},  sites(L), maxPhonons(maxPhonons)  {
  std::array<size_t, L> stateArray;
  stateArray.fill(0);   
  basis.insert({0., {0, stateArray}});
     size_t numberOfPhonons=0;
      dim++;
      // generating all states of the form 0000...0phononnumber
    while(numberOfPhonons<maxPhonons)
    {numberOfPhonons++;
     stateArray.fill(0);
     stateArray[L-1]= numberOfPhonons;		 
     do{  
       BosonState<L> newState{stateArray};
       basis.insert({newState.GetId(), {dim, newState}});
     dim++;

       }while(std::next_permutation(stateArray.begin(), stateArray.end()));
     
      
   		  // Permutating over the rest of the states
  
           StateIt it1= stateArray.begin();
        while(stateArray[0]<stateArray[L-1])
        {
         if(*it1==*(it1+1))                    // if the number of phonons at size i ar
                                    //equal to the number on site i+1 it set it to
                                    // zero and iterate one forward
        {*it1=0;
    	   it1++; }
              else                                      // if the number is smaller, i increase it
    		                                       // by one and go back to site 0
        { *it1+=1; 
         it1=stateArray.begin();
    do{
       BosonState<L> newState{stateArray};
       basis.insert({newState.GetId(), {dim, newState}});
           dim++;
              }while(std::next_permutation(stateArray.begin(), stateArray.end()));
    }   
      }
 }
      
   
}

template<size_t L>
BosonBasis<L>::BosonBasis(size_t numberOfParticles): numberOfParticles(numberOfParticles), dim{0} {
  std::array<size_t, L> statearray;
  statearray.fill(0);   

     size_t numberofbosons=0;

      // generating all states of the form 0000...0bosonnumber
    while(numberofbosons<numberOfParticles)
    {numberofbosons++;
     statearray.fill(0);
     statearray[L-1]= numberofbosons;		 
     do{
       if(static_cast<size_t>(std::accumulate(statearray.begin(), statearray.end(), 0))==numberOfParticles) {
       BosonState<L> newstate{statearray};
       basis.insert({newstate.GetId(),{dim, newstate}});
     dim++;
	}

       }while(std::next_permutation(statearray.begin(), statearray.end()));    
     
    		  // Permutating over the rest of the states

           StateIt it1= statearray.begin();
        while(statearray[0]<statearray[L-1])
        {
         if(*it1==*(it1+1))                    // if the number of bosons at size i ar
                                    //equal to the number on site i+1 it set it to
                                    // zero and iterate one forward
        {*it1=0;
    	   it1++; }
              else                                      // if the number is smaller, i increase it
    		                                       // by one and go back to site 0
        { *it1+=1; 
         it1=statearray.begin();
    do{
      if(static_cast<size_t>(std::accumulate(statearray.begin(), statearray.end(), 0))==numberOfParticles) {
      BosonState<L> newstate{statearray};
      basis.insert({newstate.GetId(), {dim, newstate}});
       
           dim++;
       }
              }while(std::next_permutation(statearray.begin(), statearray.end()));
    }   
      }
 }
   
}




