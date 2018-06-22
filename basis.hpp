#pragma once
/*
indexing: the site number is counted from left to right,
such that electronstat[0] returns 1 for the state 100                                                            
 
 */
#include<iostream>
#include<vector>
#include<bitset>
#include<array>
#include<fstream>
#include<iomanip>
#include<cassert>
#include<cmath>
#include <algorithm> 
#include<iterator>
#include<numeric>
#include<tuple>
#include<set>
#include<type_traits>
#include"accesfunctions.hpp"
// computing the factorial
template<typename State>
std::tuple<double, size_t, State> getState(double id)
{
  State aState(0);
  return {id, 0, aState};
}
template<typename State>
struct CompareState{
     bool operator()(const State& lhs, const State& rhs)
     {
       return std::get<0>(lhs)< std::get<0>(rhs);
     }
     
   };
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
    using LatticeIt=typename std::array<size_t, L>::iterator;
   ElectronState(): state(0), sites(L)  {};
   ElectronState( std::bitset<L> state): state(state), sites(L)  {};
   ElectronState(size_t i): state(i), sites(L)  {
 	 

	 
	  assert(i<static_cast<size_t>(std::pow(2, L)));
   };
   std::bitset<L> state;

const size_t operator[] ( size_t i) const{
      return static_cast<size_t>(state[L-1-i]);
    }
 inline double GetId()
   {
     return static_cast<double> (state.to_ulong());
       
       }
   friend std::ostream& operator<<(std::ostream& os, const ElectronState& state) 
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
    using LatticeIt=typename std::array<size_t, L>::iterator;
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
      const  size_t& operator[] ( size_t i) const
    {
        return state[i];
    }


       friend std::ostream& operator<<(std::ostream& os, const  BosonState& state)
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
  using Lattice=ElectronState<L>;
  using BasisType= std::tuple< double, size_t, Lattice>;
  using LatticeIt=typename Lattice::LatticeIt;
  using Basis=std::set<BasisType, CompareState<BasisType>>;
  using BasisIt= typename  Basis::iterator;
 

   
  ElectronBasis(size_t numberOfParticles): sites(L), dim(0) {
    for(size_t i=0; i<static_cast<size_t>(std::pow(2, L)); i++)
        {

	   
   	  ElectronState<L> newState(i);
            if(newState.state.count()==numberOfParticles)
             {
	      
   	      basis.insert({newState.GetId(), dim, newState});
	      dim++;    
             }
        
   }
  }
  ElectronBasis():   sites(L), dim(0){
     for(size_t i=0; i<static_cast<size_t>(std::pow(2, L)); i++)
        {
   	  ElectronState<L> newState(i);
            
             
   	     
	  basis.insert({newState.GetId(), i, newState});
	       dim++;    
             }
        
   }
  ElectronBasis(ElectronState<L> state): sites(L), dim(0){
     
       
   	  
            
             
   	     
   	      basis.insert({state.GetId(), 0, state});
	       dim++;    
     
        
   }
  ElectronBasis(const ElectronBasis& e)= default;
  friend std::ostream& operator<<(std::ostream& os,  const ElectronBasis& electronBasis)
   {
     assert(electronBasis.dim==electronBasis.basis.size());
     BasisIt it1= electronBasis.basis.begin();
   
     while( it1!=electronBasis.basis.end())
       {
	 os<< "id"<< '\t' << std::get<toBasisType(BasisInfoField::id)>(*it1) << '\t'<<"position"<< '\t'<<std::get<toBasisType(BasisInfoField::position)>(*it1)<<  '\t'<< "state"<< '\t' << std::get<toBasisType(BasisInfoField::state)>(*it1);
     it1++;
	  }

     return os;
   }
    const BasisType operator[] (double id) {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return *basis.find(aState);
  }
   size_t particlesAt(double id, size_t site)
  {
    Lattice aLattice;
    BasisType aState=*basis.find({id, 0, aLattice});
    return std::get<toBasisType(BasisInfoField::state)>(aState)[site];
    
  }
  BasisIt begin() {return basis.begin();}
  BasisIt end() {return basis.end();}
  BasisIt find(double id)
  {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return basis.find(aState);
  }
  Basis basis;
  
  size_t sites;
   size_t dim;
  //  const size_t numberofparticles;
  
};

template<size_t L>
  struct PhononBasis
{
  using Lattice= BosonState<L>;  
  using BasisType= std::tuple<double, size_t, Lattice>;
  using Basis=std::set<BasisType, CompareState<BasisType>>;
  using LatticeIt= typename Lattice::LatticeIt;
  using BasisIt= typename  Basis::iterator;
 
  PhononBasis(size_t maxPhonons);
  PhononBasis(): sites(L), maxPhonons(0) {};
  Basis basis;
  size_t dim;
  const size_t sites;
  const size_t maxPhonons; // maximum number of phonons on one site
  
  const BasisType operator[] (double id) {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return *basis.find(aState);
  }
   size_t particlesAt(double id, size_t site)
  {
    Lattice aLattice;
    BasisType aState=*basis.find({id, 0, aLattice});
    return std::get<toBasisType(BasisInfoField::state)>(aState)[site];
    
  }
  BasisIt find(double id)
  {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return basis.find(aState);
  }
  BasisIt begin() {return basis.begin();}
  BasisIt end() {return basis.end();}
  friend std::ostream& operator<<(std::ostream& os,  PhononBasis& phononBasis)
   {
     assert(phononBasis.dim==phononBasis.basis.size());
     BasisIt it1= phononBasis.basis.begin();
   
     while( it1!=phononBasis.basis.end())
       {
	 os<< "id"<< '\t' << std::get<toBasisType(BasisInfoField::id)>(*it1) << '\t'<<"position"<< '\t'<<std::get<toBasisType(BasisInfoField::position)>(*it1)<<  '\t'<< "state"<< '\t' << std::get<toBasisType(BasisInfoField::state)>(*it1);
     it1++;
	  }

     return os;
   }

   
  };

template<size_t L>
  struct BosonBasis
{
  using Lattice= BosonState<L>;
  using BasisType= std::tuple<double, size_t, Lattice>;
  using Basis=std::set<BasisType, CompareState<BasisType>>;
  using LatticeIt=typename Lattice::LatticeIt;  
  using BasisIt= typename  Basis::iterator;

  BosonBasis(size_t maxBosons);
  Basis basis;
  const size_t numberOfParticles; // maximum number of bosons on one site
  size_t dim;
   const BasisType operator[] (double id) {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return *basis.find(aState);
  }
  BasisIt find(double id)
  {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return basis.find(aState);
  }
  BasisIt begin() {return basis.begin();}
  BasisIt end() {return basis.end();}
  size_t particlesAt(double id, size_t site)
  {
    Lattice aLattice;
    BasisType aState=*basis.find({id, 0, aLattice});
    return std::get<toBasisType(BasisInfoField::state)>(aState)[site];
    
  }
  friend std::ostream& operator<<(std::ostream& os,  BosonBasis& bosonBasis)
   {
     assert(bosonBasis.dim==bosonBasis.basis.size());
     BasisIt it1= bosonBasis.basis.begin();
   
     while( it1!=bosonBasis.basis.end())
       {
os<< "id"<< '\t' << std::get<toBasisType(BasisInfoField::id)>(*it1) << '\t'<<"position"<< '\t'<<std::get<toBasisType(BasisInfoField::position)>(*it1)<<  '\t'<< "state"<< '\t' << std::get<toBasisType(BasisInfoField::state)>(*it1);
     it1++;
	  }

     return os;
   }

   
  };


template<class LeftBasis, class RightBasis>
 struct TensorProduct{
  using BasisType= std::tuple<size_t, double, double>;
  using Basis=std::set<BasisType, CompareState<BasisType>>;
  using BasisIt= typename Basis::iterator;
  // using LLattice= typename  LeftBasis::Lattice;
  // using RLattice= typename RightBasis::Lattice;
  
   
 
  
   TensorProduct( LeftBasis&  lbasis,   RightBasis& rbasis ) : sites(lbasis.sites), dim(lbasis.dim*rbasis.dim), dimfirst(lbasis.dim), dimsecond(rbasis.dim), lbasis(lbasis), rbasis(rbasis)  {
          auto it1= lbasis.basis.begin();

	
		 
       		while(it1!=lbasis.basis.end())
       		  {
       		    auto  it2= rbasis.basis.begin();
      		    
 		    while(it2!=rbasis.basis.end()){
		      //	 BasisIds 
 		
		      basis.insert({std::get<toBasisType(BasisInfoField::position)>(*it2)*lbasis.dim + std::get<toBasisType(BasisInfoField::position)>(*it1), std::get<toBasisType(BasisInfoField::id)>(*it1), std::get<toBasisType(BasisInfoField::id)>(*it2)});

 				   ++it2;
				
	
       		  }
 		 ++it1;
 		 
 		  }
	
 		assert(dim==basis.size());
		
       }
   const size_t sites;    
   size_t dim;
   size_t dimfirst;
   size_t dimsecond;
   LeftBasis lbasis;
   RightBasis rbasis;
   Basis basis;
   BasisIt begin() {return basis.begin();}
   BasisIt end() {return basis.end();}
 BasisIt find(double id)
  {
    
    BasisType aState(id, 0, 0);
    return basis.find(aState);
  }  
     friend std::ostream& operator<<(std::ostream& os,  TensorProduct& tensorProduct)
    {
      assert(tensorProduct.dim==tensorProduct.basis.size());
      
      BasisIt it1= tensorProduct.basis.begin();
   
      while( it1!=tensorProduct.basis.end())
        {
	  os<< "id and position "<< '\t' << std::get<toBasisType(TPBasisInfoField::position)>(*it1) << '\t' << "left basis id" << '\t'<<std::get<toBasisType(TPBasisInfoField::lId)>(*it1)<< '\t' <<"right basis id" <<'\t' <<std::get<toBasisType(TPBasisInfoField::rId)>(*it1) << '\n';
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
  basis.insert({0, 0, stateArray});
  size_t numberOfPhonons=0;
  dim++;
      // generating all states of the form 0000...0phononnumber
    while(numberOfPhonons<maxPhonons)
    {
      numberOfPhonons++;
      stateArray.fill(0);
      stateArray[L-1]= numberOfPhonons;		 
      do
	{  
	Lattice newState{stateArray};
	basis.insert({  newState.GetId(), dim,  newState});
	dim++;

       }
      while(std::next_permutation(stateArray.begin(), stateArray.end()));
     
      
   		  // Permutating over the rest of the states
  
      LatticeIt it1= stateArray.begin();
      while(stateArray[0]<stateArray[L-1])
        {
	  if(*it1==*(it1+1))                    // if the number of phonons at size i ar
                                    //equal to the number on site i+1 it set it to
                                    // zero and iterate one forward
        {
	  *it1=0;
	  it1++;
	}
	  else                                      // if the number is smaller, i increase it
    		                                       // by one and go back to site 0
        {
	  *it1+=1; 
	  it1=stateArray.begin();
    do
      {
	Lattice newState{stateArray};
	basis.insert({ newState.GetId(), dim,   newState});
	dim++;
      }
    while(std::next_permutation(stateArray.begin(), stateArray.end()));
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
    {
      numberofbosons++;
      statearray.fill(0);
      statearray[L-1]= numberofbosons;		 
     do
       {
	 if(static_cast<size_t>(std::accumulate(statearray.begin(), statearray.end(), 0))==numberOfParticles) {
	   Lattice newState{statearray};
	   basis.insert({newState.GetId(), dim,  newState});
	   dim++;
	}

       }while(std::next_permutation(statearray.begin(), statearray.end()));    
     
    		  // Permutating over the rest of the states

           LatticeIt it1= statearray.begin();
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
     Lattice newState{statearray};
      basis.insert({newState.GetId(), dim, newState});
       
           dim++;
       }
              }while(std::next_permutation(statearray.begin(), statearray.end()));
    }   
      }
 }
   
}




