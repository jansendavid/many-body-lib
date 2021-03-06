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
#include <algorithm> 
#include<iterator>
#include<numeric>
#include<tuple>
#include<set>
#include<type_traits>
#include"accesfunctions.hpp"
#include"numerics.hpp"
// computing the factorial
namespace Many_Body{
// template<typename State>
// std::tuple<size_t, size_t, State> getState(size_t id)
// {
//   State aState(0);
//   return {id, 0, aState};
// }
template<typename State>
struct CompareState{
     bool operator()(const State& lhs, const State& rhs) const
     {
       return std::get<0>(lhs)< std::get<0>(rhs);
     }
     
   };
  struct isFermion{
    bool isFer{1};
  };
    struct isBoson{
    bool isFer{0};
    };
  struct stateBase{
    stateBase(): sites(0), maxParNr{0}  {};


    stateBase(std::vector<size_t> state, size_t MaxParNr): sites(state.size()), maxParNr( MaxParNr)  {

         id=makeId(state);
   };
    stateBase(size_t sites, size_t maxParNr): sites(sites), maxParNr(maxParNr), id(0) {};
        stateBase(size_t sites): sites(sites), maxParNr(0), id(0) {};


    size_t sites{0};
    size_t maxParNr{0};
    size_t id{0};

    
    // member functions
    virtual std::vector<size_t> makeStateVec() const
   {
     std::vector<size_t> state(sites, 0);
          std::vector<size_t>::reverse_iterator rit = state.rbegin();
     	  size_t IdCop=id;
     for (; rit!= state.rend(); ++rit)
       {
	 
     	 *rit=IdCop%(maxParNr+1);
     	 IdCop/=(maxParNr+1);
       }
	  return state;

   }
    void setId(size_t Id){
      id=Id;
    }
     inline size_t GetId()
   {

     return id;
              }
      virtual inline size_t makeId(std::vector<size_t>& state)
   {
     size_t sum=0;
     size_t mult=1;


      std::vector<size_t>::reverse_iterator rit = state.rbegin();

  for (; rit!= state.rend(); ++rit)
    {

      sum+=(*rit)*mult;
      mult*=(maxParNr+1);
       }
  
     return sum;
     
              }
          inline void setPartNr(size_t site, size_t nr)
    {
      assert(nr<=maxParNr && nr>=0);
      std::vector<size_t> state=makeStateVec();



        state[site]=nr;


	id=makeId(state);

    }
    inline void switchPartNr(size_t site_1, size_t nr_1, size_t site_2, size_t nr_2)
    {

      std::vector<size_t> state=makeStateVec();




        state[site_1]=nr_1;
	        state[site_2]=nr_2;



	id=makeId(state);

    }
       inline size_t Count() const
   {
     std::vector<size_t> state=makeStateVec();
     size_t sum = std::accumulate(state.begin(), state.end(), 0);
     return sum;
   }

    
    auto operator[] ( size_t i) const{
  std::vector<size_t> state=makeStateVec();
  return (state[i]);
    }
       friend std::ostream& operator<<(std::ostream& os, const  stateBase& state) 
   {
     std::vector<size_t> stateVec=state.makeStateVec();
     for(size_t i=0; i<state.sites; i++)
       { os<< stateVec[i] << std::setw(3);
   	   }
     os << '\n';
     return os;
   }
  };

  struct ElectronState: public stateBase
 {
 using LatticeIt=typename std::vector<size_t>::iterator;
 using Const_LatticeIt=typename std::vector<size_t>::const_iterator;
// using stateBase::Count;
   //    using stateBase::setPartNr;
   ElectronState(): stateBase()  {};
   ElectronState(std::vector<size_t> state): stateBase(state, 1)  {
         
   };

   ElectronState( size_t sites, size_t stateInt):  stateBase(sites, 1) {
    id=stateInt;

    }
    


   void flip(size_t i)
   {
     std::vector<size_t> state=makeStateVec();

     state[i]=(state[i]+1)%2;


     id=makeId(state);

   }


 };
  struct OneElectronState: public stateBase
 {
 using LatticeIt=typename std::vector<size_t>::iterator;
 using Const_LatticeIt=typename std::vector<size_t>::const_iterator;
// using stateBase::Count;
   //    using stateBase::setPartNr;
   OneElectronState(): stateBase()  {};
   OneElectronState(std::vector<size_t> state): stateBase(state.size())  {
     maxParNr=1;

     if(std::accumulate(state.begin(), state.end(),0)>1)
       {std::cout<< "error"<<'\n';}
     size_t n=0;
     
for(long int i=sites-1; i>=0; i--)
        {
	 //	 std::cout<< i<< "x "<<std::endl;
	 if(state[i]==1){
	   id=n;
	 }
	 n++;
       }
   };

   OneElectronState( size_t sites, size_t stateInt): stateBase(sites) {
    id=stateInt;
    maxParNr=1;

    }
    
std::vector<size_t> makeStateVec() const override
   {
     std::vector<size_t> state(sites, 0);
          

	  state[sites-1-id]=1;
	  return state;

   }
   inline void setPartNr(size_t site, size_t nr)=delete;
       inline size_t makeId(std::vector<size_t>& state) override
   {

     size_t sum=0;
        //   if(std::accumulate(state.begin(), state.end(),0)>1)
       // {std::cout<< "x error"<<'\n';}
     size_t n=0;

     for(int i=state.size()-1; i>=0; i--)
       {

	 if(state[i]==1){
	   sum=n;
	 }
	 n++;
       }
     
     return sum;
     
              }

   void flip(size_t i)
   {
     std::vector<size_t> state=makeStateVec();

     state[i]=(state[i]+1)%2;


     id=makeId(state);

   }


 };

  struct BosonState: public stateBase
  {

    using LatticeIt=typename std::vector<size_t>::iterator;
        using Const_LatticeIt=typename std::vector<size_t>::const_iterator;
    BosonState(size_t sites, size_t BosonNr): stateBase(sites, BosonNr) {};
    BosonState(): stateBase() {};
   
    BosonState( std::vector<size_t>& state, size_t BosonNr): stateBase(state, BosonNr) {  

      };

       friend std::ostream& operator<<(std::ostream& os, const  BosonState& aState)
   {

     std::vector<size_t> state=aState.makeStateVec();
     for(size_t i=0; i<aState.sites; i++)
       {
	 os<< aState[i]  << std::setw(3);
   	   }
     os << '\n';
     return os;
   }

  };

struct ElectronBasis
{
  using Lattice=ElectronState;
  using BasisType= std::tuple< size_t, size_t, Lattice>;
  using LatticeIt=typename Lattice::LatticeIt;
  using Const_LatticeIt=typename Lattice::Const_LatticeIt;
  using Basis=std::set<BasisType, CompareState<BasisType>>;
  using BasisIt= typename  Basis::iterator;
  using Const_BasisIt= typename  Basis::const_iterator;

   
  ElectronBasis( size_t sites, size_t numberOfParticles): totalmaxPar(numberOfParticles), sites(sites), dim(0) {
    for(size_t i=0; i<static_cast<size_t>(std::pow(2, sites)); i++)
        {

	   
   	  ElectronState newState( sites, i);
            if(newState.Count()==numberOfParticles)
             {
	      
   	      basis.insert({newState.GetId(), dim, newState});
	      dim++;    
             }
        
   }
  }
  ElectronBasis( size_t sites, std::vector<int> spaces): totalmaxPar(*std::max_element(spaces.begin(), spaces.end())),
 sites(sites), dim(0) {
    for(size_t i=0; i<static_cast<size_t>(std::pow(2, sites)); i++)
        {

	   
   	  ElectronState newState( sites, i);
	  if(std::count(spaces.begin(), spaces.end(), newState.Count()))
             {
	      
   	      basis.insert({newState.GetId(), dim, newState});
	      dim++;    
             }
        
   }
  }
  ElectronBasis(size_t sites ): totalmaxPar(sites), sites(sites), dim(0){
    
      for(size_t i=0; i<static_cast<size_t>(std::pow(2, sites)); i++)
        {
   	  ElectronState newState( sites, i);
	  //	  std::cout<< newState<<std::endl;
             
   	     
	  basis.insert({newState.GetId(), i, newState});
	       dim++;
	       
             }
     
        
    }
  ElectronBasis(ElectronState state): sites(state.sites), totalmaxPar(state.Count()), dim(0){             
   	     
   	      basis.insert({state.GetId(), 0, state});
	       dim++;    
            
   }
  void append(size_t numberOfParticles)
  {
    // use this function to append states to an already exisiting basis 
    for(size_t i=0; i<static_cast<size_t>(std::pow(2, sites)); i++)
        {

	   
   	  ElectronState newState( sites, i);
            if(newState.Count()==numberOfParticles)
             {
	      
   	      basis.insert({newState.GetId(), dim, newState});
	      dim++;    
             }
  }
    totalmaxPar=std::max(totalmaxPar, numberOfParticles);
    
  }
  void insert(Lattice& newState) 
   	      {
		size_t i=newState.GetId();
		if (basis.count({i, dim, newState}) > 0)
	{
		std::cout << "'already exist" << std::endl;
		return;
	}
	else
	{
 basis.insert({i, dim, newState});
	      dim++;
	      return;}
	      }
   ElectronBasis(const ElectronBasis& e)= default;
   friend std::ostream& operator<<(std::ostream& os,  const ElectronBasis& electronBasis)
    {
      assert(electronBasis.dim==electronBasis.basis.size());
      Const_BasisIt it1= electronBasis.begin();
   
     while( it1!=electronBasis.basis.end())
       {
	 os<< "id"<< '\t' << std::get<toBasisType(BasisInfoField::id)>(*it1) << '\t'<<"position"<< '\t'<<std::get<toBasisType(BasisInfoField::position)>(*it1)<<  '\t'<< "state"<< '\t' << std::get<toBasisType(BasisInfoField::state)>(*it1);
     it1++;
	  }

     return os;
   }
    const BasisType operator[] (size_t id) {
      Lattice aLattice(sites, maxParticles);
    BasisType aState(id, 0, aLattice);
    return *basis.find(aState);
  }
    
   size_t particlesAt(size_t id, size_t site) const
  {
    Lattice aLattice(sites, maxParticles);
    BasisType aState=*basis.find({id, 0, aLattice});
    return std::get<toBasisType(BasisInfoField::state)>(aState)[site];
    
  }
     void flip(size_t id, size_t site) const
  {
    Lattice aLattice(sites, maxParticles);
    BasisType aState=*basis.find({id, 0, aLattice});
    return std::get<toBasisType(BasisInfoField::state)>(aState).flip(site);
    
  }
  BasisIt begin() {return basis.begin();}
  BasisIt end() {return basis.end();}
  Const_BasisIt begin() const {return basis.begin();}
  Const_BasisIt end() const {return basis.end();}
  BasisIt find(size_t id)
  {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return basis.find(aState);
  }
    Const_BasisIt find(size_t id) const
  {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return basis.find(aState);
  }
   Basis basis;
  const size_t maxParticles=1;
  size_t totalmaxPar=0;
   size_t sites;
   size_t dim;
//   //  const size_t numberofparticles;
  
 };

struct OneElectronBasis
{
  using Lattice=OneElectronState;
  using BasisType= std::tuple< size_t, size_t, Lattice>;
  using LatticeIt=typename Lattice::LatticeIt;
  using Const_LatticeIt=typename Lattice::Const_LatticeIt;
  using Basis=std::set<BasisType, CompareState<BasisType>>;
  using BasisIt= typename  Basis::iterator;
  using Const_BasisIt= typename  Basis::const_iterator;
   OneElectronBasis(size_t sites ):   sites(sites), dim(0){
     size_t n=0;

     for(long int i=0; i<sites; i++)
        {

	  std::vector<size_t> v(sites, 0);
	  v[sites-i-1]=1;


	     	  OneElectronState newState(v);


   	     
	  basis.insert({i, i, newState});
	  n++;
	  dim++;    
             }
        
    }

   OneElectronBasis(const OneElectronBasis& e)= default;
   friend std::ostream& operator<<(std::ostream& os,  const OneElectronBasis& electronBasis)
    {
      assert(electronBasis.dim==electronBasis.basis.size());
      Const_BasisIt it1= electronBasis.begin();
   
     while( it1!=electronBasis.basis.end())
       {
	 os<< "id"<< '\t' << std::get<toBasisType(BasisInfoField::id)>(*it1) << '\t'<<"position"<< '\t'<<std::get<toBasisType(BasisInfoField::position)>(*it1)<<  '\t'<< "state"<< '\t' << std::get<toBasisType(BasisInfoField::state)>(*it1);
     it1++;
	  }

     return os;
   }
    const BasisType operator[] (size_t id) {
      Lattice aLattice(sites, maxParticles);
    BasisType aState(id, 0, aLattice);
    return *basis.find(aState);
  }
    
   size_t particlesAt(size_t id, size_t site) const
  {
    Lattice aLattice(sites, maxParticles);
    BasisType aState=*basis.find({id, 0, aLattice});
    return std::get<toBasisType(BasisInfoField::state)>(aState)[site];
    
  }
     void flip(size_t id, size_t site) const
  {
    Lattice aLattice(sites, maxParticles);
    BasisType aState=*basis.find({id, 0, aLattice});
    return std::get<toBasisType(BasisInfoField::state)>(aState).flip(site);
    
  }
  BasisIt begin() {return basis.begin();}
  BasisIt end() {return basis.end();}
  Const_BasisIt begin() const {return basis.begin();}
  Const_BasisIt end() const {return basis.end();}
  BasisIt find(size_t id)
  {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return basis.find(aState);
  }
    Const_BasisIt find(size_t id) const
  {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return basis.find(aState);
  }
   Basis basis;
  const size_t maxParticles=1;
   size_t sites;
   size_t dim;
//   //  const size_t numberofparticles;
  
 };


   struct PhononBasis
 {
   using Lattice= BosonState;  
  using BasisType= std::tuple<size_t, size_t, Lattice>;
  using Basis=std::set<BasisType, CompareState<BasisType>>;
  using LatticeIt= typename Lattice::LatticeIt;
  using Const_LatticeIt= typename Lattice::Const_LatticeIt;
  using BasisIt= typename  Basis::iterator;
  using Const_BasisIt= typename  Basis::const_iterator;

   PhononBasis( size_t sites, size_t maxPhonons);
    PhononBasis(size_t sites): sites(sites), maxParticles(0) {};

   Basis basis;
   size_t dim;
   const size_t sites;
   const size_t maxParticles;
  
   const BasisType operator[] (size_t id) {
     Lattice aLattice(sites, maxParticles);
    BasisType aState(id, 0, aLattice);
    return *basis.find(aState);
  }
   size_t particlesAt(size_t id, size_t site)  const
  {

    Lattice aLattice;
    BasisType aState=*basis.find({id, 0, aLattice});


    return std::get<toBasisType(BasisInfoField::state)>(aState)[site];
    
  }
  
  BasisIt find(size_t id)
  {
    Lattice aLattice(sites, maxParticles);
    BasisType aState(id, maxParticles, aLattice);
    return basis.find(aState);
  }
    Const_BasisIt find(size_t id) const
  {
    Lattice aLattice;
    BasisType aState(id, 0, aLattice);
    return basis.find(aState);
  }
  BasisIt begin() {return basis.begin();}
  BasisIt end() {return basis.end();}
  Const_BasisIt begin() const {return basis.begin();}
  Const_BasisIt end() const {return basis.end();}
  friend std::ostream& operator<<(std::ostream& os,  PhononBasis& phononBasis)
   {
     assert(phononBasis.dim==phononBasis.basis.size());
     Const_BasisIt it1= phononBasis.begin();
   
     while( it1!=phononBasis.basis.end())
       {
 	 os<< "id"<< '\t' << std::get<toBasisType(BasisInfoField::id)>(*it1) << '\t'<<"position"<< '\t'<<std::get<toBasisType(BasisInfoField::position)>(*it1)<<  '\t'<< "state"<< '\t' << std::get<toBasisType(BasisInfoField::state)>(*it1);
     it1++;
 	  }

     return os;
   }

   
   };



template<class LeftBasis, class RightBasis>
 struct TensorProduct{
  using LB= LeftBasis;
  using RB= RightBasis;
  using BasisType= std::tuple<size_t, size_t, size_t>;
  using Basis=std::set<BasisType, CompareState<BasisType>>;
  using BasisIt= typename Basis::iterator;
  using Const_BasisIt= typename Basis::const_iterator;
  using LeftBasisIt= typename LeftBasis::BasisIt;
  using RightBasisIt= typename RightBasis::BasisIt;
  using Const_LeftBasisIt= typename LeftBasis::Const_BasisIt;
  using Const_RightBasisIt= typename RightBasis::Const_BasisIt;
  using RightLattice= typename RightBasis::Lattice;
  using LeftLattice= typename LeftBasis::Lattice;
  
  
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
  Const_BasisIt begin() const {return basis.begin();}
  Const_BasisIt end() const {return basis.end();}
  BasisIt find(size_t id)
  {
    
    BasisType aState(id, 0, 0);
    return basis.find(aState);
  }
    Const_BasisIt find(size_t id) const
  {
    
    BasisType aState(id, 0, 0);
    return basis.find(aState);
  }  
     friend std::ostream& operator<<(std::ostream& os,  TensorProduct& tensorProduct)
    {
      assert(tensorProduct.dim==tensorProduct.basis.size());
      
      Const_BasisIt it1= tensorProduct.begin();
   
      while( it1!=tensorProduct.basis.end())
        {
	  os<< "id and position "<< '\t' << std::get<toBasisType(TPBasisInfoField::position)>(*it1) << '\t' << "left basis id" << '\t'<<std::get<toBasisType(TPBasisInfoField::lId)>(*it1)<< '\t' <<"right basis id" <<'\t' <<std::get<toBasisType(TPBasisInfoField::rId)>(*it1) << '\n';
      it1++;
      	  }

      return os;
    }


  };


//   // MEMBER FUNCTION:


  PhononBasis::PhononBasis( size_t sites, size_t maxPhonons): dim{0},  sites(sites), maxParticles(maxPhonons)  {
     std::vector<size_t> stateArray(sites, 0);
     basis.insert({0, 0, {stateArray, maxPhonons}});
  size_t numberOfPhonons=0;

  dim++;
//       // generating all states of the form 0000...0phononnumber
    while(numberOfPhonons<maxPhonons)
    {
      numberOfPhonons++;
      std::fill(stateArray.begin(), stateArray.end(), 0);
      stateArray[sites-1]= numberOfPhonons;		 
      do
    	{
	  
	  Lattice newState{stateArray, maxPhonons};
	  //	  std::cout<< newState<<std::endl;

	//	std::cout<< newState;
    	basis.insert({  newState.GetId(), dim,  newState});
    	dim++;

       }
      while(std::next_permutation(stateArray.begin(), stateArray.end()));
     
      
    // 		  // Permutating over the rest of the states
  
      LatticeIt it1= stateArray.begin();
      while(stateArray[0]<stateArray[sites-1])
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

    	Lattice newState{stateArray, maxPhonons};

    	basis.insert({ newState.GetId(), dim,   newState});
    	dim++;
      }
    while(std::next_permutation(stateArray.begin(), stateArray.end()));
    }
    	}
 }

   
 }



template<typename TotalBasis>
auto generateMomBasis(TotalBasis& totalBasis, int k)->TotalBasis
{
  TotalBasis newBasis;
  int sites=totalBasis.sites;

  for(auto tpState : totalBasis)
    {
      //      const Lattice state=GetLattice(tpState);
      for(size_t i=0; i<sites; i++)
	{
	  		   
			   // Lattice state2=state;
			   // size_t newStateId=Many_Body::Translate(i, state2);
			   // if(newStateId ==tpState.id)
			   //   {
			   //    if(newBasis.find(newStateId)!=newBasis.end())
			   //   {
			   //     // auto it2= totalBasis.find(newStateId);
			   //     // const Lattice newState=GetLattice(*it2);
			   //     // int max=			       
			   //   } 
			   //   }
			   
			   
			   
			   
			   
			   
	}
    }
}

}
