#pragma once

namespace Many_Body
{


enum class BasisInfoField{ id, position, state};
enum class TPBasisInfoField{ position, lId, rId};

template<typename E>
constexpr auto toBasisType(E enumerator) noexcept
{
  return static_cast<typename std::underlying_type_t<E>>(enumerator);
}
 template<typename T>
 size_t Position( const T& state)
 {
   
   return static_cast<size_t>(std::get<toBasisType(BasisInfoField::position)>(state));
 }

 size_t Position( const std::tuple<size_t, size_t, size_t>& state)
 {
   
   return static_cast<size_t>(std::get<toBasisType(TPBasisInfoField::position)>(state));
 }
template<typename T>
 size_t Id( const T& state)
 {
   
   return std::get<toBasisType(BasisInfoField::id)>(state);
 }
size_t Id( const std::tuple<size_t, size_t, size_t>& state)
 {
   
   return std::get<toBasisType(TPBasisInfoField::position)>(state);
 }
size_t LeftId( const std::tuple<size_t, size_t, size_t>& state)
 {
   
   return std::get<toBasisType(TPBasisInfoField::lId)>(state);
 }
size_t RightId( const std::tuple<size_t, size_t, size_t>& state)
 {
   
   return std::get<toBasisType(TPBasisInfoField::rId)>(state);
 }

template< typename T>
auto GetLattice( const T& state)
 {
   
   return std::get<toBasisType(BasisInfoField::state)>(state);
 }
   template<class State, class T=double>
   T Translate(const size_t translation, const State& state)     // translates to the left
     {
       
       State newState(state);
       
       for(size_t t=0; t<translation; t++)
	 {
	   State tempState(newState);
       //       std::cout << newstate;
       for(size_t i=0; i<state.sites; i++)
	 {
	   size_t j=(i+1)%(state.sites);
	   newState[i]=tempState[j];
	 }
       
	 }

        return newState.GetId();
       
       
     }

}
