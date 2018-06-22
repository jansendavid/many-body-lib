#pragma once

enum class BasisInfoField{ id, position, state};
enum class TPBasisInfoField{ position, lId, rId};

template<typename E>
constexpr auto toBasisType(E enumerator) noexcept
{
  return static_cast<typename std::underlying_type_t<E>>(enumerator);
}
 template<typename T>
 size_t position( T state)
 {
   
   return static_cast<size_t>(std::get<toBasisType(BasisInfoField::position)>(state));
 }

 size_t position( std::tuple<size_t, double, double> state)
 {
   
   return static_cast<size_t>(std::get<toBasisType(TPBasisInfoField::position)>(state));
 }
template<typename T>
 double id( T state)
 {
   
   return std::get<toBasisType(BasisInfoField::id)>(state);
 }
double id( std::tuple<size_t, double, double> state)
 {
   
   return std::get<toBasisType(TPBasisInfoField::position)>(state);
 }
double leftId( std::tuple<size_t, double, double> state)
 {
   
   return std::get<toBasisType(TPBasisInfoField::lId)>(state);
 }
double rightId( std::tuple<size_t, double, double> state)
 {
   
   return std::get<toBasisType(TPBasisInfoField::rId)>(state);
 }

template< typename T>
auto& state( const T& state)
 {
   
   return std::get<toBasisType(BasisInfoField::state)>(state);
 }
