#pragma once
#include<fstream>
#include<iostream>
#include<string>
#include<utility>
#include <iomanip>
#include <vector>
#include <Eigen/Dense> 
namespace Many_Body{
  using namespace Eigen;
  template<typename T>
  void ToFile(const T& x, const T& y, const std::string s1, const size_t N, bool Verbose=false)
{
  // if(x.n_elem!= y.n_elem)
  //   {
  //     std::cout << "Vector dimensions are different " << std::endl;
  //   }

    std::fstream F;
    F.open(s1, std::ios::out| std::ios::app);
    for(size_t i=0; i<N; i++ )
    {
      F<<std::fixed <<  std::setprecision(15)<< x[i] << '\t' << std::setprecision(15) << y[i] << '\n';
    }
  F.close();
  if(Verbose)
    {
    std::cout << "yout data should be stored in the dirctory as " << s1  << std::endl;
    }
return;
}
    template<typename T>
    void ToFile(const T& x, const std::string s1, const size_t N, bool Verbose=false)
{
  // if(x.n_elem!= y.n_elem)
  //   {
  //     std::cout << "Vector dimensions are different " << std::endl;
  //   }

    std::fstream F;
    F.open(s1, std::ios::out| std::ios::app);
    for(size_t i=0; i<N; i++ )
    {
      F<<std::fixed <<  std::setprecision(15)<< x[i] << '\n';
    }
  F.close();
if(Verbose)
    {
  std::cout << "yout data should be stored in the dirctory as " << s1  << std::endl;
    }
return;
}
    template<typename T>
  size_t get_size(const T& v)
  {
    return v.size();
  }

  size_t get_size(const Eigen::VectorXd& v)
  {
    return v.cols()*v.rows();
  }
  // size_t get_size(const Eigen::MatrixXd& v)
  // {
  //   return v.cols()*v.rows();
  // }
    size_t get_size(const Eigen::MatrixXcd& v)
  {
    return v.cols()*v.rows();
  }
  template<typename T>
  void bin_read(const std::string filename, T& X )
  {
        using value_type=typename T::Scalar;
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	in.read((char*) X.data(), get_size(X)*sizeof(value_type));
	return;
  }
   template<typename T>
   void bin_read(const std::string filename, std::vector<T>& X )
  {


	std::ifstream in(filename, std::ios::in | std::ios::binary);

	in.read((char*) X.data(), X.size()*sizeof(T));
	return;
  }
  template<typename T>
  void bin_write(const std::string& filename, const T& X)
  {
    using value_type=typename T::Scalar;
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    if(!out){
      std::cout<< "can not open file";
    }
    out.write((char*) X.data(), get_size(X)*sizeof(value_type));
    out.close();
    return;
  }
  template<typename T>
  void bin_write(const std::string& filename, const std::vector<T> X)
  {
    //    using value_type=typename T::Scalar;
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    if(!out){
      std::cout<< "can not open file";
    }
    out.write((char*) X.data(), get_size(X)*sizeof(T));
    out.close();
    return;
  }
}
