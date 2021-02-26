#include"files.hpp"
#include<iostream>
#include<eigen3/Eigen/Dense>
#include<string>
#include<vector>
#include<cmath>
#include <cstdlib>
#include <numeric>
#include <iterator>
#include <algorithm>
using namespace Many_Body;
template <typename Container, typename T = typename std::decay<decltype(*std::begin(std::declval<Container>()))>::type>
T variance(Container && c)
{
    auto b = std::begin(c), e = std::end(c);
    auto size = std::distance(b, e);


    auto sum = std::accumulate(b, e, T());
    auto mean = sum / size;
    T accum = T();
    for (const auto d : c)
        accum += (d - mean) * (d - mean);
    return accum / (size );
}
void absav(size_t N, std::string filenameobs, std::string filenameev, double dE, double av)
{
  
  Eigen::MatrixXcd O(N, N);
  Eigen::VectorXd E(N);
  bin_read(filenameobs, O);
  bin_read(filenameev, E);
  std::cout<< E.rows() << "  " << E.cols() << std::endl;
  std::cout<< O.rows() << "  " << O.cols() << std::endl;
  E=E/av;

  std::vector<double> omega;
  std::vector<double> offdfunc;
  std::vector<double> offd;
  std::string w=std::to_string(dE);
  
  w=w.substr(0, 6);
  for(size_t i=0; i<N-1; i++)
    {
      for(size_t j=i+1; j<N; j++)
    {
      
      if(std::abs<double>(std::abs<double>((E[i]+E[j]))-1.)<dE)
  	{
  	  offd.push_back(std::abs<double>(O(i, j)));
  	  offd.push_back(std::abs<double>(O(j, i)));

  	      omega.push_back((E[i]-E[j]));

  		omega.push_back((E[j]-E[i]));
	      
	      
  	}
    }
    }
    
  bin_write("w"+w+"omega" +filenameobs, omega);
  //  bin_write("w"+w+"offdfunc"+filenameobs, offdfunc);
  bin_write("w"+w+"offd" +filenameobs, offd);
    double average = std::accumulate( offd.begin(), offd.end(), 0.0)/offd.size();
  std::cout<< average<<std::endl;
    
}
void varrat(size_t N, std::string filenameobs, size_t matSize)
{
  
  Eigen::MatrixXcd O(N, N);

  bin_read(filenameobs, O);
  
  std::vector<double> res;
  std::vector<double> step;
  // Eigen::VectorXcd v()
  std::string w=std::to_string(matSize);
  

  for(size_t i=0; i<(N-matSize); i++)
    {
      std::vector<double> diagreal;
  std::vector<double> offdreal;
  std::vector<double> offdim;

      Eigen::MatrixXcd V=O.block(i, i, matSize, matSize);      
      for(size_t l=0; l<V.rows(); l++)
     {      diagreal.push_back(real(V(l, l)));
     
       for (size_t m = l+1; m < V.cols();  m++)
 	{
	  
	  
 offdreal.push_back(real(V(l, m)));
   	  offdreal.push_back(real(V(m, l)));

 	  offdim.push_back(imag(V(m, l)));
 	  offdim.push_back(imag(V(l, m)));
  	}
     }
      // std::cout<< offdim.size()<< std::endl;
      // std::cout<< offdreal.size()<< std::endl;
      // std::cout<< diagreal.size()<< std::endl;
       step.push_back(i);
       res.push_back(variance(diagreal)/(variance(offdreal)+variance(offdim)));

       
          }
   // for(int i=0; i<step.size(); i++)
   //   { std::cout<< step[i] << "  " << res[i] << std::endl;}
  bin_write("S"+w+"steps" +filenameobs, step);
  bin_write("S"+w+"res"+filenameobs, res);
  // // bin_write("w"+w+"offd" +filenameobs, offd);
      double average = std::accumulate( res.begin(), res.end(), 0.0)/res.size();

   std::cout<< average<<std::endl;
    
}
int main(int argc, char *argv[])
{

  std::string filenameobs=argv[1];

  size_t M=std::atoi(argv[2]);
  size_t L=std::atoi(argv[3]);
  //  double dE=std::atof(argv[5]);
   size_t S=std::atoi(argv[4]);
  size_t N=std::pow((M+1), L);
  double av=M*0.5*L;
  //  absav( N,  filenameobs, filenameev, dE, av);
   varrat( N,  filenameobs, S);
  return 0;
}
