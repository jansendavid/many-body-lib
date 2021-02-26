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
using namespace std;
using namespace Many_Body;
int main(int argc, char *argv[])
{
  int M=3;
  double omega=0.5;
  int L=8;
  double av=L*M*omega*0.5;
  int N=std::pow(M+1, L);

  
  string filename="EdblHL"+ std::to_string(L) +"M" + std::to_string(M) +"t0_1.0gam0.7omg0.5k0.bin";
  vector<double> E1(N);
   bin_read(filename, E1);
   // std::move ( E.begin(), E.end(), Etot.begin() );
 filename="EdblHL"+ std::to_string(L) +"M" + std::to_string(M) +"t0_1.0gam0.7omg0.5k1.bin";
   vector<double>E2(N);
  bin_read(filename, E2);
   vector<double>E3(N);
  bin_read(filename, E3);
  filename="EdblHL"+ std::to_string(L) +"M" + std::to_string(M) +"t0_1.0gam0.7omg0.5k2.bin";
   vector<double>E22(N);
  bin_read(filename, E22);
   vector<double>E33(N);
  bin_read(filename, E33);
  // filename="EdblHL"+ std::to_string(L) +"M" + std::to_string(M) +"t0_1.0gam0.7omg0.5k3.bin";
  // vector<double>E222(N);
  //  bin_read(filename, E222);
  //  vector<double>E333(N);
  // bin_read(filename, E333);

  // filename="EdblHL"+ std::to_string(L) +"M" + std::to_string(M) +"t0_1.0gam0.7omg0.5k4.bin";
  //    vector<double> E4(N);
  //   bin_read(filename, E4);

    
       vector<double> Etot;
       //Etot.reserve( E1.size() + E2.size()  + E3.size() + E4.size()  );
 
       //Etot.insert( Etot.end(), E1.begin(), E1.end() );
        Etot.insert( Etot.end(), E2.begin(), E2.end() );
	    // Etot.insert( Etot.end(), E3.begin(), E3.end() );
	    //        Etot.insert( Etot.end(), E22.begin(), E22.end() );
	    // Etot.insert( Etot.end(), E33.begin(), E33.end() );
	    //Etot.insert( Etot.end(), E222.begin(), E222.end() );
		     //	     Etot.insert( Etot.end(), E333.begin(), E333.end());
		     //  Etot.insert( Etot.end(), E4.begin(), E4.end() );
		std::cout<< Etot.size()<<std::endl;
		std::for_each(Etot.begin(), Etot.end(),  [&av](double& i){i=i/av;});
    auto s = std::accumulate(Etot.begin(), Etot.end(), 0.0);
    std::cout<<Etot.size()-N<< "av "<< s/Etot.size()<<std::endl;
    double Tmin=0.1;
    double Tmax=500;
    double Einput=0.6;
    std::cout<< Einput<<std::endl;
    double expval=0;
   double Err=1E-14;
    int Trials=100;
    double T{0.0};
    for(int i=0; i<Trials; i++)
      {
	
    	T=(Tmax+Tmin)/2;
	//  T=2.17743;
    double Z{0};
    double sum{0};  
    for(auto&  en: Etot)
      {
	
	double ex=std::exp(-en/T);
	Z+=ex;
	sum+=ex*en;
      }
     expval= sum/Z;
     if(std::abs(Einput-expval)<Err ){
       std::cout<< "breaked "<<std::endl; break;}
      
     if(expval>Einput){ Tmax=T;}
    else{Tmin=T;}
     //  std::cout<<"Exp "<< expval << " Ein "<< Einput << " diff "<<  std::abs(Einput-expval) << " for T "<< T<< " Tmax "<< Tmax << " tmin "<< Tmin<< std::endl;
      
       }
    std::cout<< "Final result "<<expval << " for T "<< T<<std::endl;
    //std::cout<< "Final result "<<std::abs(Einput-expval) << " for T "<< T<<std::endl;
     double Z{0};
    double entropy{0};
    //T=100;
    //T=5;
    for(auto&  en: Etot)
      {
	
	double ex=std::exp(-en/T);

		Z+=ex;
	
      }
        for(auto&  en: Etot)
      {
	
	double ex=(1./Z)*std::exp(-en/T);

	entropy-=ex*std::log(ex);
	
      }
    std::cout<< " entropy " << entropy<<std::endl;
    //    std::cout<< " Z " << Z<<std::endl;
    //    std::cout<< " real S " <<std::log(Etot.size())<<std::endl;
  return 0;
}

