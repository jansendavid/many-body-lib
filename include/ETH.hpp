#pragma once
#include"files.hpp"
#include"Eigen/Dense"
#include<vector>

namespace Many_Body{
template<typename Container>  
double temp( Container& Energies, double Einput, double Tmin, double Tmax, int Trials=100, double Err=1E-9 )
{
  double T{0};
    double expval{0};
    if(Tmin>0)
      {
   for(int i=0; i<Trials; i++)
      {
	
    	T=(Tmax+Tmin)/2;
	//  T=2.17743;
    double Z{0};
    double sum{0};  
    for(auto&  en: Energies)
      {
	
	double ex=std::exp(-(en-Energies[0])/T);
	Z+=ex;
	sum+=ex*en;
      }
     expval= sum/Z;
     if(std::abs(Einput-expval)<Err ){
       std::cout<< "breaked "<<std::endl; break;}
      
     if(expval>Einput){ Tmax=T;}
    else{Tmin=T;}
        std::cout<<"Exp "<< expval << " Ein "<< Einput << " diff "<<  std::abs(Einput-expval) << " for T "<< T<< " Tmax "<< Tmax << " tmin "<< Tmin<< std::endl;
      }}
       else
      {
   for(int i=0; i<Trials; i++)
      {
	
    	T=(Tmax+Tmin)/2;
	//  T=2.17743;
    double Z{0};
    double sum{0};  
    for(auto&  en: Energies)
      {
	
	double ex=std::exp(-en/T);
	Z+=ex;
	sum+=ex*en;
      }
     expval= sum/Z;
     if(std::abs(Einput-expval)<Err ){
       std::cout<< "breaked "<<std::endl; break;}
      
     if(expval>Einput){ Tmin=T;}
    else{Tmax=T;}
        std::cout<<"Exp "<< expval << " Ein "<< Einput << " diff "<<  std::abs(Einput-expval) << " for T "<< T<< " Tmax "<< Tmax << " tmin "<< Tmin<< std::endl;
      }}
       std::cout<< "Final result "<<expval << " for T "<< T<<std::endl;
   return T;
}
  template<typename Container>  
  double expvalCan( Container& Energies, Container& Observabel,  double T)
{
    double expval{0};
    double Z{0};
    double sum{0};
    double beta=1./T;
    for(int i=0; i< Energies.size(); i++)
      {
	double exponent=-(Energies[i]-Energies[0]);
	exponent*=beta;
	double ex=std::exp(exponent);
	Z+=ex;
	sum+=ex*Observabel[i];
      }
     expval= sum/Z;
   return expval;
}
    template<typename Container>  
  double entropy( Container& Energies,   double T)
{
    double expval{0};
    double Z{0};
    double ent{0};
        for(int i=0; i< Energies.size(); i++)
      {
	
	double ex=std::exp(-Energies[i]/T);
	Z+=ex;
      }
    for(int i=0; i< Energies.size(); i++)
      {
	
	double ex=(1./Z)*std::exp(-Energies[i]/T);

	ent-=ex*std::log(ex);
      }
    
   return ent;
}
}
