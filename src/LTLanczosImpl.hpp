#ifndef LTLANCZOSIMPL_H
#define LTLANCZOSIMPL_H
#include"LanczosHelper.hpp"
using VecD=Eigen::VectorXd;
using VecC=Eigen::VectorXcd;
using MatD=Eigen::MatrixXd;
using MatC=Eigen::MatrixXcd;
template<bool Fast, typename MatrixHam,typename Container, typename TrafoMat>
void LTiterate(int lanczosDim,  std::vector<double>& beta, MatD& obs, Container& observable, VecD& eigenVals,  VecD& Z, TrafoMat& T1, VecC& iniState,  MatrixHam& hamiltonian){
   for(int j=0; j<lanczosDim; j++)
   		    {
		      VecC qvecj;
		         if constexpr(!Fast)
				     {
      qvecj=Many_Body::lanczTrafo(T1.col(j), iniState, lanczosDim, hamiltonian);
				     }
			  else{
			 qvecj=T1.col(j);
			}
			  auto link=iniState.adjoint()*qvecj;
		      for(int l=0; l<lanczosDim; l++)
		    {
		     // auto link=iniState.adjoint()*Q.col(j);
		     VecC qvecl;
		       if constexpr(!Fast)
				     {
				       
				       qvecl=Many_Body::lanczTrafo(T1.col(l), iniState, lanczosDim, hamiltonian);
				     }
		       else{
			
			 qvecl=T1.col(l);
			
		       }
		       for(size_t i=0; i<beta.size(); i++)
			 {
 double exponent=(eigenVals[j]+eigenVals[l])/2-eigenVals[0];
		      exponent*=-beta[i];
		      auto exp=std::exp(exponent);
		      auto right=(qvecl).adjoint()*iniState;
		      if(l==j)
			{
			  Z(i)+=real(right(0, 0)*link(0, 0))*exp;
			}
		      
		       for(int m=0; m<observable.size(); m++)
		       	{
			  auto expval=((qvecj).adjoint())*observable[m]*qvecl;

			  obs(i,m)+=real((link(0, 0)*right(0,0)*expval(0, 0)))*exp;
			}
			 }
			}
		    }
}
template<bool Fast, typename MatrixHam,typename Container, typename TrafoMat>
void LTiterate(int lanczosDim,  double beta, VecD& obs, Container& observable, VecD& eigenVals,  double& Z, TrafoMat& T1, VecC& iniState,  MatrixHam& hamiltonian){
   for(int j=0; j<lanczosDim; j++)
   		    {
		      VecC qvecj;
		         if constexpr(!Fast)
				     {
      qvecj=Many_Body::lanczTrafo(T1.col(j), iniState, lanczosDim, hamiltonian);
				     }
			  else{
			 qvecj=T1.col(j);
			
			
		       }
			  auto link=iniState.adjoint()*qvecj;
		      for(int l=0; l<lanczosDim; l++)
		    {


		      double exponent=(eigenVals[j]+eigenVals[l])/2-eigenVals[0];
		      exponent*=-beta;
		      auto exp=std::exp(exponent);
		      
		      // auto link=iniState.adjoint()*Q.col(j);
		     
		     VecC qvecl;
		       if constexpr(!Fast)
				     {
				       
				       qvecl=Many_Body::lanczTrafo(T1.col(l), iniState, lanczosDim, hamiltonian);
				     }
		       else{
			
			 qvecl=T1.col(l);
			
		       }

		      auto right=(qvecl).adjoint()*iniState;
		      if(l==j)
			{
			  Z+=real(right(0, 0)*link(0, 0))*exp;
			}
		      
		       for(int i=0; i<obs.size(); i++)
		       	{
			  auto expval=((qvecj).adjoint())*observable[i]*qvecl;

			  obs[i]+=real((link(0, 0)*right(0,0)*expval(0, 0)))*exp;
			}
			}
		    }

		     
  
  
}
template<typename Matrix, typename Container>
auto  calculate_lanczLT_slow(Matrix& hamiltonian, Container& observable, double beta, int lanczosDim=20, double err=1E-10)->std::tuple<std::vector<double>, double>
{
   std::random_device rd;
std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
std::uniform_real_distribution<float> dis(-1.0, 1.0);

 
 Eigen::VectorXcd iniState=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double>(dis(gen), dis(gen));});
  auto Z=std::complex<double>{0, 0};
 
    auto A=std::complex<double>{0, 0};
      

   iniState/=iniState.norm();
   //Eigen::MatrixXcd Q(hamiltonian.rows(), lanczosDim);
   auto [S, eigenVals]=estimateLD(lanczosDim, hamiltonian, iniState, err);
   std::cout<< "lancos vectors vas "<< lanczosDim<<std::endl;
     //    Q=Q*S;
    
    std::vector<double> obs(observable.size(), 0);
   		 
  LTiterate<false>(lanczosDim,  beta, obs, observable, eigenVals,  Z, S,  iniState,  hamiltonian);
		      
   			      //		
   		    
		  return {obs, real(Z)};
}
template<typename Matrix, typename Container>
auto  calculate_lanczLT_slow(Matrix& hamiltonian, Container& observable, std::vector<double> beta, int lanczosDim=20, double err=1E-10)->std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
{
   std::random_device rd;
std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
std::uniform_real_distribution<double> dis(-1.0, 1.0);

 
 Eigen::VectorXcd iniState=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double>(dis(gen), dis(gen));});

    
   iniState/=iniState.norm();
   //Eigen::MatrixXcd Q(hamiltonian.rows(), lanczosDim);
   auto [S, eigenVals]=estimateLD(lanczosDim, hamiltonian, iniState, err);
   std::cout<< "lancos vectors vas "<< lanczosDim<<std::endl;
  
   Eigen::MatrixXd obs=Eigen::MatrixXd::Zero( beta.size(),observable.size());
   Eigen::VectorXd z=Eigen::VectorXd::Zero(beta.size());

   	       LTiterate<false>(lanczosDim, beta, obs, observable, eigenVals,  z, S, iniState, hamiltonian);

   		    
		  return {obs, z};
}
template<typename Matrix, typename Container>
auto  calculate_lanczLT_fast(Matrix& hamiltonian, Container& observable, double beta, int lanczosDim=20, double err=1E-10)->std::tuple<std::vector<double>, double>
{
   std::random_device rd;
std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
std::uniform_real_distribution<float> dis(-1.0, 1.0);

 
 Eigen::VectorXcd iniState=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double>(dis(gen), dis(gen));});
  auto Z=std::complex<double>{0, 0};
 
    auto A=std::complex<double>{0, 0};
      

   iniState/=iniState.norm();
   //Eigen::MatrixXcd Q(hamiltonian.rows(), lanczosDim);
   auto [Q, S, eigenVals]=estimateLDwithQ(lanczosDim, hamiltonian, iniState, err);
   std::cout<< "lancos vectors vas "<< lanczosDim<<std::endl;
         Q=Q*S;
    
    std::vector<double> obs(observable.size(), 0);
   		 
  LTiterate<true>(lanczosDim,  beta, obs, observable, eigenVals, Z, Q,  iniState,  hamiltonian);
		      
   			      //		
   		    
		  return {obs, real(Z)};
}
template<typename Matrix, typename Container>
auto  calculate_lanczLT_fast(Matrix& hamiltonian, Container& observable, std::vector<double> beta, int lanczosDim=20, double err=1E-10)->std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
{
   std::random_device rd;
std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
std::uniform_real_distribution<double> dis(-1.0, 1.0);

 
 Eigen::VectorXcd iniState=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double>(dis(gen), dis(gen));});

    
   iniState/=iniState.norm();
  
   auto [Q, S, eigenVals]=estimateLDwithQ(lanczosDim, hamiltonian, iniState, err);
   std::cout<< "lancos vectors vas "<< lanczosDim<<std::endl;
         Q=Q*S;
   Eigen::MatrixXd obs=Eigen::MatrixXd::Zero( beta.size(),observable.size());
   Eigen::VectorXd z=Eigen::VectorXd::Zero(beta.size());

    
   	       LTiterate<true>(lanczosDim, beta, obs, observable, eigenVals,  z, Q, iniState, hamiltonian);

   		    
		  return {obs, z};
}
template<typename Matrix, typename Container>
std::vector<double> LTLM(Matrix& hamiltonian, Container& observable,  double T, int runs,     size_t lanczosDim=20)
{

   double beta=1./T;
   std::vector<double> obstot(observable.size(), 0);
       auto Ztot{ 0.0};
          for(int r=0; r<runs; r++)
	    {
        auto tup=calculate_lanczLT(hamiltonian, observable, beta, lanczosDim);
  auto obs=std::get<0>(tup);
         auto Z=std::get<1>(tup);
	 //	std::transform(obstot.begin(), obstot.end(), obs.begin(), obs.end(), std::plus<double>());
	 for(int k=0; k<obs.size(); k++)
	   {
	     obstot[k]+=obs[k];
	   }

        Ztot+=Z;  
     }
   std::for_each(obstot.begin(), obstot.end(),[&Ztot](double &val) {val/=Ztot;});
  return obstot;
     }

#endif /* LTLANCZOSIMPL_H */
