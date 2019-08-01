 #ifndef FTLANCZOSIMPL_H
#define FTLANCZOSIMPL_H
#include"LanczosHelper.hpp"
using VecD=Eigen::VectorXd;
using VecC=Eigen::VectorXcd;
using MatD=Eigen::MatrixXd;
template<bool Fast, typename Vector, typename MatrixHam, typename MatrixTrafo, typename OBS,typename Container>
void FTiterate(int lanczosDim, int lanczosDim2, double beta, OBS&& obs, Container& observable, VecD& eigenVals, VecD& eigenVals2, double& Z, MatrixTrafo& T1, MatrixTrafo& T2, Vector& iniState1, Vector& iniState2, MatrixHam& hamiltonian){
  for(int j=0; j<std::min(lanczosDim,lanczosDim2) ; j++)
   		    {
		      double exponent=-(beta*(eigenVals[j]-eigenVals[0]));
		       auto exp=std::exp(exponent);
		       VecC qvec1;
		       VecC qvec2;
		       if constexpr(!Fast)
				     {
				       qvec1=Many_Body::lanczTrafo(T1.col(j), iniState1, lanczosDim, hamiltonian);
				       qvec2=Many_Body::lanczTrafo(T2.col(j), iniState2, lanczosDim2, hamiltonian);
				     }
		       else{
			 qvec1=T1.col(j);
			 qvec2=T2.col(j);
		       }
		       //		       auto link=iniState.adjoint()*Q.col(j);

		       // should be able to vectorize thi
		       
		       //		       auto link=iniState.adjoint()*Q.col(j);
		       auto link1=iniState1.adjoint()*qvec1;
		       auto link2=iniState2.adjoint()*qvec2;
		       Z+=std::abs(link2(0, 0))*std::abs(link2(0, 0))*exp;
		       for(int i=0; i<obs.size(); i++)
		       	{
			 auto expval=qvec1.adjoint()*observable[i]*iniState1;
		       obs[i]+=real((link1(0, 0)*expval(0, 0)))*exp;
			}
   		    }
  if(lanczosDim!=lanczosDim2){
for(int j=std::min(lanczosDim,lanczosDim2); j<std::max(lanczosDim,lanczosDim2); j++)
   		    {
		          double exponent=-(beta*(eigenVals[j]-eigenVals[0]));
		       auto exp=std::exp(exponent);
		      if(lanczosDim<lanczosDim2)
			{
			  VecC qvec1;
			  
		       if constexpr(!Fast)
				     {
				       qvec1=Many_Body::lanczTrafo(T1.col(j), iniState1, lanczosDim, hamiltonian);
			
				     }
		       else{
			 qvec1=T1.col(j);
			
		       }
		       auto link1=iniState1.adjoint()*qvec1;
		       
		       for(int i=0; i<obs.size(); i++)
		       	{
			 auto expval=qvec1.adjoint()*observable[i]*iniState1;
		       obs[i]+=real((link1(0, 0)*expval(0, 0)))*exp;
			}
			}
		      else{
			 VecC qvec2;
			   if constexpr(!Fast)
				     {
				       qvec2=Many_Body::lanczTrafo(T2.col(j), iniState2, lanczosDim2, hamiltonian);
				     }
		       else{
			 qvec2=T2.col(j);
		       }
		       
		       //		       auto link=iniState.adjoint()*Q.col(j);
		       auto link2=iniState2.adjoint()*qvec2;
		       Z+=std::abs(link2(0, 0))*std::abs(link2(0, 0))*exp;
		     }

   		    }
  }
  
}
template<typename Matrix, typename Container>
auto  calculate_lanczFT_slow(Matrix& hamiltonian, Container& observable, double beta, int lanczosDim=20, double err=1E-10)->std::tuple<std::vector<double>, double>
{
   std::random_device rd;
std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
 std::uniform_real_distribution<double> dis(-1.0, 1.0);
auto Z=std::complex<double>{0, 0};
    auto A=std::complex<double>{0, 0};
 
 Eigen::VectorXcd iniState=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double>(dis(gen), dis(gen));});

iniState/=iniState.norm();

 Eigen::VectorXcd iniState2=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double> (dis(gen), dis(gen));});
   iniState2/=iniState2.norm();
   size_t lanczosDim2=lanczosDim;
   auto [S, eigenVals]=estimateLD(lanczosDim, hamiltonian, iniState, err);
   auto [S2, eigenVals2]=estimateLD(lanczosDim2, hamiltonian, iniState2, err);    
 std::cout<< "lancos vectors vas "<< lanczosDim << " and " << lanczosDim2<<std::endl;
    std::vector<double> obs(observable.size(), 0);
     FTiterate<false>(lanczosDim, lanczosDim2, beta, obs, observable, eigenVals, eigenVals2, Z, S, S2, iniState, iniState2, hamiltonian);
   		 return {obs, real(Z)};
}

template<typename Matrix, typename Container>
auto  calculate_lanczFT_slow(Matrix& hamiltonian, Container& observable, std::vector<double> beta, int lanczosDim=20, double err=1E-10)->std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
{
    std::random_device rd;
std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
  std::uniform_real_distribution<double> dis(-1.0, 1.0);

 
  Eigen::VectorXcd iniState=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double>(dis(gen), dis(gen));});

 iniState/=iniState.norm();

  Eigen::VectorXcd iniState2=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double> (dis(gen), dis(gen));});
    iniState2/=iniState2.norm();
   int lanczosDim2=lanczosDim;
    auto [S, eigenVals]=estimateLD(lanczosDim, hamiltonian, iniState, err);
    auto [S2, eigenVals2]=estimateLD(lanczosDim2, hamiltonian, iniState2, err);
  std::cout<< "lancos vectors vas "<< lanczosDim << " and " << lanczosDim2<<std::endl;
    Eigen::MatrixXd obs=Eigen::MatrixXd::Zero( beta.size(),observable.size());
    Eigen::VectorXd z=Eigen::VectorXd::Zero(beta.size());
   for(size_t i=0; i<beta.size(); i++)
     {
       FTiterate<false>(lanczosDim, lanczosDim2, beta[i], obs.row(i), observable, eigenVals, eigenVals2, z(i), S, S2, iniState, iniState2, hamiltonian);

    }
		  return {obs, z};
}
template<typename Matrix, typename Container>
auto  calculate_lanczFT_fast(Matrix& hamiltonian, Container& observable, double beta, int lanczosDim=20, double err=1E-10)->std::tuple<std::vector<double>, double>
{
   std::random_device rd;
std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
 std::uniform_real_distribution<double> dis(-1.0, 1.0);
auto Z=std::complex<double>{0, 0};
    auto A=std::complex<double>{0, 0};
 
 Eigen::VectorXcd iniState=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double>(dis(gen), dis(gen));});

iniState/=iniState.norm();

 Eigen::VectorXcd iniState2=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double> (dis(gen), dis(gen));});
   iniState2/=iniState2.norm();
   size_t lanczosDim2=lanczosDim;
   auto [Q1, S, eigenVals]=estimateLDwithQ(lanczosDim, hamiltonian, iniState, err);
   auto [Q2, S2, eigenVals2]=estimateLDwithQ(lanczosDim2, hamiltonian, iniState2, err);
    Q1=Q1*S;
   Q2=Q2*S2;
 std::cout<< "lancos vectors vas "<< lanczosDim << " and " << lanczosDim2<<std::endl;
    std::vector<double> obs(observable.size(), 0);
    FTiterate<true>(lanczosDim, lanczosDim2, beta, obs, observable, eigenVals, eigenVals2, Z, Q1, Q2, iniState, iniState2, hamiltonian);
   		 return {obs, real(Z)};
}

template<typename Matrix, typename Container>
auto  calculate_lanczFT_fast(Matrix& hamiltonian, Container& observable, std::vector<double> beta, int lanczosDim=20, double err=1E-10)->std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
{
    std::random_device rd;
std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
  std::uniform_real_distribution<double> dis(-1.0, 1.0);

 
  Eigen::VectorXcd iniState=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double>(dis(gen), dis(gen));});

 iniState/=iniState.norm();

  Eigen::VectorXcd iniState2=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return std::complex<double> (dis(gen), dis(gen));});
    iniState2/=iniState2.norm();
   int lanczosDim2=lanczosDim;
   auto [Q1, S, eigenVals]=estimateLDwithQ(lanczosDim, hamiltonian, iniState, err);
   auto [Q2, S2, eigenVals2]=estimateLDwithQ(lanczosDim2, hamiltonian, iniState2, err);
   Q1=Q1*S;
   Q2=Q2*S2;
  std::cout<< "lancos vectors vas "<< lanczosDim << " and " << lanczosDim2<<std::endl;
    Eigen::MatrixXd obs=Eigen::MatrixXd::Zero( beta.size(),observable.size());
    Eigen::VectorXd z=Eigen::VectorXd::Zero(beta.size());
   for(size_t i=0; i<beta.size(); i++)
     {
       FTiterate<true>(lanczosDim, lanczosDim2, beta[i], obs.row(i), observable, eigenVals, eigenVals2, z(i), Q1, Q2, iniState, iniState2, hamiltonian);

    }
		  return {obs, z};
}
template<typename Matrix, typename Container>
std::vector<double> FTLM(Matrix& hamiltonian, Container& observable, double T, int runs,     size_t lanczosDim=20)
{
 

   double beta=1./T;
 std::vector<double> obstot(observable.size(), 0);
 

       auto Ztot{ 0.0};

   for(int r=0; r<runs; r++)
     {
        auto tup=calculate_lanczFT_fast(hamiltonian, observable, beta, lanczosDim);
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


#endif /* FTLANCZOSIMPL_H */
