 #ifndef FTLANCZOS_H
#define FTLANCZOS_H
#include"basis.hpp"
#include"lanzcos.hpp"
#include"diag.h"
#include<random>
// #include<boost/mpi.hpp>
// #include<boost/mpi/communicator.hpp>
// #include<boost/mpi/environment.hpp>
//namespace mpi=boost::mpi;
template<typename Matrix, typename Container>
auto  calculate_lanczFT(Matrix& hamiltonian, Container& observable, double beta, size_t lanczosDim=20)->std::tuple<std::vector<double>, double>
{
   std::random_device rd;
std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
std::uniform_real_distribution<float> dis(-1.0, 1.0);

 
   Eigen::VectorXcd iniState=Eigen::VectorXcd::NullaryExpr(hamiltonian.rows(),[&](){return dis(gen);});
   iniState/=iniState.norm();
   auto Z=std::complex<double>{0, 0};
    auto A=std::complex<double>{0, 0};

    Many_Body::TriDiagMat tri=Many_Body::Lanczos(hamiltonian, iniState, lanczosDim);
    Eigen::MatrixXd S(lanczosDim, lanczosDim);
    Eigen::VectorXd eigenVals(lanczosDim);
    
     Many_Body::diag(tri, S, eigenVals);
    std::vector<double> obs(observable.size(), 0);
   		  for(int j=0; j<lanczosDim; j++)
   		    {
		       auto exp=std::exp(-beta*(eigenVals[j]));
		       auto qvec=Many_Body::lanczTrafo(S.col(j), iniState, lanczosDim, hamiltonian);
		       //		       auto link=iniState.adjoint()*Q.col(j);
		       auto link=iniState.adjoint()*qvec;
		       Z+=std::abs(link(0, 0))*std::abs(link(0, 0))*exp;
		      
		       for(int i=0; i<obs.size(); i++)
		       	{
			  
   		       auto expval=qvec.adjoint()*observable[i]*iniState;
		       obs[i]+=real((link(0, 0)*expval(0, 0)))*exp;
			}

		      
   			      //		
   		    }
		  return {obs, real(Z)};
}
template<typename Matrix, typename Container>
auto  calculate_lanczLT(Matrix& hamiltonian, Container& observable, double beta, size_t lanczosDim=20)->std::tuple<std::vector<double>, double>
{
  auto Z=std::complex<double>{0, 0};
 
    auto A=std::complex<double>{0, 0};
       Eigen::VectorXcd iniState=Eigen::VectorXcd::Random(hamiltonian.rows());

   iniState/=iniState.norm();
   //Eigen::MatrixXcd Q(hamiltonian.rows(), lanczosDim);


    Many_Body::TriDiagMat tri=Many_Body::Lanczos(hamiltonian, iniState, lanczosDim);
    Eigen::MatrixXd S(lanczosDim, lanczosDim);
    Eigen::VectorXd eigenVals(lanczosDim);
    
     Many_Body::diag(tri, S, eigenVals);

     //    Q=Q*S;
    
    std::vector<double> obs(observable.size(), 0);
   		  for(int j=0; j<lanczosDim; j++)
   		    {
		      for(int l=0; l<lanczosDim; l++)
		    {
		      auto exp=std::exp(-beta*(eigenVals[j]+eigenVals[l])/2);
		      
		      // auto link=iniState.adjoint()*Q.col(j);
		      auto qvecj=Many_Body::lanczTrafo(S.col(j), iniState, lanczosDim, hamiltonian);
		      auto link=iniState.adjoint()*qvecj;
		      auto qvecl=Many_Body::lanczTrafo(S.col(l), iniState, lanczosDim, hamiltonian);
		      //		       auto right=(Q.col(l)).adjoint()*iniState;
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

		      
   			      //		
   		    
		  return {obs, real(Z)};
}
template<typename Matrix, typename Container>
std::vector<double> FTLM(Matrix& hamiltonian, Container& observable, double T, int runs,     size_t lanczosDim=20)
{
 

   double beta=1./T;
 std::vector<double> obstot(observable.size(), 0);
 

       auto Ztot{ 0.0};

   for(int r=0; r<runs; r++)
     {
        auto tup=calculate_lanczFT(hamiltonian, observable, beta, lanczosDim);
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

// void paraLM(){
//   mpi::environment env;
//   mpi::communicator world;
//   std::cout<< " i am process # "<< world.rank() << std::endl;
// };

//  size_t runs=10;
// // make initial vector
//  auto iniState=Eigen::VectorXcd::Random(hamiltonian.rows());
//     Eigen::MatrixXcd Q(hamiltonian.rows(), lanczosDim);
// // determine a better way to deccide on Lanczos dim
// Many_Body::TriDiagMat tri1=Many_Body::Lanczos(hamiltonian, iniState, lanczosDim, Q);



#endif /* FTLANCZOS_H */
