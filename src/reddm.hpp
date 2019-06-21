 #ifndef REDDM_H
#define REDDM_H
#include"basis.hpp"
#include <Eigen/Dense>
#include"diag.h"
namespace Many_Body{
  template<typename T>
  using MatrixD =Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  template<typename T>
  using VectorD =Eigen::Matrix<T,  Eigen::Dynamic, 1>;
  template<typename T>
  inline T conj(T a)
  {return a;}
  
    template<typename T>
  inline std::complex< T > conj(std::complex< T > a)
  { std::complex< T > b(real(a), -imag(a));
    return b;}
  
  template<typename TotalBasis, typename Matrix, typename Vector>
  std::vector<Eigen::VectorXd> makeThermalRDMTP(Matrix& H, Vector& ev, TotalBasis& TP, double T, bool isDiag=false, int site=0)
{
  double beta=1./T;
  if(!isDiag)
    {
  diagMat(H, ev);
    }
  std::cout<< "GS "<< ev[0]<<std::endl;
  Eigen::MatrixXd rho=Eigen::MatrixXd::Zero(H.rows(), H.rows());
  double Z{0};
  //  std::cout<< ev << std::endl;
  for(int i=0; i<rho.rows(); i++)
    {
      rho(i, i)=std::exp(-beta*ev(i));
      Z+=std::exp(-beta*ev(i));
    }
  rho/=Z;
  //  std::cout<< rho<< std::endl;
  
  rho=rho*H.adjoint();
  rho=H*rho;
  std::vector<Eigen::VectorXd> evs;
  auto RDM= makeRedDMTP(TP, TP.rbasis, site, rho);
  for(auto& i : RDM){

  Eigen::VectorXd evRDM(i.rows());
  diagMat(i, evRDM);
  evs.push_back(evRDM);
  
  }
  
  return evs;
}
template<typename Basis, typename T>
void makeRedDM(Basis& basis, int site, MatrixD<T>& rho )
{
  size_t dim=basis.maxParticles+1;
  MatrixD<T> DM=MatrixD<T>::Zero(dim, dim);
  for(auto it= basis.begin(); it!=basis.end(); ++it)
    {
      //      std::cout<< " loop 1 "<<std::endl;
      auto currentID1 =Id(*it);
      auto p1=basis.particlesAt(currentID1, site);
      DM(p1, p1)+=rho(Position(*it), Position(*it));
      auto itx=it;
      itx++;
      for(auto it2= itx; it2!=basis.end(); ++it2)
    {
      //            std::cout<< " loop 2 "<<std::endl;
      
      auto currentID2 =Id(*it2);
       auto state1=GetLattice(*it).makeStateVec();
       auto state2=GetLattice(*it2).makeStateVec();
       auto p2=basis.particlesAt(currentID2, site);
       state1.erase(state1.begin()+site);
       state2.erase(state2.begin()+site);
       if(state1==state2)
	 {
      DM(p1, p2)+=rho(Position(*it), Position(*it2));
      DM(p2, p1)+=rho(Position(*it2), Position(*it));
      // std::cout << currentID1 << " part a 1 " << p1 << std::endl;
      //       std::cout << currentID2 << " part a 2 " << p2 << std::endl;
	 }
    }
    }
  
  std::cout<< DM << std::endl;
}
  template<typename Basis, typename T>
  void makeRedDM(Basis& basis, int site, VectorD<T>& psi)
{
  size_t dim=basis.maxParticles+1;
  
  MatrixD<T> DM=MatrixD<T>::Zero(dim, dim);
  for(auto it= basis.begin(); it!=basis.end(); ++it)
    {
      //      std::cout<< " loop 1 "<<std::endl;
      auto currentID1 =Id(*it);
      auto p1=basis.particlesAt(currentID1, site);
      DM(p1, p1)+=Many_Body::conj(static_cast<T>((psi(Position(*it)))))*(psi(Position(*it)));
	      //rho(Position(*it), Position(*it));
      auto itx=it;
      itx++;
      for(auto it2= itx; it2!=basis.end(); ++it2)
    {
      //            std::cout<< " loop 2 "<<std::endl;
      
      auto currentID2 =Id(*it2);
       auto state1=GetLattice(*it).makeStateVec();
       auto state2=GetLattice(*it2).makeStateVec();
       auto p2=basis.particlesAt(currentID2, site);
       state1.erase(state1.begin()+site);
       state2.erase(state2.begin()+site);
       if(state1==state2)
	 {
	   DM(p1, p2)+=Many_Body::conj(static_cast<T>((psi(Position(*it2)))))*(psi(Position(*it)));
	//rho(Position(*it), Position(*it2));

	   DM(p2, p1)+=Many_Body::conj(static_cast<T>((psi(Position(*it)))))*(psi(Position(*it2)));
	//rho(Position(*it2), Position(*it));
      // std::cout << currentID1 << " part a 1 " << p1 << std::endl;
      //       std::cout << currentID2 << " part a 2 " << p2 << std::endl;
	 }
    }
    }
  
  std::cout<< DM << std::endl;
}

template<typename TotalBasis, class SubBasis, typename T>
std::vector<MatrixD<T>> makeRedDMTP(const TotalBasis& totalBasis, const SubBasis& subBasis, int site, MatrixD<T>& rho )
{
  // making reduced DM for each degree of freedom in subasis L   
       using SubBasisIt= typename SubBasis::BasisIt;
       using TpBasisIt= typename TotalBasis::BasisIt;
       using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
       using RightBasisIt= typename TotalBasis::RightBasisIt;
       using LeftBasis= typename TotalBasis::LB;     
       using RightBasis= typename TotalBasis::RB;  
       using Lattice=typename SubBasis::Lattice;
       std::vector<MatrixD<T>> mats(totalBasis.lbasis.maxParticles+1);
       size_t dim=(totalBasis.rbasis.maxParticles+1);

       std::fill(mats.begin(), mats.end(), MatrixD<T>::Zero(dim, dim));

       
       //  std::cout<< dim << std::endl;

  for(auto it= totalBasis.begin(); it!=totalBasis.end(); ++it)
    {
      auto itL=totalBasis.lbasis.find(LeftId(*it));
      auto itR=totalBasis.rbasis.find(RightId(*it));
      
      //      std::cout<< " loop 1 "<<std::endl;
       auto currentID1L =Id(*itL);
       auto currentID1R =Id(*itR);
       auto p1L=totalBasis.lbasis.particlesAt(currentID1L, site);
       auto p1R=totalBasis.rbasis.particlesAt(currentID1R, site);
       // unsafe way to determine the position in the matrix
       // upper first fix electron then iterate boson
       size_t newPos= p1R;
       mats[p1L](newPos, newPos)+=rho(Position(*it), Position(*it));
       //       if(std::abs(rho(Position(*it), Position(*it)))>0.001){std::cout<< " her "<< p1L<< "  " << p1R << std::endl;}
       auto itx=it;
       itx++;
       for(auto it2= itx; it2!=totalBasis.end(); ++it2)
     {
         auto it2L=totalBasis.lbasis.find(LeftId(*it2));
      auto it2R=totalBasis.rbasis.find(RightId(*it2));
      
      //      std::cout<< " loop 1 "<<std::endl; // re think this one
       auto currentID2L =Id(*it2L);
       auto currentID2R =Id(*it2R);
       auto p2L=totalBasis.lbasis.particlesAt(currentID2L, site);
       auto p2R=totalBasis.rbasis.particlesAt(currentID2R, site);
       size_t newPos2=p2R;
       // making the delta function
       auto state1L=GetLattice(*itL).makeStateVec();
       auto state1R=GetLattice(*itR).makeStateVec();
       auto state2L=GetLattice(*it2L).makeStateVec();
       auto state2R=GetLattice(*it2R).makeStateVec();
       state1L.erase(state1L.begin()+site);
       state2L.erase(state2L.begin()+site);
       state1R.erase(state1R.begin()+site);
       state2R.erase(state2R.begin()+site);
       if(state1L==state2L and state1R==state2R)
	 {
	   if(p2L!=p1L){std::cout << "here "<< std::endl;}
       mats[p2L](newPos, newPos2)+=rho(Position(*it), Position(*it2));
       mats[p2L](newPos2, newPos)+=rho(Position(*it2), Position(*it));
	 }
         }
}
  return mats;
}

  template<typename TotalBasis, class SubBasis, typename T>
std::vector<MatrixD<T>> makeRedDMTP(const TotalBasis& totalBasis, const SubBasis& subBasis, int site, VectorD<T>& psi )
{
  // making reduced DM for each degree of freedom in subasis L   
       using SubBasisIt= typename SubBasis::BasisIt;
       using TpBasisIt= typename TotalBasis::BasisIt;
       using LeftBasisIt= typename TotalBasis::LeftBasisIt;     
       using RightBasisIt= typename TotalBasis::RightBasisIt;
       using LeftBasis= typename TotalBasis::LB;     
       using RightBasis= typename TotalBasis::RB;  
       using Lattice=typename SubBasis::Lattice;
       std::vector<MatrixD<T>> mats(totalBasis.lbasis.maxParticles+1);
       size_t dim=(totalBasis.rbasis.maxParticles+1);

       std::fill(mats.begin(), mats.end(), MatrixD<T>::Zero(dim, dim));

       
       //  std::cout<< dim << std::endl;

  for(auto it= totalBasis.begin(); it!=totalBasis.end(); ++it)
    {
      auto itL=totalBasis.lbasis.find(LeftId(*it));
      auto itR=totalBasis.rbasis.find(RightId(*it));
      
      //      std::cout<< " loop 1 "<<std::endl;
       auto currentID1L =Id(*itL);
       auto currentID1R =Id(*itR);
       auto p1L=totalBasis.lbasis.particlesAt(currentID1L, site);
       auto p1R=totalBasis.rbasis.particlesAt(currentID1R, site);
       // unsafe way to determine the position in the matrix
       // upper first fix electron then iterate boson
       size_t newPos= p1R;
       mats[p1L](newPos, newPos)+=Many_Body::conj(static_cast<T>((psi(Position(*it)))))*(psi(Position(*it)));

       //       if(std::abs(rho(Position(*it), Position(*it)))>0.001){std::cout<< " her "<< p1L<< "  " << p1R << std::endl;}
       auto itx=it;
       itx++;
       for(auto it2= itx; it2!=totalBasis.end(); ++it2)
     {
         auto it2L=totalBasis.lbasis.find(LeftId(*it2));
      auto it2R=totalBasis.rbasis.find(RightId(*it2));
      
      //      std::cout<< " loop 1 "<<std::endl;
       auto currentID2L =Id(*it2L);
       auto currentID2R =Id(*it2R);
       auto p2L=totalBasis.lbasis.particlesAt(currentID2L, site);
       auto p2R=totalBasis.rbasis.particlesAt(currentID2R, site);
       size_t newPos2=p2R;
       // making the delta function
       auto state1L=GetLattice(*itL).makeStateVec();
       auto state1R=GetLattice(*itR).makeStateVec();
       auto state2L=GetLattice(*it2L).makeStateVec();
       auto state2R=GetLattice(*it2R).makeStateVec();
       state1L.erase(state1L.begin()+site);
       state2L.erase(state2L.begin()+site);
       state1R.erase(state1R.begin()+site);
       state2R.erase(state2R.begin()+site);
       if(state1L==state2L and state1R==state2R)
	 {
	   mats[p2L](newPos, newPos2)+=Many_Body::conj(static_cast<T>((psi(Position(*it2)))))*(psi(Position(*it)));

	 //rho(Position(*it), Position(*it2));
	   mats[p2L](newPos2, newPos)+=Many_Body::conj(static_cast<T>((psi(Position(*it)))))*(psi(Position(*it2)));
	 //rho(Position(*it2), Position(*it));
       	   


	 }
         }
}
  return mats;
}
  
}


#endif /* REDDM_H */
