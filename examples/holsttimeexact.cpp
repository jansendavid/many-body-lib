#include<iostream>
#include"basis.hpp"
#include"operators.hpp"
#include"diag.h"
#include"tpoperators.hpp"
#include"files.hpp"
#include"timeev.hpp"
#include"ETH.hpp"
#include <boost/program_options.hpp>
int main(int argc, char *argv[])
{
  using namespace Eigen;
using namespace std;
using namespace Many_Body;
  using Fermi_HubbardBasis= TensorProduct<ElectronBasis, ElectronBasis>;
   using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
  using Mat= Operators::Mat;
  using boost::program_options::value;
   size_t M{};
  size_t L{};
  double t0{};
  double omega{};
  double gamma{};
  double dt{};
  double tot{};
  bool PB{};
  try
  {
    boost::program_options::options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("L", value(&L)->default_value(4), "L")
      ("M,m", value(&M)->default_value(2), "M")
      ("t", value(&t0)->default_value(1.), "t0")
      ("gam", value(&gamma)->default_value(1.), "gamma")
      ("omg", value(&omega)->default_value(1.), "omega")
    ("dt", value(&dt)->default_value(0.1), "dt")
      ("tot", value(&tot)->default_value(1.), "tot")
    ("pb", value(&PB)->default_value(true), "PB");
  


    boost::program_options::variables_map vm;
    boost::program_options::store(parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    if (vm.count("help"))
      {std::cout << desc << '\n'; return 0;}
    else{
      if (vm.count("L"))
      {      std::cout << "L: " << vm["L"].as<size_t>() << '\n';
	
      }
     if (vm.count("M"))
      {
	std::cout << "M: " << vm["M"].as<size_t>() << '\n';
	
      }
      if (vm.count("t"))
      {
	std::cout << "t0: " << vm["t"].as<double>() << '\n';	
      }
       if (vm.count("omg"))
      {
	std::cout << "omega: " << vm["omg"].as<double>() << '\n';
      }
       if (vm.count("gam"))
      {
	std::cout << "gamma: " << vm["gam"].as<double>() << '\n';
      }
       if (vm.count("dt"))
      {
	std::cout << "dt: " << vm["dt"].as<double>() << '\n';
      }
       if (vm.count("tot"))
      {
	std::cout << "total time: " << vm["tot"].as<double>() << '\n';
      }
       if (vm.count("pb"))
      {
	std::cout << "PB: " << vm["pb"].as<bool>() << '\n';
      }
    }
  }
  catch (const boost::program_options::error &ex)
  {
    std::cerr << ex.what() << '\n';
    return 0;
  }

  //std::vector<int> ee(L, 0);
  // ee[L-1]=1;

  ElectronState estate1(L, 1);
  ElectronBasis e( L, 1);
  
  //std::cout<< e<<std::endl;
std::vector<size_t> es(L, 0);
      es[0]=1;
      ElectronState aa(es);
  PhononBasis ph(L, M);
  //               std::cout<< ph<<std::endl;
      HolsteinBasis TP(e, ph);
      std::vector<size_t> state(L, 0);
  std::fill(state.begin(), state.end(), 1);
  BosonState b2(state, M);
 
     auto it3=TP.lbasis.find(aa.GetId());
     auto it4=TP.rbasis.find(b2.GetId());
     size_t StateNr= Position(*it4)*TP.lbasis.dim +Position(*it3);
     //        std::cout<< TP << std::endl;
     //    std::cout<< b2.GetId()<<std::endl;
     // 	std::cout<< b2.GetId()<<std::endl;
     // std::cout<<"state nr "<< StateNr<<std::endl;                
         Eigen::VectorXcd inistate(TP.dim);

 	 Mat O=Operators::NumberOperator(TP, ph, 1,  PB);
   inistate.setZero();

   inistate[StateNr]=1;
   	          Eigen::VectorXcd i0=inistate;
      //   std::cout<< e << std::endl;
		  
           std::cout<< TP.dim << std::endl;     
      Mat E1=Operators::EKinOperatorL(TP, e, t0, PB);
      Mat Ebdag=Operators::NBosonCOperator(TP, ph, gamma, PB);
      Mat Eb=Operators::NBosonDOperator(TP, ph, gamma, PB);
      Mat Eph=Operators::NumberOperator(TP, ph, omega,  PB);
      
      //Mat E=Operators::NumberOperatore(TP, e, 1, false);
      //    std::cout<< HH << std::xbendl;
      Eigen::VectorXd eigenVals(TP.dim);
      Mat H=E1+Eph  +Ebdag + Eb;
      //    
 
    Eigen::MatrixXd HH=Eigen::MatrixXd(H);
    //        std::cout<< HH<<std::endl;
        Eigen::MatrixXd N=Eigen::MatrixXd(Eph);
     Many_Body::diagMat(HH, eigenVals);
         Eigen::VectorXd energy(TP.dim);
    	   std::vector<double> obs(TP.dim);
	   std::vector<double> obs2(TP.dim);
	   std::vector<double> ensvec(TP.dim);
	   Eigen::MatrixXd PHD2=HH.adjoint()*N.selfadjointView<Lower>()*HH;

     Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals, dt);
   Eigen::MatrixXcd cEVec=HH.cast<std::complex<double>>();
    Eigen::VectorXcd newIn=inistate;
   int i=0;

   // O=H;
        while(i*dt<tot)
       {

       	 TimeEv::timeev_exact(newIn, cEVec, evExp);
  	 
   	  	 std::complex<double> c=(newIn.adjoint()*(O*newIn))(0);
		 //	 std::complex<double> c=(inistate.adjoint()*(newIn))(0);
  	 
   	 i++;

   			// 	outputVals(i)=real(c);
   			// 	outputVals2(i)=real(c2);
        		// outputTime(i)=i*dt;
   	 std::cout<< std::setprecision(8)<<real(c)<< " dt "<< i*dt<< std::endl;
	 //	  	 BOOST_CHECK(std::abs(real(c2)-real(c))<Many_Body::err);
	     	             }

 //std::cout<< MatrixXd(OBS2) << std::endl;
 // Eigen::MatrixXd OBS22=HH.adjoint()*(OBS2.selfadjointView<Lower>())*HH;
 // Eigen::VectorXd diagobs2=OBS22.diagonal();
	// std::cout<< en(0) << std::endl;
 // std::cout<<diagobs2(0)<<std::endl; 
    std::string filename=".bin";
    //     bin_write("E"+filename, en);
  return 0;
}
 
