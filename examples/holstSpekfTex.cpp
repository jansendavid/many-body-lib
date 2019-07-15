#include<iostream>
#define EIGEN_USE_MKL_ALL
#include"basis.hpp"
#include"operators.hpp"
#include"diag.h"
#include"tpoperators.hpp"
#include"files.hpp"
#include"timeev.hpp"
#include"ETH.hpp"
#include<iomanip>
#include <boost/program_options.hpp>

using namespace boost::program_options;
int main(int argc, char *argv[])
{
  using namespace Eigen;
using namespace std;
using namespace Many_Body;
using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
   // std::vector<size_t> ee(L, 0);
 using Mat= Operators::Mat;
  size_t M{};
  size_t L{};
  double t0{};
  double omega{};
  double gamma{};
  double T{};
  bool PB{};
  try
  {
    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("L", value(&L)->default_value(4), "L")
      ("M", value(&M)->default_value(2), "M")
      ("t", value(&t0)->default_value(1.), "t0")
      ("gam", value(&gamma)->default_value(1.), "gamma")
      ("omg", value(&omega)->default_value(1.), "omega")
      ("T", value(&T)->default_value(1.), "T")
    ("pb", value(&PB)->default_value(true), "PB");
  


    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help"))
      {std::cout << desc << '\n'; return 0;}
    else{
      if (vm.count("L"))
      {      std::cout << "L: " << vm["L"].as<size_t>() << '\n';
	
      }
     if (vm.count("M,m"))
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
         if (vm.count("T"))
      {
	std::cout << "T: " << vm["T"].as<double>() << '\n';
      }
              if (vm.count("pb"))
      {
	std::cout << "PB: " << vm["pb"].as<bool>() << '\n';
      }
    }
  }
  catch (const error &ex)
  {
    std::cerr << ex.what() << '\n';
    return 0;
  }
   ElectronBasis e( L, 1);
   //   std::cout<< e<<std::endl;
      ElectronState e2( L, 0);
   //   std::cout<< e<<std::endl;
  
  PhononBasis ph(L, M);
  //  std::cout<< ph<<std::endl;
  HolsteinBasis TP(e, ph);
  //HolsteinBasis TP2(e2, ph);
  //  std::cout<< TP<<std::endl;
  e.insert(e2);
  std::cout<< e<< std::endl;
  //  std::cout<< TP2<< std::endl;
  
  std::cout<<"total dim "<< TP.dim << std::endl;
  std::cout<<std::endl;
        Mat E1=Operators::EKinOperatorL(TP, e, t0, PB);
      Mat Ebdag=Operators::NBosonCOperator(TP, ph, gamma, PB);
      Mat Eb=Operators::NBosonDOperator(TP, ph, gamma, PB);
      Mat Eph=Operators::NumberOperator(TP, ph, omega,  PB);



      Mat H=E1  +Ebdag + Eb+ Eph;
          Eigen::MatrixXcd HH=Eigen::MatrixXcd(H);


	     Eigen::VectorXd ev=Eigen::VectorXd(TP.dim);
	     diagMat(HH, ev);
	     return 0;
}
 
