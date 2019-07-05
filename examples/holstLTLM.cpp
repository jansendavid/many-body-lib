#include<iostream>
#define EIGEN_USE_MKL_ALL
#include <Eigen/Eigenvalues> 
#include"numerics.hpp"
#include"reddm.hpp"
#include"tpoperators.hpp"
#include "files.hpp"
#include"FTLanczos.hpp"
using namespace Many_Body;
    using Mat= Operators::Mat;
int main(int argc, char *argv[])
{
    size_t M{};
  size_t L{};
  double t0{};
  double omega{};
  double gamma{};
   bool PB{};
    double mean= 0.5*omega*L*M;
  double beta=22.;

  try
  {
    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("L", value(&L)->default_value(4), "L")
      ("M", value(&M)->default_value(2), "M")
      ("t", value(&t0)->default_value(1.), "t0")
      ("gam", value(&gamma)->default_value(1.), "gamma")
      ("omg", value(&omega)->default_value(1.), "omega");
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

     

  
  double T=1./beta;
   using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
        PhononBasis g2{ 2, 1};
  ElectronBasis e( L, 1);

  
  PhononBasis ph(L, M);

  HolsteinBasis TP(e, ph);


        Mat E1=Operators::EKinOperatorL(TP, e, t0, PB);
       Mat Ebdag=Operators::BosonCOperator(TP, ph, gamma, PB);
       std::cout<< "dim "<< TP.dim<< std::endl;
       Mat Eb=Operators::BosonDOperator(TP, ph, gamma, PB);
       Mat Eph=Operators::NumberOperator(TP, ph, omega,  PB);
      
      //Mat E=Operators::NumberOperatore(TP, e, 1, PB);
      //    std::cout<< HH << std::endl;
      Eigen::VectorXd eigenVals(TP.dim);
       	Mat H=E1+Eph +Ebdag + Eb;
	auto O=Operators::NumberOperator(TP, ph, omega,  PB);
	std::vector<Mat> v{H, O};
       //       auto HH=Eigen::MatrixXd(H);

      auto ev=Eigen::VectorXd(H.rows());
      //    diagMat(HH, ev);


	auto o=LTLM(H, v, T, 10);


  
	           std::cout<< "for beta/mean = " << beta << std::endl;
       for(auto& l: o)
   	{std::cout<< l <<std::endl; }
  
  return 0;
}

