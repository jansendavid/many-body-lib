
#include<iostream>
#include <Eigen/Eigenvalues> 
#include"numerics.hpp"
#include"reddm.hpp"
#include"tpoperators.hpp"
#include "files.hpp"
#include <boost/program_options.hpp>

using namespace boost::program_options;
using namespace Many_Body;
int main(int argc, char *argv[])
{
  using Mat= Operators::Mat;
  size_t M{};
  size_t L{};
  double t0{};
  double omega{};
  double gamma{};
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

    
using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
   // std::vector<size_t> ee(L, 0);
   //    ee[L-1]=1;
   //    ElectronState aa(ee);
   // ElectronBasis e(aa);
   ElectronBasis e( L, 1);
   //   std::cout<< e<<std::endl;
  
  PhononBasis ph(L, M);
  //  std::cout<< ph<<std::endl;
  HolsteinBasis TP(e, ph);
  //  std::cout<< TP<<std::endl;

  std::cout<<"total dim "<< TP.dim << std::endl;
  std::cout<<std::endl;
        Mat E1=Operators::EKinOperatorL(TP, e, t0, PB);
      Mat Ebdag=Operators::NBosonCOperator(TP, ph, gamma, PB);
      Mat Eb=Operators::NBosonDOperator(TP, ph, gamma, PB);
      Mat Eph=Operators::NumberOperator(TP, ph, omega,  PB);


 //      Eigen::VectorXd eigenVals(TP.dim);
      Mat H=E1  +Ebdag + Eb+ Eph;
      //+phMOM;
      //
      //

	//
      //+
      //+ phKNN phK+;
      
          Eigen::MatrixXcd HH=Eigen::MatrixXcd(H);
	   std::cout<< HH << std::endl;
	  // for(int i=0; i<HH.rows(); i++)
	  //   {
	  //    for(int j=0; j<HH.rows(); j++)
	  //   {
	  //     if(std::abs(HH(i, j))>0.00001){
	  //     	if(std::abs(HH(j, i)-HH(i, j))>0.000001){ std::cout<< "err "<<std::endl;}
	  //     		std::cout<< "went from  "<< i << " to "<< j << "  " <<HH(i, j)<< std::endl;
	  //     }
	  //   }
	  //   }
	  // std::cout<< HH << std::endl;
	     Eigen::VectorXd ev=Eigen::VectorXd(TP.dim);
	     diagMat(HH, ev);
	     std::cout<<"GS "<< std::setprecision(8)<< ev[0]<< std::endl;
	     std::cout<< " nect "<< std::endl;
	     std::cout<< HH.col(0).adjoint()*Eph*HH.col(0);
   // 	     Eigen::MatrixXd MDF2=HH.adjoint()*(O)*HH;
   // Eigen::VectorXd v=MDF2.diagonal();
   // std::cout<< v << std::endl;
   // std::cout <<"mean "<< v.mean()<< std::endl;
 // 	  bool isDiag=false;
 // 	 //  std::vector<double> Tr={0.0100, 0.1000, 1.000, 2.0000, 5.0000, 10.0000};
 // 	 //  for(auto t: Tr){
 // 	 //      std::string sT=std::string(std::to_string(t)).substr(0,6);
 // auto t=0.5;
 // 	       auto optModes=makeThermalRDMTP(HH,ev,  TP, t, isDiag, 0);
 // 	 //  //	  std::cout << "sum of all eigenvalues "<< optModes.sum()<< std::endl;
 // 	 //      isDiag=true;
 // 	  int n=0;
	 
 // 	  for(auto& l : optModes)
 // 	    {
 // 	      	      std::string sn=std::string(std::to_string(n)).substr(0,1);
 // 		      //std::string filename="OML"+std::to_string(L)+"M"+std::to_string(M)+"t0_"+"1.0"+"gam"+sgam+"omg"+ somg+"T"+ sT+ "esec" +sn  + ".bin";
 // 		 std::cout << n << std::endl;
 // 		 std::cout<< "value "<< " fot T = "<< t <<std::endl;
 // 		 	      std::cout<< l<<std::endl;
 // 			      std::cout << " and sum "<< l.sum()<< std::endl;
 // 			      //bin_write(filename, l);
 // 		 n++;
 // 	    }
	  
	  
	  // }
  //   std::cout<< TP<< std::endl;
  // Eigen::VectorXd v1=Eigen::VectorXd::Zero(TP.dim);
  // v1[1]=1/std::sqrt(2);
  // Eigen::VectorXd v2=Eigen::VectorXd::Zero(TP.dim);
  //   v2[5]=1/std::sqrt(2);
  // 	Eigen::VectorXd V1=v1 +v2;
  // 	Eigen::MatrixXd V2=V1.transpose();
  // 	// std::cout<< V1<<std::endl;
  // 	// std::cout<< V2<<std::endl;	
  // 	MatrixXd M=V1*(V2);
  // 	// std::cout<< M<<std::endl;
  // 	// 	std::cout<< g2<<std::endl;
  // 	//	makeRedDM(g2, 0, M);
  // 	makeRedDMTP(TP, ph,  0, M);
  return 0;
}
