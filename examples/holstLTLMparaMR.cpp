#include<iostream>
#include <Eigen/Eigenvalues> 
#include"numerics.hpp"
#include"reddm.hpp"
#include"tpoperators.hpp"
#include "files.hpp"
#include"FTLanczos.hpp"
#include<boost/mpi.hpp>
#include<boost/mpi/communicator.hpp>
#include<boost/mpi/environment.hpp>
#include <boost/program_options.hpp>
using namespace Many_Body;
using Mat= Operators::Mat;
namespace mpi = boost::mpi;
using namespace boost::program_options;
int main(int argc, char *argv[])
{


  int M{};
  int L{};
  double t0{};
  double omega{};
  double gamma{};
  bool PB{0};
   int runs={};
      int rep={};
     int start={};
   int Ldim={};
    double T{};
    double err{};
  std::string sLdim={};
  std::string sruns={};
    std::string srep={};
std::string sT{};
std::string sM{};
std::string sL{};
std::string st0{};
std::string somega{};
std::string sgamma{};
std::string sPB{};
 std::string serr{};
 std::string filename="LTLM";
   

  try
  {
    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("L", value(&L)->default_value(4), "L")
      ("M", value(&M)->default_value(2), "M")
      ("r", value(&runs)->default_value(20), "r")
      ("rep", value(&rep)->default_value(1), "rep")
      ("start", value(&start)->default_value(0), "start")
      ("Ld", value(&Ldim)->default_value(20), "Ld")
      ("t0", value(&t0)->default_value(1.), "t0")
      ("gam", value(&gamma)->default_value(1.), "gamma")
      ("omg", value(&omega)->default_value(1.), "omega")
      ("T", value(&T)->default_value(1.), "T")
      ("pb", value(&PB)->default_value(false), "PB")
      ("err", boost::program_options::value(&err)->default_value(1E-9), "err");
    boost::program_options::variables_map vm;
    boost::program_options::store(parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
 if (vm.count("help"))
      {std::cout << desc << '\n'; return 0;}
    else{
      if (vm.count("L"))
      {
  	std::cout << "L: " << vm["L"].as<int>() << '\n';
	sL="L"+std::to_string(vm["L"].as<int>());
	filename+=sL;
	
      }
     if (vm.count("M"))
      {
  	std::cout << "M: " << vm["M"].as<int>() << '\n';
	sM="M"+std::to_string(vm["M"].as<int>());
	filename+=sM;
	
      }
     if (vm.count("r"))
      {
  	std::cout << "runs: " << vm["r"].as<int>() << '\n';
	sruns="r"+std::to_string(vm["r"].as<int>());
	filename+=sruns;
	
      }
     if (vm.count("rep"))
      {
  	std::cout << "rep: " << vm["rep"].as<int>() << '\n';
	srep="rep"+std::to_string(vm["rep"].as<int>());
	filename+=srep;
	
      }
if (vm.count("Ld"))
      {
  	std::cout << "lanczos dim: " << vm["Ld"].as<int>() << '\n';

	
      }
      if (vm.count("t"))
      {
      	std::cout << "t0: " << vm["t"].as<double>() << '\n';
      	      	st0="t0"+std::to_string(vm["t0"].as<double>()).substr(0, 3);
      	filename+=st0;
      }
       if (vm.count("omg"))
      {
      	std::cout << "omega: " << vm["omg"].as<double>() << '\n';
      	      	somega="omg"+std::to_string(vm["omg"].as<double>()).substr(0, 3);
      	filename+=somega;
      }
       if (vm.count("gam"))
      {
      	std::cout << "gamma: " << vm["gam"].as<double>() << '\n';
      	sgamma="gam"+std::to_string(vm["gam"].as<double>()).substr(0, 3);
      	filename+=sgamma;
      }
       if (vm.count("T"))
      {
      	std::cout << "T: " << vm["T"].as<double>() << '\n';
      }if (vm.count("pb"))
      {
      	std::cout << "PB: " << vm["pb"].as<bool>() << '\n';
      	sPB="PB"+std::to_string(vm["pb"].as<bool>());
      	filename+=sPB;
      }
      		     if (vm.count("err"))
      {      std::cout << "error: " << vm["err"].as<double>() << '\n';
      	 std::stringstream ss;
      	 ss<<vm["err"].as<double>();
      	 serr="err"+ss.str();
      	filename+=serr;
      }
    }
  }
  catch (const error &ex)
  {
    std::cerr << ex.what() << '\n';
    return 0;
  }
  filename+=".bin";
     using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
 


     std::vector<double> beta;
     std::vector<double> Tem;
     for(int i=1; i<11; i++)
       {
	  beta.push_back(1./(0.1*i));
	  Tem.push_back((0.1*i));
       }
     
     
    Mat H;
    Mat N;
    Mat EK;
    Mat X;
  std::vector<Mat> obs;
  {
    
    ElectronBasis e( L, 1);
  PhononBasis ph(L, M);
  HolsteinBasis TP(e, ph);
  std::cout<< TP.dim<<std::endl;
        Mat E1=Operators::EKinOperatorL(TP, e, t0,PB);
       Mat Ebdag=Operators::NBosonCOperator(TP, ph, gamma, PB);
       Mat Eb=Operators::NBosonDOperator(TP, ph, gamma, PB);
      Mat Eph=Operators::NumberOperator(TP, ph, omega,  PB);

      Eigen::VectorXd eigenVals(TP.dim);
       H=E1+Eb+Eph+Ebdag;
        N=Eph/omega;
	EK=E1;
	X=(Ebdag+Eb)/gamma;
       obs.push_back(H);
       obs.push_back(N);
       obs.push_back(E1);
       obs.push_back(X);
  }

  mpi::environment env;
  mpi::communicator world;

    for(int l=start; l<rep; l++)
      {
	if(l%world.size()!=world.rank())
	  {
	    continue;
	  }
	else{
	  std::string sl=std::to_string(l);
     Eigen::MatrixXd Astot=Eigen::MatrixXd::Zero(beta.size(), obs.size());
Eigen::VectorXd Zstot=Eigen::VectorXd::Zero(beta.size());
#pragma omp parallel
  {

#pragma omp for 
 for(int i=0; i<runs; i++)
   {
     std::cout<< " om num thread " <<omp_get_thread_num()<<" from process " << world.rank() <<std::endl;
     auto [Observables, SUMs]=calculate_lanczLT_fast(obs[0], obs, beta, Ldim, err);

	#pragma omp critical
	{
	  Astot+=Observables;
	  Zstot+=SUMs;
	}
	//	std::cout << " got meanZ "<< SUMs.mean() << std::endl;
    
   }
  }

  std::cout<< "Astot  "<<std::endl<< Astot<< std::endl;
  for(size_t i=0; i<beta.size(); i++)
    {
      Astot.row(i)/=Zstot(i);
  	   std::cout<<" T "<< 1./beta[i] << "  "<<Astot(i, 0)<<" SUM "<< Astot(i, 1)+Astot(i, 2)+Astot(i, 3)*gamma<<std::endl;
    }
  bin_write("ER"+sl+filename, Eigen::VectorXd(Astot.col(0)));
  bin_write("NphR"+sl+filename,  Eigen::VectorXd(Astot.col(1)));
  bin_write("EKR"+sl+filename, Eigen::VectorXd(Astot.col(2)));
  bin_write("nXR"+sl+filename, Eigen::VectorXd(Astot.col(3)));
  bin_write("tempR"+sl+filename, Tem);

      }
      }
  return 0;
}
