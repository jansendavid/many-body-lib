#include<iostream>
//#define EIGEN_USE_MKL_ALL
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
  using Fermi_HubbardBasis= TensorProduct<ElectronBasis, ElectronBasis>;
   using HolsteinBasis= TensorProduct<ElectronBasis, PhononBasis>;
  using Mat= Operators::Mat;
   using boost::program_options::value;
   int M{};
  int L{};
    double T{};
  double t0{};
  double omega{};
  double gamma{};
  double cutoff{};
  bool PB{};


  std::string sM{};
  std::string sL{};
  std::string star{};
  std::string sT{};
  std::string st0{};
  std::string somega{};
  std::string sgamma{};
  std::string starget{};
  std::string scutoff{};
  std::string sPB{};
  std::string filename="exT";
  try
  {
    boost::program_options::options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("L", value(&L)->default_value(4), "L")
      ("T", value(&T)->default_value(0.1), "T")
      ("M,m", value(&M)->default_value(2), "M")
      ("t0", value(&t0)->default_value(1.), "t0")
      ("gam", value(&gamma)->default_value(1.), "gamma")
      ("omg", value(&omega)->default_value(1.), "omega")
      ("pb", value(&PB)->default_value(0), "pb");

      

  


    boost::program_options::variables_map vm;
    boost::program_options::store(parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

  if (vm.count("help"))
      {std::cout << desc << '\n'; return 0;}
    else{
      if (vm.count("L"))
      {      std::cout << "L: " << vm["L"].as<int>() << '\n';
      	sL="L"+std::to_string(vm["L"].as<int>());
	filename+=sL;
      }
         if (vm.count("M"))
      {      std::cout << "M: " << vm["M"].as<int>() << '\n';
      	sM="M"+std::to_string(vm["M"].as<int>());
	filename+=sM;
      }
      	 if (vm.count("t0"))
      {      std::cout << "t0: " << vm["t0"].as<double>() << '\n';
      	st0="t0"+std::to_string(vm["t0"].as<double>()).substr(0, 3);
      	filename+=st0;
      }
      	 	 if (vm.count("omg"))
      {      std::cout << "omega: " << vm["omg"].as<double>() << '\n';
      	somega="omg"+std::to_string(vm["omg"].as<double>()).substr(0, 3);
      		filename+=somega;
      }
		 
      		 if (vm.count("gam"))
      {      std::cout << "gamma: " << vm["gam"].as<double>() << '\n';
      	sgamma="gam"+std::to_string(vm["gam"].as<double>()).substr(0, 3);
      		filename+=sgamma;
      }
      // 		 if (vm.count("T"))
      // {      std::cout << "T: " << vm["T"].as<double>() << '\n';
      // 	sT="T"+std::to_string(vm["T"].as<double>()).substr(0, 3);
      // 	filename+=sT;
      // }
		 if (vm.count("pb"))
      {      std::cout << "PB: " << vm["pb"].as<bool>() << '\n';
      	sPB="PB"+std::to_string(vm["pb"].as<bool>());
      		filename+=sPB;
      }
      }
    }
  catch (const boost::program_options::error &ex)
  {
    std::cerr << ex.what() << '\n';
    return 0;
  }
  filename+=".bin";
  double mean=0.5*L*M*omega;
  double dt=0.1;
  ElectronState estate1(L, 1);
  ElectronBasis e( L, 1);
  std::cout<< e<<std::endl;
  
  PhononBasis ph(L, M);
  //                 std::cout<< ph<<std::endl;
      HolsteinBasis TP(e, ph);
      std::cout<< TP.dim << std::endl;             
         Eigen::VectorXcd inistate(TP.dim);

   	 Mat O=Operators::NumberOperator(TP, ph, 1,  PB);
   inistate.setZero();

   inistate[0]=1;
   	          Eigen::VectorXcd i0=inistate;
      
   		  std::cout<< TP.dim << std::endl;     
      Mat E1=Operators::EKinOperatorL(TP, e, t0,PB);
       Mat Ebdag=Operators::NBosonCOperator(TP, ph, gamma, PB);
       Mat Eb=Operators::NBosonDOperator(TP, ph, gamma, PB);
      Mat Eph=Operators::NumberOperator(TP, ph, omega,  PB);

      Eigen::VectorXd eigenVals(TP.dim);
      Mat H=E1+Eb+Eph+Ebdag;
      

    Eigen::MatrixXd HH=Eigen::MatrixXd(H);
    //         std::cout<< HH<<std::endl;
        Eigen::MatrixXd N=Eigen::MatrixXd(Eph);
	  Many_Body::diagMat(HH, eigenVals);
     std::cout<<"MIN E "<<std::setprecision(15)<< eigenVals(0)<< std::endl;
     std::cout<< "MAX "<<eigenVals(eigenVals.size()-1)<< std::endl;
    		      std::cout<<endl<<eigenVals.mean()-mean<<std::endl;
         Eigen::VectorXd energy(TP.dim);
    	   std::vector<double> obs(TP.dim);
    	   std::vector<double> obs2(TP.dim);
    	   std::vector<double> ensvec(TP.dim);
	     
	   Eigen::MatrixXd PHD2=HH.transpose()*N*HH;
	    
	    Eigen::VectorXd v=PHD2.diagonal();
	    std::cout<< "GS "<<eigenVals(0)<<std::endl;
    	   for(int i=0; i<energy.size(); i++)
    {
      //obs[i]=real(obstot[i]);
      
      ensvec[i]=eigenVals(i);
      obs[i]=v(i);

      //std::cout<< ensvec[i]<<'\n';
    }// 	   	    std::cout<< "gS "<<   eigenVals[0] <<std::endl;
    		    std::cout<< "mean "<<   eigenVals.mean() <<std::endl;
    	   std::vector<double> Ovec;

   	    	   std::vector<double> Evec;
   	    	   //		   std::vector<double> Tr={0.01, 0.05, 0.1, 0.15, 0.2, 0.25};
   	    	   		   std::vector<double> Tr;
   	    	   for(int i=1; i<11; i++)
   	    	     {
   	    	       Tr.push_back(i*0.1);
		       
   	    	     }

		   
   	    	  	  for(auto t: Tr){
	     
   	    	 // //	      std::cout<< l<<std::endl;
   	         // bin_write(filename, l);
   	    	 // n++;
   	    		    std::cout<< " at T "<< t<<std::endl;
   	    		    std::cout<<std::setprecision(15)<<expvalCan(ensvec, ensvec,  t)<<std::endl;
   	    	      Ovec.push_back(expvalCan(ensvec, obs,  t));
   	    	      Evec.push_back(expvalCan(ensvec, ensvec,  t)); 

   	    }
 
   			  int pb=PB;
   		  	  std::string Hs="xH";
   		  	  	  std::string Ts="xT";
   		  		  std::string phds="PHD";
	  
		      bin_write("H"+filename, Evec);
   		   	         bin_write("Nph"+filename, Ovec);
   		   		 	         bin_write("temp"+filename, Tr);
  return 0;
}
 
