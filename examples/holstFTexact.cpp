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
      ("gam", value(&gamma)->default_value(std::sqrt(2)), "gamma")
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
   auto J=Operators::CurOperatorL(TP, e, t0, PB);    
Eigen::MatrixXcd Jcp=J*std::complex<double>(0,1);      

    Eigen::MatrixXd HH=Eigen::MatrixXd(H);
   

    //Eigen::MatrixXcd HHcp=HH;
    Eigen::MatrixXcd HHcp=Eigen::MatrixXd(HH);
             Eigen::MatrixXcd M2=HHcp*Jcp;
             Eigen::MatrixXcd M1=Jcp*HHcp;
	     Eigen::MatrixXcd COMM=M1-M2;
	     //                  std::cout<< M1<<std::endl;
	     //           std::cout<< COMM<<std::endl;
	     //   std::cout<< M2<<std::endl;
	       COMM*=Jcp;
		  
		  //           std::cout<< HHcp<<std::endl;
		  //           std::cout<< Jcp<<std::endl;
        Eigen::MatrixXd N=Eigen::MatrixXd(Eph);
	  Many_Body::diagMat(HH, eigenVals);
     std::cout<<"MIN E "<<std::setprecision(15)<< eigenVals(0)<< std::endl;
     std::cout<< "MAX "<<eigenVals(eigenVals.size()-1)<< std::endl;
    		      std::cout<<endl<<eigenVals.mean()-mean<<std::endl;
         Eigen::VectorXd energy(TP.dim);
    	   std::vector<double> obs(TP.dim);
    	   std::vector<double> obs2(TP.dim);
	   std::vector<double> obs3(TP.dim);
	   std::vector<double> obs4(TP.dim);
	   std::vector<double> obs5(TP.dim);
	   std::vector<double> obs6(TP.dim);
	   //std::vector<double> obs(TP.dim);
    	   std::vector<double> ensvec(TP.dim);
	     
	   Eigen::MatrixXd PHD2=HH.adjoint()*N*HH;
	   Eigen::MatrixXd EKIN=HH.adjoint()*E1*HH;
	   Eigen::MatrixXd EC=HH.adjoint()*(Ebdag+Eb)*HH;
	   Eigen::MatrixXcd Jmat=HH.adjoint()*(Jcp*HH);
	   Eigen::MatrixXcd JJmat=Jmat*Jmat;
	   Eigen::MatrixXcd COMMmat=HH.adjoint()*(COMM*HH);
	   //	   std::complex<double> J_in=(newIn.adjoint()*(Jcp*newIn))(0);
	   // std::complex<double> JJsq_in=(newIn.adjoint()*(Jcp*Jcp*newIn))(0);
	    Eigen::VectorXd v=PHD2.diagonal();
	    Eigen::VectorXd v2=EKIN.diagonal();
	    Eigen::VectorXd v3=EC.diagonal();
	    Eigen::VectorXd v4=(JJmat.diagonal().real());
	     Eigen::VectorXd v5=(Jmat.diagonal().real());
	     Eigen::VectorXd v6=(COMMmat.diagonal().real());
	    std::cout<< "GS "<<eigenVals(0)<<std::endl;
    	   for(int i=0; i<energy.size(); i++)
    {
      //obs[i]=real(obstot[i]);
      
      ensvec[i]=eigenVals(i);
      obs[i]=v(i);
      obs2[i]=v2(i);
      obs3[i]=v3(i);
      obs4[i]=v4(i);
      obs5[i]=v5(i);
      obs6[i]=v6(i);
      std::cout<< COMMmat(i,i)<<std::endl;
      if(gamma!=0)
	{	obs3[i]/=gamma;}

      //std::cout<< ensvec[i]<<'\n';
    }// 	   	    std::cout<< "gS "<<   eigenVals[0] <<std::endl;
    		    std::cout<< "mean "<<   eigenVals.mean() <<std::endl;
    	   std::vector<double> Ovec;
	    std::vector<double> Ovec2;
	     std::vector<double> Ovec3;
	     std::vector<double> Ovec4;
	     	     std::vector<double> Ovec5;

   	    	   std::vector<double> Evec;
   	    	   //		   std::vector<double> Tr={0.01, 0.05, 0.1, 0.15, 0.2, 0.25};
   	    	   		   std::vector<double> Tr;
   	    	   for(int i=1; i<22; i++)
   	    	     {
   	    	       Tr.push_back(i*0.05);
		       
   	    	     }

		   
   	    	  	  for(auto t: Tr){
	     
   	    	 // //	      std::cout<< l<<std::endl;
   	         // bin_write(filename, l);
   	    	 // n++;
			    double e1=expvalCan(ensvec, ensvec,  t);
			    double o=expvalCan(ensvec, obs,  t);
			    double o2=expvalCan(ensvec, obs2,  t);
			    double o3=expvalCan(ensvec, obs3,  t);
			    double o4=expvalCan(ensvec, obs4,  t);
			     double o5=expvalCan(ensvec, obs5,  t);
			     double o6=expvalCan(ensvec, obs6,  t);
   	    		    std::cout<< " at T "<< t<<std::endl;
			    std::cout<<std::setprecision(15)<<e1 <<"  "<< o+o2+o3*gamma <<" Nph "<<o << " JJ " <<o4 << "o 6 2 "<< o6<< std::endl;
			    Ovec.push_back(o);
			    Ovec2.push_back(o2);
			    Ovec3.push_back(o3);
			    Ovec4.push_back(o4);
			    			    Ovec5.push_back(o5);
			    Evec.push_back(e1);

   	    }
 
   			  int pb=PB;
   		  	  std::string Hs="H";
			  std::string Ts="T";
			  std::string phds="PHD";
	  
		      bin_write("E"+filename, Evec);
		      bin_write("Nph"+filename, Ovec);
		      bin_write("EK"+filename, Ovec2);
		      bin_write("nX"+filename, Ovec3);
		      bin_write("JJ"+filename, Ovec4);
		      bin_write("J"+filename, Ovec5);
		      bin_write("temp"+filename, Tr);
  return 0;
}
 
