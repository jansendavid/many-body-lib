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
  int i{};
  int j{};
    double T{};
      double tot{};
        double dt{};
  double t0{};
  double omega{};
  double gamma{};
  double cutoff{};
  bool PB{};


  std::string sM{};
  std::string sL{};
  std::string si{};
  std::string sj{};
  std::string star{};
  std::string sT{};
    std::string stot{};
      std::string sdt{};
  std::string st0{};
  std::string somega{};
  std::string sgamma{};
  std::string starget{};
  std::string scutoff{};
  std::string sPB{};
  std::string filename="SPexT";
  try
  {
    boost::program_options::options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("L", value(&L)->default_value(4), "L")
      ("i", value(&i)->default_value(0), "i")
      ("j", value(&j)->default_value(1), "j")
      ("T", value(&T)->default_value(0.1), "T")
      ("M,m", value(&M)->default_value(2), "M")
      ("tot", value(&tot)->default_value(1.), "tot")
      ("dt", value(&dt)->default_value(0.1), "dt")
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
	   if (vm.count("i"))
      {      std::cout << "i: " << i << '\n';
      	si="i"+std::to_string(i);
	filename+=si;
      }
	   if (vm.count("j"))
      {      std::cout << "j: " << j << '\n';
      	sj="j"+std::to_string(j);
	filename+=sj;
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
		      		 if (vm.count("tot"))
      {      std::cout << "tot: " << vm["tot"].as<double>() << '\n';
      	stot="tot"+std::to_string(vm["tot"].as<double>()).substr(0, 3);
      		filename+=stot;
      }
				 if (vm.count("dt"))
      {      std::cout << "dt: " << vm["dt"].as<double>() << '\n';
      	sdt="dt"+std::to_string(vm["dt"].as<double>()).substr(0, 3);
      		filename+=sdt;
      }

      		 if (vm.count("T"))
      {      std::cout << "T: " << T << '\n';
      	sT="T"+std::to_string(vm["T"].as<double>()).substr(0, 3);
      	filename+=sT;
      }
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

  // declaring the tensor products with 0 e, 1, e, 2 e
    ElectronBasis e( L, 0);
    e.append(1);
     e.append(2);
     // ElectronBasis e1( L, 1);
    //e0.basis.insert(e1.begin(), e1.end());
      std::cout<< e<<std::endl;
      
      //ElectronBasis e2( L, 2);
  

   PhononBasis ph(L, M);
 HolsteinBasis TP(e, ph);
//       HolsteinBasis TP1(e1, ph);
//         HolsteinBasis TP2(e2, ph);


//std::cout<< TP << std::endl;
 std::cout<< TP.dim << std::endl;
// 		  //std::cout<< ph << std::endl;
// 		  //    std::cout<< e1 << std::endl;
       Eigen::VectorXd eigenVals(TP.dim);
//       Eigen::VectorXd eigenVals2(TP2.dim);
//       Eigen::VectorXd eigenVals0(TP0.dim);
//       Mat H0=Operators::HolsteinHam(TP0, t0, omega, gamma, PB);
       Mat H=Operators::HolsteinHam(TP, t0, omega, gamma, PB);
//       Mat H2=Operators::HolsteinHam(TP2, t0, omega, gamma, PB);
//     Eigen::MatrixXd HH0=Eigen::MatrixXd(H0);
     Eigen::MatrixXd HH=Eigen::MatrixXd(H);
//     Eigen::MatrixXd HH2=Eigen::MatrixXd(H2);
//     //         std::cout<< HH<<std::endl;
//     //  Eigen::MatrixXd N=Eigen::MatrixXd(Eph);
//     Many_Body::diagMat(HH0, eigenVals0);
 	  Many_Body::diagMat(HH, eigenVals);
//        Many_Body::diagMat(HH2, eigenVals2);
//        // 
  auto Cdag=Operators::CdagOperator(TP, i,  PB);
  auto C=Operators::COperator(TP, j,  PB);
  
  std::cout<<"T "<<(Eigen::MatrixXd(Cdag)*Eigen::MatrixXd(C)).trace()<< std::endl;
//        auto C1=Operators::CdagOperatorCN(TP1, TP0,i,  PB);
//        std::cout<<"MIN E0 "<<std::setprecision(15)<< eigenVals0(0)<< std::endl;
      std::cout<<"MIN E1 "<<std::setprecision(15)<< eigenVals(0)<< std::endl;
//      std::cout<<"MIN E2 "<<std::setprecision(15)<< eigenVals2(0)<< std::endl;
//      std::cout<< " dim T0 "<< TP0.dim<<std::endl;
//           std::cout<< " dim T1 "<< TP1.dim<<std::endl;
// 	       std::cout<< " dim T2 "<< TP2.dim<<std::endl;

//      // std::cout<< Eigen::MatrixXd(Cdag) << std::endl;
//      // std::cout<< std::endl;
//      // std::cout<<TP0<<std::endl;
//      //      std::cout<<TP1<<std::endl;
// 	       std::cout<< "mat "<< std::endl;
// 	       std::cout<< Eigen::MatrixXd(C1)*Eigen::MatrixXd(C1).transpose()<< std::endl;
//       std::cout<< "endl "<<std::endl;
//       std::cout<< Eigen::MatrixXd(C2) << std::endl;
//      // std::cout<< "\n";
//      //  std::cout<< Eigen::MatrixXd(Operators::CdagOperator(TP1, TP0,1,  PB)) << std::endl;
      auto N_e=Operators::NumberOperator(TP, e, 1, PB);
        int n=0;
        Eigen::MatrixXcd evExpbeta=Eigen::MatrixXcd::Zero(eigenVals.rows(), eigenVals.rows());
        double Z{0};
	auto it=TP.basis.begin();
       for(int k=0; k<eigenVals.rows(); k++)
	 {
	   auto it2=TP.lbasis.find(LeftId(*it));
	   auto state=GetLattice(*it2);
 	   double ex{0};
	    // if(state.Count()==1)
	    //   {
	   ex=std::exp(-(eigenVals(k))/T);
	   //    }
	   evExpbeta(k, k)=std::complex<double>{ex, 0};
	Z+=ex;
	it++;
	//	sum+=ex*en;
	 }
       
//        //=TimeEv::EigenvalExponent(eigenVals1, -Many_Body::im*1./T);
    
//      std::vector<double> time;
//      std::vector<std::complex<double>> cdagc;
//      std::vector<std::complex<double>> ccdag;
     auto HT=HH.transpose();
     auto CdagTR=HT*Cdag*HH;
     auto CTR=HT*C*HH;
     std::cout<<"T "<<(Eigen::MatrixXd(CdagTR)*Eigen::MatrixXd(CTR)).trace()/Z<< std::endl;
     //   std::cout<<"T2 "<<HT*HH<< std::endl; Eigen::MatrixXd(CTR)*Eigen::MatrixXd(CdagTR)).trace()/7
     auto H2T=HH2.transpose();
       auto H0T=HH0.transpose();
        while(n*dt<tot)
        	{
     	  	  

        Eigen::MatrixXcd evExpm=TimeEv::EigenvalExponent(eigenVals, n*dt);
           Eigen::MatrixXcd evExpp=TimeEv::EigenvalExponent(eigenVals, -n*dt);
	   std::cout<< "Z " << Z << std::endl;

   	 auto CDAGC=(evExpbeta*CdagTR*evExpp*CTR*evExpm).trace()/Z;

	 auto CCDAG= (evExpbeta*evExpp*CTR*evExpm*CdagTR).trace()/Z;
;
   		 std::cout<<n*dt << "\t"<< CCDAG<< " and  "<<CDAGC<< "diff "<< CCDAG+CDAGC<< std::endl;
  		      //	std::cout<<n*dt << "\t"<<  " and  "<<CCDAG<< "  "<< "  sum  " <<Z<<std::endl;
       			    // 	    time.push_back(n*dt);
     			    // 	    ccdag.push_back(CCDAG);
     			    // cdagc.push_back(CDAGC);
  n++;

  	 	}
  // 	std::cout<< ccdag.size() << "  "<< cdagc.size()<< std::endl;
  //     bin_write("time"+filename, time);
  //     bin_write("CCDAG"+filename, ccdag);
  //     bin_write("CDAGC"+filename, cdagc);
  
 
    // 			  int pb=PB;
    // 		  	  std::string Hs="H";
    // 			  std::string Ts="T";
    // 			  std::string phds="PHD";
	  
    // 		      bin_write("E"+filename, Evec);
    // 		      bin_write("Nph"+filename, Ovec);
    // 		      bin_write("EK"+filename, Ovec2);
    // 		      bin_write("nX"+filename, Ovec3);
    // 		      bin_write("temp"+filename, Tr);
  return 0;
}
 
