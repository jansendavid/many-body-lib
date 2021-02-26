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
  std::string dirname="";
  std::string filename={};
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
      ("dir", value(&dirname)->default_value(""), "dirname")
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
	filename+="L"+std::to_string(vm["L"].as<size_t>());
      }
     if (vm.count("M"))
      {
	std::cout << "M: " << vm["M"].as<size_t>() << '\n';
	filename+="M"+std::to_string(vm["M"].as<size_t>());
      }
      if (vm.count("t"))
      {
	std::cout << "t0: " << vm["t"].as<double>() << '\n';
	filename+="t0"+std::to_string(vm["t"].as<double>()).substr(0, 3);
      }
       if (vm.count("omg"))
      {
	std::cout << "omega: " << vm["omg"].as<double>() << '\n';
	filename+="omega"+std::to_string(vm["omg"].as<double>()).substr(0, 3);
      }
       if (vm.count("gam"))
      {
	std::cout << "gamma: " << vm["gam"].as<double>() << '\n';
	filename+="gam"+std::to_string(vm["gam"].as<double>()).substr(0, 4);
      }
       if (vm.count("dt"))
      {
	std::cout << "dt: " << vm["dt"].as<double>() << '\n';
	filename+="dt"+std::to_string(vm["dt"].as<double>()).substr(0, 4);
      }
       if (vm.count("tot"))
      {
	std::cout << "total time: " << vm["tot"].as<double>() << '\n';
	filename+="tot"+std::to_string(vm["tot"].as<double>()).substr(0, 3);
      }
       if (vm.count("pb"))
      {
	std::cout << "PB: " << vm["pb"].as<bool>() << '\n';
	filename+="PB"+std::to_string(vm["pb"].as<bool>());
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
  
  //  std::cout<< e<<std::endl;
std::vector<size_t> es(L, 0);
      es[0]=1;
      ElectronState aa(es);
  PhononBasis ph(L, M);
  //       std::cout<< ph<<std::endl;
      HolsteinBasis TP(e, ph);
      std::vector<size_t> state(L, 0);
      //std::fill(state.begin(), state.end(), 1);
  BosonState b2(state, M);
 
     auto it3=TP.lbasis.find(aa.GetId());
     auto it4=TP.rbasis.find(b2.GetId());
     size_t StateNr= Position(*it4)*TP.lbasis.dim +Position(*it3);
     //       std::cout<< TP << std::endl;
     //    std::cout<< b2.GetId()<<std::endl;
     // 	std::cout<< b2.GetId()<<std::endl;
      std::cout<<"state nr "<< StateNr<<std::endl;                
         Eigen::VectorXcd inistate(TP.dim);
	 double para=1;
 	 Mat O0=Operators::NumberOperatore_1(TP, e, para, 0,  PB);
	 Mat O1=Operators::NumberOperatore_1(TP, e, para, 1,  PB);
	 Mat Nph0=Operators::NumberOperatorph_1(TP, ph, para, 0,  PB);
	 Mat Nph1=Operators::NumberOperatorph_1(TP, ph, para, 1,  PB);
	 Mat X0=std::sqrt(1./2)*(Operators::BosonDOperator_1(TP, ph, para,0,  PB)+Operators::BosonCOperator_1(TP, ph, para,0,  PB));
	 Mat X1=std::sqrt(1./2)*(Operators::BosonDOperator_1(TP, ph, para,1,  PB)+Operators::BosonCOperator_1(TP, ph, para,1,  PB));
	 Mat P0=std::sqrt(1./2)*(Operators::BosonCOperator_1(TP, ph, para,0,  PB)-Operators::BosonDOperator_1(TP, ph, para,0,  PB));
	 Mat P1=std::sqrt(1./2)*(Operators::BosonCOperator_1(TP, ph, para,1,  PB)-Operators::BosonDOperator_1(TP, ph, para,1,  PB));
	 Mat Q1=X0-X1;
	 Mat QP1=P0-P1;
	 Mat Q2=Q1*(X0-X1);
	 Mat QP2=QP1*(P0-P1);
	 Mat Q3=Q2*(X0-X1);
	 Mat QP3=QP2*(P0-P1);
	 Mat Q4=Q3*(X0-X1);
	 Mat QP4=QP3*(P0-P1);
	 //	 
   inistate.setZero();

   inistate[StateNr]=1;
   	          Eigen::VectorXcd i0=inistate;
		  std::cout<<"dim " <<e.dim << std::endl;
		  
           std::cout<< TP.dim << std::endl;     
      Mat EK=Operators::EKinOperatorL(TP, e, t0, PB);
      Mat Ebdag=Operators::NBosonCOperator(TP, ph, gamma, PB);
      Mat Eb=Operators::NBosonDOperator(TP, ph, gamma, PB);
      Mat Eph=Operators::NumberOperator(TP, ph, omega,  PB);
      
      //Mat E=Operators::NumberOperatore(TP, e, 1, false);
      //    std::cout<< HH << std::xbendl;
      Eigen::VectorXd eigenVals(TP.dim);
      Mat H=EK+Eph  +Ebdag + Eb;
      // making one Hamiltonian with infinite chem pot on one site
      // double mu=100000;
      // Eigen::VectorXd eigenVals_mit_pot(TP.dim);
      // Mat H_mit_pot=H+Operators::NumberOperatore_1(TP, e, mu, 1,  PB);
      //    
 
    Eigen::MatrixXd HH=Eigen::MatrixXd(H);
    //    Eigen::MatrixXd HH_mit_pot=Eigen::MatrixXd(H_mit_pot);
     Eigen::MatrixXd H_cop=Eigen::MatrixXd(H);
     // Eigen::MatrixXd H_mit_pot_cop=Eigen::MatrixXd(H_mit_pot);
    //        std::cout<< HH<<std::endl;
        Eigen::MatrixXd N=Eigen::MatrixXd(Eph);
	
     Many_Body::diagMat(HH, eigenVals);
     // Many_Body::diagMat(HH_mit_pot, eigenVals_mit_pot);
         Eigen::VectorXd energy(TP.dim);
	 std::vector<double> ek_vec;
    	   std::vector<double> n0_vec;
	   std::vector<double> n1_vec;
	   std::vector<double> nph0_vec;
	   std::vector<double> nph1_vec;
	   std::vector<double> x0_vec;
	   std::vector<double> x1_vec;
	   std::vector<double> p0_vec;
	   std::vector<double> p1_vec;
	   std::vector<double> q1_vec;
	   std::vector<double> qp1_vec;
	   std::vector<double> q2_vec;
	   std::vector<double> qp2_vec;
	   std::vector<double> q3_vec;
	   std::vector<double> qp3_vec;
	   std::vector<double> q4_vec;
	   std::vector<double> qp4_vec;
	   std::vector<double> ensvec_vec;
	   //	   Eigen::MatrixXd PHD2=HH.adjoint()*N.selfadjointView<Lower>()*HH;
	   std::cout<< "eigenVals(0) "<<eigenVals(0)<<std::endl;
	   //	    std::cout<< "eigenValseigenVals_mit_pot(0) "<<eigenVals_mit_pot(0)<<std::endl;
     Eigen::MatrixXcd evExp=TimeEv::EigenvalExponent(eigenVals, dt);
   Eigen::MatrixXcd cEVec=HH.cast<std::complex<double>>();
   // Eigen::MatrixXcd cEVec_mit_pot=HH_mit_pot.cast<std::complex<double>>();
   Eigen::VectorXcd newIn=inistate;
     //cEVec_mit_pot.col(0);
   //
     
   
      //
   int i=0;

   // O=H;
   std::complex<double> E_in=(newIn.adjoint()*(H_cop*newIn))(0);
      std::complex<double> n_in=(newIn.adjoint()*(O0*newIn))(0);
           std::complex<double> N_in=(newIn.adjoint()*(Nph0*newIn))(0);
   std::cout<< " init energy "<<E_in<<std::endl;
   std::cout<< " init n site 0 "<<n_in<<std::endl;
   std::cout<< " init N site 0 "<<N_in<<std::endl;
        while(i*dt<tot)
       {

       	 
  	   	 std::complex<double> ek=(newIn.adjoint()*(EK*newIn))(0);
   	  	 std::complex<double> n0=(newIn.adjoint()*(O0*newIn))(0);
 		 std::complex<double> nph0=(newIn.adjoint()*(Nph0*newIn))(0);
 		 std::complex<double> nph1=(newIn.adjoint()*(Nph1*newIn))(0);
 		 std::complex<double> e0=(newIn.adjoint()*(H_cop*newIn))(0);
 		 std::complex<double> n1=(newIn.adjoint()*(O1*newIn))(0);
 		 std::complex<double> x0=(newIn.adjoint()*(X0*newIn))(0);
 		 std::complex<double> x1=(newIn.adjoint()*(X1*newIn))(0);
 		 std::complex<double> p0=std::complex<double>(0,1)*(newIn.adjoint()*(P0*newIn))(0);
 		 std::complex<double> p1=std::complex<double>(0,1)*(newIn.adjoint()*(P1*newIn))(0);
		 std::complex<double> q1=(newIn.adjoint()*(Q1*newIn))(0);
 		 std::complex<double> qp1=std::complex<double>(0,1)*(newIn.adjoint()*(QP1*newIn))(0);
		 std::complex<double> q2=(newIn.adjoint()*(Q2*newIn))(0);
		 std::complex<double> qp2=std::complex<double>(0,1)*std::complex<double>(0,1)*(newIn.adjoint()*(QP2*newIn))(0);
		 std::complex<double> q3=(newIn.adjoint()*(Q3*newIn))(0);
		 std::complex<double> qp3=std::complex<double>(0,1)*std::complex<double>(0,1)*std::complex<double>(0,1)*(newIn.adjoint()*(QP3*newIn))(0);
		 std::complex<double> q4=(newIn.adjoint()*(Q4*newIn))(0);
		 std::complex<double> qp4=std::complex<double>(0,1)*std::complex<double>(0,1)*std::complex<double>(0,1)*std::complex<double>(0,1)*(newIn.adjoint()*(QP4*newIn))(0);
		 ek_vec.push_back(real(ek));	 
 		 n0_vec.push_back(real(n0));
 		 n1_vec.push_back(real(n1));
 		 nph0_vec.push_back(real(nph0));
 		 nph1_vec.push_back(real(nph1));
 		 x0_vec.push_back(real(x0));
 		 x1_vec.push_back(real(x1));
 		 p0_vec.push_back(real(p0));
 		 p1_vec.push_back(real(p1));
		 q1_vec.push_back(real(q1));
		 qp1_vec.push_back(real(qp1));
		 q2_vec.push_back(real(q2));
		 qp2_vec.push_back(real(qp2));
		 q3_vec.push_back(real(q3));
		 qp3_vec.push_back(real(qp3));
		 q4_vec.push_back(real(q4));
		 qp4_vec.push_back(real(qp4));
 		 //	 std::complex<double> c=(inistate.adjoint()*(newIn))(0);
  

   			// 	outputVals(i)=real(c);
   			// 	outputVals2(i)=real(c2);
        		// outputTime(i)=i*dt;
 		 std::cout<< std::setprecision(8)<< " ene "<< E_in-e0<< " partcle "<<n0 +n1<<"  " << x0 <<x1 << "  " << p0 << "  "  << p1 <<"  "<< " dt "<< i*dt<< std::endl;
 	 //	  	 BOOST_CHECK(std::abs(real(c2)-real(c))<Many_Body::err);
 	 TimeEv::timeev_exact(newIn, cEVec, evExp);
   	 i++;
 	     	             }

 //std::cout<< MatrixXd(OBS2) << std::endl;
 // Eigen::MatrixXd OBS22=HH.adjoint()*(OBS2.selfadjointView<Lower>())*HH;
 // Eigen::VectorXd diagobs2=OBS22.diagonal();
 	std::cout<< n0_vec.size()<< "  i "<< i << std::endl;
 // std::cout<<diagobs2(0)<<std::endl; 
   filename+=".bin";
     bin_write(dirname+"/ek"+filename, ek_vec);
         bin_write(dirname+"/n0"+filename, n0_vec);
 	   bin_write(dirname+"/n1"+filename, n1_vec);
 	   bin_write(dirname+"/nph0"+filename, nph0_vec);
 	   bin_write(dirname+"/nph1"+filename, nph1_vec);
 	     bin_write(dirname+"/x0"+filename, x0_vec);
 	     bin_write(dirname+"/x1"+filename, x1_vec);
 	       bin_write(dirname+"/p0"+filename, p0_vec);
 	         bin_write(dirname+"/p1"+filename, p1_vec);
		 bin_write(dirname+"/q1"+filename, q1_vec);
		 bin_write(dirname+"/qp1"+filename, qp1_vec);
		 bin_write(dirname+"/q2"+filename, q2_vec);
		 bin_write(dirname+"/qp2"+filename, qp2_vec);
		 bin_write(dirname+"/q3"+filename, q3_vec);
		 bin_write(dirname+"/qp3"+filename, qp3_vec);
		 bin_write(dirname+"/q4"+filename, q4_vec);
		 bin_write(dirname+"/qp4"+filename, qp4_vec);
  return 0;
}
 
