#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
//#define EIGEN_USE_MKL_ALL 
#include <boost/test/unit_test.hpp>
#include"basis.hpp"
#include"operators.hpp"
#include"tpoperators.hpp"
#include"diag.h"
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;
#include <Eigen/Eigenvalues>
using namespace Many_Body;
BOOST_AUTO_TEST_SUITE(basistesting)

BOOST_AUTO_TEST_SUITE(operatortesting)
BOOST_AUTO_TEST_CASE(operatoroperations)
{
  using namespace Eigen;
    using namespace std;
  MatrixXd A = MatrixXd::Random(1000,1000);
  VectorXd esx(1000);
  SelfAdjointEigenSolver<MatrixXd> es(A);
  //diagMat(A, esx );
   cout << "The eigenvalues of A are:" << endl << es.eigenvalues().rows()<< endl;
  //cout << "The eigenvalues of A are:" << endl << esx.rows()<< endl;



}      
}
// BOOST_AUTO_TEST_CASE(electrondimension)
// {
// }
BOOST_AUTO_TEST_SUITE_END()
// EOF


