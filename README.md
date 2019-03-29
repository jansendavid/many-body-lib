# many-body-lib
This is a libray written to do exact diagoanlization of 1d qunatum lattice Hamiltonians which is often used in ED and MBL studies. Typical Hamiltonians which can be studied are the Heisenberg chain, the Holstein model and Bose and Fermi Hubbard model as well as many others which might will only require minor implementations for the user. The code currently do not incoorporate symmetries, but both U1 and translation symmtry can be easily implemented. The code also includes time evolution with the Lanczos algorithms and offersinterface for different eigensolvers depeneding on the needs of the user. 
This code relies heavily on Eigen and is most efficient when combined with MKL. 
For example for doing calculations on the Fermi-Hubbard mode in the position basis we need the tensorproduct basis of two the spin up and spin down system
To do this, 
 ```C++
 #include"Many_Body.hpp"
using namespace Many_Body;
using namespace Eigen;
int main()
{
// declaring the fermi-hubbard basis as the tensor product of two electron basises 
using Fermi_HubbardBasis= TensorProduct<ElectronBasis, ElectronBasis>;
using Mat= Operators::Mat; // mat is a real eigen sparse matrix
size_t L=10;
size_t electrons=int(L/2); // define the number of electrons we want
 ElectronBasis e1(electrons, L); // defining two electronbasisies, each with L/2 electrons
 ElectronBasis e2(electrons, L);
 Fermi_HubbardBasis TP(e1, e2); // this is now the tensor product basis
 // If we now want to make a Hamiltonian with periodic boundary condidtion
    Mat E1=Operators::EKinOperatorL(TP, e, t1);
    Mat E2=Operators::EKinOperatorR(TP, e2, t2);
    Mat C=Operators::CalculateCouplungOperator(TP, e2, u); // this is u*sum_i n_{i, up}n_{i, down}
    Mat H=E1+E2+C;
    // and now to diagonalize it
    Eigen::MatrixXd HH=Eigen::MatrixXd(H); // since we want the whole spectrum and our
    //eigensolvers only work for complete matricies
    VectorXd en(TP.dim); // to store eigen vetcors
    diagMat(HH, en); // HH now contains the eigen vectors an en the eigenvalues, call MKL routine
    // different eieensolvers are implemented
}
```
