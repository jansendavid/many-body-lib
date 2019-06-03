 #ifndef FTLANCZOS_H
#define FTLANCZOS_H
#include"basis.hpp"
#include"lanzcos.hpp"
#include"diag.h"

void FTLM()
{}
template<typename SparseMatrix>
void LTLM(SparseMatrix& hamiltonian, SparseMatrix& observable)
{
size_t lanczosDim{10};
 size_t runs=10;
// make initial vector
 auto iniState=Eigen::VectorXcd::Random(hamiltonian.rows());
    Eigen::MatrixXcd Q(hamiltonian.rows(), lanczosDim);
// determine a better way to deccide on Lanczos dim
Many_Body::TriDiagMat tri1=Many_Body::Lanczos(hamiltonian, iniState, lanczosDim, Q);

}

#endif /* FTLANCZOS_H */
