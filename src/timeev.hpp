#pragma once
#include "timeev.hpp"
namespace TimeEv{
    using Many_Body::im;
  template <typename State, typename Matrix >
  State timeev_exact(State& initialState, const Matrix& hamiltonian, const Matrix&  observable, State& eigenVals, Eigen::VectorXd output, double dt, size_t numberOfSteps)
  {
    Matrix hamiltonianAdj=hamiltonian.adjoint();
    Eigen::VectorXcd expEigenVal=exp(-eigenVals*im*dt);
    for (size_t i = 0; i < N; ++i)
      {
	initialState= hamiltonian*expEigenValexps*hamAdj*initialState;
	output[i]=initialState.adjoint*observable*initialState;
      }
    return output; 
  }
}
