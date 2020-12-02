#pragma once
#include<vector>
#include<complex>
#include<iostream>
#include<algorithm>
#include"mkl_lapacke.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace Many_Body{
  

  void diag(std::complex<double>* aH, double* ev, size_t M);
   void diagzheev(std::complex<double>* aH, double* ev, size_t M);
  void diag(double* aH, double* ev, size_t M);
void diagzheevr(std::complex<double>* aH, std::complex<double>* z, double* ev, size_t M);
  void diagOnlyEv(std::complex<double>* aH, double* ev, size_t M);
  void diagOnlyEv(double* aH, double* ev, size_t M);

  template<typename T>
    void diagMat(T& aH, Eigen::VectorXd& ev)
  {
    diag(aH.data(), ev.data(), aH.cols());
  }

    template<typename T>
    void diagMatzheev(T& aH, Eigen::VectorXd& ev)
  {
    diagzheev(aH.data(), ev.data(), aH.cols());
  }
    template<typename T>
    void diagMatzheevr(T& aH, Eigen::VectorXd& ev)
  {
    T EVecs=aH;
    diagzheevr(aH.data(), EVecs.data(), ev.data(), aH.cols());
    aH=EVecs;
  }
   template<typename T>
    void diagMatOnlyEv(T& aH, Eigen::VectorXd& ev)
  {
    diagOnlyEv(aH.data(), ev.data(), aH.cols());
  }

   
       


};
