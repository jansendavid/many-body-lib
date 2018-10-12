#pragma once
#include<vector>
#include<complex>
#include<iostream>
#include"mkl_lapacke.h"
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

namespace Many_Body{
  

  void diag(std::complex<double>* aH, double* ev, size_t M);
   void diagzheev(std::complex<double>* aH, double* ev, size_t M);
  void diag(double* aH, double* ev, size_t M);
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
    void diagMatOnlyEv(T& aH, Eigen::VectorXd& ev)
  {
    diagOnlyEv(aH.data(), ev.data(), aH.cols());
  }

   
       


};
