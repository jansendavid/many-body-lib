# many-body-lib
This is a libray written to do exact diagoanlization of 1d qunatum lattice Hamiltonians which is often used in ED and MBL studies. Typical Hamiltonians which can be studied are the Heisenberg chain, the Holstein model and Bose and Fermi Hubbard model as well as many others which might will only require minor implementations for the user. The code currently do not incoorporate symmetries, but both U1 and translation symmtry can be easily implemented. The code also includes time evolution with the Lanczos algorithms and offersinterface for different eigensolvers depeneding on the needs of the user. 
This code relies heavily on Eigen and is most efficient when combined with MKL. 

 
