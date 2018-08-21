#pragma once
namespace Many_Body{
  template<typename T>
  void ToFile(const T& x, const T& y, const std::string s1, const size_t N)
{
  // if(x.n_elem!= y.n_elem)
  //   {
  //     std::cout << "Vector dimensions are different " << std::endl;
  //   }

    std::fstream F;
    F.open(s1, std::ios::out| std::ios::app);
    for(size_t i=0; i<N; i++ )
    {
      F<<std::fixed <<  std::setprecision(15)<< x[i] << '\t' << std::setprecision(15) << y[i] << '\n';
    }
  F.close();
  std::cout << "yout data should be stored in the dirctory as " << s1  << std::endl;
return;
}
    template<typename T>
    void ToFile(const T& x, const std::string s1, const size_t N)
{
  // if(x.n_elem!= y.n_elem)
  //   {
  //     std::cout << "Vector dimensions are different " << std::endl;
  //   }

    std::fstream F;
    F.open(s1, std::ios::out| std::ios::app);
    for(size_t i=0; i<N; i++ )
    {
      F<<std::fixed <<  std::setprecision(15)<< x[i] << '\n';
    }
  F.close();
  std::cout << "yout data should be stored in the dirctory as " << s1  << std::endl;
return;
}
}
