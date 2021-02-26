#pragma once
#include<cmath>
#include<complex>
namespace Many_Body{
const double pi = M_PI;
  static std::complex<double> im{0,1};
    static constexpr double err= 0.000000001;
size_t Factorial(const size_t n)
{
  if(n>20)
    {
      std::cout << "to large number " << std::endl;
}
   size_t i, x = 1;
   for (i = 1; i <= n ; i++)
   {
      x *= i;
   }
  return x;
  }
double PrimeNumber(const size_t n)
{
      size_t check,c=0;
    for(size_t i=2;i<=1000;i++)
      {
        check=0;

        for(size_t j=2;j<=i/2;j++)
        {
            if(i%j==0)
            {
              check=1;
               break;
            }
        }

      if(check==0)
        c++;

          if(c==n)
         {
          return i;
          break;
         }
     }
    return 0;
}
}

