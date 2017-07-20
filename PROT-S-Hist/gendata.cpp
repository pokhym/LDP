#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <vector>

#define SCALAR 111

using namespace std;

int main()
{
  ofstream myfile;
  myfile.open("data.txt");

  const int rows=100000;
  const int cols=1;

  std::default_random_engine generator;
  // std::normal_distribution<double> distribution(58.0,5.0);
  std::bernoulli_distribution distribution(0.5);

  for (int i=0; i<rows; ++i) {
      for(int j=0; j<cols; ++j){
          if(j!=(cols-1)){
              double number = distribution(generator);
              myfile<<number*SCALAR<<",";
          }
          else{
              double number = distribution(generator);
              myfile<<number*SCALAR;
          }
      }
      myfile<<std::endl;
  }

  myfile.close();

  return 0;
}
