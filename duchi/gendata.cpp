#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <vector>

using namespace std;

int main()
{
  ofstream myfile;
  myfile.open("data.txt");

  const int rows=10000;
  const int cols=100;

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(58.0,5.0);

  for (int i=0; i<rows; ++i) {
      for(int j=0; j<cols; ++j){
          if(j!=(cols-1)){
              double number = distribution(generator);
              myfile<<number<<",";
          }
          else{
              double number = distribution(generator);
              myfile<<number;
          }
      }
      myfile<<std::endl;
  }

  myfile.close();

  return 0;
}
