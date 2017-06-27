#include <iostream>
#include <cmath>
#include <vector>
#include <random>

using namespace std;

long long factorial(int i){
    return tgamma(i+1);
}


int main(){
    default_random_engine gen;
    uniform_int_distribution<int> asdf(0,1);

    for(int i=0; i<10000; i++){
        cout<<asdf(gen)<<endl;
    }

    return 0;
}
