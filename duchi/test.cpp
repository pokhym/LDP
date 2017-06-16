#include <iostream>
#include <cmath>

using namespace std;

long long factorial(int i){
    return tgamma(i+1);
}


int main(){
    srand(time(0));
    cout<<2%2<<endl;
    cout<<1%2<<endl;
    cout<<0%2<<endl;
    cout<<(double)rand()/RAND_MAX<<endl;
    return 0;
}

