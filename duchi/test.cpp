#include <iostream>
#include <cmath>
#include <vector>

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

    for(int i=0; i<10000; i++){
        int n=rand()%(20/2)+20/2+1;
        if(n<11){
            cout<<n<<endl;
        }
    }

    vector<double> asdf(20, 0.0);
    cout<<asdf.size()<<endl;
    for(int i=0; i<asdf.size(); i++){
        cout<<asdf[i]<<" ";
    }
    cout<<endl;
    return 0;
}
