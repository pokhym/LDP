#include <iostream>
#include <iomanip>
#include <vector>
#include "data.h"
#include "util.h"

using namespace std;

void printall(dataSet &data){
    int m=data.get_m();
    int n=data.get_n();

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            cout<<setw(8)<<data.get_dataMtx(i,j);
        }
        cout<<endl;
    }
}


int main(int argc, const char* argv[]){
    if(argc!=2){
        cout<<"usage: ./test path_to_data"<<endl;
        return 0;
    }
    dataSet *test_data=new(dataSet);
    
    cout<<"begin parsing data"<<endl;
    
    // parse data
    int lol;
    lol=parseData(argv[1], *test_data);
    if(lol==-1){
        cout<<"failed to parse data you might want to check the path"<<endl;
        return -1;
    }

    printall(*test_data);

    cout<<"parse end"<<endl;

    int m=test_data->get_m();
    int n=test_data->get_n();
    cout<<"m (users): "<<m<<" n (attributes): "<<n<<endl;

    // calculate real averages
    vector<double> avg_real(test_data->get_n(), 0.0);
    vector<int> asdf(test_data->get_n(), 0);

    for(int i=0; i<test_data->get_m(); i++){
        for(int j=0; j<test_data->get_n(); j++){
            if(test_data->get_dataMtx(i,j)!=0){
                avg_real[j]=avg_real[j]+test_data->get_dataMtx(i, j);
                asdf[j]++;
            }
        }
    }

    for(int i=0; i<test_data->get_n(); i++){
        avg_real[i]=avg_real[i]/asdf[i];
    }
 
    // data normalization
    cout<<"normalizing columns"<<endl;
    pair<vector<double>, vector<double>> scales=normalizeNeg1toPos1(*test_data);

    printall(*test_data);

    // perturb data
    vector<double> perturb_scales=tuplePerturbation(*test_data,0.5);
    printall(*test_data);
    
    /*
    // calculate estimated average
    vector<double> avg(test_data->get_n(), 0.0);
    vector<int> counters(test_data->get_n(), 0);

    for(int j=0; j<test_data->get_n(); j++){
        for(int i=0; i<test_data->get_m(); i++){
            if(test_data->get_dataMtx(i,j)!=0){
                avg[j]=avg[j]+test_data->get_dataMtx(i,j);
                counters[j]++;
            }
        }
    }

    for(int i=0; i<test_data->get_n(); i++){
        avg[i]=avg[i]/counters[i];
    }   
    
    // undo normalizaion
    undoNorm(scales, avg);
    
    // print real averages and estimated averages
    for(int i=0; i<test_data->get_n(); i++){
         cout<<"avg_real: "<<avg_real[i]<<setw(10)<<" avg_est: "<<avg[i]<<endl;   
    }
    */

    delete test_data;
    return 0;
}
