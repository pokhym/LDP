#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "data.h"
#include "util.h"

using namespace std;

void printall(dataSet &data){
    // int m=data.get_m();
    // int n=data.get_n();
    //
    // for(int i=0; i<m; i++){
    //     for(int j=0; j<n; j++){
    //         cout<<setw(12)<<data.get_dataMtx(i,j);
    //     }
    //     cout<<endl;
    //     cout<<"next row"<<endl;
    // }
}


int main(int argc, const char* argv[]){
    if(argc!=4){
        cout<<"usage: ./test path_to_data initial_m initial_n"<<endl;
        return 0;
    }
    dataSet *test_data=new(dataSet);

    cout<<"begin parsing data"<<endl;

    // parse data
    int lol;
    cout<<stoi(argv[2])<<endl;
    cout<<stoi(argv[3])<<endl;
    lol=parseData(argv[1], *test_data, stoi(argv[2]), stoi(argv[3]));
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
    vector<double> outlier={0};
    pair<vector<double>, vector<double>> scales=normalizeNeg1toPos1(*test_data, outlier);

    cout<<"done normalizing"<<endl;

    printall(*test_data);
    cout<<endl;

    // perturb data
    vector<double> perturb_scales=tuplePerturbation(*test_data,0.5);
    printall(*test_data);

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
    double rel_error=0.0;
    for(int i=0; i<test_data->get_n(); i++){
         cout<<"avg_real: "<<avg_real[i]<<setw(10)<<" avg_est: "<<avg[i]<<endl;
         rel_error+=abs((avg[i]-avg_real[i])/avg_real[i]);
    }
    cout<<"percent releative error: "<<rel_error/test_data->get_n()<<endl;

    delete test_data;
    return 0;
}
