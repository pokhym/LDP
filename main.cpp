#include <iostream>
#include <iomanip>
#include <vector>
#include "data.h"
#include "util.h"

using namespace std;

int main(){
    dataSet *test_data=new(dataSet);
    
    cout<<"begin parsing data"<<endl;

    parseData("../data-Gaussian.csv", *test_data);
   // parseData("testdata.txt", *test_data);

    cout<<"parse end"<<endl;

    int m=test_data->get_m();
    int n=test_data->get_n();
    cout<<"m: "<<m<<" n: "<<n<<endl;
    /*
     *for(int i=0; i<test_data->get_m(); i++){
     *    for(int j=0; j<test_data->get_n(); j++){
     *        cout<<setw(5)<<test_data->get_dataMtx(i,j);
     *    }
     *    cout<<endl;
     *}
     *cout<<endl;
     */

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
    cout<<"asdf"<<endl;
    for(int i=0; i<test_data->get_n(); i++){
        cout<<asdf[i]<<" ";
        avg_real[i]=avg_real[i]/asdf[i];
        // cout<<avg_real[i]<<endl;
    }
 

    cout<<"normalizing columns"<<endl;

    vector<double> outlier={0};
    pair<vector<double>, vector<double>> scales=normalizeNeg1toPos1(*test_data, outlier);
    for(int i=0; i<test_data->get_m(); i++){
        for(int j=0; j<test_data->get_n(); j++){
            if(test_data->get_dataMtx(i,j)>1 || test_data->get_dataMtx(i,j)<-1)
                cout<<"error"<<endl;
        }
    }




    vector<double> perturb_scales=tuplePerturbation(*test_data,0.5);
    
    int number_non_zero=0;
    for(int i=0; i<test_data->get_m(); i++){
        for(int j=0; j<test_data->get_n(); j++){
            if(test_data->get_dataMtx(i,j)!=0)
                number_non_zero++;
        }
    }
    cout<<"number of nonzero entries in matrix: "<<number_non_zero<<endl;
    
    //cout<<"rand column selected: "<<column+1<<endl;
    //int count=0;
    //for(int i=0; i<test_data->get_m(); i++){
        //for(int j=0; j<test_data->get_n(); j++){
            //if(test_data->get_dataMtx(i,j)!=0)
                //count++;
        //}
        //cout<<count<<endl;
        //count=0;
    //}
    
    int count=0;
    for(int j=0; j<test_data->get_n(); j++){
        for(int i=0; i<test_data->get_m(); i++){
            if(test_data->get_dataMtx(i,j)!=0)
                count++;
        }
        cout<<count<<endl;
        count=0;
    }

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

    undoNorm(scales, avg);

    for(int i=0; i<test_data->get_n(); i++){
         cout<<"avg_real: "<<avg_real[i]<<" avg_est: "<<avg[i]<<" perturb scales: "<<perturb_scales[i]<<endl;   
    }
    delete test_data;
    return 0;
}
