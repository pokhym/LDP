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

void test_basic_randomizer(){
    cout<<"testing basic randomizer"<<endl;
    // generate random vector
    vector<double> test_vec(144, 0.0);
    srand(0);
    for(int i=0; i<(int)test_vec.size(); i++){
        int idx_sign=rand()%2;
        if(idx_sign==1){
            test_vec[i]=1.0/sqrt(test_vec.size());
        }
        else{
            test_vec[i]=-1.0/sqrt(test_vec.size());
        }
    }

    // print should be +-1/12
    for(int i=0; i<(int)test_vec.size(); i++){
        if(i%12==11){
            cout<<setw(12)<<test_vec[i]<<endl;
        }
        else{
            cout<<setw(12)<<test_vec[i]<<" ";
        }
    }

    // call the randomizer function
    //R(test_vec, 0.5);

    cout<<endl;

    for(int i=0; i<(int)test_vec.size(); i++){
        if(i%12==11){
            cout<<setw(12)<<test_vec[i]<<endl;
        }
        else{
            cout<<setw(12)<<test_vec[i]<<" ";
        }
    }
    cout<<"test randomizer function end"<<endl;
}

int test_code(){
    for(int val=0; val<20; val++){
        vector<int> enc;
        double dec;
        enc=code(val);
        dec=decode(enc);
        cout<<val<<" "<<dec<<endl;
        cout<<endl;
    }

    return 0;
}

int test_PROT_PP_S_Hist_PP(dataSet &data, double epsilon){
    int real=0;
    int i=0;
    while(1){
        if(data.get_dataMtx(i,0)!=0){
            real=data.get_dataMtx(i, 0);
            break;
        }
        i++;
    }
    pair<double, double> ret;
    int n=100;
    double sum=0;
    double freq_sum=0.0;
    for(int j=0; j<5; j++){
        if(j==0){epsilon=0.05;}
        if(j==1){epsilon=0.1;}
        if(j==2){epsilon=0.5;}
        if(j==3){epsilon=1;}
        if(j==4){epsilon=5;}
        sum=0.0; freq_sum=0.0;
        for(int i=0; i<n; i++){
            ret=PROT_PP_S_Hist_pp(data, epsilon);
            sum=sum+ret.first;
            freq_sum=freq_sum+ret.second;
            if(i%20==0){
                cout<<"iteration: "<<i<<" desired value: "<<real<<" desired frequency: 0.3"<<endl;
                cout<<"m: "<<data.get_m()<<" n: "<<data.get_n()<<" epsilon: "<<epsilon<<endl;
                cout<<"heavy hitter running average: "<<(double)sum/(i+1)<<endl;
                cout<<"freq running average: "<<(double)freq_sum/(i+1)<<endl;
                cout<<endl;
            }
        }
        cout<<"avgerage final: "<<(double)sum/n<<endl;
    }

    return 0;
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

    cout<<"parse end"<<endl;

    int m=test_data->get_m();
    int n=test_data->get_n();
    cout<<"m (users): "<<m<<" n (attributes): "<<n<<endl;

    cout<<endl;

    // test_basic_randomizer();
    // test_code();
    test_PROT_PP_S_Hist_PP(*test_data, 0.05);

    delete test_data;

    return 0;
}

// int main(int argc, const char* argv[]){
//     if(argc!=4){
//         cout<<"usage: ./test path_to_data initial_m initial_n"<<endl;
//         return 0;
//     }
//     dataSet *test_data=new(dataSet);
//
//     cout<<"begin parsing data"<<endl;
//
//     // parse data
//     int lol;
//     cout<<stoi(argv[2])<<endl;
//     cout<<stoi(argv[3])<<endl;
//     lol=parseData(argv[1], *test_data, stoi(argv[2]), stoi(argv[3]));
//     if(lol==-1){
//         cout<<"failed to parse data you might want to check the path"<<endl;
//         return -1;
//     }
//
//     printall(*test_data);
//
//     cout<<"parse end"<<endl;
//
//     int m=test_data->get_m();
//     int n=test_data->get_n();
//     cout<<"m (users): "<<m<<" n (attributes): "<<n<<endl;
//
//     // calculate real averages
//     vector<double> avg_real(test_data->get_n(), 0.0);
//     vector<int> asdf(test_data->get_n(), 0);
//
//     for(int i=0; i<test_data->get_m(); i++){
//         for(int j=0; j<test_data->get_n(); j++){
//             if(test_data->get_dataMtx(i,j)!=0){
//                 avg_real[j]=avg_real[j]+test_data->get_dataMtx(i, j);
//                 asdf[j]++;
//             }
//         }
//     }
//
//     for(int i=0; i<test_data->get_n(); i++){
//         avg_real[i]=avg_real[i]/asdf[i];
//     }
//
//     // data normalization
//     cout<<"normalizing columns"<<endl;
//     vector<double> outlier={0};
//     pair<vector<double>, vector<double>> scales=normalizeNeg1toPos1(*test_data, outlier);
//
//     cout<<"done normalizing"<<endl;
//
//     printall(*test_data);
//     cout<<endl;
//
//     // perturb data
//     vector<double> perturb_scales=tuplePerturbation(*test_data,0.5);
//     printall(*test_data);
//
//     // calculate estimated average
//     vector<double> avg(test_data->get_n(), 0.0);
//     vector<int> counters(test_data->get_n(), 0);
//
//     for(int j=0; j<test_data->get_n(); j++){
//         for(int i=0; i<test_data->get_m(); i++){
//             if(test_data->get_dataMtx(i,j)!=0){
//                 avg[j]=avg[j]+test_data->get_dataMtx(i,j);
//                 counters[j]++;
//             }
//         }
//     }
//
//     for(int i=0; i<test_data->get_n(); i++){
//         avg[i]=avg[i]/counters[i];
//     }
//
//     // undo normalizaion
//     undoNorm(scales, avg);
//
//     // print real averages and estimated averages
//     double rel_error=0.0;
//     for(int i=0; i<test_data->get_n(); i++){
//          cout<<"avg_real: "<<avg_real[i]<<setw(10)<<" avg_est: "<<avg[i]<<endl;
//          rel_error+=abs((avg[i]-avg_real[i])/avg_real[i]);
//     }
//     cout<<"percent releative error: "<<rel_error/test_data->get_n()<<endl;
//
//     delete test_data;
//     return 0;
// }
