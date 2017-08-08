#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "data.h"
#include "util.h"

using namespace std;
#define LARGENUMBER 999999999

void printall(dataSet &data){
    int m=data.get_m();
    int n=data.get_n();

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            cout<<setw(12)<<data.get_dataMtx(i,j);
        }
        cout<<endl;
        cout<<"next row"<<endl;
    }
}

double MSE(vector<double> avg_real, vector<double> avg_est){
    int n=avg_real.size();
    double total=0.0;
    for(int i=0; i<n; i++){
        cout<<"total MSE: "<<total<<endl;
        total=total+(avg_real[i]-avg_est[i])*(avg_real[i]-avg_est[i]);
    }
    return (double)total/n;
}

void test_bassily_smith(dataSet &data, double epsilon, double beta, int k){

}


int main(int argc, const char* argv[]){
    if(argc<4 || argc>5){
        cout<<"usage: ./test path_to_data initial_m initial_n"<<endl;
	cout<<"usage: ./test path_to_data initial_m initial_n epsilon"<<endl;
        cout<<argc<<endl;
        return 0;
    }

    dataSet *test_data=new(dataSet);
    // parse data
    int lol;
    // cout<<stoi(argv[2])<<endl;
    // cout<<stoi(argv[3])<<endl;
    cout<<"Parse data"<<endl;
    lol=parseData(argv[1], *test_data, stoi(argv[2]), stoi(argv[3]));
    if(lol==-1){
        cout<<"failed to parse data you might want to check the path"<<endl;
        return -1;
    }

    cout<<"parse end"<<endl;

    int m=test_data->get_m();
    int n=test_data->get_n();

    cout<<"m: "<<m<<endl;
    cout<<"n: "<<n<<endl;

    cout<<endl;

    vector<double> res=bassily_smith(*test_data, 0.8, 0.3, 1000000);
    double sum=0;
    for(int i=0; i<(int)res.size(); i++){
        cout<<res[i]<<endl;
        sum+=res[i];
    }
    cout<<"sum: "<<sum<<endl;

    // cout<<endl;
    //
    // vector<double> res1=bassily_smith(*test_data, 5, 0.8, 100000);
    // sum=0;
    // for(int i=0; i<(int)res.size(); i++){
    //     //cout<<res[i]<<endl;
    //     sum+=res1[i];
    // }
    // cout<<"sum: "<<sum<<endl;
    //
    // cout<<endl;
    //
    // vector<double> res2=bassily_smith(*test_data, 2.5, 0.3, 100000);
    // sum=0;
    // for(int i=0; i<(int)res.size(); i++){
    //     //cout<<res[i]<<endl;
    //     sum+=res2[i];
    // }
    // cout<<"sum: "<<sum<<endl;
    //
    // cout<<endl;
    //
    // vector<double> res3=bassily_smith(*test_data, 2.5, 0.8, 100000);
    // sum=0;
    // for(int i=0; i<(int)res.size(); i++){
    //     //cout<<res[i]<<endl;
    //     sum+=res3[i];
    // }
    // cout<<"sum: "<<sum<<endl;
    //
    // cout<<endl;
    //
    // vector<double> res4=bassily_smith(*test_data, 5, 0.3, 100000);
    // sum=0;
    // for(int i=0; i<(int)res.size(); i++){
    //     //cout<<res[i]<<endl;
    //     sum+=res4[i];
    // }
    // cout<<"sum: "<<sum<<endl;
    //
    // cout<<endl;
    //
    // vector<double> res5=bassily_smith(*test_data, 5, 0.8, 100000);
    // sum=0;
    // for(int i=0; i<(int)res.size(); i++){
    //     //cout<<res[i]<<endl;
    //     sum+=res5[i];
    // }
    // cout<<"sum: "<<sum<<endl;
    //
    // cout<<endl;

    return 0;
}
